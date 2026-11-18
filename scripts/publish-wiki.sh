#!/usr/bin/env bash
#
# Publish the Markdown files in wiki/ to the GitHub wiki.
#
#   scripts/publish-wiki.sh --dry-run    # show what would change, change nothing
#   scripts/publish-wiki.sh              # prompt, then commit and push
#   scripts/publish-wiki.sh --yes        # commit and push without prompting
#
# GitHub keeps a wiki in a separate git repository. It does not exist until the
# first page is created in the web UI, so this script cannot initialize it.
# If the clone fails, create any page at
# https://github.com/MikeKazz/creTarget/wiki and run the script again.
#
# The script creates and pushes a backup/<UTC timestamp> branch pointing at the
# wiki's previous HEAD before it pushes. The command to restore that snapshot is
# printed at the end, so a bad sync is easy to undo.
set -euo pipefail

REPO_ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$REPO_ROOT"

WIKI_SRC="$REPO_ROOT/wiki"
REMOTE="${CRETARGET_WIKI_REMOTE:-https://github.com/MikeKazz/creTarget.wiki.git}"
DRY_RUN=0
ASSUME_YES=0

usage() {
    cat <<'EOF'
usage: scripts/publish-wiki.sh [--dry-run] [--yes] [--remote URL]

  --dry-run    Clone the wiki, copy the pages, show the diff, then stop.
  --yes, -y    Do not prompt before committing and pushing.
  --remote URL Wiki remote (default: $CRETARGET_WIKI_REMOTE or the
               MikeKazz/creTarget wiki).
  --help, -h   Show this help.
EOF
}

while [ "$#" -gt 0 ]; do
    case "$1" in
        --dry-run) DRY_RUN=1 ;;
        --yes|-y)  ASSUME_YES=1 ;;
        --remote)
            [ "$#" -ge 2 ] || { echo "ERROR: --remote needs a URL" >&2; exit 2; }
            REMOTE="$2"
            shift
            ;;
        --help|-h) usage; exit 0 ;;
        *) echo "ERROR: unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
    shift
done

if [ ! -d "$WIKI_SRC" ]; then
    echo "ERROR: $WIKI_SRC does not exist" >&2
    exit 2
fi
if [ ! -f "$WIKI_SRC/Home.md" ]; then
    echo "ERROR: $WIKI_SRC/Home.md is missing" >&2
    exit 2
fi

tmp=$(mktemp -d)
cleanup() { rm -rf "$tmp"; }
trap cleanup EXIT

if ! git clone --quiet "$REMOTE" "$tmp/wiki" 2>"$tmp/clone.err"; then
    cat "$tmp/clone.err" >&2
    echo >&2
    echo "ERROR: could not clone the wiki repository." >&2
    echo "GitHub does not create <repo>.wiki.git until the first page exists." >&2
    echo "Create any page at https://github.com/MikeKazz/creTarget/wiki," >&2
    echo "then run this script again." >&2
    exit 2
fi

cd "$tmp/wiki"
branch=$(git symbolic-ref --short HEAD)
if before=$(git rev-parse --verify HEAD 2>/dev/null); then
    :
else
    before=""
fi

# The wiki repository is flat: page names are file names at its root.
for file in "$WIKI_SRC"/*.md; do
    cp "$file" .
done

git add -A .

if git diff --cached --quiet; then
    echo "The wiki already matches $WIKI_SRC; nothing to publish."
    exit 0
fi

echo "Changes for $REMOTE (branch $branch):"
git status --short
echo
git diff --cached --stat
echo

if [ "$DRY_RUN" -eq 1 ]; then
    echo "Dry run: nothing was committed or pushed."
    echo "Run without --dry-run to publish these changes."
    exit 0
fi

if [ "$ASSUME_YES" -ne 1 ]; then
    printf 'Commit and push these changes? [y/N] '
    read -r answer
    case "$answer" in
        y|Y|yes|YES) ;;
        *) echo "Aborted; nothing was pushed."; exit 1 ;;
    esac
fi

source_rev=$(git -C "$REPO_ROOT" rev-parse --short HEAD 2>/dev/null || echo "unknown")
backup=""
if [ -n "$before" ]; then
    backup="backup/$(date -u +%Y%m%d-%H%M%S)"
    git branch "$backup" "$before"
    git push --quiet origin "$backup"
    echo "Backup branch pushed: $backup"
fi

git commit --quiet -m "Sync wiki from wiki/ (source ${source_rev})"
git push --quiet origin "HEAD:${branch}"

echo "Published to $REMOTE (branch $branch)."
echo
if [ -n "$before" ]; then
    echo "To undo this sync:"
    echo "  git clone $REMOTE /tmp/cretarget-wiki-undo"
    echo "  cd /tmp/cretarget-wiki-undo"
    echo "  git reset --hard ${before}"
    echo "  git push --force origin ${branch}"
    echo
    echo "The same snapshot is on branch ${backup} in the wiki repository."
else
    echo "The wiki had no commits before this sync, so there is no previous"
    echo "state to restore. The new commit is ${branch}'s first."
fi
