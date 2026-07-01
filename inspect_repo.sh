#!/bin/bash
# inspect_repo.sh
#
# Run this BEFORE `git init` / first commit. Checks for the kinds of things
# that are annoying or risky to fix after they're already in git history:
#   - files too large to commit (should go to a GitHub Release instead)
#   - stray/leftover artifacts (e.g. the {modules incident)
#   - bin/ scripts missing the executable bit
#   - local_paths.env present without a corresponding .example template
#   - local_paths.env not in .gitignore (real cluster paths shouldn't be
#     force-pushed to a repo that's going public later)
#   - results/ or other run-output directories that shouldn't be committed
#   - empty directories git won't track anyway (informational only)
#
# Exits non-zero if anything CRITICAL is found (large files, missing
# .gitignore entries). Warnings are printed but don't fail the run.

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

LARGE_FILE_THRESHOLD_MB=50
EXIT_CODE=0

echo "Inspecting: ${SCRIPT_DIR}"
echo ""

# ── 1. Large files (would bloat git history / exceed GitHub's 100MB hard limit) ──
echo "== Large files (>${LARGE_FILE_THRESHOLD_MB}MB) =="
FOUND_LARGE=0
while IFS= read -r -d '' f; do
    size_mb=$(du -m "$f" | cut -f1)
    echo "  CRITICAL: ${f} (${size_mb}MB) -- should be a GitHub Release asset, not committed"
    FOUND_LARGE=1
done < <(find . -type f -size +"${LARGE_FILE_THRESHOLD_MB}"M -not -path "./.git/*" -print0 2>/dev/null)
if [[ "$FOUND_LARGE" -eq 0 ]]; then
    echo "  OK -- none found"
else
    EXIT_CODE=1
fi
echo ""

# ── 2. Stray/malformed artifacts (leading {, }, or other shell-expansion debris) ──
echo "== Stray artifacts (shell-expansion leftovers, odd filenames) =="
FOUND_STRAY=0
while IFS= read -r -d '' f; do
    echo "  CRITICAL: ${f} -- looks like leftover shell-expansion debris"
    FOUND_STRAY=1
done < <(find . -not -path "./.git/*" \( -name "*{*" -o -name "*}*" \) -print0 2>/dev/null)
if [[ "$FOUND_STRAY" -eq 0 ]]; then
    echo "  OK -- none found"
else
    EXIT_CODE=1
fi
echo ""

# ── 3. bin/ scripts missing executable bit ──
echo "== bin/ executable permissions =="
FOUND_NONEXEC=0
if [[ -d bin ]]; then
    while IFS= read -r -d '' f; do
        if [[ ! -x "$f" ]]; then
            echo "  WARNING: ${f} is not executable (chmod +x)"
            FOUND_NONEXEC=1
        fi
    done < <(find bin -maxdepth 1 -type f \( -name "*.py" -o -name "*.sh" \) -print0 2>/dev/null)
fi
if [[ "$FOUND_NONEXEC" -eq 0 ]]; then
    echo "  OK -- all executable"
fi
echo ""

# ── 4. local_paths.env vs .example template ──
echo "== local_paths.env / template =="
if [[ -f local_paths.env && ! -f local_paths.env.example ]]; then
    echo "  WARNING: local_paths.env exists but no local_paths.env.example -- create the template before committing"
fi
if [[ -f local_paths.env.example ]]; then
    if grep -qE '/home/jmorrislab1|/home/[a-z]+/' local_paths.env.example 2>/dev/null; then
        echo "  WARNING: local_paths.env.example appears to contain a real absolute user path -- should be a generic placeholder"
    else
        echo "  OK -- local_paths.env.example looks generic"
    fi
fi
echo ""

# ── 5. .gitignore coverage ──
echo "== .gitignore coverage =="
if [[ ! -f .gitignore ]]; then
    echo "  CRITICAL: no .gitignore found"
    EXIT_CODE=1
else
    for pattern in "local_paths.env" "results/" "*.fai" "*.log"; do
        if ! grep -qxF "$pattern" .gitignore; then
            echo "  WARNING: .gitignore does not contain '${pattern}'"
        fi
    done
    echo "  (see .gitignore contents below for full picture)"
fi
echo ""

# ── 6. Leftover run outputs that shouldn't be committed ──
echo "== Run-output directories present (should be gitignored, not committed) =="
for d in results test_data/prep_out test_data/standalone_out test_data/launch_split_test; do
    if [[ -d "$d" ]]; then
        echo "  NOTE: ${d}/ exists locally -- fine to keep for your own testing, just confirm it's gitignored"
    fi
done
echo ""

echo "== Summary =="
if [[ "$EXIT_CODE" -eq 0 ]]; then
    echo "No CRITICAL issues found. Review WARNINGs above, then proceed with git init."
else
    echo "CRITICAL issues found above -- resolve before running 'git add .' / first commit."
fi

exit "$EXIT_CODE"
