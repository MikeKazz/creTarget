#!/usr/bin/env bash
# bump_version.sh -- set or verify the release version everywhere it appears.
#
#   ./bump_version.sh 0.2.1        # rewrite every pin, then verify
#   ./bump_version.sh --check      # verify that the pins agree with pyproject
#   ./bump_version.sh --check 0.2.1
#                                  # verify that the pins agree with 0.2.1,
#                                  # e.g. the tag that triggered a release
#
# The version lives in several files. Bumping them by hand is how a release
# ends up with a pyproject saying 0.2.1, a container tag saying 0.2.0, and a
# conda environment installing a third version. This sets all of them and then
# verifies they agree.
#
# Historical references (CHANGELOG entries, comments describing which version
# fixed what, and version-specific DOIs) are deliberately NOT touched.
set -euo pipefail

cd "$(dirname "$0")"

usage() {
    echo "usage: $0 <major.minor.patch>" >&2
    echo "       $0 --check [major.minor.patch]" >&2
    exit 1
}

# BSD sed (macOS) requires an explicit empty backup suffix; GNU sed rejects one.
# Detect once so the caller can use the same syntax on both.
if sed --version >/dev/null 2>&1; then
    sed_inplace() { sed -i "$@"; }
else
    sed_inplace() { sed -i '' "$@"; }
fi

CURRENT=$(sed -n 's/^version = "\(.*\)"/\1/p' pyproject.toml | head -1)
if [ -z "$CURRENT" ]; then
    echo "ERROR: could not read current version from pyproject.toml" >&2
    exit 1
fi

CHECK_ONLY=0
NEW=""
case "${1:-}" in
    --check|-c)
        CHECK_ONLY=1
        NEW="${2:-$CURRENT}"
        ;;
    "")
        usage
        ;;
    -h|--help)
        usage
        ;;
    *)
        NEW="$1"
        ;;
esac

if ! printf '%s' "$NEW" | grep -Eq '^[0-9]+\.[0-9]+\.[0-9]+$'; then
    echo "ERROR: '$NEW' is not a semantic version (expected MAJOR.MINOR.PATCH)" >&2
    exit 1
fi

OLD="$CURRENT"
OLD_RE="${OLD//./\\.}"

if [ "$CHECK_ONLY" -eq 1 ]; then
    if [ "$NEW" = "$OLD" ]; then
        echo "checking version references against ${NEW}"
    else
        echo "checking version references against ${NEW} (pyproject says ${OLD})"
    fi
else
    if [ "$OLD" = "$NEW" ]; then
        echo "version is already ${NEW}; verifying references"
    else
        echo "bumping ${OLD} -> ${NEW}"
        sed_inplace "s/^version = \"${OLD_RE}\"/version = \"${NEW}\"/" \
            pyproject.toml
        sed_inplace "s/^__version__ = \"${OLD_RE}\"/__version__ = \"${NEW}\"/" \
            src/cretarget/__init__.py
        sed_inplace "s/^version: ${OLD_RE}\$/version: ${NEW}/" \
            CITATION.cff
        sed_inplace "s|archive/refs/tags/v${OLD_RE}\.tar\.gz|archive/refs/tags/v${NEW}.tar.gz|" \
            assets/environment.yml
        sed_inplace "s|cretarget:${OLD_RE}|cretarget:${NEW}|g" \
            Dockerfile nextflow.config README.md
        sed_inplace "s|image\.version=\"${OLD_RE}\"|image.version=\"${NEW}\"|" \
            Dockerfile
        sed_inplace "s|^    version         = '${OLD_RE}'|    version         = '${NEW}'|" \
            nextflow.config
        sed_inplace "s|-r v${OLD_RE}|-r v${NEW}|g" \
            README.md assets/run_cretarget.sbatch
        sed_inplace "s|archive/refs/tags/v${OLD_RE}\.tar\.gz|archive/refs/tags/v${NEW}.tar.gz|g" \
            README.md

        # Documentation and wiki examples must not keep running the previous
        # release. Only the current-release patterns are rewritten, so a
        # historical note or a version-specific DOI is left alone.
        while IFS= read -r file; do
            sed_inplace \
                -e "s|cretarget:${OLD_RE}|cretarget:${NEW}|g" \
                -e "s|cretarget: ${OLD_RE}|cretarget: ${NEW}|g" \
                -e "s|creTarget ${OLD_RE}|creTarget ${NEW}|g" \
                -e "s|-r v${OLD_RE}|-r v${NEW}|g" \
                -e "s|tags/v${OLD_RE}\.tar\.gz|tags/v${NEW}.tar.gz|g" \
                -e "s|\"version\": \"${OLD_RE}\"|\"version\": \"${NEW}\"|g" \
                -e "s|\"revision\": \"v${OLD_RE}\"|\"revision\": \"v${NEW}\"|g" \
                "$file"
        done < <(grep -rlF "$OLD" docs wiki 2>/dev/null || true)
    fi
fi

echo
echo "=== verification ==="
fail=0
check() { # label expected actual
    if [ "$2" = "$3" ]; then printf '  OK   %-28s %s\n' "$1" "$3"
    else printf '  FAIL %-28s got %s, want %s\n' "$1" "$3" "$2"; fail=1; fi
}
check pyproject.toml      "$NEW" "$(sed -n 's/^version = "\(.*\)"/\1/p' pyproject.toml | head -1)"
check __init__.py         "$NEW" "$(sed -n 's/^__version__ = "\(.*\)"/\1/p' src/cretarget/__init__.py | head -1)"
check CITATION.cff        "$NEW" "$(sed -n 's/^version: //p' CITATION.cff | head -1)"
check nextflow.manifest   "$NEW" "$(sed -n "s/^    version         = '\(.*\)'/\1/p" nextflow.config | head -1)"
check container_image     "$NEW" "$(sed -n "s|.*cretarget:\([0-9.]*\)'.*|\1|p" nextflow.config | head -1)"
check environment.yml     "$NEW" "$(sed -n 's|.*refs/tags/v\([0-9.]*\)\.tar\.gz|\1|p' assets/environment.yml | head -1)"
check Dockerfile.label    "$NEW" "$(sed -n 's|.*image\.version="\([0-9.]*\)".*|\1|p' Dockerfile | head -1)"
check sbatch.revision     "$NEW" "$(sed -n 's|.*-r v\([0-9.]*\) .*|\1|p' assets/run_cretarget.sbatch | head -1)"

# Documentation and wiki examples are part of what users copy. Check every
# current-release pattern, not just the canonical pins above.
current_refs() { # pattern
    grep -rhoE -e "$1" docs wiki README.md 2>/dev/null \
        | sed -E 's/.*(cretarget:|v)([0-9]+\.[0-9]+\.[0-9]+).*/\2/' | sort -u || true
}
for pattern in 'cretarget:[0-9]+\.[0-9]+\.[0-9]+' \
               '-r v[0-9]+\.[0-9]+\.[0-9]+' \
               'tags/v[0-9]+\.[0-9]+\.[0-9]+\.tar\.gz'; do
    refs=$(current_refs "$pattern" | tr '\n' ' ')
    if [ -z "$refs" ]; then
        printf '  OK   %-28s %s\n' "docs examples ($pattern)" "none"
    elif [ "$refs" = "${NEW} " ]; then
        printf '  OK   %-28s %s\n' "docs examples ($pattern)" "${NEW}"
    else
        printf '  FAIL %-28s got %s, want %s\n' \
            "docs examples ($pattern)" "$refs" "$NEW"
        fail=1
    fi
done

if [ "$CHECK_ONLY" -eq 0 ] && [ "$OLD" != "$NEW" ]; then
    # Catch any reference the targeted edits above did not know about. A stale
    # pin in a template is worse than a stale comment: it silently runs the old
    # code. -F because the version's dots are regex wildcards otherwise:
    # "0.2.5" would then match "0" + anything + "2" + anything + "5" and report
    # files that never mentioned the old version. A check that cries wolf gets
    # ignored.
    stale=$(grep -rlnF "${OLD}" \
                --exclude-dir=.git --exclude-dir=__pycache__ \
                --exclude-dir='*.egg-info' --exclude-dir=site \
                --exclude-dir=dist --exclude-dir=build \
                --exclude-dir=.venv --exclude-dir=work --exclude-dir=results \
                --exclude-dir='results_*' \
                --exclude=CHANGELOG.md --exclude=bump_version.sh \
                --exclude=.nextflow.log --exclude='*.log' \
                . 2>/dev/null || true)
    if [ -n "$stale" ]; then
        printf '  WARN %-28s %s\n' "still mention ${OLD}:" "$(echo "$stale" | tr '\n' ' ')"
        echo "       check each one -- a historical comment or a version-specific"
        echo "       DOI is fine, a current-release pin is not."
    fi
fi

echo
if [ "$fail" -ne 0 ]; then
    if [ "$CHECK_ONLY" -eq 1 ]; then
        echo "VERSION REFERENCES DO NOT AGREE on ${NEW} -- fix before tagging." >&2
    else
        echo "SOME REFERENCES DID NOT UPDATE -- fix before tagging." >&2
    fi
    exit 1
fi
if [ "$CHECK_ONLY" -eq 1 ]; then
    echo "All version references agree on ${NEW}."
    exit 0
fi
echo "All version references agree on ${NEW}."
echo
echo "Next:"
echo "  1. move the CHANGELOG.md Unreleased entries under [${NEW}]"
echo "  2. review and commit the version bump"
echo "  3. merge to main, then tag: git tag -a v${NEW} -m \"creTarget v${NEW}\""
echo "  4. git push origin v${NEW} -- the release workflow publishes the"
echo "     GitHub Release, PyPI distribution, and container image."
