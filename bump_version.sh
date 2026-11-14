#!/usr/bin/env bash
# bump_version.sh -- set the release version everywhere it appears.
#
#   ./bump_version.sh 0.2.1
#
# The version lives in nine places across six files. Bumping them by hand is
# how a release ends up with a pyproject saying 0.2.1, a container tag saying
# 0.2.0, and a conda environment installing a third version. This sets all of
# them and then verifies they agree.
#
# Historical references (CHANGELOG entries, comments describing which version
# fixed what) are deliberately NOT touched.
set -euo pipefail

NEW="${1:-}"
if ! printf '%s' "$NEW" | grep -Eq '^[0-9]+\.[0-9]+\.[0-9]+$'; then
    echo "usage: $0 <major.minor.patch>   e.g. $0 0.2.1" >&2
    exit 1
fi

cd "$(dirname "$0")"
OLD=$(sed -n 's/^version = "\(.*\)"/\1/p' pyproject.toml)
[ -n "$OLD" ] || { echo "ERROR: could not read current version from pyproject.toml" >&2; exit 1; }
echo "bumping ${OLD} -> ${NEW}"

sed -i "s/^version = \"${OLD}\"/version = \"${NEW}\"/"                     pyproject.toml
sed -i "s/^__version__ = \"${OLD}\"/__version__ = \"${NEW}\"/"             src/cretarget/__init__.py
sed -i "s/^version: ${OLD}\$/version: ${NEW}/"                             CITATION.cff
sed -i "s|archive/refs/tags/v${OLD}\.tar\.gz|archive/refs/tags/v${NEW}.tar.gz|" assets/environment.yml
sed -i "s|cretarget:${OLD}|cretarget:${NEW}|g"                             Dockerfile nextflow.config README.md
sed -i "s|image\.version=\"${OLD}\"|image.version=\"${NEW}\"|"             Dockerfile
sed -i "s|^    version         = '${OLD}'|    version         = '${NEW}'|" nextflow.config
sed -i "s|-r v${OLD}|-r v${NEW}|g"              README.md assets/run_cretarget.sbatch

echo
echo "=== verification ==="
fail=0
check() { # label expected actual
    if [ "$2" = "$3" ]; then printf '  OK   %-28s %s\n' "$1" "$3"
    else printf '  FAIL %-28s got %s, want %s\n' "$1" "$3" "$2"; fail=1; fi
}
check pyproject.toml      "$NEW" "$(sed -n 's/^version = "\(.*\)"/\1/p' pyproject.toml)"
check __init__.py         "$NEW" "$(sed -n 's/^__version__ = "\(.*\)"/\1/p' src/cretarget/__init__.py)"
check CITATION.cff        "$NEW" "$(sed -n 's/^version: //p' CITATION.cff)"
check nextflow.manifest   "$NEW" "$(sed -n "s/^    version         = '\(.*\)'/\1/p" nextflow.config)"
check container_image     "$NEW" "$(sed -n "s|.*cretarget:\([0-9.]*\)'.*|\1|p" nextflow.config)"
check environment.yml     "$NEW" "$(sed -n 's|.*refs/tags/v\([0-9.]*\)\.tar\.gz|\1|p' assets/environment.yml)"
check Dockerfile.label    "$NEW" "$(sed -n 's|.*image\.version="\([0-9.]*\)".*|\1|p' Dockerfile)"
check sbatch.revision     "$NEW" "$(sed -n 's|.*-r v\([0-9.]*\) .*|\1|p' assets/run_cretarget.sbatch)"

# Catch any reference the targeted seds above did not know about. A stale pin
# in a template is worse than a stale comment: it silently runs the old code.
# -F because the version's dots are regex wildcards otherwise: "0.2.5" then
# matches "0" + anything + "2" + anything + "5" and reports files that never
# mentioned the old version. A check that cries wolf gets ignored.
stale=$(grep -rlnF "${OLD}" --exclude-dir=.git --exclude-dir=__pycache__ \
            --exclude-dir=*.egg-info --exclude=CHANGELOG.md \
            --exclude=bump_version.sh --exclude=.nextflow.log . 2>/dev/null || true)
if [ -n "$stale" ]; then
    printf '  WARN %-28s %s\n' "still mention ${OLD}:" "$(echo "$stale" | tr '\n' ' ')"
    echo "       check each one -- a historical comment is fine, a pin is not."
fi

echo
if [ "$fail" -ne 0 ]; then
    echo "SOME REFERENCES DID NOT UPDATE -- fix before tagging." >&2
    exit 1
fi
echo "All version references agree on ${NEW}."
echo
echo "Next:"
echo "  1. add a CHANGELOG.md entry for ${NEW}"
echo "  2. git commit -am \"Release v${NEW}\" && git push origin main"
echo "  3. docker build -t ghcr.io/mikekazz/cretarget:${NEW} . && docker push ghcr.io/mikekazz/cretarget:${NEW}"
echo "  4. git tag -a v${NEW} -m \"creTarget v${NEW}\" && git push origin v${NEW}"
echo "  5. publish a GitHub Release from the tag (Zenodo only mints a DOI on a Release)"
