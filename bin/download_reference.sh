#!/bin/bash
# download_reference.sh
#
# Downloads and prepares the hg38 reference FASTA (bgzipped, main
# chromosomes only, no alt scaffolds, no chrM) and the matching FlashFry
# SpCas9 NGG off-target index from Zenodo (DOI 10.5281/zenodo.21114489),
# verifies checksums, decompresses, and builds the .fai/.gzi indices
# needed by bedtools/samtools.
#
# Usage:
#   ./bin/download_reference.sh --outdir ~/.guidenet_reference
#
# After this completes, point local_paths.env at:
#   FASTA=<outdir>/hg38.fa.gz
#   FLASHFRY_DB=<outdir>/hg38_standard_cas9ngg_db
#
# Requires on PATH: curl or wget, sha256sum, tar, bgzip (from htslib/samtools).

set -euo pipefail

ZENODO_RECORD_ID="21114489"
BASE_URL="https://zenodo.org/records/${ZENODO_RECORD_ID}/files"

FASTA_GZ="hg38.fa.gz"
FASTA_SHA="hg38.fa.gz.sha256"
DB_TAR="hg38_standard_cas9ngg_db.tar.gz"
DB_SHA="hg38_standard_cas9ngg_db.tar.gz.sha256"

OUTDIR="${HOME}/.guidenet_reference"
SKIP_EXISTING=1

while [[ $# -gt 0 ]]; do
    case "$1" in
        --outdir) OUTDIR="$2"; shift 2 ;;
        --force) SKIP_EXISTING=0; shift ;;
        -h|--help)
            grep '^#' "$0" | sed 's/^#//'
            exit 0
            ;;
        *) echo "ERROR: unknown argument $1" >&2; exit 1 ;;
    esac
done

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

download() {
    local fname="$1"
    if [[ -f "$fname" && "$SKIP_EXISTING" -eq 1 ]]; then
        echo "  -> ${fname} already present, skipping download (use --force to re-download)" >&2
        return
    fi
    echo "  -> downloading ${fname}" >&2
    if command -v curl >/dev/null 2>&1; then
        curl -L -o "${fname}.part" "${BASE_URL}/${fname}?download=1"
    else
        wget -O "${fname}.part" "${BASE_URL}/${fname}?download=1"
    fi
    mv "${fname}.part" "${fname}"
}

verify() {
    local fname="$1" sha_fname="$2"
    echo "  -> verifying ${fname}" >&2
    local expected actual
    expected=$(awk '{print $1}' "${sha_fname}")
    actual=$(sha256sum "${fname}" | awk '{print $1}')
    if [[ "$expected" != "$actual" ]]; then
        echo "ERROR: checksum mismatch for ${fname}" >&2
        echo "  expected: ${expected}" >&2
        echo "  actual:   ${actual}" >&2
        echo "  Re-run with --force to redownload." >&2
        exit 1
    fi
}

echo "=== Downloading reference data (Zenodo DOI 10.5281/zenodo.${ZENODO_RECORD_ID}) to ${OUTDIR} ===" >&2

echo "-- FASTA --" >&2
download "${FASTA_SHA}"
download "${FASTA_GZ}"
verify "${FASTA_GZ}" "${FASTA_SHA}"

echo "-- FlashFry index --" >&2
download "${DB_SHA}"
download "${DB_TAR}"
verify "${DB_TAR}" "${DB_SHA}"

echo "=== Preparing files ===" >&2

if [[ ! -f "${FASTA_GZ}.fai" || ! -f "${FASTA_GZ}.gzi" || "$SKIP_EXISTING" -eq 0 ]]; then
    echo "  -> building .fai/.gzi index for ${FASTA_GZ} (samtools faidx)" >&2
    if ! command -v samtools >/dev/null 2>&1; then
        echo "ERROR: samtools not found on PATH -- required to index the bgzipped FASTA" >&2
        exit 1
    fi
    samtools faidx "${FASTA_GZ}"
else
    echo "  -> ${FASTA_GZ}.fai/.gzi already present, skipping" >&2
fi

DB_DIR_NAME="hg38_standard_cas9ngg_db"
if [[ ! -d "${DB_DIR_NAME}" || "$SKIP_EXISTING" -eq 0 ]]; then
    echo "  -> extracting ${DB_TAR}" >&2
    tar xzf "${DB_TAR}"
else
    echo "  -> ${DB_DIR_NAME}/ already present, skipping extraction" >&2
fi

echo "" >&2
echo "Done. Set these in local_paths.env:" >&2
echo "  FASTA=${OUTDIR}/${FASTA_GZ}" >&2
echo "  FLASHFRY_DB=${OUTDIR}/${DB_DIR_NAME}" >&2
