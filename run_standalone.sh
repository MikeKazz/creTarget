#!/bin/bash
# run_standalone.sh
#
# Standalone (non-Nextflow) driver for the crispri-guide-designer MVP logic.
#
# HOW TO RUN THIS (every time, no matter where you unzip it):
#   1. cd into the folder this script is in (e.g. the unzipped crispri-guide-designer/)
#   2. ./run_standalone.sh --snp_bed <path> --enhancer_bed <path> [options]
#
# You do NOT need to pass --fasta / --flashfry_jar / --flashfry_db / launch-dCas9
# paths -- those are read automatically from local_paths.env (sitting next to
# this script). Edit local_paths.env once if those cluster paths ever change.
#
# All relative paths below (bin/, results/, local_paths.env) resolve
# relative to THIS SCRIPT'S location, not your current directory.
#
# Usage:
#   ./run_standalone.sh \
#       [--snp_bed snps.bed] [--enhancer_bed enhancers.bed] \
#       [--fasta /path] [--flashfry_jar /path] [--flashfry_db /path] \
#       [--outdir results] [--flank_bp 50] [--snp_window_bp 150] \
#       [--blast_radius_bp 1000] [--dedup_blast_radius] \
#       [--max_mismatch 3] \
#       [--run_launch_dcas9] [--launch_chunk_size 200000] \
#       [--launch_venv_python /path] [--launch_script /path] [--launch_ckpt_dir /path]
#
# At least one of --snp_bed / --enhancer_bed is required.
# --snp_bed / --enhancer_bed paths, if relative, resolve relative to YOUR
# CURRENT DIRECTORY (not the script location).
#
# --run_launch_dcas9 adds a step (between annotate and rank) that submits two
# SLURM array jobs (enhancerFitness, promoterFitness) via sbatch, waits for
# them to finish, and merges both scores onto every guide before ranking.
# Requires a SLURM login/head node and the launch-dCas9 conda env to exist.
#
# Requires on PATH: python3 (+ pandas), bedtools, samtools, java, sbatch/squeue
# (only if --run_launch_dcas9 is used).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN_DIR="${SCRIPT_DIR}/bin"
ENV_FILE="${SCRIPT_DIR}/local_paths.env"

# ── Load default cluster paths (edit local_paths.env, not this file) ──────
FASTA=""
FLASHFRY_JAR=""
FLASHFRY_DB=""
LAUNCH_VENV_PYTHON=""
LAUNCH_SCRIPT=""
LAUNCH_CKPT_DIR=""
LAUNCH_SBATCH_PARTITION=""
LAUNCH_SBATCH_GRES=""
LAUNCH_SBATCH_MEM="16G"
LAUNCH_SBATCH_CPUS="8"
LAUNCH_SBATCH_TIME="2:00:00"
if [[ -f "$ENV_FILE" ]]; then
    # shellcheck disable=SC1090
    source "$ENV_FILE"
fi
# Default checkpoints to the ones bundled in this repo, if not overridden.
# Resolved relative to THIS SCRIPT'S location (not cwd), same as everything else.
: "${LAUNCH_CKPT_DIR:=${SCRIPT_DIR}/checkpoints/}"

# ── Other defaults ──────────────────────────────────────────────────────
OUTDIR="${SCRIPT_DIR}/results"
FLANK_BP=50
SNP_WINDOW_BP=150
BLAST_RADIUS_BP=1000
DEDUP_FLAG=""
ATAC_NARROWPEAK=""
MAX_MISMATCH=3
SNP_BED=""
ENHANCER_BED=""
SCORING_METRICS="doench2014ontarget,doench2016cfd,dangerous,hsu2013,moreno2015,jostandsantos,minot,folding"
RUN_LAUNCH_DCAS9=0
LAUNCH_CHUNK_SIZE=200000
LAUNCH_POLL_INTERVAL=30

# ── Arg parsing (CLI flags override local_paths.env / defaults above) ────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --fasta) FASTA="$2"; shift 2 ;;
        --flashfry_jar) FLASHFRY_JAR="$2"; shift 2 ;;
        --flashfry_db) FLASHFRY_DB="$2"; shift 2 ;;
        --snp_bed) SNP_BED="$2"; shift 2 ;;
        --enhancer_bed) ENHANCER_BED="$2"; shift 2 ;;
        --outdir) OUTDIR="$2"; shift 2 ;;
        --flank_bp) FLANK_BP="$2"; shift 2 ;;
        --snp_window_bp) SNP_WINDOW_BP="$2"; shift 2 ;;
        --blast_radius_bp) BLAST_RADIUS_BP="$2"; shift 2 ;;
        --dedup_blast_radius) DEDUP_FLAG="--dedup_blast_radius"; shift ;;
        --atac_narrowpeak) ATAC_NARROWPEAK="$2"; shift 2 ;;
        --max_mismatch) MAX_MISMATCH="$2"; shift 2 ;;
        --scoring_metrics) SCORING_METRICS="$2"; shift 2 ;;
        --run_launch_dcas9) RUN_LAUNCH_DCAS9=1; shift ;;
        --launch_chunk_size) LAUNCH_CHUNK_SIZE="$2"; shift 2 ;;
        --launch_venv_python) LAUNCH_VENV_PYTHON="$2"; shift 2 ;;
        --launch_script) LAUNCH_SCRIPT="$2"; shift 2 ;;
        --launch_ckpt_dir) LAUNCH_CKPT_DIR="$2"; shift 2 ;;
        -h|--help)
            grep '^#' "$0" | sed 's/^#//'
            exit 0
            ;;
        *) echo "ERROR: unknown argument $1" >&2; exit 1 ;;
    esac
done

# ── Validation ────────────────────────────────────────────────────────────
if [[ -z "$FASTA" || -z "$FLASHFRY_JAR" || -z "$FLASHFRY_DB" ]]; then
    echo "ERROR: FASTA/FLASHFRY_JAR/FLASHFRY_DB not set. Either fill them into" >&2
    echo "       ${ENV_FILE}, or pass --fasta/--flashfry_jar/--flashfry_db explicitly." >&2
    exit 1
fi
if [[ -z "$SNP_BED" && -z "$ENHANCER_BED" ]]; then
    echo "ERROR: at least one of --snp_bed / --enhancer_bed is required." >&2
    exit 1
fi
if [[ ! -f "$FASTA" ]]; then
    echo "ERROR: FASTA not found at ${FASTA} (check ${ENV_FILE})" >&2
    exit 1
fi
if [[ ! -f "${FASTA}.fai" ]]; then
    echo "No .fai found for ${FASTA} -- running samtools faidx..." >&2
    samtools faidx "${FASTA}"
fi
if [[ "$RUN_LAUNCH_DCAS9" -eq 1 ]]; then
    if [[ -z "$LAUNCH_VENV_PYTHON" || -z "$LAUNCH_SCRIPT" || -z "$LAUNCH_CKPT_DIR" ]]; then
        echo "ERROR: --run_launch_dcas9 requires LAUNCH_VENV_PYTHON/LAUNCH_SCRIPT/LAUNCH_CKPT_DIR" >&2
        echo "       set in ${ENV_FILE} or passed via --launch_venv_python/--launch_script/--launch_ckpt_dir" >&2
        exit 1
    fi
    if ! command -v sbatch >/dev/null 2>&1; then
        echo "ERROR: --run_launch_dcas9 requires sbatch/squeue on PATH (run from a SLURM login/head node)" >&2
        exit 1
    fi
fi

mkdir -p "${OUTDIR}"

echo "Using:" >&2
echo "  FASTA:        ${FASTA}" >&2
echo "  FlashFry JAR: ${FLASHFRY_JAR}" >&2
echo "  FlashFry DB:  ${FLASHFRY_DB}" >&2
echo "  Output dir:   ${OUTDIR}" >&2
if [[ "$RUN_LAUNCH_DCAS9" -eq 1 ]]; then
    echo "  launch-dCas9: ENABLED (venv=${LAUNCH_VENV_PYTHON})" >&2
fi
echo "" >&2

TOTAL_STEPS=5
if [[ "$RUN_LAUNCH_DCAS9" -eq 1 ]]; then TOTAL_STEPS=6; fi

echo "=== [1/${TOTAL_STEPS}] Preparing regions ===" >&2

FASTA_PARTS=()
META_PARTS=()

if [[ -n "$SNP_BED" ]]; then
    mkdir -p "${OUTDIR}/prep_snp"
    python3 "${BIN_DIR}/prep_regions.py" \
        --bed "${SNP_BED}" \
        --input_type snp \
        --fasta "${FASTA}" \
        --snp_window_bp "${SNP_WINDOW_BP}" \
        --link_distance_bp "${BLAST_RADIUS_BP}" \
        ${DEDUP_FLAG} \
        --outdir "${OUTDIR}/prep_snp"
    FASTA_PARTS+=("${OUTDIR}/prep_snp/regions.fasta")
    META_PARTS+=("${OUTDIR}/prep_snp/regions.metadata.tsv")
fi

if [[ -n "$ENHANCER_BED" ]]; then
    mkdir -p "${OUTDIR}/prep_enhancer"
    python3 "${BIN_DIR}/prep_regions.py" \
        --bed "${ENHANCER_BED}" \
        --input_type enhancer \
        --fasta "${FASTA}" \
        --flank_bp "${FLANK_BP}" \
        --link_distance_bp "${BLAST_RADIUS_BP}" \
        ${DEDUP_FLAG} \
        --outdir "${OUTDIR}/prep_enhancer"
    FASTA_PARTS+=("${OUTDIR}/prep_enhancer/regions.fasta")
    META_PARTS+=("${OUTDIR}/prep_enhancer/regions.metadata.tsv")
fi

cat "${FASTA_PARTS[@]}" > "${OUTDIR}/regions.combined.fasta"
head -n 1 "${META_PARTS[0]}" > "${OUTDIR}/regions.combined.metadata.tsv"
for f in "${META_PARTS[@]}"; do
    tail -n +2 "$f" >> "${OUTDIR}/regions.combined.metadata.tsv"
done

N_REGIONS=$(grep -c '^>' "${OUTDIR}/regions.combined.fasta" || true)
echo "  -> ${N_REGIONS} total regions -> ${OUTDIR}/regions.combined.fasta" >&2

echo "=== [2/${TOTAL_STEPS}] FlashFry discover ===" >&2
java -Xmx14g -jar "${FLASHFRY_JAR}" \
    discover \
    --database "${FLASHFRY_DB}" \
    --fasta "${OUTDIR}/regions.combined.fasta" \
    --output "${OUTDIR}/discover.output" \
    --maxMismatch "${MAX_MISMATCH}"

echo "=== [3/${TOTAL_STEPS}] FlashFry score ===" >&2
java -Xmx7g -jar "${FLASHFRY_JAR}" \
    score \
    --input "${OUTDIR}/discover.output" \
    --output "${OUTDIR}/scored.output" \
    --scoringMetrics "${SCORING_METRICS}" \
    --database "${FLASHFRY_DB}" \
    --maxMismatch "${MAX_MISMATCH}"

echo "=== [4/${TOTAL_STEPS}] Annotating (merge + distance/blast-radius flagging) ===" >&2
ATAC_ARGS=()
if [[ -n "$ATAC_NARROWPEAK" ]]; then
    ATAC_ARGS=(--atac_narrowpeak "$ATAC_NARROWPEAK")
fi
python3 "${BIN_DIR}/annotate.py" \
    --scored_chunks "${OUTDIR}/scored.output" \
    --region_metadata "${OUTDIR}/regions.combined.metadata.tsv" \
    --blast_radius_bp "${BLAST_RADIUS_BP}" \
    "${ATAC_ARGS[@]}" \
    --out "${OUTDIR}/guides.annotated.tsv"

ANNOTATED_FOR_RANK="${OUTDIR}/guides.annotated.tsv"

if [[ "$RUN_LAUNCH_DCAS9" -eq 1 ]]; then
    echo "=== [5/${TOTAL_STEPS}] launch-dCas9 scoring (enhancerFitness + promoterFitness) ===" >&2
    LAUNCH_DIR="${OUTDIR}/launch_dcas9"
    CHUNK_DIR="${LAUNCH_DIR}/chunks"
    ENH_RESULT_DIR="${LAUNCH_DIR}/results_enhancer"
    PRO_RESULT_DIR="${LAUNCH_DIR}/results_promoter"
    LAUNCH_LOG_DIR="${LAUNCH_DIR}/logs"
    mkdir -p "${CHUNK_DIR}" "${ENH_RESULT_DIR}" "${PRO_RESULT_DIR}" "${LAUNCH_LOG_DIR}"

    echo "  -- splitting unique spacers into chunks --" >&2
    python3 "${BIN_DIR}/launch_dcas9_split.py" \
        --annotated "${ANNOTATED_FOR_RANK}" \
        --chunk_size "${LAUNCH_CHUNK_SIZE}" \
        --outdir "${CHUNK_DIR}"

    N_CHUNKS=$(cat "${CHUNK_DIR}/n_chunks.txt")
    ARRAY_RANGE="0-$((N_CHUNKS - 1))"
    echo "  -- submitting SLURM arrays (${ARRAY_RANGE}) --" >&2

    # Build sbatch flags from local_paths.env's cluster-specific block.
    # Any left blank there are simply omitted (e.g. no partition, no GPU).
    SBATCH_FLAGS=(--parsable --array="${ARRAY_RANGE}" --mem="${LAUNCH_SBATCH_MEM}" --cpus-per-task="${LAUNCH_SBATCH_CPUS}" --time="${LAUNCH_SBATCH_TIME}")
    [[ -n "$LAUNCH_SBATCH_PARTITION" ]] && SBATCH_FLAGS+=(-p "${LAUNCH_SBATCH_PARTITION}")
    [[ -n "$LAUNCH_SBATCH_GRES" ]] && SBATCH_FLAGS+=(--gres="${LAUNCH_SBATCH_GRES}")

    ENH_JOBID=$(sbatch "${SBATCH_FLAGS[@]}" \
        -o "${LAUNCH_LOG_DIR}/enhancer_%A_%a.log" \
        -e "${LAUNCH_LOG_DIR}/enhancer_%A_%a.err" \
        --export=ALL,OUTCOME=enhancerFitness,CHUNK_DIR="${CHUNK_DIR}",RESULT_DIR="${ENH_RESULT_DIR}",LAUNCH_SCRIPT="${LAUNCH_SCRIPT}",CKPT_DIR="${LAUNCH_CKPT_DIR}",VENV_PYTHON="${LAUNCH_VENV_PYTHON}",WRAPPER_SCRIPT="${BIN_DIR}/launch_dcas9_score_chunk.py" \
        "${BIN_DIR}/launch_dcas9_array.sbatch")

    PRO_JOBID=$(sbatch "${SBATCH_FLAGS[@]}" \
        -o "${LAUNCH_LOG_DIR}/promoter_%A_%a.log" \
        -e "${LAUNCH_LOG_DIR}/promoter_%A_%a.err" \
        --export=ALL,OUTCOME=promoterFitness,CHUNK_DIR="${CHUNK_DIR}",RESULT_DIR="${PRO_RESULT_DIR}",LAUNCH_SCRIPT="${LAUNCH_SCRIPT}",CKPT_DIR="${LAUNCH_CKPT_DIR}",VENV_PYTHON="${LAUNCH_VENV_PYTHON}",WRAPPER_SCRIPT="${BIN_DIR}/launch_dcas9_score_chunk.py" \
        "${BIN_DIR}/launch_dcas9_array.sbatch")

    echo "  -- submitted enhancer job ${ENH_JOBID}, promoter job ${PRO_JOBID}; waiting --" >&2
    while true; do
        PENDING=$(squeue -j "${ENH_JOBID},${PRO_JOBID}" -h 2>/dev/null | wc -l)
        if [[ "$PENDING" -eq 0 ]]; then
            break
        fi
        echo "  ... ${PENDING} array task(s) still running, checking again in ${LAUNCH_POLL_INTERVAL}s" >&2
        sleep "${LAUNCH_POLL_INTERVAL}"
    done

    echo "  -- verifying all chunks completed --" >&2
    MISSING=0
    for i in $(seq 0 $((N_CHUNKS - 1))); do
        TASK_ID=$(printf "%04d" "$i")
        if [[ ! -s "${ENH_RESULT_DIR}/chunk_${TASK_ID}.csv" ]]; then
            echo "  MISSING enhancer chunk ${TASK_ID}" >&2
            MISSING=1
        fi
        if [[ ! -s "${PRO_RESULT_DIR}/chunk_${TASK_ID}.csv" ]]; then
            echo "  MISSING promoter chunk ${TASK_ID}" >&2
            MISSING=1
        fi
    done
    if [[ "$MISSING" -eq 1 ]]; then
        echo "ERROR: one or more launch-dCas9 chunks did not complete. Check logs in ${LAUNCH_LOG_DIR}," >&2
        echo "       then resubmit the missing indices manually (see launch_dcas9_array.sbatch) and rerun" >&2
        echo "       with the same --outdir to pick up where you left off (chunks already present are reused)." >&2
        exit 1
    fi

    echo "  -- merging scores onto guide table --" >&2
    python3 "${BIN_DIR}/launch_dcas9_merge.py" \
        --annotated "${ANNOTATED_FOR_RANK}" \
        --enhancer_dir "${ENH_RESULT_DIR}" \
        --promoter_dir "${PRO_RESULT_DIR}" \
        --out "${OUTDIR}/guides.annotated.with_launch.tsv"

    ANNOTATED_FOR_RANK="${OUTDIR}/guides.annotated.with_launch.tsv"
fi

echo "=== [${TOTAL_STEPS}/${TOTAL_STEPS}] Ranking (hard filters + weighted composite score) ===" >&2
python3 "${BIN_DIR}/rank_report.py" \
    --annotated "${ANNOTATED_FOR_RANK}" \
    --max_mismatch "${MAX_MISMATCH}" \
    --out "${OUTDIR}/guides.ranked.tsv"

echo "" >&2
echo "Done. Final output: ${OUTDIR}/guides.ranked.tsv" >&2
