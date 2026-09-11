#!/bin/bash
#SBATCH --job-name=ccre_batches
#SBATCH --cpus-per-task=40
#SBATCH --mem=110G
#SBATCH --time=48:00:00
#SBATCH -o /home/jmorrislab1/mkazemi/___ProJECTS___/Blair/gRNA_design_tool/sandbox/crispri-guide-designer/results_ccre_batches/slurm_%j.log
#SBATCH -e /home/jmorrislab1/mkazemi/___ProJECTS___/Blair/gRNA_design_tool/sandbox/crispri-guide-designer/results_ccre_batches/slurm_%j.err
#
# Runs all cCRE BED batches through `crispri-guide-designer run`, ONE AT A
# TIME (sequential, not parallel across batches -- deliberate, per the
# memory pressure seen in the chr1 pilot). Resume-safe: skips any batch
# whose guides.ranked.tsv already exists, so a resubmit after a crash
# picks up where it left off instead of redoing completed batches.
#
# Submit with: sbatch run_ccre_batches.sbatch

set -uo pipefail  # NOT -e: a single batch failing shouldn't kill the whole job;
                   # we track failures and report them at the end instead.

REPO_DIR="/home/jmorrislab1/mkazemi/___ProJECTS___/Blair/gRNA_design_tool/sandbox/crispri-guide-designer"
BATCH_DIR="${REPO_DIR}/test_data/ccre_batches"
OUT_BASE="${REPO_DIR}/results_ccre_batches"
mkdir -p "${OUT_BASE}"

cd "${REPO_DIR}"

FAILED_BATCHES=()

for batch_file in $(ls "${BATCH_DIR}"/*.bed | sort); do
    batch_name=$(basename "${batch_file}" .bed)
    outdir="${OUT_BASE}/${batch_name}"
    ranked_file="${outdir}/guides.ranked.tsv"

    if [[ -s "${ranked_file}" ]]; then
        echo "=== ${batch_name}: already completed (${ranked_file} exists), skipping ==="
        continue
    fi

    echo "=== Starting ${batch_name} at $(date) ==="
    crispri-guide-designer run \
        --enhancer-bed "${batch_file}" \
        --atac-cell-type k562 \
        --dedup-blast-radius \
        --outdir "${outdir}"
    rc=$?
    echo "=== Finished ${batch_name} at $(date), exit code ${rc} ==="

    if [[ ${rc} -ne 0 ]]; then
        echo "ERROR: ${batch_name} failed (exit ${rc})"
        FAILED_BATCHES+=("${batch_name}")
    fi
done

echo ""
echo "=== All batches processed ==="
if [[ ${#FAILED_BATCHES[@]} -gt 0 ]]; then
    echo "FAILED batches (${#FAILED_BATCHES[@]}): ${FAILED_BATCHES[*]}"
    echo "Re-run this same sbatch script to retry only the failed/incomplete batches (completed ones are skipped)."
    exit 1
else
    echo "All batches completed successfully."
fi
