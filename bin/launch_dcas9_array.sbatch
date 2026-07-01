#!/bin/bash
#SBATCH -J launch_dcas9
#
# Resource directives (partition, mem, cpus, time, gres) are NOT hardcoded
# here -- they're passed in via `sbatch` CLI flags from run_standalone.sh,
# built from the cluster-specific block in local_paths.env. This keeps this
# template identical across clusters; only local_paths.env changes.
#
# Generic launch-dCas9 scoring array job. Expects these to be set via
# `sbatch --export=ALL,VAR=value,...` (see run_standalone.sh):
#   OUTCOME        enhancerFitness | promoterFitness
#   CHUNK_DIR      directory containing chunk_XXXX.csv (from launch_dcas9_split.py)
#   RESULT_DIR     where to write chunk_XXXX.csv results for this outcome
#   LAUNCH_SCRIPT  path to launch-dCas9.py
#   CKPT_DIR       path to pretrained checkpoints
#   VENV_PYTHON    path to the launch-dCas9 conda env's python
#   WRAPPER_SCRIPT path to launch_dcas9_score_chunk.py (this repo's bin/)
#   LOG_DIR        where SBATCH -o/-e would go (informational; actual
#                   redirection is set by the sbatch --output/--error flags
#                   passed at submit time, not hardcoded here)

set -euo pipefail

TASK_ID=$(printf "%04d" "${SLURM_ARRAY_TASK_ID}")
CHUNK_CSV="${CHUNK_DIR}/chunk_${TASK_ID}.csv"
TMP_DIR="${RESULT_DIR}/tmp_${TASK_ID}"
RESULT_CSV="${RESULT_DIR}/chunk_${TASK_ID}.csv"

echo "Starting ${OUTCOME} chunk ${TASK_ID} on $(hostname) at $(date)"

"${VENV_PYTHON}" "${WRAPPER_SCRIPT}" \
    --chunk_csv "${CHUNK_CSV}" \
    --outcome "${OUTCOME}" \
    --launch_script "${LAUNCH_SCRIPT}" \
    --ckpt_dir "${CKPT_DIR}" \
    --tmp_dir "${TMP_DIR}" \
    --result_csv "${RESULT_CSV}" \
    --python_bin "${VENV_PYTHON}"

rm -rf "${TMP_DIR}"
echo "Done ${OUTCOME} chunk ${TASK_ID} at $(date)"
