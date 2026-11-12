#!/bin/bash
#
# Driver submission template for a SLURM cluster.
#
#   1. copy this file next to your inputs
#   2. edit the CONFIGURE block below
#   3. sbatch run_cretarget.sbatch
#
#SBATCH --job-name=cretarget
#SBATCH --nodes=1
#
# IMPORTANT -- the next two values are deliberately tiny.
#
# This job is only the Nextflow *driver*: it submits the real work as separate
# scheduler jobs and then waits. On a single-node cluster every resource this
# job holds is a resource the actual work cannot have. Request 8 CPUs and
# 80 GB here, as a normal analysis job would, and you reserve most of the
# machine to do nothing -- starving the jobs you are waiting on, and in the
# worst case deadlocking: waiting forever for work that can never be
# scheduled.
#
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=240:00:00
#
# No --gres here. The driver needs no GPU; the LAUNCH_DCAS9 tasks request
# their own through the site profile.
#
# Uncomment and set if your site requires them:
##SBATCH --partition=YOUR_PARTITION
##SBATCH --account=YOUR_ACCOUNT
#
#SBATCH --output=logs/cretarget_%j.out
#SBATCH --error=logs/cretarget_%j.err

set -euo pipefail

# ===========================================================================
# CONFIGURE -- everything you need to edit is in this block
# ===========================================================================

# Where the pipeline, inputs and logs live. Small; keep it backed up.
CODE_DIR="/path/to/your/project"

# Where reference data, the Nextflow work directory and results go. Large:
# budget ~300 GB of work directory and ~50 GB of output for a genome-wide run.
DATA_DIR="/path/to/large/storage"

# Reference data. `cretarget setup-reference` downloads the first and third.
# FlashFry is a separate JVM tool:
#   https://github.com/mckennalab/FlashFry/releases
REF_FASTA="${DATA_DIR}/cache/reference/hg38.fa"
FLASHFRY_JAR="${DATA_DIR}/tools/FlashFry-assembly-1.15.jar"
FLASHFRY_DB="${DATA_DIR}/cache/reference/hg38_standard_cas9ngg_db"

# Your targets: BED, or one "chr:start-end" per line.
TARGETS="${CODE_DIR}/inputs/targets.bed"

# The site profile to run under. `singlenode` is the worked example in
# conf/example-singlenode.config; `slurm` is the generic one.
PROFILE="singlenode"

# How creTarget is installed. Leave both empty if it is already on PATH.
CONDA_SH=""          # e.g. /opt/miniconda3/etc/profile.d/conda.sh
PYTHON_ENV=""        # a conda env name, or a path to a venv to activate

# ===========================================================================
# Below here you should not need to edit
# ===========================================================================

export CRETARGET_CACHE_DIR="${DATA_DIR}/cache"
export CRETARGET_CONFIG_DIR="${CODE_DIR}/.cretarget_config"
export NXF_WORK="${DATA_DIR}/work"
export NXF_HOME="${DATA_DIR}/.nextflow_home"

# The driver's own JVM. The default heap is generous for a process that mostly
# waits; capping it keeps the driver inside its 4 GB allocation.
export NXF_OPTS='-Xms500M -Xmx2G'

mkdir -p "${CODE_DIR}/logs" "${DATA_DIR}"/{cache,work,results}
cd "${CODE_DIR}"

# A batch job gets a NON-interactive shell, where `conda activate` fails
# because the conda initialize block in ~/.bashrc never runs. Sourcing
# conda.sh first is what makes activation work here.
if [ -n "${CONDA_SH}" ] && [ -r "${CONDA_SH}" ]; then
    # shellcheck disable=SC1090
    . "${CONDA_SH}"
fi
if [ -n "${PYTHON_ENV}" ]; then
    if [ -r "${PYTHON_ENV}/bin/activate" ]; then
        # shellcheck disable=SC1091
        . "${PYTHON_ENV}/bin/activate"      # a venv
    else
        conda activate "${PYTHON_ENV}"      # a conda env
    fi
fi

# Fail early and loudly rather than part-way through a long run.
command -v cretarget >/dev/null || { echo "ERROR: cretarget not on PATH" >&2; exit 1; }
command -v nextflow  >/dev/null || { echo "ERROR: nextflow not on PATH"  >&2; exit 1; }
[ -r "$TARGETS" ]              || { echo "ERROR: no targets at $TARGETS" >&2; exit 1; }
[ -r "$REF_FASTA" ]            || { echo "ERROR: no FASTA at $REF_FASTA" >&2; exit 1; }
[ -r "${REF_FASTA}.fai" ]      || { echo "ERROR: no FASTA index at ${REF_FASTA}.fai -- run: samtools faidx $REF_FASTA" >&2; exit 1; }
[ -r "$FLASHFRY_JAR" ]         || { echo "ERROR: no FlashFry jar at $FLASHFRY_JAR" >&2; exit 1; }
[ -r "$FLASHFRY_DB" ]          || { echo "ERROR: no FlashFry index at $FLASHFRY_DB" >&2; exit 1; }
[ -r "${FLASHFRY_DB}.header" ] || { echo "ERROR: no index header at ${FLASHFRY_DB}.header" >&2; exit 1; }

nextflow run MikeKazz/creTarget -r v0.2.5 \
    -profile "${PROFILE}" \
    -resume \
    --enhancer_bed "${TARGETS}" \
    --fasta        "${REF_FASTA}" \
    --flashfry_jar "${FLASHFRY_JAR}" \
    --flashfry_db  "${FLASHFRY_DB}" \
    --outdir       "${DATA_DIR}/results"

# Optional steps, appended as needed:
#
#   --snp_bed /path/to/snps.txt          # "chr:pos" per line, 1-based
#   --atac_narrowpeak /path/to/peaks.narrowPeak
#   --dedup_blast_radius
#
#   --run_launch_dcas9 \
#   --launch_script /path/to/launch-dCas9.py \
#   --launch_python /path/to/torch-env/bin/python \
#   --ckpt_dir      "${CRETARGET_CACHE_DIR}/models" \
#
#   --run_seed_finder \
#   --seed_finder_dir /path/to/seed-finder-lite \
#   --seed_finder_min_length 10 \
#   --seed_finder_tss_window 1000
