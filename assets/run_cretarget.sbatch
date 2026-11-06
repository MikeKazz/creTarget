#!/bin/bash
#SBATCH --job-name=cretarget
#SBATCH --partition=gpu
#SBATCH --nodes=1
#
# IMPORTANT -- the next two values are deliberately tiny.
#
# This job is only the Nextflow *driver*: it submits the real work as separate
# SLURM jobs and then sits waiting for them. On a single-node cluster like
# boltzmann, every resource this job holds is a resource the actual work cannot
# have. Request 8 CPUs and 80 GB here (as a normal analysis job would) and you
# reserve most of the machine to do nothing, starving the jobs you are waiting
# on -- in the worst case the driver waits forever for work that can never be
# scheduled, which is a deadlock, not a slow run.
#
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#
# No --gres here. The driver needs no GPU; the LAUNCH_DCAS9 tasks request their
# own via the boltzmann profile.
#
#SBATCH --time=240:00:00
#SBATCH --output=/home/jmorrislab1/mkazemi/___ProJECTS___/creTarget/logs/cretarget_%j.out
#SBATCH --error=/home/jmorrislab1/mkazemi/___ProJECTS___/creTarget/logs/cretarget_%j.err

set -euo pipefail

# ---------------------------------------------------------------------------
# Two filesystems, on purpose.
#
#   CODE_DIR  (jmorrislab1) -- code, inputs, logs. Small, backed up, the thing
#                              you edit and commit.
#   DATA_DIR  (jmorrislab0) -- reference genome, FlashFry index, model
#                              checkpoints, the Nextflow work directory, and
#                              results. This is where the TB-scale storage is.
#
# Sizes that justify the split: the reference + FlashFry index is ~3.2 GB, the
# Nextflow work directory runs to hundreds of GB on a genome-wide run, and the
# collated output TSV is ~63 GB for 2.35M cCREs (~132M guide rows).
# ---------------------------------------------------------------------------
CODE_DIR="/home/jmorrislab1/mkazemi/___ProJECTS___/creTarget"
DATA_DIR="/home/jmorrislab0/mkazemi/creTarget"

export CRETARGET_CACHE_DIR="${DATA_DIR}/cache"   # reference, index, checkpoints
export NXF_WORK="${DATA_DIR}/work"               # Nextflow intermediates

# Nextflow's own JVM. The default heap is generous for a driver that is mostly
# waiting; capping it keeps the driver inside its 4 GB allocation.
export NXF_OPTS='-Xms500M -Xmx2G'

# Keep Nextflow's assets/plugins cache off $HOME too -- it is small, but it
# avoids surprises if $HOME has a quota.
export NXF_HOME="${DATA_DIR}/.nextflow_home"
export CRETARGET_CONFIG_DIR="${CODE_DIR}/.cretarget_config"

mkdir -p "${CODE_DIR}/logs" "${DATA_DIR}"/{cache,work,results}

# A batch job gets a NON-interactive shell, where `conda activate` fails
# because the conda initialize block in ~/.bashrc never runs. Sourcing
# conda.sh first is what makes activation work here. Verified on boltzmann:
# bare `conda activate` FAILS, `source conda.sh && conda activate` WORKS.
CONDA_SH="/opt/miniconda3/etc/profile.d/conda.sh"
if [ -r "$CONDA_SH" ]; then
    # shellcheck disable=SC1090
    . "$CONDA_SH"
fi

# Activate whichever environment creTarget is installed into. On boltzmann
# this is a plain venv, not a conda env, so it is sourced directly.
PYENV="/home/jmorrislab1/mkazemi/myenv"
if [ -r "${PYENV}/bin/activate" ]; then
    # shellcheck disable=SC1091
    . "${PYENV}/bin/activate"
fi

# Reference data. These already exist on boltzmann, so `cretarget
# setup-reference` is unnecessary here -- it is for a fresh machine.
#
# The FASTA is used only by bedtools to extract region sequence; the FlashFry
# database is used for the off-target search. They need not be the same file,
# but they MUST be the same assembly. This pairing (hg38.fa with the
# hg38_standard index) is the one the pre-Nextflow pipeline already used.
REF_FASTA="/home/jmorrislab1/mkazemi/.guidenet_reference/hg38.fa"
FLASHFRY_JAR="/home/jmorrislab1/mkazemi/bin/FlashFry-assembly-1.15.jar"
FLASHFRY_DB="/home/jmorrislab1/mkazemi/.guidenet_reference/hg38_standard_cas9ngg_db"

# Fail early and loudly rather than part-way through a long run.
command -v cretarget >/dev/null || { echo "ERROR: cretarget not on PATH" >&2; exit 1; }
command -v nextflow  >/dev/null || { echo "ERROR: nextflow not on PATH"  >&2; exit 1; }
[ -r "$REF_FASTA" ]        || { echo "ERROR: no FASTA at $REF_FASTA" >&2; exit 1; }
[ -r "${REF_FASTA}.fai" ]  || { echo "ERROR: no FASTA index at ${REF_FASTA}.fai -- run: samtools faidx $REF_FASTA" >&2; exit 1; }
[ -r "$FLASHFRY_JAR" ]     || { echo "ERROR: no FlashFry jar at $FLASHFRY_JAR" >&2; exit 1; }
[ -r "$FLASHFRY_DB" ]      || { echo "ERROR: no FlashFry index at $FLASHFRY_DB" >&2; exit 1; }
[ -r "${FLASHFRY_DB}.header" ] || { echo "ERROR: no FlashFry index header at ${FLASHFRY_DB}.header" >&2; exit 1; }

cd "${CODE_DIR}"

# Run the pipeline from the local clone. Switch to `MikeKazz/creTarget` (and
# add `-r v0.2.0`) once you want Nextflow to pull a tagged release from GitHub
# instead -- pinning a revision is what makes a run reproducible.
nextflow run "${CODE_DIR}/creTarget" \
    -profile boltzmann \
    -resume \
    --enhancer_bed "${CODE_DIR}/inputs/enhancers_subset.txt" \
    --fasta        "${REF_FASTA}" \
    --flashfry_jar "${FLASHFRY_JAR}" \
    --flashfry_db  "${FLASHFRY_DB}" \
    --outdir       "${DATA_DIR}/results"

# To add the optional steps, append:
#
#   --run_launch_dcas9 \
#   --launch_script /path/to/launch-dCas9/launch-dCas9.py \
#   --launch_python /path/to/conda/envs/<torch-env>/bin/python \
#   --ckpt_dir      "${CRETARGET_CACHE_DIR}/models" \
#
#   --run_seed_finder \
#   --seed_finder_dir /path/to/seed-finder-lite \
#   --seed_finder_min_length 8 \
#   --seed_finder_tss_window 1000
