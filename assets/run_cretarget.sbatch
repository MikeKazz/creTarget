#!/bin/bash
#
# Driver submission template for a SLURM cluster.
#
#   1. copy this file next to your inputs
#   2. edit the CONFIGURE block below
#   3. mkdir -p logs          <-- REQUIRED before the first submission
#   4. sbatch run_cretarget.sbatch
#
# `mkdir -p logs` is step 3 and not something this script can do for you:
# slurmd opens the --output/--error files below when the job launches, long
# before the first line of the script runs. If the directory is missing the
# job dies at launch with no log to explain why. The paths are relative to
# the directory you run `sbatch` FROM.
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
#
# ADJUST: must be <= your partition's MaxTime, or sbatch rejects the job
# immediately with "Requested time limit is invalid". Check yours with:
#     sinfo -o '%P %l'
# The driver only waits on other jobs, so this needs to span the whole run.
#SBATCH --time=48:00:00
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

# The site profile to run under.
#
# `slurm` is the generic one and the right default. Do NOT use `singlenode`
# unless you have read conf/example-singlenode.config: it is a worked example
# of one specific machine and hardcodes that site's partition ('gpu'), its
# GPU request and its memory ceilings, none of which will match yours.
#
# Append a software profile if you are not using a host install, e.g.
# "slurm,singularity" or "slurm,conda".
PROFILE="slurm"

# Partition/queue and allocation code. Leave empty if your site has a usable
# default and does not require an account.
#   sinfo -s                                  # partitions
#   sacctmgr show assoc where user=$USER      # accounts
QUEUE=""             # e.g. "compute"
ACCOUNT=""           # e.g. "proj1234"

# Which pipeline to run. A local clone works on compute nodes with no
# outbound network; "MikeKazz/creTarget" makes Nextflow fetch it from GitHub
# at runtime, which needs internet access from wherever this job lands.
PIPELINE="${CODE_DIR}/creTarget"
REVISION=""          # e.g. -r v0.5.1 -- only used when fetching remotely

# If your site provides tools as environment modules, list them here. Names
# and versions differ at every site -- check `module avail java nextflow`.
SITE_MODULES=""      # e.g. "java/17.0.8 nextflow/24.10.0 singularity/3.11"

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

if [ -n "${SITE_MODULES}" ]; then
    if command -v module >/dev/null 2>&1 || [ -n "${MODULESHOME:-}" ]; then
        # Unquoted on purpose: SITE_MODULES is a list of module names.
        # shellcheck disable=SC2086
        module load ${SITE_MODULES}
    else
        echo "WARNING: SITE_MODULES is set but no 'module' command exists" >&2
    fi
fi

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

# Built as an array so an empty QUEUE/ACCOUNT contributes no argument at all
# rather than an empty string the scheduler would reject.
EXTRA_ARGS=()
[ -n "${QUEUE}" ] && EXTRA_ARGS+=(--queue "${QUEUE}")
if [ -n "${ACCOUNT}" ]; then
    # A real file, not a process substitution: Nextflow may read a -c config
    # more than once, and a fifo can only be read to the end once.
    ACCOUNT_CONFIG="${CODE_DIR}/.cretarget_account.config"
    printf "process.clusterOptions = '--account=%s'\n" "${ACCOUNT}" > "${ACCOUNT_CONFIG}"
    EXTRA_ARGS+=(-c "${ACCOUNT_CONFIG}")
fi
# shellcheck disable=SC2086
nextflow run "${PIPELINE}" ${REVISION} \
    -profile "${PROFILE}" \
    -resume \
    "${EXTRA_ARGS[@]}" \
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
#   --seed_finder_min_length 10 \
#   --seed_finder_tss_window 1000
#
# The seed screen needs its promoter reference built once, beforehand:
#   cretarget build-seed-reference --fasta "$REF_FASTA"
