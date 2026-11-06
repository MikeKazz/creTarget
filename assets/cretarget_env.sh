# cretarget_env.sh -- project-local environment, sourced deliberately.
#
#   cd /home/jmorrislab1/mkazemi/___ProJECTS___/creTarget
#   source cretarget_env.sh
#
# Nothing here is permanent. It is NOT meant for ~/.bashrc: these settings
# would then apply to every shell and every other project, which on a shared
# machine is how one project quietly breaks another. Source it when you work
# on creTarget; open a new shell and it is gone.
#
# Every variable is CRETARGET_-prefixed or already namespaced (NXF_*), so
# nothing here can collide with another project's generic CODE_DIR/DATA_DIR.
#
# The sbatch driver script sets all of this itself, in its own process, so a
# submitted job never depends on your interactive shell state.

# --- Paths ------------------------------------------------------------------
# Small and backed up: code, inputs, logs.
export CRETARGET_CODE_DIR="/home/jmorrislab1/mkazemi/___ProJECTS___/creTarget"
# Large: reference data, Nextflow work directory, results.
export CRETARGET_DATA_DIR="/home/jmorrislab0/mkazemi/creTarget"

# --- creTarget's own caches -------------------------------------------------
# Downloadable assets (model checkpoints, ATAC peaks, reference) go here
# rather than ~/.cretarget, keeping 3+ GB off the smaller filesystem.
export CRETARGET_CACHE_DIR="${CRETARGET_DATA_DIR}/cache"

# Keep the step-level config inside the project too, so `cretarget configure`
# does not write to ~/.config and cannot be confused with another project's
# settings. Unset this and the tool falls back to ~/.config/cretarget.
export CRETARGET_CONFIG_DIR="${CRETARGET_CODE_DIR}/.cretarget_config"

# --- Nextflow ---------------------------------------------------------------
# Intermediates. Hundreds of GB on a genome-wide run, so this must be on the
# large filesystem.
export NXF_WORK="${CRETARGET_DATA_DIR}/work"

# Nextflow's own assets/plugins cache. Scoped to the project so it cannot
# interact with the ~/.nextflow used by your other pipelines.
export NXF_HOME="${CRETARGET_DATA_DIR}/.nextflow_home"

# Driver JVM. Modest, because the driver mostly waits.
export NXF_OPTS='-Xms500M -Xmx2G'

# --- Conda ------------------------------------------------------------------
# `conda activate` does not work in a non-interactive shell (a batch job)
# unless conda.sh has been sourced first. Doing it here means the same
# incantation works interactively and in a script.
#
# Adjust CRETARGET_CONDA_SH if you switch conda installs. On boltzmann there
# are two: /opt/miniconda3 (on PATH, read-only base) and ~/miniforge3 (what
# CONDA_EXE points at). They are interchangeable for activating an env, but
# mixing them in one shell causes confusing failures, so pick one explicitly.
export CRETARGET_CONDA_SH="${CRETARGET_CONDA_SH:-/opt/miniconda3/etc/profile.d/conda.sh}"
if [ -r "$CRETARGET_CONDA_SH" ]; then
    # shellcheck disable=SC1090
    . "$CRETARGET_CONDA_SH"
fi

# --- Python environment -----------------------------------------------------
# Which interpreter creTarget is installed into. On boltzmann this is an
# existing venv that already has pandas, numpy, pyfaidx and pytest.
export CRETARGET_PYENV="${CRETARGET_PYENV:-/home/jmorrislab1/mkazemi/myenv}"
if [ -r "${CRETARGET_PYENV}/bin/activate" ]; then
    # shellcheck disable=SC1091
    . "${CRETARGET_PYENV}/bin/activate"
fi

# --- Reference data ---------------------------------------------------------
# Pre-existing assets on this machine; no need to re-download ~3.2 GB.
# Keep the FASTA and the FlashFry index from the SAME assembly build.
export CRETARGET_FASTA="/home/jmorrislab1/mkazemi/.guidenet_reference/hg38.fa"
export CRETARGET_FLASHFRY_JAR="/home/jmorrislab1/mkazemi/bin/FlashFry-assembly-1.15.jar"
export CRETARGET_FLASHFRY_DB="/home/jmorrislab1/mkazemi/.guidenet_reference/hg38_standard_cas9ngg_db"

# --- Optional steps ---------------------------------------------------------
export CRETARGET_LAUNCH_SCRIPT="/home/jmorrislab1/mkazemi/___ProJECTS___/BMD_STING_seq/launch-dCas9/launch-dCas9.py"
export CRETARGET_LAUNCH_PYTHON="/home/jmorrislab1/mkazemi/miniforge3/envs/launch-dCas9/bin/python"
export CRETARGET_CKPT_DIR="/home/jmorrislab1/mkazemi/___ProJECTS___/BMD_STING_seq/checkpoints"
export CRETARGET_SEED_FINDER_DIR="/home/jmorrislab1/mkazemi/___ProJECTS___/Blair/gRNA_design_tool/seed_finder/seed-finder-lite-main"

mkdir -p "${CRETARGET_DATA_DIR}"/{cache,work,results} \
         "${CRETARGET_CODE_DIR}"/{inputs,logs} 2>/dev/null

echo "creTarget environment loaded."
echo "  code    : ${CRETARGET_CODE_DIR}"
echo "  data    : ${CRETARGET_DATA_DIR}"
echo "  cache   : ${CRETARGET_CACHE_DIR}"
echo "  work    : ${NXF_WORK}"
echo "  python  : $(command -v python3 2>/dev/null || echo '<none>')"
echo "  nextflow: $(command -v nextflow 2>/dev/null || echo '<none>')"
echo
echo "Open a new shell to discard all of this."
