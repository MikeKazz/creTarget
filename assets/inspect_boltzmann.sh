#!/usr/bin/env bash
#
# inspect_boltzmann.sh -- read-only survey of the cluster before installing
# anything.
#
# SAFETY: this script is strictly read-only outside its own report file. It
# does not install, modify ~/.bashrc, create conda environments, export
# anything into your shell, or write anywhere except $REPORT (which lives in
# the project directory you pass in, default: the directory this script is
# run from). Every external command is guarded, so a missing tool produces a
# "not found" line rather than aborting the run.
#
# Usage:
#   bash inspect_boltzmann.sh                 # report into $PWD
#   bash inspect_boltzmann.sh /path/to/outdir
#
# Then paste the report, or upload the file.

set -uo pipefail   # deliberately NOT -e: a failed probe must not stop the survey

OUTDIR="${1:-$PWD}"
REPORT="${OUTDIR}/boltzmann_inspection_$(date +%Y%m%d_%H%M%S).txt"

# Paths this project cares about. Adjust if they ever change.
CODE_DIR="/home/jmorrislab1/mkazemi/___ProJECTS___/creTarget"
DATA_DIR="/home/jmorrislab0/mkazemi/creTarget"

# $USER and $HOME are not guaranteed to be set (a batch job may inherit
# neither), and this script runs under `set -u`, so resolve them defensively
# rather than letting an unbound variable abort the survey mid-way.
CURRENT_USER="${USER:-${LOGNAME:-$(id -un 2>/dev/null || echo unknown)}}"
CURRENT_HOME="${HOME:-$(getent passwd "$CURRENT_USER" 2>/dev/null | cut -d: -f6)}"
CURRENT_HOME="${CURRENT_HOME:-/tmp}"

mkdir -p "$OUTDIR" 2>/dev/null

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
section() { printf '\n==================== %s ====================\n' "$1"; }
sub()     { printf '\n--- %s ---\n' "$1"; }

# Run a command only if it exists; never let a failure kill the script.
try() {
    local cmd="$1"; shift
    if command -v "$cmd" >/dev/null 2>&1; then
        "$cmd" "$@" 2>&1 | head -40
    else
        echo "[not found: $cmd]"
    fi
}

have() { command -v "$1" >/dev/null 2>&1; }

# Walk up to the nearest existing ancestor. DATA_DIR may not exist yet, but we
# still need to know the free space and filesystem type of the volume it will
# live on -- reporting "does not exist" there would tell us nothing useful.
nearest_existing() {
    local p="$1"
    while [ -n "$p" ] && [ "$p" != "/" ] && [ ! -e "$p" ]; do
        p="$(dirname "$p")"
    done
    printf '%s' "$p"
}

# ---------------------------------------------------------------------------
# the survey
# ---------------------------------------------------------------------------
survey() {

printf 'creTarget cluster inspection\n'
printf 'generated: %s\n' "$(date -Is 2>/dev/null || date)"
printf 'script run from: %s\n' "$PWD"
printf 'report file: %s\n' "$REPORT"

# =========================================================================
section "1. HOST AND OS"
# =========================================================================
echo "hostname : $(hostname -f 2>/dev/null || hostname)"
echo "kernel   : $(uname -srm)"
echo "uptime   : $(uptime 2>/dev/null)"
sub "distribution"
if [ -r /etc/os-release ]; then grep -E '^(NAME|VERSION)=' /etc/os-release; else echo "[no /etc/os-release]"; fi
sub "locale (a non-UTF-8 locale can break tool output parsing)"
echo "LANG=${LANG:-<unset>}  LC_ALL=${LC_ALL:-<unset>}"
try locale
sub "current shell / login shell"
echo "SHELL=${SHELL:-<unset>}  BASH_VERSION=${BASH_VERSION:-<unset>}"
sub "ulimits (nofile matters: Nextflow opens many files with many concurrent tasks)"
echo "open files (soft/hard): $(ulimit -Sn)/$(ulimit -Hn)"
echo "max user processes    : $(ulimit -Su)/$(ulimit -Hu)"
echo "max locked memory     : $(ulimit -Sl)"
sub "TMPDIR (FlashFry's JVM and Nextflow both write scratch here)"
echo "TMPDIR=${TMPDIR:-<unset, defaults to /tmp>}"
df -h "${TMPDIR:-/tmp}" 2>&1 | tail -2

# =========================================================================
section "2. ENVIRONMENT MODULES"
# =========================================================================
if have module || [ -n "${MODULESHOME:-}" ] || [ -n "${LMOD_CMD:-}" ]; then
    echo "module system present (MODULESHOME=${MODULESHOME:-<unset>}, LMOD_CMD=${LMOD_CMD:-<unset>})"
    sub "currently loaded"
    module list 2>&1 | head -20
    sub "modules matching things we need (java/nextflow/bedtools/samtools/conda/singularity)"
    for m in java jdk openjdk nextflow bedtools samtools conda anaconda miniconda \
             singularity apptainer python cuda; do
        out=$(module avail "$m" 2>&1 | grep -viE '^-|^$|No module' | head -3)
        [ -n "$out" ] && printf '  %-12s -> %s\n' "$m" "$(echo "$out" | tr '\n' ' ')"
    done
else
    echo "[no environment-modules / Lmod detected]"
fi

# =========================================================================
section "3. CONDA  (the section I have the most questions about)"
# =========================================================================
sub "which conda / mamba"
echo "conda      : $(command -v conda       || echo '[not found]')"
echo "mamba      : $(command -v mamba       || echo '[not found]')"
echo "micromamba : $(command -v micromamba  || echo '[not found]')"
echo "CONDA_EXE  : ${CONDA_EXE:-<unset>}"
echo "CONDA_PREFIX (active env): ${CONDA_PREFIX:-<none active>}"
echo "CONDA_DEFAULT_ENV        : ${CONDA_DEFAULT_ENV:-<none>}"

if have conda; then
    sub "versions"
    conda --version 2>&1
    have mamba && mamba --version 2>&1

    sub "conda info"
    conda info 2>&1 | head -30

    sub "channels and channel_priority (strict priority matters for bioconda)"
    conda config --show channels channel_priority 2>&1

    sub "where environments and package caches live (and whether that is on a quota'd filesystem)"
    conda config --show envs_dirs pkgs_dirs 2>&1

    sub "existing environments"
    conda env list 2>&1

    sub "size of conda base and package cache"
    CONDA_BASE="$(conda info --base 2>/dev/null)"
    [ -n "$CONDA_BASE" ] && du -sh "$CONDA_BASE" 2>/dev/null
    for d in $(conda config --show pkgs_dirs 2>/dev/null | grep -oE '/[^ ]+'); do
        du -sh "$d" 2>/dev/null
    done

    sub "CRITICAL: does 'conda activate' work in a NON-interactive shell?"
    echo "(this is what a SLURM batch script gets; if it fails, the sbatch"
    echo " script must source conda.sh explicitly instead of calling activate)"
    # No backticks in the echoed text: inside bash -c "..." they would be
    # treated as command substitution rather than printed.
    if bash -c 'conda activate base' >/dev/null 2>&1; then
        echo "  bare 'conda activate' in non-interactive bash: WORKS"
    else
        echo "  bare 'conda activate' in non-interactive bash: FAILS (expected on most setups)"
    fi
    if [ -n "$CONDA_BASE" ] && [ -r "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
        if bash -c "source '$CONDA_BASE/etc/profile.d/conda.sh' && conda activate base" >/dev/null 2>&1; then
            echo "  source conda.sh then activate: WORKS  <-- use this form in sbatch"
        else
            echo "  source conda.sh then activate: FAILS"
        fi
        echo "  conda.sh path: $CONDA_BASE/etc/profile.d/conda.sh"
    else
        echo "  [could not find \$CONDA_BASE/etc/profile.d/conda.sh]"
    fi

    sub "is conda initialised in a shell rc file? (read-only check, nothing modified)"
    for rc in "$CURRENT_HOME/.bashrc" "$CURRENT_HOME/.bash_profile" "$CURRENT_HOME/.profile"; do
        if [ -r "$rc" ]; then
            if grep -q 'conda initialize' "$rc" 2>/dev/null; then
                echo "  $rc : contains 'conda initialize' block"
            else
                echo "  $rc : no conda initialize block"
            fi
        else
            echo "  $rc : [absent]"
        fi
    done

    sub "contents of the currently active env (looking for tools we would otherwise reinstall)"
    if [ -n "${CONDA_DEFAULT_ENV:-}" ]; then
        conda list -n "${CONDA_DEFAULT_ENV}" 2>/dev/null \
            | grep -iE '^(python|pandas|numpy|bedtools|samtools|openjdk|jdk|nextflow|torch|pytorch|pyfaidx|htslib) ' \
            || echo "  (none of the relevant packages found in ${CONDA_DEFAULT_ENV})"
    else
        echo "  [no active env]"
    fi

    sub "can conda solve against conda-forge/bioconda? (dry run, installs nothing)"
    echo "  running: conda create --dry-run -n _cretarget_probe -c conda-forge -c bioconda bedtools samtools"
    timeout 300 conda create --dry-run -n _cretarget_probe \
        -c conda-forge -c bioconda bedtools samtools 2>&1 | tail -15
    echo "  (--dry-run: no environment was created)"
else
    echo "[conda not on PATH -- check the module system section above]"
fi

# =========================================================================
section "4. PYTHON AND PIP"
# =========================================================================
echo "python3 : $(command -v python3 || echo '[not found]')"
try python3 --version
sub "pip"
echo "pip : $(command -v pip || command -v pip3 || echo '[not found]')"
try python3 -m pip --version
sub "site-packages and whether this interpreter is writable (editable install target)"
python3 - <<'PYEOF' 2>&1
import os, site, sys, sysconfig
print("prefix        :", sys.prefix)
print("base_prefix   :", sys.base_prefix)
print("in_venv/conda :", sys.prefix != sys.base_prefix)
sp = sysconfig.get_paths().get("purelib")
print("purelib       :", sp)
print("writable      :", os.access(sp, os.W_OK) if sp else "unknown")
PYEOF
sub "PYTHONPATH (a stray value here can shadow an installed package)"
echo "PYTHONPATH=${PYTHONPATH:-<unset>}"
sub "relevant packages already importable"
python3 - <<'PYEOF' 2>&1
for m in ("pandas", "numpy", "torch", "pyfaidx", "pytest"):
    try:
        mod = __import__(m)
        print(f"  {m:10s} {getattr(mod, '__version__', 'unknown')}")
    except Exception as e:
        print(f"  {m:10s} [not importable: {type(e).__name__}]")
PYEOF

# =========================================================================
section "5. JAVA  (FlashFry and Nextflow both need a JRE >= 11)"
# =========================================================================
echo "java : $(command -v java || echo '[not found]')"
echo "JAVA_HOME=${JAVA_HOME:-<unset>}"
if have java; then java -version 2>&1 | head -3; fi

# =========================================================================
section "6. NEXTFLOW"
# =========================================================================
echo "nextflow : $(command -v nextflow || echo '[not found]')"
echo "NXF_HOME=${NXF_HOME:-<unset, defaults to ~/.nextflow>}"
echo "NXF_WORK=${NXF_WORK:-<unset>}"
echo "NXF_OPTS=${NXF_OPTS:-<unset>}"
if have nextflow; then
    sub "version"
    nextflow -version 2>&1 | head -8
fi
sub "existing ~/.nextflow (size; it caches assets and plugins)"
du -sh "$CURRENT_HOME/.nextflow" 2>/dev/null || echo "[absent]"

# =========================================================================
section "7. BIOINFORMATICS TOOLS"
# =========================================================================
for t in bedtools samtools bgzip tabix; do
    if have "$t"; then
        case "$t" in
            bedtools) v=$(bedtools --version 2>&1 | head -1) ;;
            samtools) v=$(samtools --version 2>&1 | head -1) ;;
            *)        v=$("$t" --version 2>&1 | head -1) ;;
        esac
        printf '  %-10s %s  (%s)\n' "$t" "$v" "$(command -v $t)"
    else
        printf '  %-10s [not found]\n' "$t"
    fi
done

# =========================================================================
section "8. CONTAINERS"
# =========================================================================
for t in singularity apptainer docker podman; do
    if have "$t"; then
        printf '  %-12s %s\n' "$t" "$($t --version 2>&1 | head -1)"
    else
        printf '  %-12s [not found]\n' "$t"
    fi
done
echo "SINGULARITY_CACHEDIR=${SINGULARITY_CACHEDIR:-<unset>}"
echo "APPTAINER_CACHEDIR=${APPTAINER_CACHEDIR:-<unset>}"

# =========================================================================
section "9. SLURM"
# =========================================================================
sub "partitions (detail)"
try sinfo -o "%20P %5a %10l %6D %6t %8z %10m %15G %N"
sub "node detail"
try scontrol show node
sub "partition detail"
try scontrol show partition
sub "KEY CONFIG -- does SLURM treat memory as a consumable resource?"
echo "(if SelectTypeParameters lacks CR_Memory, SLURM will NOT gate on --mem,"
echo " and the maxForks values in conf/boltzmann.config are the only thing"
echo " preventing a memory overcommit)"
if have scontrol; then
    scontrol show config 2>/dev/null | grep -iE \
      'SelectType|SelectTypeParameters|DefMemPerNode|DefMemPerCPU|MaxMemPerNode|MaxMemPerCPU|MaxArraySize|MaxJobCount|MaxSubmitJobs|GresTypes|SchedulerType|AccountingStorageType|JobAcctGatherType|PriorityType|EnforcePartLimits|TaskPlugin'
else
    echo "[scontrol not found]"
fi
sub "accounting -- is sacct usable? (affects nothing in the Nextflow design, but good to know)"
if have sacct; then
    sacct --allocations --starttime now-7days --format=JobID,JobName%20,State,Elapsed,ReqMem,MaxRSS 2>&1 | head -12
else
    echo "[sacct not found]"
fi
sub "my account / QOS / limits"
try sacctmgr -n show assoc where user="$CURRENT_USER" format=Account,User,Partition,QOS,GrpTRES,MaxTRES,MaxJobs,MaxSubmit
try sacctmgr -n show qos format=Name,MaxTRESPerUser,MaxJobsPerUser,MaxWall
sub "current queue state"
try squeue -o "%.10i %.12u %.10P %.8T %.10M %.6C %.10m %R"
sub "what is running on the node right now (who else is using the memory)"
try squeue -h -o "%u %C %m %T"

# =========================================================================
section "10. GPUs"
# =========================================================================
if have nvidia-smi; then
    nvidia-smi 2>&1 | head -25
    sub "compute capability per device (launch-dCas9's torch build must support these)"
    nvidia-smi --query-gpu=index,name,compute_cap,memory.total,memory.used,utilization.gpu \
               --format=csv 2>&1
else
    echo "[nvidia-smi not found]"
fi
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-<unset>}"
sub "torch's view of the GPUs (if torch is importable in this env)"
python3 - <<'PYEOF' 2>&1
try:
    import torch
    print("torch           :", torch.__version__)
    print("cuda available  :", torch.cuda.is_available())
    print("cuda build ver  :", torch.version.cuda)
    print("device count    :", torch.cuda.device_count())
    for i in range(torch.cuda.device_count()):
        print(f"  [{i}] {torch.cuda.get_device_name(i)}  sm_{''.join(map(str, torch.cuda.get_device_capability(i)))}")
    print("arch list       :", torch.cuda.get_arch_list())
except Exception as e:
    print(f"[torch not usable here: {type(e).__name__}: {e}]")
PYEOF

# =========================================================================
section "11. STORAGE"
# =========================================================================
sub "free space (for a path that does not exist yet, its nearest existing parent)"
for d in "$CURRENT_HOME" "$CODE_DIR" "$DATA_DIR" /tmp "${TMPDIR:-/tmp}"; do
    probe="$(nearest_existing "$d")"
    if [ "$probe" = "$d" ]; then
        printf '%s\n' "$d"
    else
        printf '%s\n  [not created yet; reporting its parent %s]\n' "$d" "$probe"
    fi
    df -h "$probe" 2>&1 | tail -1
done

sub "filesystem type and device id (same device id => same filesystem)"
for d in "$CURRENT_HOME" "$CODE_DIR" "$DATA_DIR" /tmp; do
    probe="$(nearest_existing "$d")"
    [ -e "$probe" ] || continue
    printf '  %-58s type=%-10s dev=%s%s\n' "$d" \
        "$(stat -f -c '%T' "$probe" 2>/dev/null)" \
        "$(stat -c '%d' "$probe" 2>/dev/null)" \
        "$([ "$probe" = "$d" ] || printf ' (via %s)' "$probe")"
done
echo
echo "(Nextflow needs its work directory on a filesystem that supports file"
echo " locking. NFS usually does; some configurations do not. Tested below.)"

sub "mount options for the two project filesystems"
for d in "$CODE_DIR" "$DATA_DIR"; do
    probe="$(nearest_existing "$d")"
    [ -e "$probe" ] && findmnt -T "$probe" -o TARGET,SOURCE,FSTYPE,OPTIONS 2>&1 | tail -2
done

sub "quotas"
try quota -s
if have lfs; then
    sub "lustre quotas"
    for d in "$CODE_DIR" "$DATA_DIR"; do
        [ -e "$d" ] && lfs quota -h -u "$CURRENT_USER" "$d" 2>&1 | head -5
    done
fi

sub "write + file-locking test (creates and removes one small temp file per dir)"
for raw in "$CODE_DIR" "$DATA_DIR"; do
    d="$(nearest_existing "$raw")"
    if [ -d "$d" ] && [ -w "$d" ]; then
        probe="$d/.cretarget_probe.$$"
        if echo probe > "$probe" 2>/dev/null; then
            lockres="locking: unknown"
            if have flock; then
                if flock -n "$probe" true 2>/dev/null; then lockres="locking: OK"; else lockres="locking: FAILED"; fi
            fi
            echo "  $d  write: OK   $lockres"
            rm -f "$probe"
        else
            echo "  $d  write: FAILED"
        fi
    else
        echo "  $d  [missing or not writable]"
    fi
done

sub "sequential write throughput, 256 MB (writes then deletes one temp file)"
for raw in "$DATA_DIR" "$CODE_DIR"; do
    d="$(nearest_existing "$raw")"
    if [ -d "$d" ] && [ -w "$d" ]; then
        echo "  $d"
        dd if=/dev/zero of="$d/.cretarget_dd.$$" bs=1M count=256 conv=fdatasync 2>&1 | tail -1
        rm -f "$d/.cretarget_dd.$$"
    fi
done

# =========================================================================
section "12. NETWORK REACHABILITY FROM THIS HOST"
# =========================================================================
echo "(the reference genome, model checkpoints and pip packages all come from"
echo " these hosts; a proxy-only environment shows up here)"
echo "http_proxy=${http_proxy:-<unset>}  https_proxy=${https_proxy:-<unset>}"
for url in https://github.com https://zenodo.org https://pypi.org https://conda.anaconda.org; do
    if have curl; then
        code=$(curl -s -o /dev/null -w '%{http_code}' --max-time 15 -I "$url" 2>/dev/null)
        printf '  %-35s HTTP %s\n' "$url" "${code:-no-response}"
    fi
done

# =========================================================================
section "13. EXISTING PROJECT ASSETS (avoid re-downloading 3.2 GB)"
# =========================================================================
sub "any hg38 reference already on disk?"
for base in "$CURRENT_HOME" /home/jmorrislab0/mkazemi /home/jmorrislab1/mkazemi; do
    [ -d "$base" ] && find "$base" -maxdepth 4 \
        \( -name 'hg38*.fa' -o -name 'hg38*.fa.gz' -o -name 'GRCh38*.fa' \) \
        -printf '  %10s  %p\n' 2>/dev/null | head -10
done
echo "  (empty result = nothing found within 4 levels)"

sub "any FlashFry jar or index already on disk?"
for base in "$CURRENT_HOME" /home/jmorrislab0/mkazemi /home/jmorrislab1/mkazemi; do
    [ -d "$base" ] && find "$base" -maxdepth 5 \
        \( -iname 'FlashFry*.jar' -o -name '*cas9ngg_db*' \) \
        -printf '  %10s  %p\n' 2>/dev/null | head -10
done

sub "launch-dCas9 and seed-finder-lite"
for base in "$CURRENT_HOME" /home/jmorrislab0/mkazemi /home/jmorrislab1/mkazemi; do
    [ -d "$base" ] && find "$base" -maxdepth 5 \
        \( -iname 'launch-dCas9*' -o -iname 'seed-finder*' -o -iname 'targetfinder' \) \
        -printf '  %p\n' 2>/dev/null | head -20
done

sub "existing model checkpoints (.pth)"
for base in "$CURRENT_HOME" /home/jmorrislab0/mkazemi /home/jmorrislab1/mkazemi; do
    [ -d "$base" ] && find "$base" -maxdepth 6 -name 'ckpt-seq-*.pth' \
        -printf '  %10s  %p\n' 2>/dev/null | head -10
done

sub "project directory contents right now"
for d in "$CODE_DIR" "$DATA_DIR"; do
    echo "  $d:"
    if [ -d "$d" ]; then
        ls -la "$d" 2>/dev/null | head -15 | sed 's/^/    /'
    else
        echo "    [not created yet]"
    fi
done

# =========================================================================
section "14. POTENTIALLY INTERFERING ENVIRONMENT VARIABLES"
# =========================================================================
echo "(listed so we can design around them rather than overriding them)"
env | grep -iE '^(CONDA|NXF|JAVA|PYTHON|CUDA|SINGULARITY|APPTAINER|SLURM|TMPDIR|CRETARGET|LD_LIBRARY_PATH|PERL5LIB|R_LIBS)' \
    | sort | sed 's/^/  /' || echo "  (none set)"

printf '\n==================== END OF REPORT ====================\n'
}

# ---------------------------------------------------------------------------
# run it: to the report file and to the terminal at once
# ---------------------------------------------------------------------------
survey 2>&1 | tee "$REPORT"

printf '\nReport written to:\n  %s\n' "$REPORT"
printf 'Size: %s\n' "$(du -h "$REPORT" 2>/dev/null | cut -f1)"
printf '\nNothing outside that file was modified.\n'
