# Setting creTarget up on a new cluster

This walks through a first install end to end: survey the machine, install the
software, fetch the reference data, prove it works on a tiny input, then scale
up. Every step says where to run it from.

Throughout, two directories are referred to:

- `$CODE_DIR` — the pipeline, your inputs, and logs. Small, and worth backing up.
- `$DATA_DIR` — reference data, the Nextflow work directory, and results.
  Large: budget roughly 300 GB of work directory and 50 GB of output for a
  genome-wide run.

On many clusters these are the same place. On clusters with a small backed-up
home and a large scratch filesystem, they are not. Set them once:

```bash
export CODE_DIR=/path/to/your/project
export DATA_DIR=/path/to/large/storage
mkdir -p "$CODE_DIR" "$DATA_DIR"
```

---

## 1. Survey the cluster

From `$CODE_DIR`, after cloning the repository:

```bash
cd "$CODE_DIR"
git clone https://github.com/MikeKazz/creTarget.git
cd creTarget
bash assets/inspect_cluster.sh . "$CODE_DIR" "$DATA_DIR"
```

This is read-only. It writes a report to the current directory and prints a
summary. Four of its answers decide whether you need a site profile at all:

| Question | Where the answer appears | Why it matters |
|---|---|---|
| Does the scheduler account for memory when placing jobs? | `SelectTypeParameters` | If not, it will overcommit the node whatever you request, and you must cap concurrency yourself. |
| Does it *enforce* the memory it granted? | `TaskPlugin` | With `task/none`, an over-sized task is not killed cleanly — the kernel OOM killer picks a victim, possibly someone else's job. |
| Is the node shared with work the scheduler cannot see? | `AllocMem` vs `FreeMem` | An interactive box where people also run things outside the scheduler needs `maxForks`. |
| Does the work directory support file locking? | the `flock` test | Nextflow's cache requires it. NFS without locking will fail in confusing ways. |

If the scheduler accounts for and enforces memory, and the nodes are not
shared outside it, the generic `-profile slurm` is enough and you can skip
step 5.

---

## 2. Install the software

Nextflow needs Java 17 or newer. From anywhere:

```bash
java -version                      # must be >= 17
curl -s https://get.nextflow.io | bash
mkdir -p "$CODE_DIR/bin" && mv nextflow "$CODE_DIR/bin/"
export PATH="$CODE_DIR/bin:$PATH"
nextflow -version                  # must be >= 24.04.0
```

> 24.04 is a hard floor: `process.resourceLimits`, which clamps every resource
> request, was introduced there. On older versions the ceilings silently do
> nothing.

For the pipeline's own dependencies, pick one:

**Container (most portable).** Nothing to install but the engine. One image,
no dependency solve, no version skew between users. Add
`-profile singularity` (or `docker`, `apptainer`) to every run.

**Conda.** `-profile conda` builds the pinned environment in
`assets/environment.yml` on first run. If conda's solver is slow, use
`-profile micromamba` — same environment, seconds instead of an hour.

**Manual.** `pip install cretarget`, plus `bedtools`, `samtools` and a JVM on
`PATH`.

Do not add any of this to `~/.bashrc`. Put the exports in a small file you
source per session — `assets/site_env.sh.example` is a template to copy and
fill in — so that other work on a shared machine is unaffected.

---

## 3. Fetch the reference data

From `$CODE_DIR/creTarget`:

```bash
export CRETARGET_CACHE_DIR="$DATA_DIR/cache"
cretarget setup-reference --genome hg38
```

That downloads the genome FASTA and the optional model checkpoints (~3.2 GB
total) into `$CRETARGET_CACHE_DIR`. The environment variable takes precedence
over anything in `config.toml`, so set it before running if your home
directory is small.

FlashFry is a separate JVM tool with its own index, not fetched by
`setup-reference`:

```bash
mkdir -p "$DATA_DIR/tools" && cd "$DATA_DIR/tools"
curl -LO https://github.com/mckennalab/FlashFry/releases/download/1.15/FlashFry-assembly-1.15.jar

cd "$DATA_DIR/cache/reference"
samtools faidx hg38.fa                 # if not already indexed
java -Xmx32g -jar "$DATA_DIR/tools/FlashFry-assembly-1.15.jar" index \
    --tmpLocation ./tmp \
    --database hg38_standard_cas9ngg_db \
    --reference hg38.fa \
    --enzyme spcas9ngg
```

Indexing takes roughly an hour and wants 32 GB. Do it once; it is reused by
every run.

If you plan to use `--run_seed_finder`, build its promoter reference once.
This needs `bedtools` and the FASTA from above, and downloads a pinned GENCODE
annotation:

```bash
cd "$CODE_DIR/creTarget"
cretarget build-seed-reference --fasta "$DATA_DIR/cache/reference/hg38.fa"
```

It writes a TSS table, a promoter-window FASTA and a `provenance.json`
recording the annotation release into `$CRETARGET_CACHE_DIR/seed_reference`.
Takes a few minutes, and is reused by every run.

Verify:

```bash
cd "$CODE_DIR/creTarget"
cretarget doctor
```

---

## 4. Prove it works on a tiny input

From `$CODE_DIR/creTarget`:

```bash
nextflow run . -profile test,conda --outdir "$DATA_DIR/results_test"
```

The `test` profile uses a handful of regions and finishes in minutes. It is
the fastest way to find a broken Java version, an unreadable index, or a work
directory that does not support locking — before a long run finds them for
you.

Then a real but small input of your own:

```bash
head -200 /path/to/your/targets.bed > "$CODE_DIR/inputs/smoke.bed"

nextflow run . -profile slurm,conda \
    --enhancer_bed "$CODE_DIR/inputs/smoke.bed" \
    --fasta        "$DATA_DIR/cache/reference/hg38.fa" \
    --flashfry_jar "$DATA_DIR/tools/FlashFry-assembly-1.15.jar" \
    --flashfry_db  "$DATA_DIR/cache/reference/hg38_standard_cas9ngg_db" \
    --outdir       "$DATA_DIR/results_smoke"
```

Check `$DATA_DIR/results_smoke/pipeline_info/trace.txt` for the peak memory
and runtime each process actually used. Those numbers are what you put in a
site profile, if you need one.

---

## 5. Write a site profile, if the survey said you need one

Copy the template and fill it in from the trace and the survey:

```bash
cd "$CODE_DIR/creTarget"
cp conf/site-template.config conf/mysite.config
```

`conf/example-singlenode.config` is a filled-in example for the awkward case —
a single-node cluster with shared GPUs, where the scheduler accounts for
memory but does not enforce it — with the reasoning behind each number
documented inline.

Register it in the `profiles` block of `nextflow.config`:

```groovy
mysite {
    includeConfig 'conf/mysite.config'
}
```

Or keep it outside the repository entirely and pass `-c /path/to/mysite.config`,
which avoids carrying a local edit in git.

Two things to size from your trace rather than guess:

- **`PREP_REGIONS` memory** scales with the *total* number of targets, not with
  chunk size: every region is held in memory at once to flag linked
  neighbours. It also shells out to `bedtools getfasta`, whose memory does not
  show up in the Python process's own accounting — size it from the trace's
  `peak_rss`, not from a Python measurement.
- **`chunk_size`** trades scheduling overhead against the cost of a retry.
  FlashFry's runtime scales superlinearly with chunk size (measured at roughly
  n^1.3), so very large chunks are not the win they look like.

---

## 6. Scale up

Submit the driver as a batch job. `assets/run_cretarget.sbatch` is a template;
copy it next to your inputs and edit the `CONFIGURE` block at the top.

The one thing to get right: **request minimal resources for the driver job**
(1 CPU, a few GB). It only submits work and waits. On a single-node cluster,
every resource the driver holds is one the actual work cannot have — ask for
8 CPUs and 80 GB and you can deadlock, waiting forever for jobs that can never
be scheduled.

Add `-resume` to every invocation. A run that dies at hour 20 then picks up
from the last completed chunk rather than starting over.

---

## Troubleshooting

**`Process substitution ... not found`** — Nextflow stages inputs as files and
cannot stage a `<(...)` or a pipe. Write your targets to a real file first.

**Cached Singularity image not found / points at an unmounted path** — a stale
`NXF_SINGULARITY_CACHEDIR` in your environment. Override it to somewhere under
`$DATA_DIR`.

**`conda activate` fails inside a batch job** — a batch job gets a
non-interactive shell, so the conda block in `~/.bashrc` never runs. Source
`.../etc/profile.d/conda.sh` first; the sbatch template does this for you.

**Tasks rejected rather than queued** — a request escalated past what your
partition will grant. Lower `max_cpus` / `max_memory` / `max_time` in your
site profile to the largest single allocation the partition actually allows.

**Output ordering differs between runs** — it should not; this is a bug worth
reporting. The pipeline collates chunks by numeric chunk index precisely so
that two runs of the same input produce byte-identical output.
