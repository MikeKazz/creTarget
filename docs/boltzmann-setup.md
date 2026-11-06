# Running creTarget on boltzmann

Exact, copy-pasteable setup for the Morris lab machine. Every command states
the directory it runs from.

## Directory layout

Two filesystems, deliberately:

| | Path | Holds | Why |
|---|---|---|---|
| **CODE_DIR** | `/home/jmorrislab1/mkazemi/___ProJECTS___/creTarget` | git clone, input BEDs, driver script, SLURM logs | Small, backed up, the thing you edit and commit |
| **DATA_DIR** | `/home/jmorrislab0/mkazemi/creTarget` | reference genome, FlashFry index, model checkpoints, Nextflow `work/`, results | TB-scale storage |

Sizes that justify the split:

- reference genome + FlashFry index: **~3.2 GB**
- Nextflow `work/` on a genome-wide run: **~250–300 GB** (per-chunk discover,
  scored, annotated and ranked tables all live here until you delete them)
- collated output TSV for 2.35M cCREs: **~63 GB** (~132M guide rows at ~477
  bytes each, measured)

Shell variables used throughout:

These are set by `assets/cretarget_env.sh` as `CRETARGET_CODE_DIR` and
`CRETARGET_DATA_DIR`, which you source per session (step 2). They are
deliberately **not** added to `~/.bashrc`: this machine is shared and you have
other projects, so a global `CODE_DIR` would leak into all of them.

## What's unusual about this site

`boltzmann` runs SLURM, but with **exactly one node — itself**. It is both the
submit host and the only compute node:

```
PARTITION AVAIL  TIMELIMIT   NODES(A/I/O/T) NODELIST
gpu*         up   infinite          0/1/0/1 boltzmann
```

Three consequences shape the whole setup:

1. **The Nextflow driver competes with its own work.** Every CPU and GB the
   driver job holds is unavailable to the tasks it is waiting on. Submit the
   driver small — 1 CPU, 4 GB — not with the 8 CPU / 80 GB template used for
   ordinary analysis jobs. A fat driver starves the work and, in the worst
   case, deadlocks: it waits forever for jobs that can never be scheduled.
2. **Memory is the binding constraint, and it's shared.** 125 GB total with
   ~45 GB routinely in use by other people, so plan around ~75 GB. The profile
   sets explicit `maxForks` per process rather than trusting SLURM to gate on
   memory — that only works if memory is configured as a consumable resource.
   Confirm once:
   ```bash
   scontrol show config | grep -i selecttype
   ```
   If `SelectTypeParameters` does not include `CR_Memory`, the `maxForks`
   values in `conf/boltzmann.config` are the *only* thing preventing an
   overcommit.
3. **40 CPUs is really 20 physical cores**, hyperthreaded across 2 NUMA nodes.
   FlashFry is not meaningfully multithreaded here, so wide per-task CPU
   requests buy nothing and only reduce how many chunks run concurrently.

## GPUs

4 GPUs (TITAN Xp, 2× RTX 2080 Ti, TITAN X Pascal), shared with other users.

The `LAUNCH_DCAS9` step requests `--gres=gpu:1`, and **SLURM sets
`CUDA_VISIBLE_DEVICES` itself** — the pipeline never picks a device. That is
the correct behaviour on a shared box: choosing a device manually would risk
landing on a GPU allocated to someone else. `maxForks = 4` caps it at one task
per GPU.

Your torch environment reports `2.5.1+cu124` with `cuda.is_available() == True`,
which covers every card here (Pascal `sm_61` and Turing `sm_75` are both in
2.5.1's build matrix), so the driver-580 / CUDA-13 combination is not a problem.

If you ever run without SLURM (`-profile local`), nothing allocates GPUs and
every concurrent task would pile onto GPU 0. Pass `--gpu_count 4` in that case
and the pipeline round-robins devices by task index.

---

## Step 1 — Unpack and push the code

The public repo still holds v0.1.0, so the Nextflow pipeline has to be pushed
before anything can run. See `docs/release-checklist.md` item 1 for the exact
commands. After that, `$CRETARGET_CODE_DIR/creTarget` is your working clone.

## Step 2 — Load the project environment

**Nothing goes in `~/.bashrc`.** The settings below belong to this project
only; on a shared machine, exporting them globally is how one project quietly
breaks another. Source the file when you work on creTarget, and open a new
shell to discard it.

```bash
cd /home/jmorrislab1/mkazemi/___ProJECTS___/creTarget
cp creTarget/assets/cretarget_env.sh .
source cretarget_env.sh
```

It sets only `CRETARGET_*` and `NXF_*` names — no generic `CODE_DIR`/`DATA_DIR`
that could collide with another project — points `CRETARGET_CONFIG_DIR` inside
the project so `cretarget configure` does not write to your home, sources
`conda.sh`, and activates the Python environment. The sbatch driver sets all of
this itself in its own process, so a submitted job never depends on your
interactive shell.

## Step 3 — Install creTarget into the existing venv

Your `~/myenv` venv (Python 3.10.12) already has **pandas 2.2.3, numpy 1.26.4,
pyfaidx 0.9.0.3, pytest 9.1.0 and torch 2.5.1+cu121**, and its site-packages is
writable. That is everything creTarget's Python side needs, so a new conda
environment would be wasted work.

```bash
source cretarget_env.sh          # activates ~/myenv
pip install -e "$CRETARGET_CODE_DIR/creTarget"
cretarget --version              # -> creTarget 0.2.0
```

`-e` means `git pull` updates the installed command with no reinstall.

Already on the machine and found by the survey, so nothing to install:

| Tool | Version | Location |
|---|---|---|
| nextflow | 25.10.3 | `~/bin/nextflow` |
| java | openjdk 18.0.2-ea | `/usr/bin/java` |
| bedtools | v2.30.0 | `~/bin/bedtools` |
| samtools | 1.13 | `/usr/bin/samtools` |
| singularity | 4.2.0 | on PATH |

Nextflow 25.10.3 is comfortably above the pipeline's 24.04 floor, and the
pipeline is tested against it in CI.

> **If you ever do want a clean conda environment** — on another machine, or to
> isolate this from `myenv` — use `-profile conda` and Nextflow will build it
> from `assets/environment.yml`. Two notes specific to your setup: there are
> two conda installs here (`/opt/miniconda3` on PATH with a read-only base, and
> `~/miniforge3` which `CONDA_EXE` points at), so pick one explicitly rather
> than mixing them; and your `~/.condarc` lists **bioconda before
> conda-forge**, which is the reverse of what bioconda requires and a common
> cause of broken htslib/samtools solves. The pipeline's `environment.yml`
> declares the correct order itself, so it is unaffected, but it is worth
> fixing for your other projects.

## Step 4 — Reuse the reference data you already have

`cretarget setup-reference` downloads ~3.2 GB. **You do not need it** — the
survey found a complete set already on disk:

```
/home/jmorrislab1/mkazemi/.guidenet_reference/hg38.fa                        (3.1 GB)
/home/jmorrislab1/mkazemi/.guidenet_reference/hg38_standard_cas9ngg_db       (2.4 GB)
/home/jmorrislab1/mkazemi/.guidenet_reference/hg38_standard_cas9ngg_db.header
/home/jmorrislab1/mkazemi/bin/FlashFry-assembly-1.15.jar
```

`cretarget_env.sh` already points at these. One thing to verify — the survey
did not check for it, and the pipeline requires it:

```bash
ls -la "${CRETARGET_FASTA}.fai" || samtools faidx "$CRETARGET_FASTA"
```

The FASTA is used only by bedtools to extract region sequence; the FlashFry
database is used for the off-target search. They need not be the same file but
must be the same assembly, and this pairing (`hg38.fa` with the
`hg38_standard` index) is what the pre-Nextflow pipeline already used.

## Step 5 — Reuse the optional-step assets too

Also already present, so `cretarget fetch-models` is unnecessary here:

```
launch-dCas9 script : ~/___ProJECTS___/BMD_STING_seq/launch-dCas9/launch-dCas9.py
launch-dCas9 python : ~/miniforge3/envs/launch-dCas9/bin/python
checkpoints         : ~/___ProJECTS___/BMD_STING_seq/checkpoints/ckpt-seq-*.pth
seed-finder-lite    : ~/___ProJECTS___/Blair/gRNA_design_tool/seed_finder/seed-finder-lite-main
```

`cretarget_env.sh` exports all four. (The Zenodo fetch path still matters for
anyone else installing the tool; it is just redundant on this machine.)

## Step 6 — Record the paths for the step-level CLI

Only needed for `cretarget <step>` and `cretarget doctor`; Nextflow always gets
paths as explicit parameters. These write to
`$CRETARGET_CODE_DIR/.cretarget_config`, not your home directory.

```bash
cretarget configure --set reference.fasta="$CRETARGET_FASTA"
cretarget configure --set reference.flashfry_jar="$CRETARGET_FLASHFRY_JAR"
cretarget configure --set reference.flashfry_db="$CRETARGET_FLASHFRY_DB"
cretarget configure --set launch_dcas9.script="$CRETARGET_LAUNCH_SCRIPT"
cretarget configure --set launch_dcas9.python_bin="$CRETARGET_LAUNCH_PYTHON"
cretarget configure --set launch_dcas9.ckpt_dir="$CRETARGET_CKPT_DIR"
cretarget configure --set seed_finder.dir="$CRETARGET_SEED_FINDER_DIR"
cretarget configure --set cache.dir="$CRETARGET_CACHE_DIR"
```

## Step 7 — Verify

```bash
cretarget doctor --shared-path "$NXF_WORK"
```

**Every required check must say OK before you submit anything.** With the
paths above set, the optional launch-dCas9 and seed-finder checks should pass
too.

## Step 8 — Smoke test on a tiny input

Proves the full chain works on real reference data before you spend hours.

```bash
cd "$CRETARGET_CODE_DIR"
head -20 inputs/your_enhancers.txt > inputs/smoke20.txt   # or hand-write 20 lines

nextflow run "$CRETARGET_CODE_DIR/creTarget" \
    -profile boltzmann \
    --enhancer_bed "$CRETARGET_CODE_DIR/inputs/smoke20.txt" \
    --fasta        "$CRETARGET_FASTA" \
    --flashfry_jar "$CRETARGET_FLASHFRY_JAR" \
    --flashfry_db  "$CRETARGET_FLASHFRY_DB" \
    --chunk_size   10 \
    --outdir       "$CRETARGET_DATA_DIR/results/smoke20"
```

Run it in `tmux` so a dropped SSH session does not kill it. Check:

```bash
wc -l "$CRETARGET_DATA_DIR/results/smoke20/guides.ranked.tsv"     # expect ~1000+ rows
cut -f1 "$CRETARGET_DATA_DIR/results/smoke20/guides.ranked.tsv" | head -3
cat "$CRETARGET_DATA_DIR/results/smoke20/run_manifest.json"
```

If this fails, stop here — a genome-wide run will fail the same way, only
slower.

## Step 9 — Benchmark run on a subset

```bash
cd "$CRETARGET_CODE_DIR"
shuf -n 10000 inputs/your_enhancers.txt > inputs/enhancers_subset.txt
cp "$CRETARGET_CODE_DIR/creTarget/assets/run_cretarget.sbatch" "$CRETARGET_CODE_DIR/run_cretarget.sbatch"
# the template already points at inputs/enhancers_subset.txt
sbatch "$CRETARGET_CODE_DIR/run_cretarget.sbatch"
```

Monitor:

```bash
squeue -u $USER                              # driver + its task jobs
tail -f "$CRETARGET_CODE_DIR"/logs/cretarget_*.out
```

## Step 10 — Read the trace and tune

```bash
cd "$CRETARGET_DATA_DIR/results/pipeline_info"

# peak memory and runtime, per task
awk -F'\t' 'NR==1 || /FLASHFRY/ {print $4"\t"$9"\t"$11}' trace.txt | column -t
```

Then edit `$CRETARGET_CODE_DIR/creTarget/conf/boltzmann.config`:

- **`FLASHFRY_DISCOVER_SCORE` memory** is 24 GB, giving `discover` a 12 GB heap
  and `score` a 14 GB heap — the same range the pre-Nextflow pipeline used
  successfully. If peak RSS is well below that, lower it and raise `maxForks`.
- **`maxForks = 2`** on FlashFry is the main throughput lever. `2 × 24 GB = 48 GB`
  leaves comfortable headroom beside other users. On a quiet box, 3–4 is fine.
- **`chunk_size`** (default 5000) trades task count against per-task runtime.
  Fewer, larger chunks mean less scheduling overhead; the cost is a slower,
  more expensive retry when one fails.

`report.html` in the same directory shows this as a timeline, which is usually
the faster way to spot the bottleneck.

## Step 11 — Optional steps

Test these separately, on the smoke input, before combining them with a large
run. Neither has ever been executed end to end.

**launch-dCas9:**

```bash
# No fetch needed -- the checkpoints are already on this machine (step 5).
cd "$CRETARGET_CODE_DIR"
nextflow run "$CRETARGET_CODE_DIR/creTarget" \
    -profile boltzmann \
    --enhancer_bed "$CRETARGET_CODE_DIR/inputs/smoke20.txt" \
    --fasta        "$CRETARGET_FASTA" \
    --flashfry_jar "$CRETARGET_FLASHFRY_JAR" \
    --flashfry_db  "$CRETARGET_FLASHFRY_DB" \
    --chunk_size   10 \
    --run_launch_dcas9 \
    --launch_script "$CRETARGET_LAUNCH_SCRIPT" \
    --launch_python "$CRETARGET_LAUNCH_PYTHON" \
    --ckpt_dir     "$CRETARGET_CKPT_DIR" \
    --outdir       "$CRETARGET_DATA_DIR/results/smoke20_launch"
```

Confirm the new columns appeared and that a GPU was actually used:

```bash
head -1 "$CRETARGET_DATA_DIR/results/smoke20_launch/guides.ranked.tsv" | tr '\t' '\n' | grep launch
grep CUDA_VISIBLE "$NXF_WORK"/*/*/.command.err
```

**seed-finder-lite:** same shape, with `--run_seed_finder --seed_finder_dir "$CRETARGET_SEED_FINDER_DIR"
--seed_finder_min_length 8 --seed_finder_tss_window 1000`.

## Step 12 — Housekeeping after a big run

The `work/` directory is not cleaned automatically.

```bash
du -sh "$NXF_WORK"
cd "$CRETARGET_CODE_DIR" && nextflow clean -f -before <run_name>   # keeps the latest run resumable
```

Do **not** delete `work/` while you still want `-resume` to work.

## Genome-wide run: things to decide first

A 2.35M-cCRE run produces a single ~63 GB TSV. Before committing to it:

- **Does a 63 GB single file suit your downstream analysis?** Options are
  gzipping the collated output (~10–15 GB, still `zgrep`-able), splitting per
  chromosome, or writing Parquet. This is a pipeline change, not a flag —
  decide before the run, not after.
- **Budget ~300 GB for `work/`** on top of the output.
- **Run the 10k subset first** and confirm scaling is linear before
  extrapolating.
