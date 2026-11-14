# creTarget

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24148024.svg)](https://doi.org/10.5281/zenodo.24148024)
[![Container](https://img.shields.io/badge/container-ghcr.io-blue)](https://github.com/MikeKazz/creTarget/pkgs/container/cretarget)

Design, score, and rank CRISPRi guide RNAs for noncoding target regions —
SNPs/points or enhancers/intervals — in hg38, at any scale from a handful
of targets to genome-wide.

Off-target and on-target scoring comes from
[FlashFry](https://github.com/mckennalab/FlashFry), with optional CNN-based
enhancer/promoter efficacy scoring via
[launch-dCas9](https://github.com/Wancen/launch-dCas9), optional ATAC-seq
chromatin-accessibility annotation, and optional off-target/TSS-proximity
screening via seed-finder-lite.

Execution is handled by [Nextflow](https://nextflow.io), so the same command
runs on your laptop, on SLURM, PBS/PBS Pro, SGE, or LSF — the scheduler is a
`-profile` flag, not a code path.

## Features

- **Two input types, flexible formats.** SNPs (`chr:pos`, 1-based, or BED)
  and enhancers (`chr:start-end`, or BED), auto-detected per line. Malformed
  lines are skipped with a warning rather than aborting the run. Contigs are
  validated against the reference up front, so an Ensembl-vs-UCSC naming
  mismatch is reported in plain language instead of surfacing as a bedtools
  error three steps later.
- **Flat memory at any scale.** Region prep is the only global step;
  everything after it — FlashFry, annotation, ranking, optional scoring — runs
  per chunk with no cross-chunk synchronization. A 2-million-region run uses
  the same memory per task as a 3,000-region run, only with more tasks.
- **Scheduler-agnostic.** `-profile slurm`, `pbs`, `sge`, or `lsf`. Failed
  tasks retry automatically with escalated memory and time; `-resume` picks up
  where an interrupted run left off.
- **Real ranking.** Hard filters (seed score, GC content, poly-T, multi-site
  genomic match, FlashFry overflow) exclude bad guides. Passing guides get a
  weighted composite score (specificity, CFD, closest off-target distance, RNA
  folding energy, and optionally launch-dCas9 efficacy) and a per-target rank.
- **Provenance by default.** Every output directory gets a `run_manifest.json`
  (pipeline version, commit, full parameters, command line) and a
  `software_versions.yml`.

## Requirements

- [Nextflow](https://nextflow.io) **≥ 24.04** — this is a hard floor, not a
  suggestion: the pipeline uses `process.resourceLimits`, which earlier
  versions ignore, so resource ceilings would silently not apply. Tested on
  24.04, 25.10 and latest-stable in CI.
- Java **17+** (Nextflow and FlashFry both need a modern JVM). Note that
  Nextflow 23.x caps at Java 19, which is a second reason not to use it.
- Python ≥ 3.10
- `bedtools`, `samtools`
- [FlashFry](https://github.com/mckennalab/FlashFry/releases) (a JVM tool;
  download the assembly JAR separately)

A GPU is not required. The core pipeline is CPU-only; a GPU only helps the
optional `--run_launch_dcas9` step.

## Installation

```bash
pip install cretarget
```

Or with conda, which also brings `bedtools`, `samtools`, and a JRE:

```bash
conda install -c conda-forge -c bioconda cretarget
```

Check that everything is wired up:

```bash
cretarget doctor
```

## First-time setup

```bash
# 1. Reference genome + FlashFry index (~3.2GB, Zenodo DOI 10.5281/zenodo.21115562)
cretarget setup-reference

# 2. Record where things live (used by the step-level CLI; paths may use ~ and $VARS)
cretarget configure --set reference.fasta=~/.cretarget/reference/hg38.fa
cretarget configure --set reference.flashfry_jar=/path/to/FlashFry-assembly-1.15.jar
cretarget configure --set reference.flashfry_db=~/.cretarget/reference/hg38_standard_cas9ngg_db
```

> **Clusters without outbound network on compute nodes:** all downloads are
> explicit, login-node commands that populate a cache directory. Pipeline steps
> only read from that cache and never download. Point the cache somewhere
> shared with `CRETARGET_CACHE_DIR=/shared/path/.cretarget`.

## Usage

```bash
nextflow run MikeKazz/creTarget \
    -profile slurm,conda \
    --enhancer_bed enhancers.txt \
    --snp_bed snps.txt \
    --fasta ~/.cretarget/reference/hg38.fa \
    --flashfry_jar /path/to/FlashFry-assembly-1.15.jar \
    --flashfry_db ~/.cretarget/reference/hg38_standard_cas9ngg_db \
    --dedup_blast_radius \
    --outdir results
```

At least one of `--snp_bed` / `--enhancer_bed` is required; both may be given
together.

Output is `results/guides.ranked.tsv`: one row per candidate guide, annotated
with all FlashFry scores, genome coordinates, hard-filter results, and `rank`
(1 = best, per target; `exclude (<reason>)` if it failed a filter).
`--run_seed_finder` adds `results/guides.final.tsv`, the same table plus four
seed-screening columns; the ranked table is published either way.

### Common options

| Flag | Default | Description |
|---|---|---|
| `--chunk_size` | 5000 | Regions per parallel task |
| `--queue` | none | Partition/queue name, if your scheduler needs one |
| `--gpu_gres` | none | GPU request for the launch-dCas9 step, e.g. `gpu:1` |
| `--dedup_blast_radius` | off | Flag input targets within 1kb of each other (does not merge) |
| `--blast_radius_bp` | 1000 | CRISPRi silencing window, bp |
| `--flank_bp` | 50 | Enhancer-type: padding around each interval |
| `--snp_window_bp` | 150 | SNP-type: search window each side of the point |
| `--max_mismatch` | 3 | FlashFry off-target search depth |
| `--atac_narrowpeak` | none | narrowPeak file; adds `distance_to_atac_summit` |
| `--run_launch_dcas9` | off | Add launch-dCas9 efficacy scoring |
| `--run_seed_finder` | off | Add seed-finder-lite off-target/TSS columns |
| `--seed_finder_max_guides` | none | Split seed-finder into tasks of <= N guides |

Run `nextflow run MikeKazz/creTarget --help` for the full list.

### Choosing a chunk size

`--chunk_size` trades task count against per-task runtime. Smaller chunks
parallelize better but add scheduling overhead; larger chunks are more
efficient but slower to retry when one fails. The default of 5,000 regions is
a reasonable starting point. On a busy cluster, prefer fewer, larger tasks; on
an idle one, the reverse.

Nextflow throttles concurrent submissions (`executor.queueSize`, default 50 in
the cluster profiles) so that a genome-wide run does not put thousands of JVMs
onto a shared filesystem simultaneously. Raise it if your site can take it.

### Optional: launch-dCas9 on-target scoring

Requires a separate conda environment with `torch` (see
[launch-dCas9](https://github.com/Wancen/launch-dCas9)). Model checkpoints are
fetched from Zenodo on demand rather than bundled:

```bash
cretarget fetch-models

nextflow run MikeKazz/creTarget ... \
    --run_launch_dcas9 \
    --launch_script /path/to/launch-dCas9.py \
    --launch_python /path/to/conda/envs/launch-dCas9/bin/python \
    --ckpt_dir ~/.cretarget/models
```

### Optional: ATAC-seq annotation

Supply your own narrowPeak file with `--atac_narrowpeak`, or fetch one of 19
bundled hematopoietic cell types:

```bash
cretarget fetch-atac --cell-type k562
nextflow run MikeKazz/creTarget ... --atac_narrowpeak ~/.cretarget/atac_peaks/k562.narrowPeak
```

### Optional: seed-finder-lite

Asks a different question from FlashFry: not *where else could this guide
cut*, but *where else could it sit on a promoter and repress*. A perfect
PAM-proximal seed next to any NGG is enough to recruit dCas9, so this screens
every guide's seed against ±2 kb windows around all annotated TSSs. It is
orthogonal to FlashFry's cutting-oriented scores, not a substitute for them.

Requires a separate install (Zenodo DOI
[10.5281/zenodo.21599135](https://doi.org/10.5281/zenodo.21599135)).

```bash
nextflow run MikeKazz/creTarget ... \
    --run_seed_finder \
    --seed_finder_dir /path/to/seed-finder-lite \
    --seed_finder_min_length 8 \
    --seed_finder_tss_window 1000
```

`--seed_finder_min_length` and `--seed_finder_tss_window` have no defaults on
purpose: they change what the screen means, so choose them deliberately.

This step runs **once over the whole guide table**, not per chunk.
seed-finder-lite scans the promoterome a single time and then matches every
guide by dictionary lookup, so its cost barely moves with guide count --
measured at ~25 s for 20 guides and ~25 s for 2,000. Per-chunk invocation paid
that fixed scan on every chunk. The step streams the table rather than loading
it, so running it over a genome-wide table costs time, not memory.
`--seed_finder_max_guides N` splits it into tasks of at most N guides if one
task cannot finish inside its walltime; each shard repeats the scan, so keep
the shard count low.

Before the real run, creTarget asks the installed tool to recover two real
Cas9 sites taken from its own promoter reference, one on each strand. A build
whose plus-strand scan indexes guides under the reverse complement of their
seed reports only minus-strand hits -- it exits 0 and returns plausible
numbers while missing roughly half of all true hits, which no check of
creTarget's own code would catch. If your install fails this, upgrade it; the
error names the cause. `--seed_finder_skip_check` bypasses it, and
`cretarget seed-finder-check --seed-finder-dir ...` runs it on its own.

A `VERSION` file at the root of your seed-finder-lite directory is recorded in
`software_versions.yml`; without one the provenance record says `unknown`.

## How it works

1. **Prep regions** — parse and validate targets, flag linked neighbours, pad
   each region, extract sequence, encode genome coordinates into FASTA headers.
   *(global, once)*
2. **Split** — divide regions into chunks of `--chunk_size`.
3. **FlashFry discover + score** — per chunk, JVM heap derived from the task's
   actual memory allocation.
4. **Annotate** — recover genome coordinates, join region metadata, flag guides
   within the blast radius, compute ATAC summit distance. *(per chunk)*
5. **launch-dCas9** *(optional)* — CNN on-target efficacy. *(per chunk)*
6. **Rank** — hard filters, weighted composite score, per-target rank.
   *(per chunk)*
7. **Collate** — stream per-chunk tables into `guides.ranked.tsv`.
8. **seed-finder-lite** *(optional)* — promoter-proximal seed screening over
   the collated table, producing `guides.final.tsv`. *(once, not per chunk)*

Ranking is per-`region_id`, and every guide for a region lands in the same
chunk, so per-chunk ranking is exactly equivalent to ranking globally. This is
verified in CI: a chunked run produces byte-identical results to an unchunked
one.

## Step-level CLI

`cretarget <step>` runs one step on one file — useful for debugging, for small
one-off jobs, or for composing your own workflow:

```bash
cretarget prep-regions --enhancer-bed enhancers.txt --fasta hg38.fa --outdir prep
cretarget annotate --scored scored.tsv --region-metadata prep/regions.metadata.tsv --out annotated.tsv
cretarget rank --annotated annotated.tsv --out ranked.tsv
cretarget doctor
```

These read defaults from `~/.config/cretarget/config.toml`. Nextflow always
passes parameters explicitly and ignores that file.

## Repository layout

```
main.nf                     pipeline entrypoint
nextflow.config             params, profiles (local/slurm/pbs/sge/lsf, conda/docker/singularity)
conf/base.config            per-process resources and retry escalation
conf/modules.config         publishing rules
conf/limits.config          resource ceilings, applied last
conf/site-template.config   commented blank to copy for a new cluster
conf/example-singlenode.config
                            filled-in worked example (single-node SLURM)
conf/test.config            small-input smoke-test profile
workflows/cretarget.nf      the workflow DAG
modules/local/              one process per step
src/cretarget/              the Python package (steps, config, asset fetching, doctor)
assets/                     cluster survey script, sbatch template, conda env
docs/cluster-setup.md       step-by-step first install on a new cluster
tests/                      pytest suite
```

## Running on your cluster

The built-in profiles are generic. For most sites, this is all you need:

```bash
nextflow run MikeKazz/creTarget -profile slurm,conda \
    --queue short --gpu_gres gpu:1 ...
```

`-profile` takes `local`, `slurm`, `pbs`, `sge` or `lsf` for execution, and
`conda`, `mamba`, `micromamba`, `docker`, `singularity` or `apptainer` for
software.

**The container is the most portable option** -- one image, no dependency
solve, no version skew between users:

```bash
nextflow run MikeKazz/creTarget -r v0.2.6 -profile slurm,singularity ...
```

It defaults to `ghcr.io/mikekazz/cretarget:0.2.6`; override with
`--container_image`. The `conda` profile builds the same pinned environment
from `assets/environment.yml` instead; if conda's solver is slow, use
`-profile micromamba` (same environment, a solver that handles it in seconds). `--queue` and
`--gpu_gres` cover the two things that vary most between sites. The container
profiles additionally require `--container_image`; the pipeline refuses to
start with a container engine enabled and no image, rather than silently
running tasks against whatever tools happen to be on the host.

A **site profile** is only needed when your cluster has constraints the
scheduler cannot express — notably when it does *not* account for memory when
placing jobs, does not enforce the memory it granted, or shares nodes with
work it cannot see. `conf/site-template.config` is a commented blank to copy,
and it lists the specific commands that answer each question.
`conf/example-singlenode.config` is a filled-in worked example for a
single-node SLURM cluster with shared GPUs, where the scheduler accounts for
memory but does not enforce it. Its numbers are measured rather than guessed
and the reasoning is documented inline. Register your copy in the `profiles`
block of `nextflow.config`, or keep it outside the repository and pass it with
`-c /path/to/mysite.config`.

`assets/inspect_cluster.sh` is a read-only survey script that reports
everything a site profile needs — scheduler configuration, filesystem types
and locking, GPU compute capabilities, and any reference data already on disk.

`docs/cluster-setup.md` walks through a first install on a new cluster.

On a **single-node cluster** in particular, submit the Nextflow driver job
with *minimal* resources (1 CPU, a few GB). The driver only waits; every
resource it holds is one the actual work cannot have.
`assets/run_cretarget.sbatch` is a worked driver template, including the
`source conda.sh` step that batch jobs need because `conda activate` does not
work in a non-interactive shell.

## Measured performance

A full ENCODE cCRE run, on a single node (20 physical cores / 40 threads,
125 GB RAM, NFS storage) under SLURM:

| | |
|---|---|
| input | 2,348,854 cCREs (GRCh38, all classes) |
| wall time | ~24 h |
| guides produced | 92,519,254 |
| regions with >=1 usable guide | 2,337,508 (99.52%) |
| guides per region | min 1, max 205, mean 39.4 |
| passing hard filters | 67,850,895 (73.3%) |
| output table | 50 GB (~10 GB gzipped) |
| scheduler jobs | 7,049, zero failures or retries |
| peak memory per task | 2.7 GB (FlashFry), 5.2 GB (region prep) |

Configuration: `chunk_size 1000`, 10 concurrent FlashFry tasks at 6 GB each.
FlashFry runtime is superlinear in regions per chunk (measured T ~ n^1.32), so
smaller chunks are faster -- 5,000 regions per chunk costs 48 min against 6.1
min for 1,000. Concurrency, not chunk size, is the main throughput lever.

Scaling is close to linear in region count, so a 10,000-region design takes
minutes rather than hours.

## Known limitations

- FlashFry (JVM), launch-dCas9 (conda/torch), and seed-finder-lite are external
  dependencies, not bundled.
- ATAC distance is computed and reported but is not part of the ranking formula.
- hg38 is the tested reference. Other genomes should work given a matching
  FlashFry index, but are not validated.

## Citation

creTarget is archived on Zenodo. Cite
[10.5281/zenodo.24148024](https://doi.org/10.5281/zenodo.24148024) for the software in
general (this always resolves to the latest release), or
[10.5281/zenodo.24148025](https://doi.org/10.5281/zenodo.24148025) for v0.2.1
specifically when exact reproducibility matters.

If you use creTarget, please cite it (see `CITATION.cff`) along with
[FlashFry](https://github.com/mckennalab/FlashFry), and — if you use
`--run_launch_dcas9` — [launch-dCas9](https://github.com/Wancen/launch-dCas9),
whose pretrained weights this tool redistributes via Zenodo.

## License

MIT
