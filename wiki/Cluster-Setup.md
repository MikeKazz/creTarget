# Cluster setup

creTarget runs on a laptop and on a shared cluster; the scheduler is a profile,
not a separate code path. This page is the first-run checklist. The detailed
walkthrough, including a survey script and a single-node example, is in the
[full cluster guide](https://mikekazz.github.io/creTarget/cluster-setup/).

## 1. Profiles

Combine an execution profile with a software profile:

```bash
nextflow run MikeKazz/creTarget -r v0.5.1 -profile slurm,conda ...
```

| Execution | Software |
|---|---|
| `local`, `standard` | `conda`, `mamba`, `micromamba` |
| `slurm` | `docker` |
| `pbs` (PBS Pro) | `singularity` |
| `sge` | `apptainer` |
| `lsf` | — |

The container profiles use `ghcr.io/mikekazz/cretarget:0.5.1` by default. The
conda profiles build the same pinned environment from `assets/environment.yml`;
`mamba` and `micromamba` only change the solver.

## 2. Reference data on shared storage

Compute nodes usually have no outbound network. Run every fetch on a login node
with the cache pointed at shared storage:

```bash
export CRETARGET_CACHE_DIR=/shared/path/.cretarget
cretarget setup-reference
cretarget fetch-atac --cell-type k562     # optional
cretarget fetch-models                    # optional, launch-dCas9
cretarget build-seed-reference --fasta /shared/path/.cretarget/reference/hg38.fa
```

Pass absolute paths to the pipeline. It never downloads mid-run.

## 3. Survey the cluster

```bash
bash assets/inspect_cluster.sh
```

Read-only. It reports scheduler configuration, filesystem types and locking,
GPU capabilities, and reference data already on disk.

## 4. Submit the driver as a batch job

Keep the Nextflow driver request **minimal** (1 CPU, a few GB): it only waits,
and every resource it holds is one the real work cannot have. On a single-node
cluster, asking for 8 CPUs and 80 GB can deadlock while waiting for jobs that
can never be scheduled.

`assets/run_cretarget.sbatch` is a worked template. It includes the
`source .../etc/profile.d/conda.sh` step that batch jobs need, because
`conda activate` does not work in a non-interactive shell.

## 5. Resource ceilings

| Parameter | Default | Meaning |
|---|---:|---|
| `--max_cpus` | 16 | Hard ceiling after profiles and retries |
| `--max_memory` | 128 GB | Hard ceiling |
| `--max_time` | 48 h | Hard ceiling |
| `maxRetries` | 2 | Retries for OOM-killed, timed-out, or scheduler-limited tasks |

Requests escalate with `task.attempt`; a retry gets more memory, CPUs, and time
up to the ceiling. FlashFry's JVM heap is derived from the task's actual memory
allocation, so a 2× retry gets a 2× heap. If a request escalates past what a
partition will grant, the task is rejected rather than queued. Lower the
ceilings in a site profile to the largest single allocation the partition
allows.

## 6. Site profiles

The built-in profiles are generic. A site profile is needed only when the
cluster has constraints the scheduler cannot express.

- `conf/site-template.config` is a commented blank to copy.
- `conf/example-singlenode.config` is a worked single-node SLURM example with
  measured numbers and inline reasoning.
- Keep a site profile outside the repository and pass `-c /path/to/mysite.config`
  to avoid carrying a local edit in git.

## 7. Resume

Add `-resume` to every invocation. A run that dies at hour 20 picks up from the
last completed chunk. The cache is keyed by inputs and parameters, so changing
a parameter invalidates only the tasks it affects. Keep the same work directory
(`-w`) if you move the run.

## 8. Storage

For a genome-wide cCRE run:

| Item | Size |
|---|---|
| Reference (download + decompressed FASTA) | ~6.5 GB |
| Work directory | ~300 GB |
| Final output | ~50 GB (~10 GB gzipped) |

The work directory is the big one. Do not delete it until the run is complete
and the output is archived, because `-resume` depends on it.

## 9. Measured performance

A full ENCODE cCRE run on one node (20 physical cores, 125 GB RAM, NFS) under
SLURM:

| | |
|---|---|
| Input | 2,348,854 cCREs |
| Wall time | ~24 h |
| Guides | 92,519,254 |
| Regions with ≥1 usable guide | 2,337,508 (99.52%) |
| Passing hard filters | 67,850,895 (73.3%) |
| Output table | 50 GB (~10 GB gzipped) |
| Scheduler jobs | 7,049, zero failures or retries |
| Peak memory per task | 2.7 GB (FlashFry), 5.2 GB (region prep) |

Configuration: `chunk_size 1000`, 10 concurrent FlashFry tasks at 6 GB each.
Scaling is close to linear in region count; concurrency, not chunk size, is the
main throughput lever.

## Common first-run problems

| Symptom | Fix |
|---|---|
| `conda activate` fails in a batch job | Source `.../etc/profile.d/conda.sh` first |
| Singularity uses a stale cached image | Override `NXF_SINGULARITY_CACHEDIR` |
| `flock` or work-directory errors on NFS | Put the work directory on a filesystem that supports locking |
| Tasks rejected rather than queued | Lower `max_cpus`, `max_memory`, or `max_time` |
| `-resume` reuses nothing | Keep the same `-w` and input paths; parameters are part of the key |
