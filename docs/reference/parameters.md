# Pipeline parameters

Parameters are passed to the pipeline with double underscores
(`--snp_bed`). The step-level CLI uses hyphens (`--snp-bed`); see
[Step-level CLI](cli.md).

## Input

| Parameter | Default | Description |
|---|---|---|
| `--snp_bed` | none | SNP targets: `chr:pos` per line (1-based) or BED (0-based, half-open). |
| `--enhancer_bed` | none | Enhancer targets: `chr:start-end` (BED convention) or BED. |

At least one is required. Both may be supplied in the same run.

## Reference

| Parameter | Default | Description |
|---|---|---|
| `--fasta` | none | hg38 FASTA. A `.fai` alongside it is required. |
| `--flashfry_jar` | none | FlashFry assembly JAR (e.g. `FlashFry-assembly-1.15.jar`). |
| `--flashfry_db` | none | FlashFry off-target index. Its `.header` must sit alongside. |

All three are required. The pipeline validates the `.fai` and `.header` before
launching any task.

## Region preparation

| Parameter | Default | Description |
|---|---|---|
| `--flank_bp` | `50` | Padding around each enhancer interval, both sides. |
| `--snp_window_bp` | `150` | Search window each side of a SNP point. |
| `--blast_radius_bp` | `1000` | CRISPRi silencing window, used for `within_blast_radius` and linked-target flagging. |
| `--dedup_blast_radius` | `false` | Flag same-contig targets within one blast radius as `linked_regions`. Flags only; does not merge or remove. |

## FlashFry

| Parameter | Default | Description |
|---|---|---|
| `--chunk_size` | `5000` | Regions per parallel task. |
| `--max_mismatch` | `3` | Off-target search depth passed to FlashFry; also used to normalize the closest-hit distance in ranking. |
| `--scoring_metrics` | see below | Comma-separated FlashFry scoring metrics. |

The default scoring metrics are:

```text
doench2014ontarget,doench2016cfd,dangerous,hsu2013,moreno2015,jostandsantos,minot,folding
```

Change this only if you know which columns your downstream analysis expects.
The ranking step requires `JostCRISPRi_specificityscore`, `Hsu2013`,
`DoenchCFD_specificityscore`, and `basesDiffToClosestHit`.

## Optional annotation

| Parameter | Default | Description |
|---|---|---|
| `--atac_narrowpeak` | none | 10-column narrowPeak file; adds `distance_to_atac_summit`. Informational, not used in ranking. |

## Optional: launch-dCas9

| Parameter | Default | Description |
|---|---|---|
| `--run_launch_dcas9` | `false` | Add CNN-based efficacy scoring. |
| `--launch_script` | none | Path to `launch-dCas9.py`. Required when enabled. |
| `--launch_python` | none | Python from the torch environment. |
| `--ckpt_dir` | none | Checkpoint directory. Fetch first with `cretarget fetch-models`. Required when enabled. |

## Optional: seed screening

| Parameter | Default | Description |
|---|---|---|
| `--run_seed_finder` | `false` | Add promoter-proximal seed-screening columns. |
| `--seed_finder_min_length` | **none** | Minimum exact-match seed length. Required when enabled; choose deliberately. |
| `--seed_finder_tss_window` | **none** | Maximum distance from a TSS for a hit. Required when enabled; the reference extends 2000 bp. |
| `--seed_finder_promoter_fa` | cached build | Override the promoter-window FASTA. |
| `--seed_finder_tss_table` | cached build | Override the TSS table. |
| `--seed_finder_max_guides` | one task | Split the screen into tasks of at most N guides. ~1 GB RSS per million guides. |
| `--seed_finder_all_guides` | `false` | Screen guides that failed the hard filters too. |

## Output

| Parameter | Default | Description |
|---|---|---|
| `--outdir` | `results` | Output directory. |
| `--publish_dir_mode` | `copy` | Nextflow publish mode: `copy`, `symlink`, or `rellink`. |

## Software provisioning

| Parameter | Default | Description |
|---|---|---|
| `--container_image` | `ghcr.io/mikekazz/cretarget:0.5.1` | Image used by the `docker`, `singularity`, and `apptainer` profiles. Override for your own build. |

The pipeline refuses to start if a container profile is enabled and no image is
given.

## Scheduler and resources

| Parameter | Default | Description |
|---|---|---|
| `--queue` | none | Partition/queue name, when the site needs one. |
| `--gpu_gres` | none | GPU request for launch-dCas9, e.g. `gpu:1` (SLURM `--gres`). |
| `--gpu_count` | none | Local executor only: pin tasks round-robin across N GPUs. |
| `--max_cpus` | `16` | Resource ceiling after retries and profile overrides. |
| `--max_memory` | `128.GB` | Resource ceiling. |
| `--max_time` | `48.h` | Resource ceiling. |

See [Scaling and clusters](../guides/scaling-and-clusters.md) for how retries
escalate and why the ceilings matter.

## Profiles

| Profile | Purpose |
|---|---|
| `local`, `standard` | Run on the current machine. |
| `slurm` | SLURM executor. |
| `pbs` | PBS Pro executor. |
| `sge` | SGE executor. |
| `lsf` | LSF executor. |
| `singlenode` | Worked single-node SLURM example. |
| `conda` | Build the task environment from `assets/environment.yml`. |
| `mamba` | Same environment, mamba solver. |
| `micromamba` | Same environment, micromamba solver. |
| `docker` | Run tasks in the container image. |
| `singularity` | Run tasks in the container image. |
| `apptainer` | Run tasks in the container image. |
| `test` | Small resources for a smoke test. Requires reference paths. |

Combine execution and software profiles with a comma:

```bash
nextflow run MikeKazz/creTarget -profile slurm,micromamba ...
```

## Special flags

| Flag | Behavior |
|---|---|
| `--help` | Print the pipeline help and exit. |
| `--version` | Print the pipeline version and exit. |

Nextflow's own flags are also available and useful: `-resume`, `-w` (work
directory), `-c` (additional config), `-params-file`, `-with-report`,
`-with-trace`, `-with-dag`, and `-r` (revision).

## Parameter files

For a long run, keep parameters in a YAML or JSON file:

```yaml
# params.yml
snp_bed: /data/snps.txt
enhancer_bed: /data/enhancers.txt
fasta: /shared/reference/hg38.fa
flashfry_jar: /shared/tools/FlashFry-assembly-1.15.jar
flashfry_db: /shared/reference/hg38_standard_cas9ngg_db
outdir: /data/results
```

```bash
nextflow run MikeKazz/creTarget -r v0.5.1 -profile slurm,conda -params-file params.yml
```

Command-line parameters override the params file, which overrides profile
defaults. Record the final values in `run_manifest.json`.
