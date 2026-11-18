# Parameters

Pipeline parameters use double underscores (`--snp_bed`). The step-level CLI
uses hyphens (`--snp-bed`).

## Input

| Parameter | Default | Description |
|---|---|---|
| `--snp_bed` | none | SNP targets: `chr:pos` (1-based) or BED (0-based, half-open) |
| `--enhancer_bed` | none | Enhancer targets: `chr:start-end` or BED |

At least one is required; both may be used in one run.

## Reference

| Parameter | Default | Description |
|---|---|---|
| `--fasta` | none | hg38 FASTA; a `.fai` alongside it is required |
| `--flashfry_jar` | none | FlashFry assembly JAR |
| `--flashfry_db` | none | FlashFry index; its `.header` must sit alongside |

All three are required. The pipeline validates the `.fai` and `.header` before
launching any task.

## Region preparation

| Parameter | Default | Description |
|---|---|---|
| `--flank_bp` | `50` | Padding around each enhancer, both sides |
| `--snp_window_bp` | `150` | Search window each side of a SNP |
| `--blast_radius_bp` | `1000` | CRISPRi silencing window for `within_blast_radius` and link flagging |
| `--dedup_blast_radius` | `false` | Flag same-contig neighbours as `linked_regions`; does not merge or remove |

## FlashFry

| Parameter | Default | Description |
|---|---|---|
| `--chunk_size` | `5000` | Regions per parallel task |
| `--max_mismatch` | `3` | Off-target search depth; also normalizes the closest-hit distance in ranking |
| `--scoring_metrics` | see below | Comma-separated FlashFry scoring metrics |

Default metrics:

```text
doench2014ontarget,doench2016cfd,dangerous,hsu2013,moreno2015,jostandsantos,minot,folding
```

The ranking step requires `JostCRISPRi_specificityscore`, `Hsu2013`,
`DoenchCFD_specificityscore`, and `basesDiffToClosestHit`.

## Optional annotation

| Parameter | Default | Description |
|---|---|---|
| `--atac_narrowpeak` | none | 10-column narrowPeak; adds `distance_to_atac_summit` (not ranked) |

## launch-dCas9

| Parameter | Default | Description |
|---|---|---|
| `--run_launch_dcas9` | `false` | Add CNN efficacy scoring |
| `--launch_script` | none | Path to `launch-dCas9.py` |
| `--launch_python` | none | Python from the torch environment |
| `--ckpt_dir` | none | Checkpoint directory; fetch with `cretarget fetch-models` |

## Seed screening

| Parameter | Default | Description |
|---|---|---|
| `--run_seed_finder` | `false` | Add seed-screening columns |
| `--seed_finder_min_length` | **none** | Minimum exact-match seed length; required when enabled |
| `--seed_finder_tss_window` | **none** | Maximum distance from a TSS; required when enabled; reference extends 2000 bp |
| `--seed_finder_promoter_fa` | cached build | Override the promoter-window FASTA |
| `--seed_finder_tss_table` | cached build | Override the TSS table |
| `--seed_finder_max_guides` | one task | Split the screen into tasks of at most N guides (~1 GB RSS per million) |
| `--seed_finder_all_guides` | `false` | Screen guides that failed the hard filters too |

`--seed_finder_max_guides` and `--seed_finder_all_guides` are rejected when
`--run_seed_finder` is off, rather than silently ignored.

## Output

| Parameter | Default | Description |
|---|---|---|
| `--outdir` | `results` | Output directory |
| `--publish_dir_mode` | `copy` | `copy`, `symlink`, or `rellink` |

## Software provisioning

| Parameter | Default | Description |
|---|---|---|
| `--container_image` | `ghcr.io/mikekazz/cretarget:0.5.1` | Image for the `docker`, `singularity`, and `apptainer` profiles |

A container profile with no image is a hard error.

## Scheduler and resources

| Parameter | Default | Description |
|---|---|---|
| `--queue` | none | Partition/queue name, if the site needs one |
| `--gpu_gres` | none | GPU request for launch-dCas9, e.g. `gpu:1` (SLURM `--gres`) |
| `--gpu_count` | none | Local executor only: pin tasks round-robin across N GPUs |
| `--max_cpus` | `16` | Resource ceiling after retries and profile overrides |
| `--max_memory` | `128.GB` | Resource ceiling |
| `--max_time` | `48.h` | Resource ceiling |

## Profiles

| Profile | Purpose |
|---|---|
| `local`, `standard` | Run on the current machine |
| `slurm` | SLURM executor |
| `pbs` | PBS Pro executor |
| `sge` | SGE executor |
| `lsf` | LSF executor |
| `singlenode` | Worked single-node SLURM example |
| `conda` | Build the task environment from `assets/environment.yml` |
| `mamba`, `micromamba` | Same environment, faster solver |
| `docker`, `singularity`, `apptainer` | Run tasks in the container image |
| `test` | Small resources for a smoke test; reference paths still required |

Combine execution and software profiles with a comma:

```bash
nextflow run MikeKazz/creTarget -profile slurm,micromamba ...
```

## Special flags

| Flag | Behavior |
|---|---|
| `--help` | Print the pipeline help and exit |
| `--version` | Print the pipeline version and exit |

Useful Nextflow flags: `-resume`, `-w`, `-c`, `-params-file`, `-with-report`,
`-with-trace`, `-with-dag`, and `-r`.

## Parameter files

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
defaults. The final values are recorded in `run_manifest.json`.

Full reference: [Pipeline parameters](https://mikekazz.github.io/creTarget/reference/parameters/)
