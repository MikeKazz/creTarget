# Step-level CLI

The `cretarget` command runs **one step on one file**. It is not a pipeline
driver: Nextflow owns splitting, scheduling, retries, and resume. Use the CLI to
debug a step, run a small one-off job, or compose your own workflow.

```text
cretarget --version
cretarget --help
cretarget <command> --help
```

Every command writes progress to stderr and its table to the path given by
`--out` / `--outdir`. Paths read from the config file are expanded for `~` and
`$VARS`.

## Configuration

The CLI reads `~/.config/cretarget/config.toml` (override the directory with
`CRETARGET_CONFIG_DIR`). Nextflow does **not** read this file; pipeline runs
always pass parameters explicitly.

```bash
cretarget configure
cretarget configure --set reference.fasta=~/.cretarget/reference/hg38.fa
cretarget configure --set reference.flashfry_db=~/.cretarget/reference/hg38_standard_cas9ngg_db
cretarget configure --set reference.flashfry_jar=/path/to/FlashFry-assembly-1.15.jar
```

`--set` is repeatable and takes `section.key=value`. The supported sections are
`reference`, `launch_dcas9`, `seed_finder`, and `cache`.

## Data commands

### `setup-reference`

Download hg38 and the prebuilt FlashFry index (~3.2 GB), verify checksums,
decompress the FASTA, build the `.fai`, and extract the index.

```bash
cretarget setup-reference [--outdir DIR] [--force]
```

### `fetch-models`

Download and verify the launch-dCas9 checkpoints.

```bash
cretarget fetch-models [--outdir DIR] [--force]
```

### `fetch-atac`

Download one or more ATAC-seq peak files.

```bash
cretarget fetch-atac --cell-type k562 [--cell-type monocyte] [--outdir DIR] [--force]
```

`--cell-type` is required and repeatable. See
[ATAC-seq annotation](../guides/atac-seq.md) for the valid names.

### `build-seed-reference`

Build the TSS table and promoter-window FASTA for seed screening.

```bash
cretarget build-seed-reference \
    --fasta ~/.cretarget/reference/hg38.fa \
    [--gtf my.gtf] \
    [--gencode-release 47] \
    [--outdir DIR] \
    [--force]
```

`--gtf` overrides `--gencode-release`. The release, source, and checksum are
recorded in `provenance.json`.

## Processing commands

### `prep-regions`

Parse targets, validate contigs, make names unique, pad regions, and extract
sequence with `bedtools getfasta`.

```bash
cretarget prep-regions \
    --snp-bed snps.txt \
    --enhancer-bed enhancers.txt \
    --fasta hg38.fa \
    --flank-bp 50 \
    --snp-window-bp 150 \
    --blast-radius-bp 1000 \
    --dedup-blast-radius \
    --outdir prep
```

| Flag | Default | Description |
|---|---|---|
| `--snp-bed` | none | SNP input (`chr:pos` or BED). |
| `--enhancer-bed` | none | Enhancer input (`chr:start-end` or BED). |
| `--fasta` | config | Reference FASTA with `.fai`. |
| `--flank-bp` | 50 | Enhancer padding. |
| `--snp-window-bp` | 150 | SNP search window. |
| `--blast-radius-bp` | 1000 | Silencing window for link flagging. |
| `--dedup-blast-radius` | off | Flag linked targets. |
| `--outdir` | `.` | Output directory. |

Outputs: `regions.fasta`, `regions.metadata.tsv`, and an intermediate BED.

### `annotate`

Recover genome coordinates from FlashFry's `contig` field, join region
metadata, flag guides inside the blast radius, and optionally compute ATAC
summit distance.

```bash
cretarget annotate \
    --scored scored.tsv \
    --region-metadata prep/regions.metadata.tsv \
    --blast-radius-bp 1000 \
    [--atac-narrowpeak peaks.narrowPeak] \
    --out annotated.tsv
```

`--scored` accepts one or more files. Output defaults to
`guides.annotated.tsv`.

### `rank`

Apply hard filters, compute the composite score, and rank within each
`region_id`.

```bash
cretarget rank --annotated annotated.tsv --max-mismatch 3 --out ranked.tsv
```

Output defaults to `guides.ranked.tsv`. See
[Output files](../guides/outputs.md) for the ranking formula and filters.

### `score-launch-dcas9`

Add the two CNN efficacy columns.

```bash
cretarget score-launch-dcas9 \
    --annotated annotated.tsv \
    --launch-script /path/to/launch-dCas9.py \
    --ckpt-dir ~/.cretarget/models \
    [--python-bin /path/to/torch/env/bin/python] \
    [--workdir scratch] \
    --out with_launch.tsv
```

`--launch-script` falls back to `[launch_dcas9].script` and `--ckpt-dir` to the
cache. `--python-bin` defaults to the interpreter running creTarget.

### `seed-finder`

Add promoter-proximal seed-screening columns and write the hit table.

```bash
cretarget seed-finder \
    --ranked guides.ranked.tsv \
    --min-seed-length 8 \
    --tss-window 1000 \
    [--promoter-fa promoter_windows.fa.gz] \
    [--tss-table tss_table.tsv] \
    [--all-guides] \
    [--id-prefix shard1] \
    --outdir seed_finder \
    --raw-hits hits.csv \
    --out guides.final.tsv
```

`--min-seed-length` and `--tss-window` are required; there are no defaults. See
[Seed screening](../guides/seed-screening.md).

### `coords-to-bed`

Extract unique `chr:start-end` values from a table column into BED. Useful for
turning a previous results table back into pipeline input.

```bash
cretarget coords-to-bed \
    --input results.tsv \
    --column target_enhancer \
    --out regions.bed \
    [--delimiter $'\t']
```

## Diagnostics

### `doctor`

Check tools, configuration, reference data, optional cached assets, and cache
writability.

```bash
cretarget doctor [--shared-path /shared/work]
```

Exit code is 1 if any required check fails. `--shared-path` additionally checks
that a directory exists and is writable; it cannot verify that compute nodes see
it.

## Naming

The step-level CLI uses hyphens (`--snp-bed`); Nextflow uses underscores
(`--snp_bed`). They map one-to-one. The pipeline passes values explicitly; it
does not read the CLI config file.
