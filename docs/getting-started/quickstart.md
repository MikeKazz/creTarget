# Quick start

This is the shortest path from an empty directory to a ranked guide table. It
assumes a machine with Nextflow, Java 17+, `bedtools`, `samtools`, and Python
3.10+ already installed. On a cluster, read
[Scaling and clusters](../guides/scaling-and-clusters.md) first.

## 1. Install creTarget

```bash
pip install cretarget
```

If the package is not on PyPI yet:

```bash
pip install "https://github.com/MikeKazz/creTarget/archive/refs/tags/v0.5.1.tar.gz"
```

See [Installation](installation.md) for containers and conda profiles.

## 2. Fetch the reference data

This downloads the hg38 FASTA and a prebuilt FlashFry index (~3.2 GB) into
`~/.cretarget/reference`. Run it once, on a login node with network access:

```bash
cretarget setup-reference
```

It prints the exact paths. Record them for the step-level CLI:

```bash
cretarget configure --set reference.fasta=~/.cretarget/reference/hg38.fa
cretarget configure --set reference.flashfry_db=~/.cretarget/reference/hg38_standard_cas9ngg_db
```

FlashFry itself is a separate JVM tool:

```bash
curl -LO https://github.com/mckennalab/FlashFry/releases/download/1.15/FlashFry-assembly-1.15.jar
cretarget configure --set reference.flashfry_jar=$PWD/FlashFry-assembly-1.15.jar
```

Check everything:

```bash
cretarget doctor
```

## 3. Make a tiny input

Two file types are accepted, and they can be used together:

```bash
printf 'chr1:10000-10400\nchr1:50000-50300\n' > enhancers.txt
printf 'chr1:120000\nchr2:60000\n' > snps.txt
```

- Enhancer intervals use `chr:start-end` and follow the BED convention
  (0-based, half-open).
- SNP points use `chr:pos` and are 1-based.

Malformed lines are skipped with a warning. Contigs are validated against the
FASTA `.fai` before anything else runs. See
[Input formats](../guides/inputs.md) for BED, naming, and duplicate handling.

## 4. Run the pipeline

```bash
nextflow run MikeKazz/creTarget -r v0.5.1 \
    -profile local,conda \
    --enhancer_bed enhancers.txt \
    --snp_bed snps.txt \
    --fasta ~/.cretarget/reference/hg38.fa \
    --flashfry_jar "$PWD/FlashFry-assembly-1.15.jar" \
    --flashfry_db ~/.cretarget/reference/hg38_standard_cas9ngg_db \
    --outdir results
```

The `local,conda` profile runs on the current machine and builds the task
environment from `assets/environment.yml`. Use `-profile local` if the
dependencies are already on `PATH`, or `-profile local,docker` with
`--container_image ghcr.io/mikekazz/cretarget:0.5.1`.

`-r v0.5.1` pins the release. Omit it to run the current `main` branch, which is
less reproducible.

## 5. Read the result

```bash
ls results
head -2 results/guides.ranked.tsv
```

The main table is `results/guides.ranked.tsv`. One row is one candidate guide.
The `rank` column is an integer within each `region_id` (`1` is best), or
`exclude (<reason>)` when a hard filter removed the guide. `region_id` is the
target, not the guide.

Useful first checks:

```bash
# how many guides per region?
python - <<'PY'
import pandas as pd
df = pd.read_csv("results/guides.ranked.tsv", sep="\t")
print(df.groupby("region_id")["passes_hard_filters"].agg(["size", "sum"]))
print(df["exclude_reason"].value_counts(dropna=False).head())
PY
```

The `run_manifest.json` and `software_versions.yml` files in the output
directory record the exact code, parameters, and tool versions. Keep them with
the table.

## 6. What next?

- **Add optional scoring.** [launch-dCas9](../guides/launch-dcas9.md) for CNN
  on-target efficacy, [ATAC-seq](../guides/atac-seq.md) for accessibility
  annotation, and [seed screening](../guides/seed-screening.md) for
  promoter-proximal off-target repression.
- **Scale up.** [Scaling and clusters](../guides/scaling-and-clusters.md)
  covers chunk size, site profiles, SLURM/PBS/SGE/LSF, and genome-wide runs.
- **Understand the output.** [Output files](../guides/outputs.md) and
  [Output columns](../reference/output-columns.md).
- **Tune the run.** [Pipeline parameters](../reference/parameters.md).
- **Something failed.** [Troubleshooting](../guides/troubleshooting.md).
