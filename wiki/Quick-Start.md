# Quick start

From an empty directory to a ranked guide table. This assumes Nextflow, Java
17+, `bedtools`, `samtools`, and Python 3.10+ are installed.

## 1. Install

```bash
pip install cretarget
```

If the package is not on PyPI yet:

```bash
pip install "https://github.com/MikeKazz/creTarget/archive/refs/tags/v0.5.1.tar.gz"
```

## 2. Fetch the reference

```bash
cretarget setup-reference
```

This downloads the hg38 FASTA and a prebuilt FlashFry index (~3.2 GB) into
`~/.cretarget/reference`. Run it once, on a login node with network access.

```bash
curl -LO https://github.com/mckennalab/FlashFry/releases/download/1.15/FlashFry-assembly-1.15.jar
cretarget configure --set reference.fasta=~/.cretarget/reference/hg38.fa
cretarget configure --set reference.flashfry_db=~/.cretarget/reference/hg38_standard_cas9ngg_db
cretarget configure --set reference.flashfry_jar=$PWD/FlashFry-assembly-1.15.jar
cretarget doctor
```

## 3. Make a tiny input

```bash
printf 'chr1:10000-10400\nchr1:50000-50300\n' > enhancers.txt
printf 'chr1:120000\nchr2:60000\n' > snps.txt
```

- Enhancer intervals: `chr:start-end`, BED convention (0-based, half-open).
- SNP points: `chr:pos`, 1-based.

Malformed lines are skipped with a warning; contigs are validated against the
FASTA `.fai` before anything runs.

## 4. Run

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

Use `-profile local` if the dependencies are already on `PATH`, or
`-profile local,docker` with `--container_image
ghcr.io/mikekazz/cretarget:0.5.1`.

`-r v0.5.1` pins the release. Omit it to run the current `main` branch, which
is less reproducible.

## 5. Read the result

```bash
ls results
head -2 results/guides.ranked.tsv
```

The main table is `results/guides.ranked.tsv`. One row is one candidate guide.
`rank` is an integer within each `region_id` (`1` is best), or
`exclude (<reason>)` when a hard filter removed the guide. `region_id` is the
target, not the guide.

```bash
python - <<'PY'
import pandas as pd
df = pd.read_csv("results/guides.ranked.tsv", sep="\t")
print(df.groupby("region_id")["passes_hard_filters"].agg(["size", "sum"]))
print(df["exclude_reason"].value_counts(dropna=False).head())
PY
```

Keep `run_manifest.json` and `software_versions.yml` with the table: they
record the exact code, parameters, and tool versions.

## 6. What next?

- [Inputs](Inputs) and [Outputs](Outputs) — the details.
- [Parameters](Parameters) — every flag.
- [Seed screening](Seed-Screening) — promoter-proximal off-target screen.
- [Cluster setup](Cluster-Setup) — scale up.
- [Troubleshooting](Troubleshooting) — if something fails.
- [Full quick start](https://mikekazz.github.io/creTarget/getting-started/quickstart/)
