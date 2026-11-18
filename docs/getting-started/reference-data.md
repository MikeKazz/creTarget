# Reference data

creTarget does not bundle the genome, the FlashFry index, or any model
checkpoints. They are large, versioned, and have their own terms. Fetching is
always an explicit command, so compute nodes without network access never try
to download mid-run.

## What `setup-reference` does

```bash
cretarget setup-reference
```

It populates `~/.cretarget/reference` (or `$CRETARGET_CACHE_DIR/reference`)
with:

| File | Purpose |
|---|---|
| `hg38.fa` | hg38 FASTA, UCSC `chr` naming |
| `hg38.fa.fai` | samtools index |
| `hg38_standard_cas9ngg_db` | FlashFry off-target index |
| `hg38_standard_cas9ngg_db.header` | FlashFry index header; must sit next to the database |
| `*.sha256` | Checksums downloaded with the payloads |

The payloads come from Zenodo record
[10.5281/zenodo.21115562](https://doi.org/10.5281/zenodo.21115562). The command
verifies each checksum before use, writes to a `.part` file first, and
decompresses the FASTA to plain text because some bedtools builds mishandle
bgzipped input.

It requires `bgzip` (htslib), `samtools`, and enough disk space for both the
compressed download and the decompressed FASTA. Budget about 3.2 GB for the
download plus 3.2 GB for the decompressed genome, and more for the work
directory and results during a run.

## Record the paths

The step-level CLI reads `~/.config/cretarget/config.toml`:

```bash
cretarget configure --set reference.fasta=~/.cretarget/reference/hg38.fa
cretarget configure --set reference.flashfry_db=~/.cretarget/reference/hg38_standard_cas9ngg_db
```

The Nextflow pipeline does not read that file. Pass `--fasta`,
`--flashfry_jar`, and `--flashfry_db` explicitly, or put them in a params file:

```bash
nextflow run MikeKazz/creTarget -params-file params.yml ...
```

```yaml
# params.yml
fasta: /shared/reference/hg38.fa
flashfry_jar: /shared/tools/FlashFry-assembly-1.15.jar
flashfry_db: /shared/reference/hg38_standard_cas9ngg_db
```

## FlashFry

FlashFry is a separate download. Use release 1.15 or a version whose index is
compatible with your database:

```bash
curl -LO https://github.com/mckennalab/FlashFry/releases/download/1.15/FlashFry-assembly-1.15.jar
```

To build an index for another genome, follow FlashFry's instructions. The shape
of the command is:

```bash
samtools faidx mygenome.fa
mkdir -p index_tmp
java -Xmx32g -jar FlashFry-assembly-1.15.jar index \
    --tmpLocation index_tmp \
    --database mygenome_cas9ngg_db \
    --reference mygenome.fa \
    --enzyme spcas9ngg
```

Indexing a human genome takes roughly an hour and wants about 32 GB of RAM. You
need both `mygenome_cas9ngg_db` and `mygenome_cas9ngg_db.header`; the pipeline
checks for the header at launch.

## Optional: model checkpoints

For `--run_launch_dcas9`:

```bash
cretarget fetch-models
```

This downloads the two launch-dCas9 checkpoints from Zenodo record
[10.5281/zenodo.22729308](https://doi.org/10.5281/zenodo.22729308), verifies
their SHA-256 hashes, and prints the attribution. The weights remain the work
of the launch-dCas9 authors.

## Optional: ATAC-seq peaks

```bash
cretarget fetch-atac --cell-type k562
```

Available cell types are `b_cell`, `basophil`, `cd34_bone_marrow`,
`cd34_cord_blood`, `cd4_tcell`, `cd8_tcell`, `clp`, `cmp`, `erythroid`, `gmp`,
`hsc`, `k562`, `lmpp`, `megakaryocyte`, `mep`, `monocyte`, `mpp`,
`neutrophil`, and `nk_cell`. The files come from Zenodo record
[10.5281/zenodo.21268805](https://doi.org/10.5281/zenodo.21268805) and are
checked against the record's checksum manifest.

You can also pass your own 10-column narrowPeak file with
`--atac_narrowpeak`; nothing is downloaded in that case.

## Optional: promoter reference for seed screening

This one is **built**, not downloaded, from a GENCODE annotation and the genome
FASTA you already have:

```bash
cretarget build-seed-reference --fasta ~/.cretarget/reference/hg38.fa
```

It writes `tss_table.tsv`, `promoter_windows.fa.gz`, and `provenance.json` to
`~/.cretarget/seed_reference`. The default annotation is GENCODE release 47;
pin a different release with `--gencode-release N`, or use your own GTF with
`--gtf`. The release, source URL, and checksum are recorded in
`provenance.json`.

Rebuild with `--force` after changing the annotation. See
[Seed screening](../guides/seed-screening.md) for what the reference contains
and how it is used.

## Using another genome

The pipeline itself is genome-agnostic. It needs:

1. A FASTA with a `.fai` index.
2. A FlashFry index built from the same FASTA.
3. Contig names in the input that match the FASTA exactly.

hg38 is the tested reference. Other genomes should work with a matching index,
but are not validated in CI. The seed reference builder currently expects
standard human chromosomes (`chr1`–`chr22`, `chrX`, `chrY`) and a
GENCODE/Ensembl-style GTF.

## Cache layout

```
~/.cretarget/                 # or $CRETARGET_CACHE_DIR
├── reference/
│   ├── hg38.fa
│   ├── hg38.fa.fai
│   ├── hg38_standard_cas9ngg_db
│   ├── hg38_standard_cas9ngg_db.header
│   └── gencode.v47.annotation.gtf.gz
├── models/
│   ├── ckpt-seq-enhancerFitness.pth
│   └── ckpt-seq-promoterFitness.pth
├── atac_peaks/
│   └── k562.narrowPeak
└── seed_reference/
    ├── tss_table.tsv
    ├── promoter_windows.fa.gz
    └── provenance.json
```

## Offline clusters

Run every fetch command on a login node with the cache pointed at shared
storage:

```bash
export CRETARGET_CACHE_DIR=/shared/path/.cretarget
cretarget setup-reference
cretarget fetch-atac --cell-type k562
cretarget fetch-models
cretarget build-seed-reference --fasta /shared/path/.cretarget/reference/hg38.fa
```

Compute nodes then need only read access. If a file is missing, the pipeline
fails with the command to run, rather than downloading.
