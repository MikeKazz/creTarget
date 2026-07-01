# Building the hosted reference bundle

This is a one-time, maintainer-side step (not run by end users or the pipeline
itself). It produces the `genome.fa` + `flashfry_db` bundle that
`SETUP_REFERENCE` downloads automatically on a user's first run.

## 1. Prepare the FASTA: main chromosomes only, no alt scaffolds, no chrM

Start from a full hg38 reference (e.g. UCSC or Ensembl primary assembly) and
filter down to `chr1`-`chr22`, `chrX`, `chrY` only — excluding alt/random/unplaced
scaffolds and `chrM`, matching the `hg38_standard_cas9ngg_db` convention used
previously.

```bash
samtools faidx genome_full.fa
CHROMS=$(for i in $(seq 1 22) X Y; do echo chr${i}; done)
samtools faidx genome_full.fa ${CHROMS} > genome.fa
samtools faidx genome.fa
```

## 2. Build the FlashFry index

```bash
mkdir -p tmp
java -Xmx4g -jar FlashFry-assembly-1.15.jar \
    index \
    --tmpLocation ./tmp \
    --database flashfry_db \
    --reference genome.fa \
    --enzyme spcas9ngg
```

## 3. Package and version

```bash
tar czf hg38-standard-noalt-nochrM-v1.tar.gz genome.fa genome.fa.fai flashfry_db
sha256sum hg38-standard-noalt-nochrM-v1.tar.gz > hg38-standard-noalt-nochrM-v1.tar.gz.sha256
```

Upload both files to the chosen hosting location (Zenodo / GitHub Releases / S3)
under a version tag matching `params.reference_version` in `nextflow.config`
(currently `hg38-standard-noalt-nochrM-v1`).

## 4. Bump the version tag when the bundle changes

If the reference or index is ever rebuilt (new FlashFry version, different
filtering, etc.), publish it under a **new** version string and update
`params.reference_version` — never silently overwrite an existing tag, since
users' local caches key off this string.
