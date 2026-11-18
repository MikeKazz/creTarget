---
description: Design, score, and rank CRISPRi guide RNAs for noncoding target regions.
---

# creTarget

**Design, score, and rank CRISPRi guide RNAs for noncoding target regions —
SNPs, enhancers, or arbitrary intervals — in hg38, from a handful of targets to
genome-wide libraries.**

[![CI](https://github.com/MikeKazz/creTarget/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/MikeKazz/creTarget/actions/workflows/ci.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24148024.svg)](https://doi.org/10.5281/zenodo.24148024)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/MikeKazz/creTarget/blob/main/LICENSE)
[![Container](https://img.shields.io/badge/container-ghcr.io-blue)](https://github.com/MikeKazz/creTarget/pkgs/container/cretarget)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue)](https://github.com/MikeKazz/creTarget/blob/main/pyproject.toml)
[![Nextflow ≥ 24.04](https://img.shields.io/badge/Nextflow-%E2%89%A524.04-23aa62)](https://www.nextflow.io/)

<div class="hero" markdown>
creTarget is a Nextflow pipeline and Python toolkit for designing CRISPR
interference (CRISPRi) libraries against the noncoding genome. It takes SNP
points or enhancer intervals, extracts padded target windows, discovers and
scores candidate guides with FlashFry, annotates them against optional
chromatin and CNN-based efficacy data, and returns one ranked table with
provenance attached.
</div>

## Why creTarget

- **Two input types, one workflow.** SNP points (`chr:pos`) and enhancer
  intervals (`chr:start-end`) can be mixed in a single run. BED is accepted;
  malformed lines are skipped with a warning instead of aborting a long job.
- **Flat memory at any scale.** Only region preparation is global. Everything
  after it runs per chunk with no cross-chunk synchronization, so a
  two-million-region run uses the same memory per task as a three-thousand
  region run.
- **Scheduler-agnostic.** Local, SLURM, PBS/PBS Pro, SGE, and LSF are profiles,
  not code paths. `-resume` picks up where an interrupted run stopped.
- **Real ranking.** Hard filters remove risky guides; passing guides receive a
  weighted composite score and a rank within their target region.
- **Provenance by default.** Every output directory contains
  `run_manifest.json` and `software_versions.yml` with the pipeline version,
  commit, parameters, and tool versions.
- **Optional, not mandatory.** ATAC-seq annotation, launch-dCas9 CNN efficacy
  scoring, and promoter-proximal seed screening are opt-in and documented
  separately.

## At a glance

| | |
|---|---|
| **Input** | SNP points (`chr:pos`, 1-based) or enhancer intervals (`chr:start-end`), BED, one or both |
| **Reference** | hg38 + a FlashFry index; other genomes work with a matching index but are not validated |
| **Execution** | Local, SLURM, PBS/PBS Pro, SGE, LSF; conda, mamba, micromamba, Docker, Singularity, or Apptainer |
| **Output** | `guides.ranked.tsv`, optional `guides.final.tsv`, a seed hit table, and provenance |
| **Scale** | From a handful of targets to a genome-wide cCRE library |
| **Measured** | 2.35 M cCREs → 92.5 M guides in ~24 h on one 20-core node |

## Install

The Python package provides the step-level CLI. The Nextflow pipeline provides
the supported end-to-end workflow.

```bash
# Python package (PyPI once the release is published; see Installation for
# the tagged-source fallback)
pip install cretarget

# External requirements for the full pipeline
nextflow -version    # must be >= 24.04
java -version        # must be >= 17
bedtools --version
samtools --version
```

See [Installation](getting-started/installation.md) for containers, conda
profiles, and the exact version floor.

## Quick start

```bash
# 1. Fetch the hg38 FASTA + FlashFry index (~3.2 GB, login node)
cretarget setup-reference

# 2. Record the reference paths for the step-level CLI
cretarget configure --set reference.fasta=~/.cretarget/reference/hg38.fa
cretarget configure --set reference.flashfry_jar=/path/to/FlashFry-assembly-1.15.jar
cretarget configure --set reference.flashfry_db=~/.cretarget/reference/hg38_standard_cas9ngg_db

# 3. Run the pipeline
nextflow run MikeKazz/creTarget \
    -profile local,conda \
    --enhancer_bed enhancers.txt \
    --snp_bed snps.txt \
    --fasta ~/.cretarget/reference/hg38.fa \
    --flashfry_jar /path/to/FlashFry-assembly-1.15.jar \
    --flashfry_db ~/.cretarget/reference/hg38_standard_cas9ngg_db \
    --outdir results
```

At least one of `--enhancer_bed` or `--snp_bed` is required. The result is
`results/guides.ranked.tsv`, one row per candidate guide with FlashFry scores,
genome coordinates, filter status, and a per-target rank.

The full walkthrough is in [Quick start](getting-started/quickstart.md).

## Documentation map

<div class="grid cards" markdown>

- :material-rocket-launch: **New to creTarget?**

    Start with [Installation](getting-started/installation.md), then follow the
    [Quick start](getting-started/quickstart.md).

- :material-format-list-bulleted: **Inputs and outputs**

    [Input formats](guides/inputs.md), [output files](guides/outputs.md), and
    the complete [output column reference](reference/output-columns.md).

- :material-tune: **Parameters**

    Every pipeline parameter, default, and profile is in
    [Pipeline parameters](reference/parameters.md). The step-level commands are
    in [Step-level CLI](reference/cli.md).

- :material-dna: **Optional methods**

    [launch-dCas9 scoring](guides/launch-dcas9.md),
    [ATAC-seq annotation](guides/atac-seq.md), and
    [promoter-proximal seed screening](guides/seed-screening.md).

- :material-server-network: **Clusters and scale**

    [Scaling and clusters](guides/scaling-and-clusters.md), the
    [cluster walkthrough](cluster-setup.md), and the
    [architecture](project/architecture.md).

- :material-help-circle: **Something is wrong**

    [Troubleshooting](guides/troubleshooting.md) and the [FAQ](guides/faq.md).
    If neither helps, open an issue with the template.

</div>

## Citation

If you use creTarget, cite the software through its Zenodo concept DOI,
[10.5281/zenodo.24148024](https://doi.org/10.5281/zenodo.24148024), and cite the
methods you use. See [Citation](project/citation.md) for BibTeX and the
required citations for FlashFry, launch-dCas9, and seed screening.

## License

MIT. See [LICENSE](https://github.com/MikeKazz/creTarget/blob/main/LICENSE).
