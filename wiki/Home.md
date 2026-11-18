# creTarget

[![CI](https://github.com/MikeKazz/creTarget/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/MikeKazz/creTarget/actions/workflows/ci.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24148024.svg)](https://doi.org/10.5281/zenodo.24148024)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/MikeKazz/creTarget/blob/main/LICENSE)

**Design, score, and rank CRISPRi guide RNAs for noncoding target regions —
SNPs, enhancers, or arbitrary intervals — in hg38, from a handful of targets to
genome-wide libraries.**

creTarget is a Nextflow pipeline and Python toolkit. It parses SNP points or
enhancer intervals, extracts padded target windows, discovers and scores
candidate guides with [FlashFry](https://github.com/mckennalab/FlashFry),
optionally adds [launch-dCas9](https://github.com/Wancen/launch-dCas9) CNN
efficacy scores, ATAC-seq accessibility annotation, and promoter-proximal seed
screening, then returns one ranked table with provenance attached.

> **The full documentation is at
> [mikekazz.github.io/creTarget](https://mikekazz.github.io/creTarget/).**
> This wiki is a concise mirror for quick lookup.

## Install

```bash
pip install cretarget
```

If that reports no matching distribution, the release is not on PyPI yet. Use
the tagged source archive instead:

```bash
pip install "https://github.com/MikeKazz/creTarget/archive/refs/tags/v0.5.1.tar.gz"
```

The pipeline also needs Nextflow ≥ 24.04, Java ≥ 17, `bedtools`, `samtools`,
and a separate FlashFry JAR. See [Installation](Installation).

## Five-minute run

```bash
cretarget setup-reference

printf 'chr1:10000-10400\nchr1:50000-50300\n' > enhancers.txt

nextflow run MikeKazz/creTarget -r v0.5.1 \
    -profile local,conda \
    --enhancer_bed enhancers.txt \
    --fasta ~/.cretarget/reference/hg38.fa \
    --flashfry_jar /path/to/FlashFry-assembly-1.15.jar \
    --flashfry_db ~/.cretarget/reference/hg38_standard_cas9ngg_db \
    --outdir results
```

The result is `results/guides.ranked.tsv`: one row per candidate guide, with
FlashFry scores, genome coordinates, filter status, and a per-target `rank`.
See [Quick start](Quick-Start).

## Documentation map

| Page | What it covers |
|---|---|
| [Installation](Installation) | PyPI/source/container/conda, requirements, cache |
| [Quick start](Quick-Start) | A complete run from an empty directory |
| [Inputs](Inputs) | SNP, enhancer, BED, coordinate conventions |
| [Outputs](Outputs) | Output tree, ranks, filters, ranking weights |
| [Parameters](Parameters) | Every pipeline parameter and profile |
| [Seed screening](Seed-Screening) | Promoter-proximal off-target screen |
| [Cluster setup](Cluster-Setup) | First install and first run on a scheduler |
| [Troubleshooting](Troubleshooting) | Errors with known causes |
| [FAQ](FAQ) | Short answers |
| [Citation](Citation) | DOIs, BibTeX, dependency citations |
| [Contributing](Contributing) | Development setup and PR process |

## Optional features

| Feature | Flag | What it adds |
|---|---|---|
| launch-dCas9 | `--run_launch_dcas9` | CNN enhancer/promoter efficacy; the enhancer score contributes 30% to ranking |
| ATAC-seq | `--atac_narrowpeak` | `distance_to_atac_summit`; annotation only, not ranked |
| Seed screening | `--run_seed_finder` | Promoter-proximal off-target repression screen and a hit table |

`--seed_finder_min_length` and `--seed_finder_tss_window` have no defaults on
purpose: they change what the screen means, so choose them deliberately.

## Citation

Cite the software through the Zenodo concept DOI
[10.5281/zenodo.24148024](https://doi.org/10.5281/zenodo.24148024), and cite
FlashFry, launch-dCas9, or Hartman et al. 2026 when you use those methods. See
[Citation](Citation).

## License

MIT. See [LICENSE](https://github.com/MikeKazz/creTarget/blob/main/LICENSE).
