# creTarget

[![CI](https://github.com/MikeKazz/creTarget/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/MikeKazz/creTarget/actions/workflows/ci.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24148024.svg)](https://doi.org/10.5281/zenodo.24148024)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Container](https://img.shields.io/badge/container-ghcr.io-blue)](https://github.com/MikeKazz/creTarget/pkgs/container/cretarget)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue)](pyproject.toml)
[![Nextflow ≥ 24.04](https://img.shields.io/badge/Nextflow-%E2%89%A524.04-23aa62)](https://www.nextflow.io/)
[![Docs](https://img.shields.io/badge/docs-mkdocs--material-526CFE)](https://mikekazz.github.io/creTarget/)

**Design, score, and rank CRISPRi guide RNAs for noncoding target regions —
SNPs, enhancers, or arbitrary intervals — in hg38, from a handful of targets to
genome-wide libraries.**

creTarget takes SNP points or enhancer intervals, extracts padded target
windows, discovers and scores candidate guides with
[FlashFry](https://github.com/mckennalab/FlashFry), optionally adds CNN-based
efficacy scoring from [launch-dCas9](https://github.com/Wancen/launch-dCas9),
ATAC-seq accessibility annotation, and promoter-proximal seed screening, then
returns one ranked table with provenance attached.

Execution is handled by [Nextflow](https://nextflow.io), so the same command
runs on a laptop, SLURM, PBS/PBS Pro, SGE, or LSF. The scheduler is a `-profile`
flag, not a code path.

## Why creTarget

- **Two input types, one workflow.** SNP points (`chr:pos`, 1-based) and
  enhancer intervals (`chr:start-end` or BED) can be mixed in one run.
  Malformed lines are skipped with a warning; contigs are validated against the
  reference before bedtools sees them.
- **Flat memory at any scale.** Only region preparation is global. Everything
  after it runs per chunk with no cross-chunk synchronization, so a
  two-million-region run uses the same memory per task as a three-thousand
  region run.
- **Real ranking, not just scoring.** Hard filters remove guides with seed,
  GC, poly-T, multi-site, or overflow problems. Passing guides get a weighted
  composite score and a rank within their target region.
- **Optional methods are explicit.** launch-dCas9, ATAC-seq, and seed screening
  are opt-in, documented, and cite their sources.
- **Provenance by default.** Every output directory contains
  `run_manifest.json` and `software_versions.yml` with the version, commit,
  parameters, command line, and tool versions.

## Quick start

```bash
# 1. Install the Python package
pip install cretarget

# 2. Fetch hg38 + the FlashFry index (~3.2 GB, login node)
cretarget setup-reference

# 3. Record the reference paths for the step-level CLI
cretarget configure --set reference.fasta=~/.cretarget/reference/hg38.fa
cretarget configure --set reference.flashfry_jar=/path/to/FlashFry-assembly-1.15.jar
cretarget configure --set reference.flashfry_db=~/.cretarget/reference/hg38_standard_cas9ngg_db

# 4. Run
nextflow run MikeKazz/creTarget -r v0.5.1 \
    -profile local,conda \
    --enhancer_bed enhancers.txt \
    --snp_bed snps.txt \
    --fasta ~/.cretarget/reference/hg38.fa \
    --flashfry_jar /path/to/FlashFry-assembly-1.15.jar \
    --flashfry_db ~/.cretarget/reference/hg38_standard_cas9ngg_db \
    --outdir results
```

At least one of `--snp_bed` / `--enhancer_bed` is required. The main result is
`results/guides.ranked.tsv`: one row per candidate guide, with FlashFry scores,
genome coordinates, filter status, and a per-target `rank`.

> `pip install cretarget` requires the release to be on PyPI. If it is not
> there yet, install the tagged source archive:
> `pip install "https://github.com/MikeKazz/creTarget/archive/refs/tags/v0.5.1.tar.gz"`.

The full walkthrough, including input examples and output interpretation, is in
[docs/getting-started/quickstart.md](docs/getting-started/quickstart.md).

## Requirements

| Component | Minimum | Notes |
|---|---|---|
| [Nextflow](https://nextflow.io) | **24.04** | Hard floor: `process.resourceLimits` does not exist earlier, so resource ceilings would silently not apply. |
| Java | **17** | Nextflow and FlashFry both need a modern JVM. |
| Python | **3.10** | For the Python package. Tested on 3.10–3.12. |
| [bedtools](https://bedtools.readthedocs.io/) | 2.27+ | Region extraction and promoter windows. |
| [samtools](https://www.htslib.org/) | 1.10+ | FASTA indexing. |
| [FlashFry](https://github.com/mckennalab/FlashFry/releases) | 1.15 | Download the assembly JAR separately. |
| hg38 + FlashFry index | — | ~3.2 GB, fetched by `cretarget setup-reference`. |

A GPU is not required. The core pipeline is CPU-only; a GPU only helps
`--run_launch_dcas9`.

## Installation

### Python package

```bash
pip install cretarget
```

The release workflow is set up to publish tagged releases to PyPI. If that command reports
no matching distribution, use the tagged source archive above. See
[docs/getting-started/installation.md](docs/getting-started/installation.md) for
editable installs and version details.

### Container

The published image is `ghcr.io/mikekazz/cretarget:0.5.1`. It contains Python,
the package, bedtools, samtools, and a JRE. FlashFry, the genome, and the model
checkpoints are deliberately not bundled; pass them as parameters.

```bash
nextflow run MikeKazz/creTarget -r v0.5.1 \
    -profile slurm,singularity \
    --container_image ghcr.io/mikekazz/cretarget:0.5.1 \
    ...
```

### Conda

> There is no `cretarget` Bioconda package yet. The `conda`, `mamba`, and
> `micromamba` profiles build the *task environment* from
> [`assets/environment.yml`](assets/environment.yml); they do not install a
> package named `cretarget` from a channel.

```bash
nextflow run MikeKazz/creTarget -profile slurm,conda ...
```

Use `-profile micromamba` if the classic solver is slow.

### Check the installation

```bash
cretarget doctor
```

It checks external tools, configured reference paths, the `.fai` and FlashFry
header, cached models, and the seed reference.

## Output

```text
results/
├── guides.ranked.tsv
├── guides.final.tsv                  # with --run_seed_finder
├── run_manifest.json
├── software_versions.yml
├── regions/regions.metadata.tsv
├── seed_finder/<prefix>_raw_hits.csv # with --run_seed_finder
└── pipeline_info/
```

`rank` is an integer within each `region_id` (`1` is best) or
`exclude (<reason>)` for a guide that failed a hard filter. Excluded guides stay
in the table: a filtered guide is evidence, and removing rows makes "no guide
found" indistinguishable from "guides filtered out".

The complete column list, ranking weights, and filter rules are in
[docs/guides/outputs.md](docs/guides/outputs.md) and
[docs/reference/output-columns.md](docs/reference/output-columns.md).

## Optional features

| Feature | Flag | What it adds |
|---|---|---|
| [launch-dCas9](docs/guides/launch-dcas9.md) | `--run_launch_dcas9` | CNN enhancer/promoter efficacy scores; enhancer score contributes 30% to ranking |
| [ATAC-seq](docs/guides/atac-seq.md) | `--atac_narrowpeak` | `distance_to_atac_summit`; annotation only, not ranked |
| [Seed screening](docs/guides/seed-screening.md) | `--run_seed_finder` | Promoter-proximal off-target repression screen and a hit table |

`--seed_finder_min_length` and `--seed_finder_tss_window` have no defaults on
purpose. They change what the screen means, so choose them deliberately.

## How it works

1. **Prep regions** — parse and validate targets, flag linked neighbours, pad
   each region, extract sequence, encode genome coordinates into FASTA headers.
   *(global, once)*
2. **Split** — divide regions into chunks of `--chunk_size`.
3. **FlashFry discover + score** — per chunk; JVM heap derived from the task's
   actual memory allocation.
4. **Annotate** — recover genome coordinates, join region metadata, flag guides
   within the blast radius, compute ATAC summit distance. *(per chunk)*
5. **launch-dCas9** *(optional)* — CNN on-target efficacy. *(per chunk)*
6. **Rank** — hard filters, weighted composite score, per-target rank.
   *(per chunk)*
7. **Collate** — stream per-chunk tables into `guides.ranked.tsv`.
8. **Seed screening** *(optional)* — promoter-proximal seed matches over the
   collated table, producing `guides.final.tsv`. *(once, not per chunk)*

Per-chunk ranking is exactly equivalent to global ranking because every guide
for a region lands in the same chunk. CI verifies that two independent runs with
different work directories produce byte-identical output.

## Scaling

| | |
|---|---|
| Execution | `local`, `slurm`, `pbs`, `sge`, `lsf` |
| Software | `conda`, `mamba`, `micromamba`, `docker`, `singularity`, `apptainer` |
| Chunk size | `--chunk_size` regions per task (default 5,000) |
| Resume | `-resume` |
| Resource ceilings | `--max_cpus`, `--max_memory`, `--max_time` |
| Retries | Automatic, with escalated memory and time |

Measured on a single node under SLURM (20 physical cores, 125 GB RAM, NFS):

| | |
|---|---|
| Input | 2,348,854 cCREs |
| Wall time | ~24 h |
| Guides | 92,519,254 |
| Regions with ≥1 usable guide | 99.52% |
| Output table | 50 GB (~10 GB gzipped) |
| Scheduler jobs | 7,049, zero failures or retries |

See [docs/guides/scaling-and-clusters.md](docs/guides/scaling-and-clusters.md)
for the full table, storage estimates, and site-profile guidance.

## Documentation

| | |
|---|---|
| [Documentation site](https://mikekazz.github.io/creTarget/) | The rendered MkDocs site |
| [Installation](docs/getting-started/installation.md) | PyPI, source, container, conda |
| [Quick start](docs/getting-started/quickstart.md) | End-to-end run |
| [Reference data](docs/getting-started/reference-data.md) | Genome, FlashFry index, model checkpoints, TSS reference |
| [Input formats](docs/guides/inputs.md) | SNP, enhancer, BED, coordinate conventions |
| [Output files](docs/guides/outputs.md) | Tree, rank, filters, ranking weights |
| [Parameters](docs/reference/parameters.md) | Every pipeline parameter |
| [Step-level CLI](docs/reference/cli.md) | `cretarget <step>` |
| [Cluster setup](docs/cluster-setup.md) | First install on a new cluster |
| [Troubleshooting](docs/guides/troubleshooting.md) | Errors with known causes |
| [FAQ](docs/guides/faq.md) | Short answers |
| [Architecture](docs/project/architecture.md) | Design decisions and repository layout |
| [Contributing](CONTRIBUTING.md) | Development setup and PR process |
| [Releasing](docs/project/releasing.md) | Maintainer release checklist |
| [Wiki](https://github.com/MikeKazz/creTarget/wiki) | Concise guides mirrored from the docs |

## Known limitations

- FlashFry (JVM) and launch-dCas9 (conda/torch) are external dependencies, not
  bundled.
- ATAC distance is computed and reported but is not part of the ranking formula.
- hg38 is the tested reference. Other genomes work with a matching FlashFry
  index, but are not validated.
- The composite score is a ranking heuristic, not a calibrated predictor of
  editing or repression efficiency.

## Citation

creTarget is archived on Zenodo. Cite
[10.5281/zenodo.24148024](https://doi.org/10.5281/zenodo.24148024) for the
software in general (this always resolves to the latest release), or
[10.5281/zenodo.24148025](https://doi.org/10.5281/zenodo.24148025) for v0.2.1
specifically when exact reproducibility matters.

If you use creTarget, please cite it (see [`CITATION.cff`](CITATION.cff)) along
with [FlashFry](https://github.com/mckennalab/FlashFry), and — if you use
`--run_launch_dcas9` — [launch-dCas9](https://github.com/Wancen/launch-dCas9).
If you use `--run_seed_finder`, also cite Hartman et al. 2026
([doi:10.64898/2026.03.27.714658](https://doi.org/10.64898/2026.03.27.714658)).
The full citation guidance is in [docs/project/citation.md](docs/project/citation.md).

## Contributing

Contributions are welcome. Read [CONTRIBUTING.md](CONTRIBUTING.md) and the
[Code of Conduct](CODE_OF_CONDUCT.md). For security issues, follow
[SECURITY.md](SECURITY.md); do not open a public issue.

## License

MIT. See [LICENSE](LICENSE).
