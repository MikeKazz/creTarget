# Installation

creTarget has two parts:

1. **The Python package** (`cretarget`) — the step-level engine and CLI.
2. **The Nextflow pipeline** — the supported end-to-end interface, fetched from
   GitHub by Nextflow or run from a clone.

## Requirements

| Component | Minimum | Notes |
|---|---|---|
| [Nextflow](https://nextflow.io) | **24.04** | Hard floor: `process.resourceLimits` does not exist earlier. |
| Java | **17** | Nextflow and FlashFry need a modern JVM. |
| Python | **3.10** | Tested on 3.10–3.12. |
| [bedtools](https://bedtools.readthedocs.io/) | 2.27+ | Region extraction. |
| [samtools](https://www.htslib.org/) | 1.10+ | FASTA indexing. |
| [FlashFry](https://github.com/mckennalab/FlashFry/releases) | 1.15 | Separate JAR, not bundled. |
| hg38 + FlashFry index | — | ~3.2 GB, fetched by `cretarget setup-reference`. |

A GPU is **not** required. The core pipeline is CPU-only; a GPU only helps
`--run_launch_dcas9`.

## Python package

```bash
pip install cretarget
```

The release workflow is set up to publish tagged releases to PyPI. If the command reports
`No matching distribution found`, the release is not on PyPI yet; install the
tagged source archive:

```bash
pip install "https://github.com/MikeKazz/creTarget/archive/refs/tags/v0.5.1.tar.gz"
```

For development:

```bash
git clone https://github.com/MikeKazz/creTarget.git
cd creTarget
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -e ".[dev]"
```

## Container

```bash
nextflow run MikeKazz/creTarget -r v0.5.1 \
    -profile slurm,singularity \
    --container_image ghcr.io/mikekazz/cretarget:0.5.1 \
    ...
```

The image contains Python, the package, bedtools, samtools, and a JRE. It does
**not** bundle FlashFry, the genome, or model checkpoints; pass them as
parameters.

> A container profile with no image is a hard error. Running an engine with no
> image would silently use host tools while reporting a container.

## Conda

> There is no `cretarget` Bioconda package yet. The `conda`, `mamba`, and
> `micromamba` profiles build the **task environment** from
> `assets/environment.yml`; they do not install a package named `cretarget`
> from a channel.

```bash
nextflow run MikeKazz/creTarget -profile slurm,conda ...
```

Use `-profile micromamba` if the classic solver is slow.

## External tools

```bash
conda install -c conda-forge -c bioconda bedtools samtools
curl -LO https://github.com/mckennalab/FlashFry/releases/download/1.15/FlashFry-assembly-1.15.jar
curl -s https://get.nextflow.io | bash
```

## Verify

```bash
cretarget --version
cretarget --help
cretarget doctor
```

`doctor` checks external tools, configured reference paths, the `.fai` and
FlashFry header, cached models, and the seed reference. It exits nonzero when a
required piece is missing and reports optional pieces without failing.

## Cache

| Setting | Resolution order |
|---|---|
| `CRETARGET_CACHE_DIR` | Highest priority |
| `[cache].dir` in `config.toml` | Used when the environment variable is unset |
| `~/.cretarget` | Default |

The step-level CLI reads `~/.config/cretarget/config.toml` (override the
directory with `CRETARGET_CONFIG_DIR`). Nextflow does **not** read that file;
pass pipeline parameters explicitly.

On a cluster whose compute nodes have no network, set the cache to shared
storage and run every fetch command on a login node:

```bash
export CRETARGET_CACHE_DIR=/shared/path/.cretarget
cretarget setup-reference
cretarget fetch-models          # only for launch-dCas9
cretarget fetch-atac --cell-type k562
```

Pipeline tasks only read from the cache and never download.

## Next steps

- [Quick start](Quick-Start) — a complete run.
- [Cluster setup](Cluster-Setup) — first install on a scheduler.
- [Full installation guide](https://mikekazz.github.io/creTarget/getting-started/installation/)
