# Installation

creTarget has two parts:

1. **The Python package** (`cretarget`) — the step-level engine and CLI.
2. **The Nextflow pipeline** (`main.nf`, `workflows/`, `modules/`, `conf/`) —
   the supported end-to-end interface, fetched from GitHub by Nextflow or run
   from a clone.

You can install only the Python package for step-level work, but a full
pipeline run needs Nextflow and the external tools below.

## Requirements

| Component | Minimum | Notes |
|---|---|---|
| [Nextflow](https://nextflow.io) | **24.04** | A hard floor. `process.resourceLimits` does not exist earlier, so resource ceilings would silently not apply. |
| Java | **17** | Nextflow and FlashFry both need a modern JVM. Nextflow 23.x caps at Java 19, another reason to use a current release. |
| Python | **3.10** | For the Python package. Tested on 3.10–3.12. |
| [bedtools](https://bedtools.readthedocs.io/) | 2.27+ | Region extraction and promoter windows. |
| [samtools](https://www.htslib.org/) | 1.10+ | FASTA indexing. |
| [FlashFry](https://github.com/mckennalab/FlashFry/releases) | 1.15 | Download the assembly JAR separately. Not bundled. |
| hg38 FASTA + FlashFry index | — | ~3.2 GB. Fetched by `cretarget setup-reference`. |

A GPU is **not** required. The core pipeline is CPU-only; a GPU only helps the
optional `--run_launch_dcas9` step.

## Install the Python package

### From PyPI

```bash
pip install cretarget
```

The release workflow is set up to publish tagged releases to PyPI. If this command reports
`No matching distribution found`, the release is not on PyPI yet; use the
tagged source install below.

### From a tagged source archive

```bash
pip install "https://github.com/MikeKazz/creTarget/archive/refs/tags/v0.5.1.tar.gz"
```

This is the most reproducible install before a PyPI release exists. Replace the
tag with the version you want; the tag is immutable.

### From a clone, for development

```bash
git clone https://github.com/MikeKazz/creTarget.git
cd creTarget
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -e ".[dev]"
```

The editable install is what `CONTRIBUTING.md` assumes.

### Verify the CLI

```bash
cretarget --version
cretarget --help
cretarget doctor
```

`doctor` checks external tools, configured paths, and cached data. It fails on
missing required pieces and reports optional ones without failing the run.

## Container (most portable pipeline)

The published image is `ghcr.io/mikekazz/cretarget:0.5.1`. It contains Python,
the Python package, bedtools, samtools, and a JRE. It deliberately does **not**
bundle FlashFry, the genome, or model checkpoints; those are staged in by
Nextflow from paths you provide.

```bash
nextflow run MikeKazz/creTarget \
    -profile slurm,singularity \
    --container_image ghcr.io/mikekazz/cretarget:0.5.1 \
    ...
```

Use `-profile docker` where Docker is available. Use `singularity` or
`apptainer` on shared clusters. The pipeline refuses to start if a container
profile is enabled with no image, rather than silently running against host
tools.

## Conda environments

> **There is no `cretarget` Bioconda package yet.** The `conda`, `mamba`, and
> `micromamba` profiles build the *task environment* from
> [`assets/environment.yml`](https://github.com/MikeKazz/creTarget/blob/main/assets/environment.yml);
> they do not install a package called `cretarget` from a channel. Install the
> Python package with `pip`, or use the container.

```bash
# First run creates the environment; later runs reuse it
nextflow run MikeKazz/creTarget -profile slurm,conda ...
```

If conda's classic solver is slow, use `-profile micromamba` (same environment,
faster solver). The environment pins Python, pandas, numpy, bedtools, samtools,
and a JRE to the versions the test suite and benchmark used.

## Install the external tools

`bedtools` and `samtools` are available from conda-forge/bioconda or your
system package manager:

```bash
conda install -c conda-forge -c bioconda bedtools samtools
```

FlashFry is a separate download:

```bash
curl -LO https://github.com/mckennalab/FlashFry/releases/download/1.15/FlashFry-assembly-1.15.jar
```

Nextflow is a single launcher:

```bash
curl -s https://get.nextflow.io | bash
mv nextflow ~/bin/       # or anywhere on PATH
nextflow -version
```

## Cache and configuration

Reference data, ATAC peaks, and model checkpoints live in a cache directory:

| Setting | Resolution order |
|---|---|
| `CRETARGET_CACHE_DIR` | Highest priority |
| `[cache].dir` in `config.toml` | Used when the environment variable is unset |
| `~/.cretarget` | Default |

The step-level CLI reads `~/.config/cretarget/config.toml` (override the
directory with `CRETARGET_CONFIG_DIR`). Nextflow runs **do not** read that file;
pass pipeline parameters explicitly or use a params file.

On a cluster where compute nodes have no outbound network, set the cache to a
shared path and run every fetch command on a login node:

```bash
export CRETARGET_CACHE_DIR=/shared/path/.cretarget
cretarget setup-reference
cretarget fetch-models          # only if using launch-dCas9
cretarget fetch-atac --cell-type k562
```

Pipeline tasks only read from the cache and never download.

## Next steps

- [Quick start](quickstart.md) — a complete run from an empty directory.
- [Reference data](reference-data.md) — what `setup-reference` downloads and
  how to use your own genome.
- [Scaling and clusters](../guides/scaling-and-clusters.md) — running on a
  scheduler.
