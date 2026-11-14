# creTarget runtime image.
#
# Build:
#   docker build -t ghcr.io/mikekazz/cretarget:0.2.6 .
#
# Publish:
#   docker push ghcr.io/mikekazz/cretarget:0.2.6
#
# Use:
#   nextflow run MikeKazz/creTarget -profile slurm,singularity \
#       --container_image ghcr.io/mikekazz/cretarget:0.2.6 ...
#
# This is the most portable way to run the pipeline: one image, no conda solve,
# no version skew between users. It deliberately does NOT bundle FlashFry, the
# reference genome, or the model checkpoints -- those are passed as pipeline
# parameters and staged in by Nextflow, so bundling them would multiply the
# image size for no benefit and would redistribute third-party artifacts.
#
# micromamba base: small, and its solver handles the bioconda dependency set in
# seconds rather than the many minutes conda's classic solver can take.
FROM mambaorg/micromamba:2.10.0

LABEL org.opencontainers.image.title="creTarget" \
      org.opencontainers.image.description="CRISPRi guide RNA design for the noncoding genome" \
      org.opencontainers.image.source="https://github.com/MikeKazz/creTarget" \
      org.opencontainers.image.licenses="MIT" \
      org.opencontainers.image.version="0.2.6"

USER root
# git is needed to pip-install creTarget from the repository; procps provides
# `ps`, which Nextflow uses to collect per-task resource metrics. Without it
# every task's peak_rss comes back empty.
RUN apt-get update && apt-get install -y --no-install-recommends \
        git procps ca-certificates \
    && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

# Dependencies first, so editing the source does not invalidate the solve layer.
COPY --chown=$MAMBA_USER:$MAMBA_USER assets/environment.yml /tmp/environment.yml

# The pip entry in environment.yml installs creTarget from GitHub, which would
# bake in whatever `main` happened to be at build time. Strip it and install
# from the local source tree instead, so the image always matches the checkout
# it was built from.
RUN sed -i '/^  - pip:/,$d' /tmp/environment.yml \
    && micromamba install -y -n base -f /tmp/environment.yml \
    && micromamba clean --all --yes

COPY --chown=$MAMBA_USER:$MAMBA_USER . /opt/cretarget
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN pip install --no-cache-dir /opt/cretarget

ENV PATH="/opt/conda/bin:$PATH" \
    PYTHONDONTWRITEBYTECODE=1

# Fail the build rather than shipping a broken image.
RUN cretarget --version \
    && bedtools --version \
    && samtools --version | head -1 \
    && java -version 2>&1 | head -1 \
    && python -c "import pandas, numpy; print('pandas', pandas.__version__, 'numpy', numpy.__version__)"
