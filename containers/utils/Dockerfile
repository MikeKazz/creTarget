# Lightweight image used for PREP_REGIONS, ANNOTATE, and RANK_AND_REPORT.
# Just bedtools + python3 + pandas -- no Bioconductor/R needed since
# PAM-finding is delegated entirely to FlashFry.

FROM python:3.11-slim

RUN apt-get update && apt-get install -y --no-install-recommends bedtools samtools && \
    rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir pandas

WORKDIR /work
