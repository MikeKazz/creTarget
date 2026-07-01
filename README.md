# crispri-guide-designer

Design and rank CRISPRi guide RNAs for SNP/point or enhancer/interval target
regions (hg38), using FlashFry for guide discovery and off-target/on-target
scoring.

## Status

**Working MVP via `run_standalone.sh`, validated against real hg38 + real
FlashFry output on the cluster.** The Nextflow pipeline (`main.nf`) wires the
same underlying logic but hasn't been run for real yet (needs either the
hosted reference bundle + Docker images to be published, or `-profile bare`
against your existing local install). See "Nextflow (not yet run for real)"
below.

## How to run this, every time (the short version)

1. Download the zip from the chat, upload/replace it in your sandbox directory.
2. Unzip it (this overwrites the code, but **does not** touch anything you
   put outside this folder, like your own BED files):
   ```bash
   cd /home/jmorrislab1/mkazemi/___ProJECTS___/Blair/gRNA_design_tool/sandbox
   unzip -o crispri-guide-designer.zip
   cd crispri-guide-designer
   ```
3. Run it:
   ```bash
   ./run_standalone.sh --snp_bed /path/to/snps.bed --enhancer_bed /path/to/enhancers.bed --dedup_blast_radius
   ```

That's it — you do **not** need to pass `--fasta`, `--flashfry_jar`, or
`--flashfry_db` every time. Those are read automatically from
`local_paths.env` (sitting next to the script), which is pre-filled with your
real cluster paths:
```
FASTA=/home/jmorrislab1/mkazemi/bin/hg38.fa
FLASHFRY_JAR=/home/jmorrislab1/mkazemi/bin/FlashFry-assembly-1.15.jar
FLASHFRY_DB=/home/jmorrislab1/mkazemi/bin/hg38_standard_cas9ngg_db
```
If any of those paths ever change, edit `local_paths.env` (or pass
`--fasta`/`--flashfry_jar`/`--flashfry_db` on the command line to override
just for one run).

**The script also doesn't care what your current directory is** when you run
it — all of its own file lookups (the `bin/` scripts, `local_paths.env`,
default output location `results/`) are resolved relative to *where the
script itself lives*, not your working directory. `--snp_bed`/
`--enhancer_bed`, on the other hand, are resolved relative to your current
directory (or you can just use absolute paths) — so your BED files can live
anywhere, no need to copy them into this folder.

Output lands in `results/guides.ranked.tsv` inside this folder by default
(override with `--outdir`).

Only `--bed` is required. On first run, the pipeline downloads a pre-built hg38
reference FASTA + FlashFry index (main chromosomes only, no alt scaffolds, no chrM)
into `~/.crispri-guide-designer/` and reuses it on every subsequent run.

See `main.nf --help` for the full parameter list.

## Design decisions (from planning discussion)

### Nextflow (not yet run for real)

```bash
# Testing against your existing local FlashFry install (no containers needed):
nextflow run main.nf \
    --snp_bed snps.bed --enhancer_bed enhancers.bed \
    --local_fasta /home/jmorrislab1/mkazemi/bin/hg38.fa \
    --local_flashfry_jar /home/jmorrislab1/mkazemi/bin/FlashFry-assembly-1.15.jar \
    --local_flashfry_db /home/jmorrislab1/mkazemi/bin/hg38_standard_cas9ngg_db \
    --dedup_blast_radius \
    -profile bare

# Once the hosted reference bundle + Docker images are published:
nextflow run main.nf --snp_bed snps.bed --enhancer_bed enhancers.bed -profile standard,docker
```

- **PAM-finding delegated to FlashFry's `discover`** (not a custom R scanner) —
  simpler pipeline, one fewer custom code path.
- **Two input types, two BED files**: `--snp_bed` (point targets, ±150bp search
  window) and `--enhancer_bed` (interval targets, interval ±50bp flank). Either
  or both may be supplied in one run.
- **CRISPRi silencing radius (~1kb)** governs two things: (1) optional flagging
  of input targets within `blast_radius_bp` of each other (`--dedup_blast_radius`
  — flags as linked, does NOT merge/collapse them), and (2) post-hoc flagging of
  guides whose site falls outside that radius from the intended target
  (SNP position, or enhancer midpoint).
- **No `--transformPositions`**: genome coordinates are recovered by parsing the
  FASTA header (`<region_id>__<chrom>:<padded_start>-<padded_end>`, written by
  `prep_regions.py`) out of FlashFry's `contig` output column — this mirrors the
  existing `id_<n>_chr:start-end_<seq>` convention already used in production,
  adapted to region-level (not per-guide) FASTA input.
- **Scoring metrics and column names verified against real FlashFry output**
  (not just docs): `Doench2014OnTarget, DoenchCFD_maxOT,
  DoenchCFD_specificityscore, dangerous_GC, dangerous_polyT,
  dangerous_in_genome, Hsu2013, Moreno-Mateos2015OnTarget, JostCRISPRi_maxOT,
  JostCRISPRi_specificityscore, basesDiffToClosestHit, closestHitCount,
  0-1-2-3-4_mismatch, FreeEnergy20, otCount`. `JostCRISPRi_specificityscore` is
  the CRISPRi-specific off-target metric and current provisional ranking key.
- **`--maxMismatch 3`** for both `discover` and `score` (FlashFry default).
- Real sample output showed `Doench2014OnTarget = NA` when guides were scored
  with only a fixed 10bp synthetic flank (the old script's `FlashFrySEQ`
  padding). This pipeline extracts real genomic flanking sequence per region
  (default 50bp for enhancers, 150bp for SNPs) instead, which should resolve it
  — not yet confirmed against a real run.
- **Reference/index**: hosted as a versioned release bundle (see
  `assets/build_reference_index.md`), auto-downloaded and cached locally.
- **Execution**: Nextflow's native `slurm` executor replaces manual `sbatch`
  array-job scripting; `-resume` replaces manual skip-if-exists chunk tracking.
- **Containers**: one lightweight Python+bedtools image for all non-FlashFry
  steps (PAM-finding is fully delegated to FlashFry, so no Bioconductor/R stack
  is needed), one image with the FlashFry JAR baked in.

## Repo layout

```
main.nf                        Top-level workflow
nextflow.config                Default params + profile includes
conf/
  base.config                  Per-process resource defaults
  standard.config               Local execution profile
  slurm.config                  SLURM execution profile (edit queue/account for your cluster)
modules/local/
  setup_reference.nf            Cache-aware reference + FlashFry index fetch
  prep_regions.nf                BED parsing, optional blast-radius dedup, FASTA extraction
  flashfry_discover.nf           FlashFry `discover`
  flashfry_score.nf              FlashFry `score`
  annotate.nf                    Merge + map guides back to regions, distance/window flagging
  rank_report.nf                 PLACEHOLDER — ranking logic TBD
containers/
  flashfry/Dockerfile            FlashFry JAR baked in
  r-env/Dockerfile                Bioconductor R environment
assets/
  build_reference_index.md       Exact recipe for building/hosting the reference bundle
```

## Open items before this is runnable end-to-end

1. Host the reference/index bundle (Zenodo / GitHub Releases / S3) and set
   `params.reference_url_base` + implement the actual download logic in
   `SETUP_REFERENCE`.
2. Build and publish the two Docker images; replace `TODO/...` container refs
   in each module with real image URIs.
3. Implement `PREP_REGIONS` (port blast-radius dedup + windowing logic from
   `230825_design-gRNAs.R`, generalized to arbitrary BED input).
4. Implement `ANNOTATE` (merge scored chunks, join to region metadata, distance
   calc, window flagging).
5. Define and implement ranking logic in `RANK_AND_REPORT`.
6. Repeat/poly-T/duplicate-spacer hard filters from the original script still
   need a home — likely in `PREP_REGIONS` or as a new step between `FLASHFRY_SCORE`
   and `ANNOTATE`.
