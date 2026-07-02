# crispri-guide-designer

Design, score, and rank CRISPRi guide RNAs for noncoding target regions
(SNPs/points or enhancers/intervals) in hg38 — fast, at scale (validated at
~464,000 guides in ~6 minutes), with real off-target/on-target scoring via
[FlashFry](https://github.com/mckennalab/FlashFry) and optional CNN-based
enhancer/promoter efficacy scoring via
[launch-dCas9](https://github.com/Wancen/launch-dCas9).

## Features

- **Two input types, flexible formats**: SNPs (`chr:pos`, 1-based, or legacy
  BED) and enhancers (`chr:start-end`, or legacy BED) — auto-detected per
  line, malformed lines skipped with a warning rather than aborting the run.
- **Fast, memory-safe parallel execution**: FlashFry discover/score is
  chunked and run in parallel, either locally (auto-sized to your CPU count
  and available RAM) or as a SLURM array job.
- **Real ranking**: hard filters (seed score, GC content, poly-T, multi-site
  genomic match, FlashFry overflow) exclude bad guides; passing guides get a
  weighted composite score (CRISPRi/general off-target specificity, CFD
  specificity, closest off-target distance, RNA folding energy, and
  optionally launch-dCas9 on-target efficacy) and are ranked per target.
- **Optional ATAC-seq integration**: pass a narrowPeak file to annotate each
  guide with its distance to the nearest accessible-chromatin peak summit —
  one of the strongest predictors of CRISPRi efficacy in noncoding regions.
- **Zero manual reference setup**: the hg38 FASTA and FlashFry index are
  hosted on Zenodo and fetched with one command.

## Installation

Requires Python ≥3.10, [FlashFry](https://github.com/mckennalab/FlashFry)
(a JVM tool — needs Java), `bedtools`, and `samtools`.

```bash
git clone git@github.com:MikeKazz/GuideNet.git
cd GuideNet
pip install -e .
```

(Once published to PyPI: `pip install crispri-guide-designer`.)

## First-time setup

```bash
# 1. Download the hg38 reference + FlashFry index (~3.2GB total)
crispri-guide-designer setup-reference --outdir ~/.guidenet_reference

# 2. Point the tool at it + your FlashFry install
crispri-guide-designer configure --set reference.fasta=~/.guidenet_reference/hg38.fa
crispri-guide-designer configure --set reference.flashfry_jar=/path/to/FlashFry-assembly-1.15.jar
crispri-guide-designer configure --set reference.flashfry_db=~/.guidenet_reference/hg38_standard_cas9ngg_db
```

Config lives at `~/.config/crispri-guide-designer/config.toml` — set once,
used automatically on every run after. FlashFry itself isn't hosted by this
tool; download it separately from its
[GitHub releases](https://github.com/mckennalab/FlashFry/releases).

### Optional: launch-dCas9 (enhancer/promoter on-target scoring)

Requires a separate conda environment with `torch` (see
[launch-dCas9](https://github.com/Wancen/launch-dCas9)). Checkpoints ship in
this repo (`checkpoints/`, ~1MB each) — no separate download.

```bash
crispri-guide-designer configure --set launch_dcas9.venv_python=/path/to/conda/envs/launch-dCas9/bin/python
crispri-guide-designer configure --set launch_dcas9.script=/path/to/launch-dCas9.py
```

### Optional: SLURM

If running `--executor slurm`, set your cluster's resource defaults once:

```bash
crispri-guide-designer configure --set slurm.partition=your_partition
crispri-guide-designer configure --set slurm.gres=gpu:1   # blank if no GPU needed
crispri-guide-designer configure --set slurm.mem=16G
crispri-guide-designer configure --set slurm.cpus=8
crispri-guide-designer configure --set slurm.time=2:00:00
```

## Usage

```bash
crispri-guide-designer run \
    --snp-bed snps.txt \
    --enhancer-bed enhancers.txt \
    --dedup-blast-radius \
    --outdir results
```

At least one of `--snp-bed` / `--enhancer-bed` is required; both may be
given together. Input files accept either plain text (`chr:pos` per line for
SNPs, `chr:start-end` per line for enhancers) or legacy tab-separated BED.

Output: `results/guides.ranked.tsv` — one row per candidate guide, annotated
with all FlashFry scores, genome coordinates, hard-filter results, and
`rank` (1 = best, per target; `exclude (<reason>)` if it failed a filter).

### Common options

| Flag | Default | Description |
|---|---|---|
| `--dedup-blast-radius` | off | Flag input targets within 1kb of each other (does not merge) |
| `--blast-radius-bp` | 1000 | CRISPRi silencing window, bp |
| `--flank-bp` | 50 | Enhancer-type: padding added around each interval |
| `--snp-window-bp` | 150 | SNP-type: search window on each side of the point |
| `--atac-narrowpeak` | none | ENCODE-standard narrowPeak file; adds `distance_to_atac_summit` column (informational, not yet in ranking) |
| `--run-launch-dcas9` | off | Add launch-dCas9 enhancer/promoter efficacy scoring |
| `--executor` | `local` | `local` or `slurm` |
| `--n-workers` | auto | Parallel chunks; auto-sized by CPU count and available RAM if unset |

Run `crispri-guide-designer run --help` for the full list.

## How it works

1. **Prep regions** — pad each target, extract sequence from hg38, encode
   genome coordinates into FASTA headers.
2. **FlashFry discover + score** — split into chunks, run in parallel
   (local processes or SLURM array).
3. **Annotate** — merge chunks, recover genome coordinates, flag guides
   within the CRISPRi blast radius, compute ATAC summit distance if given.
4. **launch-dCas9** *(optional)* — CNN-scored enhancer/promoter on-target
   efficacy, same chunked parallel execution.
5. **Rank** — hard filters exclude bad guides; passing guides get a
   weighted composite score and a per-target rank.

## Performance

Validated on a 40-core / 125GB node: 3,517 real enhancer regions (~464,000
candidate guides) processed end-to-end in ~6 minutes — about 13× faster
than unchunked, single-process execution, with parallelism auto-capped to
stay within available memory (each FlashFry chunk loads the full off-target
index, so more workers isn't always faster — see `_default_n_workers` in
`cli.py`).

## Repository layout

```
src/crispri_guide_designer/
  cli.py               crispri-guide-designer command (run / configure / setup-reference)
  config.py            ~/.config/crispri-guide-designer/config.toml management
  executor.py           local-parallel and SLURM-array execution
  flashfry_chunked.py   chunked FlashFry orchestration
  scripts/               region prep, annotation, ranking, launch-dCas9 scoring
checkpoints/            launch-dCas9 pretrained weights (bundled, ~1MB each)
test_data/               small real-coordinate fixtures for smoke testing
run_standalone.sh        bash-only alternative entry point (no pip install needed)
inspect_repo.sh          pre-commit sanity checks (large files, stray artifacts, etc.)
```

## Known limitations

- FlashFry (JVM) and launch-dCas9 (conda/torch) are external dependencies,
  not bundled — `pip install` gets you the orchestration tool, not those.
- `--executor slurm` is implemented but not yet run-verified end-to-end
  (validated extensively with `--executor local`).
- ATAC distance is computed and reported but not yet part of the ranking
  formula.

## License

MIT
