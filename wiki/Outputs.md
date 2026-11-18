# Outputs

Everything is written below `--outdir` (default `results`). Per-chunk
intermediates stay in the Nextflow work directory.

```text
results/
├── guides.ranked.tsv                 # one row per candidate guide
├── guides.final.tsv                  # with --run_seed_finder
├── run_manifest.json
├── software_versions.yml
├── regions/
│   └── regions.metadata.tsv          # one row per input target
├── seed_finder/
│   └── <prefix>_raw_hits.csv         # with --run_seed_finder
└── pipeline_info/
    ├── report.html
    ├── timeline.html
    ├── trace.txt
    └── dag.html
```

`publish_dir_mode` defaults to `copy`. On a shared filesystem, `symlink` or
`rellink` can avoid duplicating large tables.

## `guides.ranked.tsv`

One row per candidate guide, tab-separated. Rows are sorted by `region_id`,
then passing guides before excluded guides, then by `rank` ascending. Within a
target, `rank = 1` is the best guide by the composite score.

Columns fall into five groups:

| Group | Examples |
|---|---|
| Guide and location | `contig`, `start`, `stop`, `target`, `spacer`, `orientation`, `genome_start`, `genome_end` |
| Target region | `region_id`, `chrom`, `orig_start`, `orig_end`, `input_type`, `target_ref`, `padded_start`, `padded_end`, `linked_regions` |
| FlashFry scores | `Doench2014OnTarget`, `DoenchCFD_specificityscore`, `Hsu2013`, `JostCRISPRi_specificityscore`, `basesDiffToClosestHit`, `FreeEnergy20`, `overflow`, `dangerous_*` |
| Annotation | `distance_to_target`, `within_blast_radius`, `distance_to_atac_summit` |
| Ranking | `seed_score`, `fail_*`, `exclude_reason`, `passes_hard_filters`, `composite_score`, `composite_score_max_possible`, `composite_score_normalized`, `rank` |

## Rank and hard filters

`rank` is an integer within each `region_id` for guides that pass the hard
filters. Excluded guides keep their row, but `rank` is `exclude (<reason>)` and
`passes_hard_filters` is `False`. A filtered guide is evidence; removing rows
makes "no guide found" indistinguishable from "guides filtered out".

| Reason | Filter |
|---|---|
| `seed_score` | Three or more overlapping `GG` dinucleotides in the 12 bp PAM-proximal seed |
| `dangerous_GC` | FlashFry flagged GC content |
| `dangerous_polyT` | FlashFry flagged a poly-T / Pol III terminator |
| `multisite_genomic_match` | More than one perfect match in the genome |
| `overflow` | FlashFry's off-target search overflowed |

The first failing check is reported, in that order.

## Composite score

The weights are literal percentages and are **not renormalized** when
launch-dCas9 is absent, so the score keeps its absolute meaning across runs:

| Component | Weight | Notes |
|---|---:|---|
| `launch_dCas9_enhancer` | 30% | Only when `--run_launch_dcas9` added the column |
| `JostCRISPRi_specificityscore` | 25% | |
| `Hsu2013` / 100 | 18% | Clipped to 0–1 |
| `DoenchCFD_specificityscore` | 8% | |
| `basesDiffToClosestHit` / `--max_mismatch`, capped at 1 | 5% | No off-target → 1.0 |
| `FreeEnergy20`, min–max normalized per region | 2% | Missing → 0.5 |

`composite_score_max_possible` is 0.58 without launch-dCas9 and 0.88 with it.
`composite_score_normalized` divides by that value. `launch_dCas9_promoter` is
reported but is **not** part of the ranking.

This is a ranking heuristic, not a calibrated probability. It orders guides
within a target; it does not promise that rank 1 works and rank 2 does not.

## `guides.final.tsv`

Written only with `--run_seed_finder`. Same rows and order as
`guides.ranked.tsv`, plus:

| Column | Meaning |
|---|---|
| `seedfinder_hit_count` | Non-self promoter-proximal seed hits; blank if not screened |
| `seedfinder_closest_gene` | Gene at the closest hit |
| `seedfinder_closest_tss_distance` | Signed distance to that TSS (negative upstream) |
| `seedfinder_closest_seed_length` | Seed length of the closest hit |
| `seedfinder_target_is_promoter` | Whether the guide's own target overlaps a promoter window |

A blank count means "not screened" (the guide failed a hard filter); `0` means
"screened, no hits". `--seed_finder_all_guides` screens everything.

## Seed hit table

`seed_finder/<prefix>_raw_hits.csv` has one row per seed match, including
`self_hit`, guide and hit coordinates, strand, seed length, nearest TSS and
signed distance, every TSS in the window, the PAM, the protospacer, and a UCSC
link. Self-hits are marked and excluded from `seedfinder_hit_count`.

## Provenance

`run_manifest.json` records the pipeline version, commit, dirty-tree flag,
command line, profiles, Nextflow version, work directory, user, and every
non-null parameter. `software_versions.yml` records the versions reported by
each process. Keep both with any published table.

## Interpreting a run

- **Filter before ranking.** `rank` only exists for passing guides.
- **Check `passes_hard_filters` before sorting by `composite_score`.**
- **Keep `region_id`.** It is the unit of ranking and the link to the target.
- **Look at the distribution, not one row.** A region with 40 guides and one
  with 1 are not equally well covered.
- **Do not compare scores across versions blindly.** Weights and filters are
  part of the method.

## Limitations

- FlashFry and launch-dCas9 are external dependencies, not bundled.
- `distance_to_atac_summit` is reported but not ranked.
- hg38 is the tested reference.
- The composite score is a ranking heuristic, not a validated predictor.

Full reference: [Output columns](https://mikekazz.github.io/creTarget/reference/output-columns/)
