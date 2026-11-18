# Output files

A run writes everything below `--outdir` (default `results`). Per-chunk
intermediates stay in the Nextflow work directory; the output directory contains
only the files a user is expected to read or archive.

```text
results/
├── guides.ranked.tsv                 # one row per candidate guide
├── guides.final.tsv                  # optional: ranked + seed-screening columns
├── run_manifest.json                 # code, parameters, command line
├── software_versions.yml             # versions reported by each process
├── regions/
│   └── regions.metadata.tsv          # one row per input target
├── seed_finder/
│   └── <prefix>_raw_hits.csv         # optional: one row per seed hit
└── pipeline_info/                    # Nextflow report, timeline, trace, DAG
    ├── report.html
    ├── timeline.html
    ├── trace.txt
    └── dag.html
```

`publish_dir_mode` defaults to `copy`. On a shared filesystem, set it to
`symlink` or `rellink` if you prefer references to the work directory; see
[Pipeline parameters](../reference/parameters.md).

## `guides.ranked.tsv`

The primary result. One row per candidate guide, tab-separated, with a header.
Rows are sorted by `region_id`, then passing guides before excluded guides, then
by `rank` ascending. Within a target, `rank = 1` is the best guide according to
the composite score.

The columns fall into five groups:

| Group | Examples | What it tells you |
|---|---|---|
| Guide and location | `contig`, `start`, `stop`, `target`, `spacer`, `orientation`, `genome_start`, `genome_end` | The guide sequence and its genomic span |
| Target region | `region_id`, `chrom`, `orig_start`, `orig_end`, `input_type`, `target_ref`, `padded_start`, `padded_end`, `linked_regions` | The input target and the extracted window |
| FlashFry scores | `Doench2014OnTarget`, `DoenchCFD_specificityscore`, `Hsu2013`, `JostCRISPRi_specificityscore`, `basesDiffToClosestHit`, `FreeEnergy20`, `overflow`, `dangerous_*` | On-target and off-target evidence |
| Annotation | `distance_to_target`, `within_blast_radius`, `distance_to_atac_summit` | Position relative to the target and optional chromatin data |
| Ranking | `seed_score`, `fail_*`, `exclude_reason`, `passes_hard_filters`, `composite_score`, `composite_score_max_possible`, `composite_score_normalized`, `rank` | Whether the guide passed and how it compares |

The complete column list is in
[Output columns](../reference/output-columns.md).

### `rank` and `passes_hard_filters`

`rank` is an integer within each `region_id` for guides that pass the hard
filters. Excluded guides keep their row, but `rank` is
`exclude (<reason>)` and `passes_hard_filters` is `False`. This is deliberate:
a guide that was dropped is still evidence, and silently removing rows makes it
hard to tell "no guide was found" from "a guide was filtered out".

| Reason | Filter |
|---|---|
| `seed_score` | Three or more overlapping `GG` dinucleotides in the 12 bp PAM-proximal seed |
| `dangerous_GC` | FlashFry flagged GC content |
| `dangerous_polyT` | FlashFry flagged a poly-T / Pol III terminator |
| `multisite_genomic_match` | More than one perfect match in the genome |
| `overflow` | FlashFry's off-target search overflowed its limit |

The first failing check is reported, in that order.

### `composite_score`

Passing guides receive a weighted score. The weights are literal percentages
and are **not renormalized** when launch-dCas9 is absent, so the score keeps its
absolute meaning across runs:

| Component | Weight | Notes |
|---|---:|---|
| `launch_dCas9_enhancer` | 30% | Only when `--run_launch_dcas9` added the column |
| `JostCRISPRi_specificityscore` | 25% | |
| `Hsu2013` / 100 | 18% | Clipped to 0–1 |
| `DoenchCFD_specificityscore` | 8% | |
| `basesDiffToClosestHit` / `--max_mismatch`, capped at 1 | 5% | No off-target found → 1.0 |
| `FreeEnergy20`, min–max normalized per region | 2% | Missing → 0.5 |

`composite_score_max_possible` is 0.58 without launch-dCas9 and 0.88 with it.
`composite_score_normalized` divides by that value so the two kinds of run are
visually comparable, while `composite_score` itself stays absolute.
`launch_dCas9_promoter` is reported for reference but is deliberately not part
of the ranking.

This is a **ranking heuristic**, not a calibrated probability. It orders guides
within a target; it does not tell you that rank 1 will work and rank 2 will not.
See [Known limitations](#known-limitations).

## `guides.final.tsv`

Written only with `--run_seed_finder`. It has the same rows, in the same order,
as `guides.ranked.tsv`, plus five columns:

| Column | Meaning |
|---|---|
| `seedfinder_hit_count` | Number of non-self promoter-proximal seed hits; blank if the guide was not screened |
| `seedfinder_closest_gene` | Gene at the closest hit |
| `seedfinder_closest_tss_distance` | Signed distance to that TSS (negative upstream, positive downstream) |
| `seedfinder_closest_seed_length` | Seed length of the closest hit |
| `seedfinder_target_is_promoter` | Whether the guide's own target overlaps a promoter window |

A blank count means "not screened" (the guide failed a hard filter), which is
different from `0` ("screened, no hits"). `--seed_finder_all_guides` screens
everything.

## Seed hit table

`seed_finder/<prefix>_raw_hits.csv` is written only with
`--run_seed_finder`, one row per seed match. It includes `self_hit`, the guide
and hit coordinates, strand, seed length, nearest TSS and signed distance, every
TSS in the window, the PAM and protospacer, and a UCSC browser link. The
columns are listed in [Output columns](../reference/output-columns.md#seed-hit-table).

Self-hits are marked and excluded from `seedfinder_hit_count`. They are not
off-targets: they are the guide matching its own promoter-overlapping target.

## Provenance files

`run_manifest.json` records the pipeline name and version, commit, whether the
working tree was dirty, the command line, profiles, Nextflow version, work
directory, user, and every non-null parameter. `software_versions.yml` records
the versions reported by each process. Keep both with any published table. See
[Provenance](../reference/provenance.md).

## Interpreting a run

A few habits that prevent most misreadings:

- **Filter before ranking.** `rank` only exists for passing guides. A region
  with no passing guides is a failed design, not an empty result.
- **Check `passes_hard_filters` before sorting by `composite_score`.** Excluded
  guides have a score too, but it is not meaningful for library selection.
- **Keep `region_id`.** It is the unit of ranking and the link back to the
  input target.
- **Look at the distribution, not one row.** A region with 40 guides and a
  region with 1 guide are not equally well covered.
- **Do not compare scores across versions blindly.** The weights and filters
  are part of the method; check the changelog and the manifest.

## Known limitations

- FlashFry and launch-dCas9 are external dependencies, not bundled.
- `distance_to_atac_summit` is reported but is not part of the ranking formula.
- hg38 is the tested reference. Other genomes require a matching FlashFry index
  and are not validated.
- The composite score is a ranking heuristic, not a validated predictor of
  editing or repression efficiency.
