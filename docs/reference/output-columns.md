# Output columns

Column names are the stable interface. Column **order** may change between
releases; parse by name, not position.

## `guides.ranked.tsv`

### Guide and location

| Column | Type | Meaning |
|---|---|---|
| `contig` | string | Encoded region header from `prep-regions`: `<region_id>__<chrom>:<padded_start>-<padded_end>` |
| `start` | int | FlashFry offset within the extracted window, **not** a genome coordinate |
| `stop` | int | FlashFry offset within the extracted window, **not** a genome coordinate |
| `target` | string | FlashFry target sequence (the 20 nt protospacer) |
| `context` | string | Target plus PAM context |
| `orientation` | string | Guide orientation reported by FlashFry (`+`/`-`) |
| `spacer` | string | First 20 nt of `target`; the guide sequence used for scoring and seed screening |
| `genome_chrom` | string | Chromosome recovered from `contig` |
| `window_padded_start` | int | Start of the extracted window, 0-based |
| `window_padded_end` | int | End of the extracted window, 0-based |
| `genome_start` | int | Guide start in genome coordinates |
| `genome_end` | int | Guide end in genome coordinates |
| `fasta_header` | string | Header recorded in the metadata table |
| `chrom` | string | Chromosome from the region metadata |
| `orig_start` | int | Input start (BED convention) |
| `orig_end` | int | Input end (BED convention) |
| `input_type` | string | `snp` or `enhancer` |
| `target_ref` | int | Point used for distance calculations: SNP point or enhancer midpoint |
| `padded_start` | int | Region start after padding |
| `padded_end` | int | Region end after padding |
| `linked_regions` | string | Comma-separated targets within the blast radius; empty when flagging is off |

### FlashFry scores

These are passed through from FlashFry 1.15. See FlashFry's documentation for
their definitions.

| Column | Meaning |
|---|---|
| `overflow` | `OK` or an overflow flag from the off-target search |
| `Doench2014OnTarget` | Doench 2014 on-target score |
| `DoenchCFD_maxOT` | Maximum CFD off-target score |
| `DoenchCFD_specificityscore` | Doench 2016 CFD specificity |
| `dangerous_GC` | FlashFry GC-content flag (`NONE` or a warning) |
| `dangerous_polyT` | FlashFry poly-T / Pol III terminator flag (`NONE` or a warning) |
| `dangerous_in_genome` | Perfect-match count, e.g. `IN_GENOME=1` or `NONE` |
| `Hsu2013` | Hsu 2013 specificity score |
| `Moreno-Mateos2015OnTarget` | Moreno-Mateos 2015 on-target score |
| `JostCRISPRi_maxOT` | Maximum Jost CRISPRi off-target score |
| `JostCRISPRi_specificityscore` | Jost CRISPRi specificity |
| `basesDiffToClosestHit` | Mismatches to the closest off-target hit |
| `closestHitCount` | Number of hits at that distance |
| `0-1-2-3-4_mismatch` | Off-target count by mismatch bucket |
| `FreeEnergy20` | RNA folding free energy |
| `otCount` | Off-target count |

### Annotation

| Column | Type | Meaning |
|---|---|---|
| `distance_to_target` | int | Absolute bp from the guide midpoint to `target_ref` |
| `within_blast_radius` | bool | `distance_to_target <= --blast_radius_bp` |
| `distance_to_atac_summit` | float / NaN | Absolute bp to the nearest ATAC summit; NaN when no peaks on the contig |

### Ranking

| Column | Type | Meaning |
|---|---|---|
| `seed_score` | int | Overlapping `GG` dinucleotides in the 12 bp PAM-proximal seed |
| `in_genome_count` | int | Parsed perfect-match count (`NONE` → 0) |
| `fail_seed_score` | bool | `seed_score >= 3` |
| `fail_dangerous_GC` | bool | `dangerous_GC != NONE` |
| `fail_dangerous_polyT` | bool | `dangerous_polyT != NONE` |
| `fail_multisite_genomic_match` | bool | `in_genome_count > 1` |
| `fail_overflow` | bool | `overflow != OK` |
| `exclude_reason` | string / empty | First failed filter; empty for passing guides |
| `passes_hard_filters` | bool | True when no filter failed |
| `composite_score` | float | Weighted score; see [Output files](../guides/outputs.md#composite_score) |
| `composite_score_max_possible` | float | 0.58 without launch-dCas9, 0.88 with it |
| `composite_score_normalized` | float | `composite_score / composite_score_max_possible` |
| `rank` | int / string | Per-`region_id` rank among passing guides, or `exclude (<reason>)` |

### Optional launch-dCas9

| Column | Type | Meaning |
|---|---|---|
| `launch_dCas9_enhancer` | float / NaN | Predicted enhancer fitness; contributes 30% when present |
| `launch_dCas9_promoter` | float / NaN | Predicted promoter fitness; reported, not ranked |

### Optional seed screening

Added to `guides.final.tsv`:

| Column | Type | Meaning |
|---|---|---|
| `seedfinder_hit_count` | int / blank | Non-self promoter-proximal seed hits; blank when not screened |
| `seedfinder_closest_gene` | string / blank | Gene at the closest hit |
| `seedfinder_closest_tss_distance` | int / blank | Signed distance to that TSS (negative upstream) |
| `seedfinder_closest_seed_length` | int / blank | Seed length of the closest hit |
| `seedfinder_target_is_promoter` | `yes` / `no` | Whether the guide's own target overlaps a promoter window |

## Seed hit table

`seed_finder/<prefix>_raw_hits.csv` has one row per seed hit. Coordinates are
1-based inclusive and cover the protospacer only; the PAM is reported
separately.

| Column | Meaning |
|---|---|
| `self_hit` | `yes` when the hit is the guide's own target site; excluded from the count |
| `guide_id` | Unique guide identifier (`<region>_<start>_<end>_<orientation>_<index>`, with a shard prefix when sharded) |
| `chr` | Hit chromosome |
| `start`, `end` | Hit protospacer coordinates, 1-based inclusive |
| `strand` | Hit strand (`+`/`-`) |
| `seed_length` | Consecutive PAM-proximal matching bases |
| `nearest_tss_gene` | Gene at the closest TSS |
| `tss_distance` | Signed distance to that TSS; negative upstream, positive downstream, relative to the gene's strand |
| `tss_strand` | Strand of that gene |
| `tss_genes_in_window` | Semicolon-separated genes whose promoters contain the hit, nearest first |
| `tss_distances_in_window` | Matching signed distances, same order |
| `pam` | PAM sequence as read on the guide's strand |
| `pam_start`, `pam_end` | PAM coordinates, 1-based inclusive |
| `protospacer` | Genomic protospacer sequence oriented to the guide |
| `guide_sequence` | The guide's 20 nt spacer |
| `ucsc_url` | UCSC Genome Browser link for the hit |

## `run_manifest.json`

| Section | Fields |
|---|---|
| `pipeline` | `name`, `version`, `revision`, `commit_id`, `commit_source`, `uncommitted_changes`, `home_page` |
| `run` | `run_name`, `session_id`, `start_time`, `command_line`, `profile`, `container_engine`, `nextflow_version`, `work_dir`, `launch_dir`, `user` |
| `params` | Every non-null pipeline parameter |

`commit_source` is `nextflow` when the pipeline was fetched by
`nextflow run owner/repo`, `git(projectDir)` when it was run from a local clone,
and `none` when neither is available. `uncommitted_changes` is `true`, `false`,
or `null` if git is unavailable.

## `software_versions.yml`

One block per process, with the tool versions that process reported. Blocks are
deduplicated by process name, so repeated per-chunk blocks do not repeat.

## `regions/regions.metadata.tsv`

One row per input target, with the columns described in
[Input formats](../guides/inputs.md#naming-and-downstream-columns).

## `pipeline_info/`

Nextflow's own reports: `report.html`, `timeline.html`, `trace.txt`, and
`dag.html`. They are generated by the default `nextflow.config` and are useful
for resource tuning and debugging.
