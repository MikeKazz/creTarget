# Seed screening

`--run_seed_finder` asks a different question from FlashFry. It does not ask
*where else could this guide cut*. It asks *where else could this guide sit on
a promoter and repress*.

dCas9 does not cut; it blocks transcription. A guide only needs a perfect
PAM-proximal seed next to any NGG to load onto DNA. If that happens over a
promoter, the gene can be repressed whether or not the rest of the guide
matches. FlashFry's mismatch-based scores do not capture that failure mode, so
the two screens are complementary.

## Attribution

The method is from **Hartman et al. 2026**, *Systematic identification of
seed-driven off-target effects in Perturb-seq experiments*
([doi:10.64898/2026.03.27.714658](https://doi.org/10.64898/2026.03.27.714658)).
Pipeline: [AustinHartman/perturb_seed](https://github.com/AustinHartman/perturb_seed).

The promoter-window reference design and single-pass batch architecture come
from **seed-finder-lite** (Kazemi, M.,
[doi:10.5281/zenodo.21599135](https://doi.org/10.5281/zenodo.21599135)). The
scanner in creTarget is creTarget's own code.

Cite Hartman et al. alongside creTarget when you use this feature.

## Build the promoter reference

```bash
cretarget build-seed-reference --fasta ~/.cretarget/reference/hg38.fa
```

This is built once from a GENCODE annotation and the genome FASTA; it is not
downloaded as a blob.

1. Downloads a pinned GENCODE annotation (release 47 by default; use
   `--gencode-release N`, or `--gtf` for your own file).
2. Derives a TSS table, keeping gene name and strand.
3. Extracts ±2 kb windows around every TSS, merges overlapping windows per
   chromosome, and pads each by 23 bp so a site straddling a boundary is not
   truncated.
4. Writes `tss_table.tsv`, `promoter_windows.fa.gz`, and `provenance.json` under
   `~/.cretarget/seed_reference`.

Rebuild with `--force` after changing the annotation. A reference built before
creTarget 0.5.0 has no strand column and is rejected with instructions to
rebuild; defaulting the strand would silently invert distances for every
minus-strand gene.

## Run

```bash
nextflow run MikeKazz/creTarget -r v0.5.1 \
    -profile slurm,conda \
    --enhancer_bed enhancers.txt \
    --fasta ~/.cretarget/reference/hg38.fa \
    --flashfry_jar /path/to/FlashFry-assembly-1.15.jar \
    --flashfry_db ~/.cretarget/reference/hg38_standard_cas9ngg_db \
    --run_seed_finder \
    --seed_finder_min_length 8 \
    --seed_finder_tss_window 1000 \
    --outdir results
```

`--seed_finder_min_length` and `--seed_finder_tss_window` have **no defaults on
purpose**. They change what the screen means, so choose them deliberately.

| Parameter | Meaning | Guidance |
|---|---|---|
| `--seed_finder_min_length` | Minimum exact-match PAM-proximal seed length | Hartman et al. define the seed as roughly the 10–12 PAM-proximal bases. `8` is permissive; `10`–`12` is closer to their definition; `16` is stringent. Report the value. |
| `--seed_finder_tss_window` | Maximum distance from a TSS for a hit | The reference extends 2000 bp. A larger value finds nothing extra; creTarget says so. |
| `--seed_finder_max_guides` | Split the screen into tasks of at most N guides | Memory is ~1 GB RSS per million guides. Use when memory, not time, is binding. |
| `--seed_finder_all_guides` | Screen guides that failed the hard filters | Off by default. Only passing guides go into a library. |

Only guides that pass the hard filters are screened. Excluded guides get blank
seed columns rather than a count of zero, so "looked and found nothing" stays
distinguishable from "did not look".

## Self-hits

Many cCREs are promoter-like, so a guide targeting one finds *itself* in the
promoter reference: a perfect 20/20 match that is not an off-target. In a
genome-wide cCRE run these were **73% of all hits**.

Self-hits are kept in the hit table with `self_hit=yes`, excluded from
`seedfinder_hit_count`, and surfaced as `seedfinder_target_is_promoter`.
Detection tests containment in the guide's own span rather than an exact start,
because the PAM sits at the high coordinates for a forward guide and the low
coordinates for a reverse guide.

## Why it runs once

The scan is one pass over the promoter reference and barely moves with guide
count. Running it per pipeline chunk paid that fixed pass on every chunk; at
`--chunk_size 1000`, a genome-wide run would have spent weeks re-reading the
same ~212 MB. Hits are streamed to CSV as they are found; only the seed index is
held in memory.

## Output

```text
results/guides.final.tsv
results/seed_finder/<prefix>_raw_hits.csv
```

The five summary columns are described in [Outputs](Outputs#guidesfinaltsv).

## Verify a run

```bash
bash assets/verify_seed_run.sh results/
```

Read-only and streaming, so it is safe on a 50+ GB table. It checks hits on
**both strands** (a one-sided strand column is the signature of a scanner that
indexes guides under one orientation, silently halving recall), multi-gene
hits, signed distances, self-hit accounting, and row invariants.

## Interpreting

- `0` means screened and no promoter-proximal seed match.
- Blank means not screened (failed a hard filter).
- A high count is a warning, not a verdict. It counts plausible binding sites.
- `seedfinder_target_is_promoter=yes` means the target overlaps a promoter.
  That is often expected for promoter-like cCREs.
- `tss_distance` is signed relative to the gene's strand: negative upstream.
- `tss_genes_in_window` lists every gene whose promoter contains the site,
  nearest first. A bidirectional promoter can put a site in two genes.

## Troubleshooting

| Symptom | Likely cause |
|---|---|
| `promoter-proximal seed reference is not built` | Run `cretarget build-seed-reference --fasta ...` |
| `three-column TSS table` error | Reference predates 0.5.0; rebuild with `--force` |
| `no guides to screen` | Every guide failed the hard filters; use `--seed_finder_all_guides` only if intentional |
| A window above 2000 finds nothing extra | Reference is built at ±2 kb; expected and reported |
| Out of memory on a huge table | Set `--seed_finder_max_guides` to shard |

Full guide: [Seed screening](https://mikekazz.github.io/creTarget/guides/seed-screening/)
