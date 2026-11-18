# Optional: promoter-proximal seed screening

`--run_seed_finder` asks a different question from FlashFry. It does not ask
*where else could this guide cut*. It asks *where else could this guide sit on
a promoter and repress*.

dCas9 does not cut; it blocks transcription. A guide only needs a perfect
PAM-proximal seed next to any NGG to load onto DNA. If that happens over
somebody's promoter, the gene can be repressed whether or not the rest of the
guide matches. FlashFry's mismatch-based scores do not capture that failure
mode, so the two screens are complementary and neither replaces the other.

## Attribution

The method is from **Hartman et al. 2026**, *Systematic identification of
seed-driven off-target effects in Perturb-seq experiments*
([doi:10.64898/2026.03.27.714658](https://doi.org/10.64898/2026.03.27.714658)).
Their pipeline is at
[AustinHartman/perturb_seed](https://github.com/AustinHartman/perturb_seed).

**seed-finder-lite** (Kazemi, M.,
[doi:10.5281/zenodo.21599135](https://doi.org/10.5281/zenodo.21599135)) is a
command-line implementation of that method. The merged promoter-window
reference design and the single-pass batch architecture used here come from it.
The scanner in creTarget is creTarget's own code.

If you use `--run_seed_finder`, cite Hartman et al. alongside creTarget.

## Build the promoter reference

This is built once from a GENCODE annotation and the genome FASTA you already
have; it is not downloaded as a blob:

```bash
cretarget build-seed-reference --fasta ~/.cretarget/reference/hg38.fa
```

What it does:

1. Downloads a pinned GENCODE annotation (release 47 by default; use
   `--gencode-release N` to pin another, or `--gtf` for your own file).
2. Derives a TSS table, keeping the gene name and strand. The TSS is the
   5'-most base: the start for a `+` gene, the end for a `-` gene.
3. Extracts ±2 kb windows around every TSS, merges overlapping windows per
   chromosome, and pads each merged window by 23 bp (one protospacer plus PAM)
   so a site straddling a boundary is not truncated.
4. Writes `tss_table.tsv`, `promoter_windows.fa.gz`, and `provenance.json` under
   `~/.cretarget/seed_reference`.

`provenance.json` records the release, source URL, GTF SHA-256, genome FASTA,
TSS count, window count, and the creTarget version. The promoter definition is
therefore reproducible from declared inputs rather than being a property of a
downloaded file.

Rebuild with `--force` after changing the annotation. A reference built before
creTarget 0.5.0 has no strand column and is rejected with instructions to
rebuild; defaulting the strand would silently invert distances for every
minus-strand gene.

## Run the screen

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

| Parameter | Meaning | Practical guidance |
|---|---|---|
| `--seed_finder_min_length` | Minimum exact-match PAM-proximal seed length | Hartman et al. define the seed as roughly the 10–12 PAM-proximal bases and search permissively. `8` is permissive; `10`–`12` is closer to their definition; `16` is stringent. Report the value. |
| `--seed_finder_tss_window` | Maximum distance from a TSS for a hit to count | The reference extends 2000 bp from each TSS. A larger value finds nothing extra; creTarget says so rather than silently capping. |
| `--seed_finder_max_guides` | Split the screen into tasks of at most N guides | Memory is ~1 GB of RSS per million guides. Use this when memory, not time, is the constraint. |
| `--seed_finder_all_guides` | Screen guides that failed the hard filters too | Off by default. Only passing guides go into a library. |

Only guides that pass the hard filters are screened. Excluded guides get blank
seed columns rather than a count of zero, so "we looked and found nothing" stays
distinguishable from "we did not look".

## Self-hits

Many cCREs are promoter-like, so a guide targeting one finds *itself* in the
promoter reference: a perfect 20/20 match that is not an off-target. In a
genome-wide cCRE run these were **73% of all hits**.

Self-hits are:

- kept in the hit table with `self_hit=yes`;
- excluded from `seedfinder_hit_count`;
- surfaced as `seedfinder_target_is_promoter`, which is useful on its own: it
  tells you the target region overlaps an annotated promoter.

Detection tests whether the hit falls inside the guide's own span rather than
matching an exact start, because the PAM sits at the high coordinates for a
forward guide and the low coordinates for a reverse guide.

## Why it runs once

The scan is one pass over the promoter reference and barely moves with guide
count: the per-position work is a dictionary lookup whether there are twenty
guides or two million. Running it per pipeline chunk paid that fixed pass on
every chunk. At `--chunk_size 1000`, a genome-wide run would have spent weeks
re-reading the same ~212 MB.

The table is streamed and hits are written as they are found. The only thing
held in memory is the seed index, measured at roughly **1 GB of RSS per million
guides** screened. `--seed_finder_max_guides N` splits the work into shards when
that becomes the constraint: screening 68 million guides, which would need
~68 GB in one task, is comfortable at `--seed_finder_max_guides 4000000` (~4 GB
per task). Each shard repeats the reference pass, so raise the limit until
memory is binding rather than time.

## Output

The screen adds five columns to `guides.final.tsv` and writes a hit table:

```text
results/guides.final.tsv
results/seed_finder/<prefix>_raw_hits.csv
```

The five summary columns are described in
[Output files](outputs.md#guidesfinaltsv). The hit table columns are in
[Output columns](../reference/output-columns.md#seed-hit-table).

## Verify a run

`assets/verify_seed_run.sh` checks the output after a run:

```bash
bash assets/verify_seed_run.sh results/
```

It is read-only and streams, so it is safe on a 50+ GB table. Besides the
table-shape invariants, it checks that hits appear on **both strands**: a
one-sided strand column is the signature of a scanner that indexes guides under
a single orientation, which silently halves recall while exiting 0. It also
checks for multi-gene hits and signed distances.

## Interpreting the result

- A count of `0` means screened and no promoter-proximal seed match was found.
- A blank count means the guide was not screened (it failed a hard filter).
- A high count is a warning, not a verdict. The score is not calibrated to
  repression magnitude; it counts plausible binding sites.
- `seedfinder_target_is_promoter=yes` means the target itself overlaps an
  annotated promoter. That is often expected for promoter-like cCREs.
- `tss_distance` is signed: negative upstream, positive downstream, relative to
  the gene's own strand. This matches Hartman et al.'s convention and makes the
  distance comparable with their published repression curve.
- `tss_genes_in_window` and `tss_distances_in_window` list every gene whose
  promoter the site falls in, nearest first. A bidirectional promoter can put a
  site in two genes at once.

## Troubleshooting

| Symptom | Likely cause |
|---|---|
| `promoter-proximal seed reference is not built` | Run `cretarget build-seed-reference --fasta ...` |
| `three-column TSS table` error | The reference predates 0.5.0; rebuild with `--force` |
| `no guides to screen` | Every guide failed the hard filters; use `--seed_finder_all_guides` only if that is intentional |
| A window above 2000 finds nothing extra | The reference is built at ±2 kb; this is expected and reported |
| Out of memory on a huge table | Set `--seed_finder_max_guides` to shard the screen |
