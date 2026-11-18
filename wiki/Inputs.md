# Inputs

creTarget accepts two target types, and both may be supplied in one run:

| Type | Text form | BED form | Flag |
|---|---|---|---|
| SNP / point | `chr:pos` (1-based) | `chrom start end [name]` (0-based, half-open) | `--snp_bed` |
| Enhancer / interval | `chr:start-end` (BED convention) | `chrom start end [name]` (0-based, half-open) | `--enhancer_bed` |

At least one is required. The format is detected per line, not per file, so a
file may mix text and BED rows, though keeping one convention per file is
easier to review.

## Coordinates

This is the part worth getting right before a long run.

- **`chr:pos` SNP text is 1-based.** `chr1:120000` is the 120,000th base.
  Internally it becomes `[119999, 120000)`.
- **`chr:start-end` enhancer text uses BED convention: 0-based, half-open.**
  `chr1:10000-10400` covers bases 10,001–10,400 in 1-based coordinates.
- **BED rows are always 0-based, half-open**, for both input types. A one-base
  SNP at 1-based position 120,000 is `chr1 119999 120000`.
- **`target_ref`** is the enhancer midpoint (`(start + end) // 2`, 0-based) or
  the SNP point (0-based start). It is the reference for `distance_to_target`.

## Examples

```text
# enhancers.txt
chr1:10000-10400
chr1:50000-50300
chr2:20000-20500
```

```text
# snps.txt
chr1:120000
chr2:60000
```

```text
# enhancers.bed
chr1    10000    10400    enh_1
chr1    50000    50300    enh_2
chr2    20000    20500    enh_3
```

```text
# snps.bed — one base, 0-based half-open
chr1    119999    120000    rs123456
chr2    59999     60000     rs654321
```

The fourth BED column is used as the `region_id`. If it is missing, a name is
generated. Extra columns are ignored. A BED line is recognized by a tab; spaces
alone do not make a line BED.

## What happens to the input

1. **Comments and tracks are skipped.** Blank lines and lines starting with
   `#`, `track`, or `browser` are ignored.
2. **Malformed lines are skipped with a warning.** A bad line does not abort
   the run. The number skipped is reported. A file with no usable regions is a
   hard error.
3. **Contigs are validated against the reference `.fai`.** Regions on absent
   contigs are dropped with a warning. The common cause is Ensembl-style `1`
   against UCSC-style `chr1`; creTarget says so explicitly.
4. **Duplicate names are made unique** with `_dup1`, `_dup2`, and so on.
5. **Regions are padded.** Enhancers by `--flank_bp` (default 50) per side;
   SNPs by `--snp_window_bp` (default 150) per side. Padding is clipped to the
   chromosome length.
6. **Linked targets are flagged, if requested.** `--dedup_blast_radius` marks
   same-contig targets within `--blast_radius_bp` (default 1000) in
   `linked_regions`. It flags them; it does not merge or remove them.
7. **Sequence is extracted with `bedtools getfasta`.** The FASTA header encodes
   the region name and padded coordinates.

## Naming

`region_id` identifies the target, not the guide. It is the fourth BED column
when present, otherwise a generated name or the original text line. Every guide
for a region shares the `region_id`, and ranking is per `region_id`.

The metadata table (`regions/regions.metadata.tsv`) records `region_id`,
`fasta_header`, `chrom`, `orig_start`, `orig_end`, `input_type`, `target_ref`,
`padded_start`, `padded_end`, and `linked_regions`.

## Step-level CLI

The CLI uses hyphens:

```bash
cretarget prep-regions \
    --snp-bed snps.txt \
    --enhancer-bed enhancers.txt \
    --fasta hg38.fa \
    --outdir prep
```

Nextflow uses underscores (`--snp_bed`, `--enhancer_bed`, `--fasta`).

## Pitfalls

- **Coordinates look shifted by one.** Check whether the source is 1-based
  inclusive or 0-based half-open.
- **Everything was dropped as an unknown contig.** Check the FASTA `.fai`
  naming. `chr1` and `1` are different.
- **A process substitution was passed as input.** Nextflow stages inputs by file
  name; `<(...)` does not survive. Write a real file.
- **A SNP BED row has start == end.** BED is half-open, so that is empty. Use
  `chr1 119999 120000`.
- **The wrong parser was used.** A `chr:pos` line in `--enhancer_bed` is skipped
  with a warning, and vice versa. Put each type in its own file.

Full guide: [Input formats](https://mikekazz.github.io/creTarget/guides/inputs/)
