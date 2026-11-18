# Input formats

creTarget accepts two target types, and both may be supplied in the same run:

| Type | Text form | BED form | Flag |
|---|---|---|---|
| SNP / point | `chr:pos` (1-based) | `chrom start end [name]` (0-based, half-open) | `--snp_bed` |
| Enhancer / interval | `chr:start-end` | `chrom start end [name]` (0-based, half-open) | `--enhancer_bed` |

At least one of the two is required. The format is detected per line, not per
file, so a file can contain a mix of text and BED rows (though keeping one
convention per file is easier to review).

## Coordinate conventions

This is the part worth getting right before a long run.

- **`chr:pos` SNP text** is **1-based**. `chr1:120000` is the 120,000th base.
  Internally it becomes the half-open interval `[119999, 120000)`.
- **`chr:start-end` enhancer text** uses the same numbers as BED:
  **0-based, half-open**. `chr1:10000-10400` covers bases 10,001–10,400 in
  1-based coordinates. This matches the convention used by genome browsers and
  BED files, but it is not the same as a 1-based inclusive range.
- **BED rows** are always 0-based, half-open, for both input types. A one-base
  SNP at 1-based position 120,000 is `chr1  119999  120000`.
- **`target_ref`** in the metadata is the enhancer midpoint
  (`(start + end) // 2`, 0-based) or the SNP point (0-based start). It is the
  reference position used for `distance_to_target`.

If your coordinates come from a paper or a 1-based tool, convert them before
writing the file. A quick sanity check: an enhancer interval's end should be
greater than its start, and a SNP point should be inside the chromosome length
in the FASTA `.fai`.

## Examples

### Text input

```text
# enhancers.txt — chr:start-end, BED convention
chr1:10000-10400
chr1:50000-50300
chr2:20000-20500
```

```text
# snps.txt — chr:pos, 1-based
chr1:120000
chr2:60000
```

### BED input

```text
# enhancers.bed — chrom, start, end, name
chr1    10000    10400    enh_1
chr1    50000    50300    enh_2
chr2    20000    20500    enh_3
```

```text
# snps.bed — one base, 0-based half-open
chr1    119999    120000    rs123456
chr2    59999     60000     rs654321
```

The fourth column is used as the `region_id`. If it is missing, a name is
generated. Extra columns are ignored. A BED line is recognized by the presence
of a tab; spaces alone do not make a line BED.

## What creTarget does with the input

1. **Skips comments and browser tracks.** Blank lines and lines starting with
   `#`, `track`, or `browser` are ignored.
2. **Skips malformed lines with a warning.** A bad line does not abort the run.
   At the end, the number of skipped lines is reported. A file with no usable
   regions is a hard error.
3. **Validates contigs against the reference `.fai`.** Regions on contigs
   absent from the reference are dropped with a warning. The common cause is
   Ensembl-style `1` against a UCSC-style `chr1` reference (or the reverse), and
   creTarget says so explicitly instead of letting bedtools produce a confusing
   error later.
4. **Makes duplicate names unique.** A repeated `region_id` gets `_dup1`,
   `_dup2`, and so on. Without this, the metadata join would silently keep the
   last one.
5. **Pads each region.** Enhancers are padded by `--flank_bp` (default 50) on
   both sides. SNPs are padded by `--snp_window_bp` (default 150) on both
   sides. Padding is clipped to the chromosome length.
6. **Flags linked targets, if requested.** `--dedup_blast_radius` marks
   same-contig targets within `--blast_radius_bp` (default 1000) of each other
   in the `linked_regions` column. It **flags** them; it does not merge or
   remove them. A guide against one target may silence the other, so the flag
   is a warning to review the pair.
7. **Extracts sequence with `bedtools getfasta`.** The FASTA header encodes the
   region name and padded coordinates, so genome coordinates survive into
   FlashFry's output.

## Naming and downstream columns

`region_id` is the target identifier, not the guide identifier. It is the
fourth BED column when present, otherwise a generated name such as
`enhancer_3` or the original text line (for text input). Every guide for a
region shares the same `region_id`, and ranking is per `region_id`.

The metadata table (`regions/regions.metadata.tsv`) records:

| Column | Meaning |
|---|---|
| `region_id` | Unique target name |
| `fasta_header` | Encoded header used internally |
| `chrom` | Contig |
| `orig_start`, `orig_end` | Input coordinates (BED convention) |
| `input_type` | `snp` or `enhancer` |
| `target_ref` | Point used for distance calculations |
| `padded_start`, `padded_end` | Extracted window, clipped to the contig |
| `linked_regions` | Comma-separated targets within the blast radius, if enabled |

## Inputs for the step-level CLI

The step-level CLI uses flags with hyphens:

```bash
cretarget prep-regions \
    --snp-bed snps.txt \
    --enhancer-bed enhancers.txt \
    --fasta hg38.fa \
    --outdir prep
```

Nextflow passes the same values with underscores (`--snp_bed`, `--enhancer_bed`,
`--fasta`). See [Step-level CLI](../reference/cli.md) for the full mapping.

## Common pitfalls

- **Coordinates look shifted by one.** Check whether your source is 1-based
  inclusive or 0-based half-open. BED is 0-based half-open; `chr:pos` is
  1-based.
- **Everything was dropped as an unknown contig.** Check the FASTA `.fai`
  naming. `chr1` and `1` are different contigs.
- **A process substitution was passed as input.** Nextflow stages inputs by
  file name; `<(...)` does not survive. Write a real file first.
- **A SNP BED was written as a point with start == end.** BED is half-open, so
  `chr1 120000 120000` is empty. Use `chr1 119999 120000`.
- **The wrong parser was used.** A `chr:pos` line in `--enhancer_bed` is
  reported and skipped; a `chr:start-end` line in `--snp_bed` is likewise
  reported. Put each type in its own file under the right flag.
