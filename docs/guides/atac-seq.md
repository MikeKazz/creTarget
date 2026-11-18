# Optional: ATAC-seq annotation

`--atac_narrowpeak` adds one column, `distance_to_atac_summit`: the absolute
distance in base pairs from the guide's midpoint to the nearest ATAC-seq summit
on the same contig.

This is **annotation, not ranking**. The distance is reported so you can filter
or inspect it yourself; it does not change `composite_score` or `rank`. The
rationale is that accessibility is a property of the target region and the cell
type, not a universal guide-quality term, and silently folding it into the
score would hide a biological assumption.

## Use a bundled cell type

creTarget ships with fetch presets for 19 hematopoietic cell types:

```bash
cretarget fetch-atac --cell-type k562
```

Available names:

`b_cell`, `basophil`, `cd34_bone_marrow`, `cd34_cord_blood`, `cd4_tcell`,
`cd8_tcell`, `clp`, `cmp`, `erythroid`, `gmp`, `hsc`, `k562`, `lmpp`,
`megakaryocyte`, `mep`, `monocyte`, `mpp`, `neutrophil`, `nk_cell`.

The files are downloaded from Zenodo record
[10.5281/zenodo.21268805](https://doi.org/10.5281/zenodo.21268805) into
`~/.cretarget/atac_peaks` and checked against the record's checksum manifest.
The `--cell-type` flag is repeatable if you want more than one file.

Then:

```bash
nextflow run MikeKazz/creTarget -r v0.5.1 \
    -profile slurm,conda \
    --enhancer_bed enhancers.txt \
    --fasta ~/.cretarget/reference/hg38.fa \
    --flashfry_jar /path/to/FlashFry-assembly-1.15.jar \
    --flashfry_db ~/.cretarget/reference/hg38_standard_cas9ngg_db \
    --atac_narrowpeak ~/.cretarget/atac_peaks/k562.narrowPeak \
    --outdir results
```

## Use your own narrowPeak

Pass any 10-column ENCODE-standard narrowPeak file:

```bash
--atac_narrowpeak /path/to/my_peaks.narrowPeak
```

The expected columns are:

| Column | Meaning |
|---|---|
| 1 | `chrom` |
| 2 | `start` (0-based) |
| 3 | `end` |
| 4 | `name` |
| 5 | `score` |
| 6 | `strand` |
| 7 | `signalValue` |
| 8 | `pValue` |
| 9 | `qValue` |
| 10 | `peak` — summit offset from `start`, or `-1` |

When column 10 is `-1` (no summit called), creTarget uses the peak midpoint.
Comment lines starting with `#` are ignored.

## How the distance is computed

For each guide, creTarget takes the guide midpoint
(`(genome_start + genome_end) / 2`) and finds the nearest summit on the same
contig using a vectorized search. Guides on contigs with no peaks get `NaN`,
and the step reports how many. A distance of `0` means the guide midpoint sits
exactly on a summit.

## Interpreting the result

- Small distances are evidence that the guide sits in an accessible region in
  the assayed cell type.
- Large distances do not necessarily mean inaccessible chromatin; the peak set
  is a snapshot from one cell type and one pipeline.
- A `NaN` is not zero. It means there were no peaks on that contig in the file,
  so no distance could be computed.

If you want accessibility to affect selection, filter the ranked table on
`distance_to_atac_summit` after the run and record the threshold in your
methods. Do not describe a threshold as part of creTarget's ranking unless you
have implemented and validated it.

## Bundled data attribution

The bundled peaks are derived from public ATAC-seq data. If you use them,
follow the terms of the original data source and cite the dataset and ENCODE in
addition to creTarget. The Zenodo record is a redistribution for convenience,
not a new measurement.
