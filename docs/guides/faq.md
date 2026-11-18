# FAQ

## What does creTarget actually do?

It turns target regions into candidate CRISPRi guides, scores them with
FlashFry, optionally annotates them with ATAC-seq, launch-dCas9, and
promoter-proximal seed data, then applies hard filters and ranks the survivors
within each target. The output is a table you can filter and order as a library.

It does **not** synthesize oligos, clone libraries, or run a screen for you.

## Which genomes are supported?

hg38 is the tested reference. The pipeline itself is genome-agnostic if you
provide:

1. a FASTA with a `.fai` index,
2. a FlashFry index built from that FASTA, and
3. contig names in the input that match the FASTA exactly.

Other genomes are not validated in CI. The seed reference builder currently
expects standard human chromosomes and a GENCODE/Ensembl-style GTF.

## SNP or enhancer input?

Use `--snp_bed` for points and `--enhancer_bed` for intervals. Both can be used
in one run. `--snp_bed` also accepts BED; `--enhancer_bed` also accepts BED. The
format is detected per line, not per file. See
[Input formats](inputs.md).

## Is a cluster required?

No. `-profile local` runs on a laptop or workstation. A cluster is only needed
for large libraries or for many concurrent tasks. The same command works on
SLURM, PBS/PBS Pro, SGE, and LSF by changing the profile.

## Is a GPU required?

No. The core pipeline is CPU-only. A GPU only helps `--run_launch_dcas9`; that
step also runs on CPU, just more slowly.

## Why are there two off-target screens?

FlashFry asks where else a guide could **cut**. Seed screening asks where else a
guide could **bind a promoter and repress**. dCas9 does not cut, so the second
failure mode is real and not captured by mismatch-based scores. The screens are
complementary. See [Seed screening](seed-screening.md).

## What does `rank` mean?

Within each `region_id`, passing guides are ordered by `composite_score`; `1` is
the highest score. Excluded guides have `rank = exclude (<reason>)` instead.
`rank` is a within-target ordering, not a probability that a guide will work.

## Why are excluded guides still in the table?

Because a filtered guide is evidence. Removing rows makes it impossible to tell
"no guide was found" from "guides were found but all failed a filter". Use
`passes_hard_filters` or filter on `rank`.

## Why is a seed count blank instead of zero?

A blank count means the guide was **not screened** (it failed a hard filter).
`0` means it was screened and no promoter-proximal seed match was found. The
distinction matters when you audit a library. `--seed_finder_all_guides` screens
everything.

## Why is `distance_to_atac_summit` NaN?

There were no peaks on that contig in the narrowPeak file. `NaN` means "not
computed", not "zero distance". The step reports how many guides are affected.

## Why does the composite score stop at 0.58 or 0.88?

The weights are literal percentages and are not renormalized when optional
components are absent. Without launch-dCas9, the available weights sum to 0.58;
with it, 0.88. `composite_score_normalized` divides by the applicable maximum
so runs with and without the optional score are visually comparable. See
[Output files](outputs.md#composite_score).

## How do I choose `--chunk_size`?

Default 5,000. Smaller chunks parallelize better and retry faster; larger chunks
add less scheduling overhead but make a failure more expensive. FlashFry runtime
is superlinear in chunk size, so there is a point where bigger is slower.
Measure on your cluster. See [Scaling and clusters](scaling-and-clusters.md).

## Can I use my own ATAC-seq file?

Yes. Pass a 10-column ENCODE narrowPeak file with `--atac_narrowpeak`. The
distance is annotation only; it does not enter `composite_score`. See
[ATAC-seq annotation](atac-seq.md).

## Can I use my own gene annotation for seed screening?

Yes. `cretarget build-seed-reference --fasta ... --gtf your.gtf`. The GTF must
contain `gene` records with gene names and strands. The release and checksum are
recorded in `provenance.json`.

## Is creTarget on PyPI? On Bioconda?

The release workflow is set up to publish tagged releases to PyPI. If `pip install
cretarget` reports no matching distribution, the release is not on PyPI yet;
use the tagged source install in [Installation](../getting-started/installation.md).

There is no `cretarget` Bioconda package yet. The `conda`/`mamba`/`micromamba`
Nextflow profiles build the task environment; they do not install a conda
package named `cretarget`.

## Does it support Cas12a or other PAMs?

The current pipeline is built around SpCas9 NGG guides: 20 nt spacers, a
3 bp PAM, and a FlashFry `spcas9ngg` index. Other enzymes would require a
matching FlashFry index and validation of the seed definitions. They are not
supported today.

## Can I run it offline?

Yes, after the explicit fetch steps. Run `cretarget setup-reference`,
`fetch-models`, `fetch-atac`, and `build-seed-reference` on a login node with a
shared `CRETARGET_CACHE_DIR`. Pipeline tasks only read the cache. See
[Reference data](../getting-started/reference-data.md).

## How do I cite it?

Use the Zenodo concept DOI
[10.5281/zenodo.24148024](https://doi.org/10.5281/zenodo.24148024) for the
software in general, or the version DOI when exact reproducibility matters.
Also cite FlashFry, launch-dCas9, and Hartman et al. if you use those features.
See [Citation](../project/citation.md).

## How do I report a bug?

Use the [bug report template](https://github.com/MikeKazz/creTarget/issues/new?template=bug_report.yml)
and include the version, command, `cretarget doctor` output, trace, and smallest
reproducing input. For security issues, follow
[SECURITY.md](https://github.com/MikeKazz/creTarget/blob/main/SECURITY.md).

## Does creTarget send data anywhere?

No. It reads local files and writes local files. The only network access is the
explicit fetch commands and the normal Nextflow resolution of the pipeline from
GitHub. Compute tasks do not upload data.
