# Roadmap

This is a list of candidate work, not a set of commitments. Priorities are
driven by user reports, reproducibility needs, and available time.

## Near term

- Publish the Python package to PyPI and add a Bioconda recipe.
- Keep the documentation site and GitHub wiki in sync from one source.
- Add a small, downloadable validation dataset so a new user can reproduce the
  documented output without fetching a full genome.
- Improve error messages for the most common cluster failures.

## Medium term

- Validate references other than hg38 when a matching FlashFry index is
  available.
- Add optional ranking presets, keeping the current weights as the default and
  documenting every change.
- Expand the ATAC-seq presets beyond the current hematopoietic panel.
- Provide a Nextflow `-stub-run` path for cluster plumbing tests.

## Not planned

- Bundling the reference genome, FlashFry index, or model checkpoints in the
  Python package or container image. They are large and have their own terms.
- Replacing FlashFry's off-target search with an unvalidated in-house scorer.
- A hosted web service. creTarget is designed to run where the data lives.

Ideas and corrections are welcome in the issue tracker.
