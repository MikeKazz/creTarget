# Authors and acknowledgements

## Maintainer

- **Mehrdad Kazemi** ([@MikeKazz](https://github.com/MikeKazz)) — design,
  implementation, documentation, and releases.

## Contributors

Contributors are recorded in the repository history and in release notes. If
your contribution is not visible there, open an issue so it can be corrected.

## Third-party methods and software

creTarget would not exist without the work below. If you use the corresponding
feature, cite the original work as well as creTarget.

- **FlashFry** — target discovery and on/off-target scoring.
  McKenna A, Shendure J. [mckennalab/FlashFry](https://github.com/mckennalab/FlashFry).
- **launch-dCas9** — optional CNN-based enhancer and promoter efficacy scoring.
  [Wancen/launch-dCas9](https://github.com/Wancen/launch-dCas9). Pretrained
  weights are redistributed via Zenodo with attribution.
- **Hartman et al. 2026** — promoter-proximal seed screening method.
  *Systematic identification of seed-driven off-target effects in Perturb-seq
  experiments.* [doi:10.64898/2026.03.27.714658](https://doi.org/10.64898/2026.03.27.714658).
- **seed-finder-lite** — promoter-window reference design and single-pass batch
  architecture used by the seed-screening step.
  [doi:10.5281/zenodo.21599135](https://doi.org/10.5281/zenodo.21599135).
- **Nextflow** — workflow execution and scheduler abstraction.
- **GENCODE** — gene annotation used to build the promoter reference.
- **ENCODE** — cCRE and ATAC-seq resources used in validation.

The full citation metadata is in [CITATION.cff](CITATION.cff).
