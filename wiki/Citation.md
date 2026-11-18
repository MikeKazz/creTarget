# Citation

creTarget is archived on Zenodo. Cite the software and the methods you used.

## Software

Use the **concept DOI** when readers should get the current version:

> Kazemi M. *creTarget: CRISPRi guide RNA design for the noncoding genome.*
> Zenodo. [doi:10.5281/zenodo.24148024](https://doi.org/10.5281/zenodo.24148024).

```bibtex
@software{cretarget,
  author       = {Kazemi, Mehrdad},
  title        = {creTarget: CRISPRi guide RNA design for the noncoding genome},
  year         = {2026},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.24148024},
  url          = {https://github.com/MikeKazz/creTarget},
  note         = {Concept DOI; resolves to the latest release}
}
```

Use the **version DOI** when exact reproducibility matters. The concept DOI
always resolves to the latest release; the version DOI pins one snapshot. The
version DOI for v0.2.1 is
[10.5281/zenodo.24148025](https://doi.org/10.5281/zenodo.24148025); later
releases mint their own version DOIs.

`CITATION.cff` is the machine-readable metadata behind GitHub's "Cite this
repository" button.

## FlashFry (always)

> McKenna A, Shendure J. *FlashFry: a fast and flexible tool for large-scale
> CRISPR target design.* [github.com/mckennalab/FlashFry](https://github.com/mckennalab/FlashFry).

Check the FlashFry repository for the current paper citation.

## launch-dCas9 (`--run_launch_dcas9`)

> Zhou W. *launch-dCas9.* [github.com/Wancen/launch-dCas9](https://github.com/Wancen/launch-dCas9).

The pretrained weights are redistributed via Zenodo record
[10.5281/zenodo.22729308](https://doi.org/10.5281/zenodo.22729308) with
attribution; they remain the work of the launch-dCas9 authors.

## Seed screening (`--run_seed_finder`)

> Hartman A, Blair JD, Nguyen TP, Dyson K, Bradu A, Takacsi-Nagy O,
> Santostefano K, Boade T, Bolaños M, Zhu R, Dann E, Marson A, Gitler AD,
> Satija R, Satpathy AT, Roth TL. *Systematic identification of seed-driven
> off-target effects in Perturb-seq experiments.* 2026.
> [doi:10.64898/2026.03.27.714658](https://doi.org/10.64898/2026.03.27.714658).

The promoter-window reference design and single-pass batch architecture also
come from seed-finder-lite:

> Kazemi M. *seed-finder-lite.* Zenodo.
> [doi:10.5281/zenodo.21599135](https://doi.org/10.5281/zenodo.21599135).

`--run_seed_finder` prints the Hartman citation, and the reference's
`provenance.json` records the GENCODE release and GTF checksum.

## ATAC-seq peaks

If you use the bundled peaks, cite the original data source and ENCODE in
addition to creTarget. The Zenodo record
[10.5281/zenodo.21268805](https://doi.org/10.5281/zenodo.21268805) is a
redistribution, not a new measurement.

## Reference data

The hg38 FASTA and FlashFry index are distributed via
[10.5281/zenodo.21115562](https://doi.org/10.5281/zenodo.21115562). GENCODE
annotations have their own terms; see
[gencodegenes.org](https://www.gencodegenes.org/).

## Suggested methods sentence

> CRISPRi guides were designed with creTarget vX.Y.Z
> (doi:10.5281/zenodo.24148024) using FlashFry 1.15 for target discovery and
> scoring. [Optional: Promoter-proximal seed screening followed Hartman et al.
> 2026 (doi:10.64898/2026.03.27.714658), using GENCODE release 47 and a ±2 kb
> TSS window.] [Optional: On-target efficacy was scored with launch-dCas9.]
> The exact parameters and tool versions are recorded in `run_manifest.json`
> and `software_versions.yml`.

Replace the bracketed text with the features you actually used.

Full guidance: [Citation](https://mikekazz.github.io/creTarget/project/citation/)
