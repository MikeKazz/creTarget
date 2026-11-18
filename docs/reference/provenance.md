# Provenance and reproducibility

Every output directory gets two provenance files:

- `run_manifest.json` — code, parameters, command line, and run context.
- `software_versions.yml` — versions reported by each process.

They exist because a guide table without its parameters is not a result. A year
later, "which reference, which release, which filters?" should be answerable
from the output directory alone.

## `run_manifest.json`

```json
{
  "pipeline": {
    "name": "MikeKazz/creTarget",
    "version": "0.5.1",
    "revision": "v0.5.1",
    "commit_id": "30b1ac0...",
    "commit_source": "nextflow",
    "uncommitted_changes": null,
    "home_page": "https://github.com/MikeKazz/creTarget"
  },
  "run": {
    "run_name": "elegant_curie",
    "session_id": "...",
    "start_time": "2026-11-17T...",
    "command_line": "nextflow run MikeKazz/creTarget ...",
    "profile": "slurm,conda",
    "container_engine": "none",
    "nextflow_version": "25.10.3",
    "work_dir": "/scratch/work",
    "launch_dir": "/home/user/project",
    "user": "user"
  },
  "params": {
    "snp_bed": "snps.txt",
    "enhancer_bed": "enhancers.txt",
    "chunk_size": 1000,
    "run_seed_finder": true,
    "...": "..."
  }
}
```

The exact fields are listed in
[Output columns](output-columns.md#run_manifestjson).

### What `commit_source` means

- `nextflow` — Nextflow fetched the repository, so `workflow.commitId` is
  authoritative.
- `git(projectDir)` — the pipeline ran from a local clone. creTarget asks git
  for `HEAD` and records whether the tree was dirty.
- `none` — no git metadata was available. This is a warning sign for
  reproducibility; archive the exact source tree with the results.

`uncommitted_changes` is `true`, `false`, or `null`. A `true` value means the
results came from a working tree that is not fully described by the commit ID.
For a publication, run from a clean checkout of a tagged release.

## `software_versions.yml`

One YAML block per process, deduplicated by process name:

```yaml
"PREP_REGIONS":
    cretarget: 0.5.1
    bedtools: 2.31.1
"FLASHFRY_DISCOVER_SCORE":
    java: 17.0.13
    flashfry: FlashFry-assembly-1.15.jar
"SEED_FINDER":
    cretarget: 0.5.1
    seed_reference: "62700 TSS rows"
```

Use it to record tool versions in a methods section. It does not contain the
container image digest; if you use a container, record the image tag or digest
as well.

## Determinism

The pipeline makes a specific reproducibility claim: **the same inputs and
parameters produce byte-identical `guides.ranked.tsv` across independent runs,
regardless of where the work directory lives.**

That is why chunks are collated by numeric chunk index and not by the absolute
path. `sort: true` would order by the work-directory hash and produce a
different row order between otherwise identical runs. CI checks this by running
the pipeline twice with different work directories and comparing the output.
Chunk size does not change the ranking, because every guide for a region stays
in the same chunk; changing `--chunk_size` is not a reproducibility risk, but it
is still a parameter, and it is recorded in the manifest.

What is **not** claimed:

- Bit-identical output across different FlashFry versions, container images, or
  reference builds.
- Identical output after changing a parameter, even one that looks cosmetic.
- Identical output when the input order changes, because `region_id` order and
  duplicate naming can change.

## Reproducing a run

1. Read `run_manifest.json`: version, commit, command line, profile, parameters.
2. Check out the recorded commit or pull the recorded container tag.
3. Use the same reference files. The manifest records their paths, not their
   checksums; record the Zenodo DOI or checksum alongside it.
4. Use the same Nextflow version and profile.
5. Run with the recorded command line. `-resume` is optional for a fresh
   reproduction, but useful for extending a run.
6. Compare `guides.ranked.tsv` byte-for-byte if the inputs are identical.

## Archiving a result

For a publication or a long-lived library, archive together:

- `guides.ranked.tsv` (and `guides.final.tsv` if used)
- `regions/regions.metadata.tsv`
- `seed_finder/*_raw_hits.csv` if used
- `run_manifest.json`
- `software_versions.yml`
- `pipeline_info/trace.txt`
- the input target files
- the exact reference build identifiers (Zenodo DOI, GENCODE release, GTF
  checksum)
- the container image digest, if used

Compress large tables with `gzip` or `zstd`; the genome-wide table is roughly
50 GB plain and 10 GB gzipped. Do not delete the work directory until the run is
complete and the output is archived if you may need to resume it.

## Versioning and citations

The software version is in `pyproject.toml`, `nextflow.config`, `CITATION.cff`,
the Dockerfile, and `assets/environment.yml`. `./bump_version.sh --check` fails
if they disagree. The release process is in
[Release process](../project/releasing.md).
