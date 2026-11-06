# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - Unreleased

Execution moved from a hand-written SLURM layer to Nextflow. This is a
breaking change to how the pipeline is invoked.

### Added
- `conf/site-template.config`: a commented blank site profile that names the
  specific scheduler queries behind each setting (does the scheduler account
  for memory, does it enforce it, is the node shared with work it cannot see),
  so a new cluster can be configured without reverse-engineering the example.
- `assets/inspect_boltzmann.sh`: read-only cluster survey reporting everything
  a site profile needs -- scheduler config, conda layout, filesystem types and
  locking, GPU compute capabilities, and pre-existing reference data.
- `assets/cretarget_env.sh`: project-local environment file, sourced per
  session. Uses only `CRETARGET_*`/`NXF_*` names and keeps the tool's config
  inside the project, so nothing leaks into a shared `~/.bashrc`.
- Site profile support: `--queue` and `--gpu_gres` parameters wired into the
  SLURM, PBS, SGE, and LSF profiles, plus a worked single-node SLURM example
  (`conf/boltzmann.config`), a driver submission template
  (`assets/run_cretarget.sbatch`), and `docs/boltzmann-setup.md`.
- `--gpu_count` for round-robin GPU pinning under the local executor, where no
  scheduler sets `CUDA_VISIBLE_DEVICES`. Under a scheduler the pipeline leaves
  that variable alone, so it never claims a device allocated to someone else.
- Nextflow pipeline (`nextflow run MikeKazz/creTarget`) with profiles for
  local, SLURM, PBS/PBS Pro, SGE, and LSF execution, plus conda, Docker,
  Singularity, and Apptainer software provisioning.
- Automatic retry with escalated memory and time for tasks killed by OOM or
  walltime; `-resume` to continue an interrupted run.
- `cretarget doctor` — preflight checks for external tools, configuration,
  cached assets, and filesystem writability.
- `cretarget fetch-models` and `cretarget fetch-atac` — explicit, checksum
  verified asset downloads designed for clusters whose compute nodes have no
  outbound network access.
- `run_manifest.json` and `software_versions.yml` written to every output
  directory, recording pipeline version, commit, parameters, command line,
  and tool versions.
- `composite_score_normalized` column, so runs with and without launch-dCas9
  scoring are directly comparable.
- `--version` flag.
- Test suite (86 tests) and GitHub Actions CI, including an end-to-end
  Nextflow run against a real FlashFry index.

### Changed
- **Package renamed** from `crispri_guide_designer` to `cretarget`. The
  distribution is now `cretarget` on PyPI.
- Annotation and ranking now run per chunk rather than over a single
  concatenated table, so memory is flat with respect to target count. A
  genome-wide run uses the same per-task memory as a small one.
- FlashFry JVM heap is derived from the task's actual memory allocation
  instead of hardcoded `-Xmx14g`/`-Xmx7g` values, so a retry at higher memory
  automatically gets a larger heap.
- launch-dCas9 spacer deduplication is per chunk rather than global,
  removing a cross-chunk synchronization barrier. Results are unchanged.
- Region preparation handles SNP and enhancer inputs in a single pass, so
  linked-region flagging sees every target rather than one type at a time.
- `cretarget configure --set` rewrites the config by parsing and re-emitting
  it, instead of replacing lines as text.

### Fixed
- **Declared Nextflow requirement corrected from `>=23.04.0` to `>=24.04.0`.**
  The pipeline uses `process.resourceLimits`, introduced in 24.04; on older
  versions the resource ceilings were silently ignored, so a retry could
  request more than the partition allows and be rejected outright. CI tested
  23.04, which could not have worked. The matrix is now 24.04.4, 25.10.3 and
  latest-stable.
- **The container profiles no longer silently run uncontainerised.** They
  enabled a container engine but no process declared an image, so tasks ran
  against whatever tools were on the host while the run reported a container
  engine and appeared reproducible. They now set `process.container` from
  `--container_image`, and the pipeline refuses to start with an engine
  enabled and no image.
- **`-profile conda` now builds.** `assets/environment.yml` pip-installed
  `cretarget` from PyPI, where it is not yet published, so the profile failed
  for everyone. It installs from the GitHub repository until the PyPI release.
- Configured paths containing `~` or `$VARS` are now expanded. Previously a
  literal tilde was passed to bedtools, which broke the documented quickstart.
- Duplicate region names are made unique with a `_dupN` suffix instead of
  silently colliding and then failing the metadata join.
- Input contigs are validated against the reference `.fai` before use, with a
  message identifying an Ensembl-vs-UCSC naming mismatch by name.
- Coordinates containing thousands separators (`chr1:50,000-50,300`) now parse.
- Region IDs containing underscores or colons now round-trip correctly through
  FASTA headers.
- Exclusion-reason assignment and ATAC summit distance are vectorized; both
  are verified to produce output identical to the previous implementations.
- Tarball extraction during reference setup guards against path traversal.

### Removed
- The built-in SLURM executor, chunking layer, `run_standalone.sh`,
  `local_paths.env`, and `bin/` sbatch templates. Nextflow replaces all of it.
  This removes several latent failure modes, including a job-ID parsing bug
  that reported every successful task as failed on multi-cluster SLURM sites,
  and a verification path that did the same wherever SLURM accounting was
  unavailable.
- launch-dCas9 model checkpoints are no longer bundled in the package. They
  are third-party weights and are now fetched from Zenodo on demand, with
  attribution. CI fails the build if a `.pth` reappears in the distribution.
- The `[slurm]` section of `config.toml`; scheduling is Nextflow's concern.

## [0.1.0]

Initial release: SLURM-array and local-parallel execution, FlashFry
integration, ATAC annotation, launch-dCas9 and seed-finder-lite support.
