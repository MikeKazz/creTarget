# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.3.0]

Seed screening is now part of creTarget rather than a call out to a separately
installed tool. This is a **breaking change** to the seed-screening flags.

### Removed
- `--seed_finder_dir`, `--seed_finder_python`, `--seed_finder_skip_check` and
  the `cretarget seed-finder-check` subcommand. There is no external install to
  point at, no interpreter to borrow, and no third-party build to validate at
  runtime.
- `pyfaidx` is no longer a dependency. The scan visits every promoter window
  exactly once, so it needs no random access: the reference is read as a plain
  gzip stream. That also removes the 212 MB decompression into a cache
  directory, the `.fai` index, and the BioPython requirement that came with
  reading block-gzipped FASTA.

### Added
- **A native seed scanner** (`cretarget.seed_scan`). Finds every site where a
  guide's PAM-proximal seed matches exactly next to an NGG, within a window of
  an annotated TSS, on either strand. Validated against the tool it replaces:
  on the same inputs it produced **byte-identical hit sets** -- 728 hits on a
  40-window reference, 29,338 hits for 2,000 guides against the full
  promoterome -- with every field of every hit agreeing, at slightly lower
  runtime. It is additionally tested against a brute-force implementation
  written from the biological definition, across several seed lengths and
  random fixtures.
- **`cretarget build-seed-reference`** builds the promoter reference instead of
  downloading one: a pinned GENCODE release plus the hg38 FASTA already
  required, `bedtools getfasta`, and nothing hosted. The window recipe
  reproduces seed-finder-lite's published reference exactly -- +/-2 kb around
  each TSS, merged per chromosome, then padded by 23 bp -- verified to produce
  the identical 49,497 windows on all 24 chromosomes. chrM is excluded, as in
  the original. A `provenance.json` records the annotation release, the GTF's
  SHA-256, the window geometry and the chromosome set, so the promoter
  definition is reproducible from declared inputs rather than being an opaque
  property of a blob.
- The seed screen now runs in CI. It was previously untestable there, because
  it needed a 66 MB external reference; building the reference from a synthetic
  GTF over the CI mini-genome takes under a second. CI asserts that hits appear
  on **both strands**, which is the visible symptom of the one-orientation
  defect described below, plus the table-shape invariants.
- `--seed_finder_all_guides`. By default only guides that pass the hard filters
  are screened -- they are the only ones going into a library, and excluding
  the rest cuts both the seed index and the work. Excluded guides get blank
  seed columns rather than a count of 0, keeping "we looked and found nothing"
  distinguishable from "we did not look".

### Fixed
- **Seed hits on one strand were never reported.** The external tool indexed
  every guide under `reverse_complement(guide[-min_seed_length:])` and used
  that single key for both scan branches, so the branch handling a
  plus-strand protospacer could not match anything. Only sites whose
  protospacer lay on the minus strand were ever found -- roughly half of all
  true sites -- and the tool exited 0 with plausible-looking output. Measured
  against sites taken verbatim from the tool's own reference: 0/40 plus-strand
  recovered, 40/40 minus-strand.

  The native scanner indexes both orientations and is pinned by tests that
  fail if either side regresses. **Any seed-screening results produced before
  this release are incomplete and should be regenerated.**
- `strand` in the hit table now names the strand the **protospacer** lies on,
  matching FlashFry and creTarget's own `orientation` column. The previous
  value was inverted relative to both, so the two columns in one row disagreed.
- The TSS table's header is detected rather than assumed. A headerless file
  previously lost its first gene silently.

### Changed
- The hit table gains `pam_start` and `pam_end`. The PAM sits 3' of the
  protospacer for a plus-strand hit and 5' of it for a minus-strand one, which
  is not something a reader should have to reconstruct.
- `SEED_FINDER`'s memory escalates on retry and its walltime is sized for one
  whole-table task. Guides per task is bounded by the seed index, roughly 1 GB
  per million guides.

## [0.2.6]

### Changed
- **seed-finder-lite now runs once over the whole guide table, not once per
  chunk.** Its batch mode scans the promoterome a single time and then matches
  every guide by dictionary lookup, so cost is essentially independent of guide
  count: measured at ~25 s for 20 guides and ~25 s for 2,000. Invoking it per
  chunk paid that fixed scan on every chunk -- at `--chunk_size 1000` a
  genome-wide run would have spent on the order of a month re-reading the same
  212 MB of sequence. The step now sits after collation.
  - `guides.ranked.tsv` is published on every run. `--run_seed_finder` adds
    `guides.final.tsv` rather than replacing the ranked table.
  - `--seed_finder_max_guides N` splits the work into tasks of at most N
    guides, for cases where one task cannot finish inside its walltime. Each
    shard repeats the promoterome scan, so a low shard count is the point.
    Verified to produce output byte-identical to the single-task path.
- **The step streams the guide table instead of loading it.** It previously
  read the table with `pandas.read_csv`, which was fine for one chunk and
  impossible for a collated genome-wide table of tens of GB. Both passes now
  stream, and memory tracks the number of guides that actually have a hit
  rather than the number of guides.
- `guide_id` now always ends in the row index
  (`<region>_<start>_<end>_<orientation>_<n>`). It previously appended a
  `_dupN` suffix only on collision, which needed a global de-duplication pass
  over the whole table. This changes the `guide_id` values in
  `seed_finder/*_raw_hits.csv`.
- `SEED_FINDER`'s default walltime is raised to 12 h and its CPU request
  lowered to 1, reflecting one whole-table task rather than many small ones.

### Added
- **An install check for seed-finder-lite, run before the real scan.** It asks
  the installed tool to recover two real Cas9 sites taken from its own promoter
  reference, one per strand.

  This exists because of a defect found in seed-finder-lite itself: in the
  batch scanner, guides are indexed only under
  `reverse_complement(guide[-min_seed_length:])`, and the plus-strand-PAM loop
  looks up that same key while scoring against `guide_rc`. That loop therefore
  matches nothing, and the tool reports only hits whose protospacer lies on the
  minus strand -- exiting 0, with plausible-looking output, while missing
  roughly half of all true hits. Measured against sites taken verbatim from the
  tool's own reference: 0/40 plus-strand sites recovered, 40/40 minus-strand.
  No amount of checking creTarget's own code would have caught it, and the only
  visible symptom is that every reported hit carries the same strand.

  If the check fails, the run stops and the error names the cause.
  `--seed_finder_skip_check` bypasses it; `cretarget seed-finder-check
  --seed-finder-dir ...` runs it standalone.
- `seed-finder-lite`'s version is recorded in `software_versions.yml`, read
  from a `VERSION` file at the root of its directory if present.
- `pyfaidx` added to `assets/environment.yml`. seed-finder-lite imports it to
  read its promoter reference, so `--run_seed_finder` previously failed under
  `-profile conda` and the container profiles unless `--seed_finder_python`
  pointed at an interpreter outside the environment. (biopython is not needed:
  the bundled reference is plain gzip, which seed-finder-lite decompresses
  itself.)

## [0.2.5]

### Added
- **Startup banner.** `cretarget` and `nextflow run` both open with the
  project wordmark and a schematic of a guide bound to DNA. In the CLI it is
  written to stderr, never stdout, so redirecting a step's output to a file
  stays safe; and it is printed only when stderr is a terminal, because a
  genome-wide run has thousands of tasks and a banner in every
  `.command.err` is noise rather than branding. `CRETARGET_NO_BANNER=1`
  suppresses it outright, and colour follows the
  [NO_COLOR](https://no-color.org) convention — disabled by the variable's
  presence, whatever its value.

### Fixed
- **`conf/example-singlenode.config` built an invalid SLURM request when no
  GPU was asked for.** `clusterOptions = "--gres=${params.gpu_gres}"` was
  unguarded, so running `--run_launch_dcas9` without `--gpu_gres` submitted
  the literal string `--gres=null`, which SLURM rejects. The generic `slurm`
  profile has always guarded this; the guard was lost when the site-specific
  config was genericized for 0.2.4. Without a GPU request the CNN now falls
  back to CPU, which is slower but works.

## [0.2.4]

No functional changes. The pipeline, its results and its dependencies are
identical to 0.2.3; this release makes the repository usable by someone other
than its author.

### Changed
- **Site-specific material replaced with generic equivalents.** The repository
  carried the configuration, paths and setup notes of the single cluster it was
  developed on. A reader could not tell which numbers were general and which
  were artefacts of one machine, and the documented commands referenced
  directories that exist nowhere else.
  - `conf/boltzmann.config` -> `conf/example-singlenode.config`, and the
    profile that includes it is now `singlenode`. **This renames a profile:**
    `-profile boltzmann` no longer exists. The configuration itself is
    unchanged -- every measured value and the reasoning behind it is retained,
    with only the host and path identifiers removed.
  - `assets/inspect_boltzmann.sh` -> `assets/inspect_cluster.sh`. The two
    hardcoded directories became optional arguments 2 and 3, defaulting to the
    working directory, so the survey runs anywhere without editing.
  - `assets/cretarget_env.sh` -> `assets/site_env.sh.example`, now a template
    of placeholders rather than one site's live values.
  - `assets/run_cretarget.sbatch` rewritten around a single `CONFIGURE` block
    at the top, with preflight checks that fail on a missing reference or an
    unindexed FASTA before submitting hours of work rather than during it.
  - `docs/boltzmann-setup.md` -> `docs/cluster-setup.md`, rewritten as a
    six-step first install for an arbitrary cluster: survey, install,
    reference data, smoke test, site profile, scale up. Its troubleshooting
    section covers the failures actually encountered during development --
    process substitution as an input, a stale `NXF_SINGULARITY_CACHEDIR`,
    `conda activate` in a non-interactive batch shell, requests escalating
    past a partition's limits, and output ordering.
  - Removed `docs/release-checklist.md`, which was a private worklist.

### Fixed
- `bump_version.sh` did not update the pipeline revision pinned in
  `assets/run_cretarget.sbatch`, so the driver template would have kept
  fetching the previous release after a bump. It is now updated and verified,
  and the script additionally warns about any remaining mention of the old
  version so a missed pin cannot pass silently again.

## [0.2.3]

### Fixed
- **The pipeline could not run on Nextflow 26.** Nextflow 26 stopped coercing
  command-line parameters to numbers, so `--chunk_size 1000` arrives as a
  String while the same value from a config file or profile arrives as an
  Integer. `splitFasta(by:)` rejects the String with "Value don't match: class
  java.lang.Integer" -- phrasing that reads as if Integer were the problem
  rather than what is expected. Since 26.x is what a fresh `curl
  get.nextflow.io` installs, this broke every new user. The value is now cast
  explicitly. CI gained Nextflow 26.10.5 and now passes `--chunk_size` on the
  command line: supplying it only through the test profile meant the String
  path was never exercised, which is why CI stayed green through this bug.
- **`CRETARGET_CACHE_DIR` was ignored when `[cache].dir` was configured.** A
  configured value silently won, so `CRETARGET_CACHE_DIR=/tmp/scratch cretarget
  setup-reference` downloaded 3.2 GB to the configured location instead of the
  requested one, without a word. The environment variable now takes precedence,
  which is the conventional behaviour for a per-invocation override.
- Removed a warning that fired on every run without a container profile. It was
  correct when `container_image` had no default, but once the parameter gained
  one the condition was always true.

### Added
- Unit tests for `summarize_hits` and `build_guide_ids` (the seed-finder
  aggregation, previously untested -- an end-to-end run only exercises it when
  seed-finder actually finds hits) and for `coords-to-bed`, which had no
  coverage at all. 86 tests to 115.

## [0.2.1]

Archived on Zenodo: concept DOI [10.5281/zenodo.24148024](https://doi.org/10.5281/zenodo.24148024),
this version [10.5281/zenodo.24148025](https://doi.org/10.5281/zenodo.24148025).

Release engineering and fixes found while validating the first public release.
No change to pipeline results.

### Added
- `bump_version.sh`: sets the release version in all nine places it appears
  across six files and verifies they agree. Version drift between
  `pyproject.toml`, the container tag and the conda environment's install URL
  is otherwise an easy mistake to ship.
- Published container image `ghcr.io/mikekazz/cretarget`, now the default for
  the `docker`/`singularity`/`apptainer` profiles, so `-profile singularity`
  needs no extra flags.
- `mamba` and `micromamba` profiles: the same pinned environment, solved in
  seconds rather than the many minutes conda's classic solver can take.
  Implemented as profiles rather than a parameter because a profile that reads
  `params.X` sees whatever X was when the profile block was parsed.

### Fixed
- **`-profile conda` could not build.** `assets/environment.yml` installed
  creTarget with `git+https://`, which needs the git binary, network AND
  working git credentials on whatever node builds the environment; pip runs the
  clone with `GIT_TERMINAL_PROMPT=0`, so a stale credential helper turned a
  public repo into "Authentication failed / Repository not found". It now
  installs from the release tarball over plain HTTPS -- no git, no auth, and
  immutable because it is pinned to a tag.
- **Unpinned conda environment installed untested versions.** Left unpinned it
  resolved to python 3.15 / pandas 3.1 / numpy 2.5, none of which this tool has
  been tested against, and pandas 3.x has breaking API changes. Versions are
  now pinned to the ones the test suite and published benchmark ran on, which
  also cuts the solve from an hour (or an apparent hang) to seconds. A 60-minute
  `conda.createTimeout` means a slow solve now fails with a clear error instead
  of hanging silently.
- CI could not install Nextflow: the release asset was renamed mid-series
  (`nextflow-<v>-all` through 24.04, `-dist` from 25.10), so a hardcoded name
  404s for half the matrix. It now tries one and falls back to the other.
- A test compared a CSV-round-tripped float with `==`, which passed on
  pandas 3.x and failed on the pandas 2.x that CI resolves for Python 3.10.
  Now compared with a tolerance, and `composite_score_max_possible` is rounded
  at the source so the published column cannot carry float noise.
- The run manifest recorded `commit_id: n/a` whenever the pipeline ran from a
  local clone, because `workflow.commitId` is only set when Nextflow manages
  the repository itself. It now falls back to `git rev-parse` on `projectDir`
  and also records whether the working tree was dirty -- a commit hash is
  misleading if the run used uncommitted edits.
- An unusable input path (a deleted file, a broken symlink, or a bash process
  substitution such as `<(...)`, which Nextflow stages by file-descriptor
  basename and so cannot stage at all) produced a Python traceback instead of
  a readable message.
- `PREP_REGIONS` memory raised to 12 GB. An earlier in-process measurement said
  2.22 GB, but it instrumented only the Python functions and missed the
  `bedtools getfasta` child process; the real figure at full scale is 5.20 GB.

## [0.2.0]

Execution moved from a hand-written SLURM layer to Nextflow. This is a
breaking change to how the pipeline is invoked.

### Added
- `conf/site-template.config`: a commented blank site profile that names the
  specific scheduler queries behind each setting (does the scheduler account
  for memory, does it enforce it, is the node shared with work it cannot see),
  so a new cluster can be configured without reverse-engineering the example.
- `assets/inspect_cluster.sh`: read-only cluster survey reporting everything
  a site profile needs -- scheduler config, conda layout, filesystem types and
  locking, GPU compute capabilities, and pre-existing reference data.
- `assets/site_env.sh.example`: project-local environment file, sourced per
  session. Uses only `CRETARGET_*`/`NXF_*` names and keeps the tool's config
  inside the project, so nothing leaks into a shared `~/.bashrc`.
- Site profile support: `--queue` and `--gpu_gres` parameters wired into the
  SLURM, PBS, SGE, and LSF profiles, plus a worked single-node SLURM example
  (`conf/example-singlenode.config`), a driver submission template
  (`assets/run_cretarget.sbatch`), and `docs/cluster-setup.md`.
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
- **The collated output is now reproducible.** `collectFile(sort: true)`
  ordered the per-chunk tables by absolute path, which contains each task's
  random work-directory hash, so two identical runs concatenated the chunks in
  different orders and produced byte-different output files. No rows were lost
  and every row is self-labelled, but "same input, same output" is a claim the
  pipeline has to be able to make. Chunks are now ordered by their numeric
  index (zero-padded, so chunk 2 precedes chunk 10), and CI fails if two
  independent runs of the same input differ.
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
