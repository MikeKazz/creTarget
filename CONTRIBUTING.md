# Contributing to creTarget

Thanks for taking the time to contribute. creTarget is a scientific tool, so
contributions are judged on two things at once: whether the code is sound, and
whether the result it produces is interpretable and reproducible.

This guide covers the development workflow. By participating, you agree to the
[Code of Conduct](CODE_OF_CONDUCT.md).

## Ways to contribute

- **Report a bug.** Use the [bug report template](https://github.com/MikeKazz/creTarget/issues/new?template=bug_report.yml).
  Include the version, the command line, `cretarget doctor` output, and the
  relevant part of `run_manifest.json`.
- **Improve the documentation.** Documentation fixes are first-class
  contributions. If a page is wrong, unclear, or missing a failure mode, open a
  PR or an issue.
- **Propose a feature.** Open a feature request. For changes to scoring,
  filtering, or biological assumptions, describe the evidence and the expected
  effect on existing output.
- **Add a step or an execution profile.** See [Extending the pipeline](#extending-the-pipeline).
- **Report a security issue.** Please follow [SECURITY.md](SECURITY.md) instead
  of opening a public issue.

## Development setup

Requirements: Python 3.10+, Java 17+, `bedtools`, `samtools`, and Git.

```bash
git clone https://github.com/MikeKazz/creTarget.git
cd creTarget
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -e ".[dev,docs]"
pre-commit install --install-hooks
```

The Python test suite does not need FlashFry or a real reference genome. Tests
that need `bedtools` or `samtools` are skipped automatically when those tools
are not on `PATH`.

```bash
pytest
ruff check src tests
bash bump_version.sh --check
mkdocs build --strict
shellcheck --severity=warning assets/*.sh bump_version.sh scripts/*.sh
```

The Nextflow end-to-end test needs Nextflow, Java, FlashFry, and a tiny indexed
reference. CI builds a synthetic mini-genome and a real FlashFry index for this;
see [.github/workflows/ci.yml](.github/workflows/ci.yml) for the exact recipe.

## Documentation

The canonical documentation is in `docs/` and is built with MkDocs Material.
The GitHub wiki in `wiki/` is a concise mirror for users who prefer it.

```bash
pip install -r docs/requirements.txt
mkdocs serve          # http://127.0.0.1:8000
mkdocs build --strict # what CI runs
```

When you change a parameter, an output column, or a failure mode, update the
relevant page under `docs/reference/` or `docs/guides/` in the same PR.

## Extending the pipeline

A new processing step normally touches these places:

1. `src/cretarget/steps/<name>.py` — a `run(args)` function and
   `add_arguments(parser)` for the step-level CLI.
2. `src/cretarget/cli.py` — register the subcommand.
3. `modules/local/<name>.nf` — a Nextflow process that calls the CLI and writes
   `versions.yml`.
4. `workflows/cretarget.nf` — wire the process into the DAG.
5. `conf/base.config` and `conf/modules.config` — resources and publishing
   rules.
6. `tests/` — unit and integration coverage, including the awkward cases.
7. `docs/reference/parameters.md`, `docs/guides/outputs.md`, and
   `CHANGELOG.md` — user-facing documentation.

Keep per-chunk steps independent of other chunks. The flat-memory design depends
on it: every guide for a `region_id` is processed in the same task, so ranking
is correct per chunk and no global merge is needed.

## Scientific changes

For changes to scoring, filters, seed definitions, or reference-building logic:

- State the method and cite the source in the module docstring and in the PR.
- Add a regression test with a hand-checkable example.
- Explain whether existing results change. If they do, add a **Breaking** or
  **Changed** entry to `CHANGELOG.md`.
- Do not silently change a default. Defaults are part of the method.

## Commits and pull requests

- Keep commits focused and use an imperative subject, for example
  `fix: validate contig names before bedtools`.
- Conventional prefixes (`feat:`, `fix:`, `docs:`, `test:`, `ci:`, `build:`,
  `chore:`) are welcome but not required.
- Rebase or squash as needed; the merged history should be readable.
- A PR should include tests for new behavior and documentation for new user
  surface. CI must pass before review.
- Do not commit reference genomes, FlashFry indexes, model checkpoints, or
  large result tables. `.gitignore` covers the usual names; check `git status`
  before committing.

The pull request template contains the full checklist.

## Release process

Maintainers should follow [docs/project/releasing.md](docs/project/releasing.md).
In short: update `CHANGELOG.md`, run `./bump_version.sh <version>`, merge to
`main`, tag `v<version>`, and push the tag. The release workflow builds and
publishes the artifacts.

## Questions

Use the [question issue template](https://github.com/MikeKazz/creTarget/issues/new?template=question.yml).
For bugs, include enough information for someone else to reproduce the failure.
"Does not work" is not reproducible; a command, a version, and a log are.
