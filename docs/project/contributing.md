# Contributing

The full contributor guide is
[CONTRIBUTING.md](https://github.com/MikeKazz/creTarget/blob/main/CONTRIBUTING.md).
This page is the short version for the documentation site.

## Set up a development environment

```bash
git clone https://github.com/MikeKazz/creTarget.git
cd creTarget
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -e ".[dev]"
pre-commit install --install-hooks
```

Requirements: Python 3.10+, Java 17+, `bedtools`, `samtools`, and Git. The
Python tests do not need FlashFry or a real genome.

## Run the checks

```bash
pytest
ruff check src tests
bash bump_version.sh --check
mkdocs build --strict
```

The Nextflow end-to-end test needs Nextflow, Java, FlashFry, and a tiny indexed
reference. CI builds a synthetic mini-genome and a real FlashFry index; see
`.github/workflows/ci.yml`.

## Documentation

The canonical documentation is in `docs/` and is built with MkDocs Material. The
GitHub wiki in `wiki/` is a concise mirror.

```bash
pip install -r docs/requirements.txt
mkdocs serve          # http://127.0.0.1:8000
mkdocs build --strict # what CI runs
```

When you change a parameter, an output column, or a failure mode, update the
relevant reference page in the same PR.

## Scientific changes

For scoring, filters, seed definitions, or reference-building changes:

- State the method and cite the source in the module docstring and the PR.
- Add a regression test with a hand-checkable example.
- Explain whether existing results change; if they do, add a **Breaking** or
  **Changed** entry to `CHANGELOG.md`.
- Do not silently change a default. Defaults are part of the method.

## Pull requests

- Keep commits focused and use an imperative subject.
- Include tests for new behavior and documentation for new user surface.
- Do not commit reference genomes, FlashFry indexes, model checkpoints, or
  large result tables.
- CI must pass before review.

The pull request template contains the full checklist.

## Code of conduct and security

By participating, you agree to the
[Code of Conduct](https://github.com/MikeKazz/creTarget/blob/main/CODE_OF_CONDUCT.md).
Report security issues privately according to
[SECURITY.md](https://github.com/MikeKazz/creTarget/blob/main/SECURITY.md); do
not open a public issue.
