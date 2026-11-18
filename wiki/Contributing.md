# Contributing

Contributions are welcome. The full guide is
[CONTRIBUTING.md](https://github.com/MikeKazz/creTarget/blob/main/CONTRIBUTING.md).
This page is the short version.

## Development setup

```bash
git clone https://github.com/MikeKazz/creTarget.git
cd creTarget
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -e ".[dev,docs]"
pre-commit install --install-hooks
```

Requirements: Python 3.10+, Java 17+, `bedtools`, `samtools`, and Git. The
Python tests do not need FlashFry or a real genome.

## Checks

```bash
pytest
ruff check src tests
bash bump_version.sh --check
mkdocs build --strict
shellcheck assets/*.sh bump_version.sh scripts/*.sh
```

The Nextflow end-to-end test needs Nextflow, Java, FlashFry, and a tiny indexed
reference. CI builds a synthetic mini-genome and a real FlashFry index.

## Documentation

- The canonical site is `docs/`, built with MkDocs Material.
- The GitHub wiki is `wiki/`, a concise mirror.

```bash
mkdocs serve          # http://127.0.0.1:8000
mkdocs build --strict # what CI runs
```

When you change a parameter, an output column, or a failure mode, update the
relevant reference page in the same pull request.

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

## Releasing

Maintainers should follow the
[release process](https://mikekazz.github.io/creTarget/project/releasing/).
Version bumps and checks use:

```bash
./bump_version.sh X.Y.Z
./bump_version.sh --check
```

Pushing a `v*.*.*` tag triggers the release workflow, which creates the GitHub
Release, publishes to PyPI through trusted publishing, and builds and pushes the
container image.
