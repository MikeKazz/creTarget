# Release process

This page is for maintainers. It describes how a tagged release becomes a GitHub
Release, a PyPI distribution, and a container image, and how to undo a bad
release.

## Versioning

creTarget follows [Semantic Versioning](https://semver.org/):

- **MAJOR** — incompatible changes to output columns, parameters, or the
  meaning of a score.
- **MINOR** — new functionality that is backward compatible.
- **PATCH** — bug fixes and documentation/infrastructure changes.

Tags are `vMAJOR.MINOR.PATCH`. The tag must match the
version in `pyproject.toml`, `nextflow.config`, `CITATION.cff`, the Dockerfile,
and `assets/environment.yml`. `./bump_version.sh --check` verifies this.

## One-time repository setup

These settings are not stored in the repository and must be configured once:

1. **PyPI trusted publishing.** On PyPI, add a pending/trusted publisher for
   project `cretarget`:
   - owner: `MikeKazz`
   - repository: `creTarget`
   - workflow: `release.yml`
   - environment: `pypi`

   The workflow uses OpenID Connect; there is no PyPI token in GitHub secrets.
2. **GitHub environment.** Create an environment named `pypi` in the repository
   settings if you want required reviewers or a deployment branch rule.
3. **Actions permissions.** Repository → Settings → Actions → General:
   allow GitHub Actions to create and approve pull requests only if you use
   Dependabot auto-merge; the release workflow itself needs the permissions
   declared in the workflow file.
4. **GitHub Pages.** Repository → Settings → Pages → Build and deployment →
   Source: **GitHub Actions**. The `docs.yml` workflow deploys the MkDocs site
   on every push to `main`.
5. **Container visibility.** After the first publish, set the GHCR package to
   public if you want anonymous pulls.
6. **Branch protection.** Require the CI checks on `main`. The release workflow
   assumes `main` is the default branch.
7. **Wiki.** The wiki is a separate git repository. See
   [Publishing the wiki](#publishing-the-wiki).

## Cut a release

1. **Update the changelog.** Move the `Unreleased` entries under a new
   `## [X.Y.Z] - YYYY-MM-DD` heading. Keep the format from
   [Keep a Changelog](https://keepachangelog.com/).
2. **Bump the version everywhere.**

   ```bash
   ./bump_version.sh X.Y.Z
   ```

   The script updates all version pins and then verifies they agree. It does
   not edit historical changelog entries or explain why the release exists.
3. **Run the checks.**

   ```bash
   pytest
   ruff check src tests
   mkdocs build --strict
   ./bump_version.sh --check
   ```

4. **Merge to `main`.** Open a PR, let CI pass, and merge. The release tag must
   point at a commit on `main`.
5. **Tag and push.**

   ```bash
   git switch main
   git pull --ff-only
   git tag -a vX.Y.Z -m "creTarget vX.Y.Z"
   git push origin vX.Y.Z
   ```

6. **Watch the release workflow.** It validates the tag/version match, reruns
   the tests, builds the sdist and wheel, checks the distribution, creates the
   GitHub Release with the changelog section, publishes to PyPI, and builds and
   pushes the container image to GHCR.
7. **Verify.**

   ```bash
   pip index versions cretarget
   docker manifest inspect ghcr.io/mikekazz/cretarget:X.Y.Z
   gh release view vX.Y.Z
   ```

8. **Zenodo.** If the repository is connected to Zenodo, publishing the GitHub
   Release mints a version DOI. The concept DOI
   [10.5281/zenodo.24148024](https://doi.org/10.5281/zenodo.24148024) always
   resolves to the latest version. Update `CITATION.cff` only when the version
   DOI or metadata changes.
9. **Bioconda (manual).** There is no automated Bioconda release. Once the PyPI
   sdist exists, a recipe can be submitted to
   [bioconda-recipes](https://github.com/bioconda/bioconda-recipes). Until then,
   do not claim `conda install cretarget` works.

## What the release workflow does

`.github/workflows/release.yml` runs on `v*.*.*` tags:

| Job | What it does |
|---|---|
| `validate` | Checks the tag against every version pin and the changelog, builds the sdist/wheel, runs `twine check`, asserts no model checkpoints are bundled |
| `ci` | Calls the reusable CI workflow: tests, distribution build, lint, docs build, version check, and the Nextflow smoke tests |
| `github-release` | Creates the GitHub Release with `dist/*` and the changelog section |
| `pypi` | Publishes to PyPI using trusted publishing (`environment: pypi`) |
| `container` | Builds and pushes `ghcr.io/mikekazz/cretarget` with semver tags and SBOM/provenance |

If the PyPI trusted publisher is not configured, the `pypi` job fails; the
GitHub Release and container jobs can still succeed. Configure the publisher
before the first release.

## Publishing the wiki

The GitHub wiki is a separate repository at
`https://github.com/MikeKazz/creTarget.wiki.git`. GitHub does not create it
until the first page exists. Create a placeholder `Home` page in the web UI
once, then:

```bash
bash scripts/publish-wiki.sh --dry-run   # show what would change
bash scripts/publish-wiki.sh --yes       # commit and push
```

The script creates a timestamped backup branch in the wiki repository before
pushing, so a bad sync can be reverted. It refuses to run if the wiki has not
been initialized.

## Undo a release

Nothing below is routine. Prefer a fixed patch release over rewriting public
history.

### Delete a GitHub Release and tag

If the release was just created and nothing depends on it:

```bash
gh release delete vX.Y.Z --yes
git push origin :refs/tags/vX.Y.Z
git tag -d vX.Y.Z
```

Deleting a tag does not remove the PyPI files or the container image.

### Yank a PyPI release

PyPI does not allow re-uploading the same version. **Yank** the bad release and
publish a fixed patch:

```bash
# In the PyPI web UI: Manage → Releases → Yank
```

A yanked release is still installable by exact version but is not selected by
default. Never reuse a version number.

### Remove or retag a container image

If nobody has pulled the tag, delete it in the package settings. If it may have
been consumed, publish a fixed patch and leave the old tag alone; deleting a
tag that downstream users pinned breaks reproducibility.

### Revert code

```bash
git revert <commit>
```

Use a normal revert commit; do not force-push `main`. If the release was tagged
on the bad commit, follow the patch-release path.

### Docs and wiki

- **Docs site:** re-run the `docs.yml` workflow on the last good commit, or
  revert the docs commit and push to `main`. The deployment history is in the
  Actions run.
- **Wiki:** the publish script creates a `backup/<timestamp>` branch in the wiki
  repository. Reset `master` to the backup and push, or revert the sync commit.

## Release checklist

- [ ] `CHANGELOG.md` has a dated section for the version.
- [ ] `./bump_version.sh X.Y.Z` ran and verification passed.
- [ ] `pytest`, `ruff check`, and `mkdocs build --strict` pass.
- [ ] CI is green on `main`.
- [ ] Tag `vX.Y.Z` is on a `main` commit and pushed.
- [ ] Release workflow is green.
- [ ] PyPI shows the version.
- [ ] GHCR shows the tag.
- [ ] GitHub Release notes are readable and complete.
- [ ] Zenodo has the new version DOI (if connected).
- [ ] The documentation site deployed.
- [ ] The wiki was synced, if it changed.
