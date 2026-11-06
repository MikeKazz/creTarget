# Release checklist

Ordered by what blocks what. Items 1–4 block a working install; 5–8 block the
paper's claims; 9–12 are the release itself.

---

## 1. Push v0.2.0 to GitHub — BLOCKS EVERYTHING

The public repo currently holds **v0.1.0**: the old SLURM-executor code, no
`main.nf`, no `nextflow.config`, and both `.pth` checkpoints tracked in git.
Until the rewrite is pushed, nothing else on this list can proceed.

From your laptop (or wherever you downloaded the release zip):

```bash
# 1. Get the zip onto the cluster
scp cretarget-0.2.0.zip mkazemi@boltzmann:/home/jmorrislab1/mkazemi/
```

Then on boltzmann:

Then on boltzmann. This clone becomes your working clone too, so there is only
one copy to keep track of:

```bash
export CODE_DIR=/home/jmorrislab1/mkazemi/___ProJECTS___/creTarget
mkdir -p "$CODE_DIR"

cd /home/jmorrislab1/mkazemi
unzip -q cretarget-0.2.0.zip          # -> ./cretarget-0.2.0/

# 2. Clone the existing repo (keeps history and the remote)
cd "$CODE_DIR"
git clone https://github.com/MikeKazz/creTarget.git
cd creTarget

# 3. Remove the old tracked files, keeping .git
git rm -rq .

# 4. Lay down v0.2.0
cp -r /home/jmorrislab1/mkazemi/cretarget-0.2.0/. .

# 5. Sanity-check before committing
ls main.nf nextflow.config conf/ modules/local/     # all must exist
git status --short | head -20
git ls-files --cached | grep -c '\.pth' || true     # must print 0

# 6. Commit and push
git add -A
git commit -m "Replace SLURM executor with Nextflow pipeline (v0.2.0)

Execution moves from a hand-written SLURM layer to Nextflow, with profiles
for local/SLURM/PBS/SGE/LSF. Annotation and ranking now run per chunk, so
memory is flat with respect to target count. Model checkpoints are fetched
from Zenodo rather than bundled. Adds a test suite and CI.

See CHANGELOG.md for the full list."
git push origin main
```

Step 2 of `docs/boltzmann-setup.md` ("clone the pipeline") is now already
done — `$CODE_DIR/creTarget` is your working clone. Continue from step 3.

Then verify from a clean checkout:

```bash
cd /tmp && rm -rf verify && git clone --depth 1 https://github.com/MikeKazz/creTarget.git verify
ls verify/main.nf && grep '^version' verify/pyproject.toml
```

### 1b. Optional: purge the checkpoints from git history

The two `.pth` files remain in the repo's *history* even after the commit
above removes them from the tree. They are third-party weights
(launch-dCas9), so the clean thing is to remove them entirely. The repo was
public only briefly, so a force push is low-risk here — but it rewrites every
commit hash, so do it before anyone clones.

```bash
pip install git-filter-repo
cd "$CODE_DIR/creTarget"
git filter-repo --path src/crispri_guide_designer/checkpoints --invert-paths --force
git push --force origin main
```

Judgment call: you host those weights on Zenodo yourself, so you clearly
consider redistribution-with-attribution permitted, which makes this tidiness
rather than a legal fix. Do it if it is cheap; do not let it block the paper.

---

## 2. Verify the install on boltzmann

`docs/boltzmann-setup.md`, steps 1–7, ending in a clean `cretarget doctor`.

**Why it matters:** if `doctor` is not green, every later failure is ambiguous
— you cannot tell a pipeline bug from a missing dependency.

---

## 3. Smoke test on real reference data

`docs/boltzmann-setup.md`, step 8. ~20 regions, real hg38, real FlashFry index.

**Why it matters:** every test so far has used a synthetic 320 kb mini-genome.
Real hg38 has repeats, blacklist regions, alt contigs and homopolymers that a
random sequence does not. This is the first contact between the code and the
actual data, and it costs minutes instead of hours.

---

## 4. Test the two optional paths

`docs/boltzmann-setup.md`, step 11. launch-dCas9 and seed-finder-lite.

**Why it matters:** these are the **only code paths never executed end to
end** — neither external tool was available in development. They are also the
two features most likely to be called out as the tool's novel contributions,
so a reviewer may well try them. Finding a bug here after submission is much
worse than finding it now.

---

## 5. Benchmark and re-derive the performance claims

`docs/boltzmann-setup.md`, steps 9–10. Run 10k regions, then 20k, confirm
linear scaling, then extrapolate.

**Why it matters:** the old README claimed genome-wide capability based on the
pre-Nextflow architecture and a batching script that was never committed. That
claim is currently unsupported, and I removed it rather than carry it forward.
If "runs on millions of enhancers" is a headline claim in the paper, you need
a measured number behind it, **and you must state the hardware** — a single
20-core node is a very different setting from a multi-node cluster, and
reviewers will assume the latter unless told otherwise.

Record: wall time, peak RSS per process, total CPU-hours, region count,
hardware. Put them in the paper's methods and in the README.

---

## 6. Decide the genome-wide output format

A 2.35M-cCRE run writes a single **~63 GB TSV** (~132M rows at ~477
bytes/row, measured). Options: gzip the collated output (~10–15 GB), split per
chromosome, or write Parquet.

**Why it matters:** this is the artifact a user actually consumes. A 63 GB TSV
is technically usable and practically annoying, and it is a pipeline change
rather than a flag, so it has to be decided before the definitive run rather
than after. It is also the kind of thing a reviewer flags as a usability
concern in a "resource" paper.

---

## 7. Validate against a known-good result

Pick a published CRISPRi screen with a known-good guide set, run creTarget
over the same targets, and show the top-ranked guides recover the ones that
worked.

**Why it matters:** this is the single biggest gap between "software that
runs" and "software a top journal will publish." Everything above establishes
that the tool executes correctly; none of it establishes that its rankings are
*biologically* right. A reviewer will ask what the ranking is validated
against. Having an answer is likely worth more than any remaining engineering
work on this list.

---

## 8. Write the methods text from the manifest

Every output directory contains `run_manifest.json` (pipeline version, git
commit, full parameters, command line) and `software_versions.yml`.

**Why it matters:** write the methods section *from* those files rather than
from memory, and the stated parameters cannot drift from the ones actually
used. Cite the manifest's commit hash.

---

## 9. Publish to PyPI

```bash
cd /home/jmorrislab1/mkazemi/___ProJECTS___/creTarget/creTarget
python -m build
twine check dist/*
twine upload --repository testpypi dist/*     # test first
pip install -i https://test.pypi.org/simple/ cretarget   # in a scratch env
twine upload dist/*                            # then for real
```

**Why it matters:** `pip install cretarget` is the lowest-friction install and
the one reviewers will try. Note that PyPI alone gives users a tool that will
not run until they solve `bedtools`/`samtools`/Java themselves — which is why
bioconda matters more.

---

## 10. Publish to bioconda

Recipe PR against `bioconda-recipes`.

**Why it matters:** bioconda can declare `bedtools`, `samtools` and `openjdk`
as real dependencies, so `conda install -c bioconda cretarget` produces a
working tool in one command. For a bioinformatics audience this is the
expected distribution channel, and its absence reads as an unfinished release.

---

## 11. Build a container image

An Apptainer/Singularity image pinning every dependency, published to a
registry, wired into the existing `singularity` profile.

**Why it matters:** on clusters this collapses the entire dependency story to
one file, and it is the strongest reproducibility claim available — a reviewer
can rerun your pipeline years later with the exact software stack. You already
use `-profile singularity` for nf-core work, so the workflow is familiar.

---

## 12. Cut a Zenodo-archived release

Enable the GitHub–Zenodo integration, then tag:

```bash
cd /home/jmorrislab1/mkazemi/___ProJECTS___/creTarget/creTarget
git tag -a v0.2.0 -m "creTarget v0.2.0"
git push origin v0.2.0
```

Add the resulting DOI to `CITATION.cff` and the README.

**Why it matters:** journals increasingly require a DOI'd, immutable software
artifact rather than a bare GitHub link — a repo can change or disappear after
publication. This is what the paper cites.

---

## Already done

- Nextflow pipeline with local/SLURM/PBS/SGE/LSF profiles, retry escalation,
  `-resume`
- Per-chunk annotation and ranking (flat memory), verified byte-identical to
  an unchunked run
- 86 passing tests, GitHub Actions CI including a real end-to-end FlashFry run
- Checkpoints removed from the package, fetched from Zenodo with checksum
  verification and attribution; CI fails if a `.pth` reappears
- `cretarget doctor`, `--version`, `run_manifest.json`, `CITATION.cff`,
  `CHANGELOG.md`
- Eight latent SLURM bugs eliminated with the executor they lived in
- `conf/boltzmann.config` site profile and driver submission template
