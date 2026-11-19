"""
score_launch_dcas9

Adds launch-dCas9 CNN on-target efficacy scores (enhancerFitness and
promoterFitness) to one chunk's annotated guide table.

The old pipeline split this into three SLURM-coordinated stages (split ->
array-score -> merge). Per chunk it collapses into a single step: dedup the
spacers in this chunk, score each outcome once, map back onto every row.

Deduplication is now per-chunk rather than global. A spacer appearing in two
different chunks gets scored twice, which costs a little redundant compute
but removes an entire cross-chunk synchronization barrier. The score is a
property of the sequence, so the result is identical either way.

This module runs under creTarget's own interpreter; only launch-dCas9.py
itself is invoked with the torch environment's python (--python-bin). That
way creTarget does not have to be installed inside the conda env.
"""

import argparse
import contextlib
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

ATGC_RE = re.compile(r"^[ATGCatgc]+$")
# The protospacer length launch-dCas9's checkpoints were trained on, and the
# length of a SpCas9 spacer once the 3 bp NGG PAM is removed from FlashFry's
# 23 bp target.
SPACER_LENGTH = 20
OUTCOMES = {
    "enhancerFitness": "launch_dCas9_enhancer",
    "promoterFitness": "launch_dCas9_promoter",
}


def unique_valid_spacers(spacers, expected_length=SPACER_LENGTH):
    """Unique ATGC-only spacers, all of the same length.

    Two separate requirements, both enforced here because launch-dCas9
    enforces neither safely:

    * ATGC only. Its one-hot encoder exits on anything else, so these rows
      get NaN instead.
    * A UNIFORM length. Its preprocess_seq() takes the length of the FIRST
      sequence and applies it to the whole batch (util.py: `length =
      len(data[0])`). A shorter sequence later in the batch raises IndexError
      inside a bare `except` that only prints, leaving that row's encoding
      zero-padded; a longer one is silently truncated. Either way the score
      comes back plausible and wrong rather than failing. Worse, if the
      SHORT one happens to sort first, every other spacer in the batch is
      truncated to match it.
    """
    unique = pd.Series(spacers).dropna().unique()

    valid, non_atgc, wrong_length = [], 0, []
    for spacer in unique:
        text = str(spacer)
        if not ATGC_RE.match(text):
            non_atgc += 1
        elif len(text) != expected_length:
            wrong_length.append(text)
        else:
            valid.append(text)

    if non_atgc:
        print(f"WARNING: {non_atgc} unique spacer(s) contain non-ATGC characters; "
              f"they will get NaN scores", file=sys.stderr)
    if wrong_length:
        example = ", ".join(f"{s} ({len(s)} bp)" for s in wrong_length[:3])
        print(f"WARNING: {len(wrong_length)} unique spacer(s) are not "
              f"{expected_length} bp and will get NaN scores: {example}"
              f"{' ...' if len(wrong_length) > 3 else ''}\n"
              f"         launch-dCas9 assumes one length for the whole batch, "
              f"so mixing lengths silently corrupts every score in it.",
              file=sys.stderr)
    return valid


def score_outcome(spacers, outcome, launch_script, ckpt_dir, python_bin, workdir):
    """Runs launch-dCas9.py for one outcome. Returns {spacer: prediction}."""
    workdir = Path(workdir)
    chunk_dir = workdir / f"input_{outcome}"
    result_dir = workdir / f"output_{outcome}"
    chunk_dir.mkdir(parents=True, exist_ok=True)
    result_dir.mkdir(parents=True, exist_ok=True)

    chunk_csv = chunk_dir / "spacers.csv"
    pd.DataFrame({"protospacer": spacers}).to_csv(chunk_csv, index=False)

    cmd = [
        python_bin, str(launch_script),
        "--model", "CNN",
        "--variant", "seq",
        "--outcome", outcome,
        "--model_path", str(ckpt_dir).rstrip("/") + "/",
        "--test_path", str(chunk_dir) + "/",
        "--test_filename", chunk_csv.name,
        "--result_path", str(result_dir) + "/",
    ]
    # Deliberately not touching CUDA_VISIBLE_DEVICES: if the scheduler
    # allocated a GPU it has already set visibility, and forcing a value
    # here would silently disable GPU use where it is available.
    print(f"Running: {' '.join(cmd)}", file=sys.stderr)
    proc = subprocess.run(cmd)
    if proc.returncode != 0:
        raise SystemExit(f"ERROR: launch-dCas9.py failed for outcome {outcome} "
                         f"(exit {proc.returncode})")

    pred_csv = result_dir / f"prediction-seq-{outcome}.csv"
    if not pred_csv.exists():
        raise SystemExit(f"ERROR: expected prediction output not found: {pred_csv}")

    preds = pd.read_csv(pred_csv)
    if "predict" not in preds.columns:
        raise SystemExit(f"ERROR: no 'predict' column in {pred_csv}. "
                         f"Found: {list(preds.columns)}")
    if len(preds) != len(spacers):
        raise SystemExit(
            f"ERROR: row count mismatch for {outcome}: {len(spacers)} spacers in, "
            f"{len(preds)} predictions out. Refusing to guess at the alignment."
        )

    return _join_predictions(preds, spacers, pred_csv, outcome)


def _join_predictions(preds, spacers, pred_csv, outcome):
    """Attach each prediction to its spacer BY KEY, not by row position.

    launch-dCas9 writes the input protospacer alongside its score in a
    `grna` column (runCNN.py: columns=['grna', 'predict']), so the exact
    correspondence is available and there is no reason to infer it.

    Position was the only thing tying a score to a guide before. The row
    count matching says nothing about the ORDER, so a launch-dCas9 revision
    that sorted or regrouped its output would have silently attached every
    score to the wrong spacer -- the pipeline would exit 0 and the ranking
    would be wrong. `--launch-script` points at a user-supplied clone of an
    external project, so that revision is not under creTarget's control.
    """
    key_column = next((c for c in ("grna", "protospacer") if c in preds.columns), None)
    if key_column is None:
        print(f"WARNING: {pred_csv.name} has no 'grna'/'protospacer' column, so "
              f"scores can only be matched to spacers by row order. If this "
              f"build of launch-dCas9 reorders its output, the scores will be "
              f"attached to the wrong guides. Columns: {list(preds.columns)}",
              file=sys.stderr)
        return dict(zip(spacers, preds["predict"].to_numpy()))

    returned = preds[key_column].astype(str)
    if returned.duplicated().any():
        duplicated = returned[returned.duplicated()].unique()[:3]
        raise SystemExit(
            f"ERROR: {pred_csv.name} lists the same spacer more than once "
            f"({', '.join(duplicated)}); cannot assign scores unambiguously."
        )

    values = preds["predict"].to_numpy()
    exact = dict(zip(returned, values))
    if all(s in exact for s in spacers):
        return {s: exact[s] for s in spacers}

    # Fall back to a case-insensitive match. hg38 is soft-masked, so a spacer
    # over a repeat arrives lowercase, and a tool that normalises case would
    # otherwise look like it had returned different sequences entirely.
    folded = {}
    for key, value in zip(returned, values):
        folded.setdefault(key.upper(), value)
    missing = [s for s in spacers if s.upper() not in folded]
    if missing:
        raise SystemExit(
            f"ERROR: {len(missing)} spacer(s) sent for {outcome} scoring are "
            f"absent from {pred_csv.name}, e.g. {', '.join(missing[:3])}.\n"
            f"The prediction table does not describe the guides that were "
            f"submitted; refusing to produce scores that may be misaligned."
        )
    return {s: folded[s.upper()] for s in spacers}


def run(args):
    launch_script = Path(args.launch_script)
    if not launch_script.exists():
        raise SystemExit(f"ERROR: launch-dCas9 script not found: {launch_script}")

    ckpt_dir = Path(args.ckpt_dir)
    missing_ckpts = [
        name for name in ("ckpt-seq-enhancerFitness.pth", "ckpt-seq-promoterFitness.pth")
        if not (ckpt_dir / name).exists()
    ]
    if missing_ckpts:
        raise SystemExit(
            f"ERROR: checkpoint(s) missing from {ckpt_dir}: {', '.join(missing_ckpts)}\n"
            f"Fetch them with:  cretarget fetch-models"
        )

    df = pd.read_csv(args.annotated, sep="\t")
    if "target" not in df.columns:
        raise SystemExit(f"ERROR: expected a 'target' column. Found: {list(df.columns)}")

    df["spacer"] = df["target"].astype(str).str[:20]
    spacers = unique_valid_spacers(df["spacer"])

    if not spacers:
        print("WARNING: no valid spacers in this chunk; writing NaN score columns",
              file=sys.stderr)
        for column in OUTCOMES.values():
            df[column] = pd.NA
        df.to_csv(args.out, sep="\t", index=False)
        return args.out

    print(f"Scoring {len(spacers)} unique spacer(s) from {len(df)} guide row(s)",
          file=sys.stderr)

    # Scratch goes to $TMPDIR (node-local, and auto-cleaned on SLURM sites
    # using job_container/tmpfs), not the current directory. dir=os.getcwd()
    # overrode tempfile's own TMPDIR handling, which put per-spacer CSVs on
    # the shared filesystem -- slower, and an outright failure when the
    # working directory is read-only. ExitStack removes the directory
    # afterwards; mkdtemp left one behind on every interactive run.
    # Nextflow passes --workdir . explicitly, so pipeline runs are unchanged.
    with contextlib.ExitStack() as stack:
        workdir = args.workdir or stack.enter_context(
            tempfile.TemporaryDirectory(prefix="launch_dcas9_"))
        for outcome, column in OUTCOMES.items():
            predictions = score_outcome(
                spacers, outcome, launch_script, ckpt_dir, args.python_bin, workdir)
            df[column] = df["spacer"].map(predictions)
            n_missing = int(df[column].isna().sum())
            if n_missing:
                print(f"WARNING: {n_missing} guide(s) have no {column} score "
                      f"(unscored/non-ATGC spacer)", file=sys.stderr)

    df.to_csv(args.out, sep="\t", index=False)
    print(f"Wrote {len(df)} guides with launch-dCas9 scores -> {args.out}",
          file=sys.stderr)
    return args.out


def add_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--annotated", required=True)
    parser.add_argument("--launch-script", dest="launch_script",
                        help="path to launch-dCas9.py")
    parser.add_argument("--ckpt-dir", dest="ckpt_dir",
                        help="directory holding the .pth checkpoints "
                             "(default: the `cretarget fetch-models` cache)")
    parser.add_argument("--python-bin", dest="python_bin",
                        help="python from the launch-dCas9 conda env (needs torch)")
    parser.add_argument("--workdir", default=None,
                        help="scratch directory for intermediate files "
                             "(default: a temp dir under $TMPDIR, removed "
                             "when the step finishes)")
    parser.add_argument("--out", default="guides.annotated.with_launch.tsv")
