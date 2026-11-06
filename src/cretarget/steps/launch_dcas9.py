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
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

ATGC_RE = re.compile(r"^[ATGCatgc]+$")
OUTCOMES = {
    "enhancerFitness": "launch_dCas9_enhancer",
    "promoterFitness": "launch_dCas9_promoter",
}


def unique_valid_spacers(spacers):
    """Unique ATGC-only spacers. Anything else can't be one-hot encoded by
    launch-dCas9 and would crash it, so those rows get NaN instead."""
    unique = pd.Series(spacers).dropna().unique()
    valid = [s for s in unique if ATGC_RE.match(str(s))]
    n_invalid = len(unique) - len(valid)
    if n_invalid:
        print(f"WARNING: {n_invalid} unique spacer(s) contain non-ATGC characters; "
              f"they will get NaN scores", file=sys.stderr)
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
    return dict(zip(spacers, preds["predict"].to_numpy()))


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

    workdir = args.workdir or tempfile.mkdtemp(prefix="launch_dcas9_", dir=os.getcwd())
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
                             "(default: a temp dir in the current directory)")
    parser.add_argument("--out", default="guides.annotated.with_launch.tsv")
