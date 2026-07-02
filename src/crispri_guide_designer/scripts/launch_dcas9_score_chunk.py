#!/usr/bin/env python3
"""
launch_dcas9_score_chunk.py

Scores one chunk CSV (single 'protospacer' column, already deduped and
ATGC-filtered by launch_dcas9_split.py) with launch-dCas9.py for one
outcome (enhancerFitness or promoterFitness). Simpler than the project
owner's original sbatch inline scripts since filtering/NaN-reinsertion
already happened at split time -- no invalid sequences reach this step.

MUST be invoked with the launch-dCas9 conda env's python (needs
torch/torchvision via launch-dCas9.py, and pandas for I/O here).

Writes: <result_csv> with columns [grna, predict]
"""

import argparse
import os
import subprocess
import sys
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chunk_csv", required=True)
    ap.add_argument("--outcome", required=True, choices=["enhancerFitness", "promoterFitness"])
    ap.add_argument("--launch_script", required=True)
    ap.add_argument("--ckpt_dir", required=True)
    ap.add_argument("--tmp_dir", required=True)
    ap.add_argument("--result_csv", required=True)
    ap.add_argument("--python_bin", default=sys.executable,
                     help="python interpreter to invoke launch_script with (defaults to whatever's running this script)")
    args = ap.parse_args()

    os.makedirs(args.tmp_dir, exist_ok=True)
    # NOTE: deliberately NOT forcing CUDA_VISIBLE_DEVICES here. If the job
    # was submitted with --gres=gpu:N (see local_paths.env), SLURM sets GPU
    # visibility already; forcing it to "" here would silently disable GPU
    # use even on clusters where it's available and requested.

    chunk_df = pd.read_csv(args.chunk_csv)
    if "protospacer" not in chunk_df.columns:
        sys.exit(f"ERROR: expected a 'protospacer' column in {args.chunk_csv}")

    test_filename = os.path.basename(args.chunk_csv)
    # launch-dCas9.py reads from test_path/test_filename directly; point it at
    # the chunk file's own directory so we don't need to copy it.
    test_path = os.path.dirname(os.path.abspath(args.chunk_csv)) + "/"

    cmd = [
        args.python_bin, args.launch_script,
        "--model", "CNN",
        "--variant", "seq",
        "--outcome", args.outcome,
        "--model_path", args.ckpt_dir,
        "--test_path", test_path,
        "--test_filename", test_filename,
        "--result_path", args.tmp_dir + "/",
    ]
    print(f"Running: {' '.join(cmd)}", file=sys.stderr)
    result = subprocess.run(cmd)
    if result.returncode != 0:
        sys.exit(f"ERROR: launch-dCas9.py failed for {args.chunk_csv} / {args.outcome}")

    pred_csv = os.path.join(args.tmp_dir, f"prediction-seq-{args.outcome}.csv")
    if not os.path.exists(pred_csv):
        sys.exit(f"ERROR: expected prediction output not found: {pred_csv}")

    preds_df = pd.read_csv(pred_csv)
    if "predict" not in preds_df.columns:
        sys.exit(f"ERROR: expected a 'predict' column in {pred_csv}. Found: {list(preds_df.columns)}")
    if len(preds_df) != len(chunk_df):
        sys.exit(f"ERROR: row count mismatch: {len(chunk_df)} input spacers vs {len(preds_df)} predictions")

    out = pd.DataFrame({"grna": chunk_df["protospacer"].values, "predict": preds_df["predict"].values})
    os.makedirs(os.path.dirname(os.path.abspath(args.result_csv)), exist_ok=True)
    out.to_csv(args.result_csv, index=False)
    print(f"Wrote {len(out)} predictions -> {args.result_csv}", file=sys.stderr)


if __name__ == "__main__":
    main()
