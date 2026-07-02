#!/usr/bin/env python3
"""
launch_dcas9_merge.py

Concatenates all per-chunk launch-dCas9 result CSVs (columns: grna, predict)
for both outcomes, then joins them onto the full (non-deduped) annotated
guide table by spacer sequence (many rows can share the same spacer across
different region_ids -- all get the same score, which is correct since the
score is a property of the sequence, not the genomic location).

Guides whose spacer wasn't scored (e.g. filtered for non-ATGC characters at
split time) get NaN in both new columns.
"""

import argparse
import glob
import sys
import pandas as pd


def load_predictions(result_dir):
    files = sorted(glob.glob(f"{result_dir}/chunk_*.csv"))
    if not files:
        sys.exit(f"ERROR: no chunk result files found in {result_dir}")
    dfs = [pd.read_csv(f) for f in files]
    combined = pd.concat(dfs, ignore_index=True)
    if combined["grna"].duplicated().any():
        # Expected to be unique since split_grnas deduped -- warn, don't fail.
        print(f"WARNING: duplicate grna values in {result_dir} predictions, keeping first occurrence", file=sys.stderr)
        combined = combined.drop_duplicates(subset="grna", keep="first")
    return combined.set_index("grna")["predict"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--annotated", required=True)
    ap.add_argument("--enhancer_dir", required=True)
    ap.add_argument("--promoter_dir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.annotated, sep="\t")
    if "target" not in df.columns:
        sys.exit(f"ERROR: expected a 'target' column. Found: {list(df.columns)}")
    df["spacer"] = df["target"].astype(str).str[:20]

    enh_scores = load_predictions(args.enhancer_dir)
    pro_scores = load_predictions(args.promoter_dir)

    df["launch_dCas9_enhancer"] = df["spacer"].map(enh_scores)
    df["launch_dCas9_promoter"] = df["spacer"].map(pro_scores)

    n_missing_enh = df["launch_dCas9_enhancer"].isna().sum()
    n_missing_pro = df["launch_dCas9_promoter"].isna().sum()
    if n_missing_enh or n_missing_pro:
        print(f"WARNING: {n_missing_enh} guides missing enhancer score, {n_missing_pro} missing promoter score (unscored/non-ATGC spacers)", file=sys.stderr)

    df.to_csv(args.out, sep="\t", index=False)
    print(f"Wrote {len(df)} guides with launch-dCas9 scores merged -> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
