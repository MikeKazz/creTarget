#!/usr/bin/env python3
"""
launch_dcas9_split.py

Extracts the unique, valid (ATGC-only) 20bp spacer sequences from an
annotated guide table and splits them into fixed-size chunk CSVs (one
'protospacer' column each), matching the input format launch-dCas9.py
expects. Mirrors the project owner's existing split_grnas.py, but operates
on deduped spacers (since our guide table can have the same spacer
sequence appear under multiple region_ids) rather than raw rows -- saves
redundant scoring, and the merge step joins back to every occurrence.

Writes:
  <outdir>/chunk_0000.csv, chunk_0001.csv, ...   (column: protospacer)
  <outdir>/n_chunks.txt
"""

import argparse
import re
import sys
import pandas as pd

ATGC_RE = re.compile(r'^[ATGCatgc]+$')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--annotated", required=True)
    ap.add_argument("--chunk_size", type=int, default=200000)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    import os
    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.annotated, sep="\t")
    if "target" not in df.columns:
        sys.exit(f"ERROR: expected a 'target' column. Found: {list(df.columns)}")

    spacers = df["target"].astype(str).str[:20]
    unique_spacers = spacers.dropna().unique()
    valid_mask = pd.Series(unique_spacers).apply(lambda s: bool(ATGC_RE.match(s)))
    valid_spacers = pd.Series(unique_spacers)[valid_mask].tolist()
    n_invalid = (~valid_mask).sum()
    if n_invalid:
        print(f"WARNING: {n_invalid} unique spacers contain non-ATGC characters, excluded from scoring (will get NaN)", file=sys.stderr)

    n_chunks = max(1, (len(valid_spacers) + args.chunk_size - 1) // args.chunk_size)
    for i in range(n_chunks):
        chunk = valid_spacers[i * args.chunk_size:(i + 1) * args.chunk_size]
        pd.DataFrame({"protospacer": chunk}).to_csv(
            f"{args.outdir}/chunk_{i:04d}.csv", index=False)

    with open(f"{args.outdir}/n_chunks.txt", "w") as f:
        f.write(str(n_chunks))

    print(f"{len(valid_spacers)} unique valid spacers -> {n_chunks} chunk(s) (size {args.chunk_size}) -> {args.outdir}", file=sys.stderr)
    print(f"Submit with: --array=0-{n_chunks - 1}", file=sys.stderr)


if __name__ == "__main__":
    main()
