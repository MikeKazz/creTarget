"""
seed_finder

Optional final annotation: runs the external seed-finder-lite tool
(off-target seed / TSS-proximity search, Zenodo DOI 10.5281/zenodo.21599135)
over one chunk's ranked guides and folds a per-guide summary back onto the
table.

seed-finder-lite is neither bundled nor modified. This builds its input,
invokes it, and summarizes its hit-level output -- a guide can have zero to
many hits, so the full hit table is preserved separately while the main
table stays one row per guide.
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd

SUMMARY_COLUMNS = [
    "seedfinder_hit_count",
    "seedfinder_closest_gene",
    "seedfinder_closest_tss_distance",
    "seedfinder_closest_seed_length",
]


def build_guide_ids(df):
    """region_id + coordinates + orientation, with a collision safety net.

    Orientation matters because two guides on opposite strands can share the
    same genomic span.
    """
    orientation = df["orientation"] if "orientation" in df.columns else "NA"
    base = (
        df["region_id"].astype(str) + "_"
        + df["genome_start"].astype(str) + "_"
        + df["genome_end"].astype(str) + "_"
        + pd.Series(orientation, index=df.index).astype(str)
    )
    if not base.duplicated().any():
        return base

    counts = base.groupby(base).cumcount()
    return base.where(counts == 0, base + "_dup" + counts.astype(str))


def summarize_hits(hits, guide_ids):
    all_ids = pd.DataFrame({"guide_id": guide_ids.to_numpy()})

    if hits is None or len(hits) == 0:
        summary = all_ids.copy()
        summary["seedfinder_hit_count"] = 0
        for column in SUMMARY_COLUMNS[1:]:
            summary[column] = pd.NA
        return summary

    hits = hits.copy()
    hits["_abs_tss_distance"] = hits["tss_distance"].abs()
    hit_counts = hits.groupby("guide_id").size().rename("seedfinder_hit_count")
    closest_idx = hits.groupby("guide_id")["_abs_tss_distance"].idxmin()
    closest = hits.loc[
        closest_idx, ["guide_id", "nearest_tss_gene", "tss_distance", "seed_length"]
    ].rename(columns={
        "nearest_tss_gene": "seedfinder_closest_gene",
        "tss_distance": "seedfinder_closest_tss_distance",
        "seed_length": "seedfinder_closest_seed_length",
    })

    summary = (
        all_ids
        .merge(hit_counts, on="guide_id", how="left")
        .merge(closest, on="guide_id", how="left")
    )
    summary["seedfinder_hit_count"] = (
        summary["seedfinder_hit_count"].fillna(0).astype(int)
    )
    return summary


def run(args):
    seed_finder_dir = Path(args.seed_finder_dir)
    src_dir = seed_finder_dir / "src"
    promoter_fa = seed_finder_dir / "refs" / "promoter_trim_2000_pad23.fa.gz"
    tss_table = seed_finder_dir / "refs" / "tss_table.tsv"

    for path, label in ((src_dir, "src/"), (promoter_fa, "refs/promoter fasta"),
                        (tss_table, "refs/tss_table.tsv")):
        if not path.exists():
            raise SystemExit(
                f"ERROR: seed-finder-lite {label} not found at {path}\n"
                f"Check that --seed-finder-dir points at the seed-finder-lite "
                f"directory (the one containing src/ and refs/)."
            )

    df = pd.read_csv(args.ranked, sep="\t")
    for column in ("spacer", "region_id", "genome_start", "genome_end"):
        if column not in df.columns:
            raise SystemExit(f"ERROR: expected a {column!r} column in {args.ranked}")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    guide_ids = build_guide_ids(df)
    input_tsv = outdir / "seed_finder_input.tsv"
    pd.DataFrame({"guide_id": guide_ids, "sequence": df["spacer"]}).to_csv(
        input_tsv, sep="\t", index=False)

    raw_hits = Path(args.raw_hits) if args.raw_hits else outdir / "raw_hits.csv"
    cmd = [
        args.python_bin or sys.executable, "-m", "targetfinder", "seed-find-batch",
        str(input_tsv),
        "--promoter-fa", str(promoter_fa),
        "--tss-table", str(tss_table),
        "--min-seed-length", str(args.min_seed_length),
        "--tss-window", str(args.tss_window),
        "-o", str(raw_hits),
    ]

    env = dict(os.environ)
    env["PYTHONPATH"] = str(src_dir) + (
        os.pathsep + env["PYTHONPATH"] if "PYTHONPATH" in env else "")

    print(f"Running seed-finder-lite on {len(df)} guides "
          f"(min_seed_length={args.min_seed_length}, tss_window={args.tss_window})",
          file=sys.stderr)
    log_path = outdir / "seed_finder.log"
    with open(log_path, "w") as log:
        proc = subprocess.run(cmd, env=env, stdout=log, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        raise SystemExit(f"ERROR: seed-finder-lite failed (exit {proc.returncode}), "
                         f"see {log_path}")

    if not raw_hits.exists():
        raise SystemExit(f"ERROR: expected output not found: {raw_hits}")

    # Zero hits is a legitimate result; an absent file is not.
    try:
        hits = pd.read_csv(raw_hits)
    except pd.errors.EmptyDataError:
        hits = None

    summary = summarize_hits(hits, guide_ids)
    df["_guide_id"] = guide_ids.to_numpy()
    final = df.merge(
        summary.rename(columns={"guide_id": "_guide_id"}), on="_guide_id", how="left"
    ).drop(columns=["_guide_id"])

    final.to_csv(args.out, sep="\t", index=False)
    print(f"Wrote {len(final)} guides with seed-finder-lite columns -> {args.out}",
          file=sys.stderr)
    print(f"Hit-level detail retained at {raw_hits}", file=sys.stderr)
    return args.out


def add_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--ranked", required=True)
    parser.add_argument("--seed-finder-dir", dest="seed_finder_dir",
                        help="seed-finder-lite directory (containing src/ and refs/)")
    parser.add_argument("--min-seed-length", dest="min_seed_length", type=int,
                        required=True,
                        help="minimum exact-match seed length (e.g. 8). No default "
                             "-- choose deliberately, see seed-finder-lite's docs.")
    parser.add_argument("--tss-window", dest="tss_window", type=int, required=True,
                        help="max distance (bp) from a TSS to count as a hit "
                             "(e.g. 1000). No default -- choose deliberately.")
    parser.add_argument("--python-bin", dest="python_bin", default=None,
                        help="python with pyfaidx installed (default: current)")
    parser.add_argument("--outdir", default="seed_finder")
    parser.add_argument("--raw-hits", dest="raw_hits", default=None)
    parser.add_argument("--out", default="guides.final.tsv")
