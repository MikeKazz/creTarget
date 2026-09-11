"""
seed_finder_integrate.py

Integrates the external seed-finder-lite tool (off-target seed/TSS-proximity
search, Zenodo DOI 10.5281/zenodo.21599135) into the pipeline as an optional
final step. seed-finder-lite itself is NOT bundled or modified -- this
module just builds its input, invokes it, and folds a per-guide SUMMARY of
its hit-level output back onto the main guide table. Full hit-level detail
is preserved separately, since a guide can have zero-to-many off-target
hits and the main table is one row per guide.
"""

import os
import subprocess
import sys

import pandas as pd


def build_guide_id(row):
    """Unique per-row ID: region_id + genome coordinates + orientation.
    Orientation is needed because two guides (opposite strand) can share
    the same genomic start/end span."""
    orientation = row.get("orientation", "NA") if hasattr(row, "get") else row["orientation"]
    return f"{row['region_id']}_{row['genome_start']}_{row['genome_end']}_{orientation}"


def make_unique_guide_ids(ranked_df):
    """Builds guide_ids via build_guide_id, then appends a per-duplicate
    counter to any that still collide (safety net -- region_id + coords +
    orientation is expected to be unique, but this guarantees it rather
    than assuming it)."""
    base_ids = ranked_df.apply(build_guide_id, axis=1)
    if not base_ids.duplicated().any():
        return base_ids
    dup_counter = {}
    final_ids = []
    for gid in base_ids:
        if gid not in dup_counter:
            dup_counter[gid] = 0
            final_ids.append(gid)
        else:
            dup_counter[gid] += 1
            final_ids.append(f"{gid}_dup{dup_counter[gid]}")
    return pd.Series(final_ids, index=base_ids.index)


def run_seed_finder(ranked_df, seed_finder_dir, min_seed_length, tss_window,
                     outdir, python_bin=None):
    """Runs seed-finder-lite on every guide in ranked_df. Returns
    (summary_df, raw_hits_path). summary_df has columns:
    guide_id, seedfinder_hit_count, seedfinder_closest_gene,
    seedfinder_closest_tss_distance, seedfinder_closest_seed_length.
    """
    python_bin = python_bin or sys.executable
    sf_outdir = os.path.join(outdir, "seed_finder")
    os.makedirs(sf_outdir, exist_ok=True)

    if "spacer" not in ranked_df.columns:
        raise SystemExit("ERROR: expected a 'spacer' column in the ranked guide table")

    guide_ids = make_unique_guide_ids(ranked_df)

    input_tsv = os.path.join(sf_outdir, "seed_finder_input.tsv")
    pd.DataFrame({"guide_id": guide_ids, "sequence": ranked_df["spacer"]}).to_csv(
        input_tsv, sep="\t", index=False)

    raw_hits_path = os.path.join(sf_outdir, "raw_hits.csv")
    promoter_fa = os.path.join(seed_finder_dir, "refs", "promoter_trim_2000_pad23.fa.gz")
    tss_table = os.path.join(seed_finder_dir, "refs", "tss_table.tsv")
    src_dir = os.path.join(seed_finder_dir, "src")

    cmd = [
        python_bin, "-m", "targetfinder", "seed-find-batch", input_tsv,
        "--promoter-fa", promoter_fa,
        "--tss-table", tss_table,
        "--min-seed-length", str(min_seed_length),
        "--tss-window", str(tss_window),
        "-o", raw_hits_path,
    ]
    env = dict(os.environ)
    env["PYTHONPATH"] = src_dir + (os.pathsep + env["PYTHONPATH"] if "PYTHONPATH" in env else "")

    print(f"Running seed-finder-lite on {len(ranked_df)} guides "
          f"(min_seed_length={min_seed_length}, tss_window={tss_window})...", file=sys.stderr)
    log_path = os.path.join(sf_outdir, "seed_finder.log")
    with open(log_path, "w") as logf:
        proc = subprocess.run(cmd, env=env, stdout=logf, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        raise SystemExit(f"ERROR: seed-finder-lite failed (exit {proc.returncode}), see {log_path}")

    if not os.path.exists(raw_hits_path) or os.path.getsize(raw_hits_path) == 0:
        # No hits at all is a valid outcome (0 rows) -- but if the file is
        # missing entirely, something went wrong.
        if not os.path.exists(raw_hits_path):
            raise SystemExit(f"ERROR: expected output not found: {raw_hits_path}")

    try:
        hits = pd.read_csv(raw_hits_path)
    except pd.errors.EmptyDataError:
        hits = pd.DataFrame(columns=["guide_id", "tss_distance", "nearest_tss_gene", "seed_length"])

    all_ids = pd.DataFrame({"guide_id": guide_ids})

    if len(hits) == 0:
        summary = all_ids.copy()
        summary["seedfinder_hit_count"] = 0
        summary["seedfinder_closest_gene"] = pd.NA
        summary["seedfinder_closest_tss_distance"] = pd.NA
        summary["seedfinder_closest_seed_length"] = pd.NA
    else:
        hits["_abs_tss_distance"] = hits["tss_distance"].abs()
        hit_counts = hits.groupby("guide_id").size().rename("seedfinder_hit_count")
        closest_idx = hits.groupby("guide_id")["_abs_tss_distance"].idxmin()
        closest = hits.loc[closest_idx, ["guide_id", "nearest_tss_gene", "tss_distance", "seed_length"]].rename(
            columns={
                "nearest_tss_gene": "seedfinder_closest_gene",
                "tss_distance": "seedfinder_closest_tss_distance",
                "seed_length": "seedfinder_closest_seed_length",
            })
        summary = all_ids.merge(hit_counts, on="guide_id", how="left").merge(closest, on="guide_id", how="left")
        summary["seedfinder_hit_count"] = summary["seedfinder_hit_count"].fillna(0).astype(int)

    return summary, raw_hits_path, guide_ids
