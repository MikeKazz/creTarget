#!/usr/bin/env python3
"""
rank_report.py

Implements the project owner's real ranking methodology:

HARD FILTERS (guide excluded from ranking, marked `exclude (<reason>)` in
the `rank` column, but retained in output for visibility):
  - seed_score >= 3            (see seed_score() below)
  - dangerous_GC != NONE       (FlashFry's own GC-content flag)
  - dangerous_polyT != NONE    (FlashFry's own poly-T/Pol-III-terminator flag)
  - dangerous_in_genome        (perfect-match count in genome > 1)
  - overflow != OK             (FlashFry's own off-target-search-overflow flag)

seed_score: count of overlapping "GG" dinucleotides in the 12bp PAM-proximal
end of the 20bp spacer (i.e. GGG counts as 2, GGGG as 3, ...). Verified
against project owner's examples:
  ATACTGGATGCCTGTGGCTG -> 1
  AAACCCAAGTACTATCTTTG -> 0
Guides with seed_score >= 3 carry elevated risk of widespread dCas9 binding
and are excluded.

RANKING (guides that pass all hard filters), per target region_id, weighted
composite score. Weights are the project owner's literal percentages (NOT
renormalized to sum to 1) so composite_score's scale/meaning stays stable
as more components get wired in over time -- e.g. it means the same thing
whether or not launch_dCas9_enhancer happens to be present in a given run:
  - launch_dCas9_enhancer         (launch-dCas9 CNN enhancer on-target)     30%  [only if column present]
  - JostCRISPRi_specificityscore  (CRISPRi off-target specificity)          25%
  - Hsu2013 / 100                 (general off-target specificity)          18%
  - DoenchCFD_specificityscore    (CFD specificity)                          8%
  - closest-off-target-distance   (basesDiffToClosestHit / maxMismatch, capped at 1) 5%
  - folding_norm                  (RNA folding energy, per-region min-max normalized) 2%
  (launch_dCas9_promoter, 12% in the full intended scheme, is computed and
  kept as a reference column but deliberately EXCLUDED from ranking, per
  explicit instruction.)

If launch_dCas9_enhancer is absent from the input (i.e. that step hasn't
been run), it's simply dropped from the sum -- composite_score maxes out at
0.58 instead of 0.88, and this is printed as a warning so it's never silently
ambiguous which regime a given output was ranked under.

Rank 1 = best guide, assigned locally within each region_id's passing guides
(matches project owner's existing convention: ranks are local to the
current output file/subset, not global).
"""

import argparse
import sys
import pandas as pd

# ---- Literal weights (NOT renormalized) -----------------------------------
WEIGHTS = {
    "launch_enhancer": 0.30,   # only applied if launch_dCas9_enhancer column present
    "jost_crispri": 0.25,
    "hsu2013": 0.18,
    "cfd_specificity": 0.08,
    "closest_distance": 0.05,
    "folding": 0.02,
}
MAX_POSSIBLE_WITHOUT_LAUNCH = sum(v for k, v in WEIGHTS.items() if k != "launch_enhancer")
MAX_POSSIBLE_WITH_LAUNCH = sum(WEIGHTS.values())


def seed_score(spacer):
    if not isinstance(spacer, str) or len(spacer) < 12:
        return None
    seed = spacer[-12:]
    return sum(1 for i in range(len(seed) - 1) if seed[i:i + 2] == "GG")


def parse_in_genome_count(val):
    """'IN_GENOME=1' -> 1, 'NONE' -> 0, anything unparseable -> None"""
    if not isinstance(val, str):
        return None
    if val.strip().upper() == "NONE":
        return 0
    if "=" in val:
        try:
            return int(val.split("=", 1)[1])
        except ValueError:
            return None
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--annotated", required=True)
    ap.add_argument("--out", default="guides.ranked.tsv")
    ap.add_argument("--max_mismatch", type=int, default=3,
                     help="value used for --maxMismatch in FlashFry score, for normalizing closest-off-target distance")
    args = ap.parse_args()

    df = pd.read_csv(args.annotated, sep="\t")

    required_cols = ["target", "dangerous_GC", "dangerous_polyT", "dangerous_in_genome", "overflow"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        sys.exit(f"ERROR: expected columns missing from annotated input: {missing}. Found: {list(df.columns)}")

    # ---- Derived fields ----
    df["spacer"] = df["target"].astype(str).str[:20]
    df["seed_score"] = df["spacer"].apply(seed_score)
    df["in_genome_count"] = df["dangerous_in_genome"].apply(parse_in_genome_count)

    # ---- Hard filters ----
    fail_seed      = df["seed_score"].fillna(99) >= 3
    fail_gc        = df["dangerous_GC"].astype(str).str.upper() != "NONE"
    fail_polyt     = df["dangerous_polyT"].astype(str).str.upper() != "NONE"
    fail_ingenome  = df["in_genome_count"].fillna(99) > 1
    fail_overflow  = df["overflow"].astype(str).str.upper() != "OK"

    df["fail_seed_score"] = fail_seed
    df["fail_dangerous_GC"] = fail_gc
    df["fail_dangerous_polyT"] = fail_polyt
    df["fail_multisite_genomic_match"] = fail_ingenome
    df["fail_overflow"] = fail_overflow

    def first_fail_reason(row):
        if row["fail_seed_score"]:
            return "seed_score"
        if row["fail_dangerous_GC"]:
            return "dangerous_GC"
        if row["fail_dangerous_polyT"]:
            return "dangerous_polyT"
        if row["fail_multisite_genomic_match"]:
            return "multisite_genomic_match"
        if row["fail_overflow"]:
            return "overflow"
        return None

    df["exclude_reason"] = df.apply(first_fail_reason, axis=1)
    df["passes_hard_filters"] = df["exclude_reason"].isna()

    # ---- Composite score (passing guides only) ----
    for col in ["JostCRISPRi_specificityscore", "Hsu2013", "DoenchCFD_specificityscore", "basesDiffToClosestHit"]:
        if col not in df.columns:
            sys.exit(f"ERROR: expected scoring column '{col}' not found. Found: {list(df.columns)}")

    hsu_norm = (df["Hsu2013"] / 100.0).clip(0, 1)
    closest_norm = (pd.to_numeric(df["basesDiffToClosestHit"], errors="coerce") / args.max_mismatch).clip(upper=1.0)
    closest_norm = closest_norm.fillna(1.0)  # no off-target found within search radius -> treat as maximally safe

    # Folding energy: per-region min-max normalize (less negative / higher FreeEnergy20 -> better -> higher score)
    fe = pd.to_numeric(df["FreeEnergy20"], errors="coerce") if "FreeEnergy20" in df.columns else pd.Series(0, index=df.index)
    df["_fe"] = fe

    def minmax_norm(s):
        lo, hi = s.min(), s.max()
        if pd.isna(lo) or pd.isna(hi) or hi == lo:
            return pd.Series(0.5, index=s.index)  # single guide or no variance -> neutral
        return (s - lo) / (hi - lo)

    folding_norm = df.groupby("region_id")["_fe"].transform(minmax_norm)

    composite = (
        WEIGHTS["jost_crispri"] * df["JostCRISPRi_specificityscore"].fillna(0)
        + WEIGHTS["hsu2013"] * hsu_norm.fillna(0)
        + WEIGHTS["cfd_specificity"] * df["DoenchCFD_specificityscore"].fillna(0)
        + WEIGHTS["closest_distance"] * closest_norm
        + WEIGHTS["folding"] * folding_norm
    )

    has_launch_enhancer = "launch_dCas9_enhancer" in df.columns
    if has_launch_enhancer:
        composite = composite + WEIGHTS["launch_enhancer"] * pd.to_numeric(
            df["launch_dCas9_enhancer"], errors="coerce").fillna(0)
        max_possible = MAX_POSSIBLE_WITH_LAUNCH
    else:
        max_possible = MAX_POSSIBLE_WITHOUT_LAUNCH

    df["composite_score"] = composite
    df["composite_score_max_possible"] = max_possible
    df = df.drop(columns=["_fe"])

    # ---- Local rank per region_id, among passing guides only ----
    df["rank"] = None
    passing = df[df["passes_hard_filters"]].copy()
    passing["rank"] = (
        passing.groupby("region_id")["composite_score"]
        .rank(method="first", ascending=False)
        .astype(int)
    )
    df.loc[passing.index, "rank"] = passing["rank"]
    df.loc[~df["passes_hard_filters"], "rank"] = "exclude (" + df.loc[~df["passes_hard_filters"], "exclude_reason"] + ")"

    df = df.sort_values(by=["region_id", "passes_hard_filters"], ascending=[True, False])
    # within passing guides, sort by rank ascending (best first)
    df["_sort_rank"] = pd.to_numeric(df["rank"], errors="coerce")
    df = df.sort_values(by=["region_id", "passes_hard_filters", "_sort_rank"], ascending=[True, False, True])
    df = df.drop(columns=["_sort_rank"])

    df.to_csv(args.out, sep="\t", index=False)
    n_pass = int(df["passes_hard_filters"].sum())
    print(f"Wrote {len(df)} guides ({n_pass} passing hard filters, ranked locally per region_id) -> {args.out}", file=sys.stderr)
    if has_launch_enhancer:
        print(f"composite_score includes launch_dCas9_enhancer (max possible = {max_possible:.2f}). launch_dCas9_promoter is present but NOT used in ranking, per instruction.", file=sys.stderr)
    else:
        print(f"NOTE: launch_dCas9_enhancer column not found -- composite score omits it (max possible = {max_possible:.2f}). Run the launch-dCas9 scoring step first for the full ranking.", file=sys.stderr)


if __name__ == "__main__":
    main()
