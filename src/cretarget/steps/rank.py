"""
rank

HARD FILTERS (guide kept in the output but marked `exclude (<reason>)` in
the `rank` column):
  - seed_score >= 3            (see seed_score() below)
  - dangerous_GC != NONE       (FlashFry's GC-content flag)
  - dangerous_polyT != NONE    (FlashFry's poly-T / Pol-III terminator flag)
  - dangerous_in_genome > 1    (perfect-match count in the genome)
  - overflow != OK             (FlashFry off-target-search overflow)

seed_score: count of overlapping "GG" dinucleotides in the 12bp PAM-proximal
end of the 20bp spacer (GGG counts as 2, GGGG as 3). Guides at >= 3 carry
elevated risk of widespread dCas9 binding.

RANKING, per region_id, weighted composite. Weights are literal percentages
and are NOT renormalized, so composite_score keeps the same meaning whether
or not launch-dCas9 was run:
  - launch_dCas9_enhancer        30%  [only when the column is present]
  - JostCRISPRi_specificityscore 25%
  - Hsu2013 / 100                18%
  - DoenchCFD_specificityscore    8%
  - basesDiffToClosestHit / maxMismatch, capped at 1   5%
  - folding_norm (FreeEnergy20, per-region min-max)    2%
(launch_dCas9_promoter is retained as a reference column but deliberately
excluded from ranking.)

Ranking is per-region_id and every guide for a region lives in the same
chunk, so this runs per-chunk with no global merge -- which is what keeps
memory flat at genome scale.
"""

import argparse
import sys

import numpy as np
import pandas as pd

WEIGHTS = {
    "launch_enhancer": 0.30,
    "jost_crispri": 0.25,
    "hsu2013": 0.18,
    "cfd_specificity": 0.08,
    "closest_distance": 0.05,
    "folding": 0.02,
}
# Rounded because these are reported verbatim in an output column. Summing
# binary floats can leave noise (0.5800000000000001), which both looks wrong in
# a published table and makes exact comparison against the intended value fail.
# The weights are percentages to two decimals, so 6 places is lossless here.
MAX_WITHOUT_LAUNCH = round(sum(v for k, v in WEIGHTS.items() if k != "launch_enhancer"), 6)
MAX_WITH_LAUNCH = round(sum(WEIGHTS.values()), 6)

# Order matters: the first failing check is the one reported.
FILTER_ORDER = [
    ("fail_seed_score", "seed_score"),
    ("fail_dangerous_GC", "dangerous_GC"),
    ("fail_dangerous_polyT", "dangerous_polyT"),
    ("fail_multisite_genomic_match", "multisite_genomic_match"),
    ("fail_overflow", "overflow"),
]

REQUIRED_COLUMNS = ["target", "dangerous_GC", "dangerous_polyT",
                    "dangerous_in_genome", "overflow", "region_id"]
REQUIRED_SCORES = ["JostCRISPRi_specificityscore", "Hsu2013",
                   "DoenchCFD_specificityscore", "basesDiffToClosestHit"]


def seed_score(spacer):
    """Overlapping GG count in the 12bp PAM-proximal end of the spacer."""
    if not isinstance(spacer, str) or len(spacer) < 12:
        return None
    seed = spacer[-12:].upper()
    return sum(1 for i in range(len(seed) - 1) if seed[i:i + 2] == "GG")


def parse_in_genome_count(val):
    """'IN_GENOME=1' -> 1, 'NONE' -> 0, unparseable -> None"""
    if not isinstance(val, str):
        return None
    stripped = val.strip().upper()
    if stripped == "NONE":
        return 0
    if "=" in stripped:
        try:
            return int(stripped.split("=", 1)[1])
        except ValueError:
            return None
    return None


def minmax_norm(series):
    lo, hi = series.min(), series.max()
    if pd.isna(lo) or pd.isna(hi) or hi == lo:
        return pd.Series(0.5, index=series.index)
    return (series - lo) / (hi - lo)


def first_fail_reason(df):
    """Vectorized replacement for the old df.apply(..., axis=1).

    np.select evaluates the conditions in order and takes the first match,
    which is exactly the old semantics without a Python call per row.
    """
    conditions = [df[col] for col, _ in FILTER_ORDER]
    choices = [reason for _, reason in FILTER_ORDER]
    return pd.Series(
        np.select(conditions, choices, default=None), index=df.index, dtype=object
    )


def run(args):
    df = pd.read_csv(args.annotated, sep="\t")

    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise SystemExit(f"ERROR: required column(s) missing from {args.annotated}: "
                         f"{missing}. Found: {list(df.columns)}")
    missing_scores = [c for c in REQUIRED_SCORES if c not in df.columns]
    if missing_scores:
        raise SystemExit(f"ERROR: required scoring column(s) missing: {missing_scores}. "
                         f"Found: {list(df.columns)}")

    df["spacer"] = df["target"].astype(str).str[:20]
    df["seed_score"] = df["spacer"].map(seed_score)
    df["in_genome_count"] = df["dangerous_in_genome"].map(parse_in_genome_count)

    df["fail_seed_score"] = df["seed_score"].fillna(99) >= 3
    df["fail_dangerous_GC"] = df["dangerous_GC"].astype(str).str.upper() != "NONE"
    df["fail_dangerous_polyT"] = df["dangerous_polyT"].astype(str).str.upper() != "NONE"
    df["fail_multisite_genomic_match"] = df["in_genome_count"].fillna(99) > 1
    df["fail_overflow"] = df["overflow"].astype(str).str.upper() != "OK"

    df["exclude_reason"] = first_fail_reason(df)
    df["passes_hard_filters"] = df["exclude_reason"].isna()

    hsu_norm = (pd.to_numeric(df["Hsu2013"], errors="coerce") / 100.0).clip(0, 1)

    # --- closest-hit distance ---
    raw_closest = df["basesDiffToClosestHit"]
    closest_norm = (
        pd.to_numeric(raw_closest, errors="coerce") / args.max_mismatch
    ).clip(upper=1.0)
    # "UNK" is FlashFry's sentinel for "no off-target found", not for "not
    # computed": ClosestHit.scala emits it exactly when `closest ==
    # Int.MaxValue`, i.e. when no off-target with >=1 mismatch existed, and
    # writes closestOTCount 0 alongside. Maximally safe is therefore the
    # right reading. Any OTHER unparseable value is a genuine unknown and
    # must not quietly collect full credit for safety.
    is_sentinel = raw_closest.astype(str).str.strip().str.upper().eq("UNK")
    closest_norm = closest_norm.mask(is_sentinel, 1.0)
    unresolved = closest_norm.isna()
    if unresolved.any():
        print(f"WARNING: {int(unresolved.sum())} guide(s) have an unreadable "
              f"basesDiffToClosestHit (e.g. "
              f"{raw_closest[unresolved].iloc[0]!r}); the closest-distance term "
              f"is omitted for them rather than assumed safe.", file=sys.stderr)

    if "FreeEnergy20" in df.columns:
        df["_fe"] = pd.to_numeric(df["FreeEnergy20"], errors="coerce")
    else:
        df["_fe"] = 0.0
    folding_norm = df.groupby("region_id")["_fe"].transform(minmax_norm).fillna(0.5)

    # Components are assembled as a dict so a MISSING one can be dropped from
    # the denominator as well as the numerator.
    #
    # They used to be .fillna(0)'d into a fixed-denominator sum, which scores
    # "we never measured this" identically to "we measured it and it is the
    # worst possible". That is not a small distinction here: launch-dCas9
    # deliberately returns NaN for any spacer it cannot encode, so those
    # guides silently lost 30% of the scale and were published with a
    # confident-looking mid-range score. Excluding the component instead keeps
    # composite_score_normalized meaning "of what we could measure, how good
    # is this guide".
    terms = {
        "jost_crispri": pd.to_numeric(
            df["JostCRISPRi_specificityscore"], errors="coerce").clip(0, 1),
        "hsu2013": hsu_norm,
        "cfd_specificity": pd.to_numeric(
            df["DoenchCFD_specificityscore"], errors="coerce").clip(0, 1),
        "closest_distance": closest_norm,
        "folding": folding_norm,
    }
    has_launch = "launch_dCas9_enhancer" in df.columns
    if has_launch:
        # The only component with no documented range, so it is the only one
        # that could push the sum outside [0, max_possible].
        terms["launch_enhancer"] = pd.to_numeric(
            df["launch_dCas9_enhancer"], errors="coerce").clip(0, 1)

    composite = sum(WEIGHTS[name] * value.fillna(0) for name, value in terms.items())
    # Per-row: only the weights of the components this row actually has.
    max_possible = sum(WEIGHTS[name] * value.notna() for name, value in terms.items())

    df["composite_score"] = composite
    df["composite_score_max_possible"] = max_possible.round(6)
    # Scale-free companion, and the column to compare BETWEEN guides when any
    # component may be missing. Clipped because the float sum of the weights
    # can land a half-ulp above max_possible and print as 1.0000000000000002.
    df["composite_score_normalized"] = (
        composite / max_possible.replace(0, pd.NA)).clip(0, 1)
    df["composite_components_missing"] = sum(
        value.isna().astype(int) for value in terms.values())
    df = df.drop(columns=["_fe"])

    df["rank"] = None
    passing = df[df["passes_hard_filters"]]
    if not passing.empty:
        ranks = (
            passing.groupby("region_id")["composite_score"]
            .rank(method="first", ascending=False)
            .astype(int)
        )
        df.loc[ranks.index, "rank"] = ranks

    excluded = ~df["passes_hard_filters"]
    df.loc[excluded, "rank"] = (
        "exclude (" + df.loc[excluded, "exclude_reason"].astype(str) + ")"
    )
    # Blank the scores for excluded guides too, so the three columns agree.
    # `rank` was already labelled "exclude (...)", but composite_score was
    # left populated -- and a guide failing the poly-T filter (which will not
    # transcribe at all) was published at composite_score_normalized 0.999,
    # above most guides that passed. `rank` is a string there, so sorting the
    # table by score was the natural move and surfaced exactly those rows.
    for column in ("composite_score", "composite_score_normalized",
                   "composite_score_max_possible"):
        df.loc[excluded, column] = pd.NA

    df["_sort_rank"] = pd.to_numeric(df["rank"], errors="coerce")
    df = df.sort_values(
        by=["region_id", "passes_hard_filters", "_sort_rank"],
        ascending=[True, False, True],
    ).drop(columns=["_sort_rank"])

    df.to_csv(args.out, sep="\t", index=False)

    n_pass = int(df["passes_hard_filters"].sum())
    print(f"Wrote {len(df)} guides ({n_pass} passing hard filters) -> {args.out}",
          file=sys.stderr)
    # max_possible varies per row now, so report the usual (complete) case.
    full_scale = MAX_WITH_LAUNCH if has_launch else MAX_WITHOUT_LAUNCH
    if has_launch:
        print(f"composite_score includes launch_dCas9_enhancer "
              f"(max possible {full_scale:.2f}); launch_dCas9_promoter is reported "
              f"but not ranked on.", file=sys.stderr)
    else:
        print(f"NOTE: no launch_dCas9_enhancer column -- composite_score omits it "
              f"(max possible {full_scale:.2f}).", file=sys.stderr)
    incomplete = int((df["composite_components_missing"] > 0).sum())
    if incomplete:
        print(f"NOTE: {incomplete} guide(s) are missing at least one score "
              f"component; those components are excluded from both their "
              f"composite_score and their composite_score_max_possible rather "
              f"than counted as zero. Compare such guides on "
              f"composite_score_normalized.", file=sys.stderr)
    return args.out


def add_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--annotated", required=True)
    parser.add_argument("--max-mismatch", dest="max_mismatch", type=int, default=3,
                        help="the --maxMismatch used in FlashFry score, for "
                             "normalizing closest-off-target distance")
    parser.add_argument("--out", default="guides.ranked.tsv")
