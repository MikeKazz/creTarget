"""
annotate

Recovers each guide's genome location from FlashFry's `contig` field, joins
it to its region's metadata, flags whether it falls inside the CRISPRi blast
radius, and optionally records distance to the nearest ATAC-seq summit.

Runs on ONE chunk. That is the whole point of the Nextflow rewrite: the old
version concatenated every scored chunk into a single DataFrame, which is
what actually capped the genome-wide runs -- not FlashFry. Here memory is
O(chunk), so it is flat whether you run 3,000 regions or 2.35 million.

The region metadata table is global and can be large, so it is *streamed*
and filtered down to the region_ids present in this chunk rather than
loaded whole.

Column names match real FlashFry `score` output:
  contig, start, stop, target, context, overflow, orientation,
  Doench2014OnTarget, DoenchCFD_maxOT, DoenchCFD_specificityscore,
  dangerous_GC, dangerous_polyT, dangerous_in_genome, Hsu2013,
  Moreno-Mateos2015OnTarget, JostCRISPRi_maxOT, JostCRISPRi_specificityscore,
  basesDiffToClosestHit, closestHitCount, 0-1-2-3-4_mismatch, FreeEnergy20,
  otCount

Note that FlashFry's own start/stop are offsets WITHIN the input sequence
(--transformPositions is not used); genome location comes from `contig`.
"""

import argparse
import csv
import sys

import numpy as np
import pandas as pd

METADATA_NUMERIC = ["orig_start", "orig_end", "target_ref", "padded_start", "padded_end"]


def parse_contig(contig):
    """'<region_id>__<chrom>:<start>-<end>' -> (region_id, chrom, start, end)

    rsplit on "__" because region IDs can legitimately contain underscores;
    only the final one separates the ID from the coordinates.
    """
    try:
        region_id, coord = str(contig).rsplit("__", 1)
        chrom, rng = coord.rsplit(":", 1)
        start, end = rng.split("-")
        return region_id, chrom, int(start), int(end)
    except (ValueError, IndexError, AttributeError):
        return None, None, None, None


def load_metadata_subset(metadata_path, region_ids):
    """Streams the metadata TSV, keeping only rows for `region_ids`."""
    wanted = set(region_ids)
    rows = []
    with open(metadata_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None or "region_id" not in reader.fieldnames:
            raise SystemExit(
                f"ERROR: {metadata_path} has no 'region_id' column "
                f"(found: {reader.fieldnames})"
            )
        for row in reader:
            if row["region_id"] in wanted:
                rows.append(row)

    meta = pd.DataFrame(rows, columns=reader.fieldnames)
    for col in METADATA_NUMERIC:
        if col in meta.columns:
            meta[col] = pd.to_numeric(meta[col], errors="coerce")
    return meta


def load_atac_summits(narrowpeak_path):
    """{chrom: sorted array of absolute summit positions}.

    narrowPeak column 10 is the summit offset from start, or -1 when no
    summit was called -- in which case we fall back to the peak midpoint.
    """
    cols = ["chrom", "start", "end", "name", "score", "strand",
            "signalValue", "pValue", "qValue", "peak"]
    df = pd.read_csv(narrowpeak_path, sep="\t", header=None, names=cols,
                     comment="#", usecols=range(10))

    summit = (df["start"] + df["peak"]).astype("int64")
    no_summit = df["peak"] < 0
    summit = summit.mask(no_summit, (df["start"] + df["end"]) // 2)

    summits_by_chrom = {
        chrom: np.sort(group.to_numpy())
        for chrom, group in summit.groupby(df["chrom"])
    }
    print(f"Loaded {len(df)} ATAC peaks across {len(summits_by_chrom)} contig(s) "
          f"from {narrowpeak_path}", file=sys.stderr)
    return summits_by_chrom


def nearest_summit_distances(chroms, midpoints, summits_by_chrom):
    """Vectorized nearest-summit distance, per contig.

    The old implementation did a Python-level binary search per guide. At a
    few hundred thousand guides per chunk that is the slowest thing in the
    step for no reason.
    """
    chroms = pd.Series(chroms).reset_index(drop=True)
    midpoints = pd.Series(midpoints).reset_index(drop=True).to_numpy()
    out = np.full(len(chroms), np.nan)

    for chrom, idx in chroms.groupby(chroms).groups.items():
        summits = summits_by_chrom.get(chrom)
        positions = np.asarray(idx, dtype=int)
        if summits is None or len(summits) == 0:
            continue
        mids = midpoints[positions]
        insert = np.searchsorted(summits, mids)
        right = np.where(insert < len(summits), insert, len(summits) - 1)
        left = np.where(insert > 0, insert - 1, 0)
        out[positions] = np.minimum(
            np.abs(summits[right] - mids),
            np.abs(summits[left] - mids),
        )
    return out


def run(args):
    frames = []
    for path in args.scored:
        try:
            frames.append(pd.read_csv(path, sep="\t"))
        except pd.errors.EmptyDataError:
            print(f"WARNING: {path} is empty, skipping", file=sys.stderr)
    if not frames:
        raise SystemExit("ERROR: no non-empty scored input to annotate")

    guides = pd.concat(frames, ignore_index=True) if len(frames) > 1 else frames[0]

    if "contig" not in guides.columns:
        raise SystemExit(f"ERROR: expected a 'contig' column. "
                         f"Found: {list(guides.columns)}")

    parsed = [parse_contig(c) for c in guides["contig"]]
    guides["region_id"] = [p[0] for p in parsed]
    guides["genome_chrom"] = [p[1] for p in parsed]
    guides["window_padded_start"] = [p[2] for p in parsed]
    guides["window_padded_end"] = [p[3] for p in parsed]

    n_bad = guides["region_id"].isna().sum()
    if n_bad:
        print(f"WARNING: {n_bad} row(s) had an unparseable contig field, dropping",
              file=sys.stderr)
        guides = guides.dropna(subset=["region_id"]).reset_index(drop=True)

    if guides.empty:
        raise SystemExit("ERROR: no guides left after contig parsing")

    guides["genome_start"] = guides["window_padded_start"] + guides["start"]
    guides["genome_end"] = guides["window_padded_start"] + guides["stop"]

    meta = load_metadata_subset(args.region_metadata, guides["region_id"].unique())
    if meta.empty:
        raise SystemExit(
            f"ERROR: none of this chunk's region_ids were found in "
            f"{args.region_metadata}. Are the scored input and the metadata from "
            f"the same prep-regions run?"
        )
    if meta["region_id"].duplicated().any():
        dupes = meta.loc[meta["region_id"].duplicated(), "region_id"].unique()[:5]
        raise SystemExit(
            f"ERROR: duplicate region_id(s) in {args.region_metadata}, e.g. "
            f"{list(dupes)}. Re-run prep-regions to regenerate it."
        )

    merged = guides.merge(meta, on="region_id", how="left", validate="many_to_one")

    n_unmatched = merged["target_ref"].isna().sum()
    if n_unmatched:
        print(f"WARNING: {n_unmatched} guide(s) did not match a region_id in metadata",
              file=sys.stderr)

    guide_mid = (merged["genome_start"] + merged["genome_end"]) / 2
    merged["distance_to_target"] = (guide_mid - merged["target_ref"]).abs()
    merged["within_blast_radius"] = merged["distance_to_target"] <= args.blast_radius_bp

    if args.atac_narrowpeak:
        summits = load_atac_summits(args.atac_narrowpeak)
        merged["distance_to_atac_summit"] = nearest_summit_distances(
            merged["genome_chrom"], guide_mid, summits)
        n_nan = int(merged["distance_to_atac_summit"].isna().sum())
        if n_nan:
            print(f"WARNING: {n_nan} guide(s) on contigs with no ATAC peaks -- "
                  f"distance_to_atac_summit is NaN for those", file=sys.stderr)
    else:
        merged["distance_to_atac_summit"] = np.nan

    merged.to_csv(args.out, sep="\t", index=False)
    print(f"Wrote {len(merged)} annotated guides -> {args.out}", file=sys.stderr)
    return args.out


def add_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--scored", nargs="+", required=True,
                        help="FlashFry score output file(s) for this chunk")
    parser.add_argument("--region-metadata", dest="region_metadata", required=True)
    parser.add_argument("--blast-radius-bp", dest="blast_radius_bp", type=int,
                        default=1000)
    parser.add_argument("--atac-narrowpeak", dest="atac_narrowpeak", default=None,
                        help="ENCODE-standard 10-column narrowPeak file; adds "
                             "distance_to_atac_summit (informational, not used in "
                             "ranking)")
    parser.add_argument("--out", default="guides.annotated.tsv")
