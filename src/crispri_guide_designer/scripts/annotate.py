#!/usr/bin/env python3
"""
annotate.py

Merges FlashFry `score` output chunks, recovers each guide's genome
location and source region_id by parsing the `contig` field (format:
"<region_id>__<chrom>:<padded_start>-<padded_end>", written by
prep_regions.py into the FASTA header that was fed to `discover`), joins
to region metadata, and flags whether the guide falls within the CRISPRi
blast radius of its target reference point (SNP position, or enhancer
midpoint).

Column names below match a REAL flashfry `score` output sample (not just
docs):
  contig, start, stop, target, context, overflow, orientation,
  Doench2014OnTarget, DoenchCFD_maxOT, DoenchCFD_specificityscore,
  dangerous_GC, dangerous_polyT, dangerous_in_genome, Hsu2013,
  Moreno-Mateos2015OnTarget, JostCRISPRi_maxOT, JostCRISPRi_specificityscore,
  basesDiffToClosestHit, closestHitCount, 0-1-2-3-4_mismatch, FreeEnergy20,
  otCount
Note: `start`/`stop` in FlashFry's own output are offsets WITHIN the input
sequence (not genome coordinates) since --transformPositions is not used --
genome location comes entirely from parsing `contig`.
"""

import argparse
import sys
import pandas as pd


def parse_contig(contig):
    """'<region_id>__<chrom>:<padded_start>-<padded_end>' -> (region_id, chrom, padded_start, padded_end)"""
    try:
        region_id, coord = contig.split("__", 1)
        chrom, rng = coord.split(":")
        pstart, pend = rng.split("-")
        return region_id, chrom, int(pstart), int(pend)
    except (ValueError, IndexError):
        return None, None, None, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scored_chunks", nargs="+", required=True)
    ap.add_argument("--region_metadata", required=True)
    ap.add_argument("--blast_radius_bp", type=int, default=1000)
    ap.add_argument("--out", default="guides.annotated.tsv")
    args = ap.parse_args()

    chunks = []
    for f in args.scored_chunks:
        try:
            df = pd.read_csv(f, sep="\t")
            chunks.append(df)
        except pd.errors.EmptyDataError:
            print(f"WARNING: {f} is empty, skipping", file=sys.stderr)
    if not chunks:
        sys.exit("ERROR: no non-empty scored chunks to merge")

    guides = pd.concat(chunks, ignore_index=True)

    if "contig" not in guides.columns:
        sys.exit(f"ERROR: expected a 'contig' column. Found: {list(guides.columns)}")

    parsed = guides["contig"].apply(parse_contig)
    guides["region_id"] = [p[0] for p in parsed]
    guides["genome_chrom"] = [p[1] for p in parsed]
    guides["window_padded_start"] = [p[2] for p in parsed]
    guides["window_padded_end"] = [p[3] for p in parsed]

    bad_parse = guides["region_id"].isna().sum()
    if bad_parse:
        print(f"WARNING: {bad_parse} rows had unparseable contig field, dropping", file=sys.stderr)
        guides = guides.dropna(subset=["region_id"])

    # Genome coordinates of the guide = window start + in-window offset
    # (FlashFry's start/stop are 0-based offsets within the provided sequence)
    guides["genome_start"] = guides["window_padded_start"] + guides["start"]
    guides["genome_end"] = guides["window_padded_start"] + guides["stop"]

    meta = pd.read_csv(args.region_metadata, sep="\t")

    merged = guides.merge(meta, on="region_id", how="left", validate="many_to_one")

    unmatched = merged["target_ref"].isna().sum()
    if unmatched:
        print(f"WARNING: {unmatched} guides did not match a region_id in metadata", file=sys.stderr)

    guide_mid = (merged["genome_start"] + merged["genome_end"]) / 2
    merged["distance_to_target"] = (guide_mid - merged["target_ref"]).abs()
    merged["within_blast_radius"] = merged["distance_to_target"] <= args.blast_radius_bp

    merged.to_csv(args.out, sep="\t", index=False)
    print(f"Wrote {len(merged)} annotated guides -> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
