#!/usr/bin/env python3
"""
prep_regions.py

Given a BED file of target regions (either SNP-type: single points, or
enhancer-type: intervals) and a reference FASTA, this script:
  1. Parses the BED (chrom, start, end, name).
  2. Flags pairs of regions within `link_distance_bp` of each other as
     "linked" (one guide near the shared boundary/blast radius could
     silence both) -- does NOT merge them; both stay as separate output
     regions with a `linked_regions` annotation. Only runs if
     --dedup_blast_radius is passed.
  3. Computes a search window per region:
       - SNP type:      point +/- snp_window_bp on each side
       - enhancer type: BED interval +/- flank_bp on each side
  4. Extracts the padded sequence for each region via `bedtools getfasta`.
  5. Writes FASTA headers as "<region_id>__<chrom>:<padded_start>-<padded_end>"
     -- genome coordinates are encoded directly in the header (mirroring the
     existing FlashFry workflow's id_<n>_chr:start-end_<seq> convention),
     so downstream steps can recover genome location from FlashFry's
     `contig` column without needing --transformPositions.
  6. Writes regions.metadata.tsv: region_id, chrom, orig_start, orig_end,
     input_type, target_ref (SNP position, or enhancer midpoint --
     the reference point used for blast-radius distance), padded_start,
     padded_end, linked_regions.

Requires: bedtools (getfasta) on PATH, and a faidx-indexed reference (.fai).
"""

import argparse
import csv
import subprocess
import sys
from pathlib import Path


def read_bed(bed_path, prefix):
    regions = []
    with open(bed_path) as fh:
        for i, line in enumerate(fh):
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            name = fields[3] if len(fields) > 3 and fields[3] else f"{prefix}_{i+1}"
            regions.append({"chrom": chrom, "start": start, "end": end, "name": name})
    return regions


def flag_linked_regions(regions, link_distance_bp):
    """Mark pairs of regions on the same chrom within link_distance_bp of each
    other (gap between end of one and start of the next). Does not merge --
    just annotates."""
    linked = {r["name"]: [] for r in regions}
    by_chrom = {}
    for r in regions:
        by_chrom.setdefault(r["chrom"], []).append(r)
    for chrom, regs in by_chrom.items():
        regs_sorted = sorted(regs, key=lambda r: r["start"])
        for a, b in zip(regs_sorted, regs_sorted[1:]):
            gap = b["start"] - a["end"]
            if gap < link_distance_bp:
                linked[a["name"]].append(b["name"])
                linked[b["name"]].append(a["name"])
    return linked


def compute_window(r, input_type, snp_window_bp, flank_bp, chrom_sizes):
    chrom_len = chrom_sizes.get(r["chrom"])
    if input_type == "snp":
        # treat as a point: use start (BED is 0-based half-open; a 1bp SNP
        # entry chr:pos-1:pos has its point at `start`)
        point = r["start"]
        padded_start = max(0, point - snp_window_bp)
        padded_end = point + snp_window_bp
        target_ref = point
    else:  # enhancer
        padded_start = max(0, r["start"] - flank_bp)
        padded_end = r["end"] + flank_bp
        target_ref = (r["start"] + r["end"]) // 2  # midpoint
    if chrom_len is not None:
        padded_end = min(padded_end, chrom_len)
    return padded_start, padded_end, target_ref


def load_chrom_sizes(fai_path):
    sizes = {}
    with open(fai_path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            sizes[fields[0]] = int(fields[1])
    return sizes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bed", required=True)
    ap.add_argument("--input_type", required=True, choices=["snp", "enhancer"])
    ap.add_argument("--fasta", required=True, help="reference FASTA (must have a .fai alongside it)")
    ap.add_argument("--flank_bp", type=int, default=50, help="enhancer-type padding, each side")
    ap.add_argument("--snp_window_bp", type=int, default=150, help="SNP-type search window, each side")
    ap.add_argument("--dedup_blast_radius", action="store_true")
    ap.add_argument("--link_distance_bp", type=int, default=1000)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    fai_path = args.fasta + ".fai"
    if not Path(fai_path).exists():
        sys.exit(f"ERROR: expected FASTA index at {fai_path} (run `samtools faidx {args.fasta}` first)")

    chrom_sizes = load_chrom_sizes(fai_path)
    regions = read_bed(args.bed, prefix=args.input_type)
    if not regions:
        sys.exit("ERROR: no regions parsed from BED file")

    linked = (
        flag_linked_regions(regions, args.link_distance_bp)
        if args.dedup_blast_radius
        else {r["name"]: [] for r in regions}
    )

    padded_bed_path = outdir / "regions.padded.bed"
    metadata_rows = []
    header_by_name = {}
    with open(padded_bed_path, "w") as pb:
        for r in regions:
            padded_start, padded_end, target_ref = compute_window(
                r, args.input_type, args.snp_window_bp, args.flank_bp, chrom_sizes)
            fasta_header = f"{r['name']}__{r['chrom']}:{padded_start}-{padded_end}"
            header_by_name[r["name"]] = fasta_header
            pb.write(f"{r['chrom']}\t{padded_start}\t{padded_end}\t{r['name']}\n")
            metadata_rows.append({
                "region_id": r["name"],
                "fasta_header": fasta_header,
                "chrom": r["chrom"],
                "orig_start": r["start"],
                "orig_end": r["end"],
                "input_type": args.input_type,
                "target_ref": target_ref,
                "padded_start": padded_start,
                "padded_end": padded_end,
                "linked_regions": ",".join(sorted(set(linked[r["name"]]))),
            })

    # Extract padded sequences, rewrite headers to region_id__chrom:start-end
    fasta_out = outdir / "regions.fasta"
    raw_fasta = outdir / "regions.raw.fasta"
    cmd = ["bedtools", "getfasta", "-fi", args.fasta, "-bed", str(padded_bed_path),
           "-name", "-fo", str(raw_fasta)]
    subprocess.run(cmd, check=True)

    with open(raw_fasta) as rf, open(fasta_out, "w") as fo:
        for line in rf:
            if line.startswith(">"):
                header = line[1:].strip()
                region_id = header.split("::")[0]  # bedtools -name: "region_id::chrom:start-end"
                fo.write(f">{header_by_name[region_id]}\n")
            else:
                fo.write(line)

    with open(outdir / "regions.metadata.tsv", "w", newline="") as mf:
        writer = csv.DictWriter(mf, fieldnames=[
            "region_id", "fasta_header", "chrom", "orig_start", "orig_end", "input_type",
            "target_ref", "padded_start", "padded_end", "linked_regions"], delimiter="\t")
        writer.writeheader()
        writer.writerows(metadata_rows)

    raw_fasta.unlink(missing_ok=True)
    print(f"Wrote {len(regions)} {args.input_type} regions -> {fasta_out}", file=sys.stderr)


if __name__ == "__main__":
    main()

