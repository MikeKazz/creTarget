#!/usr/bin/env python3
"""
coords_to_bed.py

Extracts unique "chr:start-end" style values from a column of a TSV/CSV
file and writes them out as a BED file (chrom, start, end, name), where
name is the original coordinate string. Useful for pulling real target
regions out of an existing results file (e.g. a `target_enhancer` column)
to use as --enhancer_bed / --snp_bed input.

Usage:
  python3 coords_to_bed.py --input file.tsv --column target_enhancer --out enhancers.bed
"""

import argparse
import csv
import re
import sys

COORD_RE = re.compile(r'^(\S+):(\d+)-(\d+)$')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--column", required=True, help="column name containing chr:start-end values")
    ap.add_argument("--out", required=True)
    ap.add_argument("--delimiter", default="\t")
    args = ap.parse_args()

    seen = set()
    rows_out = []
    with open(args.input) as f:
        reader = csv.DictReader(f, delimiter=args.delimiter)
        if args.column not in reader.fieldnames:
            sys.exit(f"ERROR: column '{args.column}' not found. Available columns: {reader.fieldnames}")
        for row in reader:
            val = row[args.column].strip()
            if val in seen or val in (".", "", "NA", "NaN"):
                continue
            seen.add(val)
            m = COORD_RE.match(val)
            if not m:
                print(f"WARNING: skipping unparseable value: {val!r}", file=sys.stderr)
                continue
            chrom, start, end = m.group(1), m.group(2), m.group(3)
            rows_out.append((chrom, start, end, val))

    with open(args.out, "w") as f:
        for chrom, start, end, name in rows_out:
            f.write(f"{chrom}\t{start}\t{end}\t{name}\n")

    print(f"Wrote {len(rows_out)} unique regions -> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
