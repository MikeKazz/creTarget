"""
coords_to_bed

Pulls unique "chr:start-end" values out of a column of a TSV/CSV and writes
them as a BED file. Handy for turning an existing results file (e.g. a
`target_enhancer` column) back into pipeline input.
"""

import argparse
import csv
import re
import sys

COORD_RE = re.compile(r"^(\S+):([\d,]+)-([\d,]+)$")
NULL_VALUES = {".", "", "NA", "NaN", "nan", "None"}


def run(args):
    seen = set()
    rows = []

    with open(args.input, newline="") as fh:
        reader = csv.DictReader(fh, delimiter=args.delimiter)
        if reader.fieldnames is None or args.column not in reader.fieldnames:
            raise SystemExit(f"ERROR: column {args.column!r} not found. "
                             f"Available: {reader.fieldnames}")
        for row in reader:
            value = (row[args.column] or "").strip()
            if value in NULL_VALUES or value in seen:
                continue
            seen.add(value)
            match = COORD_RE.match(value)
            if not match:
                print(f"WARNING: skipping unparseable value: {value!r}", file=sys.stderr)
                continue
            chrom = match.group(1)
            start = match.group(2).replace(",", "")
            end = match.group(3).replace(",", "")
            rows.append((chrom, start, end, value))

    with open(args.out, "w") as out:
        for chrom, start, end, name in rows:
            out.write(f"{chrom}\t{start}\t{end}\t{name}\n")

    print(f"Wrote {len(rows)} unique regions -> {args.out}", file=sys.stderr)
    return args.out


def add_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--input", required=True)
    parser.add_argument("--column", required=True,
                        help="column holding chr:start-end values")
    parser.add_argument("--out", required=True)
    parser.add_argument("--delimiter", default="\t")
