"""
seed_finder

Promoter-proximal seed screening: for every candidate guide, find the places
where its PAM-proximal seed matches perfectly next to an NGG within a window
around an annotated TSS. A hit means plausible off-target *transcriptional
repression*, which is a different failure mode from the off-target *cutting*
that FlashFry's scores describe.

The scan itself lives in cretarget.seed_scan, which documents the method and
credits seed-finder-lite, where the approach and the promoter-window reference
design come from. This module is the plumbing: decide which guides to screen,
run one pass over the reference, write the hit table, and fold a per-guide
summary back onto the ranked table.

WHY THIS RUNS ONCE, NOT PER CHUNK
---------------------------------
The scan's cost is one pass over the promoter reference and barely moves with
guide count -- the per-position work is a dictionary lookup whether there are
twenty guides or two million. Running it per pipeline chunk paid that fixed
pass on every chunk; at chunk_size 1000 a genome-wide run would have spent
weeks re-reading the same ~212 Mb of sequence.

The consequence is that this step sees the whole guide table rather than one
chunk, and that table can be tens of gigabytes. It is never loaded: the table
is streamed twice, hits are streamed to disk as they are found, and the only
thing held in memory is the seed index -- which is why `--seed_finder_max_guides`
exists. Budget roughly 1 GB of RAM per million guides screened.
"""

import argparse
import csv
import sys
from pathlib import Path

from .. import seed_scan

SUMMARY_COLUMNS = [
    "seedfinder_hit_count",
    "seedfinder_closest_gene",
    "seedfinder_closest_tss_distance",
    "seedfinder_closest_seed_length",
]

HIT_COLUMNS = [
    "guide_id", "chr", "start", "end", "strand", "seed_length",
    "nearest_tss_gene", "tss_distance", "pam", "pam_start", "pam_end",
    "protospacer", "guide_sequence", "ucsc_url",
]

REQUIRED_COLUMNS = ("spacer", "region_id", "genome_start", "genome_end")

# csv's default field limit is 128 KB. Nothing here should approach it, but a
# malformed line should not abort a long run with an opaque error.
csv.field_size_limit(10 * 1024 * 1024)


def build_guide_id(region_id, genome_start, genome_end, orientation, index):
    """A guide_id that is unique by construction.

    Region plus span plus orientation is almost unique but not quite -- two
    guides can genuinely share all three -- so the row index is always
    appended rather than only on collision. That also lets the id be assigned
    while streaming, with no global de-duplication pass over a table that may
    not fit in memory.
    """
    return f"{region_id}_{genome_start}_{genome_end}_{orientation}_{index}"


def _as_field(value):
    """Empty string for a genuinely absent value, but never for a falsy one --
    a tss_distance of 0 means the hit sits exactly on the TSS."""
    if value is None or value == "":
        return ""
    return str(value)


class HitSummary:
    """Per-guide aggregation, built incrementally as hits stream past.

    Holds an entry only for guides that have at least one hit, which is what
    keeps this bounded on a genome-wide run where most guides have none.
    """

    __slots__ = ("_by_guide",)

    def __init__(self):
        self._by_guide = {}

    def add(self, guide_id, gene, tss_distance, seed_length):
        try:
            key = abs(int(tss_distance))
        except (TypeError, ValueError):
            key = None

        entry = self._by_guide.get(guide_id)
        if entry is None:
            self._by_guide[guide_id] = [1, key, gene, tss_distance, seed_length]
            return
        entry[0] += 1
        # Strictly-less keeps the first hit seen on a tie, which makes the
        # choice deterministic given the scan's deterministic order.
        if key is not None and (entry[1] is None or key < entry[1]):
            entry[1:] = [key, gene, tss_distance, seed_length]

    def row_for(self, guide_id, screened=True):
        """The four summary values, in SUMMARY_COLUMNS order, always strings.

        A guide that was not screened gets empty strings rather than a count
        of 0, so "we looked and found nothing" stays distinguishable from "we
        did not look".
        """
        if not screened:
            return ["", "", "", ""]
        entry = self._by_guide.get(guide_id)
        if entry is None:
            return ["0", "", "", ""]
        return [str(entry[0]), _as_field(entry[2]), _as_field(entry[3]),
                _as_field(entry[4])]

    def __len__(self):
        return len(self._by_guide)

    @property
    def total_hits(self):
        return sum(entry[0] for entry in self._by_guide.values())


def _header_indices(header, path):
    missing = [c for c in REQUIRED_COLUMNS if c not in header]
    if missing:
        raise SystemExit(
            f"ERROR: {path} is missing expected column(s): {', '.join(missing)}")
    return {name: header.index(name) for name in REQUIRED_COLUMNS}


def _is_screened(row, filter_index):
    """Whether this guide goes into the screen.

    Guides excluded by the hard filters are not screened by default: they are
    not going into anyone's library, and leaving them out of the seed index
    cuts both the memory and the work. A table with no `passes_hard_filters`
    column is screened in full.
    """
    if filter_index is None:
        return True
    return row[filter_index].strip().lower() in ("true", "1", "yes")


def iter_guides(ranked_path, passing_only=True):
    """Stream (guide_id, spacer) for the guides to screen."""
    with open(ranked_path, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        idx = _header_indices(header, ranked_path)
        orientation_index = (header.index("orientation")
                             if "orientation" in header else None)
        filter_index = (header.index("passes_hard_filters")
                        if passing_only and "passes_hard_filters" in header
                        else None)
        for index, row in enumerate(reader):
            if not _is_screened(row, filter_index):
                continue
            yield build_guide_id(
                row[idx["region_id"]], row[idx["genome_start"]],
                row[idx["genome_end"]],
                row[orientation_index] if orientation_index is not None else "NA",
                index,
            ), row[idx["spacer"]].strip().upper()


def write_annotated(ranked_path, out_path, summary, passing_only=True):
    """Second streaming pass: the table plus the four summary columns."""
    written = 0
    with open(ranked_path, newline="") as src, open(out_path, "w", newline="") as dst:
        reader = csv.reader(src, delimiter="\t")
        header = next(reader)
        idx = _header_indices(header, ranked_path)
        orientation_index = (header.index("orientation")
                             if "orientation" in header else None)
        filter_index = (header.index("passes_hard_filters")
                        if passing_only and "passes_hard_filters" in header
                        else None)

        dst.write("\t".join(header + SUMMARY_COLUMNS) + "\n")
        for index, row in enumerate(reader):
            guide_id = build_guide_id(
                row[idx["region_id"]], row[idx["genome_start"]],
                row[idx["genome_end"]],
                row[orientation_index] if orientation_index is not None else "NA",
                index,
            )
            dst.write("\t".join(
                row + summary.row_for(guide_id, _is_screened(row, filter_index))
            ) + "\n")
            written += 1
    return written


def scan_to_file(ranked_path, promoter_fa, tss_table, min_seed_length,
                 tss_window_bp, raw_hits_path, passing_only=True,
                 progress_every=5000):
    """One pass over the reference; hits streamed to CSV, summary accumulated."""
    guides = list(iter_guides(ranked_path, passing_only=passing_only))
    if not guides:
        raise SystemExit(
            f"ERROR: no guides to screen in {ranked_path}.\n"
            f"Every row was excluded by the hard filters. Pass "
            f"--all-guides to screen them anyway."
        )
    print(f"Screening {len(guides):,} guide(s) "
          f"(min_seed_length={min_seed_length}, tss_window={tss_window_bp} bp)",
          file=sys.stderr)
    if tss_window_bp > seed_scan.PROMOTER_HALF_WIDTH:
        print(f"  note: the promoter reference only extends "
              f"{seed_scan.PROMOTER_HALF_WIDTH} bp from each TSS, so a window "
              f"of {tss_window_bp} bp finds nothing beyond that.",
              file=sys.stderr)

    tss_rows = seed_scan.load_tss_table(tss_table)
    if not tss_rows:
        raise SystemExit(f"ERROR: no TSS records read from {tss_table}")
    tss_index = seed_scan.TssIndex(tss_rows)
    print(f"  {len(tss_index):,} TSSs indexed", file=sys.stderr)

    def report(window_index):
        if window_index and window_index % progress_every == 0:
            print(f"  ... {window_index:,} promoter windows scanned",
                  file=sys.stderr)

    summary = HitSummary()
    try:
        hits = seed_scan.iter_hits(
            promoter_fa, guides, tss_index, min_seed_length, tss_window_bp,
            on_window=report)
    except seed_scan.InvalidGuide as exc:
        raise SystemExit(f"ERROR: {exc}")

    Path(raw_hits_path).parent.mkdir(parents=True, exist_ok=True)
    with open(raw_hits_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=HIT_COLUMNS,
                                extrasaction="ignore")
        writer.writeheader()
        try:
            for hit in hits:
                hit["ucsc_url"] = seed_scan.ucsc_url(
                    hit["chr"], hit["start"], hit["end"])
                writer.writerow(hit)
                summary.add(hit["guide_id"], hit["nearest_tss_gene"],
                            hit["tss_distance"], hit["seed_length"])
        except seed_scan.InvalidGuide as exc:
            raise SystemExit(f"ERROR: {exc}")

    print(f"{summary.total_hits:,} hit(s) across {len(summary):,} of "
          f"{len(guides):,} screened guides", file=sys.stderr)
    return summary


def run(args):
    from .. import config as cfg
    from .. import fetch

    config = cfg.load_config()
    promoter_fa, tss_table = fetch.require_seed_reference(
        args.promoter_fa, args.tss_table, config)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    raw_hits = Path(args.raw_hits) if args.raw_hits else outdir / "raw_hits.csv"

    summary = scan_to_file(
        args.ranked, promoter_fa, tss_table, args.min_seed_length,
        args.tss_window, raw_hits, passing_only=not args.all_guides)

    written = write_annotated(args.ranked, args.out, summary,
                             passing_only=not args.all_guides)
    print(f"Wrote {written:,} guides with seed-screening columns -> {args.out}",
          file=sys.stderr)
    print(f"Hit-level detail retained at {raw_hits}", file=sys.stderr)
    return args.out


def add_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--ranked", required=True)
    parser.add_argument("--min-seed-length", dest="min_seed_length", type=int,
                        required=True,
                        help="minimum exact-match PAM-proximal seed length "
                             "(e.g. 8). No default -- choose deliberately.")
    parser.add_argument("--tss-window", dest="tss_window", type=int, required=True,
                        help="max distance (bp) from a TSS to count as a hit "
                             "(e.g. 1000). No default -- choose deliberately.")
    parser.add_argument("--promoter-fa", dest="promoter_fa", default=None,
                        help="promoter-window FASTA (default: the built "
                             "reference in the cache)")
    parser.add_argument("--tss-table", dest="tss_table", default=None,
                        help="TSS table (default: the built reference in the cache)")
    parser.add_argument("--all-guides", dest="all_guides", action="store_true",
                        help="screen every guide, including those excluded by "
                             "the hard filters (default: only guides that pass)")
    parser.add_argument("--outdir", default="seed_finder")
    parser.add_argument("--raw-hits", dest="raw_hits", default=None)
    parser.add_argument("--out", default="guides.final.tsv")
