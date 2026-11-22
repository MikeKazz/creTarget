"""
seed_finder

Promoter-proximal seed screening: for every candidate guide, find the places
where its PAM-proximal seed matches perfectly next to an NGG within a window
around an annotated TSS. A hit means plausible off-target *transcriptional
repression*, which is a different failure mode from the off-target *cutting*
that FlashFry's scores describe.

The scan itself lives in cretarget.seed_scan, which documents the method and
credits its origin: Hartman et al. 2026 (doi:10.64898/2026.03.27.714658) for
the method, and seed-finder-lite (doi:10.5281/zenodo.21599135) for the
promoter-window reference design and the single-pass architecture. This module is the plumbing: decide which guides to screen,
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
exists. Measured at ~1 KB of RSS per guide, so budget ~1 GB per million
guides screened.
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
    # Whether the guide's OWN target overlaps a promoter window. Not an
    # off-target signal -- but it is worth knowing, and keeping it separate is
    # what lets hit_count mean off-target risk alone.
    "seedfinder_target_is_promoter",
]

HIT_COLUMNS = [
    # self_hit first so it is impossible to miss when eyeballing the table.
    "self_hit",
    "guide_id", "chr", "start", "end", "strand", "seed_length",
    # nearest_tss_gene / tss_distance describe the CLOSEST TSS; tss_distance is
    # signed (negative upstream, positive downstream, relative to that gene's
    # strand) to match the convention of Hartman et al.
    "nearest_tss_gene", "tss_distance", "tss_strand",
    # Every TSS within the window, nearest first. One site can sit in two
    # promoters at once -- a bidirectional promoter puts two genes on opposite
    # strands a few hundred bp apart -- and reporting only the closest would
    # drop the other gene silently.
    "tss_genes_in_window", "tss_distances_in_window",
    "pam", "pam_start", "pam_end",
    "protospacer", "guide_sequence", "ucsc_url",
]

REQUIRED_COLUMNS = ("spacer", "region_id", "genome_start", "genome_end")
# Needed to recognise a self-hit. Absent, screening still works and every hit
# is reported as a potential off-target.
LOCUS_COLUMN = "chrom"

# Above this fraction of screened guides with unreadable coordinates, the
# screen refuses to emit results rather than publishing inflated off-target
# counts. Mirrors the >1% ceiling used for the metadata join in annotate.py.
UNPARSEABLE_LOCUS_FATAL_FRACTION = 0.01

# csv's default field limit is 128 KB. Nothing here should approach it, but a
# malformed line should not abort a long run with an opaque error.
csv.field_size_limit(10 * 1024 * 1024)


def build_guide_id(region_id, genome_start, genome_end, orientation, index,
                   prefix=None):
    """A guide_id that is unique by construction.

    Region plus span plus orientation is almost unique but not quite -- two
    guides can genuinely share all three -- so the row index is always
    appended rather than only on collision. That also lets the id be assigned
    while streaming, with no global de-duplication pass over a table that may
    not fit in memory.

    `prefix` exists for sharded runs. The index is the row's position in the
    file being read, so with the table split across tasks it restarts at 0 in
    every shard: two genuinely duplicated guides in different shards would
    then receive the SAME id, and the merged hit table could not say which
    guide a hit belonged to. (The per-guide summary stays correct either way,
    because each shard joins against its own slice -- only the hit table is
    ambiguous.) Passing the shard's name disambiguates it. Unsharded runs
    leave it unset so ids stay short.
    """
    base = f"{region_id}_{genome_start}_{genome_end}_{orientation}_{index}"
    return f"{prefix}_{base}" if prefix else base


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

    def add(self, guide_id, gene, tss_distance, seed_length, self_hit=False):
        """Record one hit.

        A self-hit -- the guide matching its own target site, which happens
        whenever the targeted region overlaps a promoter window -- is recorded
        as a flag and deliberately NOT counted. It is a perfect full-length
        match, so counting it would make every promoter-overlapping target
        look like a maximal off-target risk. On a genome-wide cCRE library
        these were 73% of all hits.
        """
        entry = self._by_guide.get(guide_id)
        if entry is None:
            entry = [0, None, "", "", "", False]
            self._by_guide[guide_id] = entry

        if self_hit:
            entry[5] = True
            return

        try:
            # Chosen on magnitude, stored signed: "closest" means least
            # distant either way, but the sign is what makes the value
            # comparable with a published distance/repression curve.
            key = abs(int(tss_distance))
        except (TypeError, ValueError):
            key = None
        entry[0] += 1
        # Strictly-less keeps the first hit seen on a tie, which makes the
        # choice deterministic given the scan's deterministic order.
        if key is not None and (entry[1] is None or key < entry[1]):
            entry[1:5] = [key, gene, tss_distance, seed_length]

    def row_for(self, guide_id, screened=True):
        """The four summary values, in SUMMARY_COLUMNS order, always strings.

        A guide that was not screened gets empty strings rather than a count
        of 0, so "we looked and found nothing" stays distinguishable from "we
        did not look".
        """
        if not screened:
            return ["", "", "", "", ""]
        entry = self._by_guide.get(guide_id)
        if entry is None:
            return ["0", "", "", "", "no"]
        return [str(entry[0]), _as_field(entry[2]), _as_field(entry[3]),
                _as_field(entry[4]), "yes" if entry[5] else "no"]

    def __len__(self):
        return len(self._by_guide)

    @property
    def total_hits(self):
        """Off-target hits only; self-hits are excluded by construction."""
        return sum(entry[0] for entry in self._by_guide.values())

    @property
    def guides_with_hits(self):
        return sum(1 for entry in self._by_guide.values() if entry[0])

    @property
    def self_hits(self):
        return sum(1 for entry in self._by_guide.values() if entry[5])


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


def iter_guides(ranked_path, passing_only=True, id_prefix=None):
    """Stream (guide_id, spacer, locus) for the guides to screen.

    `locus` is (chrom, genome_start, genome_end) -- the guide's own target
    site, which the scanner uses to recognise a self-hit. It is None when the
    table has no `chrom` column, in which case every hit is reported as a
    potential off-target.
    """
    with open(ranked_path, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        idx = _header_indices(header, ranked_path)
        orientation_index = (header.index("orientation")
                             if "orientation" in header else None)
        filter_index = (header.index("passes_hard_filters")
                        if passing_only and "passes_hard_filters" in header
                        else None)
        chrom_index = (header.index(LOCUS_COLUMN)
                       if LOCUS_COLUMN in header else None)
        if chrom_index is None:
            print(f"WARNING: {ranked_path} has no {LOCUS_COLUMN!r} column, so a "
                  f"guide cannot be recognised when it finds its own target "
                  f"site. Self-hits will be counted as off-targets.",
                  file=sys.stderr)
        width = max(idx.values())
        unparseable_loci = 0
        screened_rows = 0
        for index, row in enumerate(reader):
            # csv.reader yields [] for a blank line, and a truncated row for a
            # table whose writer was killed mid-line. Both used to raise
            # IndexError below -- after the promoterome scan had already run,
            # which on a genome-wide table is hours in.
            if len(row) <= width:
                if row:
                    print(f"WARNING: skipping short row {index + 2} of "
                          f"{ranked_path} ({len(row)} fields)", file=sys.stderr)
                continue
            if not _is_screened(row, filter_index):
                continue
            locus = None
            if chrom_index is not None:
                try:
                    # int(float(x)) so "119869.0" parses. A table written by
                    # an older creTarget has float coordinates, and refusing
                    # to read them would disable self-hit detection for the
                    # whole run -- the exact failure this guard exists to
                    # catch. Coordinates are integral by construction, so
                    # the float round-trip is lossless below 2^53.
                    locus = (sys.intern(row[chrom_index]),
                             int(float(row[idx["genome_start"]])),
                             int(float(row[idx["genome_end"]])))
                except (ValueError, IndexError, TypeError):
                    locus = None
                    unparseable_loci += 1
                    screened_rows += 1
                    if unparseable_loci == 1:
                        print(f"WARNING: cannot parse genome coordinates for "
                              f"{row[idx['region_id']]!r} "
                              f"({row[idx['genome_start']]!r}, "
                              f"{row[idx['genome_end']]!r}); self-hit "
                              f"detection is disabled for affected guides.",
                              file=sys.stderr)
                else:
                    screened_rows += 1
            yield build_guide_id(
                row[idx["region_id"]], row[idx["genome_start"]],
                row[idx["genome_end"]],
                row[orientation_index] if orientation_index is not None else "NA",
                index, id_prefix,
            ), row[idx["spacer"]].strip().upper(), locus

    # Reported as a TOTAL, and fatal above a threshold.
    #
    # Warning once on the first bad row was too weak for the highest-severity
    # failure mode here: a null locus disables self-hit detection, self-hits
    # were 73% of all hits in a genome-wide run, and seedfinder_hit_count is
    # overstated accordingly while the run exits 0. On a sharded run that
    # single warning is one line in one of N .command.err files. This matches
    # the standard set elsewhere in this release, where >1% unmatched on the
    # metadata join is fatal.
    if unparseable_loci:
        fraction = unparseable_loci / max(screened_rows, 1)
        message = (f"{unparseable_loci:,} of {screened_rows:,} screened "
                   f"guide(s) ({fraction:.1%}) in {ranked_path} have "
                   f"unparseable genome coordinates, so they cannot be "
                   f"recognised when they find their own target site.")
        if fraction > UNPARSEABLE_LOCUS_FATAL_FRACTION:
            raise SystemExit(
                f"ERROR: {message}\n"
                f"seedfinder_hit_count would be substantially overstated for "
                f"them -- self-hits dominate the raw hit count. Refusing to "
                f"emit a screen whose off-target counts are not trustworthy.\n"
                f"Check that genome_start/genome_end are integers in the "
                f"ranked table."
            )
        print(f"WARNING: {message} Their seedfinder_hit_count will be "
              f"overstated.", file=sys.stderr)


def write_annotated(ranked_path, out_path, summary, passing_only=True,
                    id_prefix=None):
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

        # A csv.writer, to match the csv.reader above. Reading strips any
        # quoting and "\t".join put none back, so a quoted field containing a
        # tab would have been re-emitted as two columns and shifted every
        # later field in that row.
        writer = csv.writer(dst, delimiter="\t", lineterminator="\n")
        writer.writerow(header + SUMMARY_COLUMNS)
        width = max(idx.values())
        skipped = 0
        for index, row in enumerate(reader):
            # Mirrors the skip in iter_guides, so the two passes stay in step.
            if len(row) <= width:
                skipped += 1
                continue
            guide_id = build_guide_id(
                row[idx["region_id"]], row[idx["genome_start"]],
                row[idx["genome_end"]],
                row[orientation_index] if orientation_index is not None else "NA",
                index, id_prefix,
            )
            writer.writerow(
                row + summary.row_for(guide_id, _is_screened(row, filter_index)))
            written += 1

    # This pass writes the deliverable, so a row dropped here vanishes from
    # guides.final.tsv entirely. iter_guides warns about short rows; this one
    # used to drop them with no record anywhere, and nobody diffs row counts
    # between the ranked and final tables by hand on a 92M-guide run.
    if skipped:
        print(f"WARNING: {skipped} malformed row(s) in {ranked_path} are ABSENT "
              f"from {out_path}, which has {written} rows against the input's "
              f"{written + skipped}.", file=sys.stderr)
    return written


def scan_to_file(ranked_path, promoter_fa, tss_table, min_seed_length,
                 tss_window_bp, raw_hits_path, passing_only=True,
                 progress_every=5000, id_prefix=None):
    """One pass over the reference; hits streamed to CSV, summary accumulated."""
    guides = list(iter_guides(ranked_path, passing_only=passing_only,
                              id_prefix=id_prefix))
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
        # lineterminator is explicit because csv defaults to CRLF per RFC 4180.
        # That is legal CSV, but it leaves a trailing \r on the last field for
        # anyone reaching for awk or cut -- which is most people, on a table
        # this size -- so the last column silently reads as "...url\r".
        writer = csv.DictWriter(handle, fieldnames=HIT_COLUMNS,
                                extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        written = self_written = 0
        try:
            for hit in hits:
                hit["ucsc_url"] = seed_scan.ucsc_url(
                    hit["chr"], hit["start"], hit["end"])
                writer.writerow(hit)
                is_self = hit["self_hit"] == "yes"
                self_written += is_self
                written += 1
                summary.add(hit["guide_id"], hit["nearest_tss_gene"],
                            hit["tss_distance"], hit["seed_length"],
                            self_hit=is_self)
        except seed_scan.InvalidGuide as exc:
            raise SystemExit(f"ERROR: {exc}")

    print(f"{written:,} hit(s) written", file=sys.stderr)
    if self_written:
        print(f"  {self_written:,} ({100*self_written/written:.1f}%) are "
              f"SELF-hits: the guide matching its own target site, because "
              f"that target overlaps a promoter window. Recorded in the hit "
              f"table and flagged by seedfinder_target_is_promoter, but "
              f"excluded from seedfinder_hit_count -- they are not off-target "
              f"risk.", file=sys.stderr)
    print(f"  {summary.total_hits:,} off-target hit(s) across "
          f"{summary.guides_with_hits:,} of {len(guides):,} screened guides",
          file=sys.stderr)
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
        args.tss_window, raw_hits, passing_only=not args.all_guides,
        id_prefix=args.id_prefix)

    written = write_annotated(args.ranked, args.out, summary,
                             passing_only=not args.all_guides,
                             id_prefix=args.id_prefix)
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
    parser.add_argument("--id-prefix", dest="id_prefix", default=None,
                        help="prepended to every guide_id. Set this when the "
                             "table is one shard of a larger one, so ids stay "
                             "unique across shards.")
    parser.add_argument("--outdir", default="seed_finder")
    parser.add_argument("--raw-hits", dest="raw_hits", default=None)
    parser.add_argument("--out", default="guides.final.tsv")
