"""
seed_finder

Optional final annotation: runs the external seed-finder-lite batch CLI
(PAM-proximal seed / TSS-proximity search, Zenodo DOI 10.5281/zenodo.21599135)
over the ranked guide table and folds a per-guide summary back onto it.

WHY THIS RUNS ONCE, NOT PER CHUNK
---------------------------------
seed-finder-lite's batch mode scans the whole promoterome a single time and
matches every guide against it with a dictionary lookup, so its cost is
essentially independent of how many guides you hand it: measured at ~25 s for
20 guides and ~25 s for 2,000. Running it per pipeline chunk therefore pays
that fixed scan once per chunk -- at chunk_size 1000, a genome-wide run would
spend weeks re-reading the same 212 MB of sequence. One invocation over the
collated table costs one scan.

The consequence is that this step sees the entire guide table rather than a
single chunk, which can be tens of gigabytes. Nothing here loads it: the
table is streamed twice (once to build the tool's input, once to write the
annotated copy), and memory is bounded by the number of guides that actually
have hits rather than by the number of guides.

seed-finder-lite is neither bundled nor modified. This builds its input,
invokes it, and summarizes its hit-level output -- a guide can have zero to
many hits, so the full hit table is preserved separately while the main table
stays one row per guide.
"""

import argparse
import csv
import gzip
import os
import subprocess
import sys
from pathlib import Path

SUMMARY_COLUMNS = [
    "seedfinder_hit_count",
    "seedfinder_closest_gene",
    "seedfinder_closest_tss_distance",
    "seedfinder_closest_seed_length",
]

REQUIRED_COLUMNS = ("spacer", "region_id", "genome_start", "genome_end")

# csv's default field limit is 128 KB. Nothing here should approach that, but a
# malformed line would otherwise abort a long run with an opaque error.
csv.field_size_limit(10 * 1024 * 1024)


def build_guide_id(region_id, genome_start, genome_end, orientation, index):
    """A guide_id that is unique by construction.

    seed-finder-lite rejects duplicate guide_ids outright, and a duplicate
    would silently fan out the summary join besides. Region plus span plus
    orientation is almost unique but not quite -- two guides can genuinely
    share all three -- so the row index is always appended rather than only on
    collision. That also means the id can be assigned while streaming, with no
    global de-duplication pass over a table that may not fit in memory.
    """
    return f"{region_id}_{genome_start}_{genome_end}_{orientation}_{index}"


def _as_field(value):
    """Empty string for a genuinely absent value, but never for a falsy one --
    a tss_distance of 0 means the hit sits exactly on the TSS."""
    if value is None or value == "":
        return ""
    return str(value)


class HitSummary:
    """Per-guide aggregation of seed-finder hits, built incrementally.

    Holds one entry per guide that has at least one hit; guides with none are
    never stored, which is what keeps this bounded on a genome-wide run where
    most guides have no promoter-proximal seed match.
    """

    __slots__ = ("_by_guide",)

    def __init__(self):
        self._by_guide = {}

    def add(self, guide_id, gene, tss_distance, seed_length):
        """Record one hit.

        The 'closest' hit is the one with the smallest absolute TSS distance.
        tss_distance is unsigned in current seed-finder-lite output, but taking
        the absolute value keeps this correct if it ever becomes signed to
        indicate upstream/downstream.
        """
        try:
            distance = int(tss_distance)
        except (TypeError, ValueError):
            distance = None
        key = abs(distance) if distance is not None else None

        entry = self._by_guide.get(guide_id)
        if entry is None:
            self._by_guide[guide_id] = [1, key, gene, tss_distance, seed_length]
            return

        entry[0] += 1
        # Strictly-less keeps the first hit seen on a tie, matching the
        # first-occurrence behaviour of the previous pandas idxmin.
        if key is not None and (entry[1] is None or key < entry[1]):
            entry[1] = key
            entry[2] = gene
            entry[3] = tss_distance
            entry[4] = seed_length

    def row_for(self, guide_id):
        """The four summary values for one guide, in SUMMARY_COLUMNS order.

        Always strings, because these are written straight into a TSV. A guide
        with no hits gets a count of 0 and empty strings -- not an absent
        key -- so every input row keeps the same shape.
        """
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


def open_maybe_gzip(path, mode="rt"):
    """seed-finder-lite ships its promoter reference gzipped."""
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def _resolve_tool_paths(seed_finder_dir):
    seed_finder_dir = Path(seed_finder_dir)
    paths = {
        "src": seed_finder_dir / "src",
        "promoter_fa": seed_finder_dir / "refs" / "promoter_trim_2000_pad23.fa.gz",
        "tss_table": seed_finder_dir / "refs" / "tss_table.tsv",
    }
    labels = {"src": "src/", "promoter_fa": "refs/promoter fasta",
              "tss_table": "refs/tss_table.tsv"}
    for key, path in paths.items():
        if not path.exists():
            raise SystemExit(
                f"ERROR: seed-finder-lite {labels[key]} not found at {path}\n"
                f"Check that --seed-finder-dir points at the seed-finder-lite "
                f"directory (the one containing src/ and refs/)."
            )
    return paths


def _tool_env(src_dir):
    env = dict(os.environ)
    env["PYTHONPATH"] = str(src_dir) + (
        os.pathsep + env["PYTHONPATH"] if "PYTHONPATH" in env else "")
    return env


def _run_tool(python_bin, src_dir, input_tsv, promoter_fa, tss_table,
              min_seed_length, tss_window, out_csv, log_path):
    cmd = [
        python_bin or sys.executable, "-m", "targetfinder", "seed-find-batch",
        str(input_tsv),
        "--promoter-fa", str(promoter_fa),
        "--tss-table", str(tss_table),
        "--min-seed-length", str(min_seed_length),
        "--tss-window", str(tss_window),
        "-o", str(out_csv),
    ]
    with open(log_path, "w") as log:
        log.write("$ " + " ".join(cmd) + "\n\n")
        log.flush()
        proc = subprocess.run(cmd, env=_tool_env(src_dir), stdout=log,
                              stderr=subprocess.STDOUT)
    return proc.returncode, cmd


# ---------------------------------------------------------------------------
# Self-test
#
# The bug this guards against: a seed-finder-lite build whose plus-strand-PAM
# scan indexes guides under the reverse complement of their seed, so it only
# ever reports hits whose protospacer lies on the minus strand. That halves
# recall silently -- the tool exits 0 and returns plausible-looking hits -- and
# no amount of checking creTarget's own code would catch it.
#
# The test is self-contained: take two real Cas9 sites straight out of
# seed-finder-lite's own promoter reference, one in each orientation, and ask
# the tool to find them. A correct build recovers both.
# ---------------------------------------------------------------------------

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def reverse_complement(seq):
    return seq.translate(_COMPLEMENT)[::-1]


def _is_usable_probe(seq):
    """Avoid low-complexity 20-mers, which match in thousands of places and
    make the check slow without making it stronger."""
    if len(set(seq)) < 3:
        return False
    run = 1
    for a, b in zip(seq, seq[1:]):
        run = run + 1 if a == b else 1
        if run > 6:
            return False
    return True


def first_contig_sequence(promoter_fa):
    """The first record of the promoter FASTA, as a string."""
    chunks = []
    with open_maybe_gzip(promoter_fa) as handle:
        seen_header = False
        for line in handle:
            if line.startswith(">"):
                if seen_header:
                    break
                seen_header = True
                continue
            chunks.append(line.strip())
    return "".join(chunks).upper()


def find_probe_guides(sequence):
    """One plus-strand and one minus-strand real protospacer from `sequence`.

    Plus strand: a 20-mer immediately followed by NGG.
    Minus strand: CCN followed by a 20-mer, whose protospacer is that 20-mer
    reverse-complemented.
    """
    plus = minus = None
    for i in range(len(sequence) - 23):
        window = sequence[i:i + 23]
        if "N" in window:
            continue
        if plus is None and window[21] == "G" and window[22] == "G":
            if _is_usable_probe(window[:20]):
                plus = window[:20]
        if minus is None and window[0] == "C" and window[1] == "C":
            candidate = reverse_complement(window[3:23])
            if _is_usable_probe(candidate):
                minus = candidate
        if plus and minus:
            break
    return plus, minus


def run_check(args):
    """Verify a seed-finder-lite installation recovers sites in both
    orientations. Exits non-zero with an explanation if it does not."""
    paths = _resolve_tool_paths(args.seed_finder_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("Reading a probe contig from seed-finder-lite's promoter reference...",
          file=sys.stderr)
    sequence = first_contig_sequence(paths["promoter_fa"])
    plus, minus = find_probe_guides(sequence)
    if not plus or not minus:
        raise SystemExit(
            "ERROR: could not find probe sites in the promoter reference's first "
            "contig. This suggests the reference file is truncated or not the "
            "expected promoter FASTA."
        )

    input_tsv = outdir / "selftest_input.tsv"
    with open(input_tsv, "w", newline="") as handle:
        handle.write("guide_id\tsequence\n")
        handle.write(f"SELFTEST_PLUS\t{plus}\n")
        handle.write(f"SELFTEST_MINUS\t{minus}\n")

    out_csv = outdir / "selftest_hits.csv"
    log_path = outdir / "selftest.log"
    print("Running seed-finder-lite on two known sites "
          "(one promoterome scan, ~30 s)...", file=sys.stderr)
    returncode, cmd = _run_tool(
        args.python_bin, paths["src"], input_tsv, paths["promoter_fa"],
        paths["tss_table"],
        # Seed length 20 means the whole guide must match, so a recovered hit
        # is unambiguous. The TSS window is opened wide because the point is
        # strand handling, not promoter proximity.
        min_seed_length=20, tss_window=10 ** 9,
        out_csv=out_csv, log_path=log_path,
    )
    if returncode != 0:
        raise SystemExit(
            f"ERROR: seed-finder-lite self-test failed to run (exit {returncode}).\n"
            f"Command: {' '.join(cmd)}\nSee {log_path}"
        )

    expected = {"SELFTEST_PLUS": plus, "SELFTEST_MINUS": minus}
    recovered = set()
    strands = set()
    if out_csv.exists():
        with open(out_csv, newline="") as handle:
            for row in csv.DictReader(handle):
                guide_id = row.get("guide_id", "")
                if row.get("protospacer") == expected.get(guide_id):
                    recovered.add(guide_id)
                if row.get("strand"):
                    strands.add(row["strand"])

    missing = sorted(set(expected) - recovered)
    if missing:
        raise SystemExit(
            "ERROR: this seed-finder-lite installation is not finding sites in "
            "both orientations.\n\n"
            f"  plus-strand probe  {plus}: "
            f"{'FOUND' if 'SELFTEST_PLUS' in recovered else 'NOT FOUND'}\n"
            f"  minus-strand probe {minus}: "
            f"{'FOUND' if 'SELFTEST_MINUS' in recovered else 'NOT FOUND'}\n"
            f"  strand values reported: {sorted(strands) or 'none'}\n\n"
            "Both probes are real Cas9 sites taken from seed-finder-lite's own\n"
            "promoter reference, so a correct installation finds both. Missing\n"
            "one orientation means roughly half of all true seed hits will be\n"
            "invisible, with no error reported.\n\n"
            "Known cause: in the batch scanner, guides are indexed only under\n"
            "reverse_complement(guide[-min_seed_length:]), and the plus-strand-PAM\n"
            "loop looks up that same key while scoring against guide_rc. It must\n"
            "index the forward seed as well and score against the guide.\n\n"
            f"Hit table: {out_csv}\nTool log:  {log_path}\n\n"
            "Pass --seed_finder_skip_check to run anyway."
        )

    print(f"OK: both orientations recovered (strands reported: "
          f"{', '.join(sorted(strands))}).", file=sys.stderr)
    return 0


# ---------------------------------------------------------------------------
# Main step
# ---------------------------------------------------------------------------

def write_tool_input(ranked_path, input_tsv):
    """Stream the ranked table into seed-finder-lite's two-column input.

    Returns (guide_count, guide_ids_path_is_deterministic). The ids are
    regenerated identically in the second pass rather than held in memory.
    """
    count = 0
    with open(ranked_path, newline="") as src, \
            open(input_tsv, "w", newline="") as dst:
        reader = csv.DictReader(src, delimiter="\t")
        missing = [c for c in REQUIRED_COLUMNS if c not in (reader.fieldnames or [])]
        if missing:
            raise SystemExit(
                f"ERROR: {ranked_path} is missing expected column(s): "
                f"{', '.join(missing)}"
            )
        has_orientation = "orientation" in (reader.fieldnames or [])
        dst.write("guide_id\tsequence\n")
        for index, row in enumerate(reader):
            guide_id = build_guide_id(
                row["region_id"], row["genome_start"], row["genome_end"],
                row["orientation"] if has_orientation else "NA",
                index,
            )
            dst.write(f"{guide_id}\t{row['spacer']}\n")
            count += 1
    return count


def load_hit_summary(raw_hits):
    """Aggregate seed-finder-lite's hit rows into a per-guide summary."""
    summary = HitSummary()
    if not Path(raw_hits).exists():
        return summary
    with open(raw_hits, newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            return summary          # zero hits writes a header-only or empty file
        for row in reader:
            summary.add(
                row.get("guide_id", ""),
                row.get("nearest_tss_gene", ""),
                row.get("tss_distance", ""),
                row.get("seed_length", ""),
            )
    return summary


def write_annotated(ranked_path, out_path, summary):
    """Second streaming pass: copy the table with the four summary columns."""
    written = 0
    with open(ranked_path, newline="") as src, open(out_path, "w", newline="") as dst:
        reader = csv.reader(src, delimiter="\t")
        header = next(reader)
        try:
            idx = {name: header.index(name) for name in REQUIRED_COLUMNS}
        except ValueError as exc:
            raise SystemExit(f"ERROR: unexpected header in {ranked_path}: {exc}")
        orientation_idx = header.index("orientation") if "orientation" in header else None

        dst.write("\t".join(header + SUMMARY_COLUMNS) + "\n")
        for index, row in enumerate(reader):
            guide_id = build_guide_id(
                row[idx["region_id"]], row[idx["genome_start"]],
                row[idx["genome_end"]],
                row[orientation_idx] if orientation_idx is not None else "NA",
                index,
            )
            dst.write("\t".join(row + summary.row_for(guide_id)) + "\n")
            written += 1
    return written


def run(args):
    paths = _resolve_tool_paths(args.seed_finder_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    input_tsv = outdir / "seed_finder_input.tsv"
    guide_count = write_tool_input(args.ranked, input_tsv)
    print(f"Prepared {guide_count:,} guides for seed-finder-lite "
          f"(min_seed_length={args.min_seed_length}, "
          f"tss_window={args.tss_window})", file=sys.stderr)

    if guide_count == 0:
        raise SystemExit(f"ERROR: no guide rows found in {args.ranked}")

    raw_hits = Path(args.raw_hits) if args.raw_hits else outdir / "raw_hits.csv"
    log_path = outdir / "seed_finder.log"
    print("Scanning the promoterome once for all guides...", file=sys.stderr)
    returncode, cmd = _run_tool(
        args.python_bin, paths["src"], input_tsv, paths["promoter_fa"],
        paths["tss_table"], args.min_seed_length, args.tss_window,
        raw_hits, log_path,
    )
    if returncode != 0:
        raise SystemExit(
            f"ERROR: seed-finder-lite failed (exit {returncode}).\n"
            f"Command: {' '.join(cmd)}\nSee {log_path}"
        )
    if not raw_hits.exists():
        raise SystemExit(f"ERROR: expected output not found: {raw_hits}")

    summary = load_hit_summary(raw_hits)
    print(f"{summary.total_hits:,} hit(s) across {len(summary):,} of "
          f"{guide_count:,} guides", file=sys.stderr)

    written = write_annotated(args.ranked, args.out, summary)
    print(f"Wrote {written:,} guides with seed-finder-lite columns -> {args.out}",
          file=sys.stderr)
    print(f"Hit-level detail retained at {raw_hits}", file=sys.stderr)
    return args.out


def add_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--ranked", required=True)
    parser.add_argument("--seed-finder-dir", dest="seed_finder_dir",
                        help="seed-finder-lite directory (containing src/ and refs/)")
    parser.add_argument("--min-seed-length", dest="min_seed_length", type=int,
                        required=True,
                        help="minimum exact-match seed length (e.g. 8). No default "
                             "-- choose deliberately, see seed-finder-lite's docs.")
    parser.add_argument("--tss-window", dest="tss_window", type=int, required=True,
                        help="max distance (bp) from a TSS to count as a hit "
                             "(e.g. 1000). No default -- choose deliberately.")
    parser.add_argument("--python-bin", dest="python_bin", default=None,
                        help="python with pyfaidx installed (default: current)")
    parser.add_argument("--outdir", default="seed_finder")
    parser.add_argument("--raw-hits", dest="raw_hits", default=None)
    parser.add_argument("--out", default="guides.final.tsv")


def add_check_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--seed-finder-dir", dest="seed_finder_dir", required=True,
                        help="seed-finder-lite directory (containing src/ and refs/)")
    parser.add_argument("--python-bin", dest="python_bin", default=None,
                        help="python with pyfaidx installed (default: current)")
    parser.add_argument("--outdir", default="seed_finder_check")
