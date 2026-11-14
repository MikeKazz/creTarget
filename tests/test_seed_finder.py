"""Unit tests for the seed-finder-lite integration.

Two things are covered here. First, the part creTarget owns: turning
seed-finder's hit-level output (zero to many rows per guide) into one summary
row per guide, now done by streaming because this step sees the whole collated
table rather than a single chunk.

Second, and more importantly, the install check. seed-finder-lite had a defect
where its plus-strand scan indexed guides under the reverse complement of
their seed, so it reported only hits whose protospacer lay on the minus
strand -- exiting 0, with plausible numbers, while missing about half of all
true hits. `TestRunCheckCatchesTheStrandBug` stubs out a tool with exactly
that behaviour and asserts creTarget refuses to proceed.
"""

import csv
import os
import random
import textwrap
from pathlib import Path

import pytest

from cretarget.steps.seed_finder import (
    SUMMARY_COLUMNS,
    HitSummary,
    build_guide_id,
    find_probe_guides,
    first_contig_sequence,
    load_hit_summary,
    reverse_complement,
    run_check,
    write_annotated,
    write_tool_input,
    _is_usable_probe,
)


def random_dna(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


class TestBuildGuideId:
    def test_composes_from_region_coordinates_and_orientation(self):
        assert build_guide_id("reg0", 100, 123, "FWD", 0) == "reg0_100_123_FWD_0"

    def test_distinguishes_strands_at_the_same_position(self):
        a = build_guide_id("r", 100, 123, "FWD", 0)
        b = build_guide_id("r", 100, 123, "RVS", 1)
        assert a != b

    def test_index_makes_true_duplicates_unique(self):
        """Two guides can genuinely share region, span and orientation.
        seed-finder-lite rejects duplicate ids outright, so the index is
        always appended rather than only on collision -- which also means an
        id can be assigned while streaming, with no global dedup pass."""
        ids = [build_guide_id("r", 100, 123, "FWD", i) for i in range(3)]
        assert len(set(ids)) == 3


class TestHitSummaryNoHits:
    def test_unknown_guide_reports_zero_not_an_error(self):
        assert HitSummary().row_for("nope") == ["0", "", "", ""]

    def test_row_width_matches_the_column_list(self):
        assert len(HitSummary().row_for("nope")) == len(SUMMARY_COLUMNS)

    def test_guides_without_hits_are_not_stored(self):
        """This is what keeps memory bounded on a genome-wide run, where most
        guides have no promoter-proximal seed match at all."""
        summary = HitSummary()
        summary.add("g1", "GENEA", 10, 9)
        assert len(summary) == 1
        summary.row_for("g2")
        assert len(summary) == 1


class TestHitSummaryWithHits:
    def test_counts_hits_per_guide(self):
        summary = HitSummary()
        summary.add("g1", "GENEA", 100, 10)
        summary.add("g1", "GENEB", 400, 9)
        summary.add("g2", "GENEC", 50, 12)
        assert summary.row_for("g1")[0] == "2"
        assert summary.row_for("g2")[0] == "1"
        assert summary.row_for("g3")[0] == "0"
        assert summary.total_hits == 3

    def test_picks_the_closest_tss(self):
        summary = HitSummary()
        summary.add("g1", "FAR", 900, 8)
        summary.add("g1", "NEAR", 25, 11)
        summary.add("g1", "MID", 300, 9)
        count, gene, distance, seed_length = summary.row_for("g1")
        assert (gene, distance, seed_length) == ("NEAR", "25", "11")

    def test_closest_uses_absolute_distance(self):
        # Unsigned in current seed-finder-lite output, but a hit 30 bp
        # upstream is closer than one 500 bp downstream if that ever changes.
        summary = HitSummary()
        summary.add("g1", "DOWNSTREAM", 500, 8)
        summary.add("g1", "UPSTREAM", -30, 10)
        assert summary.row_for("g1")[1] == "UPSTREAM"

    def test_ties_keep_the_first_hit_seen(self):
        summary = HitSummary()
        summary.add("g1", "FIRST", 100, 9)
        summary.add("g1", "SECOND", 100, 12)
        assert summary.row_for("g1")[1] == "FIRST"

    def test_zero_distance_is_not_blanked(self):
        """A hit exactly on the TSS has distance 0, which a truthiness test
        would turn into an empty field."""
        summary = HitSummary()
        summary.add("g1", "GENEA", 0, 20)
        assert summary.row_for("g1") == ["1", "GENEA", "0", "20"]

    def test_unparseable_distance_does_not_crash(self):
        summary = HitSummary()
        summary.add("g1", "GENEA", "", 9)
        summary.add("g1", "GENEB", 40, 10)
        assert summary.row_for("g1")[0] == "2"
        # The row with a usable distance wins over the one without.
        assert summary.row_for("g1")[1] == "GENEB"


def ranked_table(tmp_path, rows, orientation=True):
    header = ["region_id", "genome_start", "genome_end", "spacer", "rank"]
    if orientation:
        header.insert(3, "orientation")
    path = tmp_path / "guides.ranked.tsv"
    with open(path, "w", newline="") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join(str(v) for v in row) + "\n")
    return path


class TestWriteToolInput:
    def test_writes_two_columns_with_unique_ids(self, tmp_path):
        ranked = ranked_table(tmp_path, [
            ("reg1", 100, 119, "FWD", "A" * 20, 1),
            ("reg1", 100, 119, "FWD", "C" * 20, 2),
        ])
        out = tmp_path / "in.tsv"
        assert write_tool_input(ranked, out) == 2
        lines = out.read_text().splitlines()
        assert lines[0] == "guide_id\tsequence"
        ids = [line.split("\t")[0] for line in lines[1:]]
        assert len(set(ids)) == 2
        assert lines[1].split("\t")[1] == "A" * 20

    def test_works_without_an_orientation_column(self, tmp_path):
        ranked = ranked_table(
            tmp_path, [("reg1", 100, 119, "G" * 20, 1)], orientation=False)
        out = tmp_path / "in.tsv"
        write_tool_input(ranked, out)
        assert "_NA_0\t" in out.read_text()

    def test_missing_required_column_is_reported_not_guessed(self, tmp_path):
        path = tmp_path / "bad.tsv"
        path.write_text("region_id\tspacer\nreg1\tAAAA\n")
        with pytest.raises(SystemExit) as excinfo:
            write_tool_input(path, tmp_path / "in.tsv")
        assert "genome_start" in str(excinfo.value)


class TestLoadHitSummary:
    def test_absent_file_is_empty_not_fatal(self, tmp_path):
        assert len(load_hit_summary(tmp_path / "nope.csv")) == 0

    def test_header_only_file_is_empty(self, tmp_path):
        path = tmp_path / "hits.csv"
        path.write_text("guide_id,nearest_tss_gene,tss_distance,seed_length\n")
        assert len(load_hit_summary(path)) == 0

    def test_completely_empty_file_is_empty(self, tmp_path):
        path = tmp_path / "hits.csv"
        path.write_text("")
        assert len(load_hit_summary(path)) == 0

    def test_aggregates_rows(self, tmp_path):
        path = tmp_path / "hits.csv"
        path.write_text(
            "guide_id,nearest_tss_gene,tss_distance,seed_length\n"
            "g1,GENEA,500,8\n"
            "g1,GENEB,20,11\n"
        )
        summary = load_hit_summary(path)
        assert summary.row_for("g1") == ["2", "GENEB", "20", "11"]


class TestWriteAnnotated:
    def test_appends_summary_columns_to_every_row(self, tmp_path):
        ranked = ranked_table(tmp_path, [
            ("reg1", 100, 119, "FWD", "A" * 20, 1),
            ("reg2", 200, 219, "RVS", "C" * 20, 1),
        ])
        summary = HitSummary()
        summary.add(build_guide_id("reg1", 100, 119, "FWD", 0), "GENEA", 42, 9)

        out = tmp_path / "final.tsv"
        assert write_annotated(ranked, out, summary) == 2
        rows = list(csv.DictReader(open(out), delimiter="\t"))
        assert [c in rows[0] for c in SUMMARY_COLUMNS] == [True] * 4
        assert rows[0]["seedfinder_hit_count"] == "1"
        assert rows[0]["seedfinder_closest_gene"] == "GENEA"
        assert rows[1]["seedfinder_hit_count"] == "0"
        assert rows[1]["seedfinder_closest_gene"] == ""

    def test_ids_match_between_the_two_streaming_passes(self, tmp_path):
        """The guide_id is regenerated in the second pass rather than held in
        memory, so the two passes must agree or the join silently misses."""
        rows = [(f"reg{i}", i, i + 19, "FWD", "A" * 20, 1) for i in range(50)]
        ranked = ranked_table(tmp_path, rows)
        tool_input = tmp_path / "in.tsv"
        write_tool_input(ranked, tool_input)
        ids_from_pass_one = [
            line.split("\t")[0] for line in tool_input.read_text().splitlines()[1:]
        ]

        summary = HitSummary()
        for guide_id in ids_from_pass_one:
            summary.add(guide_id, "GENE", 1, 20)

        out = tmp_path / "final.tsv"
        write_annotated(ranked, out, summary)
        counts = [r["seedfinder_hit_count"]
                  for r in csv.DictReader(open(out), delimiter="\t")]
        assert counts == ["1"] * 50

    def test_original_columns_are_preserved(self, tmp_path):
        ranked = ranked_table(tmp_path, [("reg1", 100, 119, "FWD", "A" * 20, 7)])
        out = tmp_path / "final.tsv"
        write_annotated(ranked, out, HitSummary())
        row = next(csv.DictReader(open(out), delimiter="\t"))
        assert row["region_id"] == "reg1" and row["rank"] == "7"


class TestProbeSelection:
    def test_rejects_low_complexity(self):
        assert not _is_usable_probe("A" * 20)
        assert not _is_usable_probe("AC" * 10)          # only 2 distinct bases
        assert not _is_usable_probe("ACGT" + "G" * 7 + "ACGTACGTA")  # long run

    def test_accepts_ordinary_sequence(self):
        assert _is_usable_probe("ACGTTGCAACGTGGCAATCG")

    def test_finds_real_sites_in_both_orientations(self):
        sequence = random_dna(3000, seed=11)
        plus, minus = find_probe_guides(sequence)
        assert plus and minus

        # The plus probe must sit immediately 5' of a real NGG.
        i = sequence.index(plus)
        assert sequence[i + 21:i + 23] == "GG"

        # The minus probe's protospacer is the reverse complement of a 20-mer
        # that follows a real CCN.
        j = sequence.index(reverse_complement(minus))
        assert sequence[j - 3:j - 1] == "CC"

    def test_reads_only_the_first_contig(self, tmp_path):
        path = tmp_path / "prom.fa"
        path.write_text(">prom|chr1|1|10\nACGTACGTAC\nGTACGTACGT\n>prom|chr2|1|4\nTTTT\n")
        assert first_contig_sequence(path) == "ACGTACGTACGTACGTACGT"


# ---------------------------------------------------------------------------
# The install check
# ---------------------------------------------------------------------------

STUB_TOOL = '''
"""A stand-in for seed-finder-lite, to exercise creTarget's install check."""
import csv, sys

def main():
    argv = sys.argv[1:]
    assert argv[0] == "seed-find-batch", argv
    input_tsv = argv[1]
    out = argv[argv.index("-o") + 1]

    guides = []
    with open(input_tsv) as handle:
        next(handle)
        for line in handle:
            if line.strip():
                gid, seq = line.rstrip("\\n").split("\\t")
                guides.append((gid, seq))

    # MODE is baked in when the stub is written.
    mode = "__MODE__"
    cols = ["guide_id", "chr", "start", "end", "strand", "seed_length",
            "nearest_tss_gene", "tss_distance", "pam", "protospacer",
            "guide_sequence", "ucsc_url"]
    with open(out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=cols)
        writer.writeheader()
        for gid, seq in guides:
            if mode == "minus_only" and gid.endswith("_PLUS"):
                continue          # the real bug: plus-strand loop finds nothing
            writer.writerow({
                "guide_id": gid, "chr": "chr1", "start": 1, "end": 20,
                "strand": "+", "seed_length": 20, "nearest_tss_gene": "GENEA",
                "tss_distance": 10, "pam": "AGG", "protospacer": seq,
                "guide_sequence": seq, "ucsc_url": "",
            })

main()
'''


def make_stub_install(tmp_path, mode):
    """A directory shaped like seed-finder-lite, with a stub tool."""
    root = tmp_path / f"sf_{mode}"
    (root / "src" / "targetfinder").mkdir(parents=True)
    (root / "refs").mkdir(parents=True)

    (root / "src" / "targetfinder" / "__init__.py").write_text("")
    (root / "src" / "targetfinder" / "__main__.py").write_text(
        STUB_TOOL.replace("__MODE__", mode))

    sequence = random_dna(3000, seed=11)
    wrapped = "\n".join(textwrap.wrap(sequence, 60))
    (root / "refs" / "promoter_trim_2000_pad23.fa.gz").write_text("")  # placeholder
    (root / "refs" / "promoter_trim_2000_pad23.fa.gz").unlink()
    # Written uncompressed under a .fa name so the helper's gzip branch is not
    # taken; _resolve_tool_paths only requires the .gz path to exist, so both
    # are provided and the .gz is the one it reads.
    import gzip
    with gzip.open(root / "refs" / "promoter_trim_2000_pad23.fa.gz", "wt") as handle:
        handle.write(f">prom|chr1|1|{len(sequence)}\n{wrapped}\n")
    (root / "refs" / "tss_table.tsv").write_text("chr\ttss\tgene_name\nchr1\t500\tGENEA\n")
    return root


class _CheckArgs:
    def __init__(self, seed_finder_dir, outdir):
        self.seed_finder_dir = str(seed_finder_dir)
        self.python_bin = None
        self.outdir = str(outdir)


class TestRunCheckCatchesTheStrandBug:
    def test_passes_when_both_orientations_are_found(self, tmp_path):
        root = make_stub_install(tmp_path, "both")
        assert run_check(_CheckArgs(root, tmp_path / "chk_ok")) == 0

    def test_fails_when_only_minus_strand_hits_come_back(self, tmp_path):
        """The regression that motivated this check: the tool succeeds, returns
        hits, and is wrong."""
        root = make_stub_install(tmp_path, "minus_only")
        with pytest.raises(SystemExit) as excinfo:
            run_check(_CheckArgs(root, tmp_path / "chk_bad"))
        message = str(excinfo.value)
        assert "both orientations" in message
        assert "NOT FOUND" in message
        assert "--seed_finder_skip_check" in message

    def test_missing_install_is_reported_clearly(self, tmp_path):
        with pytest.raises(SystemExit) as excinfo:
            run_check(_CheckArgs(tmp_path / "nothing_here", tmp_path / "chk"))
        assert "seed-finder-lite" in str(excinfo.value)
