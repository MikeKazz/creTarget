"""Tests for the seed-screening step: which guides get screened, and how
hit-level output collapses to one summary row per guide.

The scanner itself is tested in test_seed_scan.py. What is covered here is the
plumbing around it, which is where the table-shape invariants live: this step
sees the whole collated guide table rather than a single chunk, streams it
twice, and must regenerate identical guide_ids on both passes or the summary
join silently misses.
"""

import csv
import gzip

import pytest

from cretarget import seed_scan as ss
from cretarget.steps.seed_finder import (
    HIT_COLUMNS,
    SUMMARY_COLUMNS,
    HitSummary,
    build_guide_id,
    iter_guides,
    scan_to_file,
    write_annotated,
)


class TestBuildGuideId:
    def test_composes_from_region_coordinates_and_orientation(self):
        assert build_guide_id("reg0", 100, 123, "FWD", 0) == "reg0_100_123_FWD_0"

    def test_distinguishes_strands_at_the_same_position(self):
        assert (build_guide_id("r", 100, 123, "FWD", 0)
                != build_guide_id("r", 100, 123, "RVS", 1))

    def test_the_index_makes_true_duplicates_unique(self):
        """Two guides can genuinely share region, span and orientation. The
        scanner rejects duplicate ids, and a duplicate would fan out the
        summary join besides, so the index is always appended."""
        assert len({build_guide_id("r", 100, 123, "FWD", i) for i in range(3)}) == 3


class TestHitSummary:
    def test_unknown_guide_reports_zero(self):
        assert HitSummary().row_for("nope") == ["0", "", "", "", "no"]

    def test_unscreened_guide_is_blank_not_zero(self):
        """"We looked and found nothing" must stay distinguishable from "we
        did not look", because the two mean different things to a reader
        choosing guides."""
        assert HitSummary().row_for("nope", screened=False) == ["", "", "", "", ""]

    def test_row_width_matches_the_column_list(self):
        assert len(HitSummary().row_for("x")) == len(SUMMARY_COLUMNS)

    def test_guides_without_hits_are_not_stored(self):
        summary = HitSummary()
        summary.add("g1", "GENEA", 10, 9)
        summary.row_for("g2")
        assert len(summary) == 1

    def test_counts_hits_per_guide(self):
        summary = HitSummary()
        summary.add("g1", "GENEA", 100, 10)
        summary.add("g1", "GENEB", 400, 9)
        summary.add("g2", "GENEC", 50, 12)
        assert summary.row_for("g1")[0] == "2"
        assert summary.row_for("g2")[0] == "1"
        assert summary.total_hits == 3

    def test_picks_the_closest_tss(self):
        summary = HitSummary()
        summary.add("g1", "FAR", 900, 8)
        summary.add("g1", "NEAR", 25, 11)
        summary.add("g1", "MID", 300, 9)
        assert summary.row_for("g1") == ["3", "NEAR", "25", "11", "no"]

    def test_closest_uses_absolute_distance(self):
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
        assert summary.row_for("g1") == ["1", "GENEA", "0", "20", "no"]

    def test_unparseable_distance_does_not_crash(self):
        summary = HitSummary()
        summary.add("g1", "GENEA", "", 9)
        summary.add("g1", "GENEB", 40, 10)
        assert summary.row_for("g1")[0] == "2"
        assert summary.row_for("g1")[1] == "GENEB"


GUIDE_A = "GACCTGGAACTCAGTTCAGA"
GUIDE_B = "TTGACCATGGCATTCCAGGA"


def ranked_table(tmp_path, rows, orientation=True, filters=True,
                 chrom="chr1"):
    """(tmp_path is created if needed, so shard fixtures can use subdirs.)"""
    """A table shaped like guides.ranked.tsv.

    rows: (region_id, start, end, orientation, spacer, passes_hard_filters)
    """
    header = ["region_id", "chrom", "genome_start", "genome_end", "orientation",
              "spacer", "passes_hard_filters", "rank"]
    if not orientation:
        header.remove("orientation")
    if not filters:
        header.remove("passes_hard_filters")
    if chrom is None:
        header.remove("chrom")
    tmp_path.mkdir(parents=True, exist_ok=True)
    path = tmp_path / "guides.ranked.tsv"
    with open(path, "w", newline="") as handle:
        handle.write("\t".join(header) + "\n")
        for i, row in enumerate(rows):
            region, start, end, strand, spacer, passes = row
            values = {"region_id": region, "chrom": chrom, "genome_start": start,
                      "genome_end": end, "orientation": strand,
                      "spacer": spacer, "passes_hard_filters": passes,
                      "rank": i + 1}
            handle.write("\t".join(str(values[c]) for c in header) + "\n")
    return path


class TestSelfHits:
    """A guide whose target region overlaps a promoter window finds ITSELF: a
    perfect full-length match that is not an off-target. Measured at 73% of
    all hits in a genome-wide cCRE run, so counting them would make every
    promoter-overlapping target look maximally risky.
    """

    def test_a_self_hit_is_flagged_and_not_counted(self):
        summary = HitSummary()
        summary.add("g1", "GENEA", 0, 20, self_hit=True)
        count, gene, distance, seed_length, is_promoter = summary.row_for("g1")
        assert count == "0"
        assert is_promoter == "yes"
        assert summary.total_hits == 0

    def test_an_off_target_hit_is_counted(self):
        summary = HitSummary()
        summary.add("g1", "GENEA", 50, 18, self_hit=False)
        assert summary.row_for("g1")[0] == "1"
        assert summary.row_for("g1")[4] == "no"

    def test_a_guide_can_have_both(self):
        """The common real case: the target is a promoter AND there is a
        genuine off-target elsewhere. Both must survive."""
        summary = HitSummary()
        summary.add("g1", "OWN", 0, 20, self_hit=True)
        summary.add("g1", "OTHER", 300, 17, self_hit=False)
        count, gene, distance, seed_length, is_promoter = summary.row_for("g1")
        assert (count, gene, is_promoter) == ("1", "OTHER", "yes")

    def test_the_closest_gene_ignores_self_hits(self):
        """A self-hit sits ~0 bp from a TSS and would otherwise always win the
        'closest' slot, hiding the real nearest off-target gene."""
        summary = HitSummary()
        summary.add("g1", "OWN", 0, 20, self_hit=True)
        summary.add("g1", "REAL", 400, 16, self_hit=False)
        assert summary.row_for("g1")[1] == "REAL"

    def test_counts_are_reported_separately(self):
        summary = HitSummary()
        summary.add("g1", "OWN", 0, 20, self_hit=True)
        summary.add("g2", "OTHER", 10, 16, self_hit=False)
        assert summary.self_hits == 1
        assert summary.total_hits == 1
        assert summary.guides_with_hits == 1


class TestSelfHitDetection:
    def test_forward_and_reverse_are_both_recognised(self):
        """The guide's span covers protospacer plus PAM, and the PAM is at the
        high coordinates for a forward guide and the low ones for a reverse
        guide. Testing an exact start would catch one and miss the other."""
        locus = ("chr1", 9951, 9974)
        assert ss.is_self_hit(locus, "chr1", 9952, 9971)     # forward
        assert ss.is_self_hit(("chr1", 9961, 9984), "chr1", 9965, 9984)  # reverse

    def test_a_different_locus_is_not_a_self_hit(self):
        locus = ("chr1", 9951, 9974)
        assert not ss.is_self_hit(locus, "chr1", 50_000, 50_019)

    def test_a_site_starting_one_base_early_is_not_a_self_hit(self):
        """The two sides use different conventions -- the locus is 0-based
        half-open, the hit is 1-based inclusive -- and comparing them without
        converting admitted one extra base on the left. That position lies
        entirely outside the guide's own span, so a genuine off-target
        beginning one base before the guide's site was reclassified as a
        self-hit and dropped from the count. It needs a shift-tolerant
        (homopolymer or repeat) seed, which is exactly the risky case.
        """
        locus = ("chr1", 9951, 9974)          # 0-based half-open, 23 bp
        assert ss.is_self_hit(locus, "chr1", 9952, 9971)      # FWD protospacer
        assert ss.is_self_hit(locus, "chr1", 9955, 9974)      # RVS protospacer
        # One base before the guide's own span.
        assert not ss.is_self_hit(locus, "chr1", 9951, 9970)
        # One base after it.
        assert not ss.is_self_hit(locus, "chr1", 9956, 9975)

    def test_the_same_offset_on_another_chromosome_is_not_a_self_hit(self):
        locus = ("chr1", 9951, 9974)
        assert not ss.is_self_hit(locus, "chr5", 9952, 9971)

    def test_no_locus_means_nothing_is_a_self_hit(self):
        assert not ss.is_self_hit(None, "chr1", 1, 20)

    def test_the_scanner_marks_a_planted_self_hit(self):
        guide = "GACCTGGAACTCAGTTCAGA"
        left = "A" * 40
        sequence = left + guide + "AGG" + "A" * 40
        # The window starts at genomic 10,000, so the protospacer is at 10,040.
        locus = ("chr1", 10_039, 10_062)
        index = ss.build_seed_index([("g", guide, locus)], 12)
        hits = [h for h in ss.scan_window(sequence, "chr1", 9_999, index, 12)
                if h["protospacer"] == guide]
        assert len(hits) == 1
        assert hits[0]["self_hit"] == "yes"

    def test_the_scanner_does_not_mark_an_elsewhere_hit(self):
        guide = "GACCTGGAACTCAGTTCAGA"
        sequence = "A" * 40 + guide + "AGG" + "A" * 40
        index = ss.build_seed_index([("g", guide, ("chr9", 1, 24))], 12)
        hits = [h for h in ss.scan_window(sequence, "chr1", 9_999, index, 12)
                if h["protospacer"] == guide]
        assert hits and hits[0]["self_hit"] == "no"


class TestIterGuides:
    def test_screens_only_guides_that_pass_the_hard_filters(self):
        """Matches the pre-existing genome-wide workflow: excluded guides are
        not going into a library, and leaving them out cuts both the seed
        index and the work."""

    def test_excluded_guides_are_skipped_by_default(self, tmp_path):
        ranked = ranked_table(tmp_path, [
            ("r1", 1, 20, "FWD", GUIDE_A, "True"),
            ("r1", 2, 21, "FWD", GUIDE_B, "False"),
        ])
        guides = list(iter_guides(ranked))
        assert [spacer for _, spacer, _ in guides] == [GUIDE_A]

    def test_all_guides_overrides_the_filter(self, tmp_path):
        ranked = ranked_table(tmp_path, [
            ("r1", 1, 20, "FWD", GUIDE_A, "True"),
            ("r1", 2, 21, "FWD", GUIDE_B, "False"),
        ])
        guides = list(iter_guides(ranked, passing_only=False))
        assert len(guides) == 2

    def test_a_table_without_the_filter_column_is_screened_in_full(self, tmp_path):
        ranked = ranked_table(tmp_path, [
            ("r1", 1, 20, "FWD", GUIDE_A, ""),
            ("r1", 2, 21, "FWD", GUIDE_B, ""),
        ], filters=False)
        assert len(list(iter_guides(ranked))) == 2

    @pytest.mark.parametrize("value,screened", [
        ("True", True), ("true", True), ("TRUE", True), ("1", True),
        ("False", False), ("false", False), ("", False), ("nan", False),
    ])
    def test_filter_values_are_read_robustly(self, tmp_path, value, screened):
        ranked = ranked_table(tmp_path, [("r1", 1, 20, "FWD", GUIDE_A, value)])
        assert bool(list(iter_guides(ranked))) is screened

    def test_ids_are_unique_and_indexed_by_original_row(self, tmp_path):
        """The index must come from the row's position in the FULL table, not
        from its position among screened guides, or the second pass will
        regenerate different ids."""
        ranked = ranked_table(tmp_path, [
            ("r1", 1, 20, "FWD", GUIDE_A, "False"),
            ("r1", 1, 20, "FWD", GUIDE_A, "True"),
        ])
        ids = [guide_id for guide_id, _, _ in iter_guides(ranked)]
        assert ids == ["r1_1_20_FWD_1"]

    def test_duplicate_guides_collide_across_shards_without_a_prefix(self, tmp_path):
        """The row index is the position in the file being READ, so it
        restarts in every shard. Two duplicated guides landing in different
        shards would then share an id, making the merged hit table unable to
        say which guide a hit belonged to."""
        rows = [("r1", 100, 119, "FWD", GUIDE_A, "True")] * 8
        shard_a = ranked_table(tmp_path / "a", rows[:4])
        shard_b = ranked_table(tmp_path / "b", rows[4:])
        ids = ([g for g, _, _ in iter_guides(shard_a)]
               + [g for g, _, _ in iter_guides(shard_b)])
        assert len(set(ids)) < len(ids), "fixture should demonstrate the clash"

    def test_an_id_prefix_makes_shard_ids_unique(self, tmp_path):
        rows = [("r1", 100, 119, "FWD", GUIDE_A, "True")] * 8
        shard_a = ranked_table(tmp_path / "a", rows[:4])
        shard_b = ranked_table(tmp_path / "b", rows[4:])
        ids = ([g for g, _, _ in iter_guides(shard_a, id_prefix="shard1")]
               + [g for g, _, _ in iter_guides(shard_b, id_prefix="shard2")])
        assert len(set(ids)) == len(ids)
        assert ids[0].startswith("shard1_")

    def test_the_prefix_is_absent_by_default(self, tmp_path):
        ranked = ranked_table(tmp_path, [("r1", 1, 20, "FWD", GUIDE_A, "True")])
        assert list(iter_guides(ranked))[0][0] == "r1_1_20_FWD_0"

    def test_works_without_an_orientation_column(self, tmp_path):
        ranked = ranked_table(tmp_path, [("r1", 1, 20, "", GUIDE_A, "True")],
                              orientation=False)
        assert list(iter_guides(ranked))[0][0] == "r1_1_20_NA_0"

    def test_the_locus_comes_from_the_table(self, tmp_path):
        ranked = ranked_table(tmp_path, [("r1", 100, 123, "FWD", GUIDE_A, "True")],
                              chrom="chr7")
        _gid, _spacer, locus = list(iter_guides(ranked))[0]
        assert locus == ("chr7", 100, 123)

    def test_no_chrom_column_means_no_locus(self, tmp_path):
        """Screening still works; every hit is then a potential off-target."""
        ranked = ranked_table(tmp_path, [("r1", 100, 123, "FWD", GUIDE_A, "True")],
                              chrom=None)
        assert list(iter_guides(ranked))[0][2] is None

    def test_missing_required_column_is_reported(self, tmp_path):
        path = tmp_path / "bad.tsv"
        path.write_text("region_id\tspacer\nr1\tAAA\n")
        with pytest.raises(SystemExit, match="genome_start"):
            list(iter_guides(path))


class TestWriteAnnotated:
    def test_appends_summary_columns_to_every_row(self, tmp_path):
        ranked = ranked_table(tmp_path, [
            ("r1", 1, 20, "FWD", GUIDE_A, "True"),
            ("r2", 2, 21, "RVS", GUIDE_B, "True"),
        ])
        summary = HitSummary()
        summary.add(build_guide_id("r1", 1, 20, "FWD", 0), "GENEA", 42, 9)

        out = tmp_path / "final.tsv"
        assert write_annotated(ranked, out, summary) == 2
        rows = list(csv.DictReader(open(out), delimiter="\t"))
        assert all(column in rows[0] for column in SUMMARY_COLUMNS)
        assert rows[0]["seedfinder_hit_count"] == "1"
        assert rows[0]["seedfinder_closest_gene"] == "GENEA"
        assert rows[1]["seedfinder_hit_count"] == "0"

    def test_excluded_guides_get_blank_not_zero(self, tmp_path):
        ranked = ranked_table(tmp_path, [
            ("r1", 1, 20, "FWD", GUIDE_A, "True"),
            ("r2", 2, 21, "FWD", GUIDE_B, "False"),
        ])
        out = tmp_path / "final.tsv"
        write_annotated(ranked, out, HitSummary())
        rows = list(csv.DictReader(open(out), delimiter="\t"))
        assert rows[0]["seedfinder_hit_count"] == "0"
        assert rows[1]["seedfinder_hit_count"] == ""

    def test_ids_match_between_the_two_streaming_passes(self, tmp_path):
        """The id is regenerated rather than held in memory, so the passes
        must agree. A mismatch would show up as every count being 0."""
        rows = [(f"r{i}", i, i + 19, "FWD", GUIDE_A, "True") for i in range(50)]
        ranked = ranked_table(tmp_path, rows)
        summary = HitSummary()
        for guide_id, _, _ in iter_guides(ranked):
            summary.add(guide_id, "GENE", 1, 20)
        out = tmp_path / "final.tsv"
        write_annotated(ranked, out, summary)
        counts = [r["seedfinder_hit_count"]
                  for r in csv.DictReader(open(out), delimiter="\t")]
        assert counts == ["1"] * 50

    def test_the_prefix_is_used_consistently_by_both_passes(self, tmp_path):
        """If only one pass applied the prefix, every count would read 0."""
        ranked = ranked_table(tmp_path, [
            ("r1", 1, 20, "FWD", GUIDE_A, "True"),
            ("r2", 2, 21, "FWD", GUIDE_B, "True"),
        ])
        summary = HitSummary()
        for guide_id, _, _ in iter_guides(ranked, id_prefix="s7"):
            summary.add(guide_id, "GENEA", 5, 20)
        out = tmp_path / "final.tsv"
        write_annotated(ranked, out, summary, id_prefix="s7")
        counts = [r["seedfinder_hit_count"]
                  for r in csv.DictReader(open(out), delimiter="\t")]
        assert counts == ["1", "1"]

    def test_original_columns_are_preserved(self, tmp_path):
        ranked = ranked_table(tmp_path, [("r1", 1, 20, "FWD", GUIDE_A, "True")])
        out = tmp_path / "final.tsv"
        write_annotated(ranked, out, HitSummary())
        row = next(csv.DictReader(open(out), delimiter="\t"))
        assert row["region_id"] == "r1" and row["rank"] == "1"


def promoter_reference(tmp_path, sequences):
    """A promoter FASTA + TSS table pair, as build-seed-reference produces."""
    promoter = tmp_path / "promoter_windows.fa.gz"
    with gzip.open(promoter, "wt") as handle:
        for chrom, start, sequence in sequences:
            handle.write(f">prom|{chrom}|{start}|{start + len(sequence) - 1}\n")
            handle.write(sequence + "\n")
    tss = tmp_path / "tss_table.tsv"
    tss.write_text("chr\ttss\tgene_name\tstrand\nchr1\t10100\tGENEA\t+\n")
    return promoter, tss


class TestScanToFile:
    def test_writes_hits_and_returns_a_summary(self, tmp_path):
        padding = "A" * 60
        sequence = padding + GUIDE_A + "AGG" + padding
        promoter, tss = promoter_reference(tmp_path, [("chr1", 10_000, sequence)])
        ranked = ranked_table(tmp_path, [("r1", 1, 20, "FWD", GUIDE_A, "True")])

        hits_path = tmp_path / "raw_hits.csv"
        summary = scan_to_file(ranked, promoter, tss, 12, 2000, hits_path)

        rows = list(csv.DictReader(open(hits_path)))
        assert rows, "the planted site should be found"
        assert rows[0]["protospacer"] == GUIDE_A
        assert rows[0]["strand"] == "+"
        assert rows[0]["nearest_tss_gene"] == "GENEA"
        assert summary.total_hits == len(rows)

    def test_hit_table_has_the_declared_columns(self, tmp_path):
        sequence = "A" * 60 + GUIDE_A + "AGG" + "A" * 60
        promoter, tss = promoter_reference(tmp_path, [("chr1", 10_000, sequence)])
        ranked = ranked_table(tmp_path, [("r1", 1, 20, "FWD", GUIDE_A, "True")])
        hits_path = tmp_path / "raw_hits.csv"
        scan_to_file(ranked, promoter, tss, 12, 2000, hits_path)
        with open(hits_path) as handle:
            assert next(csv.reader(handle)) == HIT_COLUMNS

    def test_zero_hits_still_writes_a_header(self, tmp_path):
        """Zero hits is a legitimate result; an empty file is not, because the
        next step cannot tell it from a crash."""
        promoter, tss = promoter_reference(
            tmp_path, [("chr1", 10_000, "A" * 200)])
        ranked = ranked_table(tmp_path, [("r1", 1, 20, "FWD", GUIDE_B, "True")])
        hits_path = tmp_path / "raw_hits.csv"
        summary = scan_to_file(ranked, promoter, tss, 12, 2000, hits_path)
        assert summary.total_hits == 0
        with open(hits_path) as handle:
            assert next(csv.reader(handle)) == HIT_COLUMNS

    def test_nothing_to_screen_is_an_explicit_error(self, tmp_path):
        promoter, tss = promoter_reference(tmp_path, [("chr1", 10_000, "A" * 200)])
        ranked = ranked_table(tmp_path, [("r1", 1, 20, "FWD", GUIDE_A, "False")])
        with pytest.raises(SystemExit, match="--all-guides"):
            scan_to_file(ranked, promoter, tss, 12, 2000, tmp_path / "h.csv")

    def test_a_malformed_spacer_is_reported_clearly(self, tmp_path):
        promoter, tss = promoter_reference(tmp_path, [("chr1", 10_000, "A" * 200)])
        ranked = ranked_table(tmp_path, [("r1", 1, 20, "FWD", "ACGT", "True")])
        with pytest.raises(SystemExit, match="20 bp of ACGT"):
            scan_to_file(ranked, promoter, tss, 12, 2000, tmp_path / "h.csv")

    def test_a_pre_0_5_0_reference_is_refused(self, tmp_path):
        """A three-column TSS table has no strand, so distances could only be
        unsigned. Refusing is better than defaulting: a default would invert
        the sign for every minus-strand gene and still exit 0."""
        promoter, _ = promoter_reference(tmp_path, [("chr1", 10_000, "A" * 200)])
        old_table = tmp_path / "old_tss.tsv"
        old_table.write_text("chr\ttss\tgene_name\nchr1\t10100\tGENEA\n")
        ranked = ranked_table(tmp_path, [("r1", 1, 20, "FWD", GUIDE_A, "True")])
        with pytest.raises(SystemExit, match="build-seed-reference"):
            scan_to_file(ranked, promoter, old_table, 12, 2000, tmp_path / "h.csv")

    def test_a_shared_promoter_reports_both_genes(self, tmp_path):
        """A site between two divergent TSSs must nominate both genes, not
        just the marginally closer one."""
        padding = "A" * 40
        sequence = padding + GUIDE_A + "AGG" + padding
        promoter, _ = promoter_reference(tmp_path, [("chr1", 10_000, sequence)])
        tss = tmp_path / "bidir.tsv"
        tss.write_text("chr\ttss\tgene_name\tstrand\n"
                       "chr1\t10020\tGENE_FWD\t+\n"
                       "chr1\t10120\tGENE_REV\t-\n")
        ranked = ranked_table(tmp_path, [("r1", 1, 20, "FWD", GUIDE_A, "True")])
        hits_path = tmp_path / "raw_hits.csv"
        scan_to_file(ranked, promoter, tss, 12, 2000, hits_path)
        rows = [r for r in csv.DictReader(open(hits_path))
                if r["protospacer"] == GUIDE_A]
        assert rows
        genes = rows[0]["tss_genes_in_window"].split(";")
        assert set(genes) == {"GENE_FWD", "GENE_REV"}
        # Both distances present, and signed per the gene's own strand.
        assert len(rows[0]["tss_distances_in_window"].split(";")) == 2

    def test_an_empty_tss_table_is_reported(self, tmp_path):
        promoter, _ = promoter_reference(tmp_path, [("chr1", 10_000, "A" * 200)])
        empty = tmp_path / "empty.tsv"
        empty.write_text("chr\ttss\tgene_name\n")
        ranked = ranked_table(tmp_path, [("r1", 1, 20, "FWD", GUIDE_A, "True")])
        with pytest.raises(SystemExit, match="no TSS records"):
            scan_to_file(ranked, promoter, empty, 12, 2000, tmp_path / "h.csv")

    def test_both_orientations_are_found_end_to_end(self, tmp_path):
        """The whole-step version of the regression guarded in
        test_seed_scan.py: a plus-strand and a minus-strand site in the same
        reference must both come back."""
        padding = "A" * 40
        plus_window = padding + GUIDE_A + "AGG" + padding
        minus_window = padding + "CCT" + ss.reverse_complement(GUIDE_B) + padding
        promoter, tss = promoter_reference(tmp_path, [
            ("chr1", 10_000, plus_window),
            ("chr1", 10_000, minus_window),
        ])
        ranked = ranked_table(tmp_path, [
            ("r1", 1, 20, "FWD", GUIDE_A, "True"),
            ("r2", 2, 21, "RVS", GUIDE_B, "True"),
        ])
        hits_path = tmp_path / "raw_hits.csv"
        scan_to_file(ranked, promoter, tss, 12, 2000, hits_path)
        rows = list(csv.DictReader(open(hits_path)))
        by_strand = {r["strand"]: r for r in rows
                     if r["protospacer"] in (GUIDE_A, GUIDE_B)}
        assert set(by_strand) == {"+", "-"}
        assert by_strand["+"]["protospacer"] == GUIDE_A
        assert by_strand["-"]["protospacer"] == GUIDE_B
