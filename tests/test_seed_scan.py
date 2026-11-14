"""Tests for the promoter-proximal seed scanner.

This is the scientific core of the seed screen, so the central test is not a
hand-written expectation but a comparison against a deliberately naive
brute-force implementation written from the biological definition: a guide
hits a site when its PAM-proximal seed matches the protospacer exactly and the
protospacer is immediately 5' of an NGG on the protospacer's own strand.

The specific regression being guarded: the tool this replaces indexed every
guide under reverse_complement(seed) and used that one key for both scan
branches, so the plus-strand branch matched nothing and only minus-strand
protospacers were ever reported. It exited 0 with plausible output while
missing about half of all true sites. `TestBothOrientations` and
`TestAgainstBruteForce` both fail if that behaviour comes back.
"""

import gzip
import random

import pytest

from cretarget import seed_scan as ss


def rc(seq):
    return ss.reverse_complement(seq)


def random_dna(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def brute_force_hits(sequence, guides, min_seed_length):
    """Every (guide_id, 0-based protospacer start, strand) from first
    principles. Deliberately slow and obvious."""
    sequence = sequence.upper()
    found = set()
    for guide_id, spacer in guides:
        for i in range(len(sequence) - 22):
            # protospacer on the + strand, PAM NGG immediately 3' of it
            protospacer = sequence[i:i + 20]
            pam = sequence[i + 20:i + 23]
            if pam[1:] == "GG" and "N" not in protospacer:
                if protospacer[-min_seed_length:] == spacer[-min_seed_length:]:
                    found.add((guide_id, i, "+"))
            # protospacer on the - strand: its PAM reads NGG on the minus
            # strand, i.e. the three bases 5' of it on the reference are CCN
            if i >= 3:
                protospacer = rc(sequence[i:i + 20])
                pam = rc(sequence[i - 3:i])
                if pam[1:] == "GG" and "N" not in protospacer:
                    if protospacer[-min_seed_length:] == spacer[-min_seed_length:]:
                        found.add((guide_id, i, "-"))
    return found


def scanner_hits(sequence, guides, min_seed_length, offset=0):
    index = ss.build_seed_index(guides, min_seed_length)
    return {
        (h["guide_id"], h["start"] - 1 - offset, h["strand"])
        for h in ss.scan_window(sequence, "chrT", offset, index, min_seed_length)
    }


class TestSequenceHelpers:
    def test_reverse_complement(self):
        assert rc("ACGT") == "ACGT"
        assert rc("AAAC") == "GTTT"
        assert rc("ACGTN") == "NACGT"

    def test_reverse_complement_is_an_involution(self):
        seq = random_dna(200, seed=1)
        assert rc(rc(seq)) == seq

    def test_pam_proximal_match_length_counts_from_the_pam_end(self):
        spacer = "A" * 20
        assert ss.pam_proximal_match_length(spacer, spacer) == 20
        # A mismatch at the PAM-proximal base stops the count immediately.
        assert ss.pam_proximal_match_length("A" * 19 + "C", spacer) == 0
        # A mismatch at the far (5') end leaves 19 matching.
        assert ss.pam_proximal_match_length("C" + "A" * 19, spacer) == 19


class TestSeedIndex:
    def test_indexes_both_orientations(self):
        """The regression: one key per guide made the plus-strand branch
        unmatchable."""
        index = ss.build_seed_index([("g", "ACGT" * 5)], 8)
        seed = ("ACGT" * 5)[-8:]
        assert ("+", seed) in index
        assert ("-", rc(seed)) in index

    def test_guides_sharing_a_seed_are_both_kept(self):
        # At length 8 there are only 65,536 seeds, so collisions are expected
        # at scale and must not drop a guide.
        guides = [("a", "AAAAAAAAAAAACGTACGTA"), ("b", "CCCCCCCCCCCCCGTACGTA")]
        index = ss.build_seed_index(guides, 8)
        assert len(index[("+", "CGTACGTA")]) == 2

    def test_rejects_a_duplicate_guide_id(self):
        with pytest.raises(ss.InvalidGuide, match="duplicate"):
            ss.build_seed_index([("g", "A" * 20), ("g", "C" * 20)], 8)

    @pytest.mark.parametrize("spacer", ["A" * 19, "A" * 21, "A" * 19 + "N",
                                        "A" * 19 + "u"])
    def test_rejects_a_malformed_spacer(self, spacer):
        with pytest.raises(ss.InvalidGuide, match="20 bp of ACGT"):
            ss.build_seed_index([("g", spacer)], 8)

    @pytest.mark.parametrize("n", [0, 21])
    def test_rejects_an_impossible_seed_length(self, n):
        with pytest.raises(ValueError, match="min_seed_length"):
            ss.build_seed_index([("g", "A" * 20)], n)


class TestBothOrientations:
    """A guide must find its own site whichever strand the protospacer is on."""

    GUIDE = "GACCTGGAACTCAGTTCAGA"

    def test_finds_a_plus_strand_protospacer(self):
        left, right = random_dna(60, 2), random_dna(60, 3)
        sequence = left + self.GUIDE + "AGG" + right
        hits = list(ss.scan_window(
            sequence, "chrT", 0,
            ss.build_seed_index([("g", self.GUIDE)], 8), 8))
        exact = [h for h in hits if h["protospacer"] == self.GUIDE]
        assert len(exact) == 1
        hit = exact[0]
        assert hit["strand"] == "+"
        assert hit["seed_length"] == 20
        assert hit["pam"] == "AGG"
        assert hit["start"] == len(left) + 1          # 1-based
        assert hit["end"] == len(left) + 20
        # PAM sits 3' of the protospacer in genomic coordinates
        assert hit["pam_start"] == hit["end"] + 1
        assert hit["pam_end"] == hit["end"] + 3

    def test_finds_a_minus_strand_protospacer(self):
        left, right = random_dna(60, 4), random_dna(60, 5)
        sequence = left + "CCT" + rc(self.GUIDE) + right
        hits = list(ss.scan_window(
            sequence, "chrT", 0,
            ss.build_seed_index([("g", self.GUIDE)], 8), 8))
        exact = [h for h in hits if h["protospacer"] == self.GUIDE]
        assert len(exact) == 1
        hit = exact[0]
        assert hit["strand"] == "-"
        assert hit["seed_length"] == 20
        assert hit["pam"] == "AGG"                     # read on the minus strand
        assert hit["start"] == len(left) + 3 + 1
        # PAM lies 5' of the protospacer in genomic coordinates on this strand
        assert hit["pam_start"] == hit["start"] - 3
        assert hit["pam_end"] == hit["start"] - 1

    def test_the_offset_maps_hits_back_to_the_genome(self):
        left = random_dna(30, 6)
        sequence = left + self.GUIDE + "CGG" + random_dna(30, 7)
        index = ss.build_seed_index([("g", self.GUIDE)], 8)
        hit = next(h for h in ss.scan_window(sequence, "chr9", 1_000_000, index, 8)
                   if h["protospacer"] == self.GUIDE)
        assert hit["chr"] == "chr9"
        assert hit["start"] == 1_000_000 + len(left) + 1


class TestAgainstBruteForce:
    @pytest.mark.parametrize("seed", [11, 12, 13])
    @pytest.mark.parametrize("min_seed_length", [8, 12, 20])
    def test_matches_a_naive_implementation_on_real_sites(self, seed,
                                                          min_seed_length):
        """Guides taken from the sequence itself, so there are sites to find."""
        sequence = random_dna(4000, seed)
        guides = []
        for i in range(0, 3000, 311):
            guides.append((f"p{i}", sequence[i:i + 20]))
            guides.append((f"m{i}", rc(sequence[i:i + 20])))
        expected = brute_force_hits(sequence, guides, min_seed_length)
        assert scanner_hits(sequence, guides, min_seed_length) == expected
        assert expected, "the fixture should contain some real sites"

    def test_matches_a_naive_implementation_on_unrelated_guides(self):
        sequence = random_dna(3000, 21)
        guides = [(f"g{i}", random_dna(20, 500 + i)) for i in range(40)]
        expected = brute_force_hits(sequence, guides, 8)
        assert scanner_hits(sequence, guides, 8) == expected

    def test_a_plus_strand_site_is_actually_in_the_expectation(self):
        """Guards the guard: if brute_force_hits itself only produced minus
        strand hits, the comparison above would pass while both were wrong."""
        sequence = random_dna(4000, 11)
        guides = [(f"p{i}", sequence[i:i + 20]) for i in range(0, 3000, 311)]
        strands = {strand for _, _, strand in brute_force_hits(sequence, guides, 12)}
        assert "+" in strands


class TestAmbiguousBases:
    def test_a_window_of_ns_yields_nothing(self):
        index = ss.build_seed_index([("g", "A" * 20)], 8)
        assert list(ss.scan_window("N" * 100, "chrT", 0, index, 8)) == []

    def test_an_n_inside_the_seed_prevents_a_match(self):
        guide = "GACCTGGAACTCAGTTCAGA"
        broken = guide[:-4] + "N" + guide[-3:]
        sequence = random_dna(40, 8) + broken + "AGG" + random_dna(40, 9)
        index = ss.build_seed_index([("g", guide)], 8)
        assert [h for h in ss.scan_window(sequence, "chrT", 0, index, 8)
                if h["protospacer"] == guide] == []

    def test_a_short_window_yields_nothing(self):
        index = ss.build_seed_index([("g", "A" * 20)], 8)
        assert list(ss.scan_window("ACGT", "chrT", 0, index, 8)) == []


class TestTssIndex:
    def test_nearest_matches_a_linear_scan(self):
        rng = random.Random(99)
        rows = [("chr1", rng.randrange(1, 100_000), f"G{i}") for i in range(500)]
        rows += [("chr2", rng.randrange(1, 100_000), f"H{i}") for i in range(500)]
        index = ss.TssIndex(rows)
        for _ in range(300):
            chrom = rng.choice(["chr1", "chr2"])
            position = rng.randrange(1, 100_000)
            _gene, distance = index.nearest(chrom, position)
            expected = min(abs(position - t) for c, t, _ in rows if c == chrom)
            assert distance == expected

    def test_unknown_chromosome_is_not_an_error(self):
        index = ss.TssIndex([("chr1", 100, "A")])
        assert index.nearest("chrZ", 100) == ("", None)

    def test_reports_the_gene_name(self):
        index = ss.TssIndex([("chr1", 100, "AAA"), ("chr1", 5000, "BBB")])
        assert index.nearest("chr1", 120) == ("AAA", 20)
        assert index.nearest("chr1", 4000) == ("BBB", 1000)

    def test_counts_every_row(self):
        assert len(ss.TssIndex([("chr1", 1, "A"), ("chr2", 2, "B")])) == 2


class TestTssTable:
    def test_header_is_detected_not_assumed(self, tmp_path):
        """Skipping line 1 unconditionally silently drops a gene from a
        headerless file."""
        with_header = tmp_path / "a.tsv"
        with_header.write_text("chr\ttss\tgene_name\nchr1\t100\tAAA\nchr1\t200\tBBB\n")
        without = tmp_path / "b.tsv"
        without.write_text("chr1\t100\tAAA\nchr1\t200\tBBB\n")
        assert len(ss.load_tss_table(with_header)) == 2
        assert len(ss.load_tss_table(without)) == 2

    def test_reads_gzip(self, tmp_path):
        path = tmp_path / "t.tsv.gz"
        with gzip.open(path, "wt") as handle:
            handle.write("chr\ttss\tgene_name\nchr1\t100\tAAA\n")
        assert ss.load_tss_table(path) == [("chr1", 100, "AAA")]

    def test_blank_and_short_rows_are_skipped(self, tmp_path):
        path = tmp_path / "t.tsv"
        path.write_text("chr1\t100\tAAA\n\nchr1\t200\n")
        assert ss.load_tss_table(path) == [("chr1", 100, "AAA")]

    def test_missing_gene_name_becomes_a_placeholder(self, tmp_path):
        path = tmp_path / "t.tsv"
        path.write_text("chr1\t100\t\n")
        assert ss.load_tss_table(path) == [("chr1", 100, "?")]


class TestWindowNames:
    def test_parses_a_promoter_window_name(self):
        assert ss.parse_window_name("prom|chr16|189257|193304") == (
            "chr16", 189257, 193304)

    @pytest.mark.parametrize("name", ["", "chr1|1|2", "prom|chr1|1",
                                       "notprom|chr1|1|2"])
    def test_rejects_anything_else(self, name):
        with pytest.raises(ValueError):
            ss.parse_window_name(name)


class TestFastaStreaming:
    def test_reads_multi_line_records(self, tmp_path):
        path = tmp_path / "x.fa"
        path.write_text(">a\nACGT\nTTTT\n>b\nGG\n")
        assert list(ss.iter_fasta(path)) == [("a", "ACGTTTTT"), ("b", "GG")]

    def test_reads_gzip(self, tmp_path):
        path = tmp_path / "x.fa.gz"
        with gzip.open(path, "wt") as handle:
            handle.write(">a\nACGT\n")
        assert list(ss.iter_fasta(path)) == [("a", "ACGT")]

    def test_empty_file_yields_nothing(self, tmp_path):
        path = tmp_path / "x.fa"
        path.write_text("")
        assert list(ss.iter_fasta(path)) == []


def write_promoter_fasta(path, windows):
    with gzip.open(path, "wt") as handle:
        for chrom, start, sequence in windows:
            handle.write(f">prom|{chrom}|{start}|{start + len(sequence) - 1}\n")
            handle.write(sequence + "\n")


class TestIterHits:
    GUIDE = "GACCTGGAACTCAGTTCAGA"

    def _reference(self, tmp_path, tss_position):
        left = random_dna(100, 31)
        sequence = left + self.GUIDE + "AGG" + random_dna(100, 32)
        path = tmp_path / "prom.fa.gz"
        write_promoter_fasta(path, [("chr1", 10_000, sequence)])
        # The window name carries 1-based coordinates, so sequence[0] IS
        # position 10,000 and the protospacer starts len(left) bases later.
        site_start = 10_000 + len(left)
        index = ss.TssIndex([("chr1", tss_position, "GENEA")])
        return path, index, site_start

    def test_annotates_the_nearest_tss(self, tmp_path):
        path, tss_index, site_start = self._reference(tmp_path, 10_100)
        hits = [h for h in ss.iter_hits(path, [("g", self.GUIDE)], tss_index, 8, 2000)
                if h["protospacer"] == self.GUIDE]
        assert len(hits) == 1
        assert hits[0]["nearest_tss_gene"] == "GENEA"
        assert hits[0]["chr"] == "chr1"
        assert hits[0]["start"] == site_start

    def test_drops_hits_beyond_the_tss_window(self, tmp_path):
        path, tss_index, _ = self._reference(tmp_path, 10_100)
        inside = list(ss.iter_hits(path, [("g", self.GUIDE)], tss_index, 8, 2000))
        outside = list(ss.iter_hits(path, [("g", self.GUIDE)], tss_index, 8, 5))
        assert inside and not outside

    def test_windows_with_unparseable_names_are_skipped(self, tmp_path):
        path = tmp_path / "prom.fa.gz"
        sequence = self.GUIDE + "AGG"
        with gzip.open(path, "wt") as handle:
            handle.write(f">just_a_contig\n{sequence}\n")
        tss_index = ss.TssIndex([("chr1", 1, "GENEA")])
        assert list(ss.iter_hits(path, [("g", self.GUIDE)], tss_index, 8, 2000)) == []

    def test_reports_progress_per_window(self, tmp_path):
        path = tmp_path / "prom.fa.gz"
        write_promoter_fasta(path, [("chr1", 1, random_dna(50, i)) for i in range(4)])
        seen = []
        list(ss.iter_hits(path, [("g", self.GUIDE)], ss.TssIndex([]), 8, 2000,
                          on_window=seen.append))
        assert seen == [0, 1, 2, 3]


class TestUcscUrl:
    def test_builds_a_link(self):
        url = ss.ucsc_url("chr1", 100, 119)
        assert "db=hg38" in url and "chr1:100-119" in url

    def test_adds_a_missing_chr_prefix(self):
        assert "chr7:1-2" in ss.ucsc_url("7", 1, 2)
