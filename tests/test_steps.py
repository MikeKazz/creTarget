"""Unit tests for the pure, science-bearing functions.

These are the pieces where a silent change of behaviour would alter which
guides a user orders and synthesizes, so they are pinned explicitly --
including the two worked examples from the original project notes.
"""

import pandas as pd
import pytest

from cretarget.steps.annotate import parse_contig
from cretarget.steps.prep_regions import (
    compute_window, flag_linked_regions, make_names_unique, validate_contigs,
)
from cretarget.steps.rank import (
    FILTER_ORDER, first_fail_reason, minmax_norm, parse_in_genome_count, seed_score,
)
from cretarget.steps.launch_dcas9 import unique_valid_spacers


# --------------------------------------------------------------------------
# seed_score
# --------------------------------------------------------------------------

class TestSeedScore:
    def test_project_reference_examples(self):
        # Both verified against the original project notes; if either of
        # these changes, the hard filter has changed meaning.
        assert seed_score("ATACTGGATGCCTGTGGCTG") == 1
        assert seed_score("AAACCCAAGTACTATCTTTG") == 0

    def test_counts_only_pam_proximal_12bp(self):
        # GG in the distal 8bp must be ignored entirely.
        assert seed_score("GGGGGGGG" + "ATATATATATAT") == 0

    def test_overlapping_runs(self):
        assert seed_score("A" * 8 + "GGG" + "A" * 9) == 2   # GGG -> 2
        assert seed_score("A" * 8 + "GGGG" + "A" * 8) == 3  # GGGG -> 3

    def test_is_case_insensitive(self):
        assert seed_score("a" * 8 + "ggg" + "a" * 9) == 2

    @pytest.mark.parametrize("value", [None, "", "ACGT", 42, float("nan")])
    def test_returns_none_for_unusable_input(self, value):
        assert seed_score(value) is None


# --------------------------------------------------------------------------
# parse_in_genome_count
# --------------------------------------------------------------------------

class TestParseInGenomeCount:
    @pytest.mark.parametrize("value,expected", [
        ("NONE", 0),
        ("none", 0),
        ("IN_GENOME=1", 1),
        ("IN_GENOME=12", 12),
    ])
    def test_parses_known_forms(self, value, expected):
        assert parse_in_genome_count(value) == expected

    @pytest.mark.parametrize("value", [None, 3, "IN_GENOME=x", "garbage"])
    def test_returns_none_when_unparseable(self, value):
        assert parse_in_genome_count(value) is None


# --------------------------------------------------------------------------
# first_fail_reason
# --------------------------------------------------------------------------

class TestFirstFailReason:
    def _frame(self, **flags):
        return pd.DataFrame({col: [flags.get(col, False)] for col, _ in FILTER_ORDER})

    def test_none_when_all_pass(self):
        assert first_fail_reason(self._frame()).iloc[0] is None

    def test_reports_highest_priority_failure(self):
        both = self._frame(fail_seed_score=True, fail_overflow=True)
        assert first_fail_reason(both).iloc[0] == "seed_score"

    @pytest.mark.parametrize("col,reason", FILTER_ORDER)
    def test_each_filter_reports_itself(self, col, reason):
        assert first_fail_reason(self._frame(**{col: True})).iloc[0] == reason

    def test_matches_rowwise_reference_implementation(self):
        """Guards the vectorized rewrite against the original semantics."""
        import numpy as np
        rng = np.random.default_rng(0)
        df = pd.DataFrame({col: rng.random(2000) < 0.3 for col, _ in FILTER_ORDER})

        def reference(row):
            for col, reason in FILTER_ORDER:
                if row[col]:
                    return reason
            return None

        expected = df.apply(reference, axis=1)
        assert (first_fail_reason(df).fillna("~") == expected.fillna("~")).all()


# --------------------------------------------------------------------------
# minmax_norm
# --------------------------------------------------------------------------

class TestMinMaxNorm:
    def test_scales_to_unit_interval(self):
        result = minmax_norm(pd.Series([-10.0, -5.0, 0.0]))
        assert list(result) == [0.0, 0.5, 1.0]

    def test_neutral_when_no_variance(self):
        # A region with one guide, or with identical folding energies, must
        # not become 0 or 1 -- that would silently reward or punish it.
        assert list(minmax_norm(pd.Series([3.0, 3.0]))) == [0.5, 0.5]

    def test_neutral_when_all_nan(self):
        assert list(minmax_norm(pd.Series([float("nan")] * 3))) == [0.5, 0.5, 0.5]


# --------------------------------------------------------------------------
# parse_contig
# --------------------------------------------------------------------------

class TestParseContig:
    def test_round_trips_the_header_format(self):
        assert parse_contig("enh1__chr1:1000-2000") == ("enh1", "chr1", 1000, 2000)

    def test_region_id_may_contain_underscores_and_colons(self):
        # Region IDs are often the user's own coordinate string, so they
        # routinely contain both. rsplit is what makes this work.
        assert parse_contig("chr1:50,000-50,300__chr1:49950-50350") == (
            "chr1:50,000-50,300", "chr1", 49950, 50350)
        assert parse_contig("my_region_2__chr7:10-20") == ("my_region_2", "chr7", 10, 20)

    @pytest.mark.parametrize("value", ["nonsense", "", None, "a__b", "x__chr1:a-b"])
    def test_returns_nones_when_unparseable(self, value):
        assert parse_contig(value) == (None, None, None, None)


# --------------------------------------------------------------------------
# prep_regions helpers
# --------------------------------------------------------------------------

class TestComputeWindow:
    def test_snp_window_is_symmetric_around_the_point(self):
        region = {"chrom": "chr1", "start": 1000, "end": 1001, "input_type": "snp"}
        start, end, ref = compute_window(region, 150, 50, {"chr1": 100000})
        assert (start, end, ref) == (850, 1150, 1000)

    def test_enhancer_window_pads_and_uses_midpoint(self):
        region = {"chrom": "chr1", "start": 1000, "end": 2000, "input_type": "enhancer"}
        start, end, ref = compute_window(region, 150, 50, {"chr1": 100000})
        assert (start, end, ref) == (950, 2050, 1500)

    def test_clamps_to_contig_bounds(self):
        region = {"chrom": "chr1", "start": 10, "end": 20, "input_type": "enhancer"}
        start, end, _ = compute_window(region, 150, 50, {"chr1": 60})
        assert start == 0 and end == 60

    def test_no_clamping_when_contig_length_unknown(self):
        region = {"chrom": "chrUn", "start": 10, "end": 20, "input_type": "enhancer"}
        _, end, _ = compute_window(region, 150, 50, {})
        assert end == 70


class TestMakeNamesUnique:
    def test_suffixes_repeats_only(self):
        regions = [{"name": n} for n in ["a", "b", "a", "a"]]
        assert [r["name"] for r in make_names_unique(regions)] == [
            "a", "b", "a_dup1", "a_dup2"]

    def test_leaves_unique_names_untouched(self):
        regions = [{"name": n} for n in ["a", "b", "c"]]
        assert [r["name"] for r in make_names_unique(regions)] == ["a", "b", "c"]


class TestFlagLinkedRegions:
    def _region(self, name, chrom, start, end):
        return {"name": name, "chrom": chrom, "start": start, "end": end}

    def test_flags_near_neighbours_symmetrically(self):
        regions = [self._region("a", "chr1", 1000, 1100),
                   self._region("b", "chr1", 1500, 1600)]
        linked = flag_linked_regions(regions, 1000)
        assert linked["a"] == ["b"] and linked["b"] == ["a"]

    def test_ignores_distant_regions(self):
        regions = [self._region("a", "chr1", 1000, 1100),
                   self._region("b", "chr1", 50000, 50100)]
        linked = flag_linked_regions(regions, 1000)
        assert linked["a"] == [] and linked["b"] == []

    def test_never_links_across_contigs(self):
        regions = [self._region("a", "chr1", 1000, 1100),
                   self._region("b", "chr2", 1000, 1100)]
        linked = flag_linked_regions(regions, 1000)
        assert linked["a"] == [] and linked["b"] == []


class TestValidateContigs:
    def test_drops_unknown_contigs_and_keeps_the_rest(self):
        regions = [{"chrom": "chr1", "name": "a"}, {"chrom": "99", "name": "b"}]
        kept = validate_contigs(regions, {"chr1": 1000})
        assert [r["name"] for r in kept] == ["a"]

    def test_fails_loudly_when_nothing_matches(self):
        # The Ensembl-vs-UCSC naming mismatch. Failing here with a hint beats
        # failing inside bedtools with a cryptic message.
        regions = [{"chrom": "1", "name": "a"}, {"chrom": "2", "name": "b"}]
        with pytest.raises(SystemExit):
            validate_contigs(regions, {"chr1": 1000, "chr2": 1000})


# --------------------------------------------------------------------------
# launch-dCas9 input filtering
# --------------------------------------------------------------------------

class TestUniqueValidSpacers:
    def test_dedups_and_preserves_order(self):
        assert unique_valid_spacers(["ACGT", "TTTT", "ACGT"]) == ["ACGT", "TTTT"]

    def test_drops_non_atgc(self):
        # N-containing spacers cannot be one-hot encoded and would crash the
        # CNN; they must be filtered out and scored as NaN instead.
        assert unique_valid_spacers(["ACGT", "ACGN", "AC-T", ""]) == ["ACGT"]

    def test_accepts_lowercase(self):
        assert unique_valid_spacers(["acgt"]) == ["acgt"]
