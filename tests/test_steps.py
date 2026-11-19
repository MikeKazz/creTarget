"""Unit tests for the pure, science-bearing functions.

These are the pieces where a silent change of behaviour would alter which
guides a user orders and synthesizes, so they are pinned explicitly --
including the two worked examples from the original project notes.
"""

from pathlib import Path

import pandas as pd
import pytest

from cretarget.steps.annotate import parse_contig
from cretarget.steps.prep_regions import (
    compute_window, flag_linked_regions, make_names_unique, validate_contigs,
)
from cretarget.steps.rank import (
    FILTER_ORDER, first_fail_reason, minmax_norm, parse_in_genome_count, seed_score,
)
from cretarget.steps.launch_dcas9 import _join_predictions, unique_valid_spacers


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

    def test_finds_a_region_nested_inside_a_wider_one(self):
        # The regression: comparing only ADJACENT pairs in start order meant
        # `inner` sorted between `big` and `late`, so `big` was never compared
        # with `late` even though `late` sits entirely inside it. The bias was
        # systematic -- the wider the enhancer, the more neighbours it hid --
        # and a falsely empty linked_regions is the damaging direction, since
        # the column warns that one guide may silence two targets.
        regions = [self._region("big", "chr1", 1000, 9000),
                   self._region("inner", "chr1", 1100, 1200),
                   self._region("late", "chr1", 8000, 8100)]
        linked = flag_linked_regions(regions, 1000)
        assert sorted(linked["big"]) == ["inner", "late"]
        assert linked["late"] == ["big"]
        assert linked["inner"] == ["big"]

    def test_is_independent_of_input_order(self):
        import itertools
        regions = [self._region("big", "chr1", 1000, 9000),
                   self._region("inner", "chr1", 1100, 1200),
                   self._region("late", "chr1", 8000, 8100)]
        seen = {
            tuple(sorted((k, tuple(sorted(v)))
                         for k, v in flag_linked_regions(list(p), 1000).items()))
            for p in itertools.permutations(regions)
        }
        assert len(seen) == 1

    def test_links_regions_exactly_at_the_distance_threshold(self):
        # "within link_distance_bp" includes the boundary.
        regions = [self._region("a", "chr1", 0, 100),
                   self._region("b", "chr1", 1100, 1200)]
        assert flag_linked_regions(regions, 1000)["a"] == ["b"]


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

A20 = "ACGT" * 5          # a valid 20 bp spacer
T20 = "TTTT" * 5
N20 = "ACGN" + "ACGT" * 4  # 20 bp but contains N


class TestUniqueValidSpacers:
    def test_dedups_and_preserves_order(self):
        assert unique_valid_spacers([A20, T20, A20]) == [A20, T20]

    def test_drops_non_atgc(self):
        # N-containing spacers cannot be one-hot encoded and would crash the
        # CNN; they must be filtered out and scored as NaN instead.
        assert unique_valid_spacers([A20, N20, "AC-T", ""]) == [A20]

    def test_accepts_lowercase(self):
        # hg38 is soft-masked, so a spacer over a repeat arrives lowercase.
        assert unique_valid_spacers([A20.lower()]) == [A20.lower()]

    def test_rejects_spacers_of_the_wrong_length(self):
        # launch-dCas9's preprocess_seq() takes the length of the FIRST
        # sequence and applies it to the whole batch, so a short one is
        # zero-padded and a long one truncated -- silently, with a plausible
        # score coming back either way. Letting a 19 bp spacer through would
        # corrupt every score in the batch alongside it.
        assert unique_valid_spacers([A20, A20[:19], A20 + "A"]) == [A20]

    def test_all_wrong_length_yields_nothing(self):
        assert unique_valid_spacers(["ACGT", "TTTT"]) == []


# --------------------------------------------------------------------------
# launch-dCas9 prediction join
#
# launch-dCas9 writes the input protospacer next to its score in a `grna`
# column. Matching on that column rather than on row position is what keeps a
# score attached to the guide it was computed for: a row-count check says
# nothing about ORDER, and `--launch-script` points at a user-supplied clone
# of an external project whose revision creTarget does not control.
# --------------------------------------------------------------------------

class TestJoinPredictions:
    @staticmethod
    def _preds(grnas, values):
        return pd.DataFrame({"grna": grnas, "predict": values})

    def test_matches_on_the_key_not_the_row_order(self):
        spacers = [A20, T20]
        # Predictions returned in the OPPOSITE order to the input.
        preds = self._preds([T20, A20], [0.25, 0.75])
        result = _join_predictions(preds, spacers, Path("p.csv"), "enhancerFitness")
        assert result == {A20: 0.75, T20: 0.25}

    def test_positional_zip_would_have_been_wrong_here(self):
        # Pins the actual regression: the old implementation did
        # dict(zip(spacers, preds["predict"])), which on reordered output
        # attaches each score to the wrong spacer while exiting 0.
        spacers = [A20, T20]
        preds = self._preds([T20, A20], [0.25, 0.75])
        positional = dict(zip(spacers, preds["predict"].to_numpy()))
        keyed = _join_predictions(preds, spacers, Path("p.csv"), "enhancerFitness")
        assert positional == {A20: 0.25, T20: 0.75}   # wrong
        assert keyed != positional

    def test_case_insensitive_fallback(self):
        # Soft-masked (lowercase) input, uppercase returned.
        spacers = [A20.lower()]
        preds = self._preds([A20], [0.5])
        result = _join_predictions(preds, spacers, Path("p.csv"), "enhancerFitness")
        assert result == {A20.lower(): 0.5}

    def test_unknown_spacer_is_an_error_not_a_guess(self):
        preds = self._preds([T20], [0.5])
        with pytest.raises(SystemExit):
            _join_predictions(preds, [A20], Path("p.csv"), "enhancerFitness")

    def test_duplicate_keys_are_an_error(self):
        preds = self._preds([A20, A20], [0.1, 0.9])
        with pytest.raises(SystemExit):
            _join_predictions(preds, [A20], Path("p.csv"), "enhancerFitness")

    def test_falls_back_to_order_when_there_is_no_key_column(self, capsys):
        # An older build with no grna column: positional is all that is
        # available, but it must say so rather than pretending to be sure.
        preds = pd.DataFrame({"predict": [0.3, 0.6]})
        result = _join_predictions(preds, [A20, T20], Path("p.csv"), "enhancerFitness")
        assert result == {A20: 0.3, T20: 0.6}
        assert "row order" in capsys.readouterr().err
