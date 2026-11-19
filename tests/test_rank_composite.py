"""Unit tests for rank.run()'s composite assembly.

rank.py decides which guides a user synthesizes, and until now its scoring
block had no unit test at all -- only an end-to-end integration test whose
fixture manufactures a distinct score per row, so the missing-component and
saturation regimes that real FlashFry output actually produces were never
exercised.

The fixture here deliberately mirrors real output: a `basesDiffToClosestHit`
of "UNK" (which is what FlashFry emits when it finds no off-target) and
specificity scores at their maximum.
"""

import argparse

import pandas as pd
import pytest

from cretarget.steps import rank as rank_step
from cretarget.steps.rank import MAX_WITH_LAUNCH, MAX_WITHOUT_LAUNCH, WEIGHTS

SPACER = "ACGT" * 5


def _row(**overrides):
    """One clean, passing guide row."""
    row = {
        "region_id": "r1",
        "contig": "r1__chr1:1000-2000",
        "target": SPACER + "AGG",
        "genome_start": 1000,
        "genome_end": 1023,
        "orientation": "FWD",
        # Hard filters, all passing.
        "dangerous_GC": "NONE",
        "dangerous_polyT": "NONE",
        "dangerous_in_genome": "IN_GENOME=1",
        "overflow": "OK",
        # Score components, all at maximum.
        "JostCRISPRi_specificityscore": 1.0,
        "Hsu2013": 100.0,
        "DoenchCFD_specificityscore": 1.0,
        "basesDiffToClosestHit": "UNK",
        "FreeEnergy20": 0.0,
    }
    row.update(overrides)
    return row


def _run(tmp_path, rows, max_mismatch=3):
    tmp_path.mkdir(parents=True, exist_ok=True)
    source = tmp_path / "annotated.tsv"
    pd.DataFrame(rows).to_csv(source, sep="\t", index=False)
    out = tmp_path / "ranked.tsv"
    rank_step.run(argparse.Namespace(
        annotated=str(source), max_mismatch=max_mismatch, out=str(out)))
    return pd.read_csv(out, sep="\t")


class TestClosestHitSentinel:
    def test_unk_is_treated_as_no_offtarget_found(self, tmp_path):
        """FlashFry's ClosestHit.scala emits "UNK" exactly when `closest ==
        Int.MaxValue`, i.e. no off-target with >=1 mismatch existed, and
        writes closestOTCount 0 alongside. Maximally safe is the correct
        reading, and the component must still count toward the denominator.

        Asserted as a CONTRAST against a guide whose closest off-target is at
        0 mismatches, so the test does not depend on the folding term.
        """
        unk = _run(tmp_path / "a", [_row(basesDiffToClosestHit="UNK")]).iloc[0]
        worst = _run(tmp_path / "b", [_row(basesDiffToClosestHit=0)]).iloc[0]

        # Neither is "missing" -- both were measured.
        assert unk["composite_components_missing"] == 0
        assert worst["composite_components_missing"] == 0
        assert unk["composite_score_max_possible"] == pytest.approx(
            MAX_WITHOUT_LAUNCH)
        assert worst["composite_score_max_possible"] == pytest.approx(
            MAX_WITHOUT_LAUNCH)
        # UNK earns the full closest-distance weight; 0 mismatches earns none.
        assert unk["composite_score"] - worst["composite_score"] == pytest.approx(
            WEIGHTS["closest_distance"])

    def test_an_unreadable_value_is_missing_not_safe(self, tmp_path):
        """Anything that is NOT the documented sentinel is a genuine unknown
        and must not collect full credit for safety."""
        unreadable = _run(tmp_path / "a", [_row(basesDiffToClosestHit="???")]).iloc[0]
        unk = _run(tmp_path / "b", [_row(basesDiffToClosestHit="UNK")]).iloc[0]

        assert unreadable["composite_components_missing"] == 1
        # The closest-distance weight is dropped from the scale entirely...
        assert unreadable["composite_score_max_possible"] == pytest.approx(
            MAX_WITHOUT_LAUNCH - WEIGHTS["closest_distance"])
        # ...and from the score, so it is not silently credited as safe.
        assert unreadable["composite_score"] == pytest.approx(
            unk["composite_score"] - WEIGHTS["closest_distance"])

    def test_a_real_distance_is_scaled_by_max_mismatch(self, tmp_path):
        far = _run(tmp_path / "a", [_row(basesDiffToClosestHit=3)],
                   max_mismatch=3).iloc[0]
        near = _run(tmp_path / "b", [_row(basesDiffToClosestHit=0)],
                    max_mismatch=3).iloc[0]
        assert far["composite_score"] - near["composite_score"] == pytest.approx(
            WEIGHTS["closest_distance"])
        # The scale is unchanged either way: both were measured.
        assert far["composite_score_max_possible"] == pytest.approx(
            near["composite_score_max_possible"])


class TestMissingComponents:
    def test_a_missing_component_is_excluded_not_scored_zero(self, tmp_path):
        """The regression: .fillna(0) scored "never measured" identically to
        "measured, and the worst possible". launch-dCas9 returns NaN for any
        spacer it cannot encode, so those guides silently lost 30% of the
        scale and were published with a confident mid-range score.
        """
        rows = [
            _row(launch_dCas9_enhancer=0.4),
            _row(launch_dCas9_enhancer=None),
        ]
        df = _run(tmp_path, rows)
        scored = df[df["launch_dCas9_enhancer"].notna()].iloc[0]
        unscored = df[df["launch_dCas9_enhancer"].isna()].iloc[0]

        assert scored["composite_components_missing"] == 0
        assert unscored["composite_components_missing"] == 1

        # Its denominator excludes the component it never had...
        assert unscored["composite_score_max_possible"] == pytest.approx(
            MAX_WITH_LAUNCH - WEIGHTS["launch_enhancer"])
        assert scored["composite_score_max_possible"] == pytest.approx(
            MAX_WITH_LAUNCH)
        # ...so on the comparable scale it is NOT dragged below the guide that
        # was actually measured at a mediocre 0.4. Under the old .fillna(0)
        # the unscored guide scored strictly lower, which read as "worse"
        # rather than "unknown".
        assert (unscored["composite_score_normalized"]
                > scored["composite_score_normalized"])

    def test_normalized_never_exceeds_one_when_saturated(self, tmp_path):
        # Every component at maximum: the float sum can land a half-ulp above
        # max_possible, which printed as 1.0000000000000002 in a real run.
        df = _run(tmp_path, [_row(launch_dCas9_enhancer=1.0)])
        assert df.loc[0, "composite_score_normalized"] <= 1.0

    def test_launch_score_outside_unit_range_is_clipped(self, tmp_path):
        # launch-dCas9's fitness head has no documented range and is the only
        # unbounded input to the sum.
        df = _run(tmp_path, [_row(launch_dCas9_enhancer=5.0)])
        assert df.loc[0, "composite_score"] <= MAX_WITH_LAUNCH + 1e-9
        assert df.loc[0, "composite_score_normalized"] <= 1.0


class TestExcludedGuides:
    def test_excluded_guides_have_no_score(self, tmp_path):
        rows = [_row(), _row(dangerous_polyT="POLY_T")]
        df = _run(tmp_path, rows)
        excluded = df[~df["passes_hard_filters"]]
        assert len(excluded) == 1
        assert str(excluded.iloc[0]["rank"]).startswith("exclude (")
        for column in ("composite_score", "composite_score_normalized",
                       "composite_score_max_possible"):
            assert pd.isna(excluded.iloc[0][column]), column

    def test_passing_guides_are_ranked_from_one_per_region(self, tmp_path):
        rows = [
            _row(region_id="r1", JostCRISPRi_specificityscore=0.9),
            _row(region_id="r1", JostCRISPRi_specificityscore=0.5),
            _row(region_id="r2", JostCRISPRi_specificityscore=0.7),
        ]
        df = _run(tmp_path, rows)
        passing = df[df["passes_hard_filters"]]
        for _, group in passing.groupby("region_id"):
            assert sorted(group["rank"].astype(int)) == list(
                range(1, len(group) + 1))
