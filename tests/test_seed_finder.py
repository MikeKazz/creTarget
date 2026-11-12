"""Unit tests for the seed-finder-lite integration.

These cover the part creTarget owns: turning seed-finder's hit-level output
(zero to many rows per guide) into one summary row per guide. That logic had
no test coverage, and an end-to-end run only exercises it if seed-finder
actually finds hits -- if it finds none, the empty-result branch runs and the
aggregation below is never touched.
"""

import pandas as pd
import pytest

from cretarget.steps.seed_finder import (
    SUMMARY_COLUMNS, build_guide_ids, summarize_hits,
)


def hits_frame(rows):
    """A hit table shaped like seed-finder-lite's output."""
    return pd.DataFrame(
        rows, columns=["guide_id", "nearest_tss_gene", "tss_distance", "seed_length"]
    )


class TestSummarizeHitsNoHits:
    """Zero hits is a legitimate result, not an error."""

    @pytest.mark.parametrize("empty", [None, "empty_frame"])
    def test_every_guide_is_still_reported(self, empty):
        ids = pd.Series(["g1", "g2", "g3"])
        hits = hits_frame([]) if empty == "empty_frame" else None
        out = summarize_hits(hits, ids)
        assert list(out["guide_id"]) == ["g1", "g2", "g3"]
        assert list(out["seedfinder_hit_count"]) == [0, 0, 0]
        assert out["seedfinder_closest_gene"].isna().all()

    def test_all_summary_columns_present(self):
        out = summarize_hits(None, pd.Series(["g1"]))
        for column in SUMMARY_COLUMNS:
            assert column in out.columns


class TestSummarizeHitsWithHits:
    def test_counts_hits_per_guide(self):
        hits = hits_frame([
            ("g1", "GENEA", 100, 10),
            ("g1", "GENEB", 400, 9),
            ("g2", "GENEC", 50, 12),
        ])
        out = summarize_hits(hits, pd.Series(["g1", "g2", "g3"]))
        counts = dict(zip(out["guide_id"], out["seedfinder_hit_count"]))
        assert counts == {"g1": 2, "g2": 1, "g3": 0}

    def test_picks_the_closest_tss(self):
        hits = hits_frame([
            ("g1", "FAR", 900, 8),
            ("g1", "NEAR", 25, 11),
            ("g1", "MID", 300, 9),
        ])
        row = summarize_hits(hits, pd.Series(["g1"])).iloc[0]
        assert row["seedfinder_closest_gene"] == "NEAR"
        assert row["seedfinder_closest_tss_distance"] == 25
        assert row["seedfinder_closest_seed_length"] == 11

    def test_closest_uses_absolute_distance(self):
        # tss_distance is signed (upstream/downstream). A hit 30 bp upstream is
        # closer than one 500 bp downstream, so the comparison must be on |d|.
        hits = hits_frame([
            ("g1", "DOWNSTREAM", 500, 8),
            ("g1", "UPSTREAM", -30, 10),
        ])
        row = summarize_hits(hits, pd.Series(["g1"])).iloc[0]
        assert row["seedfinder_closest_gene"] == "UPSTREAM"
        assert row["seedfinder_closest_tss_distance"] == -30

    def test_guides_without_hits_keep_count_zero_not_nan(self):
        hits = hits_frame([("g2", "GENEA", 10, 9)])
        out = summarize_hits(hits, pd.Series(["g1", "g2"]))
        g1 = out[out["guide_id"] == "g1"].iloc[0]
        assert g1["seedfinder_hit_count"] == 0
        assert pd.isna(g1["seedfinder_closest_gene"])

    def test_hit_count_is_integer(self):
        # Merging introduces NaN for unhit guides, which makes the column
        # float unless it is cast back.
        hits = hits_frame([("g1", "GENEA", 10, 9)])
        out = summarize_hits(hits, pd.Series(["g1", "g2"]))
        assert out["seedfinder_hit_count"].dtype.kind == "i"

    def test_one_row_per_guide_even_with_many_hits(self):
        hits = hits_frame([("g1", f"GENE{i}", i * 10, 8) for i in range(1, 25)])
        out = summarize_hits(hits, pd.Series(["g1"]))
        assert len(out) == 1

    def test_input_order_is_preserved(self):
        hits = hits_frame([("g3", "GENEA", 10, 9)])
        out = summarize_hits(hits, pd.Series(["g3", "g1", "g2"]))
        assert list(out["guide_id"]) == ["g3", "g1", "g2"]


class TestBuildGuideIds:
    def _frame(self, orientation=True, n=2):
        data = {
            "region_id": [f"reg{i}" for i in range(n)],
            "genome_start": [100 + i for i in range(n)],
            "genome_end": [123 + i for i in range(n)],
        }
        if orientation:
            data["orientation"] = ["FWD"] * n
        return pd.DataFrame(data)

    def test_composes_from_region_and_coordinates(self):
        ids = build_guide_ids(self._frame(n=1))
        assert list(ids) == ["reg0_100_123_FWD"]

    def test_distinguishes_strands_at_the_same_position(self):
        df = pd.DataFrame({
            "region_id": ["r", "r"],
            "genome_start": [100, 100],
            "genome_end": [123, 123],
            "orientation": ["FWD", "RVS"],
        })
        ids = build_guide_ids(df)
        assert ids.iloc[0] != ids.iloc[1]

    def test_suffixes_true_collisions(self):
        # Same region, span and strand: without a suffix these would collide
        # and the summary join would fan out.
        df = pd.DataFrame({
            "region_id": ["r", "r", "r"],
            "genome_start": [100, 100, 100],
            "genome_end": [123, 123, 123],
            "orientation": ["FWD", "FWD", "FWD"],
        })
        ids = build_guide_ids(df)
        assert len(set(ids)) == 3
        assert ids.iloc[0] == "r_100_123_FWD"
        assert all("_dup" in i for i in ids.iloc[1:])

    def test_works_without_an_orientation_column(self):
        ids = build_guide_ids(self._frame(orientation=False, n=1))
        assert list(ids) == ["reg0_100_123_NA"]
