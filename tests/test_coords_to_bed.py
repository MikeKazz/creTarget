"""Unit tests for the coords-to-bed utility.

This subcommand had no test coverage at all. It exists to turn a column of
"chr:start-end" strings in an existing results table back into pipeline input,
so a silent mis-parse here feeds wrong coordinates into a whole design run.
"""

import argparse

import pytest

from cretarget.steps.coords_to_bed import run


def make_args(input_path, out_path, column="target_enhancer", delimiter="\t"):
    return argparse.Namespace(
        input=str(input_path), column=column, out=str(out_path), delimiter=delimiter
    )


def write_table(path, rows, header="target_enhancer", sep="\t"):
    path.write_text(header + "\n" + "\n".join(rows) + "\n")
    return path


def read_bed(path):
    return [l.rstrip("\n").split("\t") for l in open(path) if l.strip()]


class TestCoordsToBed:
    def test_converts_coordinates_to_bed(self, tmp_path):
        src = write_table(tmp_path / "in.tsv", ["chr1:1000-2000"])
        out = tmp_path / "out.bed"
        run(make_args(src, out))
        assert read_bed(out) == [["chr1", "1000", "2000", "chr1:1000-2000"]]

    def test_strips_thousands_separators(self, tmp_path):
        # Coordinates pasted from a genome browser carry commas.
        src = write_table(tmp_path / "in.tsv", ["chr1:1,000,000-1,000,500"])
        out = tmp_path / "out.bed"
        run(make_args(src, out))
        row = read_bed(out)[0]
        assert row[1] == "1000000" and row[2] == "1000500"

    def test_deduplicates_repeated_values(self, tmp_path):
        src = write_table(tmp_path / "in.tsv",
                          ["chr1:100-200", "chr1:100-200", "chr2:300-400"])
        out = tmp_path / "out.bed"
        run(make_args(src, out))
        assert len(read_bed(out)) == 2

    @pytest.mark.parametrize("null", [".", "NA", "NaN", "None", ""])
    def test_skips_null_placeholders(self, tmp_path, null):
        src = write_table(tmp_path / "in.tsv", [null, "chr1:100-200"])
        out = tmp_path / "out.bed"
        run(make_args(src, out))
        assert read_bed(out) == [["chr1", "100", "200", "chr1:100-200"]]

    def test_skips_unparseable_values_but_keeps_the_rest(self, tmp_path):
        src = write_table(tmp_path / "in.tsv",
                          ["not-a-coordinate", "chr1:100-200", "chr3:xyz"])
        out = tmp_path / "out.bed"
        run(make_args(src, out))
        assert read_bed(out) == [["chr1", "100", "200", "chr1:100-200"]]

    def test_preserves_the_original_string_as_the_name(self, tmp_path):
        # The BED name column must round-trip exactly, so results can be
        # joined back to the source table.
        src = write_table(tmp_path / "in.tsv", ["chr1:1,000-2,000"])
        out = tmp_path / "out.bed"
        run(make_args(src, out))
        assert read_bed(out)[0][3] == "chr1:1,000-2,000"

    def test_fails_loudly_on_a_missing_column(self, tmp_path):
        src = write_table(tmp_path / "in.tsv", ["chr1:100-200"])
        with pytest.raises(SystemExit) as excinfo:
            run(make_args(src, tmp_path / "out.bed", column="no_such_column"))
        assert "no_such_column" in str(excinfo.value)

    def test_supports_comma_delimited_input(self, tmp_path):
        src = tmp_path / "in.csv"
        src.write_text("id,target_enhancer\n1,chr1:100-200\n")
        out = tmp_path / "out.bed"
        run(make_args(src, out, delimiter=","))
        assert read_bed(out) == [["chr1", "100", "200", "chr1:100-200"]]

    def test_handles_contigs_with_underscores_and_dots(self, tmp_path):
        src = write_table(tmp_path / "in.tsv", ["chr1_KI270706v1_random:10-20"])
        out = tmp_path / "out.bed"
        run(make_args(src, out))
        assert read_bed(out)[0][0] == "chr1_KI270706v1_random"
