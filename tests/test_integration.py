"""End-to-end integration test for the step chain.

Runs the real `cretarget` subcommands against a synthetic mini-genome.
FlashFry is not required: its `score` output is a plain TSV with a
documented column set, so it is fabricated here. The full run including
real FlashFry is covered by the Nextflow smoke test in CI.

Skipped automatically when bedtools/samtools are unavailable.
"""

import random
import shutil
import subprocess
import sys

import pandas as pd
import pytest

pytestmark = pytest.mark.skipif(
    shutil.which("bedtools") is None or shutil.which("samtools") is None,
    reason="requires bedtools and samtools on PATH",
)

# Column set of a real FlashFry `score` output (verified against 1.15).
FLASHFRY_COLUMNS = [
    "contig", "start", "stop", "target", "context", "overflow", "orientation",
    "Doench2014OnTarget", "DoenchCFD_maxOT", "DoenchCFD_specificityscore",
    "dangerous_GC", "dangerous_polyT", "dangerous_in_genome", "Hsu2013",
    "Moreno-Mateos2015OnTarget", "JostCRISPRi_maxOT",
    "JostCRISPRi_specificityscore", "basesDiffToClosestHit", "closestHitCount",
    "0-1-2-3-4_mismatch", "FreeEnergy20", "otCount",
]


def run_cli(*args):
    result = subprocess.run(
        [sys.executable, "-m", "cretarget.cli", *map(str, args)],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, f"cretarget {args[0]} failed:\n{result.stderr}"
    return result


@pytest.fixture(scope="module")
def reference(tmp_path_factory):
    """A small two-contig synthetic genome with a .fai index."""
    ref_dir = tmp_path_factory.mktemp("reference")
    fasta = ref_dir / "mini.fa"
    rng = random.Random(7)
    with open(fasta, "w") as fh:
        for chrom, length in [("chr1", 60000), ("chr2", 30000)]:
            seq = "".join(rng.choice("ACGT") for _ in range(length))
            fh.write(f">{chrom}\n")
            for i in range(0, length, 60):
                fh.write(seq[i:i + 60] + "\n")
    subprocess.run(["samtools", "faidx", str(fasta)], check=True)
    return fasta


@pytest.fixture(scope="module")
def prepared(tmp_path_factory, reference):
    """Runs prep-regions over inputs that exercise the awkward cases."""
    work = tmp_path_factory.mktemp("prep")

    (work / "enh.txt").write_text(
        "chr1:10000-10400\n"
        "chr1:10600-11000\n"     # within a blast radius of the one above
        "chr1:20,000-20,300\n"   # thousands separators
        "chr1:10000-10400\n"     # exact duplicate
        "1:5000-5100\n"          # Ensembl-style contig, absent from reference
        "chr1:nonsense\n"        # unparseable
        "chr2:9000-9300\n"
    )
    (work / "snps.txt").write_text("chr1:30000\nchr2:15000\n")

    run_cli("prep-regions",
            "--snp-bed", work / "snps.txt",
            "--enhancer-bed", work / "enh.txt",
            "--fasta", reference,
            "--dedup-blast-radius",
            "--outdir", work / "out")
    return work / "out"


class TestPrepRegions:
    def test_keeps_only_valid_regions(self, prepared):
        meta = pd.read_csv(prepared / "regions.metadata.tsv", sep="\t")
        # 7 enhancer lines - 1 bad contig - 1 unparseable = 5, plus 2 SNPs
        assert len(meta) == 7

    def test_deduplicates_names(self, prepared):
        meta = pd.read_csv(prepared / "regions.metadata.tsv", sep="\t")
        assert not meta["region_id"].duplicated().any()
        assert "chr1:10000-10400_dup1" in set(meta["region_id"])

    def test_parses_comma_separated_coordinates(self, prepared):
        meta = pd.read_csv(prepared / "regions.metadata.tsv", sep="\t")
        row = meta[meta["region_id"] == "chr1:20,000-20,300"].iloc[0]
        assert row["orig_start"] == 20000 and row["orig_end"] == 20300

    def test_flags_linked_neighbours(self, prepared):
        meta = pd.read_csv(prepared / "regions.metadata.tsv", sep="\t")
        row = meta[meta["region_id"] == "chr1:10600-11000"].iloc[0]
        assert "chr1:10000-10400" in str(row["linked_regions"])

    def test_fasta_headers_encode_coordinates(self, prepared):
        headers = [l[1:].strip() for l in open(prepared / "regions.fasta")
                   if l.startswith(">")]
        assert all("__" in h and ":" in h.split("__")[-1] for h in headers)
        assert len(headers) == 7

    def test_snp_window_is_symmetric(self, prepared):
        meta = pd.read_csv(prepared / "regions.metadata.tsv", sep="\t")
        row = meta[meta["region_id"] == "chr1:30000"].iloc[0]
        # "chr1:30000" is 1-based, so the BED point is 29999 and the window
        # is that point +/- snp_window_bp.
        assert row["orig_start"] == 29999 and row["target_ref"] == 29999
        assert row["padded_start"] == 29999 - 150
        assert row["padded_end"] == 29999 + 150


@pytest.fixture(scope="module")
def scored(tmp_path_factory, prepared):
    """A fabricated FlashFry score table over the prepared regions."""
    meta = pd.read_csv(prepared / "regions.metadata.tsv", sep="\t")
    rows = []
    rng = random.Random(11)

    for _, region in meta.iterrows():
        for i in range(4):
            spacer = "".join(rng.choice("ACGT") for _ in range(20))
            if i == 0:
                spacer = "A" * 8 + "GGGG" + "A" * 8   # seed_score 3 -> excluded
            rows.append({
                "contig": region["fasta_header"],
                "start": 10 + i * 25, "stop": 33 + i * 25,
                "target": spacer + "CGG",
                "context": "NNNN" + spacer + "CGGNNN",
                "overflow": "OK" if i != 1 else "OVERFLOW",
                "orientation": "FWD" if i % 2 == 0 else "RVS",
                "Doench2014OnTarget": 0.5, "DoenchCFD_maxOT": 0.2,
                "DoenchCFD_specificityscore": 0.8,
                "dangerous_GC": "NONE", "dangerous_polyT": "NONE",
                "dangerous_in_genome": "IN_GENOME=1",
                "Hsu2013": 70.0, "Moreno-Mateos2015OnTarget": 0.4,
                "JostCRISPRi_maxOT": 0.3,
                "JostCRISPRi_specificityscore": 0.6 + i * 0.05,
                "basesDiffToClosestHit": 3, "closestHitCount": 1,
                "0-1-2-3-4_mismatch": "1,0,0,2,5",
                "FreeEnergy20": -5.0 - i, "otCount": 8,
            })

    path = tmp_path_factory.mktemp("flashfry") / "scored.tsv"
    pd.DataFrame(rows, columns=FLASHFRY_COLUMNS).to_csv(path, sep="\t", index=False)
    return path


@pytest.fixture(scope="module")
def annotated(tmp_path_factory, prepared, scored):
    out = tmp_path_factory.mktemp("annotate") / "annotated.tsv"
    run_cli("annotate",
            "--scored", scored,
            "--region-metadata", prepared / "regions.metadata.tsv",
            "--blast-radius-bp", 1000,
            "--out", out)
    return out


class TestAnnotate:
    def test_recovers_genome_coordinates(self, annotated, prepared):
        df = pd.read_csv(annotated, sep="\t")
        # annotate already carries padded_start through from the metadata join,
        # so use it directly rather than re-merging into a column collision.
        assert (df["genome_start"] == df["padded_start"] + df["start"]).all()
        assert (df["genome_end"] == df["padded_start"] + df["stop"]).all()

    def test_every_guide_matches_a_region(self, annotated):
        df = pd.read_csv(annotated, sep="\t")
        assert df["target_ref"].notna().all()
        assert len(df) == 28

    def test_computes_blast_radius_flag(self, annotated):
        df = pd.read_csv(annotated, sep="\t")
        assert df["within_blast_radius"].all()

    def test_atac_column_is_nan_without_peaks(self, annotated):
        df = pd.read_csv(annotated, sep="\t")
        assert df["distance_to_atac_summit"].isna().all()


@pytest.fixture(scope="module")
def ranked(tmp_path_factory, annotated):
    out = tmp_path_factory.mktemp("rank") / "ranked.tsv"
    run_cli("rank", "--annotated", annotated, "--max-mismatch", 3, "--out", out)
    return out


class TestRank:
    def test_excludes_high_seed_score_guides(self, ranked):
        df = pd.read_csv(ranked, sep="\t")
        seedy = df[df["spacer"] == "A" * 8 + "GGGG" + "A" * 8]
        assert len(seedy) > 0
        assert (seedy["exclude_reason"] == "seed_score").all()

    def test_excludes_overflow_guides(self, ranked):
        df = pd.read_csv(ranked, sep="\t")
        assert (df[df["overflow"] != "OK"]["exclude_reason"] == "overflow").all()

    def test_ranks_restart_within_each_region(self, ranked):
        df = pd.read_csv(ranked, sep="\t")
        passing = df[df["passes_hard_filters"]]
        for _, group in passing.groupby("region_id"):
            ranks = sorted(int(r) for r in group["rank"])
            assert ranks == list(range(1, len(ranks) + 1))

    def test_better_specificity_ranks_higher(self, ranked):
        df = pd.read_csv(ranked, sep="\t")
        passing = df[df["passes_hard_filters"]]
        for _, group in passing.groupby("region_id"):
            ordered = group.sort_values("rank")
            assert ordered["composite_score"].is_monotonic_decreasing

    def test_reports_max_possible_without_launch_dcas9(self, ranked):
        df = pd.read_csv(ranked, sep="\t")
        assert (df["composite_score_max_possible"] == 0.58).all()
        assert (df["composite_score"] <= 0.58 + 1e-9).all()

    def test_normalized_score_is_on_unit_scale(self, ranked):
        df = pd.read_csv(ranked, sep="\t")
        assert df["composite_score_normalized"].between(0, 1).all()

    def test_excluded_guides_are_retained_and_labelled(self, ranked):
        df = pd.read_csv(ranked, sep="\t")
        excluded = df[~df["passes_hard_filters"]]
        assert len(excluded) > 0
        assert excluded["rank"].astype(str).str.startswith("exclude (").all()
