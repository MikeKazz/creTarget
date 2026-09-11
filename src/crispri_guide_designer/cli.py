"""
cli.py

Entry point for the `cretarget` command.

Subcommands:
  run              design + score + rank guides (the main pipeline)
  configure        print the config file location / create it if missing
  setup-reference  download hg38 + FlashFry index from Zenodo

Design note: the individual pipeline steps (region prep, annotate, rank,
launch-dCas9 scoring) are NOT reimplemented here -- they're the same
proven scripts from bin/, bundled into the package under
crispri_guide_designer/scripts/, invoked as subprocesses. This keeps the
already-validated logic exactly as-is; this module's job is orchestration,
config, and chunked/parallel execution (the new capability), not
re-deriving guide-design logic.
"""

import argparse
import glob
import os
import shutil
import subprocess
import sys
from importlib.resources import files
from pathlib import Path

import pandas as pd

from . import config as cfg
from . import atac_cache
from . import seed_finder_integrate
from .executor import get_executor
from .flashfry_chunked import run_chunked_discover_score


def _script_path(name: str) -> str:
    return str(files("crispri_guide_designer.scripts").joinpath(name))


def _run_py_script(name: str, args: list[str]):
    cmd = [sys.executable, _script_path(name)] + args
    result = subprocess.run(cmd)
    if result.returncode != 0:
        raise SystemExit(f"ERROR: {name} failed (exit {result.returncode})")


def _default_n_workers(n_regions: int) -> int:
    """min(cpu_count, available_ram_gb // 3), capped at n_regions -- validated
    on real hardware (see project notes): this formula kept a 3517-region /
    ~464K-guide run comfortably within available memory while matching the
    speed of a much more aggressive, memory-risky chunk count."""
    cpu_count = os.cpu_count() or 4
    try:
        import psutil
        available_gb = psutil.virtual_memory().available / (1024 ** 3)
    except ImportError:
        # Fallback if psutil isn't installed: read /proc/meminfo directly (Linux)
        available_gb = 16.0
        try:
            with open("/proc/meminfo") as f:
                for line in f:
                    if line.startswith("MemAvailable:"):
                        available_gb = int(line.split()[1]) / (1024 ** 2)
                        break
        except FileNotFoundError:
            pass
    mem_limited = max(1, int(available_gb // 3))
    return max(1, min(cpu_count, mem_limited, n_regions))


def cmd_configure(args):
    path = cfg.ensure_config_exists()
    print(f"Config file: {path}")
    print("Edit it directly, or set individual values with:")
    print("  cretarget configure --set reference.fasta=/path/to/hg38.fa")
    if args.set:
        # Minimal manual TOML patch (avoid adding a write-capable TOML dep):
        # read as text, replace matching `key = ` lines under the right section.
        text = path.read_text()
        for kv in args.set:
            section_key, value = kv.split("=", 1)
            section, key = section_key.split(".")
            lines = text.splitlines()
            in_section = False
            found = False
            for i, line in enumerate(lines):
                if line.strip() == f"[{section}]":
                    in_section = True
                    continue
                if line.strip().startswith("[") and in_section:
                    in_section = False
                if in_section and line.strip().startswith(f"{key} "):
                    lines[i] = f'{key} = "{value}"'
                    found = True
                    break
            if not found:
                print(f"WARNING: could not find [{section}].{key} in config to update", file=sys.stderr)
            text = "\n".join(lines) + "\n"
        path.write_text(text)
        print(f"Updated {path}")


def cmd_setup_reference(args):
    # Thin Python wrapper around the same logic as bin/download_reference.sh
    import hashlib
    import urllib.request

    ZENODO_RECORD_ID = "21115562"
    BASE_URL = f"https://zenodo.org/records/{ZENODO_RECORD_ID}/files"
    outdir = Path(args.outdir).expanduser()
    outdir.mkdir(parents=True, exist_ok=True)

    files_to_get = [
        ("hg38.fa.gz.sha256", None),
        ("hg38.fa.gz", "hg38.fa.gz.sha256"),
        ("hg38_standard_cas9ngg_db.tar.gz.sha256", None),
        ("hg38_standard_cas9ngg_db.tar.gz", "hg38_standard_cas9ngg_db.tar.gz.sha256"),
    ]
    for fname, sha_name in files_to_get:
        dest = outdir / fname
        if dest.exists() and not args.force:
            print(f"  -> {fname} already present, skipping", file=sys.stderr)
            continue
        print(f"  -> downloading {fname}", file=sys.stderr)
        urllib.request.urlretrieve(f"{BASE_URL}/{fname}?download=1", str(dest))
        if sha_name:
            expected = (outdir / sha_name).read_text().split()[0]
            actual = hashlib.sha256(dest.read_bytes()).hexdigest()
            if expected != actual:
                raise SystemExit(f"ERROR: checksum mismatch for {fname}")

    fasta_gz = outdir / "hg38.fa.gz"
    fasta = outdir / "hg38.fa"
    if not fasta.exists() or args.force:
        print("  -> decompressing FASTA (bedtools compatibility)", file=sys.stderr)
        subprocess.run(["bgzip", "-d", "-k", "-f", str(fasta_gz)], check=True)
        subprocess.run(["samtools", "faidx", str(fasta)], check=True)

    db_tar = outdir / "hg38_standard_cas9ngg_db.tar.gz"
    db_file = outdir / "hg38_standard_cas9ngg_db"
    if not db_file.exists() or args.force:
        print("  -> extracting FlashFry index", file=sys.stderr)
        subprocess.run(["tar", "xzf", str(db_tar), "-C", str(outdir)], check=True)

    print(f"\nDone. Set in config ({cfg.CONFIG_PATH}):")
    print(f"  [reference] fasta = \"{fasta}\"")
    print(f"  [reference] flashfry_db = \"{db_file}\"")


def cmd_run(args):
    conf = cfg.load_config()
    fasta = cfg.require(conf, "reference", "fasta") if not args.fasta else args.fasta
    flashfry_jar = cfg.require(conf, "reference", "flashfry_jar") if not args.flashfry_jar else args.flashfry_jar
    flashfry_db = cfg.require(conf, "reference", "flashfry_db") if not args.flashfry_db else args.flashfry_db

    if not args.snp_bed and not args.enhancer_bed:
        raise SystemExit("ERROR: at least one of --snp-bed / --enhancer-bed is required")

    if args.atac_narrowpeak and args.atac_cell_type:
        raise SystemExit("ERROR: --atac-narrowpeak and --atac-cell-type are mutually exclusive, pick one")

    if args.run_seed_finder and (args.seed_finder_min_length is None or args.seed_finder_tss_window is None):
        raise SystemExit("ERROR: --run-seed-finder requires both --seed-finder-min-length and --seed-finder-tss-window")

    atac_narrowpeak_path = args.atac_narrowpeak
    if args.atac_cell_type:
        atac_narrowpeak_path = atac_cache.get_atac_path(args.atac_cell_type)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # ---- 1. Region prep (reuses bin/prep_regions.py, unmodified) ----
    fasta_parts, meta_parts = [], []
    if args.snp_bed:
        d = outdir / "prep_snp"
        d.mkdir(exist_ok=True)
        prep_args = ["--bed", args.snp_bed, "--input_type", "snp", "--fasta", fasta,
                     "--snp_window_bp", str(args.snp_window_bp),
                     "--link_distance_bp", str(args.blast_radius_bp), "--outdir", str(d)]
        if args.dedup_blast_radius:
            prep_args.append("--dedup_blast_radius")
        _run_py_script("prep_regions.py", prep_args)
        fasta_parts.append(d / "regions.fasta")
        meta_parts.append(d / "regions.metadata.tsv")

    if args.enhancer_bed:
        d = outdir / "prep_enhancer"
        d.mkdir(exist_ok=True)
        prep_args = ["--bed", args.enhancer_bed, "--input_type", "enhancer", "--fasta", fasta,
                     "--flank_bp", str(args.flank_bp),
                     "--link_distance_bp", str(args.blast_radius_bp), "--outdir", str(d)]
        if args.dedup_blast_radius:
            prep_args.append("--dedup_blast_radius")
        _run_py_script("prep_regions.py", prep_args)
        fasta_parts.append(d / "regions.fasta")
        meta_parts.append(d / "regions.metadata.tsv")

    combined_fasta = outdir / "regions.combined.fasta"
    with open(combined_fasta, "w") as out:
        for p in fasta_parts:
            out.write(Path(p).read_text())

    combined_meta = outdir / "regions.combined.metadata.tsv"
    with open(combined_meta, "w") as out:
        out.write(Path(meta_parts[0]).read_text().splitlines()[0] + "\n")  # header once
        for p in meta_parts:
            lines = Path(p).read_text().splitlines()[1:]
            out.write("\n".join(lines) + ("\n" if lines else ""))

    n_regions = sum(1 for line in open(combined_fasta) if line.startswith(">"))
    print(f"=== {n_regions} total regions ===", file=sys.stderr)

    # ---- 2. Chunked FlashFry discover+score ----
    n_workers = args.n_workers or _default_n_workers(n_regions)
    print(f"=== FlashFry (executor={args.executor}, n_workers={n_workers}) ===", file=sys.stderr)
    if args.executor == "local":
        executor = get_executor("local", n_workers=n_workers)
    else:
        slurm_conf = conf.get("slurm", {})
        executor = get_executor(
            "slurm",
            partition=slurm_conf.get("partition", ""),
            gres=slurm_conf.get("gres", ""),
            mem=slurm_conf.get("mem", "16G"),
            cpus=slurm_conf.get("cpus", "8"),
            time=slurm_conf.get("time", "2:00:00"),
        )

    scored_chunks = run_chunked_discover_score(
        fasta_path=str(combined_fasta),
        flashfry_jar=flashfry_jar,
        flashfry_db=flashfry_db,
        max_mismatch=args.max_mismatch,
        scoring_metrics=args.scoring_metrics,
        executor=executor,
        outdir=str(outdir),
        n_chunks=n_workers,
    )

    # ---- 3. Annotate ----
    annotated = outdir / "guides.annotated.tsv"
    annotate_args = [
        "--scored_chunks", *scored_chunks,
        "--region_metadata", str(combined_meta),
        "--blast_radius_bp", str(args.blast_radius_bp),
        "--out", str(annotated),
    ]
    if atac_narrowpeak_path:
        annotate_args += ["--atac_narrowpeak", atac_narrowpeak_path]
    _run_py_script("annotate.py", annotate_args)

    annotated_for_rank = annotated

    # ---- 4. Optional launch-dCas9 scoring ----
    if args.run_launch_dcas9:
        venv_python = cfg.require(conf, "launch_dcas9", "venv_python")
        launch_script = cfg.require(conf, "launch_dcas9", "script")
        ckpt_dir = conf.get("launch_dcas9", {}).get("ckpt_dir") or str(
            files("crispri_guide_designer").joinpath("checkpoints"))

        launch_dir = outdir / "launch_dcas9"
        chunk_dir = launch_dir / "chunks"
        enh_dir = launch_dir / "results_enhancer"
        pro_dir = launch_dir / "results_promoter"
        for d in (chunk_dir, enh_dir, pro_dir):
            d.mkdir(parents=True, exist_ok=True)

        _run_py_script("launch_dcas9_split.py", [
            "--annotated", str(annotated), "--chunk_size", str(args.launch_chunk_size),
            "--outdir", str(chunk_dir),
        ])
        n_chunks = int((chunk_dir / "n_chunks.txt").read_text().strip())

        commands, labels = [], []
        for outcome, result_dir in (("enhancerFitness", enh_dir), ("promoterFitness", pro_dir)):
            for i in range(n_chunks):
                chunk_csv = chunk_dir / f"chunk_{i:04d}.csv"
                result_csv = result_dir / f"chunk_{i:04d}.csv"
                cmd = (
                    f"{venv_python} {_script_path('launch_dcas9_score_chunk.py')} "
                    f"--chunk_csv {chunk_csv} --outcome {outcome} "
                    f"--launch_script {launch_script} --ckpt_dir {ckpt_dir} "
                    f"--tmp_dir {result_dir}/tmp_{i:04d} --result_csv {result_csv} "
                    f"--python_bin {venv_python}"
                )
                commands.append(cmd)
                labels.append(f"launch_{outcome}_{i:04d}")

        results = executor.run_batch(commands, logdir=str(launch_dir / "logs"), labels=labels)
        failed = [r for r in results if r.returncode != 0]
        if failed:
            raise SystemExit(f"ERROR: {len(failed)} launch-dCas9 chunk(s) failed, see {launch_dir / 'logs'}")

        merged = outdir / "guides.annotated.with_launch.tsv"
        _run_py_script("launch_dcas9_merge.py", [
            "--annotated", str(annotated), "--enhancer_dir", str(enh_dir),
            "--promoter_dir", str(pro_dir), "--out", str(merged),
        ])
        annotated_for_rank = merged

    # ---- 5. Rank ----
    ranked = outdir / "guides.ranked.tsv"
    _run_py_script("rank_report.py", [
        "--annotated", str(annotated_for_rank), "--max_mismatch", str(args.max_mismatch),
        "--out", str(ranked),
    ])

    final_output = ranked

    # ---- 6. Optional seed-finder-lite off-target annotation ----
    if args.run_seed_finder:
        seed_finder_dir = cfg.require(conf, "seed_finder", "dir")
        seed_finder_python = conf.get("seed_finder", {}).get("python_bin") or None

        ranked_df = pd.read_csv(ranked, sep="\t")
        summary, raw_hits_path, guide_ids = seed_finder_integrate.run_seed_finder(
            ranked_df, seed_finder_dir, args.seed_finder_min_length, args.seed_finder_tss_window,
            str(outdir), python_bin=seed_finder_python)

        ranked_df["_guide_id"] = guide_ids
        summary = summary.rename(columns={"guide_id": "_guide_id"})
        final_df = ranked_df.merge(summary, on="_guide_id", how="left").drop(columns=["_guide_id"])

        final_output = outdir / "guides.final.tsv"
        final_df.to_csv(final_output, sep="\t", index=False)
        print(f"seed-finder-lite: raw hit-level detail -> {raw_hits_path}", file=sys.stderr)

    print(f"\nDone. Final output: {final_output}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(prog="cretarget")
    sub = parser.add_subparsers(dest="command", required=True)

    p_run = sub.add_parser("run", help="design, score, and rank guides")
    p_run.add_argument("--snp-bed", dest="snp_bed")
    p_run.add_argument("--enhancer-bed", dest="enhancer_bed")
    p_run.add_argument("--dedup-blast-radius", dest="dedup_blast_radius", action="store_true")
    p_run.add_argument("--atac-narrowpeak", dest="atac_narrowpeak", default=None,
                        help="optional ENCODE-standard 10-column narrowPeak file (single cell type); "
                             "adds distance_to_atac_summit column, informational only (not used in ranking yet). "
                             "Mutually exclusive with --atac-cell-type.")
    p_run.add_argument("--atac-cell-type", dest="atac_cell_type", default=None,
                        choices=sorted(atac_cache.CELL_TYPES.keys()),
                        help="use bundled ATAC-seq peaks for this cell type instead of supplying your own file "
                             "(auto-downloaded and cached on first use). Mutually exclusive with --atac-narrowpeak.")
    p_run.add_argument("--blast-radius-bp", dest="blast_radius_bp", type=int, default=1000)
    p_run.add_argument("--flank-bp", dest="flank_bp", type=int, default=50)
    p_run.add_argument("--snp-window-bp", dest="snp_window_bp", type=int, default=150)
    p_run.add_argument("--max-mismatch", dest="max_mismatch", type=int, default=3)
    p_run.add_argument("--scoring-metrics", dest="scoring_metrics",
                        default="doench2014ontarget,doench2016cfd,dangerous,hsu2013,moreno2015,jostandsantos,minot,folding")
    p_run.add_argument("--outdir", default="results")
    p_run.add_argument("--executor", choices=["local", "slurm"], default="local")
    p_run.add_argument("--n-workers", dest="n_workers", type=int, default=None)
    p_run.add_argument("--run-launch-dcas9", dest="run_launch_dcas9", action="store_true")
    p_run.add_argument("--launch-chunk-size", dest="launch_chunk_size", type=int, default=200000)
    p_run.add_argument("--run-seed-finder", dest="run_seed_finder", action="store_true",
                        help="run seed-finder-lite (off-target seed/TSS-proximity search) as a final step, "
                             "adding summary columns to the output. Requires --seed-finder-min-length and "
                             "--seed-finder-tss-window, and [seed_finder].dir set in config.")
    p_run.add_argument("--seed-finder-min-length", dest="seed_finder_min_length", type=int, default=None,
                        help="minimum exact-match seed length for seed-finder-lite (e.g. 8). No default -- "
                             "you must choose this deliberately (see seed-finder-lite's own documentation).")
    p_run.add_argument("--seed-finder-tss-window", dest="seed_finder_tss_window", type=int, default=None,
                        help="max distance (bp) from a TSS to count as a hit for seed-finder-lite (e.g. 1000). "
                             "No default -- you must choose this deliberately.")
    p_run.add_argument("--fasta", default=None, help="override config's [reference].fasta")
    p_run.add_argument("--flashfry-jar", dest="flashfry_jar", default=None)
    p_run.add_argument("--flashfry-db", dest="flashfry_db", default=None)
    p_run.set_defaults(func=cmd_run)

    p_conf = sub.add_parser("configure", help="show/edit config file")
    p_conf.add_argument("--set", action="append", default=[], help="section.key=value, repeatable")
    p_conf.set_defaults(func=cmd_configure)

    p_ref = sub.add_parser("setup-reference", help="download hg38 + FlashFry index from Zenodo")
    p_ref.add_argument("--outdir", default=str(Path.home() / ".guidenet_reference"))
    p_ref.add_argument("--force", action="store_true")
    p_ref.set_defaults(func=cmd_setup_reference)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
