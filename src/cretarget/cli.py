"""
cli.py

`cretarget <step>` -- the step-level command line interface.

This is deliberately NOT a pipeline driver. Each subcommand does one step on
one chunk of data, reads its inputs from flags, and writes its outputs to
flags. Nextflow composes them into the real pipeline
(`nextflow run MikeKazz/creTarget`), and owns everything to do with
splitting, scheduling, retries, and resume.

Use these directly to debug a single step, to run something small by hand,
or to script your own composition.

Flags fall back to ~/.config/cretarget/config.toml where it makes sense, so
interactive use does not require retyping reference paths. Nextflow always
passes them explicitly.
"""

import argparse
import sys

from . import __version__
from . import config as cfg
from . import doctor as doctor_mod
from . import fetch
from .steps import annotate as annotate_step
from .steps import coords_to_bed as coords_to_bed_step
from .steps import launch_dcas9 as launch_dcas9_step
from .steps import prep_regions as prep_regions_step
from .steps import rank as rank_step
from .steps import seed_finder as seed_finder_step


def cmd_prep_regions(args):
    config = cfg.load_config()
    args.fasta = cfg.expand(args.fasta) or cfg.require_path(config, "reference", "fasta")
    prep_regions_step.run(args)


def cmd_annotate(args):
    annotate_step.run(args)


def cmd_rank(args):
    rank_step.run(args)


def cmd_score_launch_dcas9(args):
    config = cfg.load_config()
    args.launch_script = (
        cfg.expand(args.launch_script)
        or cfg.require_path(config, "launch_dcas9", "script")
    )
    args.python_bin = (
        cfg.expand(args.python_bin)
        or cfg.get_path(config, "launch_dcas9", "python_bin")
        or sys.executable
    )
    ckpt_dir = cfg.expand(args.ckpt_dir) or cfg.get_path(config, "launch_dcas9", "ckpt_dir")
    # require_models never downloads -- this may be running on a compute node
    # with no network. It fails with instructions instead.
    args.ckpt_dir = str(fetch.require_models(ckpt_dir or None, config))
    launch_dcas9_step.run(args)


def cmd_seed_finder(args):
    config = cfg.load_config()
    args.seed_finder_dir = (
        cfg.expand(args.seed_finder_dir)
        or cfg.require_path(config, "seed_finder", "dir")
    )
    args.python_bin = (
        cfg.expand(args.python_bin)
        or cfg.get_path(config, "seed_finder", "python_bin")
        or None
    )
    seed_finder_step.run(args)


def cmd_coords_to_bed(args):
    coords_to_bed_step.run(args)


def cmd_configure(args):
    if args.set:
        path = cfg.set_values(args.set)
        print(f"Updated {path}")
    else:
        path = cfg.ensure_config_exists()
        print(f"Config file: {path}\n")
        print(path.read_text())
        print("Set values with:")
        print("  cretarget configure --set reference.fasta=/path/to/hg38.fa")


def cmd_fetch_models(args):
    config = cfg.load_config()
    fetch.fetch_models(dest=args.outdir, force=args.force, config=config)


def cmd_fetch_atac(args):
    config = cfg.load_config()
    for cell_type in args.cell_type:
        fetch.fetch_atac(cell_type, dest=args.outdir, force=args.force, config=config)


def cmd_setup_reference(args):
    config = cfg.load_config()
    fetch.setup_reference(dest=args.outdir, force=args.force, config=config)


def cmd_doctor(args):
    raise SystemExit(doctor_mod.run(args))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="cretarget",
        description="CRISPRi guide RNA design for the noncoding genome "
                    "(step-level CLI; use `nextflow run MikeKazz/creTarget` "
                    "for full runs).",
    )
    parser.add_argument("--version", action="version",
                        version=f"creTarget {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("prep-regions",
                       help="parse targets, pad them, extract sequence")
    prep_regions_step.add_arguments(p)
    p.set_defaults(func=cmd_prep_regions)

    p = sub.add_parser("annotate",
                       help="recover coordinates and join region metadata")
    annotate_step.add_arguments(p)
    p.set_defaults(func=cmd_annotate)

    p = sub.add_parser("rank", help="apply hard filters and rank per target")
    rank_step.add_arguments(p)
    p.set_defaults(func=cmd_rank)

    p = sub.add_parser("score-launch-dcas9",
                       help="add launch-dCas9 CNN on-target efficacy scores")
    launch_dcas9_step.add_arguments(p)
    p.set_defaults(func=cmd_score_launch_dcas9)

    p = sub.add_parser("seed-finder",
                       help="add seed-finder-lite off-target/TSS columns")
    seed_finder_step.add_arguments(p)
    p.set_defaults(func=cmd_seed_finder)

    p = sub.add_parser("coords-to-bed",
                       help="extract chr:start-end values from a table into BED")
    coords_to_bed_step.add_arguments(p)
    p.set_defaults(func=cmd_coords_to_bed)

    p = sub.add_parser("configure", help="show or edit the config file")
    p.add_argument("--set", action="append", default=[],
                   help="section.key=value, repeatable")
    p.set_defaults(func=cmd_configure)

    p = sub.add_parser("fetch-models",
                       help="download launch-dCas9 checkpoints into the cache")
    p.add_argument("--outdir", default=None)
    p.add_argument("--force", action="store_true")
    p.set_defaults(func=cmd_fetch_models)

    p = sub.add_parser("fetch-atac", help="download ATAC-seq peaks into the cache")
    p.add_argument("--cell-type", dest="cell_type", action="append", required=True,
                   choices=fetch.ATAC_CELL_TYPES, help="repeatable")
    p.add_argument("--outdir", default=None)
    p.add_argument("--force", action="store_true")
    p.set_defaults(func=cmd_fetch_atac)

    p = sub.add_parser("setup-reference",
                       help="download hg38 + the FlashFry index (~3.2GB)")
    p.add_argument("--outdir", default=None)
    p.add_argument("--force", action="store_true")
    p.set_defaults(func=cmd_setup_reference)

    p = sub.add_parser("doctor", help="check tools, config, and cached data")
    doctor_mod.add_arguments(p)
    p.set_defaults(func=cmd_doctor)

    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
