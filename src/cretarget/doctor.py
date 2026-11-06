"""
doctor

Preflight checks. The goal is that a user who hits a problem runs
`cretarget doctor` once and sees exactly which piece is missing, instead of
filing an issue about a bedtools error message three steps downstream.

Exit code is 1 if any REQUIRED check fails, 0 otherwise. Optional
components report their status without failing the run.
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

from . import __version__
from . import config as cfg
from . import fetch

OK = "  OK      "
WARN = "  WARN    "
FAIL = "  FAIL    "


def _tool_version(name, flag="--version"):
    """Best-effort version string.

    errors="replace" is not optional: some builds of samtools emit a
    copyright line containing non-UTF-8 bytes, which otherwise raises
    UnicodeDecodeError and takes the whole doctor run down with it. A
    diagnostic command must never be the thing that crashes.
    """
    try:
        proc = subprocess.run(
            [name, flag], capture_output=True, text=True, timeout=20,
            errors="replace",
        )
        output = (proc.stdout or proc.stderr).strip().splitlines()
        return output[0].strip() if output else "(version unknown)"
    except Exception:
        return "(version unknown)"


def _check_tool(name, required, hint, version_flag="--version"):
    path = shutil.which(name)
    if path:
        print(f"{OK}{name}: {_tool_version(name, version_flag)}")
        return True
    marker = FAIL if required else WARN
    status = "required" if required else "optional"
    print(f"{marker}{name}: not found on PATH ({status}) -- {hint}")
    return not required


def _check_file(label, path, required, hint=""):
    if not path:
        marker = FAIL if required else WARN
        print(f"{marker}{label}: not configured{' -- ' + hint if hint else ''}")
        return not required
    if Path(path).exists():
        print(f"{OK}{label}: {path}")
        return True
    marker = FAIL if required else WARN
    print(f"{marker}{label}: configured but missing at {path}")
    return not required


def run(args):
    config = cfg.load_config()
    ok = True

    print(f"creTarget {__version__}")
    print(f"python    {sys.version.split()[0]} ({sys.executable})")
    print(f"config    {cfg.CONFIG_PATH}")
    print(f"cache     {cfg.cache_dir(config)}")
    print()

    print("== External tools ==")
    ok &= _check_tool("bedtools", True, "conda install -c bioconda bedtools")
    ok &= _check_tool("samtools", True, "conda install -c bioconda samtools")
    ok &= _check_tool("java", True, "needed by FlashFry and Nextflow; install a JRE >= 11",
                      version_flag="-version")
    # Nextflow takes a single-dash -version, as does java.
    _check_tool("nextflow", False, "only needed to run the full pipeline; "
                                   "see https://nextflow.io",
                version_flag="-version")
    print()

    print("== Reference data ==")
    fasta = cfg.get_path(config, "reference", "fasta")
    ok &= _check_file("reference.fasta", fasta, True,
                      "run `cretarget setup-reference`, then "
                      "`cretarget configure --set reference.fasta=...`")
    if fasta and Path(fasta).exists():
        ok &= _check_file("reference .fai index", fasta + ".fai", True,
                          f"run `samtools faidx {fasta}`")
    ok &= _check_file("reference.flashfry_jar",
                      cfg.get_path(config, "reference", "flashfry_jar"), True,
                      "download from "
                      "https://github.com/mckennalab/FlashFry/releases")
    ok &= _check_file("reference.flashfry_db",
                      cfg.get_path(config, "reference", "flashfry_db"), True,
                      "run `cretarget setup-reference`")
    print()

    print("== Optional: launch-dCas9 on-target scoring ==")
    models = fetch.models_dir(config)
    present = [n for n in fetch.MODEL_FILES if (models / n).exists()]
    if len(present) == len(fetch.MODEL_FILES):
        print(f"{OK}model checkpoints: {models}")
    else:
        print(f"{WARN}model checkpoints: {len(present)}/{len(fetch.MODEL_FILES)} "
              f"present in {models} -- run `cretarget fetch-models`")
    _check_file("launch_dcas9.script", cfg.get_path(config, "launch_dcas9", "script"),
                False, "clone https://github.com/Wancen/launch-dCas9")
    _check_file("launch_dcas9.python_bin",
                cfg.get_path(config, "launch_dcas9", "python_bin"), False,
                "python from a conda env with torch installed")
    print()

    print("== Optional: seed-finder-lite off-target screening ==")
    _check_file("seed_finder.dir", cfg.get_path(config, "seed_finder", "dir"), False,
                "see Zenodo DOI 10.5281/zenodo.21599135")
    print()

    print("== Optional: ATAC-seq peaks ==")
    atac = fetch.atac_dir(config)
    cached = sorted(p.stem for p in atac.glob("*.narrowPeak")) if atac.exists() else []
    if cached:
        print(f"{OK}cached cell types ({len(cached)}): {', '.join(cached)}")
    else:
        print(f"{WARN}no cached ATAC peaks in {atac} -- "
              f"run `cretarget fetch-atac --cell-type <name>` if you want them")
    print()

    print("== Filesystem ==")
    cache = cfg.cache_dir(config)
    try:
        cache.mkdir(parents=True, exist_ok=True)
        probe = cache / ".cretarget_write_probe"
        probe.write_text("ok")
        probe.unlink()
        print(f"{OK}cache directory is writable: {cache}")
    except OSError as exc:
        print(f"{FAIL}cache directory is not writable: {cache} ({exc})")
        ok = False

    if args.shared_path:
        shared = Path(args.shared_path)
        try:
            shared.mkdir(parents=True, exist_ok=True)
            probe = shared / ".cretarget_write_probe"
            probe.write_text("ok")
            probe.unlink()
            print(f"{OK}--shared-path is writable: {shared}")
            print(f"{WARN}note: creTarget cannot verify from here that compute nodes "
                  f"can see this path. If your cluster nodes do not share a "
                  f"filesystem with this machine, point -w/--outdir somewhere they do.")
        except OSError as exc:
            print(f"{FAIL}--shared-path is not writable: {shared} ({exc})")
            ok = False

    if os.environ.get("CRETARGET_CACHE_DIR"):
        print(f"{OK}CRETARGET_CACHE_DIR is set "
              f"({os.environ['CRETARGET_CACHE_DIR']})")
    print()

    if ok:
        print("All required checks passed.")
    else:
        print("Some REQUIRED checks failed (see FAIL lines above).")
    return 0 if ok else 1


def add_arguments(parser) -> None:
    parser.add_argument("--shared-path", dest="shared_path", default=None,
                        help="also check that this directory exists and is writable "
                             "(use your intended Nextflow work dir on a cluster)")
