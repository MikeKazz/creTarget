"""
config.py

Replaces local_paths.env now that this is a pip-installed tool (no more
"script directory" to keep a config file next to). Config lives at
~/.config/cretarget/config.toml, managed via `cretarget configure`
or manual edit.
"""

import os
from pathlib import Path

try:
    import tomllib  # Python 3.11+
except ImportError:
    import tomli as tomllib  # type: ignore

CONFIG_DIR = Path(os.environ.get("CRETARGET_CONFIG_DIR", Path.home() / ".config" / "cretarget"))
CONFIG_PATH = CONFIG_DIR / "config.toml"

DEFAULT_CONFIG = """\
# cretarget config
# Edit this file, or run `cretarget configure` for a guided setup.

[reference]
fasta = ""             # e.g. ~/.guidenet_reference/hg38.fa (see `cretarget setup-reference`)
flashfry_jar = ""       # path to FlashFry-assembly-*.jar (not hosted by us, install separately)
flashfry_db = ""        # e.g. ~/.guidenet_reference/hg38_standard_cas9ngg_db

[launch_dcas9]
venv_python = ""        # optional: path to launch-dCas9 conda env's python
script = ""             # optional: path to launch-dCas9.py
ckpt_dir = ""           # optional: defaults to the checkpoints/ bundled with this package

[slurm]
partition = ""
gres = ""
mem = "16G"
cpus = "8"
time = "2:00:00"
max_array_size = "1000"   # SLURM cluster's MaxArraySize -- check with `scontrol show config | grep MaxArraySize`
array_throttle = ""        # optional: max concurrently-running array tasks (SLURM %N syntax), blank = no throttle

[seed_finder]
dir = ""                # path to seed-finder-lite-main (the dir containing src/ and refs/)
python_bin = ""         # optional: python with pyfaidx installed; defaults to whatever runs cretarget itself
"""


def ensure_config_exists():
    CONFIG_DIR.mkdir(parents=True, exist_ok=True)
    if not CONFIG_PATH.exists():
        CONFIG_PATH.write_text(DEFAULT_CONFIG)
    return CONFIG_PATH


def load_config() -> dict:
    ensure_config_exists()
    with open(CONFIG_PATH, "rb") as f:
        return tomllib.load(f)


def require(config: dict, section: str, key: str) -> str:
    val = config.get(section, {}).get(key, "")
    if not val:
        raise SystemExit(
            f"ERROR: [{section}].{key} is not set in {CONFIG_PATH}\n"
            f"Run `cretarget configure` or edit that file directly."
        )
    return val
