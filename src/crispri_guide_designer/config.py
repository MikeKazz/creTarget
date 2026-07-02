"""
config.py

Replaces local_paths.env now that this is a pip-installed tool (no more
"script directory" to keep a config file next to). Config lives at
~/.config/crispri-guide-designer/config.toml, managed via `crispri-guide-designer configure`
or manual edit.
"""

import os
from pathlib import Path

try:
    import tomllib  # Python 3.11+
except ImportError:
    import tomli as tomllib  # type: ignore

CONFIG_DIR = Path(os.environ.get("CGD_CONFIG_DIR", Path.home() / ".config" / "crispri-guide-designer"))
CONFIG_PATH = CONFIG_DIR / "config.toml"

DEFAULT_CONFIG = """\
# crispri-guide-designer config
# Edit this file, or run `crispri-guide-designer configure` for a guided setup.

[reference]
fasta = ""             # e.g. ~/.guidenet_reference/hg38.fa (see `crispri-guide-designer setup-reference`)
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
            f"Run `crispri-guide-designer configure` or edit that file directly."
        )
    return val
