"""
config.py

Config for the step-level `cretarget <step>` CLI, at
~/.config/cretarget/config.toml (override the directory with
CRETARGET_CONFIG_DIR).

Scheduling/resource settings are deliberately NOT here any more -- Nextflow
owns execution (see nextflow.config profiles). This file only records
*where things live* on a given machine: reference data, external tools,
cache directories.

Every path read through `get_path()` / `require_path()` is passed through
os.path.expanduser + expandvars, so "~/..." and "$SCRATCH/..." in the TOML
behave the way users expect.
"""

import os
from pathlib import Path

try:
    import tomllib  # Python 3.11+
except ImportError:  # pragma: no cover - exercised only on 3.10
    import tomli as tomllib  # type: ignore

CONFIG_DIR = Path(
    os.environ.get("CRETARGET_CONFIG_DIR", Path.home() / ".config" / "cretarget")
)
CONFIG_PATH = CONFIG_DIR / "config.toml"

DEFAULT_CONFIG = """\
# creTarget config (step-level CLI).
#
# Nextflow runs do NOT read this file -- pass --fasta/--flashfry_jar/etc as
# pipeline params, or put them in a -params-file / custom -c config.
# This file is for `cretarget <step>` invocations and for defaults you get
# tired of retyping.
#
# Paths may use ~ and $ENV_VARS; both are expanded.

[reference]
fasta = ""              # e.g. ~/.cretarget/reference/hg38.fa (see `cretarget setup-reference`)
flashfry_jar = ""       # path to FlashFry-assembly-*.jar (install separately)
flashfry_db = ""        # e.g. ~/.cretarget/reference/hg38_standard_cas9ngg_db

[launch_dcas9]
python_bin = ""         # python from the launch-dCas9 conda env (needs torch)
script = ""             # path to launch-dCas9.py
ckpt_dir = ""           # optional: defaults to the cache dir used by `cretarget fetch-models`

[seed_finder]
dir = ""                # path to seed-finder-lite (the dir containing src/ and refs/)
python_bin = ""         # optional: python with pyfaidx installed; defaults to the current interpreter

[cache]
dir = ""                # optional: overrides CRETARGET_CACHE_DIR / ~/.cretarget
"""


def ensure_config_exists() -> Path:
    CONFIG_DIR.mkdir(parents=True, exist_ok=True)
    if not CONFIG_PATH.exists():
        CONFIG_PATH.write_text(DEFAULT_CONFIG)
    return CONFIG_PATH


def load_config() -> dict:
    ensure_config_exists()
    with open(CONFIG_PATH, "rb") as f:
        return tomllib.load(f)


def expand(value: str) -> str:
    """Expands ~ and $VARS. Empty/None passes through unchanged.

    This is the fix for the single most common first-run failure: the docs
    tell you to set `fasta = "~/..."`, and a literal tilde then gets handed
    to bedtools, which has no idea what to do with it.
    """
    if not value:
        return value
    return os.path.expanduser(os.path.expandvars(str(value)))


def get_path(config: dict, section: str, key: str, default: str = "") -> str:
    return expand(config.get(section, {}).get(key, "") or default)


def require_path(config: dict, section: str, key: str) -> str:
    val = get_path(config, section, key)
    if not val:
        raise SystemExit(
            f"ERROR: [{section}].{key} is not set in {CONFIG_PATH}\n"
            f"Set it with:  cretarget configure --set {section}.{key}=/path/to/thing\n"
            f"(or pass the equivalent command-line flag)."
        )
    return val


def set_values(pairs) -> Path:
    """Applies "section.key=value" strings to the config file.

    Rewritten from the old naive line-replacer, which could not add a key
    that was absent, silently mangled values containing quotes or
    backslashes, and could not represent multi-line values at all.
    We now parse, update, and re-emit the whole document.
    """
    path = ensure_config_exists()
    config = load_config()

    for pair in pairs:
        if "=" not in pair:
            raise SystemExit(f"ERROR: --set expects section.key=value, got {pair!r}")
        section_key, value = pair.split("=", 1)
        if "." not in section_key:
            raise SystemExit(
                f"ERROR: --set expects section.key=value, got {pair!r} "
                f"(missing the section, e.g. reference.fasta=...)"
            )
        section, key = section_key.split(".", 1)
        config.setdefault(section, {})[key] = value

    path.write_text(_dump_toml(config))
    return path


def _toml_escape(value: str) -> str:
    return (
        str(value)
        .replace("\\", "\\\\")
        .replace('"', '\\"')
        .replace("\n", "\\n")
        .replace("\t", "\\t")
    )


def _dump_toml(config: dict) -> str:
    """Minimal TOML writer for our flat section/string schema.

    Deliberately not a general-purpose serializer -- the config is only ever
    sections of string values, and we would rather not add a write-capable
    TOML dependency for this.
    """
    lines = ["# creTarget config -- `cretarget configure --set section.key=value`", ""]
    for section, values in config.items():
        lines.append(f"[{section}]")
        for key, value in values.items():
            lines.append(f'{key} = "{_toml_escape(value)}"')
        lines.append("")
    return "\n".join(lines)


def cache_dir(config: dict | None = None) -> Path:
    """Resolution order: $CRETARGET_CACHE_DIR -> [cache].dir -> ~/.cretarget

    The environment variable wins deliberately. It is the conventional
    per-invocation override, and someone who writes

        CRETARGET_CACHE_DIR=/tmp/scratch cretarget setup-reference

    means it. The earlier order let a configured [cache].dir silently take
    precedence, so that command downloaded 3.2 GB to the configured location
    instead -- an explicit instruction ignored without a word.

    Kept in one place because every downloadable asset (model checkpoints,
    ATAC peaks, reference genome) has to agree on it, and because clusters
    routinely need it pointed somewhere other than $HOME.
    """
    env = os.environ.get("CRETARGET_CACHE_DIR")
    if env:
        return Path(expand(env))
    if config:
        configured = get_path(config, "cache", "dir")
        if configured:
            return Path(configured)
    return Path.home() / ".cretarget"
