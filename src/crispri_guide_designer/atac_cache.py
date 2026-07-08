"""
atac_cache.py

Maps a cell-type slug (e.g. "b_cell") to its ATAC-seq narrowPeak file,
hosted on Zenodo (DOI 10.5281/zenodo.21268805). Downloads on first use into
a local cache dir, verifies against the project's own checksums.sha256
(also hosted in the same Zenodo record), and returns the local path -- used
exactly like a user-supplied --atac-narrowpeak file downstream.
"""

import os
import sys
import urllib.request
from pathlib import Path

ZENODO_RECORD_ID = "21268805"
BASE_URL = f"https://zenodo.org/records/{ZENODO_RECORD_ID}/files"

# Slug -> filename (both identical by construction, kept as a dict in case
# that ever changes). Order here also defines --help's displayed choices.
CELL_TYPES = {
    "b_cell": "b_cell.narrowPeak",
    "cd34_cord_blood": "cd34_cord_blood.narrowPeak",
    "cd34_bone_marrow": "cd34_bone_marrow.narrowPeak",
    "cd4_tcell": "cd4_tcell.narrowPeak",
    "cd8_tcell": "cd8_tcell.narrowPeak",
    "clp": "clp.narrowPeak",
    "cmp": "cmp.narrowPeak",
    "erythroid": "erythroid.narrowPeak",
    "gmp": "gmp.narrowPeak",
    "hsc": "hsc.narrowPeak",
    "k562": "k562.narrowPeak",
    "lmpp": "lmpp.narrowPeak",
    "mep": "mep.narrowPeak",
    "mpp": "mpp.narrowPeak",
    "megakaryocyte": "megakaryocyte.narrowPeak",
    "monocyte": "monocyte.narrowPeak",
    "nk_cell": "nk_cell.narrowPeak",
    "basophil": "basophil.narrowPeak",
    "neutrophil": "neutrophil.narrowPeak",
}


def default_cache_dir() -> Path:
    return Path.home() / ".guidenet_reference" / "atac_peaks"


def _load_checksums(cache_dir: Path) -> dict:
    """Downloads (or reuses cached) checksums.sha256, returns {filename: sha256}."""
    sha_path = cache_dir / "checksums.sha256"
    if not sha_path.exists():
        cache_dir.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(f"{BASE_URL}/checksums.sha256?download=1", str(sha_path))
    checksums = {}
    for line in sha_path.read_text().splitlines():
        parts = line.split()
        if len(parts) == 2:
            checksums[parts[1]] = parts[0]
    return checksums


def get_atac_path(cell_type: str, cache_dir: Path | None = None, force: bool = False) -> str:
    """Returns the local path to the given cell type's narrowPeak file,
    downloading it (and verifying checksum) on first use."""
    import hashlib

    if cell_type not in CELL_TYPES:
        valid = ", ".join(sorted(CELL_TYPES))
        raise SystemExit(f"ERROR: unknown --atac-cell-type '{cell_type}'. Valid options: {valid}")

    cache_dir = cache_dir or default_cache_dir()
    cache_dir.mkdir(parents=True, exist_ok=True)
    filename = CELL_TYPES[cell_type]
    dest = cache_dir / filename

    if dest.exists() and not force:
        return str(dest)

    checksums = _load_checksums(cache_dir)
    print(f"Downloading ATAC peaks for '{cell_type}' ({filename})...", file=sys.stderr)
    urllib.request.urlretrieve(f"{BASE_URL}/{filename}?download=1", str(dest))

    expected = checksums.get(filename)
    if expected:
        actual = hashlib.sha256(dest.read_bytes()).hexdigest()
        if actual != expected:
            dest.unlink(missing_ok=True)
            raise SystemExit(f"ERROR: checksum mismatch for {filename} (expected {expected}, got {actual})")
    else:
        print(f"WARNING: no checksum found for {filename}, skipping verification", file=sys.stderr)

    return str(dest)
