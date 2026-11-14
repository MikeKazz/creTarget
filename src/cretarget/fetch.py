"""
fetch.py

Downloads and caches the assets creTarget needs but does not ship:
model checkpoints, ATAC-seq peak files, and the hg38 reference + FlashFry
index.

Design constraint that shapes everything here: **compute nodes on many
clusters have no outbound network access.** So fetching is always an
explicit, user-invoked, login-node step (`cretarget fetch-models`,
`fetch-atac`, `setup-reference`) that populates a cache directory. Pipeline
steps only ever *read* from that cache, and if an asset is missing they
fail with a message telling the user which command to run from a networked
node -- never by silently attempting a download mid-run.
"""

import gzip
import hashlib
import json
import shutil
import subprocess
import sys
import tarfile
import urllib.error
import urllib.request
from pathlib import Path

from . import config as cfg

# --- Model checkpoints ----------------------------------------------------
# launch-dCas9 pretrained weights. These are NOT ours and are NOT bundled in
# the PyPI package; they are redistributed from our Zenodo record with
# attribution to the original authors (see MODEL_ATTRIBUTION below).
MODELS_ZENODO_RECORD = "22729308"
MODELS_DOI = "10.5281/zenodo.22729308"
MODELS_BASE_URL = f"https://zenodo.org/records/{MODELS_ZENODO_RECORD}/files"

MODEL_FILES = {
    "ckpt-seq-enhancerFitness.pth":
        "295ef0fb318e5aa1c25d74aa3e18fab0ca5fa9e733438d0b04ff6a12f04908a3",
    "ckpt-seq-promoterFitness.pth":
        "b78fbc1f06b4e39573a36695ec493458c9a816984b52e6b2fdc8e38638c0ea24",
}

MODEL_ATTRIBUTION = """\
These model checkpoints are the pretrained weights from launch-dCas9
(https://github.com/Wancen/launch-dCas9), redistributed via Zenodo
({doi}) for convenience. They are the work of the launch-dCas9
authors, not of creTarget, and remain under their original license and
terms. If you use the --run_launch_dcas9 step, please cite launch-dCas9 in
addition to creTarget.\
"""

# --- ATAC-seq peaks -------------------------------------------------------
ATAC_ZENODO_RECORD = "21268805"
ATAC_DOI = "10.5281/zenodo.21268805"
ATAC_BASE_URL = f"https://zenodo.org/records/{ATAC_ZENODO_RECORD}/files"

ATAC_CELL_TYPES = [
    "b_cell", "basophil", "cd34_bone_marrow", "cd34_cord_blood", "cd4_tcell",
    "cd8_tcell", "clp", "cmp", "erythroid", "gmp", "hsc", "k562", "lmpp",
    "megakaryocyte", "mep", "monocyte", "mpp", "neutrophil", "nk_cell",
]

# --- Reference genome + FlashFry index ------------------------------------
REFERENCE_ZENODO_RECORD = "21115562"
REFERENCE_DOI = "10.5281/zenodo.21115562"
REFERENCE_BASE_URL = f"https://zenodo.org/records/{REFERENCE_ZENODO_RECORD}/files"


class MissingAssetError(SystemExit):
    """Raised when a cached asset is absent and we refuse to auto-download."""


def models_dir(config: dict | None = None) -> Path:
    return cfg.cache_dir(config) / "models"


def atac_dir(config: dict | None = None) -> Path:
    return cfg.cache_dir(config) / "atac_peaks"


def reference_dir(config: dict | None = None) -> Path:
    return cfg.cache_dir(config) / "reference"


def sha256_of(path: Path, chunk_bytes: int = 1 << 20) -> str:
    """Streaming hash -- the reference FASTA is ~3GB, so read_bytes() is not
    an option."""
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk_bytes), b""):
            digest.update(block)
    return digest.hexdigest()


def _download(url: str, dest: Path, description: str = "") -> None:
    """Downloads to a .part file then renames, so an interrupted transfer
    never leaves something that looks like a complete, valid asset."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    label = description or dest.name
    print(f"  -> downloading {label}", file=sys.stderr)
    try:
        with urllib.request.urlopen(url) as response, open(tmp, "wb") as out:
            shutil.copyfileobj(response, out)
    except urllib.error.URLError as exc:
        tmp.unlink(missing_ok=True)
        raise SystemExit(
            f"ERROR: could not download {label}: {exc}\n"
            f"  URL: {url}\n"
            f"If this machine has no outbound network access (common on cluster\n"
            f"compute nodes), run this command from a login node or download the\n"
            f"file manually and place it in the cache directory."
        ) from exc
    tmp.replace(dest)


def _verify(path: Path, expected_sha256: str, label: str) -> None:
    actual = sha256_of(path)
    if actual != expected_sha256:
        path.unlink(missing_ok=True)
        raise SystemExit(
            f"ERROR: checksum mismatch for {label}\n"
            f"  expected: {expected_sha256}\n"
            f"  actual:   {actual}\n"
            f"The partial/corrupt file has been removed; re-run to try again."
        )


# --------------------------------------------------------------------------
# Model checkpoints
# --------------------------------------------------------------------------

def fetch_models(dest: Path | None = None, force: bool = False,
                 config: dict | None = None) -> Path:
    dest = Path(dest) if dest else models_dir(config)
    dest.mkdir(parents=True, exist_ok=True)

    print(f"=== launch-dCas9 checkpoints (Zenodo {MODELS_DOI}) -> {dest} ===",
          file=sys.stderr)
    for filename, expected in MODEL_FILES.items():
        target = dest / filename
        if target.exists() and not force:
            if sha256_of(target) == expected:
                print(f"  -> {filename} already present and verified", file=sys.stderr)
                continue
            print(f"  -> {filename} present but checksum differs, re-downloading",
                  file=sys.stderr)
        _download(f"{MODELS_BASE_URL}/{filename}?download=1", target, filename)
        _verify(target, expected, filename)

    print("", file=sys.stderr)
    print(MODEL_ATTRIBUTION.format(doi=MODELS_DOI), file=sys.stderr)
    return dest


def require_models(ckpt_dir: str | None = None, config: dict | None = None) -> Path:
    """Returns the checkpoint directory, or fails with instructions.

    Never downloads: this is called from inside pipeline steps, which may be
    running on a node with no network.
    """
    path = Path(cfg.expand(ckpt_dir)) if ckpt_dir else models_dir(config)
    missing = [name for name in MODEL_FILES if not (path / name).exists()]
    if missing:
        raise MissingAssetError(
            f"ERROR: launch-dCas9 checkpoint(s) not found in {path}: {', '.join(missing)}\n"
            f"Fetch them once, from a machine with internet access:\n"
            f"    cretarget fetch-models\n"
            f"If your compute nodes cannot reach the internet, run that on a login\n"
            f"node with a shared cache location, e.g.:\n"
            f"    CRETARGET_CACHE_DIR=/shared/path/.cretarget cretarget fetch-models"
        )
    return path


# --------------------------------------------------------------------------
# ATAC-seq peaks
# --------------------------------------------------------------------------

def _atac_checksums(dest: Path, force: bool = False) -> dict:
    sha_path = dest / "checksums.sha256"
    if not sha_path.exists() or force:
        _download(f"{ATAC_BASE_URL}/checksums.sha256?download=1", sha_path,
                  "checksums.sha256")
    checksums = {}
    for line in sha_path.read_text().splitlines():
        parts = line.split()
        if len(parts) == 2:
            checksums[parts[1]] = parts[0]
    return checksums


def fetch_atac(cell_type: str, dest: Path | None = None, force: bool = False,
               config: dict | None = None) -> Path:
    if cell_type not in ATAC_CELL_TYPES:
        raise SystemExit(
            f"ERROR: unknown ATAC cell type {cell_type!r}.\n"
            f"Valid options: {', '.join(ATAC_CELL_TYPES)}"
        )
    dest = Path(dest) if dest else atac_dir(config)
    dest.mkdir(parents=True, exist_ok=True)

    filename = f"{cell_type}.narrowPeak"
    target = dest / filename
    if target.exists() and not force:
        print(f"  -> {filename} already present", file=sys.stderr)
        return target

    checksums = _atac_checksums(dest, force=force)
    print(f"=== ATAC peaks '{cell_type}' (Zenodo {ATAC_DOI}) -> {dest} ===",
          file=sys.stderr)
    _download(f"{ATAC_BASE_URL}/{filename}?download=1", target, filename)

    expected = checksums.get(filename)
    if expected:
        _verify(target, expected, filename)
    else:
        print(f"WARNING: no checksum listed for {filename}, skipping verification",
              file=sys.stderr)
    return target


def require_atac(cell_type: str, config: dict | None = None) -> Path:
    """Cache-only lookup for use inside pipeline steps."""
    target = atac_dir(config) / f"{cell_type}.narrowPeak"
    if not target.exists():
        raise MissingAssetError(
            f"ERROR: ATAC peaks for '{cell_type}' not found at {target}\n"
            f"Fetch them once, from a machine with internet access:\n"
            f"    cretarget fetch-atac --cell-type {cell_type}"
        )
    return target


# --------------------------------------------------------------------------
# Reference genome + FlashFry index
# --------------------------------------------------------------------------

def setup_reference(dest: Path | None = None, force: bool = False,
                    config: dict | None = None) -> Path:
    dest = Path(cfg.expand(str(dest))) if dest else reference_dir(config)
    dest.mkdir(parents=True, exist_ok=True)

    print(f"=== hg38 + FlashFry index (Zenodo {REFERENCE_DOI}) -> {dest} ===",
          file=sys.stderr)
    print("    ~3.2GB download; this takes a while.", file=sys.stderr)

    pairs = [
        ("hg38.fa.gz", "hg38.fa.gz.sha256"),
        ("hg38_standard_cas9ngg_db.tar.gz", "hg38_standard_cas9ngg_db.tar.gz.sha256"),
    ]
    for payload, sha_file in pairs:
        sha_path = dest / sha_file
        if not sha_path.exists() or force:
            _download(f"{REFERENCE_BASE_URL}/{sha_file}?download=1", sha_path, sha_file)
        expected = sha_path.read_text().split()[0]

        target = dest / payload
        if target.exists() and not force:
            print(f"  -> {payload} already present, verifying", file=sys.stderr)
            if sha256_of(target) == expected:
                continue
            print(f"  -> checksum differs, re-downloading {payload}", file=sys.stderr)
        _download(f"{REFERENCE_BASE_URL}/{payload}?download=1", target, payload)
        _verify(target, expected, payload)

    fasta = dest / "hg38.fa"
    if not fasta.exists() or force:
        # Decompressed to plain text deliberately: some bedtools builds
        # mishandle bgzipped FASTA input.
        print("  -> decompressing hg38.fa.gz", file=sys.stderr)
        _require_tool("bgzip", "htslib")
        subprocess.run(["bgzip", "-d", "-k", "-f", str(dest / "hg38.fa.gz")], check=True)

    if not (dest / "hg38.fa.fai").exists() or force:
        print("  -> building .fai index (samtools faidx)", file=sys.stderr)
        _require_tool("samtools", "samtools")
        subprocess.run(["samtools", "faidx", str(fasta)], check=True)

    db_file = dest / "hg38_standard_cas9ngg_db"
    if not db_file.exists() or not (dest / "hg38_standard_cas9ngg_db.header").exists() or force:
        print("  -> extracting FlashFry index", file=sys.stderr)
        with tarfile.open(dest / "hg38_standard_cas9ngg_db.tar.gz") as tar:
            _safe_extract(tar, dest)

    if not db_file.exists() or not (dest / "hg38_standard_cas9ngg_db.header").exists():
        raise SystemExit(
            "ERROR: expected both hg38_standard_cas9ngg_db and "
            "hg38_standard_cas9ngg_db.header after extraction."
        )

    print("", file=sys.stderr)
    print("Done. Record these paths:", file=sys.stderr)
    print(f"    cretarget configure --set reference.fasta={fasta}", file=sys.stderr)
    print(f"    cretarget configure --set reference.flashfry_db={db_file}", file=sys.stderr)
    print("FlashFry itself is a separate download (JVM tool):", file=sys.stderr)
    print("    https://github.com/mckennalab/FlashFry/releases", file=sys.stderr)
    return dest


def _safe_extract(tar: tarfile.TarFile, dest: Path) -> None:
    """Guards against path-traversal entries in the tarball ("../../etc/x").

    Python 3.12 added filter="data" for exactly this, but we support 3.10+,
    so the check is explicit.
    """
    dest = dest.resolve()
    for member in tar.getmembers():
        target = (dest / member.name).resolve()
        if not str(target).startswith(str(dest)):
            raise SystemExit(
                f"ERROR: refusing to extract {member.name!r} -- it would write "
                f"outside {dest}"
            )
    tar.extractall(dest)


# --------------------------------------------------------------------------
# Promoter-proximal seed reference
#
# Built rather than hosted. The two files seed screening needs -- a TSS table
# and a FASTA of promoter windows -- are derived entirely from a GENCODE
# annotation plus the hg38 FASTA that creTarget already requires, so there is
# nothing to host, nothing to keep alive, and the exact annotation release is
# recorded in the run manifest instead of being an opaque property of a
# downloaded blob.
#
# The window recipe reproduces seed-finder-lite's published reference exactly
# (verified against it for all 24 chromosomes): take +/-2 kb around every
# annotated TSS, merge overlapping windows per chromosome, then pad each merged
# window by 23 bp -- one protospacer plus PAM -- so a site straddling a
# boundary is not truncated. chrM is excluded: at 16.6 kb, +/-2 kb windows
# would blanket most of it, and CRISPRi on mtDNA is not a meaningful screen.
# --------------------------------------------------------------------------

GENCODE_RELEASE = "44"
GENCODE_GTF_URL = (
    f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
    f"release_{GENCODE_RELEASE}/gencode.v{GENCODE_RELEASE}.annotation.gtf.gz"
)

SEED_REFERENCE_ATTRIBUTION = """\
The promoter-proximal seed screen, and the design of the promoter-window
reference it scans, are from seed-finder-lite (Kazemi, M.;
doi:10.5281/zenodo.21599135). Gene annotations are GENCODE, which is
distributed under its own terms -- see https://www.gencodegenes.org/.
"""

# Standard assembled chromosomes. Scaffolds and patches are left out: a seed
# hit on an unplaced contig is not actionable, and including them would inflate
# the scan for no gain.
SEED_REFERENCE_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


def seed_reference_dir(config: dict | None = None) -> Path:
    return cfg.cache_dir(config) / "seed_reference"


def _parse_gencode_tss(gtf_path: Path):
    """(chrom, tss, gene_name) for every gene entry.

    The TSS is the 5'-most base: the start for a + strand gene, the end for a
    - strand one. Getting that backwards would put every minus-strand gene's
    promoter at the wrong end of the gene.
    """
    rows = []
    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    with opener(gtf_path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            chrom = fields[0]
            if not chrom.startswith("chr"):
                chrom = "chr" + chrom
            start, end, strand = int(fields[3]), int(fields[4]), fields[6]
            tss = start if strand == "+" else end
            gene_name = "?"
            for item in fields[8].split(";"):
                item = item.strip()
                if item.startswith("gene_name "):
                    gene_name = item.split(" ", 1)[1].strip().strip('"')
                    break
            rows.append((chrom, tss, gene_name))
    return rows


def _promoter_windows(tss_rows, chroms, half_width, pad):
    """Merged, padded promoter windows as (chrom, start_1based, end_1based)."""
    by_chrom: dict[str, list[tuple[int, int]]] = {}
    allowed = set(chroms)
    for chrom, tss, _gene in tss_rows:
        if chrom not in allowed:
            continue
        by_chrom.setdefault(chrom, []).append(
            (max(1, tss - half_width), tss + half_width))

    windows = []
    for chrom in chroms:
        intervals = by_chrom.get(chrom)
        if not intervals:
            continue
        intervals.sort()
        current_start, current_end = intervals[0]
        merged = []
        for start, end in intervals[1:]:
            # +1 because these are inclusive: windows that abut should merge.
            if start <= current_end + 1:
                current_end = max(current_end, end)
            else:
                merged.append((current_start, current_end))
                current_start, current_end = start, end
        merged.append((current_start, current_end))
        for start, end in merged:
            windows.append((chrom, max(1, start - pad), end + pad))
    return windows


def _read_fai(fasta: Path) -> dict:
    fai = Path(str(fasta) + ".fai")
    if not fai.exists():
        raise SystemExit(
            f"ERROR: no FASTA index at {fai}\n"
            f"Run: samtools faidx {fasta}"
        )
    lengths = {}
    with open(fai) as handle:
        for line in handle:
            fields = line.split("\t")
            if len(fields) >= 2:
                lengths[fields[0]] = int(fields[1])
    return lengths


def build_seed_reference(fasta: str, gtf: str | None = None,
                         dest: Path | None = None, force: bool = False,
                         config: dict | None = None) -> Path:
    """Build the TSS table and promoter-window FASTA for seed screening."""
    from . import seed_scan

    dest = Path(dest) if dest else seed_reference_dir(config)
    dest.mkdir(parents=True, exist_ok=True)
    tss_path = dest / "tss_table.tsv"
    promoter_path = dest / "promoter_windows.fa.gz"
    provenance_path = dest / "provenance.json"

    if not force and tss_path.exists() and promoter_path.exists():
        print(f"Seed reference already present in {dest} (use --force to rebuild)",
              file=sys.stderr)
        return dest

    _require_tool("bedtools", "bedtools")
    # Resolved so the provenance record names an unambiguous path rather than
    # whatever was relative to the caller's working directory.
    fasta = Path(cfg.expand(str(fasta))).resolve()
    if not fasta.is_file():
        raise SystemExit(f"ERROR: genome FASTA not found: {fasta}")
    chrom_lengths = _read_fai(fasta)

    if gtf:
        gtf_path = Path(cfg.expand(str(gtf)))
        gtf_source = str(gtf_path)
        if not gtf_path.is_file():
            raise SystemExit(f"ERROR: GTF not found: {gtf_path}")
    else:
        gtf_path = reference_dir(config) / f"gencode.v{GENCODE_RELEASE}.annotation.gtf.gz"
        gtf_source = GENCODE_GTF_URL
        if not gtf_path.exists():
            _download(GENCODE_GTF_URL, gtf_path,
                      f"GENCODE release {GENCODE_RELEASE} annotation (~50 MB)")

    print(f"Reading gene annotations from {gtf_path.name}", file=sys.stderr)
    tss_rows = _parse_gencode_tss(gtf_path)
    if not tss_rows:
        raise SystemExit(
            f"ERROR: no gene entries found in {gtf_path}\n"
            f"Is this a GENCODE/Ensembl GTF? A GFF3 will not parse here."
        )
    tss_rows.sort(key=lambda r: (r[0], r[1]))
    with open(tss_path, "w") as handle:
        handle.write("chr\ttss\tgene_name\n")
        for chrom, tss, gene in tss_rows:
            handle.write(f"{chrom}\t{tss}\t{gene}\n")
    print(f"  {len(tss_rows):,} TSSs -> {tss_path.name}", file=sys.stderr)

    windows = _promoter_windows(
        tss_rows, SEED_REFERENCE_CHROMS,
        seed_scan.PROMOTER_HALF_WIDTH, seed_scan.PROMOTER_PAD)
    if not windows:
        raise SystemExit(
            "ERROR: no promoter windows built. The GTF's chromosome names did "
            "not match any of chr1..chr22, chrX, chrY."
        )

    # BED is 0-based half-open; the window name carries 1-based inclusive
    # coordinates, which is what the scanner parses to remap hits.
    bed_path = dest / "promoter_windows.bed"
    skipped = 0
    with open(bed_path, "w") as handle:
        for chrom, start_1, end_1 in windows:
            limit = chrom_lengths.get(chrom)
            if limit is None:
                skipped += 1
                continue
            clipped_end = min(end_1, limit)
            if clipped_end - start_1 + 1 < seed_scan.SITE_LENGTH:
                continue
            handle.write(
                f"{chrom}\t{start_1 - 1}\t{clipped_end}\t"
                f"prom|{chrom}|{start_1}|{clipped_end}\n")
    if skipped:
        print(f"  note: {skipped} window(s) on chromosomes absent from the "
              f"FASTA were skipped", file=sys.stderr)

    print(f"  extracting {len(windows):,} promoter windows with bedtools",
          file=sys.stderr)
    plain_fa = dest / "promoter_windows.fa"
    with open(plain_fa, "w") as out:
        result = subprocess.run(
            ["bedtools", "getfasta", "-fi", str(fasta), "-bed", str(bed_path),
             "-nameOnly"],
            stdout=out, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        plain_fa.unlink(missing_ok=True)
        raise SystemExit(
            f"ERROR: bedtools getfasta failed (exit {result.returncode}):\n"
            f"{result.stderr}")

    with open(plain_fa, "rb") as src, gzip.open(promoter_path, "wb") as out:
        shutil.copyfileobj(src, out)
    plain_fa.unlink()
    bed_path.unlink()

    provenance = {
        "gencode_release": None if gtf else GENCODE_RELEASE,
        "gtf_source": gtf_source,
        "gtf_sha256": sha256_of(gtf_path),
        "genome_fasta": str(fasta),
        "tss_count": len(tss_rows),
        "promoter_window_count": len(windows),
        "promoter_half_width_bp": seed_scan.PROMOTER_HALF_WIDTH,
        "promoter_pad_bp": seed_scan.PROMOTER_PAD,
        "chromosomes": SEED_REFERENCE_CHROMS,
        "cretarget_version": __import__("cretarget").__version__,
    }
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n")

    print(f"\nSeed reference ready in {dest}", file=sys.stderr)
    print(SEED_REFERENCE_ATTRIBUTION, file=sys.stderr)
    return dest


def require_seed_reference(promoter_fa: str | None = None,
                           tss_table: str | None = None,
                           config: dict | None = None):
    """Resolve the seed reference, explaining how to build it if absent."""
    if promoter_fa and tss_table:
        promoter, tss = Path(cfg.expand(promoter_fa)), Path(cfg.expand(tss_table))
    else:
        base = seed_reference_dir(config)
        promoter = Path(cfg.expand(promoter_fa)) if promoter_fa else base / "promoter_windows.fa.gz"
        tss = Path(cfg.expand(tss_table)) if tss_table else base / "tss_table.tsv"
    missing = [str(p) for p in (promoter, tss) if not p.is_file()]
    if missing:
        raise SystemExit(
            "ERROR: the promoter-proximal seed reference is not built.\n"
            f"  missing: {', '.join(missing)}\n\n"
            "Build it once with:\n"
            "  cretarget build-seed-reference --fasta /path/to/hg38.fa\n\n"
            "It needs bedtools and a GENCODE GTF (downloaded automatically, or "
            "pass --gtf)."
        )
    return promoter, tss


def _require_tool(name: str, package_hint: str) -> None:
    if shutil.which(name) is None:
        raise SystemExit(
            f"ERROR: '{name}' not found on PATH (install {package_hint}, e.g. "
            f"`conda install -c bioconda {package_hint}`)"
        )
