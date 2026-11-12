"""
prep_regions

Turns target regions (SNP points and/or enhancer intervals) into a single
padded-region FASTA plus a metadata table.

  1. Parse inputs, auto-detecting format per line.
  2. Validate contigs against the reference .fai (catches hg38-vs-Ensembl
     naming before bedtools does, with a readable message).
  3. Make region names unique (duplicates in the input used to collide
     silently and then blow up the downstream metadata join).
  4. Optionally flag regions within blast_radius_bp of each other as
     "linked" -- one guide could silence both. Flags, never merges.
  5. Pad each region and extract sequence with `bedtools getfasta`.
  6. Write FASTA headers as "<region_id>__<chrom>:<padded_start>-<padded_end>"
     so genome coordinates survive into FlashFry's `contig` column without
     needing --transformPositions.

Both input types are handled in one invocation so that linked-region
flagging sees every target, not just those of one type.
"""

import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

RANGE_RE = re.compile(r"^(\S+):([\d,]+)-([\d,]+)$")
POINT_RE = re.compile(r"^(\S+):([\d,]+)$")

METADATA_COLUMNS = [
    "region_id", "fasta_header", "chrom", "orig_start", "orig_end",
    "input_type", "target_ref", "padded_start", "padded_end", "linked_regions",
]


def _int(value: str) -> int:
    """Coordinates copy-pasted from genome browsers often carry thousands
    separators (chr1:1,000,000-1,000,500)."""
    return int(value.replace(",", ""))


def read_targets(bed_path, input_type, prefix=None):
    """Parses target regions, auto-detecting format per line:

      - tab-separated: BED (chrom, start, end[, name]), 0-based half-open
      - "chr:start-end" (enhancer files): same convention as BED start/end
      - "chr:pos" (SNP files), 1-based: becomes start=pos-1, end=pos

    Malformed lines are skipped with a warning rather than aborting the run.
    """
    prefix = prefix or input_type
    regions = []
    n_skipped = 0

    # A clear message beats a traceback. This fires for a deleted file, a
    # broken symlink, or -- the case that prompted it -- a bash process
    # substitution passed as a pipeline input: Nextflow stages inputs by
    # basename, so `<(...)` becomes a file literally named after the file
    # descriptor and nothing is staged.
    if not os.path.isfile(bed_path):
        raise SystemExit(
            f"ERROR: input file not found or not a regular file: {bed_path}\n"
            f"If you passed a process substitution such as <(...) or a pipe, "
            f"write the targets to a real file first -- they cannot be staged."
        )

    with open(bed_path) as fh:
        for i, raw in enumerate(fh):
            line = raw.rstrip("\n").strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue

            if "\t" in raw:
                fields = raw.rstrip("\n").split("\t")
                try:
                    chrom, start, end = fields[0], _int(fields[1]), _int(fields[2])
                except (ValueError, IndexError):
                    print(f"WARNING: skipping malformed BED line {i + 1} in "
                          f"{bed_path}: {line!r}", file=sys.stderr)
                    n_skipped += 1
                    continue
                name = fields[3] if len(fields) > 3 and fields[3] else f"{prefix}_{i + 1}"
            else:
                m_range = RANGE_RE.match(line)
                m_point = POINT_RE.match(line)
                if input_type == "enhancer" and m_range:
                    chrom = m_range.group(1)
                    start, end = _int(m_range.group(2)), _int(m_range.group(3))
                    name = line
                elif input_type == "snp" and m_point:
                    chrom, pos = m_point.group(1), _int(m_point.group(2))
                    start, end = pos - 1, pos
                    name = line
                elif input_type == "enhancer" and m_point:
                    print(f"WARNING: skipping line {i + 1} in {bed_path} -- looks like "
                          f"a SNP (chr:pos) in an enhancer file, expected "
                          f"chr:start-end: {line!r}", file=sys.stderr)
                    n_skipped += 1
                    continue
                elif input_type == "snp" and m_range:
                    print(f"WARNING: skipping line {i + 1} in {bed_path} -- looks like "
                          f"an interval (chr:start-end) in a SNP file, expected "
                          f"chr:pos: {line!r}", file=sys.stderr)
                    n_skipped += 1
                    continue
                else:
                    print(f"WARNING: skipping unparseable line {i + 1} in "
                          f"{bed_path}: {line!r}", file=sys.stderr)
                    n_skipped += 1
                    continue

            if end < start:
                print(f"WARNING: skipping line {i + 1} in {bed_path} -- end < start: "
                      f"{line!r}", file=sys.stderr)
                n_skipped += 1
                continue

            regions.append({
                "chrom": chrom, "start": start, "end": end,
                "name": name, "input_type": input_type,
            })

    if n_skipped:
        print(f"WARNING: skipped {n_skipped} malformed line(s) total in {bed_path}",
              file=sys.stderr)
    return regions


def load_chrom_sizes(fai_path):
    sizes = {}
    with open(fai_path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                sizes[fields[0]] = int(fields[1])
    return sizes


def validate_contigs(regions, chrom_sizes):
    """Drops regions whose contig isn't in the reference, with a diagnosis.

    The common case is Ensembl-style naming ("1") against a UCSC-style
    reference ("chr1"). bedtools would fail on this too, but with an error
    that sends people hunting in the wrong direction.
    """
    known = set(chrom_sizes)
    kept, dropped = [], []
    for region in regions:
        (kept if region["chrom"] in known else dropped).append(region)

    if dropped:
        missing = sorted({r["chrom"] for r in dropped})
        preview = ", ".join(missing[:5]) + (" ..." if len(missing) > 5 else "")
        print(f"WARNING: dropped {len(dropped)} region(s) on {len(missing)} contig(s) "
              f"absent from the reference: {preview}", file=sys.stderr)

        ucsc_style = any(c.startswith("chr") for c in known)
        plain_missing = [c for c in missing if not c.startswith("chr")]
        if ucsc_style and plain_missing:
            print("HINT: the reference uses UCSC-style contig names (chr1, chr2, ...) "
                  "but your input does not. Prefix your contigs with 'chr', or use a "
                  "reference matching your input's naming.", file=sys.stderr)
        elif not ucsc_style and any(c.startswith("chr") for c in missing):
            print("HINT: the reference uses Ensembl-style contig names (1, 2, ...) "
                  "but your input uses 'chr' prefixes. Strip the prefix, or use a "
                  "matching reference.", file=sys.stderr)

    if not kept:
        raise SystemExit(
            "ERROR: no input regions remain after contig validation -- every contig "
            "in your input is absent from the reference .fai. See the hint above."
        )
    return kept


def make_names_unique(regions):
    """Appends _dup1, _dup2, ... to repeated region names.

    Duplicates are routine in real input (the same enhancer listed twice, a
    BED whose name column isn't unique). Previously the last one silently
    won the FASTA-header lookup, and the metadata join then failed its
    many_to_one validation.
    """
    seen = defaultdict(int)
    n_renamed = 0
    for region in regions:
        name = region["name"]
        seen[name] += 1
        if seen[name] > 1:
            region["name"] = f"{name}_dup{seen[name] - 1}"
            n_renamed += 1
    if n_renamed:
        print(f"WARNING: {n_renamed} duplicate region name(s) made unique with a "
              f"_dupN suffix", file=sys.stderr)
    return regions


def flag_linked_regions(regions, link_distance_bp):
    """Marks same-contig regions within link_distance_bp of each other."""
    linked = {r["name"]: [] for r in regions}
    by_chrom = defaultdict(list)
    for region in regions:
        by_chrom[region["chrom"]].append(region)

    for regs in by_chrom.values():
        regs_sorted = sorted(regs, key=lambda r: r["start"])
        for a, b in zip(regs_sorted, regs_sorted[1:]):
            if b["start"] - a["end"] < link_distance_bp:
                linked[a["name"]].append(b["name"])
                linked[b["name"]].append(a["name"])
    return linked


def compute_window(region, snp_window_bp, flank_bp, chrom_sizes):
    if region["input_type"] == "snp":
        point = region["start"]
        padded_start = max(0, point - snp_window_bp)
        padded_end = point + snp_window_bp
        target_ref = point
    else:
        padded_start = max(0, region["start"] - flank_bp)
        padded_end = region["end"] + flank_bp
        target_ref = (region["start"] + region["end"]) // 2

    chrom_len = chrom_sizes.get(region["chrom"])
    if chrom_len is not None:
        padded_end = min(padded_end, chrom_len)
        padded_start = min(padded_start, max(0, chrom_len - 1))
    return padded_start, padded_end, target_ref


def run(args):
    if not args.snp_bed and not args.enhancer_bed:
        raise SystemExit("ERROR: at least one of --snp-bed / --enhancer-bed is required")
    if shutil.which("bedtools") is None:
        raise SystemExit("ERROR: 'bedtools' not found on PATH "
                         "(conda install -c bioconda bedtools)")

    fasta = args.fasta
    fai_path = fasta + ".fai"
    if not Path(fai_path).exists():
        raise SystemExit(
            f"ERROR: no FASTA index at {fai_path}\n"
            f"Build one with:  samtools faidx {fasta}"
        )
    chrom_sizes = load_chrom_sizes(fai_path)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    regions = []
    if args.snp_bed:
        regions += read_targets(args.snp_bed, "snp")
    if args.enhancer_bed:
        regions += read_targets(args.enhancer_bed, "enhancer")
    if not regions:
        raise SystemExit("ERROR: no usable regions parsed from the input file(s)")

    regions = validate_contigs(regions, chrom_sizes)
    regions = make_names_unique(regions)

    linked = (
        flag_linked_regions(regions, args.blast_radius_bp)
        if args.dedup_blast_radius
        else {r["name"]: [] for r in regions}
    )

    padded_bed = outdir / "regions.padded.bed"
    metadata_rows = []
    header_by_name = {}

    with open(padded_bed, "w") as bed_out:
        for region in regions:
            padded_start, padded_end, target_ref = compute_window(
                region, args.snp_window_bp, args.flank_bp, chrom_sizes)
            fasta_header = (f"{region['name']}__{region['chrom']}:"
                            f"{padded_start}-{padded_end}")
            header_by_name[region["name"]] = fasta_header
            bed_out.write(f"{region['chrom']}\t{padded_start}\t{padded_end}\t"
                          f"{region['name']}\n")
            metadata_rows.append({
                "region_id": region["name"],
                "fasta_header": fasta_header,
                "chrom": region["chrom"],
                "orig_start": region["start"],
                "orig_end": region["end"],
                "input_type": region["input_type"],
                "target_ref": target_ref,
                "padded_start": padded_start,
                "padded_end": padded_end,
                "linked_regions": ",".join(sorted(set(linked[region["name"]]))),
            })

    fasta_out = outdir / "regions.fasta"
    raw_fasta = outdir / "regions.raw.fasta"
    subprocess.run(
        ["bedtools", "getfasta", "-fi", fasta, "-bed", str(padded_bed),
         "-name", "-fo", str(raw_fasta)],
        check=True,
    )

    # bedtools >=2.27 writes "name::chrom:start-end" for -name; older builds
    # write just "name". Splitting on "::" handles both.
    with open(raw_fasta) as src, open(fasta_out, "w") as dst:
        for line in src:
            if line.startswith(">"):
                region_id = line[1:].strip().split("::")[0]
                dst.write(f">{header_by_name[region_id]}\n")
            else:
                dst.write(line)
    raw_fasta.unlink(missing_ok=True)

    metadata_path = outdir / "regions.metadata.tsv"
    with open(metadata_path, "w", newline="") as meta_out:
        writer = csv.DictWriter(meta_out, fieldnames=METADATA_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(metadata_rows)

    n_snp = sum(1 for r in regions if r["input_type"] == "snp")
    n_enh = len(regions) - n_snp
    print(f"Wrote {len(regions)} regions ({n_snp} SNP, {n_enh} enhancer) -> {fasta_out}",
          file=sys.stderr)
    return fasta_out


def add_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--snp-bed", dest="snp_bed")
    parser.add_argument("--enhancer-bed", dest="enhancer_bed")
    parser.add_argument("--fasta", help="reference FASTA (needs a .fai alongside)")
    parser.add_argument("--flank-bp", dest="flank_bp", type=int, default=50)
    parser.add_argument("--snp-window-bp", dest="snp_window_bp", type=int, default=150)
    parser.add_argument("--blast-radius-bp", dest="blast_radius_bp", type=int, default=1000)
    parser.add_argument("--dedup-blast-radius", dest="dedup_blast_radius",
                        action="store_true")
    parser.add_argument("--outdir", default=".")
