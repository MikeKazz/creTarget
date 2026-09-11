"""
flashfry_chunked.py

Splits a combined region FASTA into N chunks (by sequence count) and runs
FlashFry discover+score on each chunk independently via an executor
(local parallel processes, or a SLURM array) -- generalizing the pattern
already proven for launch-dCas9 scoring to FlashFry itself, so the whole
pipeline stays fast at hundreds-of-thousands-of-guides scale instead of
running FlashFry as one giant serial call.

Chunking is by FASTA record count (one region = one record), not by guide
count -- we don't know how many guides a region will yield until FlashFry
runs, so this is an approximation. A --n_workers of N roughly means "split
into N roughly-equal-sized groups of regions", not "exactly N guides per
chunk". Good enough in practice: FlashFry's own per-guide cost dominates,
and region sizes in a typical enhancer/SNP input are fairly uniform.
"""

import os
import sys
from pathlib import Path


def compute_slurm_chunking(n_regions: int, chunk_size: int, max_array_size: int) -> tuple:
    """Given a target chunk_size (regions per array task) and a hard cap on
    array size (max_array_size -- SLURM's MaxArraySize, cluster-configured,
    commonly 1000-10001), returns (n_chunks, effective_chunk_size).

    If the naive chunk count (ceil(n_regions / chunk_size)) would exceed
    max_array_size, chunk_size is increased just enough to bring the count
    under the cap -- this is a deliberate, visible adjustment (the caller
    should print a warning when effective_chunk_size != chunk_size), not a
    silent one.
    """
    if chunk_size < 1:
        raise ValueError("chunk_size must be >= 1")
    if max_array_size < 1:
        raise ValueError("max_array_size must be >= 1")

    def _ceil_div(a, b):
        return -(-a // b)

    n_chunks = max(1, min(_ceil_div(n_regions, chunk_size), n_regions))
    if n_chunks <= max_array_size:
        return n_chunks, chunk_size

    effective_chunk_size = max(1, _ceil_div(n_regions, max_array_size))
    n_chunks = max(1, min(_ceil_div(n_regions, effective_chunk_size), n_regions))
    return n_chunks, effective_chunk_size


def read_fasta_records(fasta_path):
    """Yields (header, sequence_lines_joined) tuples."""
    header = None
    seq_lines = []
    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_lines)
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            yield header, "".join(seq_lines)


def split_fasta(fasta_path, outdir, n_chunks):
    records = list(read_fasta_records(fasta_path))
    if not records:
        raise SystemExit(f"ERROR: no FASTA records found in {fasta_path}")
    n_chunks = max(1, min(n_chunks, len(records)))  # never more chunks than records

    Path(outdir).mkdir(parents=True, exist_ok=True)
    chunk_size = (len(records) + n_chunks - 1) // n_chunks
    chunk_paths = []
    for i in range(n_chunks):
        chunk_records = records[i * chunk_size:(i + 1) * chunk_size]
        if not chunk_records:
            continue
        chunk_path = os.path.join(outdir, f"regions_chunk_{i:04d}.fasta")
        with open(chunk_path, "w") as f:
            for header, seq in chunk_records:
                f.write(f">{header}\n{seq}\n")
        chunk_paths.append(chunk_path)

    print(f"  {len(records)} regions -> {len(chunk_paths)} chunk(s) (~{chunk_size} regions/chunk)", file=sys.stderr)
    return chunk_paths


def run_chunked_discover_score(
    fasta_path, flashfry_jar, flashfry_db, max_mismatch, scoring_metrics,
    executor, outdir, n_chunks,
):
    """Returns list of scored chunk output paths (all successful; raises if
    any chunk failed rather than silently proceeding on incomplete data)."""
    chunk_dir = os.path.join(outdir, "flashfry_chunks")
    fasta_chunks = split_fasta(fasta_path, chunk_dir, n_chunks)

    commands = []
    labels = []
    scored_paths = []
    for i, chunk_fasta in enumerate(fasta_chunks):
        discover_out = os.path.join(chunk_dir, f"discover_{i:04d}.out")
        scored_out = os.path.join(chunk_dir, f"scored_{i:04d}.out")
        scored_paths.append(scored_out)
        cmd = (
            f"java -Xmx14g -jar {flashfry_jar} discover "
            f"--database {flashfry_db} --fasta {chunk_fasta} "
            f"--output {discover_out} --maxMismatch {max_mismatch} && "
            f"java -Xmx7g -jar {flashfry_jar} score "
            f"--input {discover_out} --output {scored_out} "
            f"--scoringMetrics {scoring_metrics} --database {flashfry_db} "
            f"--maxMismatch {max_mismatch}"
        )
        commands.append(cmd)
        labels.append(f"flashfry_chunk_{i:04d}")

    print(f"  running FlashFry discover+score on {len(commands)} chunk(s)", file=sys.stderr)
    results = executor.run_batch(commands, logdir=os.path.join(chunk_dir, "logs"), labels=labels)

    failed = [r for r in results if r.returncode != 0]
    if failed:
        details = "\n".join(f"  - {r.label}: see {r.log_path}" for r in failed)
        raise SystemExit(
            f"ERROR: {len(failed)} of {len(results)} FlashFry chunk(s) failed:\n{details}\n"
            f"Not proceeding with incomplete/partial results."
        )

    missing = [p for p in scored_paths if not os.path.exists(p) or os.path.getsize(p) == 0]
    if missing:
        raise SystemExit(f"ERROR: expected scored output missing/empty for: {missing}")

    return scored_paths
