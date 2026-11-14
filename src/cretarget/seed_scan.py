"""
Promoter-proximal seed scanning.

THE QUESTION THIS ANSWERS
-------------------------
FlashFry asks "where else could this guide cut?" -- a mismatch-tolerant
alignment problem. For CRISPRi that is the wrong question: dCas9 does not cut,
it sits. A guide only needs a perfect PAM-proximal seed next to any NGG to
load dCas9 onto DNA, and if that happens over somebody's promoter the gene is
repressed whether or not the rest of the guide matches. So this asks the other
question: where does this guide's seed match perfectly, close to a
transcription start site?

Both screens are worth running. Neither substitutes for the other.

WHERE THIS COMES FROM
---------------------
The approach, and the pre-trimmed promoter reference it scans, are from
seed-finder-lite (Kazemi, M., Zenodo 10.5281/zenodo.21599135): screen seeds
against +/-2 kb windows around every annotated TSS rather than the whole
genome, which cuts the search space from 3.1 Gb to ~212 Mb. The single-pass
batch architecture -- scan the reference once, match all guides by dictionary
lookup -- is also from that tool's batch mode.

HOW IT SCALES
-------------
Cost is dominated by one pass over the promoter reference and is essentially
independent of guide count: the per-position work is a dictionary lookup
whether there are twenty guides or two million. Guides are indexed by their
seed, hits are streamed out as they are found, and the reference is read as a
gzip stream, so nothing here holds either the reference or the hit list in
memory.

A NOTE ON ORIENTATION, WHICH IS EASY TO GET WRONG
-------------------------------------------------
A guide's spacer sequence is identical to the protospacer, and the PAM sits
immediately 3' of the protospacer on the protospacer's own strand. So there
are two distinct cases, and they need *different* dictionary keys:

  * protospacer on the + strand: NGG lies at +20..+22 of the reference, and
    the guide's seed appears in the reference as-is.
  * protospacer on the - strand: the PAM reads NGG on the minus strand, which
    is CCN on the reference, and the guide's seed appears in the reference
    reverse-complemented.

Indexing guides under only one of those two keys makes the other case match
nothing. The tool this came from had exactly that defect -- every guide was
indexed under reverse_complement(seed) and both scan branches looked up that
one key -- so it silently reported only minus-strand hits and found about half
of all true sites while exiting 0. `strand` here names the strand the
PROTOSPACER lies on, matching FlashFry and creTarget's own `orientation`
column. test_seed_scan.py pins both cases.
"""

from __future__ import annotations

import gzip
from bisect import bisect_left

GUIDE_LENGTH = 20
PAM_LENGTH = 3
SITE_LENGTH = GUIDE_LENGTH + PAM_LENGTH

# How the promoter reference is built. The scan cannot see past the window
# half-width, so a --seed_finder_tss_window larger than this finds nothing
# extra; cretarget warns rather than silently capping.
PROMOTER_HALF_WIDTH = 2000
# Each merged window is padded so that a site straddling the boundary is not
# truncated: a protospacer plus its PAM is 23 bp.
PROMOTER_PAD = SITE_LENGTH

_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def open_text(path):
    """Promoter references are large enough that they are normally gzipped."""
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


# ---------------------------------------------------------------------------
# TSS lookup
# ---------------------------------------------------------------------------

def load_tss_table(path):
    """Read `chr<TAB>tss<TAB>gene_name`, tolerating a header row.

    The header is detected rather than assumed: skipping the first line
    unconditionally drops a real gene from a headerless file, and silently
    losing one TSS is the kind of error nothing downstream would reveal.
    """
    rows = []
    with open_text(path) as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            try:
                position = int(fields[1])
            except ValueError:
                continue        # header, or a malformed row
            rows.append((fields[0], position, fields[2] or "?"))
    return rows


class TssIndex:
    """Nearest-TSS lookup in O(log n) per query.

    A linear scan over every TSS for every hit costs ~1.5 ms here, which is
    four minutes of pure annotation on a 170k-hit run and scales with the
    product of the two.
    """

    __slots__ = ("_by_chrom",)

    def __init__(self, rows):
        grouped = {}
        for chrom, position, gene in rows:
            grouped.setdefault(chrom, []).append((position, gene))
        self._by_chrom = {}
        for chrom, entries in grouped.items():
            entries.sort(key=lambda e: e[0])
            self._by_chrom[chrom] = (
                [e[0] for e in entries],
                [e[1] for e in entries],
            )

    def nearest(self, chrom, position):
        """(gene_name, distance) for the closest TSS, or ("", None)."""
        found = self._by_chrom.get(chrom)
        if not found:
            return "", None
        positions, genes = found
        i = bisect_left(positions, position)
        best_gene, best_distance = "", None
        # The nearest TSS is immediately on one side of the insertion point.
        for candidate in (i - 1, i):
            if 0 <= candidate < len(positions):
                distance = abs(position - positions[candidate])
                if best_distance is None or distance < best_distance:
                    best_distance = distance
                    best_gene = genes[candidate]
        return best_gene, best_distance

    def __len__(self):
        return sum(len(p) for p, _ in self._by_chrom.values())


# ---------------------------------------------------------------------------
# Guide indexing
# ---------------------------------------------------------------------------

class InvalidGuide(ValueError):
    pass


def build_seed_index(guides, min_seed_length):
    """Map each guide's seed to the guides carrying it, per orientation.

    Returns {("+"|"-", seed): [(guide_id, spacer), ...]}.

    The "+" key is the seed as it appears in the reference when the
    protospacer is on the plus strand; the "-" key is its reverse complement,
    which is how it appears when the protospacer is on the minus strand.
    Indexing both is what makes the scan find sites in both orientations.
    Several guides can share a seed -- at length 8 there are only 65,536
    possibilities -- so the values are lists.
    """
    if not 1 <= min_seed_length <= GUIDE_LENGTH:
        raise ValueError(
            f"min_seed_length must be between 1 and {GUIDE_LENGTH}, "
            f"got {min_seed_length}")

    index = {}
    seen = set()
    for guide_id, spacer in guides:
        spacer = spacer.strip().upper()
        if len(spacer) != GUIDE_LENGTH or not set(spacer) <= set("ACGT"):
            raise InvalidGuide(
                f"guide {guide_id!r} must be {GUIDE_LENGTH} bp of ACGT, "
                f"got {spacer!r}")
        if guide_id in seen:
            raise InvalidGuide(f"duplicate guide_id: {guide_id!r}")
        seen.add(guide_id)

        seed = spacer[-min_seed_length:]
        index.setdefault(("+", seed), []).append((guide_id, spacer))
        index.setdefault(("-", reverse_complement(seed)), []).append(
            (guide_id, spacer))
    return index


def pam_proximal_match_length(protospacer, spacer):
    """Consecutive matching bases from the PAM-proximal end (index 19 down)."""
    n = 0
    for i in range(GUIDE_LENGTH - 1, -1, -1):
        if protospacer[i] != spacer[i]:
            break
        n += 1
    return n


# ---------------------------------------------------------------------------
# The scan
# ---------------------------------------------------------------------------

def scan_window(sequence, chrom, offset, seed_index, min_seed_length):
    """Yield hits in one promoter window.

    `offset` is the 0-based genomic coordinate of sequence[0]. Emitted
    coordinates are 1-based inclusive and cover the protospacer only; the PAM
    is reported separately because it lies on opposite sides of the
    protospacer in the two orientations.
    """
    sequence = sequence.upper()
    length = len(sequence)
    if length < SITE_LENGTH:
        return
    seed_start = GUIDE_LENGTH - min_seed_length

    for i in range(length - SITE_LENGTH + 1):
        # --- protospacer on the + strand: reference reads <20bp><NGG> ---
        if sequence[i + 21] == "G" and sequence[i + 22] == "G":
            candidates = seed_index.get(
                ("+", sequence[i + seed_start:i + GUIDE_LENGTH]))
            if candidates:
                protospacer = sequence[i:i + GUIDE_LENGTH]
                start = offset + i + 1
                for guide_id, spacer in candidates:
                    seed_length = pam_proximal_match_length(protospacer, spacer)
                    if seed_length < min_seed_length:
                        continue        # seed collision between guides
                    yield {
                        "guide_id": guide_id,
                        "chr": chrom,
                        "start": start,
                        "end": start + GUIDE_LENGTH - 1,
                        "strand": "+",
                        "seed_length": seed_length,
                        "pam": sequence[i + GUIDE_LENGTH:i + SITE_LENGTH],
                        "pam_start": start + GUIDE_LENGTH,
                        "pam_end": start + SITE_LENGTH - 1,
                        "protospacer": protospacer,
                        "guide_sequence": spacer,
                    }

        # --- protospacer on the - strand: reference reads <CCN><20bp> ---
        if sequence[i] == "C" and sequence[i + 1] == "C":
            candidates = seed_index.get(
                ("-", sequence[i + PAM_LENGTH:i + PAM_LENGTH + min_seed_length]))
            if candidates:
                protospacer = reverse_complement(
                    sequence[i + PAM_LENGTH:i + SITE_LENGTH])
                start = offset + i + PAM_LENGTH + 1
                for guide_id, spacer in candidates:
                    seed_length = pam_proximal_match_length(protospacer, spacer)
                    if seed_length < min_seed_length:
                        continue
                    yield {
                        "guide_id": guide_id,
                        "chr": chrom,
                        "start": start,
                        "end": start + GUIDE_LENGTH - 1,
                        "strand": "-",
                        "seed_length": seed_length,
                        "pam": reverse_complement(sequence[i:i + PAM_LENGTH]),
                        "pam_start": start - PAM_LENGTH,
                        "pam_end": start - 1,
                        "protospacer": protospacer,
                        "guide_sequence": spacer,
                    }


def parse_window_name(name):
    """`prom|chr16|189257|189267` -> ("chr16", 189257, 189267), 1-based."""
    fields = (name or "").split("|")
    if len(fields) < 4 or fields[0] != "prom":
        raise ValueError(f"unrecognized promoter window name: {name!r}")
    return fields[1], int(fields[2]), int(fields[3])


def iter_fasta(path):
    """Stream (name, sequence) from a FASTA.

    Sequential streaming on purpose: the scan visits every window once, so it
    needs no random access. That is what lets the reference stay gzipped with
    no index, no decompression to a cache directory and no pyfaidx dependency.
    """
    name = None
    chunks = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def iter_hits(promoter_fa, guides, tss_index, min_seed_length, tss_window_bp,
              on_window=None):
    """Scan the promoter reference once and yield annotated hits.

    Each hit gains `nearest_tss_gene` and `tss_distance`, and hits further
    than `tss_window_bp` from any TSS are dropped: the reference is built at a
    fixed half-width, so the window is applied here rather than baked in.
    """
    seed_index = build_seed_index(guides, min_seed_length)

    for window_index, (name, sequence) in enumerate(iter_fasta(promoter_fa)):
        try:
            chrom, start_1, _end_1 = parse_window_name(name)
        except ValueError:
            continue            # not a promoter window; skip rather than fail
        if on_window is not None:
            on_window(window_index)

        for hit in scan_window(sequence, chrom, start_1 - 1, seed_index,
                               min_seed_length):
            midpoint = (hit["start"] + hit["end"]) // 2
            gene, distance = tss_index.nearest(hit["chr"], midpoint)
            if distance is None or distance > tss_window_bp:
                continue
            hit["nearest_tss_gene"] = gene
            hit["tss_distance"] = distance
            yield hit


def ucsc_url(chrom, start, end, db="hg38"):
    prefixed = chrom if str(chrom).startswith("chr") else f"chr{chrom}"
    return (f"https://genome.ucsc.edu/cgi-bin/hgTracks"
            f"?db={db}&position={prefixed}:{start}-{end}")
