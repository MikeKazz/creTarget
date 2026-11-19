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
The method is from Hartman et al. 2026, "Systematic identification of
seed-driven off-target effects in Perturb-seq experiments"
(doi:10.64898/2026.03.27.714658): screen a guide's PAM-proximal seed for exact
matches next to an NGG within a window of an annotated TSS, and treat a match
as evidence of possible off-target transcriptional repression. Their pipeline
is at https://github.com/AustinHartman/perturb_seed.

seed-finder-lite (Kazemi, M., Zenodo 10.5281/zenodo.21599135) is a
command-line implementation of that method. The pre-trimmed promoter reference
scanned here -- +/-2 kb windows around every annotated TSS rather than the
whole genome, which cuts the search space from 3.1 Gb to ~212 Mb -- and the
single-pass batch architecture (scan the reference once, match all guides by
dictionary lookup) come from it.

The scanner below is creTarget's own code.

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

class StrandlessTssTable(SystemExit):
    """Raised for a reference built before strand was recorded."""


def load_tss_table(path):
    """Read `chr<TAB>tss<TAB>gene_name<TAB>strand`, tolerating a header row.

    The header is detected rather than assumed: skipping the first line
    unconditionally drops a real gene from a headerless file, and silently
    losing one TSS is the kind of error nothing downstream would reveal.

    A three-column table -- one built before 0.5.0 -- is a hard error rather
    than a default of "+". Defaulting would silently invert the distance sign
    for every minus-strand gene while exiting 0, which is the same shape of
    failure as the one-orientation seed index fixed in 0.3.0.
    """
    rows = []
    widths = set()
    with open_text(path) as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").rstrip("\r").split("\t")
            if len(fields) < 3:
                continue
            try:
                position = int(fields[1])
            except ValueError:
                widths.add(len(fields))     # header, or a malformed row
                continue
            widths.add(len(fields))
            if len(fields) < 4 or fields[3] not in ("+", "-"):
                raise StrandlessTssTable(
                    f"ERROR: {path} has no usable strand column.\n"
                    f"  Expected 4 columns (chr, tss, gene_name, strand); "
                    f"found {len(fields)}.\n\n"
                    f"References built before creTarget 0.5.0 did not record "
                    f"strand, so guide-to-TSS distances could only be "
                    f"unsigned. Rebuild it:\n"
                    f"  cretarget build-seed-reference --fasta /path/to/hg38.fa --force"
                )
            rows.append((fields[0], position, fields[2] or "?", fields[3]))
    return rows


class TssIndex:
    """Nearest-TSS lookup in O(log n) per query.

    A linear scan over every TSS for every hit costs ~1.5 ms here, which is
    four minutes of pure annotation on a 170k-hit run and scales with the
    product of the two.
    """

    __slots__ = ("_by_chrom",)

    def __init__(self, rows):
        """`rows` are (chrom, position, gene_name, strand) 4-tuples. The
        strand is required rather than defaulted, because a wrong strand
        silently inverts the sign of every distance for that gene."""
        grouped = {}
        for chrom, position, gene, strand in rows:
            grouped.setdefault(chrom, []).append((position, gene, strand))
        self._by_chrom = {}
        for chrom, entries in grouped.items():
            entries.sort(key=lambda e: e[0])
            self._by_chrom[chrom] = (
                [e[0] for e in entries],
                [e[1] for e in entries],
                [e[2] for e in entries],
            )

    @staticmethod
    def signed_distance(position, tss, strand):
        """Negative upstream of the TSS, positive downstream, relative to the
        gene's own strand -- the convention Hartman et al. use, so that these
        distances can be compared with their repression curve directly.

        For a minus-strand gene "downstream" runs toward lower coordinates, so
        the genomic difference is negated.
        """
        return (position - tss) if strand == "+" else (tss - position)

    def within(self, chrom, position, window):
        """Every TSS within `window` bp, as (gene, signed_distance, strand),
        nearest first.

        Returns ALL of them, not just the closest. Promoter windows routinely
        span several TSSs -- at a bidirectional promoter two genes sit on
        opposite strands a few hundred bp apart, which is roughly 10% of human
        genes -- and dCas9 loaded there can repress either. Reporting only the
        marginally closer one drops the other silently, and biases the screen
        against exactly the compact GC-rich promoters where seed matches
        concentrate.
        """
        found = self._by_chrom.get(chrom)
        if not found:
            return []
        positions, genes, strands = found
        i = bisect_left(positions, position)

        hits = []
        j = i - 1
        while j >= 0 and position - positions[j] <= window:
            hits.append(j)
            j -= 1
        j = i
        while j < len(positions) and positions[j] - position <= window:
            hits.append(j)
            j += 1

        out = [
            (genes[k], self.signed_distance(position, positions[k], strands[k]),
             strands[k])
            for k in hits
        ]
        # Nearest first, then by gene name so ties are deterministic.
        out.sort(key=lambda e: (abs(e[1]), e[0]))
        return out

    def nearest(self, chrom, position):
        """(gene_name, signed_distance) for the closest TSS, or ("", None).

        Signed, like `within` -- one convention for the class.
        """
        found = self._by_chrom.get(chrom)
        if not found:
            return "", None
        positions, genes, strands = found
        i = bisect_left(positions, position)
        best = None
        # The nearest TSS is immediately on one side of the insertion point.
        for candidate in (i - 1, i):
            if 0 <= candidate < len(positions):
                distance = abs(position - positions[candidate])
                if best is None or distance < best[0]:
                    best = (distance, genes[candidate],
                            self.signed_distance(position, positions[candidate],
                                                 strands[candidate]))
        if best is None:
            return "", None
        return best[1], best[2]

    def __len__(self):
        return sum(len(p) for p, _, _ in self._by_chrom.values())


# ---------------------------------------------------------------------------
# Guide indexing
# ---------------------------------------------------------------------------

class InvalidGuide(ValueError):
    pass


def build_seed_index(guides, min_seed_length):
    """Map each guide's seed to the guides carrying it, per orientation.

    `guides` yields (guide_id, spacer) or (guide_id, spacer, locus), where
    locus is (chrom, start, end) giving the guide's OWN target site. Supplying
    it lets the scan recognise a self-hit: a guide whose target region overlaps
    a promoter window finds itself, a perfect full-length match that is not an
    off-target at all. On a cCRE library a large share of hits are these --
    measured at 73% of all hits in a genome-wide run -- so counting them as
    off-target risk is badly misleading.

    Returns {("+"|"-", seed): [(guide_id, spacer, locus), ...]}.

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
    for entry in guides:
        if len(entry) == 3:
            guide_id, spacer, locus = entry
        else:
            (guide_id, spacer), locus = entry, None
        spacer = spacer.strip().upper()
        if len(spacer) != GUIDE_LENGTH or not set(spacer) <= set("ACGT"):
            raise InvalidGuide(
                f"guide {guide_id!r} must be {GUIDE_LENGTH} bp of ACGT, "
                f"got {spacer!r}")
        if guide_id in seen:
            raise InvalidGuide(f"duplicate guide_id: {guide_id!r}")
        seen.add(guide_id)

        # One tuple per guide, referenced from both orientation keys, so the
        # locus costs memory once rather than twice.
        record = (guide_id, spacer, locus)
        seed = spacer[-min_seed_length:]
        index.setdefault(("+", seed), []).append(record)
        index.setdefault(("-", reverse_complement(seed)), []).append(record)
    return index


def is_self_hit(locus, chrom, start, end):
    """Whether this hit IS the guide's own target site.

    The guide's recorded span covers protospacer plus PAM, and which side the
    PAM sits on depends on orientation -- so containment is tested rather than
    an exact start, which would catch forward guides and miss every reverse
    one.

    The two sides use different conventions and must be reconciled here:

        locus        0-based half-open [start, end)   (from guides.ranked.tsv)
        start/end    1-based inclusive [start, end]   (from the hit table)

    `start - 1` converts the hit to 0-based. Comparing the 1-based start
    directly against the 0-based locus start, as this used to, admitted one
    extra base on the left -- a position lying entirely OUTSIDE the guide's
    own span. A genuine off-target protospacer beginning one base before the
    guide's own site was therefore reclassified as a self-hit and dropped
    from seedfinder_hit_count. Such a site needs a shift-tolerant seed, i.e.
    a homopolymer or repeat, which is exactly the risky case.
    """
    if not locus:
        return False
    return locus[0] == chrom and (start - 1) >= locus[1] and end <= locus[2]


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
                end = start + GUIDE_LENGTH - 1
                for guide_id, spacer, locus in candidates:
                    seed_length = pam_proximal_match_length(protospacer, spacer)
                    if seed_length < min_seed_length:
                        continue        # seed collision between guides
                    yield {
                        "self_hit": "yes" if is_self_hit(locus, chrom, start, end) else "no",
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
                end = start + GUIDE_LENGTH - 1
                for guide_id, spacer, locus in candidates:
                    seed_length = pam_proximal_match_length(protospacer, spacer)
                    if seed_length < min_seed_length:
                        continue
                    yield {
                        "self_hit": "yes" if is_self_hit(locus, chrom, start, end) else "no",
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
            nearby = tss_index.within(hit["chr"], midpoint, tss_window_bp)
            if not nearby:
                continue
            # First entry is the closest; the rest are the other genes whose
            # promoters this same site sits in, which would otherwise be lost.
            hit["nearest_tss_gene"] = nearby[0][0]
            hit["tss_distance"] = nearby[0][1]
            hit["tss_strand"] = nearby[0][2]
            hit["tss_genes_in_window"] = ";".join(e[0] for e in nearby)
            hit["tss_distances_in_window"] = ";".join(str(e[1]) for e in nearby)
            yield hit


def ucsc_url(chrom, start, end, db="hg38"):
    prefixed = chrom if str(chrom).startswith("chr") else f"chr{chrom}"
    return (f"https://genome.ucsc.edu/cgi-bin/hgTracks"
            f"?db={db}&position={prefixed}:{start}-{end}")
