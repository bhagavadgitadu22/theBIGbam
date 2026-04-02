"""Convert a circular-mapped BAM (doubled reference) to SAM-spec circular BAM.

After ``mapping-per-sample --circular``, every contig in the BAM has twice
its real length because the reference was concatenated to itself.  This
module rewrites the BAM so that:

* ``@SQ`` headers carry the **original** (halved) LN.
* Read POS is kept between 1 and LN.  Reads whose alignment extends past
  LN are left as-is — this is valid per the SAM specification (§1.4):

      *"POS plus the sum of the lengths of M/=/X/D/N CIGAR operations
      may exceed LN.  Coordinates greater than LN are interpreted by
      subtracting LN so that bases at LN+1, LN+2, LN+3, … are
      considered to be mapped at positions 1, 2, 3, …; thus each
      (1-based) position p is interpreted as ((p − 1) mod LN) + 1."*

* Ghost secondary and supplementary alignments produced by the doubled
  reference are filtered out.

.. warning::

   Although the output file is SAM-spec-compliant, **most downstream
   tools do not yet implement circular coordinate handling**.  Programs
   such as IGV, samtools depth, or variant callers will likely produce
   incorrect results for reads that wrap around the origin.  Use with
   caution and verify that your downstream tool explicitly supports
   circular genomes before relying on the output.
"""

from itertools import groupby
from pathlib import Path

import pysam


# ---------------------------------------------------------------------------
# Header helpers
# ---------------------------------------------------------------------------

_CIRCULAR_CO = "theBIGbam:circular=true"
_WARNING_CO = (
    "WARNING: This BAM follows the SAM specification for circular genomes "
    "(section 1.4): POS + alignment length may exceed SQ:LN, with coordinates "
    ">LN wrapping via ((p-1) mod LN)+1. While spec-compliant, most tools do "
    "not yet implement circular coordinate handling and may produce incorrect "
    "results. Use with caution."
)


def _build_converted_header(header: pysam.AlignmentHeader) -> dict:
    """Return a header dict with halved SQ lengths and updated CO lines."""
    hd = header.to_dict()

    # Halve every SQ length
    for sq in hd.get("SQ", []):
        sq["LN"] = sq["LN"] // 2

    # Keep existing CO lines as-is (circular=true tag stays) and add warning
    new_co = list(hd.get("CO", []))
    new_co.append(_WARNING_CO)
    hd["CO"] = new_co

    return hd


# ---------------------------------------------------------------------------
# Ghost alignment filtering
# ---------------------------------------------------------------------------

def _contained_in_circular(inner_start, inner_end, outer_start, outer_end, orig_ln):
    """Check if [inner_start, inner_end) is fully contained in [outer_start, outer_end) on a circular reference."""
    inner_len = inner_end - inner_start
    outer_len = outer_end - outer_start
    if inner_len <= 0 or outer_len <= 0 or inner_len > orig_ln:
        return False
    if outer_len >= orig_ln:
        return True
    offset = (inner_start - outer_start) % orig_ln
    return offset + inner_len <= outer_len


def _filter_ghost_alignments(reads, original_lengths):
    """Filter ghost secondary/supplementary alignments from a read name group.

    A secondary or supplementary is a ghost if its aligned interval is fully
    contained within any primary's aligned interval on the same contig
    (in circular coordinates).  Ghost supplementaries are also deduplicated
    by ``pos % original_ln``.

    Returns (kept_reads, n_ghost_secondary, n_ghost_supplementary,
    mapq_fix_reads) where mapq_fix_reads is a set of read objects whose
    ghosts were removed and have no remaining real secondaries/supplementaries.
    """
    primaries = []
    secondaries = []
    supplementaries = []

    for read in reads:
        if read.is_unmapped:
            primaries.append(read)  # unmapped → keep unconditionally
        elif read.is_secondary:
            secondaries.append(read)
        elif read.is_supplementary:
            supplementaries.append(read)
        else:
            primaries.append(read)

    # Build list of (tid, ref_start, ref_end) for all primary alignments
    primary_intervals = []
    for r in primaries:
        if not r.is_unmapped:
            tid = r.reference_id
            if 0 <= tid < len(original_lengths):
                primary_intervals.append((tid, r.reference_start, r.reference_end))

    # Filter ghost secondaries (fully contained in a primary's circular interval)
    kept_sec = []
    n_ghost_sec = 0
    for r in secondaries:
        tid = r.reference_id
        if 0 <= tid < len(original_lengths):
            orig_ln = original_lengths[tid]
            is_ghost = any(
                p_tid == tid and _contained_in_circular(
                    r.reference_start, r.reference_end, p_start, p_end, orig_ln)
                for p_tid, p_start, p_end in primary_intervals
            )
            if is_ghost:
                n_ghost_sec += 1
                continue
        kept_sec.append(r)

    # Filter ghost supplementaries (fully contained in primary, or duplicate pos_mod)
    kept_supp = []
    n_ghost_supp = 0
    seen_supp_keys = set()
    for r in supplementaries:
        tid = r.reference_id
        if 0 <= tid < len(original_lengths):
            orig_ln = original_lengths[tid]
            if any(
                p_tid == tid and _contained_in_circular(
                    r.reference_start, r.reference_end, p_start, p_end, orig_ln)
                for p_tid, p_start, p_end in primary_intervals
            ):
                n_ghost_supp += 1
                continue
            pos_mod = r.reference_start % orig_ln
            key = (tid, pos_mod)
            if key in seen_supp_keys:
                n_ghost_supp += 1
                continue
            seen_supp_keys.add(key)
        kept_supp.append(r)

    # Identify primaries needing MAPQ fix: ghosts were removed and
    # no real secondaries/supplementaries remain for this read group
    mapq_fix_reads = set()
    if (n_ghost_sec > 0 or n_ghost_supp > 0) and not kept_sec and not kept_supp:
        for r in primaries:
            if not r.is_unmapped:
                mapq_fix_reads.add(id(r))

    return primaries + kept_sec + kept_supp, n_ghost_sec, n_ghost_supp, mapq_fix_reads


# ---------------------------------------------------------------------------
# Read coordinate adjustment
# ---------------------------------------------------------------------------

def _adjust_read(read, original_lengths):
    """Shift POS / MPOS so they fall within [0, original_ln).

    CIGAR is left unchanged — the alignment simply starts at the
    normalised position and may extend past LN, which the SAM spec
    allows for circular references.
    """
    if read.is_unmapped:
        return

    tid = read.reference_id
    if tid < 0 or tid >= len(original_lengths):
        return
    orig_ln = original_lengths[tid]

    # --- primary position ---
    if read.reference_start >= orig_ln:
        read.reference_start -= orig_ln

    # --- mate position ---
    if read.is_paired and not read.mate_is_unmapped:
        mtid = read.next_reference_id
        if 0 <= mtid < len(original_lengths):
            mate_orig_ln = original_lengths[mtid]
            if read.next_reference_start >= mate_orig_ln:
                read.next_reference_start -= mate_orig_ln

        # Recompute TLEN for same-reference pairs using circular-aware shortest distance
        if read.reference_id == read.next_reference_id:
            direct = abs(read.reference_start - read.next_reference_start)
            wrapped = orig_ln - direct
            shortest = min(direct, wrapped)
            # Keep sign convention: positive for leftmost read
            if read.reference_start <= read.next_reference_start:
                read.template_length = shortest
            else:
                read.template_length = -shortest


# ---------------------------------------------------------------------------
# Main conversion
# ---------------------------------------------------------------------------

def convert_circular_bam(
    input_bam: Path,
    output_bam: Path,
    threads: int = 4,
) -> None:
    """Convert a doubled-reference circular BAM to SAM-spec circular BAM.

    This function:
    1. Name-sorts the input BAM
    2. Filters ghost secondary/supplementary alignments (duplicates from
       the doubled reference)
    3. Shifts coordinates so POS is in [0, original_LN)
    4. Position-sorts and indexes the output
    """
    input_bam = Path(input_bam)
    output_bam = Path(output_bam)

    # Step 1: Name-sort the input for grouping by read name
    namesorted = str(output_bam) + ".namesort.tmp.bam"
    pysam.sort("-n", "-@", str(threads), "-o", namesorted, str(input_bam))

    total_ghost_sec = 0
    total_ghost_supp = 0
    total_mapq_fixed = 0

    try:
        with pysam.AlignmentFile(namesorted, "rb", threads=threads) as infile:
            header = infile.header

            # Original (halved) lengths, indexed by tid
            original_lengths = [sq["LN"] // 2 for sq in header.to_dict()["SQ"]]
            new_header = _build_converted_header(header)

            # Step 2+3: Filter ghosts and adjust coordinates, write unsorted output
            unsorted_out = str(output_bam) + ".unsorted.tmp.bam"
            with pysam.AlignmentFile(
                unsorted_out, "wb", header=new_header, threads=threads,
            ) as outfile:
                for _qname, group in groupby(
                    infile.fetch(until_eof=True),
                    key=lambda r: r.query_name,
                ):
                    reads = list(group)
                    kept, n_gs, n_gsup, mapq_fix = _filter_ghost_alignments(reads, original_lengths)
                    total_ghost_sec += n_gs
                    total_ghost_supp += n_gsup
                    total_mapq_fixed += len(mapq_fix)

                    for read in kept:
                        _adjust_read(read, original_lengths)
                        if id(read) in mapq_fix:
                            read.mapping_quality = 60
                        outfile.write(read)

    finally:
        # Clean up name-sorted temp
        try:
            Path(namesorted).unlink()
        except OSError:
            pass

    # Step 4: Position-sort and index
    pysam.sort("-@", str(threads), "-o", str(output_bam), unsorted_out)
    try:
        Path(unsorted_out).unlink()
    except Exception:
        pass
    pysam.index("-@", str(threads), str(output_bam))

    print(
        f"Circular BAM conversion: discarded {total_ghost_sec} ghost secondary "
        f"and {total_ghost_supp} ghost supplementary alignments, "
        f"fixed MAPQ to 60 for {total_mapq_fixed} primary alignments",
        flush=True,
    )
