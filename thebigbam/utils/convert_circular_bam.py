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
_SPEC_CO = "theBIGbam:circular=spec"
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

    # Update CO lines
    new_co = []
    for co in hd.get("CO", []):
        if _CIRCULAR_CO in co:
            new_co.append(co.replace(_CIRCULAR_CO, _SPEC_CO))
        else:
            new_co.append(co)
    new_co.append(_WARNING_CO)
    hd["CO"] = new_co

    return hd


# ---------------------------------------------------------------------------
# Ghost alignment filtering
# ---------------------------------------------------------------------------

def _filter_ghost_alignments(reads, original_lengths):
    """Filter ghost secondary/supplementary alignments from a read name group.

    Ghost secondaries arise because the mapper maps the same read to both
    copies of the doubled reference.  A secondary is a ghost if its
    ``pos % original_ln`` matches any primary's position on the same contig.

    Ghost supplementaries are rarer: a non-junction supplementary could
    appear twice (once per copy).  We discard a supplementary only if
    another supplementary of the same query name has the same
    ``pos % original_ln`` on the same contig.

    Returns (kept_reads, n_ghost_secondary, n_ghost_supplementary).
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

    # Build set of (tid, pos_mod) for all primary alignments
    primary_positions = set()
    for r in primaries:
        if not r.is_unmapped:
            tid = r.reference_id
            if 0 <= tid < len(original_lengths):
                primary_positions.add((tid, r.reference_start % original_lengths[tid]))

    # Filter ghost secondaries
    kept_sec = []
    n_ghost_sec = 0
    for r in secondaries:
        tid = r.reference_id
        if 0 <= tid < len(original_lengths):
            pos_mod = r.reference_start % original_lengths[tid]
            if (tid, pos_mod) in primary_positions:
                n_ghost_sec += 1
                continue
        kept_sec.append(r)

    # Filter ghost supplementaries (duplicate supp at same pos_mod)
    kept_supp = []
    n_ghost_supp = 0
    seen_supp_positions = set()
    for r in supplementaries:
        tid = r.reference_id
        if 0 <= tid < len(original_lengths):
            pos_mod = r.reference_start % original_lengths[tid]
            key = (tid, pos_mod)
            if key in seen_supp_positions:
                n_ghost_supp += 1
                continue
            seen_supp_positions.add(key)
        kept_supp.append(r)

    return primaries + kept_sec + kept_supp, n_ghost_sec, n_ghost_supp


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
                    kept, n_gs, n_gsup = _filter_ghost_alignments(reads, original_lengths)
                    total_ghost_sec += n_gs
                    total_ghost_supp += n_gsup

                    for read in kept:
                        _adjust_read(read, original_lengths)
                        outfile.write(read)

    finally:
        # Clean up name-sorted temp
        try:
            Path(namesorted).unlink()
        except Exception:
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
        f"and {total_ghost_supp} ghost supplementary alignments",
        flush=True,
    )
