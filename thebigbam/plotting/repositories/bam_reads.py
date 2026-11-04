"""Indexed regional BAM reads, independent of plotting models."""

from __future__ import annotations

from time import perf_counter
from typing import Iterator

import pysam

from ..models.read_candidate import ReadCandidate

BAM_THREADS = 2


def _alignment_class(read: pysam.AlignedSegment) -> str:
    if read.is_supplementary:
        return "supplementary"
    if read.is_secondary:
        return "secondary"
    return "primary"


def _candidate(read: pysam.AlignedSegment, alignment_class: str) -> ReadCandidate:
    return ReadCandidate(
        name=read.query_name or "",
        start=read.reference_start,
        end=read.reference_end,
        alignment_class=alignment_class,
        mapq=read.mapping_quality,
        cigar=read.cigarstring or "",
        mate_start=read.next_reference_start if read.is_paired and not read.mate_is_unmapped else None,
        flag=read.flag,
    )


def iter_bam_records(bam_path, contig, start, end, profiler=None, *, threads=BAM_THREADS):
    """Keep iterator/I/O timings exclusive of checks, summaries and sampling."""
    if profiler is None:
        with pysam.AlignmentFile(bam_path, "rb", require_index=True, threads=threads) as bam:
            yield from bam.fetch(contig, start, end)
        return
    with profiler.phase("bam_open"):
        bam = pysam.AlignmentFile(bam_path, "rb", require_index=True, threads=threads)
    try:
        with profiler.phase("bam_region_setup"):
            records = iter(bam.fetch(contig, start, end))
        while True:
            started = perf_counter()
            try:
                read = next(records)
            except StopIteration:
                return
            finally:
                profiler.add_time("bam_iterator", perf_counter() - started)
            yield read
    finally:
        with profiler.phase("bam_close"):
            bam.close()


def _eligible_class(read, classes):
    if read.is_unmapped or read.reference_end is None:
        return None
    alignment_class = _alignment_class(read)
    return alignment_class if alignment_class in classes else None


def fetch_read_candidates(
    bam_path: str, contig: str, start: int, end: int, classes: frozenset[str], *, profiler=None,
) -> Iterator[ReadCandidate]:
    """Yield eligible reads overlapping a 0-based, half-open reference interval."""
    for read in iter_bam_records(bam_path, contig, start, end, profiler):
        if profiler is None:
            alignment_class = _eligible_class(read, classes)
        else:
            profiler.count("bam_examined")
            with profiler.phase("bam_eligibility"):
                alignment_class = _eligible_class(read, classes)
        if alignment_class is None:
            continue
        if profiler is None:
            candidate = _candidate(read, alignment_class)
        else:
            with profiler.phase("bam_candidate"):
                candidate = _candidate(read, alignment_class)
        # Equivalent to the old nonempty get_blocks() check, without building
        # per-operation blocks for alignments that will be discarded.
        if "M" in candidate.cigar or "=" in candidate.cigar or "X" in candidate.cigar:
            yield candidate
