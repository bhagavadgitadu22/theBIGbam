"""Indexed regional BAM reads, independent of plotting models."""

from __future__ import annotations

from time import perf_counter
from typing import Iterator

import pysam

from ..models.read_candidate import ReadCandidate
from .read_matching import exact_matching_bases

BAM_THREADS = 2


def _alignment_class(read: pysam.AlignedSegment) -> str:
    if read.is_supplementary:
        return "supplementary"
    if read.is_secondary:
        return "secondary"
    return "primary"


def _candidate(read: pysam.AlignedSegment, alignment_class: str, offset, region, include_matching) -> ReadCandidate:
    matching, length = exact_matching_bases(read) if include_matching else (None, None)
    return ReadCandidate(
        name=read.query_name or "",
        start=read.reference_start,
        end=read.reference_end,
        alignment_class=alignment_class,
        mapq=read.mapping_quality,
        cigar=read.cigarstring or "",
        mate_start=read.next_reference_start if read.is_paired and not read.mate_is_unmapped else None,
        flag=read.flag,
        virtual_offset=offset,
        fetch_region=region,
        template_length=read.template_length,
        mate_same_reference=read.reference_id == read.next_reference_id,
        matching_bases=matching,
        full_read_length=length,
    )


def _located_records(bam, contig, start, end):
    # A regional iterator can skip records that do not overlap the window.
    # Its previous result's end offset is therefore not the next record's
    # start. Seek with the index once, then capture offsets before raw reads.
    first = next(bam.fetch(contig, start, end), None)
    if first is None:
        return
    yield first, None
    reference_id = first.reference_id
    while True:
        before = bam.tell()
        read = next(bam, None)
        if read is None or read.reference_id != reference_id or read.reference_start >= end:
            return
        if read.reference_end is not None and read.reference_end > start:
            yield read, before


def iter_bam_records(bam_path, contig, start, end, profiler=None, *, threads=BAM_THREADS, with_offsets=False):
    """Keep iterator/I/O timings exclusive of checks, summaries and sampling."""
    if profiler is None:
        with pysam.AlignmentFile(bam_path, "rb", require_index=True, threads=threads) as bam:
            yield from _located_records(bam, contig, start, end) if with_offsets else bam.fetch(contig, start, end)
        return
    with profiler.phase("bam_open"):
        bam = pysam.AlignmentFile(bam_path, "rb", require_index=True, threads=threads)
    try:
        with profiler.phase("bam_region_setup"):
            records = iter(_located_records(bam, contig, start, end) if with_offsets else bam.fetch(contig, start, end))
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


def _region_candidates(
    bam_path: str, contig: str, start: int, end: int, classes: frozenset[str], *, profiler=None, include_matching=False,
) -> Iterator[ReadCandidate]:
    """Yield eligible reads overlapping a 0-based, half-open reference interval."""
    for read, offset in iter_bam_records(bam_path, contig, start, end, profiler, with_offsets=True):
        if profiler is None:
            alignment_class = _eligible_class(read, classes)
        else:
            profiler.count("bam_examined")
            with profiler.phase("bam_eligibility"):
                alignment_class = _eligible_class(read, classes)
        if alignment_class is None:
            continue
        if profiler is None:
            candidate = _candidate(read, alignment_class, offset, (start, end), include_matching)
        else:
            with profiler.phase("bam_candidate"):
                candidate = _candidate(read, alignment_class, offset, (start, end), include_matching)
            if include_matching:
                profiler.count("alignment_percentage_known" if candidate.matching_bases is not None else "alignment_percentage_unknown")
        # Equivalent to the old nonempty get_blocks() check, without building
        # per-operation blocks for alignments that will be discarded.
        if "M" in candidate.cigar or "=" in candidate.cigar or "X" in candidate.cigar:
            yield candidate


def fetch_read_candidates(bam_path, contig, start, end, classes, *, profiler=None, circular_length=None, include_matching=False):
    """Include the wrapped ends of circular alignments when querying the origin."""
    seen = set()
    for candidate in _region_candidates(bam_path, contig, start, end, classes, profiler=profiler, include_matching=include_matching):
        if circular_length and candidate.end > circular_length:
            seen.add((candidate.name, candidate.start, candidate.cigar, candidate.flag))
        yield candidate
    if circular_length and end < circular_length:
        for candidate in _region_candidates(bam_path, contig, circular_length - 1, circular_length, classes, profiler=profiler, include_matching=include_matching):
            key = (candidate.name, candidate.start, candidate.cigar, candidate.flag)
            if candidate.end - circular_length > start and key not in seen:
                seen.add(key)
                yield candidate
