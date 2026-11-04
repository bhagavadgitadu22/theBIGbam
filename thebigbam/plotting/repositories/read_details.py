"""Retrieve sequence and MD tags only for selected alignment candidates."""

from time import perf_counter

import pysam

from ..shared.read_identity import read_identity
from .bam_reads import BAM_THREADS, _alignment_class


def _record_at(bam, contig, candidate):
    offset = candidate.virtual_offset
    try:
        if offset is None:
            if candidate.fetch_region is None:
                raise ValueError("Alignment has no BAM file location.")
            return next(bam.fetch(contig, *candidate.fetch_region))
        current = bam.tell()
        if current <= offset and current >> 16 == offset >> 16:
            while bam.tell() < offset:
                next(bam)
        else:
            bam.seek(offset)
        return next(bam)
    except StopIteration as error:
        raise ValueError("Selected alignment is no longer available in the BAM.") from error


def selected_read_details(path, contig, candidates, *, profiler=None):
    with pysam.AlignmentFile(path, "rb", require_index=True, threads=BAM_THREADS) as bam:
        for candidate in sorted(candidates, key=lambda item: item.virtual_offset or -1):
            before = perf_counter() if profiler is not None else 0
            record = _record_at(bam, contig, candidate)
            identity = (record.query_name or "", record.reference_start, record.reference_end,
                        record.cigarstring or "", record.is_read1, _alignment_class(record))
            if identity != read_identity(candidate) or record.flag != candidate.flag:
                raise ValueError("Selected BAM record does not match the sampled alignment.")
            sequence = record.query_sequence or ""
            md = str(record.get_tag("MD")) if record.has_tag("MD") else None
            if profiler is not None:
                profiler.add_time("read_detail_fetch", perf_counter() - before)
                profiler.count("read_details_loaded")
                profiler.count("read_sequence_bases", len(sequence))
            yield candidate, sequence, md
