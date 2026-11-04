"""Prepare a bounded read plot from an indexed BAM."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from time import perf_counter

from ..models.read_plot import ReadPlotData, ReadPlotRequest
from ..repositories.bam_reads import fetch_read_candidates
from ..shared.data_cache import SessionDataCache
from ..shared.read_diagnostics import timed_bam_reads
from ..shared.read_identity import read_identity
from .read_discrepancies import add_read_discrepancies
from .read_event_coordinates import project_events
from .read_filtering import filtered_candidates
from .read_sampling import build_read_pool, sample_reads
from .read_summaries import summarize_candidate

REUSABLE_READ_CAPACITY = 10_000


def _fingerprint(path: str) -> tuple[int, int]:
    stat = Path(path).stat()
    return stat.st_size, stat.st_mtime_ns


def prepare_read_plot(request: ReadPlotRequest, cache: SessionDataCache | None = None) -> ReadPlotData:
    if not request.classes:
        return ReadPlotData((), 0, "Select a read alignment class in Variables to show reads.")
    if request.max_reads < 1:
        return ReadPlotData((), 0, "Number of subsampled reads must be at least 1.")
    if request.event_variables and abs(request.end - request.start) > request.max_discrepancy_window:
        if request.profiler is not None:
            request.profiler.count("discrepancies_skipped_window")
        request = replace(request, event_variables=frozenset())
    if request.members:
        if request.members == ((request.contig, request.contig_length, 0),):
            return prepare_read_plot(replace(request, members=()), cache)
        return _prepare_members(request, cache)
    start = max(0, request.start - 1)
    end = min(request.contig_length, request.end)
    if end <= start:
        return ReadPlotData((), 0, "Select a valid contig interval to show reads.")
    region_key = (request.bam_path, _fingerprint(request.bam_path), request.contig,
                  start, end, tuple(sorted(request.classes)), request.circular, request.filters)
    reference_identity = id(getattr(request.reference_source, "conn", request.reference_source))
    key = ("read_plot", *region_key, request.max_reads, tuple(sorted(request.event_variables)), reference_identity)
    if cache is not None:
        found, value = cache.get(key)
        if request.profiler is not None:
            request.profiler.count("cache_hits" if found else "cache_misses")
        if found:
            return value
    profiler = request.profiler
    pool = None
    pool_key = ("read_candidates", *region_key)
    if cache is not None:
        found, cached_pool = cache.get(pool_key)
        if found and (cached_pool.capacity >= request.max_reads or len(cached_pool.candidates) == cached_pool.fetched):
            pool = cached_pool
        if profiler is not None:
            profiler.count("candidate_cache_hits" if pool is not None else "candidate_cache_misses")
    started = perf_counter() if profiler is not None else 0
    scan_before = profiler.seconds.get("bam_scan_including_candidates", 0) if profiler is not None else 0
    filter_before = profiler.seconds.get("read_filtering", 0) if profiler is not None else 0
    try:
        if pool is None:
            source = fetch_read_candidates(
                request.bam_path, request.contig, start, end, request.classes,
                **({"profiler": profiler} if profiler is not None else {}),
                **({"circular_length": request.contig_length} if request.circular else {}),
                **({"include_matching": True} if any(rule.column == "Alignment percentage"
                    for section in request.filters.sections for rule in section.predicates) else {}),
            )
            if profiler is not None:
                source = timed_bam_reads(source, profiler)
            if request.filters.active:
                source = filtered_candidates(source, request.filters, request.contig_length, request.circular, profiler)
            pool = build_read_pool(
                source,
                start, end, max(REUSABLE_READ_CAPACITY, request.max_reads) if cache is not None else request.max_reads,
            )
            if cache is not None:
                size = pool.retained_bytes()
                while size > cache.max_retained_bytes and pool.capacity > request.max_reads:
                    capacity = max(request.max_reads, pool.capacity // 2)
                    smaller = build_read_pool(pool.candidates, start, end, capacity)
                    pool = replace(smaller, fetched=pool.fetched)
                    size = pool.retained_bytes()
                retained = cache.put(pool_key, pool, retained_bytes=size)
                if profiler is not None:
                    profiler.count("candidate_cache_stores" if retained else "candidate_cache_oversize")
        candidates, _ = sample_reads(pool.candidates, start, end, request.max_reads)
        count = pool.fetched
    except (OSError, ValueError) as error:
        if profiler is not None:
            profiler.count("bam_errors")
        return ReadPlotData((), 0, f"Could not read indexed BAM region: {error}")
    finally:
        if profiler is not None:
            scan_elapsed = profiler.seconds.get("bam_scan_including_candidates", 0) - scan_before
            filter_elapsed = profiler.seconds.get("read_filtering", 0) - filter_before
            profiler.add_time("sampling", max(0, perf_counter() - started - scan_elapsed - filter_elapsed))
    summary_started = perf_counter() if profiler is not None else 0
    reads = tuple(summarize_candidate(read) for read in candidates)
    if profiler is not None:
        profiler.add_time("bam_summary", perf_counter() - summary_started)
        profiler.count("summaries_built", len(reads))
        profiler.count("candidate_pool_size", len(pool.candidates))
        profiler.count("candidate_pool_capacity", pool.capacity)
    reads, warnings, details_loaded = add_read_discrepancies(reads, candidates, request, cache)
    empty_message = "No reads match the read filters in this interval." if request.filters.active else "No selected reads overlap this interval."
    result = ReadPlotData(reads, count, None if reads else empty_message, warnings, details_loaded)
    if cache is not None:
        if details_loaded:
            cache.put(key, result)
        else:
            cache.discard(pool_key)
    return result


def _prepare_members(request, cache):
    """Fetch visible member intervals, then pack reads in concatenated coordinates."""
    alignments, messages, warnings = [], [], []
    fetched, details_loaded = 0, True
    for name, length, offset in request.members:
        start, end = max(1, request.start - offset), min(length, request.end - offset)
        if end < start:
            continue
        data = prepare_read_plot(replace(
            request, members=(), contig=name, contig_length=length, start=start, end=end,
        ), cache)
        fetched += data.fetched
        details_loaded &= data.details_loaded
        warnings.extend(data.warnings)
        if data.message:
            messages.append(f"{name}: {data.message}")
        started = perf_counter() if request.profiler is not None else 0
        for read in data.alignments:
            # Split circular alignments before adding the member offset, so
            # a contig-end crossing never spills into the next MAG member.
            left = read.start
            while left < read.end:
                local = left % length
                width = min(read.end - left, length - local)
                alignments.append(replace(
                    read, start=local + offset, end=local + width + offset,
                    blocks=((local + offset, local + width + offset),),
                    mate_start=read.mate_start + offset if read.mate_start is not None else None,
                    events=project_events(read.events, left, left + width, length, offset, read_start=read.start),
                    alignment_id=(name, *(read.alignment_id or read_identity(read))),
                ))
                left += width
        if request.profiler is not None:
            request.profiler.add_time("member_coordinates", perf_counter() - started)
    started = perf_counter() if request.profiler is not None else 0
    reads, _ = sample_reads(alignments, request.start - 1, request.end, request.max_reads)
    if request.profiler is not None:
        request.profiler.add_time("member_resampling", perf_counter() - started)
    return ReadPlotData(reads, fetched, "; ".join(messages) if messages else None, tuple(dict.fromkeys(warnings)), details_loaded)
