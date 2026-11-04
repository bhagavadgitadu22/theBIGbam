"""Enrich sampled alignments without retaining BAM records or full sequences."""

from dataclasses import replace
from time import perf_counter

from ..repositories.read_details import selected_read_details
from ..shared.read_identity import read_identity
from .read_event_extraction import extract_read_events
from .read_reference import ReadReference

MAX_READ_EVENTS = 100_000


def add_read_discrepancies(alignments, candidates, request, cache):
    if not request.event_variables or not candidates:
        return alignments, (), True
    profiler, before = request.profiler, perf_counter()
    reference = ReadReference(request, cache)
    events, remaining, unavailable, truncated = {}, MAX_READ_EVENTS, 0, False
    try:
        for read, sequence, md in selected_read_details(request.bam_path, request.contig, candidates, profiler=profiler):
            current, missing, limited = extract_read_events(read, sequence, md, reference, request, remaining)
            events[read_identity(read)] = current
            remaining -= len(current)
            unavailable += int(missing)
            truncated |= limited
            if remaining <= 0:
                truncated = True
                break
    except (OSError, ValueError, StopIteration) as error:
        return alignments, (f"Could not load read discrepancies: {error}",), False
    finally:
        if profiler is not None:
            profiler.add_time("read_discrepancies", perf_counter() - before)
            profiler.count("discrepancy_events", MAX_READ_EVENTS - remaining)
    warnings = []
    if unavailable:
        warnings.append(f"Discrepancy bases unavailable for {unavailable:,} alignments without sequence, MD tags or reference bases.")
    if truncated:
        warnings.append(f"Discrepancy display limited to {MAX_READ_EVENTS:,} events; zoom in and click APPLY for more detail.")
    return tuple(replace(read, events=events.get(read_identity(read), ())) for read in alignments), tuple(warnings), True
