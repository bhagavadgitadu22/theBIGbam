"""Opt-in timings for preparing and drawing a read track on APPLY."""

import json
from contextlib import contextmanager
from dataclasses import dataclass, field
from time import perf_counter


@dataclass
class ReadProfiler:
    operation: str = "apply"
    seconds: dict[str, float] = field(default_factory=dict)
    counts: dict[str, int] = field(default_factory=dict)
    started: float = field(default_factory=perf_counter)

    def add_time(self, name, elapsed):
        self.seconds[name] = self.seconds.get(name, 0.0) + elapsed

    def count(self, name, amount=1):
        self.counts[name] = self.counts.get(name, 0) + amount

    @contextmanager
    def phase(self, name):
        started = perf_counter()
        try:
            yield
        finally:
            self.add_time(name, perf_counter() - started)

    def report(self, request, data):
        seconds = dict(self.seconds)
        if "bam_scan_including_candidates" in seconds:
            measured = sum(seconds.get(name, 0) for name in (
                "bam_open", "bam_region_setup", "bam_iterator", "bam_eligibility", "bam_candidate", "bam_close",
            ))
            # Generator bookkeeping, counters and timing overhead remain here.
            seconds["bam_scan_other"] = max(0, seconds["bam_scan_including_candidates"] - measured)
        print("[read-profile] " + json.dumps({
            "operation": self.operation, "sample": request.sample_name, "contig": request.contig,
            "start": request.start, "end": request.end, "limit": request.max_reads,
            "elapsed_seconds": perf_counter() - self.started,
            "seconds": seconds, "counts": self.counts,
            "eligible_in_region": data.fetched, "plotted_alignments": len(data.alignments),
        }, sort_keys=True), flush=True)


def timed_bam_reads(reads, profiler):
    """Measure generator work without counting time spent by its consumer."""
    iterator = iter(reads)
    while True:
        started = perf_counter()
        try:
            read = next(iterator)
        except StopIteration:
            return
        finally:
            profiler.add_time("bam_scan_including_candidates", perf_counter() - started)
        profiler.count("bam_eligible")
        yield read
