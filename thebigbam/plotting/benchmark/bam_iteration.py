"""Compare raw regional BAM iteration across paths and decompression threads."""

import argparse
import json
from pathlib import Path
from statistics import median
from time import perf_counter

import pysam

from ..repositories.bam_reads import iter_bam_records


def benchmark_iteration(paths, *, threads=(1, 2, 4), repeats=3, contig=None, start=1, end=None):
    """Alternate thread-count order; never use the application's read cache."""
    if repeats < 1 or start < 1 or not threads or min(threads) < 1:
        raise ValueError("Repeats, start and thread counts must be positive.")
    for path in paths:
        path = Path(path).resolve()
        with pysam.AlignmentFile(str(path), "rb", require_index=True) as bam:
            reference = contig or bam.references[0]
            stop = min(end if end is not None else bam.get_reference_length(reference), bam.get_reference_length(reference))
        if stop < start:
            raise ValueError("End must be at least start.")
        timings = {count: [] for count in threads}
        expected = None
        for trial in range(repeats):
            for count in threads if trial % 2 == 0 else tuple(reversed(threads)):
                before = perf_counter()
                records = sum(1 for _ in iter_bam_records(str(path), reference, start - 1, stop, threads=count))
                timings[count].append(perf_counter() - before)
                if expected is not None and records != expected:
                    raise ValueError("Record count changed between benchmark passes.")
                expected = records
        yield {
            "path": str(path), "bam_bytes": path.stat().st_size, "contig": reference,
            "start": start, "end": stop, "records": expected,
            "application_cache": False, "filesystem_cache": "uncontrolled",
            "threads": {count: {"seconds": values, "median_seconds": median(values)}
                        for count, values in timings.items()},
        }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bam_paths", nargs="+", help="Indexed BAM paths, including copies on different storage.")
    parser.add_argument("--threads", nargs="+", type=int, default=[1, 2, 4])
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--contig")
    parser.add_argument("--start", type=int, default=1)
    parser.add_argument("--end", type=int)
    args = parser.parse_args()
    for result in benchmark_iteration(args.bam_paths, threads=args.threads, repeats=args.repeats,
                                      contig=args.contig, start=args.start, end=args.end):
        print(json.dumps(result, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
