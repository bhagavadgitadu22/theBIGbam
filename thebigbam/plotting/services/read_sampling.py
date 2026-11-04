"""Deterministic, spatially balanced sampling and reusable candidate reservoirs."""

import heapq
from collections import deque
from dataclasses import dataclass
from sys import getsizeof

from ..shared.read_identity import read_identity
from ..shared.read_identity import read_score as _hash_score


def _score(read):
    cached = getattr(read, "sampling_score", None)
    return _hash_score(read) if cached is None else cached


def _retain(heap, item, limit):
    if len(heap) < limit:
        heapq.heappush(heap, item)
    elif item[0] > heap[0][0]:
        heapq.heapreplace(heap, item)


def _heaps(reads, start, end, limit, bin_size):
    bin_count = max(1, (end - start + bin_size - 1) // bin_size)
    per_bin = max(1, (limit + bin_count - 1) // bin_count)
    bins, occupied, global_heap = {}, [], []
    count = 0
    for read in reads:
        count += 1
        item = (-_score(read), count, read)
        _retain(global_heap, item, limit)
        bucket = min(bin_count - 1, max(0, (read.start - start) // bin_size))
        if bucket not in bins:
            # Round-robin only visits the first `limit` occupied bins when
            # there are more bins than places. Later bins need no storage.
            if len(bins) >= limit:
                if bucket >= -occupied[0]:
                    continue
                del bins[-heapq.heapreplace(occupied, -bucket)]
            else:
                heapq.heappush(occupied, -bucket)
            bins[bucket] = []
        _retain(bins[bucket], item, per_bin)
    return bins, global_heap, count


def _ranked(heap):
    return sorted(heap, key=lambda item: (-item[0], read_identity(item[2])))


def sample_reads(reads, start, end, limit, *, bin_size=1000):
    """Keep round-robin selection, hashing each candidate once."""
    if limit < 1:
        return (), 0
    bins, global_heap, count = _heaps(reads, start, end, limit, bin_size)
    queues = [deque(_ranked(bins[index])) for index in sorted(bins)]
    chosen = {}
    while queues and len(chosen) < limit:
        remaining = []
        for queue in queues:
            if len(chosen) >= limit:
                break
            read = queue.popleft()[2]
            chosen[read_identity(read)] = read
            if queue:
                remaining.append(queue)
        queues = remaining
    for _, _, read in _ranked(global_heap):
        if len(chosen) >= limit:
            break
        chosen.setdefault(read_identity(read), read)
    return tuple(sorted(chosen.values(), key=lambda read: (read.start, read.end, read.name, read.is_read1))), count


@dataclass(frozen=True)
class ReadPool:
    candidates: tuple
    fetched: int
    capacity: int

    def retained_bytes(self):
        return getsizeof(self) + getsizeof(vars(self)) + getsizeof(self.candidates) + sum(
            getsizeof(read) + getsizeof(vars(read)) + sum(getsizeof(value) for value in vars(read).values())
            for read in self.candidates
        )


def build_read_pool(reads, start, end, capacity):
    """Retain sufficient candidates to reproduce any sample up to capacity."""
    bins, global_heap, count = _heaps(reads, start, end, capacity, 1000)
    retained = {serial: read for _, serial, read in global_heap}
    for heap in bins.values():
        retained.update((serial, read) for _, serial, read in heap)
    return ReadPool(tuple(read for _, read in sorted(retained.items())), count, capacity)
