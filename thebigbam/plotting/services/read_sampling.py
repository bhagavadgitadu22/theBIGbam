"""Uniform deterministic sampling and reusable candidate reservoirs."""

import heapq
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


def _top_heap(reads, limit):
    """Keep the lowest stable identity hashes from the full candidate stream."""
    heap = []
    count = 0
    for read in reads:
        count += 1
        item = (-_score(read), count, read)
        _retain(heap, item, limit)
    return heap, count


def _ranked(heap):
    return sorted(heap, key=lambda item: (-item[0], read_identity(item[2])))


def sample_reads(reads, limit):
    """Select by stable hash across the full population, then order by position."""
    if limit < 1:
        return (), 0
    heap, count = _top_heap(reads, limit)
    selected = (read for _, _, read in _ranked(heap))
    return tuple(sorted(selected, key=lambda read: (read.start, read.end, read.name, read.is_read1))), count


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


def build_read_pool(reads, capacity):
    """Retain sufficient candidates to reproduce any sample up to capacity."""
    heap, count = _top_heap(reads, capacity)
    return ReadPool(tuple(read for _, _, read in _ranked(heap)), count, capacity)
