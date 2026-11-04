"""Pack sampled reads into ordered, non-overlapping display rows."""

from __future__ import annotations

import heapq
from dataclasses import dataclass

from ..models.read_plot import ReadAlignment


@dataclass(frozen=True)
class PlacedRead:
    read: ReadAlignment
    row: int


@dataclass(frozen=True)
class PairConnector:
    left: int
    right: int
    row: int


@dataclass(frozen=True)
class ReadLayout:
    reads: tuple[PlacedRead, ...]
    connectors: tuple[PairConnector, ...]
    row_count: int


def _groups(reads: tuple[ReadAlignment, ...], link_pairs: bool):
    if not link_pairs:
        return [((read,), read.start, read.end) for read in reads]
    grouped: dict[tuple[str, str], list[ReadAlignment]] = {}
    for read in reads:
        if read.is_paired:
            grouped.setdefault((read.name, read.alignment_class), []).append(read)
    used: set[int] = set()
    result = []
    for candidates in grouped.values():
        first = next((read for read in candidates if read.is_read1), None)
        second = next((read for read in candidates if not read.is_read1), None)
        if first is not None and second is not None and first is not second:
            result.append(((first, second), min(first.start, second.start), max(first.end, second.end)))
            used.update((id(first), id(second)))
    result.extend(((read,), read.start, read.end) for read in reads if id(read) not in used)
    return result


def pack_read_rows(reads: tuple[ReadAlignment, ...], link_pairs: bool) -> ReadLayout:
    """Place the earliest genomic start first, reusing rows after their end."""
    groups = sorted(_groups(reads, link_pairs), key=lambda item: (item[1], item[2], item[0][0].name))
    active: list[tuple[int, int]] = []
    free_rows: list[int] = []
    placed: list[PlacedRead] = []
    connectors: list[PairConnector] = []
    row_count = 0
    for members, start, end in groups:
        while active and active[0][0] + 2 <= start:
            _, free = heapq.heappop(active)
            heapq.heappush(free_rows, free)
        if free_rows:
            row = heapq.heappop(free_rows)
        else:
            row = row_count
            row_count += 1
        heapq.heappush(active, (end, row))
        placed.extend(PlacedRead(read, row) for read in members)
        if len(members) == 2:
            left, right = sorted(members, key=lambda read: read.start)
            if left.end < right.start:
                connectors.append(PairConnector(left.end, right.start, row))
    return ReadLayout(tuple(placed), tuple(connectors), row_count)
