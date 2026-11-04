"""Pack sampled reads into ordered, non-overlapping display rows."""

from __future__ import annotations

import heapq
from bisect import bisect_left, insort
from dataclasses import dataclass

from ..models.read_plot import ReadAlignment
from .read_layout_groups import packing_groups, pair_connector


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


def _linear_rows(groups):
    active: list[tuple[int, int]] = []
    free_rows: list[int] = []
    row_count = 0
    for group in groups:
        start, end = group.intervals[0]
        while active and active[0][0] + 2 <= start:
            _, free = heapq.heappop(active)
            heapq.heappush(free_rows, free)
        if free_rows:
            row = heapq.heappop(free_rows)
        else:
            row = row_count
            row_count += 1
        heapq.heappush(active, (end, row))
        yield group, row


def _fits(occupied, intervals):
    for start, end in intervals:
        index = bisect_left(occupied, (start,))
        if index and occupied[index - 1][1] + 2 > start:
            return False
        if index < len(occupied) and end + 2 > occupied[index][0]:
            return False
    return True


def _circular_rows(groups):
    rows = []
    for group in groups:
        row = next((index for index, occupied in enumerate(rows) if _fits(occupied, group.intervals)), len(rows))
        if row == len(rows):
            rows.append([])
        for interval in group.intervals:
            insort(rows[row], interval)
        yield group, row


def pack_read_rows(reads: tuple[ReadAlignment, ...], link_pairs: bool, contig_length: int = 0) -> ReadLayout:
    """Reserve every displayed interval; only explicitly linked mates may overlap."""
    groups = packing_groups(reads, link_pairs, contig_length)
    assignments = _linear_rows(groups) if all(len(group.intervals) == 1 for group in groups) else _circular_rows(groups)
    placed, connectors, row_count = [], [], 0
    for group, row in assignments:
        row_count = max(row_count, row + 1)
        placed.extend(PlacedRead(read, row) for read in group.reads)
        connector = pair_connector(group, contig_length)
        if connector is not None and connector[0] < connector[1]:
            connectors.append(PairConnector(*connector, row))
    return ReadLayout(tuple(placed), tuple(connectors), row_count)
