"""Group circular pieces and optionally linked mates for row placement."""

from collections import defaultdict
from dataclasses import dataclass

from .read_coordinates import merge_intervals, wrap_interval


@dataclass(frozen=True)
class PackingGroup:
    reads: tuple
    intervals: tuple[tuple[int, int], ...]
    linked: bool = False


def read_intervals(read, length):
    return wrap_interval(read.start, read.end, length)


def packing_groups(reads, link_pairs, length):
    pieces = defaultdict(list)
    for read in reads:
        pieces[read.alignment_id if read.alignment_id is not None else id(read)].append(read)
    alignments = list(pieces.values())
    grouped, used = [], set()
    if link_pairs:
        pairs = defaultdict(list)
        for index, alignment in enumerate(alignments):
            first = alignment[0]
            if first.is_paired:
                pairs[(first.name, first.alignment_class)].append(index)
        for indices in pairs.values():
            first = next((i for i in indices if alignments[i][0].is_read1), None)
            second = next((i for i in indices if not alignments[i][0].is_read1), None)
            if first is not None and second is not None:
                grouped.append(((*alignments[first], *alignments[second]), True))
                used.update((first, second))
    grouped.extend((tuple(alignment), False) for index, alignment in enumerate(alignments) if index not in used)
    result = []
    for members, linked in grouped:
        intervals = merge_intervals(span for read in members for span in read_intervals(read, length))
        if not intervals:
            continue
        if linked:
            # Reserve the connector as well; overlapping mates are deliberately
            # allowed inside this group, but unrelated reads must stay clear.
            intervals = ((intervals[0][0], intervals[-1][1]),)
        result.append(PackingGroup(tuple(members), intervals, linked))
    return sorted(result, key=lambda group: (group.intervals[0][0], min(read.start for read in group.reads),
                                            group.intervals[-1][1], group.reads[0].name))


def pair_connector(group, length):
    if not group.linked:
        return None
    first = merge_intervals(span for read in group.reads if read.is_read1 for span in read_intervals(read, length))
    second = merge_intervals(span for read in group.reads if not read.is_read1 for span in read_intervals(read, length))
    gaps = []
    for a, b in first:
        for c, d in second:
            if a < d and c < b:
                return None  # Overlapping linked mates need no connecting gap.
            gaps.append((b, c) if b <= c else (d, a))
    return min(gaps, key=lambda gap: gap[1] - gap[0]) if gaps else None
