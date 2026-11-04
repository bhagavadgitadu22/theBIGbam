"""Split read events at a circular origin and place them in MAG coordinates."""

from dataclasses import replace


def project_events(events, left, right, length, offset=0, *, read_start=None):
    result = []
    for event in events:
        if event.start == event.end:
            if event.kind in {"right_clippings", "insertions"}:
                belongs = left < event.start <= right or (event.kind == "insertions" and event.start == left == read_start)
            else:
                belongs = left <= event.start < right
            if belongs:
                position = left % length + event.start - left + offset
                result.append(replace(event, start=position, end=position))
            continue
        lo, hi = max(left, event.start), min(right, event.end)
        if lo < hi:
            position = left % length + lo - left + offset
            a, b = lo - event.start, hi - event.start
            result.append(replace(event, start=position, end=position + hi - lo,
                                  sequence=event.sequence[a:b], reference=event.reference[a:b]))
    return tuple(result)


def circular_events(read, length):
    if length < 1:
        return read.events
    result, left = [], read.start
    while left < read.end:
        right = min(read.end, left + length - left % length)
        result.extend(project_events(read.events, left, right, length, read_start=read.start))
        left = right
    return tuple(result)
