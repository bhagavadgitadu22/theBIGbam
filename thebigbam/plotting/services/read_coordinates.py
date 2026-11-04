"""The reference intervals actually drawn for a circular alignment."""


def wrap_interval(start, end, length):
    if length < 1:
        return ((start, end),) if start < end else ()
    intervals = []
    while start < end:
        local = start % length
        width = min(end - start, length - local)
        intervals.append((local, local + width))
        start += width
    return tuple(intervals)


def merge_intervals(intervals):
    result = []
    for start, end in sorted(intervals):
        if result and start <= result[-1][1]:
            result[-1] = (result[-1][0], max(end, result[-1][1]))
        else:
            result.append((start, end))
    return tuple(result)
