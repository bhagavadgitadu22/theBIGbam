"""Prepare the visible portion of a MAG overview track."""

from __future__ import annotations

from bisect import bisect_left, bisect_right
from dataclasses import dataclass


@dataclass(frozen=True)
class VisibleMagTrackData:
    boundaries: tuple[float, ...]
    segments: tuple[tuple[str, int, int, int], ...]
    dot_xs: tuple[float, ...] = ()
    dot_colors: tuple[str, ...] = ()


def visible_mag_track_data(members, visible_start, visible_end, track_dots=None):
    members = tuple(members)
    if not members:
        return VisibleMagTrackData((), ())
    offsets = [int(member[2]) for member in members]
    ends = [int(offset + length) for _name, length, offset in members]
    first = bisect_right(ends, visible_start)
    last = bisect_left(offsets, visible_end, lo=first)
    visible = tuple(
        (str(name), int(length), int(offset), int(offset + length)) for name, length, offset in members[first:last]
    )
    boundary_values = []
    for _name, _length, start, end in visible:
        if visible_start <= start <= visible_end:
            boundary_values.append(start)
        if visible_start <= end <= visible_end:
            boundary_values.append(end)
    boundaries = tuple(dict.fromkeys(boundary_values))
    dot_xs, dot_colors = (), ()
    if track_dots:
        pairs = [
            (float(x), str(color))
            for x, color in zip(track_dots.get("xs", ()), track_dots.get("colors", ()))
            if visible_start <= x <= visible_end
        ]
        dot_xs = tuple(x for x, _color in pairs)
        dot_colors = tuple(color for _x, color in pairs)
    return VisibleMagTrackData(boundaries, visible, dot_xs, dot_colors)
