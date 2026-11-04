"""Rows reserve all circular pieces, allowing overlap only for linked mates."""

from dataclasses import replace

from thebigbam.plotting.models.read_plot import ReadAlignment
from thebigbam.plotting.services.read_coordinates import wrap_interval
from thebigbam.plotting.services.read_layout import pack_read_rows


def _read(name, start, end, **kwargs):
    return replace(ReadAlignment(name, start, end, ((start, end),), True, "primary", 60,
                                  f"{end-start}M", None, False, False), **kwargs)


def _assert_clear(layout, length, linked=False):
    for index, first in enumerate(layout.reads):
        for second in layout.reads[index + 1:]:
            if first.row != second.row:
                continue
            a, b = first.read, second.read
            if a.alignment_id is not None and a.alignment_id == b.alignment_id:
                continue
            if linked and a.is_paired and b.is_paired and a.name == b.name and a.is_read1 != b.is_read1:
                continue
            assert all(hi <= lo2 or hi2 <= lo for lo, hi in wrap_interval(a.start, a.end, length)
                       for lo2, hi2 in wrap_interval(b.start, b.end, length))


def test_wrapped_read_reserves_both_ends_and_reuses_its_middle_gap():
    reads = (_read("circular", 95, 110), _read("origin", 5, 15),
             _read("middle", 30, 40), _read("tail", 92, 97))
    layout = pack_read_rows(reads, False, 100)
    rows = {item.read.name: item.row for item in layout.reads}
    assert rows["circular"] != rows["origin"]
    assert rows["circular"] != rows["tail"]
    assert rows["middle"] == rows["circular"]
    assert layout.row_count == 2
    _assert_clear(layout, 100)


def test_multiple_origin_crossings_reserve_the_entire_contig():
    layout = pack_read_rows((_read("long", 95, 210), _read("middle", 30, 40)), False, 100)
    assert layout.row_count == 2
    _assert_clear(layout, 100)


def test_only_linked_mates_are_allowed_to_overlap():
    reads = (_read("pair", 10, 40, is_paired=True, is_read1=True),
             _read("pair", 20, 35, is_paired=True), _read("unrelated", 25, 30))
    unlinked = pack_read_rows(reads, False, 100)
    assert unlinked.row_count == 3
    _assert_clear(unlinked, 100)
    linked = pack_read_rows(reads, True, 100)
    mate_rows = {placed.row for placed in linked.reads if placed.read.name == "pair"}
    assert len(mate_rows) == 1 and linked.row_count == 2
    assert linked.connectors == ()
    _assert_clear(linked, 100, linked=True)


def test_split_mag_alignment_pieces_stay_on_one_nonoverlapping_row():
    key = ("chr", "wrapped", 95, 110)
    reads = (_read("wrapped", 595, 600, alignment_id=key),
             _read("wrapped", 500, 510, alignment_id=key),
             _read("origin", 502, 507), _read("middle", 530, 550))
    layout = pack_read_rows(reads, False, 600)
    assert len({placed.row for placed in layout.reads if placed.read.name == "wrapped"}) == 1
    assert layout.row_count == 2
    _assert_clear(layout, 600)
