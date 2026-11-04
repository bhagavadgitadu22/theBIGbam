"""Exact matching percentages include clips and reject unknown match counts."""

import pysam
import pytest

from thebigbam.plotting.models.filters import FilterExpression, FilterPredicate, FilterSection
from thebigbam.plotting.models.read_candidate import ReadCandidate
from thebigbam.plotting.models.read_plot import ReadPlotRequest
from thebigbam.plotting.repositories.read_matching import exact_matching_bases
from thebigbam.plotting.services.read_filtering import compile_read_filters, validate_read_filters
from thebigbam.plotting.services.read_plots import prepare_read_plot


@pytest.mark.parametrize("cigar,nm,md,expected", [
    ("5H10S75=5X5I", None, None, (75, 100)),
    ("10S85M5I", 8, None, (82, 100)),
    ("5S80M5=5X5I", 13, None, (82, 100)),
    ("2S8M2I2D3M5H", None, "4A3^GT3", (10, 20)),
    ("2S8M2I2D3M5H", 5, None, (10, 20)),
    ("8M", None, None, (None, 8)),
    ("8M", "bad", None, (None, 8)),
    ("8M", -1, None, (None, 8)),
    ("8M", None, "10", (None, 8)),
    ("8M", 100, "4A3", (7, 8)),
])
def test_exact_matching_counts(cigar, nm, md, expected):
    read = pysam.AlignedSegment()
    read.reference_start = 0
    read.cigarstring = cigar
    if nm is not None:
        read.set_tag("NM", nm)
    if md is not None:
        read.set_tag("MD", md)
    assert exact_matching_bases(read) == expected


def _expression(operator, value):
    return FilterExpression((FilterSection((FilterPredicate("Reads", "Alignment percentage", operator, value),)),))


def test_percentage_filter_uses_full_read_length_and_unknowns_do_not_match():
    read = ReadCandidate("read", 0, 90, "10S90M", "primary", 0, 60, None, matching_bases=89)
    assert read.alignment_percentage == 89
    assert compile_read_filters(_expression(">", 88.5), 1000, False)(read)
    assert compile_read_filters(_expression("=", 89), 1000, False)(read)
    unknown = ReadCandidate("unknown", 0, 90, "10S90M", "primary", 0, 60, None)
    assert unknown.alignment_percentage is None
    assert not compile_read_filters(_expression("!=", 89), 1000, False)(unknown)
    with pytest.raises(ValueError):
        validate_read_filters(_expression(">", 101))


def test_percentage_is_calculated_before_subsampling(tmp_path):
    path = tmp_path / "matching.bam"
    header = {"HD": {"SO": "coordinate"}, "SQ": [{"SN": "chr", "LN": 1000}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for index, edits in enumerate((0, 1, 20, None)):
            read = pysam.AlignedSegment()
            read.query_name, read.query_sequence = f"r{index}", "A" * 100
            read.reference_id, read.reference_start = 0, index * 10
            read.cigarstring = "10S90M"
            if edits is not None:
                read.set_tag("NM", edits)
            bam.write(read)
    pysam.index(str(path))
    request = ReadPlotRequest(str(path), "chr", 1000, 1, 200, 1, False, frozenset({"primary"}),
                              filters=_expression(">", 85))
    data = prepare_read_plot(request)
    assert data.fetched == 2 and len(data.alignments) == 1
    assert data.alignments[0].name in {"r0", "r1"}
