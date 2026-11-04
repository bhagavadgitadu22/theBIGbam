"""Read predicates run before spatial sampling and preserve AND/OR grouping."""

from dataclasses import replace

import pytest

from thebigbam.plotting.models.filters import FilterExpression, FilterPredicate, FilterSection
from thebigbam.plotting.models.read_candidate import ReadCandidate
from thebigbam.plotting.models.read_plot import ReadPlotRequest
from thebigbam.plotting.services.read_filtering import compile_read_filters, validate_read_filters
from thebigbam.plotting.services.read_plots import prepare_read_plot
from thebigbam.plotting.shared.data_cache import SessionDataCache
from thebigbam.plotting.shared.read_diagnostics import ReadProfiler


def _read(**kwargs):
    return replace(ReadCandidate("read", 10, 20, "2S6M1I1D3M2H", "primary", 33, 42, 40,
                                  template_length=-50), **kwargs)


def _rule(field, operator, value):
    return FilterPredicate("Reads", field, operator, value)


def _expression(*rules, connectors=()):
    return FilterExpression((FilterSection(rules, connectors),))


@pytest.mark.parametrize("field,value", [
    ("Read length", 14), ("Alignment length", 10), ("Alignment start", 11),
    ("Alignment end", 20), ("MAPQ", 42), ("Insert sizes", 50),
])
def test_numeric_read_fields(field, value):
    assert compile_read_filters(_expression(_rule(field, "=", value)), 100, False)(_read())
    assert compile_read_filters(_expression(_rule(field, ">", value - 1)), 100, False)(_read())
    assert not compile_read_filters(_expression(_rule(field, "!=", value)), 100, False)(_read())


def test_alignment_includes_is_one_based_and_wraps_circular_origin():
    read = _read(start=95, end=105, cigar="10M")
    for point in (96, 100, 1, 5):
        assert compile_read_filters(_expression(_rule("Alignment includes", "=", point)), 100, True)(read)
    assert compile_read_filters(_expression(_rule("Alignment includes", "!=", 50)), 100, True)(read)
    assert not compile_read_filters(_expression(_rule("Alignment includes", "=", 95)), 100, True)(read)
    assert compile_read_filters(_expression(_rule("Alignment end", "=", 105)), 100, True)(read)
    assert compile_read_filters(_expression(_rule("Alignment includes", "=", 17)), 100, False)(_read())


def test_boolean_flags_and_circular_pair_orientation():
    assert compile_read_filters(_expression(_rule("Has mate mapped", "=", "TRUE")), 100, False)(_read())
    assert compile_read_filters(_expression(_rule("Non-inward pair", "=", "FALSE")), 100, False)(_read())
    assert compile_read_filters(_expression(_rule("Non-inward pair", "=", "TRUE")), 100, False)(_read(flag=1))
    assert not compile_read_filters(_expression(_rule("Has mate mapped", "=", "TRUE")), 100, False)(_read(flag=9))
    wrapped = _read(start=95, end=105, mate_start=5)
    assert compile_read_filters(_expression(_rule("Circular", "=", "TRUE")), 100, True)(wrapped)
    assert not compile_read_filters(_expression(_rule("Circular", "=", "TRUE")), 100, False)(wrapped)
    inward = _expression(_rule("Non-inward pair", "=", "FALSE"))
    assert compile_read_filters(inward, 100, True)(wrapped)
    assert not compile_read_filters(inward, 100, False)(wrapped)


def test_unknown_insert_size_does_not_match_numeric_inequality():
    expression = _expression(_rule("Insert sizes", "!=", 10))
    for read in (_read(template_length=0), _read(flag=0), _read(mate_same_reference=False)):
        assert not compile_read_filters(expression, 100, False)(read)


def test_block_and_row_connectors_follow_existing_left_to_right_grouping():
    a, b, c = _rule("MAPQ", ">", 40), _rule("Read length", "<", 10), _rule("Alignment start", ">", 20)
    first = FilterSection((a, b, c), ("OR", "AND"))
    expression = FilterExpression((first, FilterSection((_rule("Insert sizes", "=", 50),))), ("OR",))
    assert not compile_read_filters(FilterExpression((first,)), 100, False)(_read())
    assert compile_read_filters(expression, 100, False)(_read())
    assert compile_read_filters(FilterExpression(()), 100, False)(_read())


@pytest.mark.parametrize("rule", [
    _rule("MAPQ", "=", float("nan")), _rule("MAPQ", "=", float("inf")),
    _rule("Alignment includes", ">", 5), _rule("Alignment includes", "=", 1.5),
    _rule("Alignment includes", "=", 0), _rule("Has mate mapped", "=", "maybe"),
])
def test_invalid_read_filters_are_rejected(rule):
    with pytest.raises(ValueError):
        validate_read_filters(_expression(rule))


def test_filters_change_cache_identity_and_run_before_subsampling(monkeypatch):
    calls = []
    reads = tuple(_read(name=f"r{i}", mapq=60 if i >= 95 else 0) for i in range(100))

    def fetch(*args, **kwargs):
        calls.append(1)
        return iter(reads)

    monkeypatch.setattr("thebigbam.plotting.services.read_plots.fetch_read_candidates", fetch)
    monkeypatch.setattr("thebigbam.plotting.services.read_plots._fingerprint", lambda _: (1, 1))
    request = ReadPlotRequest("fake.bam", "chr", 100, 1, 100, 1, False, frozenset({"primary"}))
    cache = SessionDataCache()
    prepare_read_plot(request, cache)
    filtered = replace(request, max_reads=3, filters=_expression(_rule("MAPQ", ">", 20)), profiler=ReadProfiler())
    data = prepare_read_plot(filtered, cache)
    assert len(calls) == 2
    assert data.fetched == 5 and len(data.alignments) == 3
    assert all(read.mapq == 60 for read in data.alignments)
    assert filtered.profiler.counts["read_filter_examined"] == 100
    assert filtered.profiler.counts["read_filter_matched"] == 5
    larger = replace(filtered, max_reads=10, profiler=ReadProfiler())
    assert len(prepare_read_plot(larger, cache).alignments) == 5
    assert len(calls) == 2 and larger.profiler.counts["candidate_cache_hits"] == 1
