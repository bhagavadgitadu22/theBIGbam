from types import SimpleNamespace

from thebigbam.plotting.renderers.filter_distributions import FilterVisualizations


def test_numeric_visualization_skips_non_numeric_columns():
    visualizations = FilterVisualizations(
        SimpleNamespace(),
        {"Sample": {"columns": {"Group": {"type": "text"}}}},
        False,
        "",
        lambda: None,
    )
    row = {}

    assert visualizations.build_numeric_histogram(row, "Sample", "Group", object()) is None
    assert row["histogram_pane"] is None
    assert row["threshold_span"] is None


def test_text_visualization_handles_empty_distribution():
    service = SimpleNamespace(value_counts=lambda *_args: [])
    visualizations = FilterVisualizations(service, {}, False, "", lambda: None)
    row = {}

    assert visualizations.build_text_treemap(row, "Sample", "Group", {}, object()) is None
    assert row["treemap_pane"] is None
