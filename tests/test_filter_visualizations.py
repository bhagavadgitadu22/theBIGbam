from thebigbam.plotting.controls import filter_visualizations
from thebigbam.plotting.controls.filter_visualizations import FilterVisualizations


def test_numeric_visualization_skips_non_numeric_columns():
    visualizations = FilterVisualizations(
        "example.db",
        {"Sample": {"columns": {"Group": {"type": "text"}}}},
        False,
        "",
        lambda: None,
    )
    row = {}

    assert visualizations.build_numeric_histogram(row, "Sample", "Group", object()) is None
    assert row["histogram_pane"] is None
    assert row["threshold_span"] is None


def test_text_visualization_handles_empty_distribution(monkeypatch):
    monkeypatch.setattr(filter_visualizations, "resolve_value_counts", lambda *args, **kwargs: [])
    visualizations = FilterVisualizations("example.db", {}, False, "", lambda: None)
    row = {}

    assert visualizations.build_text_treemap(row, "Sample", "Group", {}, object()) is None
    assert row["treemap_pane"] is None
