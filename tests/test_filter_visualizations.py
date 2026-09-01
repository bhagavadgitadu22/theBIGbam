from types import SimpleNamespace

import numpy as np
import panel as pn

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


def test_distribution_scale_records_explicit_axis_state(monkeypatch):
    actions = []
    monkeypatch.setattr(
        "thebigbam.plotting.renderers.filter_distributions.curdoc",
        lambda: SimpleNamespace(add_next_tick_callback=lambda _callback: None),
    )
    service = SimpleNamespace(
        histogram_bins=lambda *_args, **_kwargs: SimpleNamespace(
            edges=np.array([1.0, 2.0, 3.0]),
            counts=np.array([2, 3]),
            approximate=False,
            sampled_rows=5,
        ),
        column_null_stats=lambda *_args: None,
    )
    visualizations = FilterVisualizations(
        service,
        {"Contig": {"columns": {"Length": {"type": "numeric", "is_bool": False}}}},
        False,
        "",
        lambda: None,
        record_action=lambda action, details: actions.append((action, details)),
        occurrence_for=lambda _row: 2,
    )
    row = {"hist_container": pn.Column(), "loading_gen": 0}
    result = visualizations.build_numeric_histogram(
        row,
        "Contig",
        "Length",
        SimpleNamespace(step=1, value=1),
    )

    result[0].objects[0].objects[0].clicks += 1

    assert actions == [
        (
            "filter_distribution_scale",
            {
                "category": "Contig",
                "column": "Length",
                "occurrence": 2,
                "axis": "x",
                "enabled": True,
            },
        )
    ]
