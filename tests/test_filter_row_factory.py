from types import SimpleNamespace

import pytest
from bokeh.models import InlineStyleSheet

from thebigbam.plotting.controls.filter_rows import FilterRowFactory
from thebigbam.plotting.controls.searchable_select import SearchableSelect


class Visualizations:
    def build_numeric_histogram(self, *args, **kwargs):
        return None

    def build_text_treemap(self, *args, **kwargs):
        return None


def test_filter_row_factory_builds_row_from_cached_column_options():
    metadata = {
        "Contig": {
            "columns": {"Contig_length": {"type": "numeric", "is_bool": False}},
        }
    }
    factory = FilterRowFactory(
        metadata_service=SimpleNamespace(),
        filtering_metadata=metadata,
        columns={"Contig": [("Contig_length", "Contig length")]},
        raw_columns={"Contig": ["Contig_length"]},
        visualizations=Visualizations(),
        refresh=lambda: None,
        stylesheet=InlineStyleSheet(css=":host { color: red; }"),
        enable_timing=False,
    )
    factory.attach_controller(SimpleNamespace(count_rows=lambda: 1, sections=[]))

    result = factory.create_row({}, initial_category="Contig", initial_column="Contig_length", initial_operator=">")

    assert result["category_select"].value == "Contig"
    assert result["subcategory_select"].value == "Contig_length"
    assert result["comparison_select"].value == ">"
    assert result["input_ref"]["widget"].value is None
    assert result["minus_btn"].stylesheets == [":host { color: red; }"]


def test_lookup_records_semantic_target(monkeypatch):
    callbacks = []
    events = []
    monkeypatch.setattr(
        "thebigbam.plotting.controls.filter_rows.curdoc",
        lambda: SimpleNamespace(add_next_tick_callback=callbacks.append),
    )
    metadata = {"Coverage": {"columns": {"Coverage_mean": {"type": "numeric", "is_bool": False}}}}

    def record(action, details):
        events.append((action, details))

    factory = FilterRowFactory(
        metadata_service=SimpleNamespace(),
        filtering_metadata=metadata,
        columns={"Coverage": [("Coverage_mean", "Coverage mean")]},
        raw_columns={"Coverage": ["Coverage_mean"]},
        visualizations=Visualizations(),
        refresh=lambda: None,
        stylesheet=InlineStyleSheet(css=""),
        enable_timing=False,
        record_action=record,
    )
    factory.attach_controller(SimpleNamespace(count_rows=lambda: 1, sections=[]))
    row = factory.create_row({}, initial_category="Coverage", initial_column="Coverage_mean")

    factory._toggle_distribution(row)
    callbacks.pop(0)()

    assert events == [
        (
            "filter_lookup",
            {"category": "Coverage", "column": "Coverage_mean", "opening": True},
        ),
    ]


def test_semantic_lookup_targets_current_duplicate_occurrence(monkeypatch):
    callbacks = []
    monkeypatch.setattr(
        "thebigbam.plotting.controls.filter_rows.curdoc",
        lambda: SimpleNamespace(add_next_tick_callback=callbacks.append),
    )
    metadata = {"Termini": {"columns": {"Packaging_mechanism": {"type": "text"}}}}
    factory = FilterRowFactory(
        metadata_service=SimpleNamespace(),
        filtering_metadata=metadata,
        columns={"Termini": [("Packaging_mechanism", "Packaging mechanism")]},
        raw_columns={"Termini": ["Packaging_mechanism"]},
        visualizations=Visualizations(),
        refresh=lambda: None,
        stylesheet=InlineStyleSheet(css=""),
        enable_timing=False,
    )
    first = factory.create_row({}, initial_category="Termini", initial_column="Packaging_mechanism")
    second = factory.create_row({}, initial_category="Termini", initial_column="Packaging_mechanism")
    factory.attach_controller(SimpleNamespace(sections=[{"rows": [first, second]}]))

    factory.set_distribution("Termini", "Packaging_mechanism", occurrence=2, opening=True)

    assert first["loading_gen"] == 0
    assert second["loading_gen"] == 1
    assert len(callbacks) == 1
    with pytest.raises(ValueError, match="occurrence 3 does not exist"):
        factory.set_distribution("Termini", "Packaging_mechanism", occurrence=3)


def test_category_change_updates_shared_metric_immediately_and_defers_one_refresh(monkeypatch):
    callbacks = []
    monkeypatch.setattr(
        "thebigbam.plotting.controls.filter_rows.curdoc",
        lambda: SimpleNamespace(add_next_tick_callback=callbacks.append),
    )
    metadata = {
        "Coverage": {"columns": {"Status": {"type": "numeric", "is_bool": False}}},
        "MAG coverage": {"columns": {"Status": {"type": "text", "is_bool": False}}},
    }
    refreshes = []
    factory = FilterRowFactory(
        metadata_service=SimpleNamespace(
            search_distinct_values=lambda category, column, search, limit: ["high", "low"]
        ),
        filtering_metadata=metadata,
        columns={"Coverage": [("Status", "Status")], "MAG coverage": [("Status", "Status")]},
        raw_columns={"Coverage": ["Status"], "MAG coverage": ["Status"]},
        visualizations=Visualizations(),
        refresh=lambda: refreshes.append(True),
        stylesheet=InlineStyleSheet(css=""),
        enable_timing=False,
    )
    factory.attach_controller(SimpleNamespace(count_rows=lambda: 1, sections=[]))
    row = factory.create_row({}, initial_category="Coverage", initial_column="Status")

    row["category_select"].value = "MAG coverage"

    assert row["subcategory_select"].options == [("Status", "Status")]
    assert row["comparison_select"].options == ["=", "!=", "has", "has not"]
    assert isinstance(row["input_ref"]["widget"], SearchableSelect)
    assert row["input_ref"]["widget"].placeholder == "Search..."
    assert row["input_ref"]["widget"].server_search
    assert row["input_ref"]["widget"].min_search_chars == 2
    assert row["input_ref"]["widget"].options == []
    assert refreshes == []

    while callbacks:
        callbacks.pop(0)()

    assert refreshes == [True]
    row["input_ref"]["widget"].search_query = "h"
    row["input_ref"]["widget"].search_nonce += 1
    assert row["input_ref"]["widget"].options == []
    row["input_ref"]["widget"].search_query = "hi"
    row["input_ref"]["widget"].search_nonce += 1
    assert row["input_ref"]["widget"].options == ["high", "low"]
