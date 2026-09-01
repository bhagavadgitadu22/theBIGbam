import pytest

from thebigbam.plotting.models.filters import FilterExpression, FilterPredicate, FilterSection
from thebigbam.plotting.settings.history import HistoryEntry, SessionHistory, describe_history_entry
from thebigbam.plotting.settings.persistence import serialize_filter_expression


def _expression(value=10):
    return FilterExpression(
        (FilterSection((FilterPredicate("Contig", "Length", ">", value),)),)
    )


def test_committed_filter_expression_uses_settings_schema():
    assert serialize_filter_expression(_expression()) == [
        {
            "section_and_or": None,
            "rows": [
                {
                    "category": "Contig",
                    "column": "Length",
                    "operator": ">",
                    "value": 10,
                    "row_and_or": None,
                }
            ],
        }
    ]


def test_history_is_chronological_but_can_be_split_by_action():
    history = SessionHistory("example.db")
    first = history.append("apply_filters", {"filtering": []})
    second = history.append("apply_plot", {"selection": {"contig": "c1"}})

    assert history.entries == (first, second)
    assert history.for_action("apply_filters") == (first,)
    assert history.for_action("apply_plot") == (second,)
    document = history.document()
    assert document["_meta"]["format"] == "thebigbam-session-history"
    assert [entry["action"] for entry in document["entries"]] == ["apply_filters", "apply_plot"]


def test_history_copies_settings_and_bounds_each_action():
    history = SessionHistory("example.db", limit_per_action=1)
    settings = {"selection": {"contig": "c1"}}
    history.append("apply_plot", settings)
    settings["selection"]["contig"] = "changed"
    latest = history.append("apply_plot", {"selection": {"contig": "c2"}})

    assert history.entries == (latest,)
    assert latest.settings["selection"]["contig"] == "c2"


def test_history_rejects_unknown_actions():
    with pytest.raises(ValueError, match="unsupported history action"):
        SessionHistory("example.db").append("widget_change", {})


def test_removing_history_entry_preserves_other_entries_and_sequence():
    history = SessionHistory("example.db")
    removed = history.append("apply_filters", {"filtering": []})
    retained = history.append("apply_plot", {"selection": {"contig": "c1"}})

    assert history.remove(removed.id) is True
    assert history.remove(removed.id) is False
    assert history.entries == (retained,)
    assert retained.sequence == 2
    assert [entry["id"] for entry in history.document()["entries"]] == [retained.id]


def test_filter_history_description_lists_only_useful_condition_details():
    entry = HistoryEntry(
        "id",
        1,
        "apply_filters",
        "2026-09-01T14:32:08+00:00",
        {
            "filtering": [
                {
                    "rows": [
                        {
                            "category": "Contig",
                            "column": "Length",
                            "operator": ">",
                            "value": 10000,
                        }
                    ]
                }
            ]
        },
    )

    assert describe_history_entry(entry) == (
        "1 condition:\n"
        "• Contig · Length > 10000"
    )


def test_plot_history_description_uses_scope_subject_variables_and_filters():
    entry = HistoryEntry(
        "id",
        2,
        "apply_plot",
        "2026-09-01T14:32:08+00:00",
        {
            "view_mode": {"mag_or_contig": 1, "one_or_all_samples": 0},
            "selection": {"contig": "NODE_42", "sample": "sample-A"},
            "contig": {"position_range": {"from": 1, "to": 50000}},
            "variables": {
                "Coverage": {"selected_one": ["Coverage", "GC content"]},
            },
            "filtering": [],
        },
    )

    description = describe_history_entry(entry)
    assert "Mode: One sample" in description
    assert "Contig: NODE_42" in description
    assert "Sample: sample-A" in description
    assert "Range: 1–50000" in description
    assert "Variables: Coverage, GC content" in description
    assert "Filters: No active filters" in description
    assert "Generated plot" not in description
    assert "Time:" not in description


def test_description_truncates_long_variable_lists():
    entry = HistoryEntry(
        "id",
        3,
        "apply_plot",
        "2026-09-01T14:32:08+00:00",
        {
            "view_mode": {"mag_or_contig": 0, "one_or_all_samples": 1},
            "selection": {"mag": "MAG-12"},
            "variables": {"module": {"selected_all": [f"v{i}" for i in range(7)]}},
        },
    )

    description = describe_history_entry(entry)
    assert "Mode: All samples" in description
    assert "MAG: MAG-12" in description
    assert "Variables: v0, v1, v2, v3, v4 (+2 more)" in description
