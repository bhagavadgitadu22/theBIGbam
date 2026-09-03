import pytest

from thebigbam.plotting.models.filters import FilterExpression, FilterPredicate, FilterSection
from thebigbam.plotting.settings.history import (
    HistoryEntry,
    SessionHistory,
    entries_from_session_document,
    history_description_lines,
    history_diff_lines,
)
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


def test_git_like_description_reports_changes_and_struck_removals():
    first = HistoryEntry(
        "first",
        1,
        "apply_plot",
        "time",
        {"selection": {"contig": "c1", "sample": "s1"}},
    )
    second = HistoryEntry(
        "second",
        2,
        "apply_plot",
        "time",
        {"selection": {"contig": "c2"}},
    )

    lines = history_diff_lines(first, second)

    assert [line.text for line in lines] == [
        "Selection · Contig: c1 → c2",
        "Selection · Sample: s1",
    ]
    assert [line.removed for line in lines] == [False, True]


def test_complete_description_uses_canonical_settings_and_ignores_metadata():
    entry = HistoryEntry(
        "id",
        1,
        "apply_plot",
        "time",
        {"_meta": {"saved_at": "volatile"}, "selection": {"mag": "MAG-1"}},
    )

    assert [line.text for line in history_description_lines(entry)] == [
        "Selection · MAG: MAG-1"
    ]


def test_removing_history_entry_preserves_other_entries_and_sequence():
    history = SessionHistory("example.db")
    removed = history.append("apply_filters", {"filtering": []})
    retained = history.append("apply_plot", {"selection": {"contig": "c1"}})

    assert history.remove(removed.id) is True
    assert history.remove(removed.id) is False
    assert history.entries == (retained,)
    assert retained.sequence == 2
    assert [entry["id"] for entry in history.document()["entries"]] == [retained.id]


def test_saved_session_entries_can_be_loaded_with_identity_and_order_preserved():
    history = SessionHistory("example.db")
    first = history.append("apply_filters", {"filtering": []})
    second = history.append("apply_plot", {"selection": {"contig": "c1"}})

    assert entries_from_session_document(history.document()) == (first, second)


@pytest.mark.parametrize(
    "change, message",
    [
        (lambda document: document["_meta"].update(version=999), "version"),
        (lambda document: document.update(entries={}), "must be a list"),
        (lambda document: document["entries"][0].update(action="unknown"), "unsupported action"),
        (lambda document: document["entries"][0].update(settings=[]), "invalid settings"),
    ],
)
def test_invalid_saved_sessions_are_rejected(change, message):
    history = SessionHistory("example.db")
    history.append("apply_filters", {"filtering": []})
    document = history.document()
    change(document)

    with pytest.raises(ValueError, match=message):
        entries_from_session_document(document)
