import pytest

from thebigbam.plotting.models.filters import FilterExpression, FilterPredicate, FilterSection
from thebigbam.plotting.settings.history import (
    HistoryEntry,
    SessionHistory,
    entries_from_session_document,
)
from thebigbam.plotting.settings.history_descriptions import (
    HistoryDescriptionContext,
    canonical_history_description_lines,
    diff_description_lines,
    history_description_lines,
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

    lines = diff_description_lines(
        canonical_history_description_lines(first),
        canonical_history_description_lines(second),
    )

    assert [line.text for line in lines] == [
        "Selection · Contig: c2",
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


@pytest.mark.parametrize(
    ("view_mode", "expected"),
    [
        (
            {"mag_or_contig": 0, "one_or_all_samples": 0},
            [
                "View mode · MAG/contig view: MAG",
                "View mode · Sample mode: One sample",
            ],
        ),
        (
            {"mag_or_contig": 1, "one_or_all_samples": 1},
            [
                "View mode · MAG/contig view: Contig",
                "View mode · Sample mode: All samples",
            ],
        ),
    ],
)
def test_history_description_translates_encoded_view_modes(view_mode, expected):
    entry = HistoryEntry("id", 1, "apply_plot", "time", {"view_mode": view_mode})

    assert [line.text for line in history_description_lines(entry)] == expected


def test_history_description_orders_selection_and_omits_inactive_controls():
    entry = HistoryEntry(
        "id",
        1,
        "apply_plot",
        "time",
        {
            "selection": {
                "sample": "GL34_UP",
                "contig": None,
                "mag": "337R_metabat_172",
            },
            "contig": {
                "feature_widgets": {"sequence": False, "translated_sequence": None},
                "coloring": {
                    "custom_color_rows": [
                        {
                            "qualifier": "activity",
                            "operator": "=",
                            "value": "",
                            "color": "#cccccc",
                        }
                    ]
                },
            },
        },
    )

    assert [line.text for line in history_description_lines(entry)] == [
        "Selection · MAG: 337R_metabat_172",
        "Selection · Sample: GL34_UP",
    ]


def test_history_description_compacts_position_and_coloring_rules():
    entry = HistoryEntry(
        "id",
        1,
        "apply_plot",
        "time",
        {
            "contig": {
                "position_range": {"from": 1, "to": 3254835},
                "coloring": {
                    "custom_color_rows": [
                        {
                            "qualifier": "activity",
                            "operator": "has",
                            "value": "Defense",
                            "color": "#abcdef",
                        }
                    ],
                    "mag_track_color_rows": [
                        {
                            "qualifier": "Coverage_mean",
                            "operator": ">",
                            "value": 10,
                            "color": "#123456",
                        }
                    ],
                },
            }
        },
    )

    assert [line.text for line in history_description_lines(entry)] == [
        "Contig · Position range: 1–3254835",
        "Contig · Coloring · #1: activity has Defense · #abcdef",
        "Contig · MAG coloring · #1: Coverage_mean > 10 · #123456",
    ]


def test_plotting_defaults_are_hidden_but_reverting_to_one_is_a_change():
    defaults = HistoryEntry(
        "defaults",
        1,
        "apply_plot",
        "time",
        {
            "view_mode": {"mag_or_contig": 1, "one_or_all_samples": 1},
            "plotting_params": {
                "max_genemap_window": 100_000,
                "sample_params": {"max_samples": 20, "same_y_scale": False},
            }
        },
    )
    customized = HistoryEntry(
        "customized",
        2,
        "apply_plot",
        "time",
        {
            "view_mode": {"mag_or_contig": 1, "one_or_all_samples": 1},
            "plotting_params": {
                "max_genemap_window": 50_000,
                "sample_params": {"max_samples": 40, "same_y_scale": True},
            }
        },
    )

    assert [line.text for line in history_description_lines(defaults)] == [
        "View mode · MAG/contig view: Contig",
        "View mode · Sample mode: All samples",
    ]
    default_lines = canonical_history_description_lines(defaults)
    customized_lines = canonical_history_description_lines(customized)
    assert [line.text for line in diff_description_lines((), default_lines)] == [
        "View mode · MAG/contig view: Contig",
        "View mode · Sample mode: All samples",
    ]
    assert [line.text for line in diff_description_lines(default_lines, customized_lines)] == [
        "Plotting parameters · Max genemap window: 50000",
        "Plotting parameters · Sample params · Max samples: 40",
        "Plotting parameters · Sample params · Same y scale: Enabled",
    ]
    assert [line.text for line in diff_description_lines(customized_lines, default_lines)] == [
        "Plotting parameters · Max genemap window: 100000",
        "Plotting parameters · Sample params · Max samples: 20",
        "Plotting parameters · Sample params · Same y scale: Disabled",
    ]


def test_plot_history_description_includes_compact_committed_filters():
    settings = {
        "selection": {"contig": "c1"},
        "filtering": [
            {
                "rows": [
                    {
                        "category": "Contig",
                        "column": "Length",
                        "operator": ">",
                        "value": 1000,
                    }
                ]
            }
        ],
    }
    entry = HistoryEntry("plot", 1, "apply_plot", "time", settings)

    assert [line.text for line in history_description_lines(entry)] == [
        "Selection · Contig: c1",
        "Filter #1.1: Contig · Length > 1000",
    ]
    assert entry.settings["filtering"] == settings["filtering"]


def test_plot_diff_reports_a_filter_only_change():
    first = HistoryEntry(
        "first",
        1,
        "apply_plot",
        "time",
        {"selection": {"contig": "c1"}, "filtering": []},
    )
    second = HistoryEntry(
        "second",
        2,
        "apply_plot",
        "time",
        {
            "selection": {"contig": "c1"},
            "filtering": [
                {
                    "rows": [
                        {
                            "category": "Annotations",
                            "column": "activity",
                            "operator": "has",
                            "value": "Defense",
                        }
                    ]
                }
            ],
        },
    )

    assert [
        line.text
        for line in diff_description_lines(
            canonical_history_description_lines(first),
            canonical_history_description_lines(second),
        )
    ] == ["Filter #1.1: Annotations · activity has Defense"]


def test_history_description_always_capitalizes_mag_acronym():
    entry = HistoryEntry(
        "plot",
        1,
        "apply_plot",
        "time",
        {
            "contig": {
                "coloring": {"apply_annotation_rules_to_mag_track": True}
            }
        },
    )

    assert [line.text for line in history_description_lines(entry)] == [
        "Contig · Coloring · Apply annotation rules to MAG track: Enabled"
    ]


def test_plot_description_uses_visible_parameter_groups_and_compact_variables():
    context = HistoryDescriptionContext(
        mag_order=("Coverage", "Coverage_mean"),
        sample_order=("Sample", "Sample_name"),
    )
    one_sample = HistoryEntry(
        "one",
        1,
        "apply_plot",
        "time",
        {
            "view_mode": {"mag_or_contig": 0, "one_or_all_samples": 0},
            "variables": {
                "Coverage": {
                    "module_enabled": True,
                    "selected_one": ["Primary alignments", "Alignments by strand"],
                    "selected_all": ["Other alignments"],
                }
            },
            "plotting_params": {
                "sample_params": {
                    "max_samples": 20,
                    "order_category": "Sample",
                    "order_metric": "Sample_name",
                    "order_direction": 0,
                },
                "mag_params": {
                    "category": "Coverage",
                    "metric": "Coverage_mean",
                    "direction": 1,
                    "sort_sample": "unused-in-one-sample",
                    "max_dots": 1000,
                },
            },
        },
    )
    partial = HistoryEntry(
        "partial",
        2,
        "apply_plot",
        "time",
        {
            "view_mode": {"mag_or_contig": 1, "one_or_all_samples": 0},
            "variables": {
                "Coverage": {
                    "module_enabled": False,
                    "selected_one": ["Primary alignments", "Alignments by strand"],
                    "selected_all": ["Other alignments"],
                }
            },
        },
    )

    one_lines = [line.text for line in history_description_lines(one_sample, context)]
    assert one_lines == [
        "View mode · MAG/contig view: MAG",
        "View mode · Sample mode: One sample",
        "Variables · Coverage: All",
    ]
    assert [line.text for line in history_description_lines(partial, context)] == [
        "View mode · MAG/contig view: Contig",
        "View mode · Sample mode: One sample",
        "Variables · Coverage: Primary alignments, Alignments by strand",
    ]


def test_plot_description_combines_non_default_visible_ordering():
    context = HistoryDescriptionContext(
        mag_order=("Contig", "Contig_length"),
        sample_order=("Sample", "Sample_name"),
    )
    mag_entry = HistoryEntry(
        "mag",
        1,
        "apply_plot",
        "time",
        {
            "view_mode": {"mag_or_contig": 0, "one_or_all_samples": 0},
            "plotting_params": {
                "mag_params": {
                    "category": "Coverage",
                    "metric": "Coverage_mean",
                    "direction": 0,
                    "max_dots": 1000,
                }
            },
        },
    )

    assert [line.text for line in history_description_lines(mag_entry, context)][-1] == (
        "Plotting parameters · MAG order: Coverage · Coverage_mean · ↑"
    )


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
