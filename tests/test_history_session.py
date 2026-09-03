from types import SimpleNamespace

import panel as pn
from bokeh.events import ButtonClick
from bokeh.models import Div

from thebigbam.plotting.application.history_session import build_history_session
from thebigbam.plotting.settings.history import HistoryEntry


def test_history_row_has_restore_and_coloring_style_minus_actions():
    restored = []
    session = build_history_session(
        db_path="example.db",
        restore_entry=restored.append,
        stylesheet="",
    )
    entry = session.append(
        "apply_plot", {"selection": {"contig": "c1"}}, apply_step=23
    )

    history_row = session.plot_entries.objects[0]
    description_cell = history_row.objects[0]
    assert description_cell.margin == 0
    assert description_cell.css_classes == ["history-description"]
    assert [button.name for button in history_row.objects[1:]] == ["←", "\u2212"]
    assert [button.width for button in history_row.objects[1:]] == [None, 30]
    assert [button.sizing_mode for button in history_row.objects[1:]] == [None, "fixed"]
    assert [button.css_classes for button in history_row.objects[1:]] == [
        ["history-restore"],
        ["history-action"],
    ]
    assert "control-row" in history_row.css_classes
    assert history_row.margin == 0
    assert history_row.objects[1].description.content == "Restore this state"
    assert history_row.objects[1].description.position == "left"
    assert history_row.objects[2].description.content == "Remove from history"
    assert history_row.objects[2].description.position == "left"
    assert entry in session.history.entries


def test_history_remove_refreshes_ui_and_export_and_respects_busy_state():
    busy = {"value": True}
    session = build_history_session(
        db_path="example.db",
        restore_entry=lambda _entry: None,
        can_mutate=lambda: not busy["value"],
        stylesheet="",
    )
    entry = session.append("apply_filters", {"filtering": []})

    session._remove(entry.id)
    assert entry in session.history.entries

    busy["value"] = False
    session._remove(entry.id)
    assert entry not in session.history.entries
    assert all(item["id"] != entry.id for item in session.history.document()["entries"])
    assert session.filter_entries.objects[0].object == "_No entries yet._"


def test_history_save_style_and_collapsible_sections_match_sidebars():
    session = build_history_session(
        db_path="example.db",
        restore_entry=lambda _entry: None,
        stylesheet="",
        toggle_stylesheet="",
    )

    save = session.drawer.objects[5].objects[0].objects[0]
    assert session.drawer.width == 250
    assert session.drawer.css_classes == ["history-drawer", "sidebar-content"]
    assert save.css_classes == ["action-primary", "apply-btn", "save-session-btn"]
    assert save.button_type == "primary"
    assert save.sizing_mode != "stretch_width"
    assert save.margin == (5, 0, 0, 0)
    assert save.description.content == "Save all history of applied actions (filters and plots)"
    assert save.description.position == "left"
    assert session.filter_entries.visible is True
    assert session.plot_entries.visible is True
    assert session.filter_entries.css_classes == ["history-entry-list"]
    assert session.plot_entries.css_classes == ["history-entry-list"]
    assert session.filter_entries.stylesheets == [""]
    assert session.plot_entries.stylesheets == [""]
    assert session.filter_toggle.label == "▼"
    assert session.plot_toggle.label == "▼"
    assert session.filter_short_descriptions.value is True
    assert session.plot_short_descriptions.value is True

    filter_section = session.drawer.objects[1]
    divider = session.drawer.objects[2]
    plot_section = session.drawer.objects[3]
    assert divider.object.styles["background-color"] == "#333"
    for section in (filter_section, plot_section):
        header = section.objects[0].object
        assert "section-header" in header.css_classes
        assert "section-title" in header.children[1].css_classes
        assert header.children[0].width == header.children[2].width
        assert section.objects[1].name == "Use git-like short descriptions"
        assert section.objects[1].css_classes == ["history-description-mode"]
    assert session.drawer.objects[4].height == 10
    assert session.drawer.objects[4].css_classes == ["history-save-top-space"]
    assert session.drawer.objects[6].height == 20
    assert session.drawer.objects[6].css_classes == ["history-bottom-space"]


def test_history_save_row_can_host_settings_and_session_actions_together():
    settings_save = SimpleNamespace(
        button=pn.widgets.Button(
            name="SAVE SETTINGS",
            description="Save current settings only, including unapplied ones",
            margin=(5, 0, 0, 0),
        ),
        confirmation=Div(text="", visible=False),
    )
    session = build_history_session(
        db_path="example.db",
        restore_entry=lambda _entry: None,
        stylesheet="",
        settings_save_controls=settings_save,
    )

    save_group = session.drawer.objects[5]
    assert isinstance(save_group, pn.Column)
    save_rows = save_group.objects[:2]
    assert [row.objects[0].name for row in save_rows] == ["SAVE SESSION", "SAVE SETTINGS"]
    assert save_rows[1].objects[0].description == (
        "Save current settings only, including unapplied ones"
    )
    assert all(row.objects[0].sizing_mode != "stretch_width" for row in save_rows)
    assert all(row.objects[0].margin == (5, 0, 0, 0) for row in save_rows)
    assert all(row.stylesheets == [""] for row in save_rows)
    confirmation_carrier = save_group.objects[2]
    assert (confirmation_carrier.width, confirmation_carrier.height) == (0, 0)
    assert confirmation_carrier.css_classes == ["save-confirmation-carrier"]


def test_history_sections_fold_independently_and_stay_folded_on_refresh():
    session = build_history_session(
        db_path="example.db",
        restore_entry=lambda _entry: None,
        stylesheet="",
    )

    session.filter_toggle._trigger_event(ButtonClick(model=session.filter_toggle))
    session.append("apply_filters", {"filtering": []})

    assert session.filter_entries.visible is False
    assert session.filter_short_descriptions.visible is False
    assert session.filter_toggle.label == "▶"
    assert session.plot_entries.visible is True
    assert session.plot_toggle.label == "▼"


def test_history_description_uses_safe_full_text_title_without_html_injection():
    session = build_history_session(
        db_path="example.db",
        restore_entry=lambda _entry: None,
        stylesheet="",
    )
    session.append(
        "apply_plot",
        {
            "view_mode": {"mag_or_contig": 1, "one_or_all_samples": 0},
            "selection": {"contig": "<script>alert(1)</script>"},
        },
    )

    line_html = session.plot_entries.objects[0].objects[0].objects[2].object
    assert 'class="history-description-line"' in line_html
    assert "&lt;script&gt;alert(1)&lt;/script&gt;" in line_html
    assert "<script>" not in line_html


def test_append_updates_only_affected_section_and_preserves_row_identity():
    session = build_history_session(
        db_path="example.db",
        restore_entry=lambda _entry: None,
        stylesheet="",
    )
    first_plot = session.append("apply_plot", {"selection": {"contig": "c1"}})
    first_plot_row = session.rows_by_id[first_plot.id]
    plot_objects = session.plot_entries.objects

    session.append("apply_filters", {"filtering": []})

    assert session.rows_by_id[first_plot.id] is first_plot_row
    assert session.plot_entries.objects is plot_objects

    second_plot = session.append("apply_plot", {"selection": {"contig": "c2"}})
    assert session.plot_entries.objects == [first_plot_row, session.rows_by_id[second_plot.id]]
    assert session.rows_by_id[first_plot.id] is first_plot_row


def test_description_mode_switches_between_diff_and_complete_and_remove_recomputes():
    session = build_history_session(
        db_path="example.db",
        restore_entry=lambda _entry: None,
        stylesheet="",
    )
    first = session.append("apply_plot", {"selection": {"contig": "c1"}})
    second = session.append("apply_plot", {"selection": {"contig": "c2"}})

    second_line = session.rows_by_id[second.id].objects[0].objects[0].object
    assert "c1 → c2" in second_line

    session.plot_short_descriptions.value = False
    second_line = session.rows_by_id[second.id].objects[0].objects[0].object
    assert "c1 → c2" not in second_line
    assert "c2" in second_line

    session.plot_short_descriptions.value = True
    session._remove(first.id)
    second_line = session.rows_by_id[second.id].objects[0].objects[0].object
    assert "c1 → c2" not in second_line
    assert "c2" in second_line


def test_initial_history_is_rebuilt_oldest_first_without_creating_new_entries():
    entries = (
        HistoryEntry("filters", 1, "apply_filters", "2026-09-01T10:00:00+00:00", {"filtering": []}),
        HistoryEntry(
            "plot-old", 2, "apply_plot", "2026-09-01T10:01:00+00:00", {"selection": {"contig": "c1"}}
        ),
        HistoryEntry(
            "plot-new", 3, "apply_plot", "2026-09-01T10:02:00+00:00", {"selection": {"contig": "c2"}}
        ),
    )

    session = build_history_session(
        db_path="example.db",
        restore_entry=lambda _entry: None,
        stylesheet="",
        initial_entries=entries,
    )

    assert session.history.entries == entries
    assert session.filter_entries.objects == [session.rows_by_id["filters"]]
    assert session.plot_entries.objects == [
        session.rows_by_id["plot-old"],
        session.rows_by_id["plot-new"],
    ]
    appended = session.append("apply_filters", {"filtering": []})
    assert appended.sequence == 4


def test_retention_removes_only_expired_row_and_keeps_models_bounded():
    session = build_history_session(
        db_path="example.db",
        restore_entry=lambda _entry: None,
        stylesheet="",
    )
    session.history.limit_per_action = 2
    first = session.append("apply_plot", {"selection": {"contig": "c1"}})
    second = session.append("apply_plot", {"selection": {"contig": "c2"}})
    second_row = session.rows_by_id[second.id]
    third = session.append("apply_plot", {"selection": {"contig": "c3"}})

    assert first.id not in session.rows_by_id
    assert session.rows_by_id[second.id] is second_row
    assert session.plot_entries.objects == [second_row, session.rows_by_id[third.id]]
    assert len(session.history.for_action("apply_plot")) == 2


def test_restore_and_remove_emit_explicit_scenario_actions():
    actions = []
    restored = []
    session = build_history_session(
        db_path="example.db",
        restore_entry=restored.append,
        record_action=lambda action, details: actions.append((action, details)),
        stylesheet="",
    )
    entry = session.append(
        "apply_plot", {"selection": {"contig": "c1"}}, apply_step=23
    )

    session._restore(entry)
    session._remove(entry.id)

    assert restored == [entry]
    assert actions == [
        (
            "restore_history",
            {
                "history_sequence": entry.sequence,
                "history_action": "apply_plot",
                "apply_step": 23,
                "settings": entry.settings,
            },
        ),
        (
            "remove_history",
            {
                "history_sequence": entry.sequence,
                "history_action": "apply_plot",
                "apply_step": 23,
            },
        ),
    ]


def test_saved_session_history_actions_keep_explicit_pre_scenario_origin():
    actions = []
    entry = HistoryEntry(
        "saved-plot",
        4,
        "apply_plot",
        "2026-09-01T10:00:00+00:00",
        {"selection": {"contig": "c1"}},
    )
    session = build_history_session(
        db_path="example.db",
        restore_entry=lambda _entry: None,
        record_action=lambda action, details: actions.append((action, details)),
        stylesheet="",
        initial_entries=(entry,),
    )

    session._restore(entry)
    session._remove(entry.id)

    assert [details["apply_step"] for _action, details in actions] == [None, None]


def test_save_session_persists_and_records_explicit_action_only_when_clicked(monkeypatch, tmp_path):
    actions = []
    saved = []
    monkeypatch.setattr(
        "thebigbam.plotting.application.history_session.save_session_document",
        lambda document, db_path: saved.append((document, db_path)) or tmp_path / "saved.json",
    )
    session = build_history_session(
        db_path="example.db",
        restore_entry=lambda _entry: None,
        record_action=lambda action, details: actions.append((action, details)),
        stylesheet="",
    )
    save_row = session.drawer.objects[5].objects[0]
    save = save_row.objects[0]
    confirmation = session.drawer.objects[5].objects[1].objects[0].object

    assert actions == []
    save.clicks += 1

    assert saved[0][0]["_meta"]["format"] == "thebigbam-session-history"
    assert saved[0][1] == "example.db"
    assert confirmation.text == "1"
    assert actions == [("save_session", {})]
