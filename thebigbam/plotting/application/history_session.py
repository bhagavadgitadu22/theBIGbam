"""Panel projection and restore hooks for per-browser-session history."""

from __future__ import annotations

import html
import os
from dataclasses import dataclass
from typing import Any, Callable

import panel as pn
from bokeh.models import Button as BokehButton
from bokeh.models import CustomJS, Div

from ..settings.controls import save_confirmation_js
from ..settings.history import HistoryEntry, SessionHistory, describe_history_entry
from ..settings.persistence import save_session_document
from ..shared.styles import (
    panel_control_row,
    panel_stylesheet,
    right_panel_tooltip,
    section_header,
)
from .layout import separator


@dataclass(frozen=True)
class HistorySession:
    history: SessionHistory
    drawer: Any
    filter_entries: Any
    plot_entries: Any
    filter_toggle: Any
    plot_toggle: Any
    restore_entry: Callable[[HistoryEntry], None]
    can_mutate: Callable[[], bool]
    record_action: Callable[[str, dict[str, Any]], None]
    save: Callable[..., Any]
    rows_by_id: dict[str, Any]
    actions_by_id: dict[str, str]

    def append(self, action: str, settings: dict[str, Any]) -> HistoryEntry:
        previous_ids = {item.id for item in self.history.entries}
        entry = self.history.append(action, settings)
        retained_ids = {item.id for item in self.history.entries}
        for expired_id in previous_ids - retained_ids:
            self._remove_row(expired_id)
        row = self._entry_row(entry)
        self.rows_by_id[entry.id] = row
        self.actions_by_id[entry.id] = action
        content = self._content_for(action)
        existing = [] if self._is_empty(content) else list(content.objects)
        content.objects = [row, *existing]
        return entry

    @staticmethod
    def _empty_row() -> Any:
        return pn.pane.Markdown("_No entries yet._", margin=(0, 5))

    @staticmethod
    def _is_empty(content: Any) -> bool:
        return len(content.objects) == 1 and getattr(content.objects[0], "object", None) == "_No entries yet._"

    def _content_for(self, action: str) -> Any:
        return self.filter_entries if action == "apply_filters" else self.plot_entries

    def _entry_row(self, entry: HistoryEntry) -> Any:
        description = describe_history_entry(entry)
        time = pn.widgets.Button(
            name=html.escape(entry.created_at[11:19]),
            description=right_panel_tooltip(description),
            sizing_mode="stretch_width",
            height=30,
            margin=0,
            css_classes=["history-time"],
            stylesheets=list(self.drawer.stylesheets),
        )
        button_stylesheets = list(self.drawer.stylesheets)
        restore = pn.widgets.Button(
            name="←",
            description=right_panel_tooltip("Restore this state"),
            height=30,
            margin=0,
            css_classes=["history-restore"],
            stylesheets=button_stylesheets,
        )
        restore.on_click(lambda _event, item=entry: self._restore(item))
        remove = pn.widgets.Button(
            name="\u2212",
            description=right_panel_tooltip("Remove from history"),
            width=30,
            height=30,
            sizing_mode="fixed",
            margin=0,
            css_classes=["history-action"],
            stylesheets=button_stylesheets,
        )
        remove.on_click(lambda _event, entry_id=entry.id: self._remove(entry_id))
        return panel_control_row(
            time,
            restore,
            remove,
            css_classes=["history-entry"],
            sizing_mode="stretch_width",
            margin=0,
        )

    def _remove(self, entry_id: str, *, record: bool = True) -> None:
        if not self.can_mutate():
            return
        entry = next((item for item in self.history.entries if item.id == entry_id), None)
        if entry is not None and self.history.remove(entry_id):
            self._remove_row(entry_id)
            if record:
                self.record_action(
                    "remove_history",
                    {"history_sequence": entry.sequence, "history_action": entry.action},
                )

    def _restore(self, entry: HistoryEntry) -> None:
        if not self.can_mutate():
            return
        self.record_action(
            "restore_history",
            {
                "history_sequence": entry.sequence,
                "history_action": entry.action,
                "settings": entry.settings,
            },
        )
        self.restore_entry(entry)

    def remove_sequence(self, sequence: int, action: str) -> bool:
        entry = next(
            (item for item in self.history.entries if item.sequence == sequence and item.action == action),
            None,
        )
        if entry is None:
            return False
        self._remove(entry.id, record=False)
        return entry.id not in self.rows_by_id

    def _remove_row(self, entry_id: str) -> None:
        row = self.rows_by_id.pop(entry_id, None)
        action = self.actions_by_id.pop(entry_id, None)
        if row is None or action is None:
            return
        content = self._content_for(action)
        retained = [item for item in content.objects if item is not row]
        content.objects = retained or [self._empty_row()]

def build_history_session(
    *,
    db_path: str,
    restore_entry: Callable[[HistoryEntry], None],
    can_mutate: Callable[[], bool] = lambda: True,
    record_action: Callable[[str, dict[str, Any]], None] = lambda _action, _details: None,
    stylesheet: Any,
    toggle_stylesheet: Any | None = None,
    settings_save_controls: Any | None = None,
) -> HistorySession:
    history = SessionHistory(os.path.basename(db_path))
    filter_entries = pn.Column(
        HistorySession._empty_row(),
        sizing_mode="stretch_width",
        css_classes=["history-entry-list"],
        stylesheets=[stylesheet],
    )
    plot_entries = pn.Column(
        HistorySession._empty_row(),
        sizing_mode="stretch_width",
        css_classes=["history-entry-list"],
        stylesheets=[stylesheet],
    )

    confirmation = Div(text="", visible=False)
    confirmation.js_on_change("text", CustomJS(code=save_confirmation_js("save-session-btn")))

    def save_session(_event=None, *, record=True):
        path = save_session_document(history.document(), db_path)
        print(f"[history] Saved to {path}", flush=True)
        confirmation.text = "1"
        if record:
            record_action("save_session", {})
        return path

    save_button = pn.widgets.Button(
        name="SAVE SESSION",
        description=right_panel_tooltip(
            "Save all history of applied actions (filters and plots)"
        ),
        button_type="primary",
        align="center",
        css_classes=["action-primary", "apply-btn", "save-session-btn"],
        margin=(5, 0, 0, 0),
        stylesheets=[stylesheet],
    )
    save_button.on_click(save_session)
    session_save_row = pn.Row(
        save_button,
        sizing_mode="stretch_width",
        align="center",
        css_classes=["action-row", "sidebar-actions"],
        stylesheets=[panel_stylesheet(stylesheet)],
    )
    save_rows = [session_save_row]
    confirmations = [confirmation]
    if settings_save_controls is not None:
        save_rows.append(
            pn.Row(
                settings_save_controls.button,
                sizing_mode="stretch_width",
                align="center",
                css_classes=["action-row", "sidebar-actions"],
                stylesheets=[panel_stylesheet(stylesheet)],
            )
        )
        confirmations.append(settings_save_controls.confirmation)
    confirmation_carrier = pn.Column(
        *confirmations,
        width=0,
        height=0,
        sizing_mode="fixed",
        margin=0,
        styles={"overflow": "hidden"},
        css_classes=["save-confirmation-carrier"],
    )
    save_group = pn.Column(
        *save_rows,
        confirmation_carrier,
        sizing_mode="stretch_width",
        margin=0,
    )

    def collapsible_section(title: str, content: Any) -> tuple[Any, BokehButton]:
        toggle = BokehButton(
            label="▼",
            width=20,
            height=20,
            button_type="primary",
            align="center",
            margin=0,
            stylesheets=[toggle_stylesheet if toggle_stylesheet is not None else stylesheet],
        )
        header = section_header(title, toggle=toggle, css_classes=["history-section-header"])

        def toggle_content() -> None:
            content.visible = not content.visible
            toggle.label = "▼" if content.visible else "▶"

        toggle.on_click(lambda _event: toggle_content())
        return pn.Column(header, content, sizing_mode="stretch_width"), toggle

    filter_section, filter_toggle = collapsible_section("Filter applications", filter_entries)
    plot_section, plot_toggle = collapsible_section("Plot applications", plot_entries)
    drawer = pn.Column(
        pn.pane.Markdown("## History", margin=(5, 8)),
        filter_section,
        separator(),
        plot_section,
        pn.Spacer(height=10, sizing_mode="fixed", css_classes=["history-save-top-space"]),
        save_group,
        pn.Spacer(height=20, sizing_mode="fixed", css_classes=["history-bottom-space"]),
        width=250,
        sizing_mode="stretch_height",
        css_classes=["history-drawer", "sidebar-content"],
        stylesheets=[stylesheet],
    )
    session = HistorySession(
        history,
        drawer,
        filter_entries,
        plot_entries,
        filter_toggle,
        plot_toggle,
        restore_entry,
        can_mutate,
        record_action,
        save_session,
        {},
        {},
    )
    return session
