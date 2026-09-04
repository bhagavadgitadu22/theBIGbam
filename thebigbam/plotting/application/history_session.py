"""Panel projection and restore hooks for per-browser-session history."""

from __future__ import annotations

import html
import os
from dataclasses import dataclass
from typing import Any, Callable, Iterable

import panel as pn
from bokeh.models import Button as BokehButton
from bokeh.models import CustomJS, Div

from ..settings.controls import save_confirmation_js
from ..settings.history import HistoryEntry, SessionHistory
from ..settings.history_descriptions import (
    HistoryDescriptionContext,
    HistoryDescriptionLine,
    canonical_history_description_lines,
    diff_description_lines,
)
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
    filter_short_descriptions: Any
    plot_short_descriptions: Any
    restore_entry: Callable[[HistoryEntry], None]
    can_mutate: Callable[[], bool]
    record_action: Callable[[str, dict[str, Any]], None]
    save: Callable[..., Any]
    rows_by_id: dict[str, Any]
    actions_by_id: dict[str, str]
    apply_steps_by_id: dict[str, int]
    description_lines_by_id: dict[str, tuple[HistoryDescriptionLine, ...]]
    description_context: HistoryDescriptionContext

    def append(
        self, action: str, settings: dict[str, Any], *, apply_step: int | None = None
    ) -> HistoryEntry:
        previous_ids = {item.id for item in self.history.entries}
        entry = self.history.append(action, settings)
        self.description_lines_by_id[entry.id] = canonical_history_description_lines(
            entry, self.description_context
        )
        retained_ids = {item.id for item in self.history.entries}
        expired = previous_ids - retained_ids
        for expired_id in expired:
            self._remove_row(expired_id, refresh=False)
        row = self._entry_row(entry)
        self.rows_by_id[entry.id] = row
        self.actions_by_id[entry.id] = action
        if apply_step is not None:
            self.apply_steps_by_id[entry.id] = apply_step
        last_index = len(self.history.for_action(action)) - 1
        changed = {last_index}
        if expired:
            changed.add(0)
        self._refresh_action(action, changed)
        return entry

    def restore(self, entries: Iterable[HistoryEntry]) -> None:
        """Install saved entries without recording them as new actions."""
        self.history.restore(entries)
        self.rows_by_id.clear()
        self.actions_by_id.clear()
        self.apply_steps_by_id.clear()
        self.description_lines_by_id.clear()
        for entry in self.history.entries:
            self.description_lines_by_id[entry.id] = canonical_history_description_lines(
                entry, self.description_context
            )
            row = self._entry_row(entry)
            self.rows_by_id[entry.id] = row
            self.actions_by_id[entry.id] = entry.action
        self._refresh_action("apply_filters")
        self._refresh_action("apply_plot")

    @staticmethod
    def _empty_row() -> Any:
        return pn.pane.Markdown("_No entries yet._", margin=(0, 5))

    def _content_for(self, action: str) -> Any:
        return self.filter_entries if action == "apply_filters" else self.plot_entries

    def _entry_row(self, entry: HistoryEntry) -> Any:
        description = self._description_cell(entry, None)
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
            description,
            restore,
            remove,
            css_classes=["history-entry"],
            sizing_mode="stretch_width",
            margin=0,
        )

    def _short_descriptions(self, action: str) -> bool:
        checkbox = (
            self.filter_short_descriptions
            if action == "apply_filters"
            else self.plot_short_descriptions
        )
        return bool(checkbox.value)

    def _description_lines(self, entry: HistoryEntry) -> tuple[HistoryDescriptionLine, ...]:
        lines = self.description_lines_by_id.get(entry.id)
        if lines is None:
            lines = canonical_history_description_lines(entry, self.description_context)
            self.description_lines_by_id[entry.id] = lines
        return lines

    def _description_cell(self, entry: HistoryEntry, previous: HistoryEntry | None) -> Any:
        current_lines = self._description_lines(entry)
        lines = (
            diff_description_lines(
                self._description_lines(previous) if previous is not None else (),
                current_lines,
            )
            if self._short_descriptions(entry.action)
            else tuple(line for line in current_lines if not line.default)
        )
        rendered_lines = []
        for line in lines:
            escaped = html.escape(line.text)
            displayed = f"<s>{escaped}</s>" if line.removed else escaped
            rendered_lines.append(
                f'<div class="history-description-line" title="{html.escape(line.text, quote=True)}">'
                f"{displayed}</div>"
            )
        return pn.pane.HTML(
            "".join(rendered_lines),
            sizing_mode="stretch_width",
            margin=0,
            css_classes=["history-description"],
            stylesheets=list(self.drawer.stylesheets),
        )

    def _refresh_action(self, action: str, changed_indices: set[int] | None = None) -> None:
        entries = list(self.history.for_action(action))
        previous = None
        rows = []
        for entry in entries:
            row = self.rows_by_id.get(entry.id)
            if row is None:
                row = self._entry_row(entry)
                self.rows_by_id[entry.id] = row
                self.actions_by_id[entry.id] = action
            if changed_indices is None or len(rows) in changed_indices:
                row.objects = [self._description_cell(entry, previous), *row.objects[1:]]
            rows.append(row)
            previous = entry
        self._content_for(action).objects = rows or [self._empty_row()]

    def _remove(self, entry_id: str, *, record: bool = True) -> None:
        if not self.can_mutate():
            return
        entry = next((item for item in self.history.entries if item.id == entry_id), None)
        action_entries = list(self.history.for_action(entry.action)) if entry is not None else []
        removed_index = action_entries.index(entry) if entry is not None else -1
        apply_step = self.apply_steps_by_id.get(entry_id)
        if entry is not None and self.history.remove(entry_id):
            self._remove_row(entry_id, refresh=False)
            changed = {removed_index} if removed_index < len(self.history.for_action(entry.action)) else set()
            self._refresh_action(entry.action, changed)
            if record:
                self.record_action(
                    "remove_history",
                    {
                        "history_sequence": entry.sequence,
                        "history_action": entry.action,
                        "apply_step": apply_step,
                    },
                )

    def _restore(self, entry: HistoryEntry) -> None:
        if not self.can_mutate():
            return
        self.record_action(
            "restore_history",
            {
                "history_sequence": entry.sequence,
                "history_action": entry.action,
                "apply_step": self.apply_steps_by_id.get(entry.id),
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

    def _remove_row(self, entry_id: str, *, refresh: bool = True) -> None:
        row = self.rows_by_id.pop(entry_id, None)
        action = self.actions_by_id.pop(entry_id, None)
        self.apply_steps_by_id.pop(entry_id, None)
        self.description_lines_by_id.pop(entry_id, None)
        if row is None or action is None:
            return
        if refresh:
            self._refresh_action(action)

def build_history_session(
    *,
    db_path: str,
    restore_entry: Callable[[HistoryEntry], None],
    can_mutate: Callable[[], bool] = lambda: True,
    record_action: Callable[[str, dict[str, Any]], None] = lambda _action, _details: None,
    stylesheet: Any,
    toggle_stylesheet: Any | None = None,
    settings_save_controls: Any | None = None,
    initial_entries: Iterable[HistoryEntry] = (),
    description_context: HistoryDescriptionContext | None = None,
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
    filter_short_descriptions = pn.widgets.Checkbox(
        name="Use git-like short descriptions",
        value=True,
        margin=(0, 5, 4, 5),
        css_classes=["history-description-mode"],
    )
    plot_short_descriptions = pn.widgets.Checkbox(
        name="Use git-like short descriptions",
        value=True,
        margin=(0, 5, 4, 5),
        css_classes=["history-description-mode"],
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

    def collapsible_section(title: str, mode: Any, content: Any) -> tuple[Any, BokehButton]:
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
            mode.visible = content.visible
            toggle.label = "▼" if content.visible else "▶"

        toggle.on_click(lambda _event: toggle_content())
        return pn.Column(header, mode, content, sizing_mode="stretch_width"), toggle

    filter_section, filter_toggle = collapsible_section(
        "Filter applications", filter_short_descriptions, filter_entries
    )
    plot_section, plot_toggle = collapsible_section(
        "Plot applications", plot_short_descriptions, plot_entries
    )
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
        filter_short_descriptions,
        plot_short_descriptions,
        restore_entry,
        can_mutate,
        record_action,
        save_session,
        {},
        {},
        {},
        {},
        description_context or HistoryDescriptionContext(),
    )
    filter_short_descriptions.param.watch(
        lambda _event: session._refresh_action("apply_filters"), "value"
    )
    plot_short_descriptions.param.watch(
        lambda _event: session._refresh_action("apply_plot"), "value"
    )
    session.restore(initial_entries)
    return session
