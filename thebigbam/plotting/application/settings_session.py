"""Session-level assembly for settings collection, saving, and restoration."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import panel as pn

from ..settings.collection import SettingsCollector
from ..settings.controls import SettingsSaveControls
from ..settings.restoration import restore_settings
from ..shared.styles import panel_stylesheet
from .wiring import make_settings_collector_bindings, make_settings_restore_bindings


@dataclass(frozen=True)
class SettingsSession:
    buttons_row: Any
    confirmation: Any
    collector: SettingsCollector
    restore_bindings: Any
    filter_projection: Any
    save_controls: Any

    def restore(self, settings: dict[str, Any]) -> None:
        restore_settings(settings, self.restore_bindings)

    def collect_applied(self) -> dict[str, Any]:
        return self.collector.collect_applied(self.filter_projection.applied_expression())


def build_settings_session(
    *,
    db_path: str,
    widgets: dict[str, Any],
    sample_scope: Any,
    genome: Any,
    parameters: Any,
    filtering: Any,
    filtering_metadata: Any,
    create_query_row: Any,
    apply_button: Any,
    filter_projection: Any,
    stylesheet: Any,
    record_action: Any | None = None,
) -> SettingsSession:
    collector = SettingsCollector(
        make_settings_collector_bindings(db_path, widgets, sample_scope, genome, parameters, filtering)
    )
    # SAVE SETTINGS preserves the complete workspace, including filter edits
    # which have not yet been committed with APPLY FILTERS. Action history uses
    # collect_applied() below because it represents work that actually ran.
    save_controls = SettingsSaveControls(db_path, collector.collect, stylesheet, record_action)
    buttons = pn.Row(
        apply_button,
        align="center",
        css_classes=["action-row", "sidebar-actions"],
        stylesheets=[panel_stylesheet(stylesheet)],
    )
    restore_bindings = make_settings_restore_bindings(
        widgets=widgets,
        sample_scope=sample_scope,
        genome=genome,
        parameters=parameters,
        filtering=filtering,
        filtering_metadata=filtering_metadata,
        create_query_row=create_query_row,
    )
    return SettingsSession(buttons, save_controls.confirmation, collector, restore_bindings, filter_projection, save_controls)
