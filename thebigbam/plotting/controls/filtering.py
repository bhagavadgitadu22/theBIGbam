"""Dynamic filtering section presentation state."""

from __future__ import annotations

from typing import Any, Callable

import panel as pn
from bokeh.models.widgets import Select

from ..shared.styles import panel_stylesheet


class FilterSectionController:
    """Own OR sections, row connectors, and their Panel containers."""

    def __init__(
        self,
        create_row: Callable[[dict[str, Any]], dict[str, Any]],
        on_change: Callable[[], None],
        stylesheet: Any,
        add_stylesheet: Any,
    ) -> None:
        self._create_row = create_row
        self._on_change = on_change
        self._stylesheet = stylesheet
        self.sections: list[dict[str, Any]] = []
        self.inter_section_selects: list[Any] = []
        self.content = pn.Column(sizing_mode="stretch_width")
        self._global_add = pn.widgets.Button(
            name="+ Add AND/OR",
            margin=(10, 0, 5, 0),
            button_type="primary",
            stylesheets=[panel_stylesheet(add_stylesheet)],
        )
        self._global_add.on_click(self._add_section)
        self._tail_widget: Any = self._global_add
        self.sections.append(self.create_section())
        self.rebuild_content()

    def count_rows(self) -> int:
        return sum(len(section["rows"]) for section in self.sections)

    def rebuild_section(self, section: dict[str, Any]) -> None:
        children = []
        rows = section["rows"]
        for index, row in enumerate(rows):
            row["and_div"] = row.get("and_div")
            children.append(row["row_wrapper"])
            connector_items = [row["minus_btn"]]
            if index < len(rows) - 1:
                following = rows[index + 1]
                connector = following["and_div"]
                if connector is None:
                    connector = Select(options=["AND", "OR"], value="AND", height=30, margin=(2, 0, 2, 0))
                    connector.on_change("value", lambda attr, old, new: self._on_change())
                    following["and_div"] = connector
                connector_items.append(connector)
            else:
                connector_items.append(section["add_and_btn"])
            children.append(pn.Row(*connector_items, sizing_mode="stretch_width", height=30, margin=(3, 0, 6, 0)))
        section["column"].objects = children

    def create_section(self) -> dict[str, Any]:
        section = {
            "column": pn.Column(
                sizing_mode="stretch_width",
                styles={"border-left": "3px solid #00b17c", "padding-left": "10px", "margin-left": "5px"},
            ),
            "rows": [],
            "add_and_btn": pn.widgets.Button(
                name="+ Add AND/OR",
                height=30,
                margin=(2, 0, 2, 0),
                button_type="success",
                stylesheets=[panel_stylesheet(self._stylesheet)],
            ),
        }

        def add_row(event: Any) -> None:
            section["rows"].append(self._create_row(section))
            self.rebuild_section(section)
            self._on_change()

        section["add_and_btn"].on_click(add_row)
        section["rows"].append(self._create_row(section))
        self.rebuild_section(section)
        return section

    def rebuild_content(self) -> None:
        old_values = [widget.value for widget in self.inter_section_selects]
        self.inter_section_selects.clear()
        children = []
        for index, section in enumerate(self.sections):
            if index:
                value = old_values[index - 1] if index - 1 < len(old_values) else "AND"
                connector = Select(options=["AND", "OR"], value=value, margin=(2, 0, 2, 0))
                connector.on_change("value", lambda attr, old, new: self._on_change())
                self.inter_section_selects.append(connector)
                children.append(connector)
            children.append(section["column"])
        children.append(self._tail_widget)
        self.content.objects = children

    def replace_sections(self, sections: list[dict[str, Any]]) -> None:
        self.sections[:] = sections
        self.rebuild_content()

    def _add_section(self, event: Any) -> None:
        self.sections.append(self.create_section())
        self.rebuild_content()
        self._on_change()
