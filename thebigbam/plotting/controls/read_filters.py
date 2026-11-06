"""Read filter blocks reuse Filtering's controller, expression and serialization."""

from bokeh.models import Select

from ..models.read_filters import BOOLEAN_READ_FILTERS, READ_FILTER_FIELDS, SINGLE_VALUE_READ_FILTERS
from ..settings.persistence import serialize_filter_sections
from ..shared.styles import ACTION_CONTROL_MARGIN
from .filter_projection import expression_from_widgets
from .filtering import FilterSectionController
from .read_filter_rows import ReadFilterRowFactory


class ReadFilterControls:
    def __init__(self, stylesheet="", circular_available=False):
        fields = [field for field in READ_FILTER_FIELDS if field != "Circular" or circular_available]
        self.factory = ReadFilterRowFactory(fields, stylesheet)
        self.controller = FilterSectionController(self.factory.create_row, lambda: None, stylesheet, stylesheet)
        self.factory.controller = self.controller
        self.content = self.controller.content

    def expression(self):
        return expression_from_widgets(self.controller.sections, self.controller.inter_section_selects)

    def serialize(self):
        return serialize_filter_sections(self.controller.sections, self.controller.inter_section_selects)

    def restore(self, saved):
        sections, connectors = [], []
        for block in saved or []:
            section = self.controller.create_section()
            section["rows"].clear()
            for item in block.get("rows", []):
                column = item.get("column")
                if column == "Alignment percentage":
                    column = "Read % aligned"
                if column not in self.factory.fields:
                    print(f"[settings] Warning: skipped unavailable read filter '{column}'.", flush=True)
                    continue
                if column in SINGLE_VALUE_READ_FILTERS and item.get("value") in (None, ""):
                    continue
                row = self.factory.create_row(section, column, item.get("operator"))
                value = item.get("value")
                if value is not None:
                    if column in BOOLEAN_READ_FILTERS and isinstance(value, bool):
                        value = "TRUE" if value else "FALSE"
                    row["input_ref"]["widget"].value = value
                if section["rows"]:
                    row["and_div"] = Select(options=["AND", "OR"], value=item.get("row_and_or") or "AND",
                                             height=30, margin=ACTION_CONTROL_MARGIN)
                section["rows"].append(row)
            if section["rows"]:
                self.controller.rebuild_section(section)
                sections.append(section)
                connectors.append(block.get("section_and_or") or "AND")
        self.controller.replace_sections(sections or [self.controller.create_section()])
        for widget, value in zip(self.controller.inter_section_selects, connectors[1:]):
            widget.value = value
