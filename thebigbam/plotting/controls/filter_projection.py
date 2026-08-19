"""Translate filtering widgets into domain expressions and evaluated results."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from typing import Any

from ..models.filters import FilterExpression, FilterPredicate, FilterSection
from ..services.filtering import FilterExpressionService


def expression_from_widgets(
    sections: Sequence[dict[str, Any]], inter_section_selects: Sequence[Any]
) -> FilterExpression:
    """Project the mutable filtering UI into an immutable domain expression."""
    projected_sections: list[FilterSection] = []
    section_connectors: list[str] = []
    for section_index, section_data in enumerate(sections):
        predicates: list[FilterPredicate] = []
        connectors: list[str] = []
        for row_data in section_data["rows"]:
            input_ref = row_data["input_ref"]
            value = input_ref["widget"].value
            if value is None or value == "":
                continue
            if input_ref["is_panel"] and isinstance(value, str) and not value.strip():
                continue
            if predicates:
                connector = row_data.get("and_div")
                connectors.append("OR" if connector is not None and connector.value == "OR" else "AND")
            predicates.append(
                FilterPredicate(
                    row_data["category_select"].value,
                    row_data["subcategory_select"].value,
                    row_data["comparison_select"].value,
                    value,
                )
            )
        if not predicates:
            continue
        if projected_sections:
            connector_index = section_index - 1
            connector = (
                inter_section_selects[connector_index].value if connector_index < len(inter_section_selects) else "AND"
            )
            section_connectors.append("OR" if connector == "OR" else "AND")
        projected_sections.append(FilterSection(tuple(predicates), tuple(connectors)))
    return FilterExpression(tuple(projected_sections), tuple(section_connectors))


class FilterWidgetProjection:
    """Evaluate the current widget expression and expose the legacy result shape."""

    def __init__(
        self,
        service: FilterExpressionService,
        sections: Sequence[dict[str, Any]],
        inter_section_selects: Sequence[Any],
        *,
        enable_timing: bool = False,
        set_operation: Callable[[str], None] = lambda _operation: None,
    ) -> None:
        self.service = service
        self.sections = sections
        self.inter_section_selects = inter_section_selects
        self.enable_timing = enable_timing
        self.set_operation = set_operation
        self._valid = False
        self._result: dict[str, Any] | None = None

    def invalidate(self) -> None:
        self._valid = False
        self.service.invalidate()

    def evaluate(self) -> dict[str, Any] | None:
        if self._valid:
            return self._result
        self.set_operation("filtering_pairs")
        try:
            result = self.service.evaluate(expression_from_widgets(self.sections, self.inter_section_selects))
            if result is None:
                self._result = None
            else:
                self._result = {
                    "sql": result.compiled.sql,
                    "params": list(result.compiled.parameters),
                    "count_pairs": result.pair_count,
                    "count_contigs": result.contig_count,
                    "count_samples": result.sample_count,
                    "count_mag_pairs": result.mag_pair_count,
                }
                if self.enable_timing:
                    print(
                        f"[timing] Filter SQL built; count_pairs={result.pair_count} "
                        f"contigs={result.contig_count} samples={result.sample_count}",
                        flush=True,
                    )
            self._valid = True
            return self._result
        finally:
            self.set_operation("idle")
