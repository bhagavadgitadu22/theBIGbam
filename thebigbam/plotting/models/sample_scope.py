"""Pure ONE/ALL sample-scope decisions and a minimal UI projection adapter."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable

from ..application.state import SampleScope


@dataclass(frozen=True)
class SampleScopeVisibility:
    """Visibility state for controls affected by a sample-scope transition."""

    sample_section: bool
    variables_one: bool
    variables_all: bool
    sample_params_header: bool
    sample_params_content: bool
    mag_sort_sample_row: bool


def visibility_for_sample_scope(
    scope: SampleScope,
    *,
    mag_sort_category: str,
) -> SampleScopeVisibility:
    """Return UI visibility without importing Bokeh, Panel, or DuckDB."""
    is_all = SampleScope(scope) is SampleScope.ALL
    return SampleScopeVisibility(
        sample_section=not is_all,
        variables_one=not is_all,
        variables_all=is_all,
        sample_params_header=is_all,
        sample_params_content=False,
        mag_sort_sample_row=is_all and mag_sort_category != "Contig",
    )


def apply_changed_attributes(assignments: Iterable[tuple[Any, str, Any]]) -> int:
    """Assign only changed attributes and return the number of mutations."""
    changed = 0
    for target, attribute, value in assignments:
        if getattr(target, attribute) != value:
            setattr(target, attribute, value)
            changed += 1
    return changed
