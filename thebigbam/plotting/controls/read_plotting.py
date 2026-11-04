"""Read variable toggles and the reserved read-filter section."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from bokeh.models import CheckboxButtonGroup, CheckboxGroup, Div

from ..shared.styles import section_separator


@dataclass(frozen=True)
class ReadPlotControls:
    section: Any
    add_reads: CheckboxButtonGroup
    all_reads: CheckboxButtonGroup
    max_reads: Any
    link_pairs: CheckboxGroup
    plot_height: Any
    read_height: Any


def build_read_plot_controls(
    *,
    max_reads: Any,
    link_pairs: CheckboxGroup,
    plot_height: Any,
    read_height: Any,
) -> ReadPlotControls:
    """Create the variable toggle and reserve a section for read filters."""
    import panel as pn

    add_reads = CheckboxButtonGroup(labels=["Read plot"], active=[], sizing_mode="stretch_width")
    all_reads = CheckboxButtonGroup(labels=["Read plot"], active=[], sizing_mode="stretch_width", visible=False)
    section = pn.Column(
        section_separator(),
        Div(text="<span style='font-size: 1.2em;'><b>Read plotting</b></span>"),
        add_reads,
        all_reads,
        sizing_mode="stretch_width",
        margin=0,
    )

    return ReadPlotControls(section, add_reads, all_reads, max_reads, link_pairs, plot_height, read_height)
