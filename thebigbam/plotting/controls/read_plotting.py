"""Read variable toggles and the shared AND/OR read-filter editor."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from bokeh.models import CheckboxButtonGroup, CheckboxGroup, Div

from ..shared.styles import section_separator
from .collapsible import build_collapsible
from .read_filters import ReadFilterControls


@dataclass(frozen=True)
class ReadPlotControls:
    section: Any
    add_reads: CheckboxButtonGroup
    all_reads: CheckboxButtonGroup
    max_reads: Any
    link_pairs: CheckboxGroup
    plot_height: Any
    read_height: Any
    filters: ReadFilterControls | None = None
    max_discrepancy_window: Any = None
    filter_header: Any = None
    filter_content: Any = None


def build_read_plot_controls(
    *,
    max_reads: Any,
    link_pairs: CheckboxGroup,
    plot_height: Any,
    read_height: Any,
    stylesheet: Any = "",
    circular_available: bool = False,
    max_discrepancy_window: Any = None,
    toggle_stylesheet: Any = "",
    make_toggle_callback: Any = None,
) -> ReadPlotControls:
    """Create the read toggles and alignment filtering controls."""
    import panel as pn

    add_reads = CheckboxButtonGroup(labels=["Read plot"], active=[], sizing_mode="stretch_width")
    all_reads = CheckboxButtonGroup(labels=["Read plot"], active=[], sizing_mode="stretch_width", visible=False)
    filters = ReadFilterControls(stylesheet, circular_available)
    filter_title = Div(text='<span title="Read length includes soft/hard clips. Alignment length counts aligned query bases (M/I/=/X). '
                 'Alignment (%) is 100 × exact matching bases / full read length, including clips. '
                 'Unknown percentages without sufficient CIGAR/MD/NM information do not match numeric comparisons. '
                 'Coordinates are 1-based and local to the contig, including in MAG view. Alignment includes uses the alignment span. '
                 'Alignment end is the BAM end coordinate, which may exceed contig length for circular reads. '
                 'Insert sizes use absolute SAM TLEN; unknown sizes do not match numeric filters. '
                 'Non-inward pair tests the orientation of mapped, same-contig mates, not the proper-pair flag. '
                 'Circular means an alignment crosses the origin in a circular-mapping sample.">Filter reads:</span>')
    filter_header, filter_content = build_collapsible(
        filter_title, [filters.content], toggle_stylesheet, make_toggle_callback, slug="filter-reads",
    )
    section = pn.Column(
        section_separator(),
        Div(text="<span style='font-size: 1.2em;'><b>Read plotting</b></span>"),
        add_reads,
        all_reads,
        filter_header,
        filter_content,
        sizing_mode="stretch_width",
        margin=0,
    )

    return ReadPlotControls(section, add_reads, all_reads, max_reads, link_pairs, plot_height, read_height,
                            filters, max_discrepancy_window, filter_header, filter_content)
