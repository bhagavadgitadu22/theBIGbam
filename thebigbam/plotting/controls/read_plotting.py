"""Read plot selection, filtering, and discrepancy display controls."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from bokeh.layouts import row
from bokeh.models import CheckboxButtonGroup, CheckboxGroup, Div, InlineStyleSheet
from bokeh.models.widgets import MultiChoice

from ..shared.styles import FIELD_LABEL_MARGIN, section_separator
from .collapsible import build_collapsible
from .read_color_rules import ReadColorRuleControls
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
    discrepancies: Any = None
    customise_header: Any = None
    customise_content: Any = None
    color_rules: Any = None


def build_read_plot_controls(
    *,
    max_reads: Any,
    plot_height: Any,
    read_height: Any,
    stylesheet: Any = "",
    circular_available: bool = False,
    max_discrepancy_window: Any = None,
    strand_colors: Any = None,
    toggle_stylesheet: Any = "",
    make_toggle_callback: Any = None,
) -> ReadPlotControls:
    """Create the read toggles and alignment filtering controls."""
    import panel as pn

    add_reads = CheckboxButtonGroup(labels=["Read plot"], active=[], sizing_mode="stretch_width")
    all_reads = CheckboxButtonGroup(labels=["Read plot"], active=[], sizing_mode="stretch_width", visible=False)
    link_pairs = CheckboxGroup(labels=["Link paired reads"], active=[], sizing_mode="stretch_width")
    filters = ReadFilterControls(stylesheet, circular_available)
    color_rules = ReadColorRuleControls(
        stylesheet, circular_available=circular_available, strand_colors=strand_colors,
    )
    filter_title = Div(text='<span title="These rules select alignments to include in the plot. The initial Alignment class = Primary rule shows primary alignments; add OR rules to include secondary or supplementary alignments. Read length includes soft/hard clips. Alignment length counts aligned query bases (M/I/=/X). '
                 'Read % aligned is 100 × exact matching bases / full read length, including clips. '
                 'Unknown percentages without sufficient CIGAR/MD/NM information do not match numeric comparisons. '
                 'Coordinates are 1-based and local to the contig, including in MAG view. Alignment includes uses the alignment span. '
                 'Alignment end is the BAM end coordinate, which may exceed contig length for circular reads. '
                 'Alignment class identifies primary, secondary, and supplementary records. '
                 'Insert sizes use absolute SAM TLEN; unknown sizes do not match numeric filters. '
                 'Non-inward pair tests the orientation of mapped, same-contig mates, not the proper-pair flag. '
                 'Circular means an alignment crosses the origin in a circular-mapping sample.">Select alignments to show</span>')
    filter_header, filter_content = build_collapsible(
        filter_title, [filters.content], toggle_stylesheet, make_toggle_callback,
        slug="filter-reads", initially_open=True,
    )
    discrepancies_label = Div(text="Discrepancies to plot:", margin=FIELD_LABEL_MARGIN)
    discrepancies = MultiChoice(
        title="",
        options=["Clippings", "Insertions", "Deletions", "Mismatches"],
        value=[],
        placeholder="Choose discrepancies to plot",
        sizing_mode="stretch_width",
        margin=FIELD_LABEL_MARGIN,
        stylesheets=[InlineStyleSheet(css=":host { background-color: white; }")],
    )
    customise_header, customise_content = build_collapsible(
        color_rules.title,
        [color_rules.content, discrepancies_label, discrepancies,
         row(link_pairs, sizing_mode="stretch_width")],
        toggle_stylesheet,
        make_toggle_callback,
        slug="customise-reads",
    )
    section = pn.Column(
        section_separator(),
        Div(text="<span style='font-size: 1.2em;'><b>Read plotting</b></span>"),
        filter_header,
        filter_content,
        customise_header,
        customise_content,
        sizing_mode="stretch_width",
        margin=0,
    )

    return ReadPlotControls(section, add_reads, all_reads, max_reads, link_pairs, plot_height, read_height,
                            filters, max_discrepancy_window, filter_header, filter_content,
                            discrepancies, customise_header, customise_content, color_rules)
