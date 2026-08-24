"""Construction of the application's core selection controls."""

from __future__ import annotations

from bokeh.models import CheckboxButtonGroup, CheckboxGroup, RadioButtonGroup, Tooltip

from ..models.preload import PreloadedPlotData
from .searchable_select import SearchableSelect


def build_controls(preloaded: PreloadedPlotData):
    """Create core selection widgets without performing database queries."""
    mags = preloaded.mags
    contigs = preloaded.contigs
    samples = preloaded.samples
    mag_select = SearchableSelect(
        value=mags[0] if len(mags) == 1 else "",
        options=list(mags),
        placeholder="Type to search MAGs...",
        server_search=True,
        sizing_mode="stretch_width",
        margin=(0, 5, 0, 5),
        visible=preloaded.has_mags,
    )
    view_radio = RadioButtonGroup(
        labels=["MAG VIEW", "CONTIG VIEW"],
        active=1,
        visible=preloaded.has_mags,
        sizing_mode="stretch_width",
        margin=(0, 5, 10, 5),
    )
    contig_select = SearchableSelect(
        value=contigs[0] if len(contigs) == 1 else "",
        options=list(contigs),
        placeholder="Type to search contigs...",
        server_search=True,
        sizing_mode="stretch_width",
        margin=(0, 5, 0, 5),
    )
    sample_select = SearchableSelect(
        value=samples[0] if len(samples) == 1 else "",
        options=list(samples),
        placeholder="Type to search samples...",
        server_search=True,
        sizing_mode="stretch_width",
        margin=(0, 5, 0, 5),
    )
    module_widgets_one = []
    variables_widgets_one = []
    variables_widgets_all = []
    helps_widgets = []
    for module, variables, help_text in zip(
        preloaded.module_names,
        preloaded.module_variables,
        preloaded.module_helps,
    ):
        module_widgets_one.append(CheckboxGroup(labels=[module], active=[]))
        variables_widgets_one.append(
            CheckboxButtonGroup(labels=list(variables), active=[], sizing_mode="stretch_width", orientation="vertical")
        )
        variables_widgets_all.append(
            CheckboxButtonGroup(labels=list(variables), active=[], sizing_mode="stretch_width", orientation="vertical")
        )
        helps_widgets.append(Tooltip(content=help_text, position="right") if help_text else None)
    return {
        "sample_select": sample_select,
        "contig_select": contig_select,
        "mag_select": mag_select,
        "view_radio": view_radio,
        "contig_name_to_id": preloaded.contig_name_to_id,
        "contig_id_to_name": preloaded.contig_id_to_name,
        "sample_name_to_id": preloaded.sample_name_to_id,
        "sample_id_to_name": preloaded.sample_id_to_name,
        "mag_to_contigs": preloaded.mag_to_contigs,
        "contig_to_mag": preloaded.contig_to_mag,
        "mag_to_sample_ids": preloaded.mag_to_sample_ids,
        "mag_to_contig_offsets": preloaded.mag_to_contig_offsets,
        "mags": mags,
        "has_mags": preloaded.has_mags,
        "module_names": preloaded.module_names,
        "module_widgets_one": module_widgets_one,
        "helps_widgets": helps_widgets,
        "variables_widgets_one": variables_widgets_one,
        "variables_widgets_all": variables_widgets_all,
        "contigs": contigs,
        "contig_lengths": preloaded.contig_lengths,
        "samples": samples,
        "custom_contig_subplots": preloaded.custom_contig_subplots,
        "annotation_types": preloaded.annotation_types,
        "has_samples": preloaded.has_samples,
    }
