"""Construction of genome annotation, sequence, and position controls."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Mapping

import panel as pn
from bokeh.layouts import row
from bokeh.models import Div, InlineStyleSheet
from bokeh.models.widgets import (
    Button,
    CheckboxButtonGroup,
    CheckboxGroup,
    HelpButton,
    MultiChoice,
    TextInput,
)

from .color_rules import build_color_rule_controls


@dataclass(frozen=True)
class GenomeControls:
    content: Any
    combined_features: Any | None
    feature_types: Any | None
    plot_isoforms: Any | None
    feature_label: Any | None
    sequence: Any | None
    translated_sequence: Any | None
    genome_master: Any | None
    from_position: Any
    to_position: Any
    apply_annotation_rules: Any | None
    color_qualifier_options: list[str]
    custom_color_rows: list[dict[str, Any]]
    mag_color_rows: list[dict[str, Any]]
    create_color_row: Any
    rebuild_color_rows: Any
    rebuild_mag_color_rows: Any


def _build_collapsible(title, children, toggle_stylesheet, make_toggle_callback, *, margin=(0, 0, 5, 0)):
    toggle = Button(
        label="▶",
        width=20,
        height=20,
        button_type="primary",
        align="center",
        margin=0,
        stylesheets=[toggle_stylesheet],
    )
    toggle.styles = {"padding": "0px", "line-height": "20px"}
    header = row(toggle, Div(text=title, align="center"), sizing_mode="stretch_width", align="center")
    content = pn.Column(*children, visible=False, sizing_mode="stretch_width", margin=margin)
    toggle.on_click(make_toggle_callback(toggle, content))
    return header, content


def _wire_master_checkbox(master, features, interaction_lock) -> None:
    lock = {"locked": False}

    def master_changed(attr, old, new):
        if lock["locked"] or interaction_lock.get("locked", False):
            return
        lock["locked"] = True
        try:
            features.active = list(range(len(features.labels))) if 0 in master.active else []
        finally:
            lock["locked"] = False

    def features_changed(attr, old, new):
        if lock["locked"] or interaction_lock.get("locked", False):
            return
        lock["locked"] = True
        try:
            master.active = [0] if features.labels and len(features.active) == len(features.labels) else []
        finally:
            lock["locked"] = False

    master.on_change("active", master_changed)
    features.on_change("active", features_changed)


def build_genome_controls(
    metadata_service: Any,
    color_templates: Mapping[str, list[dict[str, Any]]],
    genome_capabilities: Any,
    widgets: Mapping[str, Any],
    filtering_metadata: Mapping[str, Any],
    genome_checkbox: Any | None,
    genome_index: int | None,
    stylesheet: Any,
    toggle_stylesheet: Any,
    make_toggle_callback: Callable[[Any, Any], Callable[..., None]],
    enable_timing: bool,
    interaction_lock: Mapping[str, bool],
) -> GenomeControls:
    genome_cbg_one = genome_checkbox
    genome_index_one = genome_index
    ## Build Genome module controls (placed in Contigs section, shared between views)
    genome_section = None
    combined_features_cbg = None

    # Feature type filter (MultiChoice) - only show if annotation types exist
    # Gene map is plotted if at least one feature type is selected
    feature_type_multichoice = None
    if widgets["annotation_types"]:
        # Initially select only CDS if available, otherwise nothing
        initial_value = ["CDS"] if "CDS" in widgets["annotation_types"] else []
        multichoice_stylesheet = InlineStyleSheet(css=":host { background-color: white; }")
        feature_type_multichoice = MultiChoice(
            options=list(widgets["annotation_types"]),
            value=initial_value,
            placeholder="Choose feature types to plot",
            sizing_mode="stretch_width",
            stylesheets=[multichoice_stylesheet],
        )

    color_controls = build_color_rule_controls(
        metadata_service,
        filtering_metadata,
        color_templates,
        stylesheet,
        enable_timing,
    )
    template_select = color_controls.template_select
    color_qualifier_options = color_controls.qualifier_options
    custom_color_rows = color_controls.custom_rows
    custom_color_column = color_controls.custom_column
    mag_track_color_rows = color_controls.mag_rows
    mag_track_color_column = color_controls.mag_column
    feature_label_select = color_controls.feature_label_select
    create_color_row = color_controls.create_row
    rebuild_color_rows = color_controls.rebuild_custom
    rebuild_mag_track_color_rows = color_controls.rebuild_mag
    # Plot isoforms checkbox - only show if at least one locus_tag appears more than once
    plot_isoforms_cbg = None
    if genome_capabilities.has_isoforms:
        plot_isoforms_cbg = CheckboxGroup(
            labels=["Plot isoforms"],
            active=[],  # Unchecked by default
        )

    # Build combined labels: Genome features (without Gene map) + Custom contig features
    combined_labels = []
    if genome_cbg_one is not None:
        combined_labels.extend(genome_cbg_one.labels)  # Already without Gene map
    if widgets["custom_contig_subplots"]:
        combined_labels.extend(widgets["custom_contig_subplots"])

    if combined_labels:
        combined_features_cbg = CheckboxButtonGroup(
            labels=combined_labels, active=[], sizing_mode="stretch_width", orientation="vertical"
        )

    # Add Genome section to contig_content
    below_contig_children = []

    # Create position range inputs
    from_position_input = TextInput(
        value="1", placeholder="Start position", sizing_mode="stretch_width", margin=(0, 0, 0, 0)
    )
    to_position_input = TextInput(
        value="", placeholder="End position", sizing_mode="stretch_width", margin=(0, 0, 0, 0)
    )

    position_label_from = Div(text="From", width=40, margin=(5, 0, 5, 5))
    position_label_to = Div(text="to", width=25, margin=(5, 0, 5, 5))

    # Create Reset button to reset position inputs
    position_reset_button = Button(label="Reset", stylesheets=[stylesheet], margin=(0, 5, 0, 5))

    def reset_position_inputs():
        from_position_input.value = "1"
        if widgets["has_mags"] and widgets["view_radio"].active == 0:
            # MAG view: clear contig selection and show full MAG extent
            widgets["contig_select"].value = ""
            selected_mag = widgets["mag_select"].value
            if selected_mag:
                total = sum(
                    widgets["contig_lengths"].get(c, 0) for c in widgets["mag_to_contigs"].get(selected_mag, [])
                )
                to_position_input.value = str(total)
            else:
                to_position_input.value = ""
        elif widgets["contig_select"].value and widgets["contig_select"].value in widgets["contig_lengths"]:
            to_position_input.value = str(widgets["contig_lengths"][widgets["contig_select"].value])
        else:
            to_position_input.value = ""

    position_reset_button.on_click(lambda event: reset_position_inputs())

    position_row = row(
        position_label_from,
        from_position_input,
        position_label_to,
        to_position_input,
        position_reset_button,
        sizing_mode="stretch_width",
        margin=(10, 0, 5, 0),
    )

    below_contig_children.append(position_row)

    sequence_cbg = None
    sequence_row = None
    if genome_capabilities.has_sequence:
        sequence_cbg = CheckboxGroup(labels=["Plot sequence"], active=[])
        sequence_row = row(sequence_cbg, sizing_mode="stretch_width")

    translated_sequence_cbg = None
    translated_sequence_row = None
    if genome_capabilities.has_translation:
        translated_sequence_cbg = CheckboxGroup(labels=["Plot translated sequence"], active=[])
        translated_sequence_row = row(translated_sequence_cbg, sizing_mode="stretch_width")

    genome_master_cbg = None
    apply_annotation_rules_cbg = None
    if feature_type_multichoice is not None or combined_features_cbg is not None:
        # --- Subsection 1: "Genomic annotations to plot" (collapsible) ---
        # Covers the feature type multichoice, isoform toggle, and both
        # sequence/translated-sequence rows. Coloring and labelling widgets
        # live in the separate "Customise genomic annotations" subsection.
        annotations_children = []
        if feature_type_multichoice is not None:
            annotations_children.append(feature_type_multichoice)
        if plot_isoforms_cbg is not None:
            annotations_children.append(plot_isoforms_cbg)
        if sequence_row is not None:
            annotations_children.append(sequence_row)
        if translated_sequence_row is not None:
            annotations_children.append(translated_sequence_row)
        annotations_header, annotations_content = _build_collapsible(
            "Genomic annotations to plot",
            annotations_children,
            toggle_stylesheet,
            make_toggle_callback,
            margin=(0, 0, 0, 0),
        )

        # --- Subsection 2: "Customise genomic annotations" (collapsible) ---
        # Independent sibling of the annotations-to-plot section. Bundles
        # the coloring rules (template + custom color rows) and the label
        # dropdown. Starts collapsed.
        customise_header = None
        customise_content = None
        has_customise = template_select is not None or color_qualifier_options or feature_label_select is not None
        if has_customise:
            customise_children = []
            if feature_label_select is not None:
                customise_children.append(feature_label_select)
            if template_select is not None:
                customise_children.append(template_select)
            if template_select is not None or color_qualifier_options:
                customise_children.append(Div(text="Color annotations with:"))
            if color_qualifier_options:
                customise_children.append(custom_color_column)
            customise_header, customise_content = _build_collapsible(
                "Customise genomic annotations",
                customise_children,
                toggle_stylesheet,
                make_toggle_callback,
            )

        # --- MAG track coloring section (only in MAG databases with annotation columns) ---
        mag_track_section_header = None
        mag_track_section_content = None
        apply_annotation_rules_cbg = None
        if widgets["has_mags"] and color_qualifier_options:
            apply_annotation_rules_cbg = CheckboxGroup(
                labels=["Apply genomic annotations coloring rules"],
                active=[],
                margin=(4, 0, 4, 0),
            )
            mag_track_section_header, mag_track_section_content = _build_collapsible(
                "Customise MAG track",
                [Div(text="Color annotations on MAG track with:"), mag_track_color_column, apply_annotation_rules_cbg],
                toggle_stylesheet,
                make_toggle_callback,
            )

        # --- Subsection 3: "Other genomic features to plot" (collapsible) ---
        other_features_toggle_btn = Button(
            label="▶",
            width=20,
            height=20,
            button_type="primary",
            align="center",
            margin=0,
            stylesheets=[toggle_stylesheet],
        )
        other_features_toggle_btn.styles = {"padding": "0px", "line-height": "20px"}
        # Master checkbox whose label doubles as the section title — same
        # pattern Variables modules use (CheckboxGroup with [module_name]).
        # Toggling it checks/unchecks every feature in combined_features_cbg
        # at once; bidirectional sync is wired below.
        genome_master_cbg = CheckboxGroup(labels=["Other genomic features to plot"], active=[])

        genome_help_tooltip = widgets["helps_widgets"][genome_index_one] if genome_index_one is not None else None
        if genome_help_tooltip is not None:
            help_btn = HelpButton(
                tooltip=genome_help_tooltip,
                width=20,
                height=20,
                align="center",
                button_type="light",
                stylesheets=[toggle_stylesheet],
            )
            other_features_header = row(
                other_features_toggle_btn,
                genome_master_cbg,
                help_btn,
                sizing_mode="stretch_width",
                align="center",
                styles={"overflow": "hidden"},
            )
        else:
            other_features_header = row(
                other_features_toggle_btn, genome_master_cbg, sizing_mode="stretch_width", align="center"
            )

        other_features_content = None
        if combined_features_cbg is not None:
            combined_features_cbg.visible = True
            other_features_content = pn.Column(
                combined_features_cbg, visible=False, sizing_mode="stretch_width", margin=(0, 0, 0, 0)
            )
            other_features_toggle_btn.on_click(make_toggle_callback(other_features_toggle_btn, other_features_content))

            # Bidirectional sync: master ⇄ individual feature toggles.
            # Mirrors make_module_callback / make_variable_callback used by
            # the Variables subsections at ~line 2025.
            _wire_master_checkbox(genome_master_cbg, combined_features_cbg, interaction_lock)

        # Assemble the full genome section from the independent collapsibles.
        genome_section_children = [annotations_header, annotations_content]
        if mag_track_section_header is not None:
            genome_section_children.append(mag_track_section_header)
            genome_section_children.append(mag_track_section_content)
        if customise_header is not None:
            genome_section_children.append(customise_header)
            if customise_content is not None:
                genome_section_children.append(customise_content)
        genome_section_children.append(other_features_header)
        if other_features_content is not None:
            genome_section_children.append(other_features_content)
        genome_section = pn.Column(
            *genome_section_children, visible=True, sizing_mode="stretch_width", margin=(0, 0, 0, 0)
        )

    if genome_section is not None:
        below_contig_children = list(below_contig_children) + [genome_section]

    # Fallback: if no genome section exists, append sequence/translated rows directly
    if genome_section is None and sequence_row is not None:
        below_contig_children.append(sequence_row)
    if genome_section is None and translated_sequence_row is not None:
        below_contig_children.append(translated_sequence_row)

    below_contig_content = pn.Column(
        *below_contig_children, visible=True, sizing_mode="stretch_width", margin=(0, 0, 0, 0)
    )

    # Initialize position inputs if contig is pre-filled
    if widgets["contig_select"].value and widgets["contig_select"].value in widgets["contig_lengths"]:
        to_position_input.value = str(widgets["contig_lengths"][widgets["contig_select"].value])

    return GenomeControls(
        content=below_contig_content,
        combined_features=combined_features_cbg,
        feature_types=feature_type_multichoice,
        plot_isoforms=plot_isoforms_cbg,
        feature_label=feature_label_select,
        sequence=sequence_cbg,
        translated_sequence=translated_sequence_cbg,
        genome_master=genome_master_cbg,
        from_position=from_position_input,
        to_position=to_position_input,
        apply_annotation_rules=apply_annotation_rules_cbg,
        color_qualifier_options=color_qualifier_options,
        custom_color_rows=custom_color_rows,
        mag_color_rows=mag_track_color_rows,
        create_color_row=create_color_row,
        rebuild_color_rows=rebuild_color_rows,
        rebuild_mag_color_rows=rebuild_mag_track_color_rows,
    )
