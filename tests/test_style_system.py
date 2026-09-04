import re
from pathlib import Path

from bokeh.models import Button, Div, Spacer

from thebigbam.plotting.shared.styles import (
    CONTROL_GAP_PX,
    bokeh_control_row,
    panel_control_row,
    section_header,
)

STATIC = Path("thebigbam/static")


def test_sidebar_design_tokens_and_semantic_classes_share_one_stylesheet():
    css = (STATIC / "bokeh_styles.css").read_text(encoding="utf-8")

    for contract in (
        "--sidebar-gutter",
        "--control-gap",
        "--drawer-tab-width",
        "--panel-background",
        ".sidebar-content",
        ".sidebar-busy",
        ".control-row",
        ".sidebar-field",
        ".nested-section",
        ".action-row",
        ".sidebar-actions",
        ":host(.action-primary)",
        ":host(.action-add)",
        ":host(.action-muted)",
        ":host(.mode-switch)",
        ":host(.drawer-toggle-left)",
        ":host(.drawer-toggle-right)",
    ):
        assert contract in css

    assert ":root,\n:host {" in css
    assert re.search(r"var\([^)]*,", css) is None

    assert not (STATIC / "pink_buttons.css").exists()
    assert not (STATIC / "grey_buttons.css").exists()

    assert "border-radius: 0 6px 6px 0" in css
    assert "border-radius: 6px 0 0 6px" in css
    assert ":host(.drawer-toggle) .bk-btn" in css
    assert "background: var(--color-pink) !important" in css
    assert "color: white !important" in css
    page_rule = css.split("html, body {", 1)[1].split("}", 1)[0]
    assert "margin: 0" in page_rule
    assert "overflow-x: hidden" in page_rule
    main_layout_rule = css.split(".main-layout {", 1)[1].split("}", 1)[0]
    assert "overflow-x: hidden" in main_layout_rule
    assert "max-width: 100vw" in main_layout_rule
    plot_area_rule = css.split(":host(.plot-area) {", 1)[1].split("}", 1)[0]
    assert "min-width: 0" in plot_area_rule
    assert "overflow-x: hidden" in plot_area_rule
    rail_rule = css.split(":host(.panel-resize-rail) {", 1)[1].split("}", 1)[0]
    assert "width: 8px" in rail_rule
    assert "min-width: 8px" in rail_rule
    assert "z-index: 20" in rail_rule
    shell_rule = css.split(":host(.panel-shell) {", 1)[1].split("}", 1)[0]
    assert "position: relative" in shell_rule
    assert "overflow: visible" in shell_rule
    responsive_rule = css.split(":host(.responsive-control-row) {", 1)[1].split("}", 1)[0]
    assert "flex-wrap: wrap" in responsive_rule
    control_rule = css.split(":host(.control-row) {", 1)[1].split("}", 1)[0]
    assert "column-gap: var(--control-gap) !important" in control_rule
    assert "gap: var(--control-gap) !important" in control_rule
    assert ".control-row > .bk-Row" in css
    nested_control_rule = css.split(":host(.control-row) .bk-Row {", 1)[1].split("}", 1)[0]
    assert "column-gap: var(--control-gap) !important" in nested_control_rule
    assert "gap: var(--control-gap) !important" in nested_control_rule
    assert ":host(.filtering-row) {" in css
    assert ":host(.coloring-row) {" in css
    assert ":host(.filtering-row) .bk-Row" not in css
    assert ":host(.coloring-row) .bk-Row" not in css
    filtering_rule = css.split(":host(.filtering-row) {", 1)[1].split("}", 1)[0]
    coloring_rule = css.split(":host(.coloring-row) {", 1)[1].split("}", 1)[0]
    assert "minmax(0, 2.5fr) 42px" in filtering_rule
    assert "padding-right: var(--control-row-end-inset)" in filtering_rule
    assert "padding-right: var(--control-row-end-inset)" in coloring_rule
    assert ":host(.filtering-row) > *" in css
    assert ":host(.coloring-row) > *" in css
    lookup_rule = css.split(":host(.filter-lookup) {", 1)[1].split("}", 1)[0]
    lookup_button_rule = css.split(":host(.filter-lookup) .bk-btn {", 1)[1].split("}", 1)[0]
    assert "width: 42px !important" in lookup_rule
    assert "width: 100% !important" in lookup_button_rule
    assert "min-width: 0 !important" in lookup_button_rule
    filter_action_rule = css.split(":host(.filter-action-row) {", 1)[1].split("}", 1)[0]
    color_picker_rule = css.split(":host(.color-picker-field) {", 1)[1].split("}", 1)[0]
    color_input_rule = css.split(':host(.color-picker-field) input[type="color"] {', 1)[1].split("}", 1)[0]
    assert "align-items: center" in filter_action_rule
    assert "min-height: 30px" in filter_action_rule
    assert "height: 30px !important" in color_picker_rule
    assert "height: 30px !important" in color_input_rule
    assert "padding: 0 !important" in color_input_rule
    for field_class in (
        "filter-category",
        "filter-metric",
        "filter-operator",
        "filter-value",
        "color-qualifier",
        "color-operator",
        "color-value",
    ):
        field_rule = css.split(f":host(.{field_class}) {{", 1)[1].split("}", 1)[0]
        assert "flex:" in field_rule
        assert "min-width: 0" in field_rule
    nested_rule = css.split(":host(.nested-section) {", 1)[1].split("}", 1)[0]
    assert "width: calc(100% - var(--space-sm))" in nested_rule
    assert "max-width: calc(100% - var(--space-sm))" in nested_rule
    color_list_rule = css.split(":host(.color-rule-list) {", 1)[1].split("}", 1)[0]
    assert "overflow-x: hidden" in color_list_rule
    assert ":host(.history-entry-list)" in css
    assert "max-height: 500px" in css
    history_drawer_rule = css.split(":host(.history-drawer) {", 1)[1].split("}", 1)[0]
    assert "padding:" not in history_drawer_rule
    sidebar_rule = css.split(":host(.sidebar-content) {", 1)[1].split("}", 1)[0]
    assert "padding: var(--sidebar-gutter)" in sidebar_rule
    history_rule = css.split(":host(.history-entry) {", 1)[1].split("}", 1)[0]
    assert "border-left" not in history_rule
    assert "border: 1px solid #ccc" in history_rule
    assert "border-radius: 5px" in history_rule
    assert "padding: 5px" in history_rule
    description_rule = css.split(":host(.history-description) {", 1)[1].split("}", 1)[0]
    description_line_rule = css.split(":host(.history-description) .history-description-line {", 1)[1].split("}", 1)[0]
    assert "overflow-x: hidden" in description_rule
    assert "overflow-y: visible" in description_rule
    assert "text-overflow: ellipsis" in description_line_rule
    assert "white-space: nowrap" in description_line_rule


def test_control_row_constructors_apply_the_shared_spacing_contract():
    bokeh_row = bokeh_control_row(Div(), Div(), css_classes=["position-row"])
    panel_row = panel_control_row(Div(), Div(), css_classes=["history-entry"])

    assert bokeh_row.spacing == CONTROL_GAP_PX
    assert bokeh_row.css_classes == ["position-row", "control-row"]
    assert panel_row.css_classes == ["history-entry", "control-row"]
    assert panel_row.styles == {"column-gap": "4px", "gap": "4px"}


def test_collapsible_section_header_centers_title_with_balanced_toggle_slot():
    toggle = Button(label="▼", width=20, height=20)
    header = section_header("Filter applications", toggle=toggle)

    assert header.children[0] is toggle
    assert header.children[1].text == (
        "<span style='font-size: 1.2em;'><b>Filter applications</b></span>"
    )
    assert "section-title" in header.children[1].css_classes
    assert isinstance(header.children[2], Spacer)
    assert header.children[2].width == toggle.width
    assert "section-header" in header.css_classes
