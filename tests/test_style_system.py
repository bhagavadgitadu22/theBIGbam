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
    assert "margin-left: -4px" in rail_rule
    assert "margin-right: -4px" in rail_rule
    assert "z-index: 20" in rail_rule
    assert ":host(.history-entry-list)" in css
    assert "max-height: 252px" in css
    history_drawer_rule = css.split(":host(.history-drawer) {", 1)[1].split("}", 1)[0]
    assert "padding:" not in history_drawer_rule
    sidebar_rule = css.split(":host(.sidebar-content) {", 1)[1].split("}", 1)[0]
    assert "padding: var(--sidebar-gutter)" in sidebar_rule
    history_rule = css.split(":host(.history-entry) {", 1)[1].split("}", 1)[0]
    assert "border-left" not in history_rule
    assert "border: 1px solid #ccc" in history_rule
    assert "border-radius: 5px" in history_rule
    assert "padding: 5px" in history_rule
    history_time_rule = css.split(":host(.history-time) {", 1)[1].split("}", 1)[0]
    assert "height: 30px" in history_time_rule
    assert "padding: 0" in history_time_rule
    assert "width: 100%" in history_time_rule
    history_button_rule = css.split(":host(.history-time) .bk-btn {", 1)[1].split("}", 1)[0]
    assert "align-items: center" in history_button_rule
    assert "font-weight: bold" in history_button_rule


def test_control_row_constructors_apply_the_shared_spacing_contract():
    bokeh_row = bokeh_control_row(Div(), Div(), css_classes=["position-row"])
    panel_row = panel_control_row(Div(), Div(), css_classes=["history-entry"])

    assert bokeh_row.spacing == CONTROL_GAP_PX
    assert bokeh_row.css_classes == ["position-row", "control-row"]
    assert panel_row.css_classes == ["history-entry", "control-row"]


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
