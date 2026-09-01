"""Construction of the complete filtering panel and its projection boundary."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Mapping

from bokeh.layouts import row
from bokeh.models import Div, Tooltip
from bokeh.models.widgets import Button, HelpButton

from ..renderers.filter_distributions import FilterVisualizations
from .filter_projection import FilterWidgetProjection
from .filter_rows import FilterRowFactory
from .filtering import FilterSectionController
from .parameter_options import ParameterOptionCatalog


@dataclass(frozen=True)
class FilterPanel:
    header: Any
    content: Any
    controller: FilterSectionController
    projection: FilterWidgetProjection
    options: ParameterOptionCatalog
    create_query_row: Callable[..., Any]
    set_distribution: Callable[..., None]
    set_distribution_scale: Callable[..., None]


def build_filter_panel(
    *,
    preloaded: Any,
    widgets: dict[str, Any],
    expression_service: Any,
    metadata_service: Any,
    refresh: Callable[[], None],
    make_toggle_callback: Callable[[Any, Any], Callable[[], None]],
    stylesheet: Any,
    toggle_stylesheet: Any,
    button_stylesheet: Any,
    muted_button_stylesheet: Any,
    enable_timing: bool,
    set_operation: Callable[[str], None],
    record_action: Callable[[str, Mapping[str, Any]], Any] | None = None,
    header: Any | None = None,
    toggle: Any | None = None,
) -> FilterPanel:
    if header is None or toggle is None:
        toggle = Button(
            label="▼",
            width=20,
            height=20,
            button_type="primary",
            align="center",
            margin=0,
            stylesheets=[toggle_stylesheet],
        )
        title = Div(text="<b>Filtering</b>", align="center")
        help_button = HelpButton(
            tooltip=Tooltip(content="Filtering rows are independent from each other.", position="right"),
            width=20,
            height=20,
            align="center",
            button_type="light",
            stylesheets=[toggle_stylesheet],
        )
        header = row(toggle, title, help_button, sizing_mode="stretch_width", align="center")

    metadata = preloaded["filtering_metadata"]
    options = ParameterOptionCatalog.from_filtering_metadata(metadata)
    visualizations = FilterVisualizations(
        metadata_service,
        metadata,
        enable_timing,
        muted_button_stylesheet,
        refresh,
        preloaded.filter_encode,
        set_operation,
        record_action,
    )
    row_factory = FilterRowFactory(
        metadata_service=metadata_service,
        filtering_metadata=metadata,
        columns=options.columns,
        raw_columns=options.columns_raw,
        visualizations=visualizations,
        refresh=refresh,
        stylesheet=stylesheet,
        enable_timing=enable_timing,
        filter_encode=preloaded.filter_encode,
        record_action=record_action,
    )
    controller = FilterSectionController(row_factory.create_row, refresh, stylesheet, button_stylesheet)
    row_factory.attach_controller(controller)
    projection = FilterWidgetProjection(
        expression_service,
        controller.sections,
        controller.inter_section_selects,
        enable_timing=enable_timing,
        set_operation=set_operation,
    )
    toggle.on_click(make_toggle_callback(toggle, controller.content))
    return FilterPanel(
        header,
        controller.content,
        controller,
        projection,
        options,
        row_factory.create_row,
        row_factory.set_distribution,
        row_factory.set_distribution_scale,
    )
