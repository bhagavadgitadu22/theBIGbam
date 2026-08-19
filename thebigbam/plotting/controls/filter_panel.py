"""Construction of the complete filtering panel and its projection boundary."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable

from bokeh.layouts import row
from bokeh.models import Div, Tooltip
from bokeh.models.widgets import Button, HelpButton

from ..repositories.filtering import FilteringRepository
from ..services.filter_query import FilterQueryBuilder
from ..services.filtering import FilterExpressionService
from .filter_projection import FilterWidgetProjection
from .filter_rows import FilterRowFactory
from .filter_visualizations import FilterVisualizations
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


def build_filter_panel(
    *,
    db_path: str,
    preloaded: Any,
    widgets: dict[str, Any],
    repository: FilteringRepository,
    refresh: Callable[[], None],
    make_toggle_callback: Callable[[Any, Any], Callable[[], None]],
    stylesheet: Any,
    toggle_stylesheet: Any,
    button_stylesheet: Any,
    grey_button_stylesheet: Any,
    enable_timing: bool,
    set_operation: Callable[[str], None],
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
    expression_service = FilterExpressionService(
        repository,
        FilterQueryBuilder(metadata, preloaded.filter_encode, has_samples=widgets["has_samples"]),
        has_mags=widgets["has_mags"],
    )
    visualizations = FilterVisualizations(
        db_path,
        metadata,
        enable_timing,
        grey_button_stylesheet,
        refresh,
        preloaded.filter_encode,
    )
    row_factory = FilterRowFactory(
        db_path=db_path,
        filtering_metadata=metadata,
        columns=options.columns,
        raw_columns=options.columns_raw,
        visualizations=visualizations,
        refresh=refresh,
        stylesheet=stylesheet,
        enable_timing=enable_timing,
        filter_encode=preloaded.filter_encode,
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
    return FilterPanel(header, controller.content, controller, projection, options, row_factory.create_row)
