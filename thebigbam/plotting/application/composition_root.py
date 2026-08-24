"""Panel/Bokeh plotting application construction."""

import os
import time
from dataclasses import dataclass
from typing import Any, Callable

import duckdb
import panel as pn
from bokeh.io import curdoc
from bokeh.layouts import row
from bokeh.models import Div, InlineStyleSheet, Tooltip
from bokeh.models.callbacks import CustomJS
from bokeh.models.widgets import (
    Button,
    HelpButton,
    RadioButtonGroup,
)

from ..controls.base import build_controls
from ..controls.filter_panel import build_filter_panel
from ..controls.genome import build_genome_controls
from ..controls.plot_parameters import build_plot_parameter_controls
from ..controls.variables import attach_module_synchronization, build_variable_controls
from ..models.preload import PreloadedPlotData
from ..models.session import PlotSessionContext
from ..repositories.color_templates import ColorTemplateRepository
from ..repositories.filter_metadata import FilterMetadataRepository
from ..repositories.filtering import FilteringRepository
from ..repositories.genome_controls import GenomeControlRepository
from ..repositories.preload import PreloadRepository
from ..services.filter_metadata import FilterMetadataService
from ..services.filter_query import FilterQueryBuilder
from ..services.filtering import FilterExpressionService, FilteringAvailabilityService
from ..shared.data_cache import SessionDataCache

# Import the plotting function from the repo
from ..shared.diagnostics import PlotDiagnostics
from ..shared.lifecycle import PlotLifecycle
from ..shared.paths import static_directory
from ..shared.timing import BrowserTimingRelay, estimate_grid_data_size, rss_mb, start_rss_watchdog
from .availability import AvailabilityBindings, AvailabilityController
from .availability_facade import AvailabilityFacade
from .layout import LayoutParts, separator
from .output_controls import build_output_controls
from .scope_transition import ScopeTransitionBindings, ScopeTransitionController
from .session_assembly import finalize_session, inspectable_application
from .session_callbacks import SessionCallbacks
from .settings_session import build_settings_session
from .state import PlotController
from .subject import SubjectBindings, SubjectController

_STATIC_DIR = static_directory()
_CSS_CACHE = {}


def _get_cached_stylesheet(filename):
    if filename not in _CSS_CACHE:
        with open(os.path.join(_STATIC_DIR, filename)) as f:
            _CSS_CACHE[filename] = f.read()
    return InlineStyleSheet(css=_CSS_CACHE[filename])


_get_rss_mb = rss_mb
_estimate_grid_data_size = estimate_grid_data_size

_TIMING = None
_current_op: str = "idle"


def _set_current_operation(operation: str) -> None:
    global _current_op
    _current_op = operation


@dataclass
class PlotApplication:
    layout: Any
    context: PlotSessionContext
    controller: PlotController
    widgets: dict[str, Any]
    diagnostics: PlotDiagnostics
    apply: Callable[[], None]


def _start_rss_watchdog(interval: int = 5) -> None:
    start_rss_watchdog(lambda: _current_op, interval)


def preload_db_data(db_path, enable_timing=False):
    """Run all expensive DB queries once at startup. Returns a dict of pure data."""
    started = time.perf_counter()
    repository = PreloadRepository(db_path)
    preloaded = repository.load()
    if enable_timing:
        elapsed = time.perf_counter() - started
        print(
            f"[timing] Preload repository: {elapsed:.3f}s "
            f"({repository.query_count} queries, {_get_rss_mb():.0f} MB RSS){_TIMING.tag(elapsed)}",
            flush=True,
        )
    return preloaded


def create_layout(
    db_path, preloaded: PreloadedPlotData, enable_timing=False, initial_settings=None, _return_application=False
):
    """Create and return the application layout for Panel serve."""
    global _current_op
    _current_op = "session_init"
    diagnostics = PlotDiagnostics(enabled=enable_timing)
    data_cache = SessionDataCache()
    if enable_timing:
        print(f"[timing] RSS at session init start: {_get_rss_mb():.0f} MB", flush=True)
        _TIMING.start_phase("Session init")

    ### Event functions
    ## Helper function to create collapsible section toggle callbacks
    def make_toggle_callback(btn, content):
        def callback():
            content.visible = not content.visible
            if content.visible:
                btn.label = "▼"
            else:
                btn.label = "▶"

        return callback

    _contig_to_mag = {}  # populated after widgets is built (see below)
    filter_projection_ref = {}

    def get_filtering_filtered_pairs():
        return filter_projection_ref["projection"].evaluate()

    availability_facade = AvailabilityFacade()
    update_widget_completions = availability_facade.update_widget

    def push_search_completions(widget, completions):
        """Update completions for a server-side search response.

        Unlike update_widget_completions, always force-fires the 'options'
        watcher (via param.trigger) even when the resulting list happens to
        be identical to what's already loaded — the frontend's pending
        Tom Select load() callback is resolved from that event, and Param's
        onlychanged watcher would otherwise never fire for a no-op update,
        leaving the search box stuck in its loading state.
        """
        availability_facade.push_search(widget, completions)

    _compute_contig_completions = availability_facade.compute_contigs
    refresh_contig_options_unlocked = availability_facade.refresh_contigs
    _compute_sample_completions = availability_facade.compute_samples
    _sort_sample_completions_for = availability_facade.sort_samples_for_mag
    refresh_sample_options_unlocked = availability_facade.refresh_samples
    _compute_mag_completions = availability_facade.compute_mags
    refresh_mag_options_unlocked = availability_facade.refresh_mags
    update_section_titles = availability_facade.update_titles

    session_callbacks = SessionCallbacks(lambda callback: curdoc().add_next_tick_callback(callback))
    apply_clicked = session_callbacks.apply_clicked
    _do_apply = session_callbacks.do_apply
    peruse_clicked = session_callbacks.show_summary

    # Track current plot state for x-range preservation across APPLY clicks
    plot_lifecycle = PlotLifecycle()
    current_plot_state = plot_lifecycle.state

    timing_relay = BrowserTimingRelay(enable_timing)
    _send_timing_ping = timing_relay.send

    def _report_download_timing(label, elapsed):
        if enable_timing:
            print(f"[timing] {label}: {elapsed:.3f}s{_TIMING.tag(elapsed)}", flush=True)

    ### Creating all DOM elements
    # Open DuckDB database connection to build widgets depending on data
    if enable_timing:
        t_init = time.perf_counter()
    conn = duckdb.connect(db_path, read_only=True)
    filtering_repository = FilteringRepository(conn)
    availability_service = FilteringAvailabilityService(filtering_repository)
    widgets = build_controls(preloaded)

    _sync_from_to_for_selected_contig = session_callbacks.sync_selected_contig_position

    if widgets["has_mags"]:
        for _mag, _ctgs in widgets["mag_to_contigs"].items():
            for _c in _ctgs:
                _contig_to_mag[_c] = _mag
    if enable_timing:
        _step = time.perf_counter() - t_init
        print(f"[timing] Session: widget creation: {_step:.3f}s{_TIMING.tag(_step)}", flush=True)

    _subplot_to_varnames = preloaded["subplot_to_varnames"]
    _encoding_by_feature = preloaded.get("encoding_by_feature")
    _total_coverage_count = preloaded.get("total_coverage_count", 0)

    if enable_timing:
        t_ui = time.perf_counter()
        t_section = time.perf_counter()

    # Load cached CSS stylesheets
    stylesheet = _get_cached_stylesheet("bokeh_styles.css")
    pink_buttons_stylesheet = _get_cached_stylesheet("pink_buttons.css")
    widgets["view_radio"].stylesheets = [pink_buttons_stylesheet]
    toggle_stylesheet = _get_cached_stylesheet("toggle_styles.css")
    grey_buttons_stylesheet = _get_cached_stylesheet("grey_buttons.css")

    # Create main elements
    ## Views section
    logo_url_local = "/assets/LOGO.png"
    logo_url_remote = "https://raw.githubusercontent.com/bhagavadgitadu22/theBIGbam/master/thebigbam/static/LOGO.png"
    logo = Div(
        text=f"""<img src="{logo_url_local}" onerror="this.onerror=function(){{this.style.display='none'}}; this.src='{logo_url_remote}'" style="width:100%; max-width:800px; padding: 0 25%;">"""
    )
    views = RadioButtonGroup(
        labels=["ONE SAMPLE", "ALL SAMPLES"],
        active=0,
        sizing_mode="stretch_width",
        stylesheets=[pink_buttons_stylesheet],
    )
    if enable_timing:
        views.js_on_change(
            "active",
            CustomJS(code="window.__thebigbam_view_change_started = performance.now();"),
        )

    # Global lock for toggles when enforcing "All samples" view (single-variable mode)
    global_toggle_lock = {"locked": False}

    ## Build Filtering section (dynamic query builder with AND/OR logic)
    filtering_toggle_btn = Button(
        label="▼", width=20, height=20, button_type="primary", align="center", margin=0, stylesheets=[toggle_stylesheet]
    )
    filtering_title = Div(text="<b>Filtering</b>", align="center")
    _filter_help_tooltip = Tooltip(
        content=(
            "Filtering rows are independent from each other. "
            'For example, "Number of samples" is based on the total number of presences '
            "of a genome in the database, and is not recomputed based on the other filtering rows"
        ),
        position="right",
    )
    _filter_help_btn = HelpButton(
        tooltip=_filter_help_tooltip,
        width=20,
        height=20,
        align="center",
        button_type="light",
        stylesheets=[toggle_stylesheet],
    )
    filtering_header = row(
        filtering_toggle_btn, filtering_title, _filter_help_btn, sizing_mode="stretch_width", align="center"
    )

    filtering_metadata = preloaded["filtering_metadata"]
    filter_metadata_service = FilterMetadataService(
        FilterMetadataRepository(db_path, filtering_metadata, enable_timing=enable_timing)
    )

    filter_ui = {}

    def mark_filters_dirty():
        """Mark draft widget changes without querying availability."""
        projection = filter_projection_ref.get("projection")
        if projection is None or not filter_ui:
            return
        dirty = projection.has_pending_changes()
        filter_ui["apply"].disabled = not dirty
        filter_ui["apply"].button_type = "primary" if dirty else "default"

    def apply_filter_changes():
        """Commit draft filters and atomically refresh availability choices."""
        global _current_op
        _current_op = "filter_change"
        if enable_timing:
            print(f"[timing] RSS at filter_change start: {_get_rss_mb():.0f} MB", flush=True)
        doc = curdoc()
        doc.hold("combine")
        filter_projection_ref["projection"].apply()
        availability_service.invalidate()
        data_cache.invalidate("filter_change")
        filter_projection_ref["projection"].invalidate()
        global_toggle_lock["locked"] = True
        try:
            refresh_contig_options_unlocked()
            refresh_sample_options_unlocked()
            refresh_mag_options_unlocked()
            update_section_titles()
        finally:
            global_toggle_lock["locked"] = False
            doc.unhold()

    filter_panel = build_filter_panel(
        preloaded=preloaded,
        widgets=widgets,
        expression_service=FilterExpressionService(
            filtering_repository,
            FilterQueryBuilder(filtering_metadata, preloaded.filter_encode, has_samples=widgets["has_samples"]),
            has_mags=widgets["has_mags"],
        ),
        metadata_service=filter_metadata_service,
        refresh=mark_filters_dirty,
        make_toggle_callback=make_toggle_callback,
        stylesheet=stylesheet,
        toggle_stylesheet=toggle_stylesheet,
        button_stylesheet=stylesheet,
        grey_button_stylesheet=grey_buttons_stylesheet,
        enable_timing=enable_timing,
        set_operation=_set_current_operation,
        header=filtering_header,
        toggle=filtering_toggle_btn,
    )
    filtering_header = filter_panel.header
    filtering_controller = filter_panel.controller
    filter_projection_ref["projection"] = filter_panel.projection
    parameter_options = filter_panel.options
    create_query_row = filter_panel.create_query_row

    filter_apply_button = Button(
        label="APPLY",
        button_type="default",
        disabled=True,
        width=120,
        align="center",
        stylesheets=[
            stylesheet,
            pink_buttons_stylesheet,
            InlineStyleSheet(
                css="""
                :host .bk-btn:disabled {
                    background: #d0d0d0 !important;
                    border-color: #b0b0b0 !important;
                    color: #666 !important;
                }
                """
            ),
        ],
        css_classes=["apply-btn"],
        margin=(8, 0, 3, 0),
    )
    filter_ui.update({"apply": filter_apply_button})

    def _set_availability_controls_disabled(disabled):
        for key in ("contig_select", "sample_select", "mag_select"):
            control = widgets.get(key)
            if control is not None and hasattr(control, "disabled"):
                control.disabled = disabled

    def _apply_filters_clicked():
        if not filter_panel.projection.has_pending_changes():
            mark_filters_dirty()
            return
        filter_apply_button.disabled = True
        _set_availability_controls_disabled(True)

        def _run():
            try:
                apply_filter_changes()
            except Exception:
                filter_apply_button.disabled = False
                filter_apply_button.button_type = "primary"
                raise
            finally:
                _set_availability_controls_disabled(False)
            filter_apply_button.button_type = "default"

        curdoc().add_next_tick_callback(_run)

    filter_apply_button.on_click(lambda _event: _apply_filters_clicked())
    filter_apply_row = pn.Row(
        pn.Spacer(sizing_mode="stretch_width"),
        filter_apply_button,
        pn.Spacer(sizing_mode="stretch_width"),
        sizing_mode="stretch_width",
        margin=0,
    )
    filtering_content = pn.Column(
        filter_panel.content,
        filter_apply_row,
        sizing_mode="stretch_width",
    )

    if enable_timing:
        _step = time.perf_counter() - t_section
        print(f"[timing]   Filtering section: {_step:.3f}s{_TIMING.tag(_step)}", flush=True)
        t_section = time.perf_counter()

    ## Build Sample section
    sample_title = Div(text="<b>Samples</b>")

    ## Build Contig section
    # Keep original full lists so we can restore when filters are off
    orig_contigs = list(widgets["contigs"])
    orig_samples = list(widgets["samples"])

    def _get_filtered_samples_for_contig(contig):
        contig_id = widgets["contig_name_to_id"].get(contig)
        if contig_id is None:
            return []
        result = list(filtering_repository.samples_for_contig_unbounded(contig_id))
        filtering = get_filtering_filtered_pairs()
        if filtering is not None:
            allowed = set(
                filtering_repository.sample_names_for_filter(
                    filtering["sql"],
                    filtering["params"],
                )
            )
            result = [sample for sample in result if sample in allowed]
        return result

    def _get_filtered_samples_for_mag(mag_name):
        sample_ids = widgets["mag_to_sample_ids"].get(mag_name, set())
        result = [sample for sample in orig_samples if widgets["sample_name_to_id"].get(sample) in sample_ids]
        filtering = get_filtering_filtered_pairs()
        if filtering is not None:
            allowed = set(
                filtering_repository.filtered_samples(
                    filtering["sql"],
                    filtering["params"],
                    mag_name=mag_name,
                )
            )
            result = [sample for sample in result if sample in allowed]
        return result

    # MAGs section header (visible only in MAG-mode databases)
    mag_title = Div(
        text=f"<span style='font-size: 1.2em;'><b>MAGs</b></span> ({len(widgets['mags'])} available)",
        align="center",
        visible=widgets["has_mags"],
    )
    mag_header = row(
        mag_title, sizing_mode="stretch_width", align="center", margin=(0, 0, 0, 0), visible=widgets["has_mags"]
    )
    separator_mags = Div(
        text="",
        height=2,
        sizing_mode="stretch_width",
        styles={"background-color": "#333", "margin-top": "10px", "margin-bottom": "10px"},
        visible=widgets["has_mags"],
    )

    # Create Contigs section header — no outer collapse; the inner
    # "Genomic annotations" and "Other genomic features" subsections each
    # have their own toggle below.
    contig_title = Div(text="<b>Contigs</b>", align="center")
    contig_header = row(contig_title, sizing_mode="stretch_width", align="center", margin=(0, 0, 0, 0))

    if enable_timing:
        _step = time.perf_counter() - t_section
        print(f"[timing]   Sample/Contig/MAG sections: {_step:.3f}s{_TIMING.tag(_step)}", flush=True)
        t_section = time.perf_counter()

    variable_controls = build_variable_controls(widgets, toggle_stylesheet, make_toggle_callback)
    variables_section_one = variable_controls.one_sample
    variables_section_all = variable_controls.all_samples
    genome_cbg_one = variable_controls.genome_checkbox
    genome_index_one = variable_controls.genome_index

    genome_controls = build_genome_controls(
        metadata_service=filter_metadata_service,
        color_templates=ColorTemplateRepository(conn).load(),
        genome_capabilities=GenomeControlRepository(conn).capabilities(),
        widgets=widgets,
        filtering_metadata=filtering_metadata,
        genome_checkbox=genome_cbg_one,
        genome_index=genome_index_one,
        stylesheet=stylesheet,
        toggle_stylesheet=toggle_stylesheet,
        make_toggle_callback=make_toggle_callback,
        enable_timing=enable_timing,
        interaction_lock=global_toggle_lock,
    )
    below_contig_content = genome_controls.content
    combined_features_cbg = genome_controls.combined_features
    from_position_input = genome_controls.from_position
    to_position_input = genome_controls.to_position

    attach_module_synchronization(widgets, global_toggle_lock)

    if enable_timing:
        _step = time.perf_counter() - t_section
        print(f"[timing]   Variables + Genome/Annotations sections: {_step:.3f}s{_TIMING.tag(_step)}", flush=True)
        t_section = time.perf_counter()

    plot_parameters = build_plot_parameter_controls(
        widgets=widgets,
        sample_scope=views,
        original_samples=orig_samples,
        parameter_options=parameter_options,
        toggle_stylesheet=toggle_stylesheet,
        stylesheet=stylesheet,
        make_toggle_callback=make_toggle_callback,
        compute_sample_completions=_compute_sample_completions,
        push_search_completions=push_search_completions,
        sort_sample_completions=_sort_sample_completions_for,
        interaction_lock=global_toggle_lock,
    )
    separator_plotting_params = plot_parameters.separator
    plotting_params_header = plot_parameters.header
    plotting_params_content = plot_parameters.content
    mag_params_category_select = plot_parameters.mag_category
    mag_params_sort_sample_select = plot_parameters.mag_sort_sample
    mag_params_sort_sample_row = plot_parameters.mag_sort_sample_row
    _mag_sort_category_sources = plot_parameters.mag_category_sources
    sample_params_header = plot_parameters.sample_header
    sample_params_content = plot_parameters.sample_content
    sample_order_category_select = plot_parameters.sample_order_category
    _sample_sort_category_sources = plot_parameters.sample_category_sources
    _sample_contig_categories = plot_parameters.sample_contig_categories
    _sample_mag_categories = plot_parameters.sample_mag_categories
    _sample_sort_current_categories = plot_parameters.sample_current_categories

    availability_controller = AvailabilityController(
        AvailabilityBindings(
            availability_service=availability_service,
            filtering_pairs=get_filtering_filtered_pairs,
            original_contigs=orig_contigs,
            original_samples=orig_samples,
            sample_scope=views,
            widgets=widgets,
            sort_sample_select=mag_params_sort_sample_select,
            update_completions=update_widget_completions,
            total_coverage_count=_total_coverage_count,
            filtering_title=filtering_title,
            contig_title=contig_title,
            sample_title=sample_title,
            mag_title=mag_title,
        )
    )
    availability_facade.attach(availability_controller)

    subject_controller = SubjectController(
        SubjectBindings(
            widgets=widgets,
            interaction_lock=global_toggle_lock,
            compute_contigs=_compute_contig_completions,
            compute_samples=_compute_sample_completions,
            compute_mags=_compute_mag_completions,
            push_completions=push_search_completions,
            refresh_contigs=refresh_contig_options_unlocked,
            refresh_samples=refresh_sample_options_unlocked,
            refresh_mags=refresh_mag_options_unlocked,
            update_titles=update_section_titles,
            from_position=from_position_input,
            to_position=to_position_input,
            sample_order_category=sample_order_category_select,
            sample_contig_categories=_sample_contig_categories,
            sample_mag_categories=_sample_mag_categories,
            sample_current_categories=_sample_sort_current_categories,
        )
    )
    session_callbacks.attach_subject(subject_controller)
    subject_controller.attach()

    ## Create final Apply and Peruse data buttons
    apply_button = Button(
        label="APPLY",
        align="center",
        button_type="primary",
        stylesheets=[pink_buttons_stylesheet],
        css_classes=["apply-btn"],
        margin=(5, 0, 0, 0),
    )
    apply_button.on_click(apply_clicked)

    # Peruse button will be positioned in the plot area, styled to match toolbar
    output_controls = build_output_controls(
        db_path=db_path,
        connection=conn,
        widgets=widgets,
        sample_scope=views,
        filtered_samples=_get_filtered_samples_for_contig,
        combined_features=combined_features_cbg,
        subplot_variables=_subplot_to_varnames,
        from_position=from_position_input,
        to_position=to_position_input,
        stylesheet=pink_buttons_stylesheet,
        enable_timing=enable_timing,
        timing=_TIMING,
        report_timing=_report_download_timing if enable_timing else None,
        show_summary=peruse_clicked,
    )
    peruse_button = output_controls.peruse_button

    # Hidden carrier to send summary HTML to client browser via JS
    summary_carrier = output_controls.summary_carrier
    session_callbacks.attach_summary(output_controls.summary_controller)

    command_hint_pane = output_controls.command_hint
    download_widgets = output_controls.download_widgets

    settings_session = build_settings_session(
        db_path=db_path,
        widgets=widgets,
        sample_scope=views,
        genome=genome_controls,
        parameters=plot_parameters,
        filtering=filtering_controller,
        filtering_metadata=filtering_metadata,
        create_query_row=create_query_row,
        apply_button=apply_button,
        stylesheet=stylesheet,
    )
    buttons_row = settings_session.buttons_row

    # sample_section must exist before settings restoration runs, since restoring the
    # ALL SAMPLES scope fires on_view_change, which sets sample_section.visible.
    separator_samples = separator()
    sample_section = pn.Column(
        separator_samples,
        sample_title,
        widgets["sample_select"],
        sizing_mode="stretch_width",
        margin=0,
    )

    def refresh_scope_availability():
        # Preserve the selected subject when an old scope snapshot contains a
        # sample which is no longer valid under the active filter.
        refresh_sample_options_unlocked()
        refresh_mag_options_unlocked()
        refresh_contig_options_unlocked()
        update_section_titles()

    scope_transition = ScopeTransitionController(
        ScopeTransitionBindings(
            widgets=widgets,
            sample_scope=views,
            sample_section=sample_section,
            variables_one=variables_section_one,
            variables_all=variables_section_all,
            sample_parameters_header=sample_params_header,
            sample_parameters_content=sample_params_content,
            mag_sort_sample_row=mag_params_sort_sample_row,
            mag_sort_category=mag_params_category_select,
            interaction_lock=global_toggle_lock,
            diagnostics=diagnostics,
            enable_timing=enable_timing,
            timing=_TIMING,
            send_timing_ping=_send_timing_ping,
            set_operation=_set_current_operation,
            update_completions=update_widget_completions,
            refresh_availability=refresh_scope_availability,
        )
    )
    session_callbacks.attach_scope(scope_transition)
    controller = scope_transition.plot_controller
    for variable_group in widgets["variables_widgets_all"]:
        variable_group.on_change("active", scope_transition.variable_callback(variable_group))
    views.on_change("active", scope_transition.view_changed)

    if initial_settings:
        settings_session.restore(initial_settings)
        apply_filter_changes()

    ## Initialize section titles with counts
    update_section_titles()

    if enable_timing:
        _step = time.perf_counter() - t_section
        print(f"[timing]   Plotting params + layout assembly: {_step:.3f}s{_TIMING.tag(_step)}", flush=True)

    peruse_button.visible = False  # Initially hidden
    finalized = finalize_session(
        parts=LayoutParts(
            logo=logo,
            sample_scope=views,
            filtering_header=filtering_header,
            filtering_content=filtering_content,
            mag_separator=separator_mags,
            mag_header=mag_header,
            view_radio=widgets["view_radio"],
            mag_select=widgets["mag_select"],
            contig_header=contig_header,
            contig_select=widgets["contig_select"],
            below_contig_content=below_contig_content,
            sample_section=sample_section,
            variables_one=variables_section_one,
            variables_all=variables_section_all,
            plotting_separator=separator_plotting_params,
            plotting_header=plotting_params_header,
            plotting_content=plotting_params_content,
            buttons=buttons_row,
        ),
        has_samples=widgets["has_samples"],
        summary_carrier=summary_carrier,
        stylesheet=stylesheet,
        timing_models=(timing_relay.ping, timing_relay.ack) if enable_timing else (),
        apply_arguments=dict(
            db_path=db_path,
            connection=conn,
            widgets=widgets,
            sample_scope=views,
            genome=genome_controls,
            parameters=plot_parameters,
            filtered_samples=_get_filtered_samples_for_contig,
            filtered_mag_samples=_get_filtered_samples_for_mag,
            encoding_by_feature=_encoding_by_feature,
            send_timing_ping=_send_timing_ping,
            sync_position=_sync_from_to_for_selected_contig,
            command_hint=command_hint_pane,
            plot_state=current_plot_state,
            diagnostics=diagnostics,
            data_cache=data_cache,
            buttons={**download_widgets, "peruse": peruse_button},
            plot_lifecycle=plot_lifecycle,
            enable_timing=enable_timing,
            timing=_TIMING,
            set_operation=_set_current_operation,
        ),
    )
    session_callbacks.attach_placeholder(finalized.placeholder)
    session_callbacks.attach_apply(finalized.apply_controller)
    layout = finalized.layout

    if enable_timing:
        _step = time.perf_counter() - t_ui
        print(f"[timing] Session: UI construction: {_step:.3f}s{_TIMING.tag(0)}", flush=True)
        _step = time.perf_counter() - t_init
        print(f"[timing] Session ready (total: {_step:.3f}s){_TIMING.tag(0)}", flush=True)
        _TIMING.summary("Session init")

    if _return_application:
        return inspectable_application(
            PlotApplication,
            session=finalized,
            connection=conn,
            preloaded=preloaded,
            widgets={**widgets, "sample_scope": views},
            diagnostics=diagnostics,
            plot_state=current_plot_state,
            plot_controller=controller,
            apply=_do_apply,
        )
    return layout


def create_application(db_path, preloaded, enable_timing=False, initial_settings=None):
    """Build an inspectable application object for tests and diagnostics."""
    return create_layout(
        db_path,
        preloaded,
        enable_timing=enable_timing,
        initial_settings=initial_settings,
        _return_application=True,
    )
