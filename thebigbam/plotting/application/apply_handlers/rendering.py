"""Mode-specific rendering handlers for Apply requests."""

from __future__ import annotations

import time

from ...models.composition import (
    CompositionDisplay,
    MagCompositionRequest,
    MagOrdering,
    SingleSampleCompositionRequest,
    TrackSelection,
)
from ...shared.timing import rss_mb
from ..apply_inputs import ApplyInputs
from ..apply_modes import all_sample_features, one_sample_features, track_visibility, without_gene_map
from ..apply_pipeline import (
    PlotPresenter,
    PlotResult,
)
from ..subject import sample_order_category_for_view
from .bindings import ApplyBindings


def warm_plot_pipeline_imports() -> None:
    """Load heavy rendering pipelines before the first interactive Apply."""
    from .. import all_samples_pipeline, sample_mag_pipeline

    # Keep explicit references so static analysis and future refactors cannot
    # accidentally turn this warm-up into a no-op import.
    _ = (
        all_samples_pipeline.build_all_samples_plot,
        sample_mag_pipeline.build_mag_plot,
        sample_mag_pipeline.build_single_sample_plot,
    )


class ApplyRenderEngine:
    """Shared rendering implementation used by typed mode handlers."""

    def __init__(self, bindings: ApplyBindings, plot_presenter: PlotPresenter) -> None:
        self.bindings = bindings
        self.plot_presenter = plot_presenter

    def render_mag(self, inputs: ApplyInputs, t_apply_start: float) -> PlotResult:
        bindings = self.bindings
        _encoding_by_feature = bindings._encoding_by_feature
        _get_filtered_samples_for_mag = bindings._get_filtered_samples_for_mag
        _mag_sort_category_sources = bindings._mag_sort_category_sources
        _sample_sort_category_sources = bindings._sample_sort_category_sources
        _send_timing_ping = bindings._send_timing_ping
        _sync_from_to_for_selected_contig = bindings._sync_from_to_for_selected_contig
        combined_features_cbg = bindings.combined_features_cbg
        conn = bindings.conn
        enable_timing = bindings.enable_timing
        from_position_input = bindings.from_position_input
        mag_params_category_select = bindings.mag_params_category_select
        mag_params_direction = bindings.mag_params_direction
        mag_params_metric_select = bindings.mag_params_metric_select
        mag_params_sort_sample_select = bindings.mag_params_sort_sample_select
        max_samples_input = bindings.max_samples_input
        max_sequence_window_input = bindings.max_sequence_window_input
        sample_order_category_select = bindings.sample_order_category_select
        sample_order_direction = bindings.sample_order_direction
        sample_order_metric_select = bindings.sample_order_metric_select
        sequence_cbg = bindings.sequence_cbg
        to_position_input = bindings.to_position_input
        translated_sequence_cbg = bindings.translated_sequence_cbg
        widgets = bindings.widgets
        timing = bindings.timing
        set_operation = bindings.set_operation
        has_samples = inputs.has_samples
        is_all = inputs.is_all
        sample = inputs.sample
        selected_feature_types = inputs.selected_feature_types
        genbank_path = inputs.genbank_path
        plot_isoforms = inputs.plot_isoforms
        feature_label_key = inputs.feature_label_key
        custom_colors = inputs.custom_colors
        mag_track_colors = inputs.mag_track_colors
        max_track_dots = inputs.max_track_dots
        active_variables_widgets = inputs.active_variables_widgets
        max_genemap_window = inputs.max_genemap_window
        same_y_scale = inputs.same_y_scale
        subplot_size = inputs.subplot_size
        genemap_size = inputs.genemap_size
        sequence_size = inputs.sequence_size
        translated_sequence_size = inputs.translated_sequence_size
        max_binning = inputs.max_binning
        min_coverage_freq = inputs.min_coverage_freq

        active_mag = widgets["mag_select"].value
        from ...shared.timing import ApplyProfiler

        profiler = ApplyProfiler("mag_all" if is_all else "mag_one", active_mag, sample, enabled=enable_timing)

        # Parse position inputs (None → use full MAG extent inside plotting function)
        xstart = int(from_position_input.value) if from_position_input.value.strip() else 1
        xend = int(to_position_input.value) if to_position_input.value.strip() else None

        mag_window = (xend - xstart) if xend is not None else float("inf")

        # Sequence / translated-sequence visibility — same threshold logic
        # as the contig-view path (see below).
        mag_max_sequence_window = int(max_sequence_window_input.value)
        visibility = track_visibility(
            mag_window,
            max_genemap_window=max_genemap_window,
            max_sequence_window=mag_max_sequence_window,
            sequence_requested=sequence_cbg is not None and 0 in sequence_cbg.active,
            translation_requested=translated_sequence_cbg is not None and 0 in translated_sequence_cbg.active,
        )
        if sequence_cbg is not None and 0 in sequence_cbg.active:
            if not visibility.sequence:
                print(
                    f"Warning: Sequence will not be plotted for regions larger than {mag_max_sequence_window} bp.",
                    flush=True,
                )
        if translated_sequence_cbg is not None and 0 in translated_sequence_cbg.active:
            if not visibility.translation:
                print(
                    f"Warning: Translated sequence will not be plotted for regions larger than {mag_max_sequence_window} bp.",
                    flush=True,
                )

        # Collect requested features
        if is_all:
            # ALL SAMPLES mode: pick one sample-level variable (same logic as contig-view all-samples)
            mag_requested_features, selected_var = all_sample_features(combined_features_cbg, active_variables_widgets)
            if selected_var:
                mag_requested_features.append(selected_var)
            mag_requested_features = [f for f in mag_requested_features if f != "Gene map"]
            # Compute filtered samples for this MAG (union across all its contigs)
            mag_allowed_samples = set(_get_filtered_samples_for_mag(active_mag))
        else:
            # ONE SAMPLE mode: collect all selected features (same logic as the One-sample Contig path below)
            mag_requested_features = one_sample_features(combined_features_cbg, active_variables_widgets)
            # Gene map is handled via genbank_path, not as a feature name
            mag_requested_features = [f for f in mag_requested_features if f != "Gene map"]
            mag_allowed_samples = None

        # Read MAG contig sort parameters
        _sort_cat = mag_params_category_select.value
        _sort_metric = mag_params_metric_select.value
        _sort_ascending = mag_params_direction.active == 0
        _sort_source = _mag_sort_category_sources.get(_sort_cat)
        _sort_sample_name = None
        if _sort_cat != "Contig" and _sort_source:
            if is_all:
                _sort_sample_name = mag_params_sort_sample_select.value
                # Required selection is validated by ApplyRequestBuilder.
            else:
                _sort_sample_name = sample

        # Read sample order parameters
        _sample_sort_cat = sample_order_category_for_view(sample_order_category_select.value, mag_view=True)
        _sample_sort_metric = sample_order_metric_select.value
        _sample_sort_ascending = sample_order_direction.active == 0
        _sample_sort_source = _sample_sort_category_sources.get(_sample_sort_cat)

        # If from/to track a selected contig's cached position, pass it as
        # focus_contig so the MAG composer can derive the new
        # position from the already-computed sorted members list — no extra
        # DB query needed.
        _sel_contig = widgets["contig_select"].value
        focus_contig = None
        if _sel_contig and active_mag:
            _cached_off = widgets["mag_to_contig_offsets"].get(active_mag, {})
            if _sel_contig in _cached_off:
                _old_off = _cached_off[_sel_contig]
                _sel_clen = widgets["contig_lengths"].get(_sel_contig, 0)
                if xstart == _old_off + 1 and xend == _old_off + _sel_clen:
                    focus_contig = _sel_contig
                    # mag_window stays correct: contig length is unchanged by reordering

        print(
            f"[start_bokeh_server] MAG view: mag={active_mag}, is_all={is_all}, sample={sample}, "
            f"sort={_sort_cat}/{_sort_metric}/{'asc' if _sort_ascending else 'desc'}, "
            f"features={mag_requested_features}",
            flush=True,
        )
        if enable_timing:
            t_params = time.perf_counter()
            _step = t_params - t_apply_start
            print(f"[timing] Parameter parsing: {_step:.3f}s{timing.tag(_step)}", flush=True)

        range_snapshot = self.plot_presenter.capture_range(active_mag, sample, is_all, xstart, xend)

        set_operation("apply/generate_plot_mag")
        if enable_timing:
            print(f"[timing] RSS before MAG composition: {rss_mb():.0f} MB", flush=True)
            t_plot = time.perf_counter()
        request = MagCompositionRequest(
            mag_name=active_mag,
            sample_name=sample or None,
            xstart=xstart,
            xend=xend,
            genbank_path=genbank_path if visibility.gene_map else None,
            tracks=TrackSelection(
                features=tuple(mag_requested_features),
                feature_types=tuple(selected_feature_types or ()),
                plot_isoforms=plot_isoforms,
                plot_sequence=visibility.sequence,
                plot_translation=visibility.translation,
                label_key=feature_label_key,
                color_rules=tuple(custom_colors or ()),
            ),
            display=CompositionDisplay(
                subplot_height=subplot_size,
                genemap_height=genemap_size,
                sequence_height=sequence_size,
                translation_height=translated_sequence_size,
                same_y_scale=same_y_scale,
                max_base_resolution=max_binning,
                max_genemap_window=max_genemap_window,
                max_sequence_window=mag_max_sequence_window,
                min_relative_value=min_coverage_freq,
            ),
            mag_track_colors=tuple(mag_track_colors or ()),
            max_track_dots=max_track_dots,
            is_all=is_all,
            allowed_samples=frozenset(mag_allowed_samples) if mag_allowed_samples is not None else None,
            max_samples=int(max_samples_input.value),
            enable_timing=enable_timing,
            ordering=MagOrdering(_sort_source, _sort_metric, _sort_ascending, _sort_sample_name),
            encoding_by_feature=_encoding_by_feature,
            data_cache=bindings.data_cache,
            profiler=profiler,
            sample_ordering=MagOrdering(
                _sample_sort_source,
                _sample_sort_metric,
                _sample_sort_ascending,
            ),
            focus_contig=focus_contig,
        )
        from ..sample_mag_pipeline import build_mag_plot

        grid, mag_meta = build_mag_plot(conn, request)
        if enable_timing:
            _step = time.perf_counter() - t_plot
            print(
                f"[timing] MAG composition (DB queries + plotting): {_step:.3f}s{timing.tag(_step)}",
                flush=True,
            )
            print(f"[timing] Bokeh model count in grid: {len(grid.references())}{timing.tag(0)}", flush=True)
            _overview = mag_meta.get("overview", {})
            print(
                f"[timing] MAG overview: boundaries={_overview.get('boundaries', 0)}, "
                f"visible_segments={_overview.get('visible_segments', 0)}, "
                f"data_sources={_overview.get('data_sources', 0)}",
                flush=True,
            )

        self.plot_presenter.install_range(
            grid,
            range_snapshot,
            subject=active_mag,
            sample=sample,
            is_all=is_all,
            start=xstart,
            end=xend,
        )

        actual_offsets = mag_meta.get("contig_offsets")
        if actual_offsets is not None:
            widgets["mag_to_contig_offsets"][active_mag] = actual_offsets
            _sync_from_to_for_selected_contig()

        _dots_shown = mag_meta.get("dots_shown", 0)
        _dots_total = mag_meta.get("dots_total", 0)
        if _dots_total > _dots_shown > 0:
            _dots_html = f"<i style='color:#555;font-size:0.85em;margin-left:10px'>{_dots_shown:,}/{_dots_total:,} dots only could be plotted on the MAG track</i>"
        else:
            _dots_html = ""
        toolbar_row = self.plot_presenter.mag_toolbar(_dots_html, has_samples)
        n_contigs = len(widgets["mag_to_contigs"].get(active_mag, []))
        return PlotResult(
            grid=grid,
            toolbar_row=toolbar_row,
            operation="apply/panel_update_mag",
            description=(
                f"MAG view for '{active_mag}', sample='{sample}', is_all={is_all}, "
                f"features={len(mag_requested_features)}, contigs={n_contigs}"
            ),
            total_label="MAG view",
            profiler=profiler,
        )

    def render_contig(self, inputs: ApplyInputs, t_apply_start: float) -> PlotResult:
        bindings = self.bindings
        _encoding_by_feature = bindings._encoding_by_feature
        _get_filtered_samples_for_contig = bindings._get_filtered_samples_for_contig
        _sample_sort_category_sources = bindings._sample_sort_category_sources
        _send_timing_ping = bindings._send_timing_ping
        combined_features_cbg = bindings.combined_features_cbg
        conn = bindings.conn
        enable_timing = bindings.enable_timing
        from_position_input = bindings.from_position_input
        max_samples_input = bindings.max_samples_input
        max_sequence_window_input = bindings.max_sequence_window_input
        sample_order_category_select = bindings.sample_order_category_select
        sample_order_direction = bindings.sample_order_direction
        sample_order_metric_select = bindings.sample_order_metric_select
        sequence_cbg = bindings.sequence_cbg
        to_position_input = bindings.to_position_input
        translated_sequence_cbg = bindings.translated_sequence_cbg
        widgets = bindings.widgets
        timing = bindings.timing
        set_operation = bindings.set_operation
        contig = inputs.contig
        has_samples = inputs.has_samples
        is_all = inputs.is_all
        sample = inputs.sample
        selected_feature_types = inputs.selected_feature_types
        genbank_path = inputs.genbank_path
        plot_isoforms = inputs.plot_isoforms
        feature_label_key = inputs.feature_label_key
        custom_colors = inputs.custom_colors
        active_variables_widgets = inputs.active_variables_widgets
        max_genemap_window = inputs.max_genemap_window
        same_y_scale = inputs.same_y_scale
        subplot_size = inputs.subplot_size
        genemap_size = inputs.genemap_size
        sequence_size = inputs.sequence_size
        translated_sequence_size = inputs.translated_sequence_size
        max_binning = inputs.max_binning
        min_coverage_freq = inputs.min_coverage_freq
        from ...shared.timing import ApplyProfiler

        profiler = ApplyProfiler("contig_all" if is_all else "contig_one", contig, sample, enabled=enable_timing)

        # Get contig length for validation
        contig_length = widgets["contig_lengths"].get(contig, 0)

        # Parse position inputs
        xstart = int(from_position_input.value) if from_position_input.value.strip() else 1
        xend = int(to_position_input.value) if to_position_input.value.strip() else contig_length

        max_seq_window = int(max_sequence_window_input.value)
        visibility = track_visibility(
            xend - xstart,
            max_genemap_window=max_genemap_window,
            max_sequence_window=max_seq_window,
            sequence_requested=sequence_cbg is not None and 0 in sequence_cbg.active,
            translation_requested=translated_sequence_cbg is not None and 0 in translated_sequence_cbg.active,
        )
        if not visibility.gene_map:
            print(
                f"Warning: Genome map will not be plotted for regions larger than {max_genemap_window} bp.",
                flush=True,
            )

        if sequence_cbg is not None and 0 in sequence_cbg.active:
            if not visibility.sequence:
                print(
                    f"Warning: Sequence will not be plotted for regions larger than {max_seq_window} bp.",
                    flush=True,
                )

        if translated_sequence_cbg is not None and 0 in translated_sequence_cbg.active:
            if not visibility.translation:
                print(
                    f"Warning: Translated sequence will not be plotted for regions larger than {max_seq_window} bp.",
                    flush=True,
                )

        # Check whether to preserve x-range from previous plot
        range_snapshot = self.plot_presenter.capture_range(contig, sample, is_all, xstart, xend)

        if enable_timing:
            t_params = time.perf_counter()
            _step = t_params - t_apply_start
            print(f"[timing] Parameter parsing: {_step:.3f}s{timing.tag(_step)}", flush=True)

        grid = None
        if is_all:
            # All-samples view: require exactly one variable selected from non-Genome modules
            # Genome module is shared (in Contigs section):
            #   - Gene map (index 0) is handled via genbank_path
            #   - Other Genome features (Repeats, etc.) are passed via genome_features parameter
            genome_features, selected_var = all_sample_features(combined_features_cbg, active_variables_widgets)

            filtered_samples = _get_filtered_samples_for_contig(contig)

            _sample_sort_cat = sample_order_category_for_view(sample_order_category_select.value, mag_view=False)
            _sample_sort_metric = sample_order_metric_select.value
            _sample_sort_ascending = sample_order_direction.active == 0
            _sample_sort_source = _sample_sort_category_sources.get(_sample_sort_cat)

            print(
                f"[start_bokeh_server] Generating plot for all samples with variable={selected_var}, contig={contig}, genome_features={genome_features}, filtered_samples={len(filtered_samples)}",
                flush=True,
            )
            genome_features = without_gene_map(genome_features, visibility.gene_map)
            set_operation("apply/generate_plot_all_samples")
            if enable_timing:
                print(f"[timing] RSS before compose_all_samples_plot: {rss_mb():.0f} MB", flush=True)
                t_plot = time.perf_counter()
            from ..all_samples_pipeline import build_all_samples_plot

            grid = build_all_samples_plot(
                conn,
                selected_var,
                contig,
                xstart=xstart,
                xend=xend,
                genbank_path=genbank_path,
                genome_features=genome_features if genome_features else None,
                allowed_samples=set(filtered_samples),
                feature_types=selected_feature_types,
                plot_isoforms=plot_isoforms,
                plot_sequence=visibility.sequence,
                plot_translated_sequence=visibility.translation,
                same_y_scale=same_y_scale,
                subplot_size=subplot_size,
                genemap_size=genemap_size,
                sequence_size=sequence_size,
                translated_sequence_size=translated_sequence_size,
                order_by_column=_sample_sort_metric,
                order_by_source=_sample_sort_source,
                order_ascending=_sample_sort_ascending,
                max_base_resolution=max_binning,
                max_genemap_window=max_genemap_window,
                min_relative_value=min_coverage_freq,
                feature_label_key=feature_label_key,
                custom_colors=custom_colors if custom_colors else None,
                max_samples=int(max_samples_input.value),
                enable_timing=enable_timing,
                encoding_by_feature=_encoding_by_feature,
                data_cache=bindings.data_cache,
                profiler=profiler,
            )
            if enable_timing:
                _step = time.perf_counter() - t_plot
                print(
                    f"[timing] compose_all_samples_plot (DB queries + plotting): {_step:.3f}s{timing.tag(_step)}",
                    flush=True,
                )
                print(f"[timing] Bokeh model count in grid: {len(grid.references())}{timing.tag(0)}", flush=True)
        else:
            # One-sample view: collect possibly-many requested features and call per-sample plot
            requested_features = one_sample_features(combined_features_cbg, active_variables_widgets)

            print(
                f"[start_bokeh_server] Generating plot for sample={sample}, contig={contig}, features={requested_features}",
                flush=True,
            )
            # Remove "Gene map" from requested_features if window too large
            requested_features = without_gene_map(requested_features, visibility.gene_map)
            # MAG view is handled by the early path above; MAG track is never shown in Contig view.
            if enable_timing:
                t_plot = time.perf_counter()
            request = SingleSampleCompositionRequest(
                contig_name=contig,
                sample_name=sample or None,
                xstart=xstart,
                xend=xend,
                genbank_path=genbank_path,
                tracks=TrackSelection(
                    features=tuple(requested_features),
                    feature_types=tuple(selected_feature_types or ()),
                    plot_isoforms=plot_isoforms,
                    plot_sequence=visibility.sequence,
                    plot_translation=visibility.translation,
                    label_key=feature_label_key,
                    color_rules=tuple(custom_colors or ()),
                ),
                display=CompositionDisplay(
                    subplot_height=subplot_size,
                    genemap_height=genemap_size,
                    sequence_height=sequence_size,
                    translation_height=translated_sequence_size,
                    max_base_resolution=max_binning,
                    max_genemap_window=max_genemap_window,
                    max_sequence_window=int(max_sequence_window_input.value),
                    min_relative_value=min_coverage_freq,
                ),
                enable_timing=enable_timing,
                encoding_by_feature=_encoding_by_feature,
                data_cache=bindings.data_cache,
                profiler=profiler,
            )
            from ..sample_mag_pipeline import build_single_sample_plot

            grid = build_single_sample_plot(conn, request)
            if enable_timing:
                _step = time.perf_counter() - t_plot
                print(
                    f"[timing] single-sample composition (DB queries + plotting): {_step:.3f}s{timing.tag(_step)}",
                    flush=True,
                )
                print(f"[timing] Bokeh model count in grid: {len(grid.references())}{timing.tag(0)}", flush=True)

        # Restore preserved x-range and update state
        self.plot_presenter.install_range(
            grid,
            range_snapshot,
            subject=contig,
            sample=sample,
            is_all=is_all,
            start=xstart,
            end=xend,
        )

        # Create toolbar-style row with buttons positioned top-right
        # Use Panel Row to mix Bokeh and Panel widgets
        toolbar_row = self.plot_presenter.contig_toolbar(has_samples, widgets["has_mags"])

        # Display the plot
        view_label = "all samples" if is_all else f"one sample ({sample})"
        return PlotResult(
            grid=grid,
            toolbar_row=toolbar_row,
            operation="apply/panel_update",
            description=(
                f"contig='{contig}', view={view_label}, "
                f"features={len(requested_features) if not is_all else '1 variable + genome'}"
            ),
            total_label="all samples" if is_all else "one sample",
            profiler=profiler,
        )
