import time

from bokeh.models import Range1d

from ..composers.layout import assemble_grid
from ..models.composition import (
    MagCompositionRequest,
    SingleSampleCompositionRequest,
)
from ..renderers.composition import apply_per_feature_y_ranges, apply_primary_relative_y_range
from ..renderers.features import render_feature_tracks
from ..renderers.mag_tracks import MagTrackRenderer
from ..repositories.composition import CompositionRepository
from ..repositories.features import FeatureRepository
from ..repositories.mag_overview import MagOverviewRepository
from ..services.composition import CompositionDataService
from ..services.mag_overview import MagOverviewService
from ..shared.defaults import DEFAULT_GENEMAP_WINDOW, DEFAULT_SEQUENCE_WINDOW
from ..shared.timing import PipelineTimings, profile_phase
from .genome_track_pipeline import GenomeTrackPipeline


def build_single_sample_plot(conn, request: SingleSampleCompositionRequest):
    """Generate a single-sample plot from an immutable request."""
    contig_name = request.contig_name
    sample_name = request.sample_name
    xstart, xend = request.xstart, request.xend
    feature_types = request.tracks.feature_types
    plot_isoforms = request.tracks.plot_isoforms
    plot_sequence = request.tracks.plot_sequence
    plot_translated_sequence = request.tracks.plot_translation
    feature_label_key = request.tracks.label_key
    custom_colors = request.tracks.color_rules
    subplot_size = request.display.subplot_height
    genemap_size = request.display.genemap_height
    sequence_size = request.display.sequence_height
    translated_sequence_size = request.display.translation_height
    same_y_scale = request.display.same_y_scale
    max_genemap_window = request.display.max_genemap_window
    max_sequence_window = request.display.max_sequence_window
    genbank_path = request.genbank_path
    mag_name = request.mag_name
    enable_timing = request.enable_timing
    phase_timings = PipelineTimings()
    repository = CompositionRepository(conn)
    feature_repository = FeatureRepository(conn)
    profiler = request.profiler

    # Get contig characteristics
    with phase_timings.phase("repository"):
        with profile_phase(profiler, "repository"):
            contig = repository.contig(contig_name)
    contig_id, locus_name, locus_size = contig.contig_id, contig.name, contig.length
    print(f"Locus {locus_name} validated ({locus_size} bp)", flush=True)

    # --- Main gene annotation plot (only if genbank provided and window <= 100kb) ---
    shared_xrange = Range1d(0, locus_size)
    if xstart is not None and xend is not None:
        shared_xrange.start = xstart
        shared_xrange.end = xend
    genome_tracks = GenomeTrackPipeline(conn)
    genome_members = ((int(contig_id), locus_name, int(locus_size), 0),)

    _genemap_threshold = max_genemap_window if max_genemap_window is not None else DEFAULT_GENEMAP_WINDOW
    annotation_fig = None
    if genbank_path and xstart is not None and xend is not None and (xend - xstart) <= _genemap_threshold:
        annotation_fig = genome_tracks.gene_map(
            locus_name,
            locus_size,
            genome_members,
            xstart,
            xend,
            genemap_size if genemap_size is not None else subplot_size,
            shared_xrange,
            feature_types=feature_types,
            plot_isoforms=plot_isoforms,
            label_key=feature_label_key,
            color_rules=custom_colors,
        )

    # Get sample characteristics (optional ÃƒÂ¢Ã¢â€šÂ¬Ã¢â‚¬Å“ contig-level features work without a sample)
    with phase_timings.phase("repository"):
        with profile_phase(profiler, "repository"):
            sample_id, sample_name = repository.sample(sample_name)
    if sample_id is not None:
        print(f"Sample {sample_name} validated.", flush=True)

    # --- Add one subplot per feature requested ---
    # Requested features are variables like 'coverage', 'reads_starts', etc.
    subplots = []

    # --- Add sequence subplot right after annotation (top of data tracks) ---
    _seq_threshold = max_sequence_window if max_sequence_window is not None else DEFAULT_SEQUENCE_WINDOW
    if plot_sequence and xstart is not None and xend is not None and (xend - xstart) <= _seq_threshold:
        _seq_height = sequence_size if sequence_size is not None else subplot_size // 2
        seq_subplot = genome_tracks.nucleotides(genome_members, xstart, xend, _seq_height, shared_xrange)
        if seq_subplot:
            subplots.append(seq_subplot)
    elif plot_sequence and xstart is not None and xend is not None and (xend - xstart) > _seq_threshold:
        print(f"Sequence not plotted: window > {_seq_threshold} bp", flush=True)

    # --- Add translated sequence subplot right after DNA sequence ---
    if plot_translated_sequence and xstart is not None and xend is not None and (xend - xstart) <= _seq_threshold:
        _trans_height = (
            translated_sequence_size
            if translated_sequence_size is not None
            else (sequence_size if sequence_size is not None else subplot_size // 2)
        )
        trans_subplot = genome_tracks.translation(genome_members, xstart, xend, _trans_height, shared_xrange)
        if trans_subplot:
            subplots.append(trans_subplot)
    elif plot_translated_sequence and xstart is not None and xend is not None and (xend - xstart) > _seq_threshold:
        print(f"Translated sequence not plotted: window > {_seq_threshold} bp", flush=True)

    t_features = time.perf_counter()
    data_service = CompositionDataService(feature_repository, repository, cache=request.data_cache, profiler=profiler)
    with phase_timings.phase("service"):
        with profile_phase(profiler, "service_transformation"):
            prepared = data_service.for_contig(request, contig_id, sample_id)
    feature_subplots = []
    render_phase = profiler.phase("track_rendering") if profiler is not None else phase_timings.phase("renderer")
    with render_phase:
        for feature in (*prepared.contig, *prepared.sample):
            subplot = render_feature_tracks(feature.series, subplot_size, shared_xrange, show_tooltips=True)
            if subplot is None:
                continue
            subplots.append(subplot)
            if same_y_scale and feature in prepared.sample:
                maximum = max((value for item in feature.series for value in item["y"]), default=0)
                feature_subplots.append((feature.name, subplot, maximum))

    # Post-process y-ranges for same_y_scale (per-sample view)
    if same_y_scale and sample_id is not None and feature_subplots:
        # Find primary_reads max y
        primary_max = 0
        for fname, fig, my in feature_subplots:
            if fname == "Primary alignments":
                primary_max = max(primary_max, my)

        # If Primary alignments not plotted, fetch its data to find the max
        if primary_max == 0:
            try:
                primary_max = data_service.primary_maximum(request, contig_id, sample_id)
            except Exception:
                pass

        apply_primary_relative_y_range(feature_subplots, primary_max)

    if enable_timing:
        print(
            f"[timing]   feature subplots ({len(subplots)} plots): {time.perf_counter() - t_features:.3f}s", flush=True
        )

    # --- Optional MAG track (only when a MAG is selected in MAG-mode DBs) ---
    mag_fig = None
    if mag_name:
        try:
            mag_fig = MagTrackRenderer().render(repository.mag_contigs(mag_name), 30)
        except Exception as e:
            print(f"Error building MAG track for '{mag_name}': {e}", flush=True)
            mag_fig = None

    top_plots = [mag_fig] if mag_fig is not None else []
    grid_phase = profiler.phase("grid_assembly") if profiler is not None else phase_timings.phase("renderer")
    with grid_phase:
        grid = assemble_grid(top_plots + [annotation_fig] + subplots, enable_timing=enable_timing)
    if enable_timing:
        print(f"[timing] single-sample phases: {phase_timings.format()}", flush=True)
    return grid


def build_mag_plot(conn, request: MagCompositionRequest):
    """Generate a concatenated MAG plot from an immutable request."""
    mag_name = request.mag_name
    sample_name = request.sample_name
    xstart, xend = request.xstart, request.xend
    feature_types = request.tracks.feature_types
    plot_isoforms = request.tracks.plot_isoforms
    plot_sequence = request.tracks.plot_sequence
    plot_translated_sequence = request.tracks.plot_translation
    feature_label_key = request.tracks.label_key
    custom_colors = request.tracks.color_rules
    subplot_size = request.display.subplot_height
    genemap_size = request.display.genemap_height
    sequence_size = request.display.sequence_height
    translated_sequence_size = request.display.translation_height
    same_y_scale = request.display.same_y_scale
    max_genemap_window = request.display.max_genemap_window
    max_sequence_window = request.display.max_sequence_window
    genbank_path = request.genbank_path
    mag_track_colors = request.mag_track_colors
    max_track_dots = request.max_track_dots
    is_all = request.is_all
    enable_timing = request.enable_timing

    phase_timings = PipelineTimings()
    repository = CompositionRepository(conn)
    feature_repository = FeatureRepository(conn)
    profiler = request.profiler
    with phase_timings.phase("repository"):
        with profile_phase(profiler, "mag_context"):
            context = repository.mag(request)
    members = context.members
    mag_members = context.feature_members
    total_len = context.total_length
    xstart, xend = context.xstart, context.xend
    genome_tracks = GenomeTrackPipeline(conn)
    with profile_phase(profiler, "genome_member_preparation"):
        genome_members = context.genome_members
    print(f"MAG {mag_name}: {len(members)} contigs, {total_len} bp total", flush=True)

    shared_xrange = Range1d(
        xstart if xstart is not None else 1,
        xend if xend is not None else total_len,
    )

    # --- MAG track at top (same height as nucleotide sequence track) ---
    _seq_height = sequence_size if sequence_size is not None else subplot_size // 2
    track_dots = None
    if mag_track_colors:
        try:
            with profile_phase(profiler, "mag_overview_retrieval"):
                dot_xs, dot_colors, dot_total = MagOverviewService(MagOverviewRepository(conn)).annotation_dots(
                    mag_members, mag_track_colors, max_track_dots
                )
            if dot_xs:
                track_dots = {"xs": dot_xs, "colors": dot_colors, "total": dot_total}
        except Exception as e:
            print(f"[mag_track] Dot computation error: {e}", flush=True)
    try:
        with profile_phase(profiler, "mag_overview_rendering"):
            mag_fig = MagTrackRenderer().render(members, _seq_height, shared_xrange, track_dots)
    except Exception as e:
        print(f"Error building MAG track for '{mag_name}': {e}", flush=True)
        mag_fig = None

    # --- Combined gene map (only when window is small enough) ---
    _genemap_threshold = max_genemap_window if max_genemap_window is not None else DEFAULT_GENEMAP_WINDOW
    annotation_fig = None
    if genbank_path and xstart is not None and xend is not None and (xend - xstart) <= _genemap_threshold:
        _gm_size = genemap_size if genemap_size is not None else subplot_size
        try:
            with profile_phase(profiler, "gene_map"):
                annotation_fig = genome_tracks.gene_map(
                    mag_name,
                    total_len,
                    genome_members,
                    xstart,
                    xend,
                    _gm_size,
                    shared_xrange,
                    feature_types=feature_types,
                    plot_isoforms=plot_isoforms,
                    label_key=feature_label_key,
                    color_rules=custom_colors,
                )
        except Exception as e:
            print(f"Error building combined gene map for MAG '{mag_name}': {e}", flush=True)

    # --- Resolve sample ---
    with phase_timings.phase("repository"):
        with profile_phase(profiler, "repository"):
            sample_id, sample_name = repository.sample(sample_name if not is_all else None)

    t_features = time.perf_counter()
    subplots = []

    # --- Sequence subplot (stitched across member contigs) ---
    _seq_threshold = max_sequence_window if max_sequence_window is not None else DEFAULT_SEQUENCE_WINDOW
    if plot_sequence and xstart is not None and xend is not None and (xend - xstart) <= _seq_threshold:
        with profile_phase(profiler, "nucleotide_track"):
            seq_subplot = genome_tracks.nucleotides(genome_members, xstart, xend, _seq_height, shared_xrange)
        if seq_subplot:
            subplots.append(seq_subplot)
    elif plot_sequence and xstart is not None and xend is not None and (xend - xstart) > _seq_threshold:
        print(f"Sequence not plotted: window > {_seq_threshold} bp", flush=True)

    # --- Translated-sequence subplot (CDS across member contigs) ---
    if plot_translated_sequence and xstart is not None and xend is not None and (xend - xstart) <= _seq_threshold:
        _trans_height = (
            translated_sequence_size
            if translated_sequence_size is not None
            else (sequence_size if sequence_size is not None else subplot_size // 2)
        )
        with profile_phase(profiler, "translation_track"):
            trans_subplot = genome_tracks.translation(genome_members, xstart, xend, _trans_height, shared_xrange)
        if trans_subplot:
            subplots.append(trans_subplot)
    elif plot_translated_sequence and xstart is not None and xend is not None and (xend - xstart) > _seq_threshold:
        print(f"Translated sequence not plotted: window > {_seq_threshold} bp", flush=True)

    with phase_timings.phase("service"):
        with profile_phase(profiler, "service_transformation"):
            prepared = CompositionDataService(
                feature_repository, repository, cache=request.data_cache, profiler=profiler
            ).for_mag(request, context, sample_id)
    feature_subplot_info = []
    render_phase = profiler.phase("track_rendering") if profiler is not None else phase_timings.phase("renderer")
    with render_phase:
        for feature in (*prepared.contig, *prepared.sample):
            subplot = render_feature_tracks(feature.series, subplot_size, shared_xrange, show_tooltips=True)
            if subplot is None:
                continue
            subplots.append(subplot)
            if same_y_scale and feature.sample_name is not None:
                maximum = max(
                    (value for item in feature.series for value in item.get("y", []) if value is not None),
                    default=0,
                )
                feature_subplot_info.append((feature.name, subplot, maximum))
    if same_y_scale and feature_subplot_info:
        apply_per_feature_y_ranges(feature_subplot_info)

    if enable_timing:
        print(
            f"[timing]   MAG feature subplots ({len(subplots)} plots): {time.perf_counter() - t_features:.3f}s",
            flush=True,
        )

    top_plots = [mag_fig] if mag_fig is not None else []
    all_plots = top_plots + ([annotation_fig] if annotation_fig is not None else []) + subplots
    grid_phase = profiler.phase("grid_assembly") if profiler is not None else phase_timings.phase("renderer")
    with grid_phase:
        grid = assemble_grid(all_plots, empty_message="No plots to display for MAG view", enable_timing=enable_timing)
    actual_contig_offsets = {name: off for name, _clen, off in members}
    mag_meta = {
        "dots_shown": len(track_dots["xs"]) if track_dots else 0,
        "dots_total": track_dots.get("total", 0) if track_dots else 0,
        "contig_offsets": actual_contig_offsets,
        "overview": mag_fig.tags[-1] if mag_fig is not None and mag_fig.tags else {},
        "phase_seconds": dict(phase_timings.seconds),
    }
    if enable_timing:
        print(f"[timing] MAG phases: {phase_timings.format()}", flush=True)
    return grid, mag_meta
