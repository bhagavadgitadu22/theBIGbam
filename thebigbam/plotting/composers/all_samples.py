"""Application-level composition for the ALL SAMPLES contig plot."""

from __future__ import annotations

import time

from ..models.plots import AllSamplesPlotRequest, DisplayOptions, GenomicRegion, SampleOrdering
from ..renderers.all_samples import AllSamplesRenderer
from ..repositories.all_samples import AllSamplesRepository
from ..services.all_samples import AllSamplesDataService
from ..shared.defaults import DEFAULT_GENEMAP_WINDOW
from ..shared.timing import PipelineTimings, profile_phase
from .genome_tracks import GenomeTrackComposer


def _expand_genome_features(features):
    expanded = []
    for feature in features or ():
        if feature.lower().strip() in {"repeats", "repeat", "direct repeats", "inverted repeats"}:
            expanded.extend(("Repeat count", "Max repeat identity"))
        else:
            expanded.append(feature)
    return tuple(dict.fromkeys(expanded))


def compose_all_samples_plot(
    conn,
    variable,
    contig_name,
    xstart=None,
    xend=None,
    subplot_size=130,
    genbank_path=None,
    genome_features=None,
    allowed_samples=None,
    feature_types=None,
    plot_isoforms=True,
    plot_sequence=False,
    plot_translated_sequence=False,
    same_y_scale=False,
    genemap_size=None,
    sequence_size=None,
    translated_sequence_size=None,
    order_by_column=None,
    order_by_source=None,
    order_ascending=True,
    max_base_resolution=None,
    max_genemap_window=None,
    min_relative_value=0.0,
    feature_label_key=None,
    custom_colors=None,
    max_samples=None,
    enable_timing=False,
    encoding_by_feature=None,
    data_cache=None,
    profiler=None,
):
    """Generate an all-samples plot through repository, service and renderer layers."""
    started = time.perf_counter()
    phase_timings = PipelineTimings()
    repository = AllSamplesRepository(conn)
    with phase_timings.phase("repository"):
        with profile_phase(profiler, "repository"):
            contig_id, locus_name, locus_size = repository.resolve_contig(contig_name)
    start = int(xstart if xstart is not None else 1)
    end = int(xend if xend is not None else locus_size)
    request = AllSamplesPlotRequest(
        contig_name=contig_name,
        variable=variable,
        region=GenomicRegion(start, end, max_base_resolution or 10_000),
        allowed_samples=frozenset(allowed_samples) if allowed_samples is not None else None,
        ordering=SampleOrdering(order_by_source, order_by_column, order_ascending, max_samples),
        display=DisplayOptions(subplot_size, same_y_scale, min_relative_value),
        genome_features=_expand_genome_features(genome_features),
        data_cache=data_cache,
        profiler=profiler,
    )
    service = AllSamplesDataService(repository)
    with phase_timings.phase("service"):
        with profile_phase(profiler, "service_transformation"):
            plot_data = service.load(request)
    if profiler is not None:
        for name, seconds in service.phase_seconds.items():
            profiler.seconds[name] = profiler.seconds.get(name, 0.0) + seconds
    renderer = AllSamplesRenderer()
    render_phase = profiler.phase("track_rendering") if profiler is not None else phase_timings.phase("renderer")
    with render_phase:
        shared_xrange, genome_figures, sample_figures = renderer.render_figures(plot_data, subplot_size, same_y_scale)

    genome_track_composer = GenomeTrackComposer(conn)
    members = ((contig_id, contig_name, locus_size, 0),)
    supplemental_figures = []
    threshold = max_genemap_window if max_genemap_window is not None else DEFAULT_GENEMAP_WINDOW
    if genbank_path and end - start <= threshold:
        annotation = genome_track_composer.gene_map(
            locus_name,
            locus_size,
            members,
            start,
            end,
            genemap_size if genemap_size is not None else subplot_size,
            shared_xrange,
            feature_types=feature_types,
            plot_isoforms=plot_isoforms,
            label_key=feature_label_key,
            color_rules=custom_colors,
        )
        if annotation is not None:
            supplemental_figures.append(annotation)
    sequence_height = sequence_size if sequence_size is not None else subplot_size // 2
    if plot_sequence:
        sequence = genome_track_composer.nucleotides(members, start, end, sequence_height, shared_xrange)
        if sequence is not None:
            supplemental_figures.append(sequence)
    if plot_translated_sequence:
        translated = genome_track_composer.translation(
            members,
            start,
            end,
            translated_sequence_size if translated_sequence_size is not None else sequence_height,
            shared_xrange,
        )
        if translated is not None:
            supplemental_figures.append(translated)

    grid_phase = profiler.phase("grid_assembly") if profiler is not None else phase_timings.phase("renderer")
    with grid_phase:
        grid = renderer.compose([*supplemental_figures, *genome_figures, *sample_figures])
    if enable_timing:
        phases = ", ".join(f"{name}={seconds:.3f}s" for name, seconds in service.phase_seconds.items())
        print(
            f"[timing] all-samples pipeline: total={time.perf_counter() - started:.3f}s "
            f"queries={repository.query_count} {phase_timings.format()} service_detail={phases} "
            f"genome_tracks={genome_track_composer.diagnostics}",
            flush=True,
        )
    return grid
