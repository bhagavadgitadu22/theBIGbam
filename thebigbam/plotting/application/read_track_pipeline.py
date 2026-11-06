"""Compose BAM tracks for ordered samples in contig or concatenated MAG coordinates."""

from ..models.read_plot import ReadPlotRequest, ReadTrackConfig
from ..renderers.read_alignments import render_read_track
from ..services.read_plots import prepare_read_plot
from ..shared.read_diagnostics import ReadProfiler
from ..shared.read_colors import READ_CLASS_COLORS

ALL_READ_CLASSES = frozenset(("primary", "secondary", "supplementary"))


def read_track_config(inputs, widgets, *, diagnose=False, reference_source=None):
    if not inputs.read_plot.enabled:
        return None
    colors = {**READ_CLASS_COLORS, **widgets.get("read_strand_colors", {})}
    return ReadTrackConfig(inputs.read_plot, widgets["bam_catalog"], ALL_READ_CLASSES,
                           colors, diagnose, reference_source)


def read_plot_request(config, sample, members, x_range):
    """Build the same immutable BAM request for single and multi-sample views."""
    if config is None or not members:
        return None
    bam = config.catalog.by_sample.get(sample)
    if bam is None or not bam.available:
        return None
    options = config.options
    total_length = max(offset + length for _, length, offset in members)
    start, end = (x_range.start, x_range.end) if hasattr(x_range, "start") else x_range
    return ReadPlotRequest(
        str(bam.path), members[0][0], total_length, int(start), int(end),
        options.max_reads, options.link_pairs, config.classes, config.colors,
        plot_height=options.plot_height, read_height=options.read_height,
        members=tuple(members), sample_name=sample,
        profiler=ReadProfiler() if config.diagnose else None,
        event_variables=options.event_variables,
        filters=options.filters,
        max_discrepancy_window=options.max_discrepancy_window,
        color_rules=options.color_rules,
        reference_source=config.reference_source,
        circular=bool(
            bam.circular
            or (config.reference_source is not None and config.reference_source.circular_sample(sample))
        ),
    )


def build_read_tracks(config, samples, members, x_range, cache):
    if config is None:
        return []
    figures = []
    for sample in samples:
        request = read_plot_request(config, sample, members, x_range)
        if request is None:
            continue
        data = prepare_read_plot(request, cache)
        plot, scrollbar = render_read_track(data, request, x_range)
        figures.append((plot, scrollbar))
    return figures
