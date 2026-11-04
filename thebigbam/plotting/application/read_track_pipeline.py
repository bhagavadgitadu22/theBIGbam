"""Compose BAM tracks for ordered samples in contig or concatenated MAG coordinates."""

from html import escape

from bokeh.models import Div

from ..models.read_plot import ReadPlotRequest, ReadTrackConfig
from ..renderers.read_alignments import render_read_track
from ..services.read_plots import prepare_read_plot
from ..services.read_selection import read_classes_for_subplots
from ..shared.read_diagnostics import ReadProfiler


def read_track_config(inputs, widgets, *, diagnose=False):
    if not inputs.read_plot.enabled:
        return None
    classes = read_classes_for_subplots(inputs.read_plot.selected_subplots, widgets["subplot_to_varnames"])
    return ReadTrackConfig(inputs.read_plot, widgets["bam_catalog"], classes, widgets["read_colors"], diagnose)


def build_read_tracks(config, samples, members, x_range, cache):
    if config is None:
        return []
    total_length = max(offset + length for _, length, offset in members)
    figures = []
    for sample in samples:
        bam = config.catalog.by_sample.get(sample)
        if bam is None or not bam.available:
            warning = config.catalog.warning_for_sample(sample) or "No usable indexed BAM."
            figures.append(Div(text=escape(f"{sample}: {warning}"), sizing_mode="stretch_width"))
            continue
        options = config.options
        request = ReadPlotRequest(
            str(bam.path), members[0][0], total_length, int(x_range.start), int(x_range.end),
            options.max_reads, options.link_pairs, config.classes, config.colors,
            plot_height=options.plot_height, read_height=options.read_height,
            members=tuple(members), sample_name=sample,
            profiler=ReadProfiler() if config.diagnose else None,
        )
        data = prepare_read_plot(request, cache)
        plot, scrollbar = render_read_track(data, request, x_range)
        figures.append((plot, scrollbar))
    return figures
