"""Bokeh rendering for already-loaded ALL SAMPLES plot data."""

from __future__ import annotations

from bokeh.layouts import gridplot
from bokeh.models import Range1d

from ..models.plots import AllSamplesPlotData
from .features import render_feature_tracks


def _legacy_feature_dict(series):
    data = {key: list(value) for key, value in series.data.items()}
    data.update(
        type=series.style.plot_type,
        color=series.style.color,
        alpha=series.style.alpha,
        fill_alpha=series.style.fill_alpha,
        size=series.style.size,
        title=series.style.title,
        is_relative_scaled=series.is_relative_scaled,
        has_stats="mean" in data,
        has_sequences="sequence" in data,
    )
    return data


class AllSamplesRenderer:
    """Render plain data without a database connection or SQL."""

    def render_figures(self, plot_data: AllSamplesPlotData, subplot_height: int, same_y_scale: bool):
        shared_xrange = Range1d(plot_data.region.start, plot_data.region.end)
        genome_figures = []
        for track in plot_data.genome_tracks:
            figure = render_feature_tracks(
                [_legacy_feature_dict(series) for series in track.series],
                subplot_height,
                shared_xrange,
                sample_title=track.feature_name,
            )
            if figure is not None:
                genome_figures.append(figure)

        sample_figures = []
        for track in plot_data.sample_tracks:
            figure = render_feature_tracks(
                [_legacy_feature_dict(series) for series in track.series],
                subplot_height,
                shared_xrange,
                sample_title=track.sample_name,
            )
            if figure is not None:
                if same_y_scale and plot_data.y_max > 0:
                    figure.y_range = Range1d(plot_data.y_min, plot_data.y_max)
                sample_figures.append(figure)
        return shared_xrange, genome_figures, sample_figures

    @staticmethod
    def compose(figures):
        figures = list(figures)
        if not figures:
            raise ValueError("No plots to display")
        return gridplot([[figure] for figure in figures], merge_tools=True, sizing_mode="stretch_width")
