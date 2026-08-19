"""Bokeh rendering for already-decoded feature series."""

from __future__ import annotations

from bokeh.models import (
    BoxZoomTool,
    ColumnDataSource,
    HoverTool,
    NumeralTickFormatter,
    Range1d,
    WheelZoomTool,
)
from bokeh.plotting import figure


def _source_data(item):
    data = {"x": item["x"], "y": item["y"]}
    optional = (
        "first_pos",
        "last_pos",
        "width",
        "linked_start",
        "linked_end",
        "length",
        "repeat_positions",
        "partner_contig",
        "mean",
        "median",
        "std",
        "sequence",
        "sequence_prevalence",
        "codon_category",
        "codon_change",
        "aa_change",
    )
    for key in optional:
        if key in item:
            data[key] = item[key]
    if "sequence_prevalence" in item:
        data["sequence_prevalence_str"] = [
            f"{value / 100.0:.2f}" if value is not None else "" for value in item["sequence_prevalence"]
        ]
    return data


def _tooltips(items):
    if any(item.get("is_duplication", False) for item in items):
        return [
            ("First position", "@first_pos{0,0}"),
            ("Last position", "@last_pos{0,0}"),
            ("Linked start", "@linked_start{0,0}"),
            ("Linked end", "@linked_end{0,0}"),
            ("Length", "@length{0,0}"),
            ("Identity", "@y{0.01}%"),
        ]
    if any("partner_contig" in item for item in items):
        return [("Position", "@x{0,0}"), ("Identity (%)", "@y{0.01}"), ("Closest contig:", "@partner_contig")]
    if any("repeat_positions" in item for item in items):
        return [("Position", "@x{0,0}"), ("Value", "@y{0.00}"), ("Repeat positions", "@repeat_positions")]
    if any("width" in item and any(width != 1 for width in item["width"]) for item in items):
        return [("First position", "@first_pos{0,0}"), ("Last position", "@last_pos{0,0}"), ("Value", "@y{0.00}")]
    if any(item.get("has_stats", False) for item in items):
        has_mean = any(any(value is not None for value in item.get("mean", ())) for item in items)
        tips = [("Position", "@x{0,0}"), ("Value", "@y{0.00}")]
        if has_mean:
            tips.append(("Mean", "@mean{0.00}"))
        tips.append(("Median" if has_mean else "Median clipping", "@median{0.00}"))
        if has_mean:
            tips.append(("Std", "@std{0.00}"))
        if any(item.get("has_sequences", False) for item in items):
            tips.extend((("Main variant", "@sequence"), ("Prevalence", "@sequence_prevalence_str")))
        return tips
    if any(item.get("has_sequences", False) for item in items):
        tips = [
            ("Position", "@x{0,0}"),
            ("Value", "@y{0.00}"),
            ("Main variant", "@sequence"),
            ("Prevalence", "@sequence_prevalence_str"),
        ]
        if any("codon_category" in item for item in items):
            tips.extend((("Category", "@codon_category"), ("Codon", "@codon_change"), ("Amino acid", "@aa_change")))
        return tips
    return [("Position", "@x{0,0}"), ("Value", "@y{0.00}")]


class FeatureRenderer:
    """Render legacy-shaped, decoded feature dictionaries without database access."""

    def render(self, series, height, x_range, sample_title=None, show_tooltips=True):
        series = list(series or ())
        if not series or not any(any(value != 0 for value in item["y"]) for item in series):
            return None
        plot = figure(height=height, x_range=x_range, tools="xpan,reset,save")
        plot.xaxis.formatter = NumeralTickFormatter(format="0,0")
        for item in series:
            title = f"{sample_title} {item['title']}" if sample_title else item["title"]
            source_data = _source_data(item)
            source = ColumnDataSource(source_data)
            if item["type"] == "curve":
                plot.varea(
                    x="x",
                    y1=0,
                    y2="y",
                    source=source,
                    fill_color=item["color"],
                    fill_alpha=item["fill_alpha"],
                    legend_label=title,
                )
                plot.line(
                    x="x",
                    y="y",
                    source=source,
                    line_color=item["color"],
                    line_alpha=item["alpha"],
                    line_width=item["size"],
                    legend_label=title,
                )
            elif item["type"] == "bars":
                width = (
                    "width"
                    if "width" in source_data and any(value != 1 for value in source_data["width"])
                    else item["size"]
                )
                plot.vbar(
                    x="x",
                    bottom=0,
                    top="y",
                    source=source,
                    color=item["color"],
                    alpha=item["alpha"],
                    width=width,
                    legend_label=title,
                )
        if show_tooltips:
            plot.add_tools(HoverTool(tooltips=_tooltips(series), mode="vline"))
        plot.toolbar.logo = None
        plot.xgrid.visible = False
        plot.y_range.start = min(0, *(value for item in series for value in item["y"] if value is not None))
        if all(item.get("is_relative_scaled", False) for item in series):
            plot.y_range = Range1d(plot.y_range.start, 1)
        plot.ygrid.grid_line_alpha = 0.2
        plot.yaxis.axis_label = None
        plot.outline_line_color = None
        plot.min_border_left = 40
        plot.min_border_right = 10
        plot.legend.location = "top_left"
        if len(series) > 1:
            plot.legend.click_policy = "hide"
        wheel = WheelZoomTool(dimensions="width")
        plot.add_tools(wheel, BoxZoomTool(dimensions="width"))
        plot.toolbar.active_scroll = wheel
        return plot


def render_feature_tracks(series, height, x_range, sample_title=None, show_tooltips=True):
    """Functional entry point used by plot composers."""
    return FeatureRenderer().render(series, height, x_range, sample_title, show_tooltips)
