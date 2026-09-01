"""Bokeh rendering for prepared filtering distributions."""

from __future__ import annotations

from typing import Any, Callable, Mapping

import panel as pn
from bokeh.events import Tap
from bokeh.io import curdoc
from bokeh.models import ColumnDataSource, CustomJSTickFormatter, HoverTool, Range1d, Span, TapTool, TextInput
from bokeh.models.callbacks import CustomJS
from bokeh.plotting import figure as bk_figure

from ..shared.styles import panel_stylesheet


class FilterVisualizations:
    def __init__(
        self,
        metadata_service: Any,
        filtering_metadata: Mapping[str, Any],
        enable_timing: bool,
        muted_buttons_stylesheet: Any,
        refresh_on_filter_change: Callable[[], None],
        filter_encode: Mapping[str, float] | None = None,
        set_operation: Callable[[str], None] | None = None,
        record_action: Callable[[str, Mapping[str, Any]], Any] | None = None,
        occurrence_for: Callable[[Any], int] | None = None,
    ) -> None:
        self.metadata_service = metadata_service
        self.filtering_metadata = filtering_metadata
        self.enable_timing = enable_timing
        self.muted_buttons_stylesheet = muted_buttons_stylesheet
        self.refresh_on_filter_change = refresh_on_filter_change
        self.filter_encode = filter_encode or {}
        self.set_operation = set_operation or (lambda _operation: None)
        self.record_action = record_action
        self.occurrence_for = occurrence_for or (lambda _row: 1)

    def build_numeric_histogram(self, row_data, category, col_name, spinner, log_mode=False, log_y=False):
        """Build a histogram with draggable threshold bar for a numeric column."""
        import numpy as np

        col_info = self.filtering_metadata.get(category, {}).get("columns", {}).get(col_name, {})
        if col_info.get("type") != "numeric" or col_info.get("is_bool"):
            row_data["histogram_pane"] = None
            row_data["histogram_fig"] = None
            row_data["threshold_span"] = None
            return None

        scale = self.filter_encode.get(col_name)
        self.set_operation("filter_histogram")
        try:
            bin_result = self.metadata_service.histogram_bins(
                category, col_name, n_bins=50, log_mode=log_mode, scale=scale
            )
        finally:
            self.set_operation("idle")
        if bin_result is None:
            row_data["histogram_pane"] = None
            row_data["histogram_fig"] = None
            row_data["threshold_span"] = None
            return None

        edges, counts = bin_result.edges, bin_result.counts
        row_data["log_mode"] = log_mode
        row_data["log_y"] = log_y

        if not scale:
            import math as _m

            bin_width = float(edges[-1] - edges[0]) / max(len(edges) - 1, 1)
            if bin_width > 0 and bin_width < 1:
                spinner.step = 10 ** _m.floor(_m.log10(bin_width))
            elif bin_width >= 1:
                spinner.step = 1
            else:
                spinner.step = 0.001

        total_count = int(sum(counts))
        mids = ((edges[:-1] + edges[1:]) / 2).tolist()
        pcts = [(c / total_count * 100 if total_count else 0) for c in counts.tolist()]
        if log_y:
            counts = np.array([np.log10(c) if c > 0 else 0 for c in counts])
        max_count = float(max(counts)) if len(counts) else 1
        if log_mode:
            display_min = [f"{10**v:.4g}" for v in edges[:-1].tolist()]
            display_max = [f"{10**v:.4g}" for v in edges[1:].tolist()]
        else:
            display_min = [f"{v:.4g}" for v in edges[:-1].tolist()]
            display_max = [f"{v:.4g}" for v in edges[1:].tolist()]
        hist_source = ColumnDataSource(
            data=dict(
                left=edges[:-1].tolist(),
                right=edges[1:].tolist(),
                top=counts.tolist(),
                bottom=[0] * len(counts),
                mid=mids,
                min_val=display_min,
                max_val=display_max,
                pct=[f"{p:.1f}%" for p in pcts],
            )
        )

        fig = bk_figure(
            height=60,
            sizing_mode="stretch_width",
            toolbar_location=None,
            tools="",
            outline_line_color=None,
            min_border_left=5,
            min_border_right=5,
            min_border_top=2,
            min_border_bottom=15,
            y_range=Range1d(0, max_count * 1.15),
        )
        fig.quad(
            left="left",
            right="right",
            top="top",
            bottom="bottom",
            source=hist_source,
            fill_color="#00b17c",
            fill_alpha=0.6,
            line_color="white",
            line_width=0.5,
            nonselection_fill_color="#00b17c",
            nonselection_fill_alpha=0.6,
            nonselection_line_color="white",
            nonselection_line_width=0.5,
            selection_fill_color="#00b17c",
            selection_fill_alpha=0.6,
            selection_line_color="white",
            selection_line_width=0.5,
        )
        fig.xaxis.minor_tick_line_color = None
        fig.xaxis.major_label_text_font_size = "7pt"
        if log_mode:
            fig.xaxis.formatter = CustomJSTickFormatter(code="return (Math.pow(10, tick)).toPrecision(3);")
        fig.yaxis.visible = False
        fig.xgrid.grid_line_color = None
        fig.ygrid.grid_line_color = None

        fig.add_tools(HoverTool(tooltips=[("Min", "@min_val"), ("Max", "@max_val"), ("%", "@pct")], mode="vline"))

        if log_mode and spinner.value and spinner.value > 0:
            init_loc = np.log10(spinner.value)
        elif spinner.value:
            init_loc = spinner.value
        else:
            init_loc = float("-inf")
        threshold_span = Span(location=init_loc, dimension="height", line_color="black", line_width=2)
        fig.add_layout(threshold_span)

        is_log = [1] if log_mode else [0]
        fig.add_tools(TapTool())
        click_js = CustomJS(
            args=dict(span=threshold_span, spinner=spinner, source=hist_source, is_log=is_log),
            code="""
            const x = cb_obj.x;
            if (x === undefined || isNaN(x)) return;
            const step = spinner.step || 1;
            const decimals = step < 1 ? Math.ceil(-Math.log10(step)) : 0;
            const lefts = source.data.left;
            const rights = source.data.right;
            const mids = source.data.mid;
            let val = x;
            for (let i = 0; i < lefts.length; i++) {
                if (x >= lefts[i] && x <= rights[i]) {
                    val = mids[i];
                    break;
                }
            }
            span.location = val;
            let real_val = is_log[0] ? Math.pow(10, val) : val;
            const rounded = parseFloat((Math.round(real_val / step) * step).toFixed(decimals));
            spinner.value = rounded;
        """,
        )
        fig.js_on_event(Tap, click_js)

        hist_pane = pn.pane.Bokeh(fig, sizing_mode="stretch_width", margin=(0, 0, 0, 0))

        def _rebuild_histogram(new_log_x, new_log_y):

            hist_container = row_data["hist_container"]
            hist_container.loading = True
            row_data["loading_gen"] += 1
            gen = row_data["loading_gen"]

            def _do_rebuild():
                try:
                    if row_data["loading_gen"] != gen:
                        return
                    row_data["histogram_pane"] = None
                    row_data["histogram_fig"] = None
                    row_data["threshold_span"] = None
                    result = self.build_numeric_histogram(
                        row_data, category, col_name, spinner, log_mode=new_log_x, log_y=new_log_y
                    )
                    if result:
                        hist_container.objects = result
                finally:
                    # A newer rebuild owns the loading indicator when its
                    # generation has advanced; only its callback may clear it.
                    if row_data["loading_gen"] == gen:
                        hist_container.loading = False

            curdoc().add_next_tick_callback(_do_rebuild)

        log_x_btn = pn.widgets.Button(
            name="log x",
            height=30,
            margin=(0, 2, 3, 0),
            button_type="primary" if log_mode else "default",
            description="Only positive values will be considered for the plot",
            stylesheets=[panel_stylesheet(self.muted_buttons_stylesheet)],
            css_classes=["action-muted"],
        )
        log_y_btn = pn.widgets.Button(
            name="log y",
            height=30,
            margin=(0, 0, 3, 2),
            button_type="primary" if log_y else "default",
            stylesheets=[panel_stylesheet(self.muted_buttons_stylesheet)],
            css_classes=["action-muted"],
        )

        def _on_log_x(event):
            new_log = not row_data.get("log_mode", False)
            if self.record_action is not None:
                self.record_action(
                    "filter_distribution_scale",
                    {
                        "category": category,
                        "column": col_name,
                        "occurrence": self.occurrence_for(row_data),
                        "axis": "x",
                        "enabled": new_log,
                    },
                )
            _rebuild_histogram(new_log, row_data.get("log_y", False))

        def _on_log_y(event):
            new_log_y = not row_data.get("log_y", False)
            if self.record_action is not None:
                self.record_action(
                    "filter_distribution_scale",
                    {
                        "category": category,
                        "column": col_name,
                        "occurrence": self.occurrence_for(row_data),
                        "axis": "y",
                        "enabled": new_log_y,
                    },
                )
            _rebuild_histogram(row_data.get("log_mode", False), new_log_y)

        log_x_btn.on_click(_on_log_x)
        log_y_btn.on_click(_on_log_y)
        row_data["set_histogram_scale"] = lambda axis, enabled: _rebuild_histogram(
            enabled if axis == "x" else row_data.get("log_mode", False),
            enabled if axis == "y" else row_data.get("log_y", False),
        )

        pane = pn.Column(
            pn.Row(log_x_btn, log_y_btn, margin=0),
            hist_pane,
            sizing_mode="stretch_width",
            margin=0,
            styles={"background": "white", "border-radius": "5px", "padding": "5px"},
        )
        row_data["histogram_pane"] = pane
        row_data["histogram_fig"] = fig
        row_data["threshold_span"] = threshold_span
        null_stats = self.metadata_service.column_null_stats(category, col_name)
        if null_stats is not None:
            non_null_count, total_possible = null_stats
            approximation = (
                f" Approximate distribution — {bin_result.sampled_rows:,} sampled values."
                if bin_result.approximate else ""
            )
            label_html = (
                f'<span style="font-size:11px; color:#666; font-style:italic;">'
                f"Used {non_null_count:,} times (out of {total_possible:,} possible).{approximation}</span>"
            )
            dist_label = pn.pane.HTML(label_html, sizing_mode="stretch_width", margin=(2, 5, 0, 5))
            return [dist_label, pane]
        if bin_result.approximate:
            label_html = (
                '<span style="font-size:11px; color:#666; font-style:italic;">'
                f"Approximate distribution — {bin_result.sampled_rows:,} sampled values.</span>"
            )
            return [pn.pane.HTML(label_html, sizing_mode="stretch_width", margin=(2, 5, 0, 5)), pane]
        return [pane]

    def build_text_treemap(self, row_data, category, col_name, input_ref, hist_container):
        """Build a 1D treemap showing value distribution for a text column."""
        from bokeh.palettes import Category20_20

        value_counts = self.metadata_service.value_counts(category, col_name)
        if not value_counts:
            row_data["treemap_pane"] = None
            return None

        total = sum(c for _, c in value_counts)
        if total == 0:
            row_data["treemap_pane"] = None
            return None

        filtered = [(v, c) for v, c in value_counts if (c / total * 100) >= 1.0]
        if not filtered:
            row_data["treemap_pane"] = None
            return None

        lefts, rights, labels, values, pcts, colors = [], [], [], [], [], []
        cursor = 0.0
        for i, (val, cnt) in enumerate(filtered):
            pct = cnt / total * 100
            lefts.append(cursor)
            rights.append(cursor + pct)
            labels.append(str(val))
            values.append(str(val))
            pcts.append(f"{pct:.1f}%")
            colors.append(Category20_20[i % 20])
            cursor += pct

        other_count = total - sum(c for _, c in filtered)
        if other_count > 0:
            other_pct = other_count / total * 100
            lefts.append(cursor)
            rights.append(cursor + other_pct)
            labels.append("Other")
            values.append("")
            pcts.append(f"{other_pct:.1f}%")
            colors.append("#cccccc")

        source = ColumnDataSource(
            data=dict(
                left=lefts,
                right=rights,
                top=[1] * len(lefts),
                bottom=[0] * len(lefts),
                label=labels,
                value=values,
                pct=pcts,
                color=colors,
            )
        )

        fig = bk_figure(
            height=40,
            sizing_mode="stretch_width",
            toolbar_location=None,
            tools="",
            outline_line_color=None,
            x_range=Range1d(0, 100),
            y_range=Range1d(0, 1),
            min_border_left=0,
            min_border_right=0,
            min_border_top=0,
            min_border_bottom=0,
        )
        fig.quad(
            left="left",
            right="right",
            top="top",
            bottom="bottom",
            source=source,
            fill_color="color",
            line_color="white",
            line_width=1,
        )
        fig.xaxis.visible = False
        fig.yaxis.visible = False
        fig.xgrid.grid_line_color = None
        fig.ygrid.grid_line_color = None

        fig.add_tools(HoverTool(tooltips=[("Value", "@label"), ("%", "@pct")]), TapTool())

        bridge = TextInput(value="", visible=False)
        bridge_pane = pn.pane.Bokeh(bridge, height=0, sizing_mode="fixed", margin=0)

        tap_cb = CustomJS(
            args=dict(source=source, bridge=bridge),
            code="""
            const x = cb_obj.x;
            if (x === undefined || isNaN(x)) return;
            const data = source.data;
            for (let i = 0; i < data.left.length; i++) {
                if (x >= data.left[i] && x <= data.right[i]) {
                    if (data.value[i]) {
                        bridge.value = data.value[i];
                    }
                    break;
                }
            }
        """,
        )
        fig.js_on_event(Tap, tap_cb)

        def on_bridge_change(attr, old, new):
            if new:
                widget = input_ref["widget"]
                if hasattr(widget, "value"):
                    widget.value = new
                bridge.value = ""
                self.refresh_on_filter_change()

        bridge.on_change("value", on_bridge_change)

        pane = pn.pane.Bokeh(fig, sizing_mode="stretch_width", margin=(0, 0, 0, 0))
        row_data["treemap_pane"] = pane
        row_data["bridge_input"] = bridge
        null_stats = self.metadata_service.column_null_stats(category, col_name)
        if null_stats is not None:
            non_null_count, total_possible = null_stats
            label_html = (
                f'<span style="font-size:11px; color:#666; font-style:italic;">'
                f"Used {non_null_count:,} times (out of {total_possible:,} possible)</span>"
            )
            dist_label = pn.pane.HTML(label_html, sizing_mode="stretch_width", margin=(2, 5, 0, 5))
            return [dist_label, pane, bridge_pane]
        return [pane, bridge_pane]
