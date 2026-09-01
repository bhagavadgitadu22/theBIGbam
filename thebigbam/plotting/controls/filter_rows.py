"""Factory for dynamic filtering query rows."""

from __future__ import annotations

from typing import Any, Callable, Mapping

import panel as pn
from bokeh.io import curdoc
from bokeh.models.widgets import Select, Spinner, TextInput

from ..renderers.filter_distributions import FilterVisualizations
from ..shared.styles import panel_stylesheet
from .distinct_value_select import build_distinct_value_select


class FilterRowFactory:
    def __init__(
        self,
        metadata_service: Any,
        filtering_metadata: Mapping[str, Any],
        columns: Mapping[str, Any],
        raw_columns: Mapping[str, Any],
        visualizations: FilterVisualizations,
        refresh: Callable[[], None],
        stylesheet: Any,
        enable_timing: bool,
        filter_encode: Mapping[str, float] | None = None,
        record_action: Callable[[str, Mapping[str, Any]], Any] | None = None,
    ) -> None:
        self.metadata_service = metadata_service
        self.filtering_metadata = filtering_metadata
        self.columns = columns
        self.raw_columns = raw_columns
        self.visualizations = visualizations
        self.refresh = refresh
        self.stylesheet = stylesheet
        self.enable_timing = enable_timing
        self.filter_encode = filter_encode or {}
        self.record_action = record_action
        self.controller: Any | None = None
        self._row_sequence = 0

    def attach_controller(self, controller: Any) -> None:
        self.controller = controller
        self.visualizations.occurrence_for = self._occurrence_for

    def _occurrence_for(self, target: Any) -> int:
        matches = [
            row_data
            for section in self.controller.sections
            for row_data in section["rows"]
            if row_data["category_select"].value == target["category_select"].value
            and row_data["subcategory_select"].value == target["subcategory_select"].value
        ]
        return matches.index(target) + 1

    def _remove_row(self, row_data) -> None:
        row_data["loading_gen"] += 1
        if self.controller.count_rows() <= 1:
            return
        for section in self.controller.sections:
            if row_data not in section["rows"]:
                continue
            section["rows"].remove(row_data)
            self.controller.rebuild_section(section)
            if not section["rows"]:
                self.controller.sections.remove(section)
                self.controller.rebuild_content()
            self.refresh()
            return

    def _toggle_distribution(self, row_data, *, record_action: bool = True) -> None:
        container = row_data["hist_container"]
        row_data["loading_gen"] += 1
        category = row_data["category_select"].value
        column = row_data["subcategory_select"].value
        if record_action and self.record_action is not None:
            self.record_action(
                "filter_lookup",
                {"category": category, "column": column, "opening": not bool(container.objects)},
            )
        if container.objects:
            container.objects = []
            for key in ("histogram_pane", "histogram_fig", "threshold_span", "treemap_pane", "bridge_input"):
                row_data[key] = None
            return
        generation = row_data["loading_gen"]

        def load_distribution():
            if row_data["loading_gen"] != generation:
                return
            info = self.filtering_metadata.get(category, {}).get("columns", {}).get(column, {})
            if info.get("type") == "numeric" and not info.get("is_bool"):
                result = self.visualizations.build_numeric_histogram(
                    row_data, category, column, row_data["input_ref"]["widget"]
                )
            elif info.get("type") == "text":
                result = self.visualizations.build_text_treemap(
                    row_data, category, column, row_data["input_ref"], container
                )
            else:
                stats = self.metadata_service.column_null_stats(category, column)
                result = None
                if stats is not None:
                    non_null_count, total_possible = stats
                    result = [
                        pn.pane.HTML(
                            '<span style="font-size:11px; color:#666; font-style:italic;">'
                            f"Used {non_null_count:,} times (out of {total_possible:,} possible)</span>",
                            sizing_mode="stretch_width",
                            margin=(2, 5, 0, 5),
                        )
                    ]
            if result and row_data["loading_gen"] == generation:
                container.objects = result

        curdoc().add_next_tick_callback(load_distribution)

    def set_distribution(
        self,
        category: str,
        column: str,
        *,
        occurrence: int = 1,
        opening: bool = True,
    ) -> None:
        """Set a current row's distribution state using semantic identity."""
        if not isinstance(occurrence, int) or isinstance(occurrence, bool) or occurrence < 1:
            raise ValueError("filter lookup occurrence must be a positive integer")
        if self.controller is None:
            raise RuntimeError("filter row controller is not attached")
        matches = [
            row_data
            for section in self.controller.sections
            for row_data in section["rows"]
            if row_data["category_select"].value == category
            and row_data["subcategory_select"].value == column
        ]
        if occurrence > len(matches):
            raise ValueError(
                f"filter lookup target {category}.{column} occurrence {occurrence} does not exist; "
                f"found {len(matches)} matching current rows"
            )
        row_data = matches[occurrence - 1]
        is_open = bool(row_data["hist_container"].objects)
        if is_open != opening:
            self._toggle_distribution(row_data, record_action=False)

    def set_distribution_scale(
        self,
        category: str,
        column: str,
        *,
        occurrence: int,
        axis: str,
        enabled: bool,
    ) -> None:
        matches = [
            row_data
            for section in self.controller.sections
            for row_data in section["rows"]
            if row_data["category_select"].value == category
            and row_data["subcategory_select"].value == column
        ]
        if occurrence < 1 or occurrence > len(matches):
            raise ValueError("filter distribution scale target does not exist")
        setter = matches[occurrence - 1].get("set_histogram_scale")
        if setter is None:
            raise ValueError("filter distribution is not open")
        setter(axis, enabled)

    def create_row(self, section_data, initial_category=None, initial_column=None, initial_operator=None):
        """Create a single query row with cascading selects, comparison, dynamic input and remove button.

        initial_category/initial_column/initial_operator let a caller (e.g. the
        settings restore path) build a row for a specific saved column/operator
        synchronously, using the same construction path as the default first
        row instead of simulating UI changes through the (partly deferred)
        update_subcategories/update_input_widget callbacks.
        """
        row_index = self._row_sequence
        self._row_sequence += 1

        # Get categories from metadata
        categories = list(self.filtering_metadata.keys())
        if not categories:
            categories = ["No data"]

        if initial_category is None or initial_category not in categories:
            initial_category = categories[0]
        initial_columns_raw = list(self.filtering_metadata.get(initial_category, {}).get("columns", {}).keys())
        if not initial_columns_raw:
            initial_columns_raw = ["No columns"]
        initial_columns = [(c, c.replace("_", " ").replace("percentage", "(%)")) for c in initial_columns_raw]
        if initial_column is None or initial_column not in initial_columns_raw:
            initial_column = initial_columns_raw[0]

        # Determine initial column type
        initial_col_info = self.filtering_metadata.get(initial_category, {}).get("columns", {}).get(initial_column, {})
        initial_is_text = initial_col_info.get("type") == "text"

        # First level select (categories)
        category_select = Select(name=f"benchmark-filter-{row_index}-category", options=categories, value=initial_category, width=70, margin=0)

        # Second level select (columns)
        subcategory_select = Select(
            name=f"benchmark-filter-{row_index}-metric",
            options=initial_columns, value=initial_column, sizing_mode="stretch_width", margin=0
        )

        # Comparison operator select - "=" and "!=" for text, all operators for numeric
        _default_operator = "=" if initial_is_text else ">"
        _initial_ops = ["=", "!=", "has", "has not"] if initial_is_text else ["=", ">", "<", "!="]
        if initial_operator is not None and initial_operator in _initial_ops:
            _default_operator = initial_operator
        comparison_select = Select(name=f"benchmark-filter-{row_index}-operator", options=_initial_ops, value=_default_operator, width=50, margin=0)

        # Container for the dynamic input widget
        input_container = pn.Column(width=90, margin=0)

        def make_searchable_input(category, column, *, placeholder="Search..."):
            return build_distinct_value_select(
                self.metadata_service,
                category,
                column,
                name=f"benchmark-filter-{row_index}-value",
                on_value=self.refresh,
                width=90,
            )

        # Create initial input widget based on column type
        if initial_is_text and _default_operator in ("has", "has not"):
            initial_input = TextInput(value="", placeholder="Search...", width=90, margin=(0, 2, 0, 0))
            input_container.objects = [initial_input]
            initial_input.on_change("value", lambda attr, old, new: self.refresh())
            initial_is_panel = False
        elif initial_is_text:
            initial_input = make_searchable_input(initial_category, initial_column)
            input_container.objects = [initial_input]
            initial_is_panel = True
        else:
            _enc_scale = self.filter_encode.get(initial_column)
            _step = 1.0 / _enc_scale if _enc_scale else 1
            initial_input = Spinner(name=f"benchmark-filter-{row_index}-value", value=None, step=_step, placeholder="Value...", width=90, margin=(0, 2, 0, 0))
            input_container.objects = [initial_input]

            # Add callback for Bokeh Spinner (also syncs histogram threshold)
            def _on_initial_spinner(attr, old, new):
                if row_data.get("threshold_span") is not None and new is not None:
                    import math as _m

                    if row_data.get("log_mode") and new > 0:
                        row_data["threshold_span"].location = _m.log10(new)
                    else:
                        row_data["threshold_span"].location = new
                self.refresh()

            initial_input.on_change("value", _on_initial_spinner)
            initial_is_panel = False

        # Store reference to current input widget (for later retrieval)
        current_input_ref = {"widget": initial_input, "is_panel": initial_is_panel}
        refresh_state = {"pending": False}
        input_update_state = {"active": False}

        def schedule_refresh():
            """Let the browser render control changes before expensive filtering."""
            if refresh_state["pending"]:
                return
            refresh_state["pending"] = True

            def _run_refresh():
                refresh_state["pending"] = False
                self.refresh()

            curdoc().add_next_tick_callback(_run_refresh)

        def update_input_widget(col_name):
            """Update the input widget based on column type.

            Synchronous part: update comparison/input widgets from metadata
            (no DB call) and clear the inset immediately.  Deferred part:
            resolve distinct values and rebuild the inset via next-tick.
            """

            category = category_select.value
            col_info = self.filtering_metadata.get(category, {}).get("columns", {}).get(col_name, {})
            is_text = col_info.get("type") == "text"

            # --- synchronous: update comparison + input widget (no DB) ---
            input_update_state["active"] = True
            try:
                if is_text:
                    comparison_select.options = ["=", "!=", "has", "has not"]
                    if comparison_select.value not in comparison_select.options:
                        comparison_select.value = "="
                else:
                    comparison_select.options = ["=", ">", "<", "!="]
                    if comparison_select.value not in comparison_select.options:
                        comparison_select.value = "="
            finally:
                input_update_state["active"] = False

            if is_text:
                current_op = comparison_select.value
                if current_op in ("has", "has not"):
                    new_input = TextInput(name=f"benchmark-filter-{row_index}-value", value="", placeholder="Search...", width=90, margin=(0, 2, 0, 0))
                    input_container.objects = [new_input]
                    current_input_ref["widget"] = new_input
                    current_input_ref["is_panel"] = False
                    new_input.on_change("value", lambda attr, old, new: self.refresh())
                else:
                    searchable_input = make_searchable_input(category, col_name)
                    input_container.objects = [searchable_input]
                    current_input_ref["widget"] = searchable_input
                    current_input_ref["is_panel"] = True
            else:
                _enc_scale = self.filter_encode.get(col_name)
                _step = 1.0 / _enc_scale if _enc_scale else 1
                new_input = Spinner(name=f"benchmark-filter-{row_index}-value", value=None, step=_step, placeholder="Value...", width=90, margin=(0, 2, 0, 0))
                input_container.objects = [new_input]
                current_input_ref["widget"] = new_input
                current_input_ref["is_panel"] = False

                def _on_spinner_change(attr, old, new):
                    if row_data.get("threshold_span") is not None and new is not None:
                        import math as _m

                        if row_data.get("log_mode") and new > 0:
                            row_data["threshold_span"].location = _m.log10(new)
                        else:
                            row_data["threshold_span"].location = new
                    self.refresh()

                new_input.on_change("value", _on_spinner_change)

            # --- synchronous: clear inset immediately, show spinner if rebuilding ---
            had_inset = bool(hist_container.objects)
            if had_inset:
                # Cancel any loading overlay owned by a superseded histogram
                # rebuild before replacing its contents with this row spinner.
                hist_container.loading = False
                _spinner_html = pn.pane.HTML(
                    '<div style="display:flex;align-items:center;gap:8px;padding:8px 0;">'
                    '<div style="width:18px;height:18px;border:2px solid #ddd;border-top:2px solid #888;'
                    'border-radius:50%;animation:spin 0.8s linear infinite;"></div>'
                    '<span style="font-size:11px;color:#888;font-style:italic;">Loading...</span>'
                    "</div>"
                    "<style>@keyframes spin{to{transform:rotate(360deg)}}</style>",
                    sizing_mode="stretch_width",
                    margin=(2, 5, 0, 5),
                )
                hist_container.objects = [_spinner_html]
                row_data["histogram_pane"] = None
                row_data["histogram_fig"] = None
                row_data["threshold_span"] = None
                row_data["treemap_pane"] = None
                row_data["bridge_input"] = None

            schedule_refresh()

            # --- deferred: DB queries for distinct values + inset rebuild ---
            row_data["loading_gen"] += 1
            gen = row_data["loading_gen"]

            def _deferred_update():
                if (
                    row_data["loading_gen"] != gen
                    or category_select.value != category
                    or subcategory_select.value != col_name
                ):
                    return
                if had_inset and row_data["loading_gen"] == gen:
                    if is_text:
                        result = self.visualizations.build_text_treemap(
                            row_data, category, col_name, current_input_ref, hist_container
                        )
                        if result:
                            hist_container.objects = result
                        else:
                            hist_container.objects = []
                    elif not col_info.get("is_bool"):
                        result = self.visualizations.build_numeric_histogram(
                            row_data, category, col_name, current_input_ref["widget"]
                        )
                        if result:
                            hist_container.objects = result
                        else:
                            hist_container.objects = []
                    else:
                        hist_container.objects = []

            curdoc().add_next_tick_callback(_deferred_update)

        def update_subcategories(attr, old, new):
            """Update column options when category changes.

            Uses the startup-time self.columns cache — no self.filtering_metadata
            re-derivation on every change.  Clamps the selected value to the first
            valid column so stale/invalid selections are never left behind.
            """
            raw_columns = self.raw_columns.get(new, [])
            columns = raw_columns if raw_columns else ["No columns"]
            new_options = self.columns.get(new) or [(c, c) for c in columns]
            if subcategory_select.options != new_options:
                subcategory_select.options = new_options
            # Clamp: reset to first column when current value is not valid for new category.
            if subcategory_select.value not in set(columns):
                subcategory_select.value = columns[0]
                # The value callback owns the input rebuild in this branch.
                return
            # The same metric name can exist in both categories, so its value
            # callback will not fire. Rebuild it explicitly with new metadata.
            update_input_widget(subcategory_select.value)

        def update_input_on_column_change(attr, old, new):
            """Update input widget when column changes."""
            update_input_widget(new)

        def update_input_on_operator_change(attr, old, new):
            """Swap input widget between TextInput and SearchableSelect based on operator."""

            if input_update_state["active"]:
                return

            category = category_select.value
            col_name = subcategory_select.value
            col_info = self.filtering_metadata.get(category, {}).get("columns", {}).get(col_name, {})
            is_text = col_info.get("type") == "text"

            if is_text:
                if new in ("has", "has not"):
                    new_input = TextInput(name=f"benchmark-filter-{row_index}-value", value="", placeholder="Search...", width=90, margin=(0, 2, 0, 0))
                    input_container.objects = [new_input]
                    current_input_ref["widget"] = new_input
                    current_input_ref["is_panel"] = False
                    new_input.on_change("value", lambda attr, old, new: self.refresh())
                else:
                    searchable_input = make_searchable_input(category, col_name)
                    input_container.objects = [searchable_input]
                    current_input_ref["widget"] = searchable_input
                    current_input_ref["is_panel"] = True

            schedule_refresh()

        category_select.on_change("value", update_subcategories)
        subcategory_select.on_change("value", update_input_on_column_change)
        comparison_select.on_change("value", update_input_on_operator_change)

        dist_toggle = pn.widgets.Button(
            name="\U0001f50d",
            width=40,
            height=30,
            margin=(0, 10, 0, 0),
            description="See distribution of values",
            stylesheets=[panel_stylesheet(self.stylesheet)],
            css_classes=[f"benchmark-filter-{row_index}-lookup"],
        )

        # Remove button (Panel button for proper dynamic event handling)
        minus_btn = pn.widgets.Button(
            name="−",
            width=30,
            height=30,
            margin=(2, 5, 0, 0),
            stylesheets=[panel_stylesheet(self.stylesheet)],
        )

        query_row = pn.Row(
            category_select,
            subcategory_select,
            comparison_select,
            input_container,
            dist_toggle,
            sizing_mode="stretch_width",
            margin=(3, 0, 3, 0),
            css_classes=["control-row"],
            stylesheets=[panel_stylesheet(self.stylesheet)],
        )

        hist_container = pn.Column(sizing_mode="stretch_width", margin=(0, 5, 0, 0))
        row_wrapper = pn.Column(query_row, hist_container, sizing_mode="stretch_width", margin=(0, 0, 2, 0))

        # Store reference to this row
        row_data = {
            "query_row": query_row,
            "row_wrapper": row_wrapper,
            "hist_container": hist_container,
            "dist_toggle": dist_toggle,
            "category_select": category_select,
            "subcategory_select": subcategory_select,
            "comparison_select": comparison_select,
            "input_ref": current_input_ref,
            "minus_btn": minus_btn,
            "and_div": None,
            "histogram_pane": None,
            "histogram_fig": None,
            "threshold_span": None,
            "treemap_pane": None,
            "bridge_input": None,
            "log_mode": False,
            "log_y": False,
            "loading_gen": 0,
        }

        minus_btn.on_click(lambda _event: self._remove_row(row_data))
        dist_toggle.on_click(lambda _event: self._toggle_distribution(row_data))

        return row_data
