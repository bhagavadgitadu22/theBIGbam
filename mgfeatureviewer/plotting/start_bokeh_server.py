import argparse
import os
import sqlite3
import traceback

from bokeh.layouts import column, row
from bokeh.models import Div, InlineStyleSheet, Tooltip
from bokeh.models.callbacks import CustomJS
from bokeh.models.widgets import AutocompleteInput, CheckboxGroup, HelpButton, Button, RadioButtonGroup, CheckboxButtonGroup, RangeSlider, Select, TextInput
from bokeh.models.plots import GridPlot

# Import the plotting function from the repo
from .plotting_data_per_sample import generate_bokeh_plot_per_sample
from .plotting_data_all_samples import generate_bokeh_plot_all_samples

def build_controls(conn):
    """Query DB and return widgets and helper mappings."""
    cur = conn.cursor()

    # Widget Selector for Contigs (autocomplete with max 20 suggestions)
    cur.execute("SELECT Contig_name, Contig_length FROM Contig ORDER BY Contig_name")
    rows = cur.fetchall()
    contigs = [r[0] for r in rows]
    contig_lengths = {r[0]: r[1] for r in rows}  # Dictionary mapping contig_name -> length
    
    # If only one contig, pre-fill and disable the field
    contig_select = AutocompleteInput(value=contigs[0] if len(contigs) == 1 else "", 
                                      completions=contigs, 
                                      min_characters=0,
                                      case_sensitive=False,
                                      restrict=False,
                                      max_completions=20,
                                      placeholder="Type to search contigs...",
                                      sizing_mode="stretch_width",
                                      disabled=len(contigs) == 1)
    if len(contigs) == 1:
        contig_select.styles = {'background-color': '#e0e0e0'}

    # Widget Selector for Samples (autocomplete with max 20 suggestions)
    cur.execute("SELECT Sample_name FROM Sample ORDER BY Sample_name")
    samples = [r[0] for r in cur.fetchall()]
    
    # If only one sample, pre-fill and disable the field
    sample_select = AutocompleteInput(value=samples[0] if len(samples) == 1 else "", 
                                      completions=samples,
                                      min_characters=0,
                                      case_sensitive=False,
                                      restrict=False,
                                      max_completions=20,
                                      placeholder="Type to search samples...",
                                      sizing_mode="stretch_width",
                                      disabled=len(samples) == 1)
    if len(samples) == 1:
        sample_select.styles = {'background-color': '#e0e0e0'}
    
    # Build presence mappings: sample -> contigs and contig -> samples
    cur.execute("""
    SELECT Contig.Contig_name, Sample.Sample_name FROM Presences
      JOIN Contig ON Presences.Contig_id = Contig.Contig_id
      JOIN Sample ON Presences.Sample_id = Sample.Sample_id
    """)
    sample_to_contigs = {}
    contig_to_samples = {}
    for contig_name, sample_name in cur.fetchall():
        sample_to_contigs.setdefault(sample_name, set()).add(contig_name)
        contig_to_samples.setdefault(contig_name, set()).add(sample_name)

    # Modules and variables
    cur.execute("SELECT DISTINCT Variable_name FROM Variable")
    variables = [r[0] for r in cur.fetchall()]

    cur.execute("SELECT DISTINCT Module FROM Variable")
    modules = [r[0] for r in cur.fetchall()]

    # For each module get variables
    module_names = []
    module_widgets = []
    variables_widgets = []
    helps_widgets = []
    for module in modules:
        module_names.append(module)
        cur.execute("SELECT DISTINCT Subplot, Help FROM Variable WHERE Module=?", (module,))
        variables_checkbox = [r[0] for r in cur.fetchall()]

        if len(variables_checkbox) > 1:
            module_checkbox = CheckboxGroup(labels=[module], active=[])
            module_widgets.append(module_checkbox)
        else:
            module_widgets.append(None)

        # use a CheckboxButtonGroup for selecting individual variables in this module
        cbg = CheckboxButtonGroup(labels=variables_checkbox, active=[], sizing_mode="stretch_width", orientation="vertical")
        variables_widgets.append(cbg)

        # Consolidate help texts for the module into a single HelpButton attached to module title
        combined_help = ""
        cur.execute("SELECT DISTINCT Subplot, Title, Help FROM Variable WHERE Module=?", (module,))
        records = cur.fetchall()
        for subplot, title, help_text in records:
            if help_text is None or not help_text.strip():
                continue  # Skip empty or None help texts
            combined_help += f"{title} ({subplot} subplot): {help_text}\n"

        if combined_help:
            tooltip = Tooltip(content=combined_help, position="right")
        else:
            tooltip = None
        helps_widgets.append(tooltip)

    widgets = {
        'sample_select': sample_select,
        'contig_select': contig_select,
        'sample_to_contigs': sample_to_contigs,
        'contig_to_samples': contig_to_samples,
        'module_names': module_names,
        'module_widgets': module_widgets,
        'helps_widgets': helps_widgets,
        'variables_widgets': variables_widgets,
        'contigs': contigs,
        'contig_lengths': contig_lengths,
        'samples': samples,
        'variables': variables,
        'contig_originally_disabled': len(contigs) == 1,
        'sample_originally_disabled': len(samples) == 1
    }
    return widgets

def modify_doc_factory(db_path):
    """Return a modify_doc(doc) function to be used by Bokeh server application."""

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
    
    ## Helper functions to refresh completions based on filters
    def query_summary_filtered_items(entity_type, var_name, comparison, threshold):
        """Query Summary table for filtered contig or sample names.
        
        Args:
            entity_type: "Contig" or "Sample"
            var_name: Variable name to filter by
            comparison: ">" or "<"
            threshold: Row count threshold
        
        Returns:
            Set of matching entity names
        """
        operator = ">" if comparison == ">" else "<"
        query = f"""
            SELECT DISTINCT {entity_type}.{entity_type}_name 
            FROM Summary
            JOIN {entity_type} ON Summary.{entity_type}_id = {entity_type}.{entity_type}_id
            JOIN Variable ON Summary.Variable_id = Variable.Variable_id
            WHERE Variable.Variable_name = ? AND Summary.Row_count {operator} ?
        """
        cur = conn.cursor()
        cur.execute(query, (var_name, threshold))
        return {row[0] for row in cur.fetchall()}
    
    def get_variable_filtered_contigs():
        """Apply variable-based filters to get allowed contigs."""
        if not variable_filter_rows:
            return set(orig_contigs)
        
        allowed_contigs = set(orig_contigs)
        
        for filter_row in variable_filter_rows:
            var_name = filter_row.children[0].value
            comparison = filter_row.children[1].value
            threshold_str = filter_row.children[2].value
            
            # Skip invalid filters - only apply if variable name is valid and exists
            if not var_name or not threshold_str or var_name not in widgets['variables']:
                continue
            
            try:
                threshold = int(threshold_str)
                if threshold < 0:
                    continue
            except ValueError:
                continue
            
            matching = query_summary_filtered_items("Contig", var_name, comparison, threshold)
            allowed_contigs &= matching  # AND logic: all filters must match
        
        return allowed_contigs
    
    def get_variable_filtered_samples():
        """Apply variable-based filters to get allowed samples."""
        if not variable_filter_rows:
            return set(orig_samples)
        
        allowed_samples = set(orig_samples)
        
        for filter_row in variable_filter_rows:
            var_name = filter_row.children[0].value
            comparison = filter_row.children[1].value
            threshold_str = filter_row.children[2].value
            
            # Skip invalid filters - only apply if variable name is valid and exists
            if not var_name or not threshold_str or var_name not in widgets['variables']:
                continue
            
            try:
                threshold = int(threshold_str)
                if threshold < 0:
                    continue
            except ValueError:
                continue
            
            matching = query_summary_filtered_items("Sample", var_name, comparison, threshold)
            allowed_samples &= matching  # AND logic: all filters must match
        
        return allowed_samples
    
    def update_widget_completions(widget, completions, originally_disabled=False):
        """Update widget completions. Auto-fill and disable when filtering reduces to 1 option."""
        widget.completions = completions
        
        if len(completions) == 1:
            # Only one option after filtering: fill it, disable, and gray background
            widget.value = completions[0]
            widget.disabled = True
            widget.styles = {'background-color': '#e0e0e0'}
        else:
            # Multiple options: re-enable (unless originally disabled), restore normal background
            if not originally_disabled:
                widget.disabled = False
                widget.styles = {}
            # Clear value if it's not in the new completions
            if widget.value not in completions:
                widget.value = ""
    
    def refresh_contig_options():
        # Start with presence filter if active
        if 0 in filter_contigs.active and views.active == 0:
            sel_sample = widgets['sample_select'].value
            allowed = widgets['sample_to_contigs'].get(sel_sample, set())
            completions = [c for c in orig_contigs if c in allowed]
        else:
            completions = list(orig_contigs)

        # Apply length filter
        if length_slider is not None:
            min_length, max_length = length_slider.value
            completions = [c for c in completions if min_length <= widgets['contig_lengths'].get(c, 0) <= max_length]
        
        # Apply variable-based filters (only in "Per sample" view)
        if views.active == 0:
            var_allowed = get_variable_filtered_contigs()
            completions = [c for c in completions if c in var_allowed]

        update_widget_completions(widgets['contig_select'], completions, widgets['contig_originally_disabled'])

    def refresh_sample_options():
        # Start with presence filter if active
        if 0 in filter_samples.active and views.active == 0:
            sel_contig = widgets['contig_select'].value
            allowed = widgets['contig_to_samples'].get(sel_contig, set())
            completions = [s for s in orig_samples if s in allowed]
        else:
            completions = list(orig_samples)
        
        # Apply variable-based filters (only in "Per sample" view)
        if views.active == 0:
            var_allowed = get_variable_filtered_samples()
            completions = [s for s in completions if s in var_allowed]

        update_widget_completions(widgets['sample_select'], completions, widgets['sample_originally_disabled'])

    ## Views function
    # Enforce single-variable selection when in "All samples" view
    def make_global_variable_callback(cbg):
        def callback(attr, old, new):
            if global_toggle_lock['locked']:
                return
            # Only enforce in All samples mode
            if views.active != 1:
                return

            # Determine which index was most-recently changed.
            # Use set difference to find an added index (preferred).
            sel_index = None
            old_set = set(old) if old else set()
            new_set = set(new) if new else set()
            added = new_set - old_set

            if added:
                # pick the (one) newly added index
                sel_index = next(iter(added))
            elif new:
                # no clear addition, fall back to last element
                sel_index = new[-1]

            global_toggle_lock['locked'] = True
            for other in widgets['variables_widgets']:
                if other is cbg:
                    if sel_index is None:
                        other.active = []
                    else:
                        other.active = [sel_index]
                else:
                    other.active = []
            global_toggle_lock['locked'] = False
        return callback
    
    # Views (One sample / All samples) callback: show/hide sample-related controls
    def on_view_change(attr, old, new):
        is_all = (new == 1)  # True means All samples
        
        # Lock callbacks during view change to prevent cascading updates
        global_toggle_lock['locked'] = True
        
        # Toggle Sample section
        separator_samples.visible = not is_all
        sample_header.visible = not is_all
        sample_content.visible = not is_all
        filter_samples.visible = not is_all
        widgets['sample_select'].visible = not is_all

        # Contig section: keep visible but hide filter checkbox and per-variable filtering in All samples mode
        # separator_contigs always visible
        # contig_header always visible
        # contig_content always visible
        filter_contigs.visible = not is_all
        per_variable_header.visible = not is_all
        variable_filters_column.visible = not is_all
        
        if is_all:
            filter_contigs.active = []
            filter_samples.active = []
            # Hide module checkboxes and clear all variable selections
            for mw in widgets['module_widgets']:
                if mw is not None:
                    mw.visible = False
                    mw.active = []
            for cbg in widgets['variables_widgets']:
                cbg.active = []
        else:
            # Show module checkboxes
            for mw in widgets['module_widgets']:
                if mw is not None:
                    mw.visible = True
        
        # Toggle header visibility: checkbox-header in One-sample, title-header in All-samples
        for i, mw in enumerate(widgets['module_widgets']):
            hdr_cb = module_header_with_checkbox[i]
            hdr_title = module_header_with_title[i]
            
            if mw is None:
                # Single-variable module: always show title header
                if hdr_cb:
                    hdr_cb.visible = False
                if hdr_title:
                    hdr_title.visible = True
            else:
                # Module with checkbox: toggle between headers
                hdr_cb.visible = not is_all
                hdr_title.visible = is_all
        
        # Unlock callbacks
        global_toggle_lock['locked'] = False
        
        # Refresh options after view change
        # In "All samples": clears presence and per-variable filters, keeps only length filter
        # In "Per sample": applies all filters including presence and per-variable
        refresh_contig_options()
        if not is_all:
            refresh_sample_options()

    def create_variable_filter_row():
        """Create a new row of variable filter widgets."""
        var_input = AutocompleteInput(
            completions=widgets['variables'],
            placeholder="Select variable...",
            min_characters=0,
            case_sensitive=False,
            restrict=False,
            max_completions=20,
            sizing_mode="stretch_width"
        )
        
        comparison_select = Select(
            options=[">", "<"],
            value=">",
            width=45
        )
        
        threshold_input = TextInput(
            value="0",
            placeholder="Threshold",
            width=60
        )
        
        plus_btn = Button(label="+", width=30, height=30)
        # len(variable_filter_rows)>0 to make - button invisible initially
        minus_btn = Button(label="−", width=30, height=30, visible=len(variable_filter_rows)>0)
        
        filter_row = row(var_input, comparison_select, threshold_input, plus_btn, minus_btn, sizing_mode="stretch_width")
        
        def add_row_callback():
            new_row = create_variable_filter_row()
            variable_filter_rows.append(new_row)
            variable_filters_column.children = list(variable_filters_column.children) + [new_row]
            # Make minus buttons visible when there's more than one row
            for row_widget in variable_filter_rows:
                row_widget.children[-1].visible = True  # Last child is minus button
        
        def remove_row_callback():
            if filter_row in variable_filter_rows:
                variable_filter_rows.remove(filter_row)
                variable_filters_column.children = [r for r in variable_filters_column.children if r != filter_row]
                # Hide minus buttons if only one row remains
                if len(variable_filter_rows) == 1:
                    variable_filter_rows[0].children[-1].visible = False
                    
                # Refresh options after removing filter
                refresh_contig_options()
                refresh_sample_options()
        
        plus_btn.on_click(add_row_callback)
        minus_btn.on_click(remove_row_callback)
        
        # Create a shared callback that refreshes both contig and sample options
        def refresh_on_filter_change(attr, old, new):
            refresh_contig_options()
            refresh_sample_options()
        
        # Attach to all three inputs
        var_input.on_change('value', refresh_on_filter_change)
        comparison_select.on_change('value', refresh_on_filter_change)
        threshold_input.on_change('value', refresh_on_filter_change)

        return filter_row
    
    ## Apply button function
    def apply_clicked():
        try:
            sample = widgets['sample_select'].value
            contig = widgets['contig_select'].value

            # Check if gene map should be shown
            genbank_path = db_path if (0 in show_genemap.active) else None

            # Build requested_features list
            requested_features = []
            for cbg in widgets['variables_widgets']:
                # cbg.active is a list of indices
                for idx in cbg.active:
                    requested_features.append(cbg.labels[idx])

            # Save current positions of the plot for restoration after re-plot
            fig = main_placeholder.children[0]
            
            xstart = None
            xend = None
            if isinstance(fig, GridPlot):
                subplot = fig.children[0][0]
                xstart = subplot.x_range.start
                xend = subplot.x_range.end

            if views.active == 1:
                # All-samples view: require exactly one variable selected across all modules
                selected_var = None
                for cbg in widgets['variables_widgets']:
                    if cbg.active:
                        # pick the last selected index in this group
                        selected_var = cbg.labels[cbg.active[-1]]
                        break

                if not selected_var:
                    raise ValueError("When in 'All samples' view you must select one variable to plot.")

                print(f"[start_bokeh_server] Generating plot for all samples with variable={selected_var}, contig={contig}")
                grid = generate_bokeh_plot_all_samples(conn, selected_var, contig, xstart=xstart, xend=xend, genbank_path=genbank_path)
            else:
                # One-sample view: collect possibly-many requested features and call per-sample plot
                requested_features = []
                for cbg in widgets['variables_widgets']:
                    for idx in cbg.active:
                        requested_features.append(cbg.labels[idx])

                print(f"[start_bokeh_server] Generating plot for sample={sample}, contig={contig}, features={requested_features}")
                grid = generate_bokeh_plot_per_sample(conn, requested_features, contig, sample, xstart=xstart, xend=xend, genbank_path=genbank_path)

            main_placeholder.children = [grid]

        except Exception as e:
            tb = traceback.format_exc()
            print(f"[start_bokeh_server] Exception: {tb}", flush=True)
            main_placeholder.children = [Div(text=f"<pre>Error building plot:\n{tb}</pre>")]



    ### Creating all DOM elements
    # Open SQLite database connection to build widgets depending on data
    conn = sqlite3.connect(db_path)
    widgets = build_controls(conn)

    # Load the CSS
    css_path = os.path.join(os.path.dirname(__file__), "..", "..", "static", "bokeh_styles.css")
    with open(css_path) as f:
        css_text = f.read()
    stylesheet = InlineStyleSheet(css=css_text)

    # Create main elements
    ## Views section
    views = RadioButtonGroup(labels=["One sample", "All samples"], active=0, sizing_mode="stretch_width", button_type="primary")

    # Global lock for toggles when enforcing "All samples" view (single-variable mode)
    global_toggle_lock = {'locked': False}
    # Attach global variable callbacks to all CheckboxButtonGroups
    for cbg in widgets['variables_widgets']:
        cbg.on_change('active', make_global_variable_callback(cbg))
    views.on_change('active', on_view_change)


    ## Build Sample section
    sample_toggle_btn = Button(label="▼", width=20, height=20, button_type="primary", align="center", margin=0, stylesheets=[stylesheet])
    sample_title = Div(text="<b>Samples</b>", align="center")
    sample_header = row(sample_toggle_btn, sample_title, sizing_mode="stretch_width", align="center")

    filter_samples = CheckboxGroup(labels=["Only show samples present with selected contig"], active=[])
    filter_samples.on_change('active', lambda attr, old, new: refresh_sample_options())

    sample_children = [filter_samples]
    sample_content = column(
        *sample_children,
        visible=True, sizing_mode="stretch_width"
    )

    sample_toggle_btn.on_click(make_toggle_callback(sample_toggle_btn, sample_content))
    widgets['sample_select'].on_change('value', lambda attr, old, new: refresh_contig_options())


    ## Build Contig section
    # Keep original full lists so we can restore when filters are off
    orig_contigs = list(widgets['contigs'])
    orig_samples = list(widgets['samples'])

    filter_contigs = CheckboxGroup(labels=["Only show contigs present with selected sample"], active=[])
    filter_contigs.on_change('active', lambda attr, old, new: refresh_contig_options())

    # Length filter slider (only if multiple contigs)
    min_len = min(widgets['contig_lengths'].values())
    max_len = max(widgets['contig_lengths'].values())
    length_slider = None
    if len(widgets["contigs"]) > 1:
        length_slider = RangeSlider(start=min_len, end=max_len, value=(min_len, max_len), step=1, title="Contig length")
        length_slider.on_change('value', lambda attr, old, new: refresh_contig_options())

    # Create collapsible Filtering section
    contig_toggle_btn = Button(label="▼", width=20, height=20, button_type="primary", align="center", margin=0, stylesheets=[stylesheet])
    contig_toggle_btn.styles = {'padding': '0px', 'line-height': '20px'}
    contig_title = Div(text="<b>Contigs</b>", align="center")
    contig_header = row(contig_toggle_btn, contig_title, sizing_mode="stretch_width", align="center")
    contig_children = [filter_contigs]
    if length_slider is not None:
        contig_children.append(length_slider)
    
    # Add "Per variable" filtering subsection
    combined_help = "Only suggest contigs where more (>) or less (<) than a specified number of points exist for the variable in the selected sample"
    tooltip = Tooltip(content=combined_help, position="right")
    help_per_variable = HelpButton(tooltip=tooltip, width=20, height=20, align="center", button_type="light", stylesheets=[stylesheet])
    per_variable_title = Div(text="Per variable:")
    per_variable_header = row(per_variable_title, help_per_variable, sizing_mode="stretch_width")

    # Store all variable filter rows for dynamic management (already initialized above)
    variable_filters_column = column(sizing_mode="stretch_width")
    
    # Create initial filter row
    variable_filter_rows = []
    initial_row = create_variable_filter_row()
    variable_filter_rows.append(initial_row)
    variable_filters_column.children = [initial_row]
    
    contig_children.append(per_variable_header)
    contig_children.append(variable_filters_column)
    contig_content = column(
        *contig_children,
        visible=True, sizing_mode="stretch_width"
    )
    contig_toggle_btn.on_click(make_toggle_callback(contig_toggle_btn, contig_content))

    widgets['contig_select'].on_change('value', lambda attr, old, new: refresh_sample_options())


    ## Build Variables section
    variables_title = Div(text="<b>Variables</b>")
    show_genemap = CheckboxGroup(labels=["Show gene map"], active=[0])
    
    # Append variable selectors. For modules that have a module-checkbox widget we
    # show either the checkbox (One-sample) or the plain module title (All-samples)
    # We create a container that holds both variants and toggle visibility
    module_header_with_checkbox = []
    module_header_with_title = []
    # Store toggle buttons and content containers for collapsible sections
    module_toggles = []
    module_contents = []

    controls_variables = []
    for i, module_widget in enumerate(widgets['module_widgets']):
        module_name = widgets['module_names'][i]
        
        help_tooltip = widgets['helps_widgets'][i]

        # Create a small toggle button for collapsible section (just the arrow)
        # Start with modules folded (collapsed)
        toggle_btn = Button(label="▶", width=20, height=20, button_type="primary", align="center", margin=0, stylesheets=[stylesheet])
        toggle_btn.styles = {'padding': '0px', 'line-height': '20px'}
        module_toggles.append(toggle_btn)

        # Build two header variants: one with the checkbox (shows module name as label)
        # and one with the plain title (used when checkbox is hidden). Both may include the help button.
        if module_widget is not None:
            # Create separate title div for the title-only header
            module_title_div = Div(text=f"{module_name}", align="center")
            
            if help_tooltip is not None:
                # Need separate help buttons for each header to avoid "already in doc" error
                help_btn_cb = HelpButton(tooltip=help_tooltip, width=20, height=20, align="center", button_type="light", stylesheets=[stylesheet])
                help_btn_title = HelpButton(tooltip=help_tooltip, width=20, height=20, align="center", button_type="light", stylesheets=[stylesheet])
                hdr_cb = row(toggle_btn, module_widget, help_btn_cb, sizing_mode="stretch_width", align="center")
                hdr_title = row(toggle_btn, module_title_div, help_btn_title, sizing_mode="stretch_width", align="center")
            else:
                hdr_cb = row(toggle_btn, module_widget, sizing_mode="stretch_width", align="center")
                hdr_title = row(toggle_btn, module_title_div, sizing_mode="stretch_width", align="center")

            # Default: show the checkbox header, hide the plain title header
            hdr_cb.visible = True
            hdr_title.visible = False
            controls_variables.append(hdr_cb)
            controls_variables.append(hdr_title)
            
            module_header_with_checkbox.append(hdr_cb)
            module_header_with_title.append(hdr_title)
        else:
            # No module checkbox exists: show plain title (with help if available)
            module_title_div = Div(text=f"{module_name}", align="center")
            if help_btn is not None:
                hdr = row(toggle_btn, module_title_div, help_btn, sizing_mode="stretch_width", align="center")
            else:
                hdr = row(toggle_btn, module_title_div, sizing_mode="stretch_width", align="center")
            hdr.visible = True
            controls_variables.append(hdr)
            
            module_header_with_checkbox.append(None)
            module_header_with_title.append(hdr)

        # Add the module's CheckboxButtonGroup for variables (this will be collapsible)
        # Start with modules folded (collapsed)
        cbg = widgets['variables_widgets'][i]
        cbg.visible = False
        module_contents.append(cbg)
        controls_variables.append(cbg)

    # Add callbacks for collapsible sections
    for i, toggle_btn in enumerate(module_toggles):
        content = module_contents[i]
        toggle_btn.on_click(make_toggle_callback(toggle_btn, content))


    ### Attach callbacks
    for i, mc in enumerate(widgets['module_widgets']):
        if mc is None:
            continue

        toggles = widgets['variables_widgets'][i]
        lock = {"locked": False}  # per-module lock

        # Module → toggles (CheckboxButtonGroup)
        def make_module_callback(mc, toggles, lock):
            def callback(attr, old, new):
                if lock.get("locked", False):
                    return
                lock["locked"] = True
                module_on = 0 in mc.active
                if module_on:
                    toggles.active = list(range(len(toggles.labels)))
                else:
                    toggles.active = []
                lock["locked"] = False
            return callback

        mc.on_change("active", make_module_callback(mc, toggles, lock))

        # Variable → module (update module checkbox only)
        def make_variable_callback(mc, toggles, lock):
            def callback(attr, old, new):
                if lock.get("locked", False):
                    return
                total = len(toggles.labels)
                active_count = len(toggles.active)

                lock["locked"] = True
                if active_count == total and total > 0:
                    mc.active = [0]
                else:
                    mc.active = []
                lock["locked"] = False
            return callback

        toggles.on_change("active", make_variable_callback(mc, toggles, lock))


    ## Create final Apply button
    apply_button = Button(label="Apply", button_type="primary", align="center")
    apply_button.on_click(lambda: apply_clicked())


    ## Put together all DOM elements
    # Create visual separators (horizontal lines) using background color
    separator_samples = Div(text="", height=1, width=350, styles={'background-color': '#333', 'margin': '10px 0'})
    separator_contigs = Div(text="", height=1, width=350, styles={'background-color': '#333', 'margin': '10px 0'})
    separator_variables = Div(text="", height=1, width=350, styles={'background-color': '#333', 'margin': '10px 0'})

    controls_children = [views, separator_samples, sample_header, sample_content, widgets['sample_select'], separator_contigs, contig_header, contig_content, widgets['contig_select'], separator_variables, variables_title, show_genemap] + controls_variables + [apply_button]

    controls_column = column(*controls_children, width=350, sizing_mode="stretch_height", spacing=0)
    controls_column.css_classes = ["left-col"]

    main_placeholder = column(Div(text="<i>No plot yet. Select one sample, one contig and at least one variable in \"One sample\" mode or one contig and one variable in \"All samples\" mode and click Apply.</i>"), sizing_mode="stretch_both")

    # Wrap everything in a Flex container
    layout = row(controls_column, main_placeholder, sizing_mode="stretch_both", spacing = 0)
    layout.stylesheets = [stylesheet]

    def modify_doc(doc):
        doc.add_root(layout)

    return modify_doc

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", required=True, help="Path to sqlite DB")
    parser.add_argument("--port", type=int, default=5006, help="Port to serve Bokeh app")
    args = parser.parse_args()

    modify_doc = modify_doc_factory(args.db)

    # Start bokeh server programmatically
    try:
        from bokeh.server.server import Server
        from bokeh.application import Application
        from bokeh.application.handlers.function import FunctionHandler
    except Exception as e:
        print("Bokeh server components not installed or unavailable:", e)
        print("You can run this file with `bokeh serve --show start_bokeh_server.py --args --db <DB>` as alternative.")
        return

    app = Application(FunctionHandler(modify_doc))
    server = Server({'/': app}, port=args.port)
    server.start()
    print(f"Bokeh server running at http://localhost:{args.port}/")
    server.io_loop.start()

def add_serve_args(parser):
    parser.add_argument("--db", required=True, help="Path to sqlite DB")
    parser.add_argument("--port", type=int, default=5006, help="Port to serve Bokeh app")

def run_serve(args):
    modify_doc = modify_doc_factory(args.db)

    try:
        from bokeh.server.server import Server
        from bokeh.application import Application
        from bokeh.application.handlers.function import FunctionHandler
    except Exception as e:
        print("Bokeh server components not installed or unavailable:", e)
        print("You can run this file with `bokeh serve --show start_bokeh_server.py --args --db <DB>` as alternative.")
        return 1

    app = Application(FunctionHandler(modify_doc))
    server = Server({'/': app}, port=args.port)
    server.start()
    print(f"Bokeh server running at http://localhost:{args.port}/")
    server.io_loop.start()
    return 0

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    add_serve_args(parser)
    args = parser.parse_args()
    raise SystemExit(run_serve(args))