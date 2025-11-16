import argparse
import os
import sqlite3
import traceback

from bokeh.layouts import column, row
from bokeh.models import Div, InlineStyleSheet, Tooltip
from bokeh.models.widgets import Select, CheckboxGroup, HelpButton, Button, RadioButtonGroup, CheckboxButtonGroup

# Import the plotting function from the repo
from plotting_data_per_sample import generate_bokeh_plot_per_sample
from plotting_data_all_samples import generate_bokeh_plot_all_samples
from plotting_data_all_samples import generate_bokeh_plot_all_samples

def build_controls(conn):
    """Query DB and return widgets and helper mappings."""
    cur = conn.cursor()

    # Widget Selector for Contigs
    cur.execute("SELECT Contig_name FROM Contig ORDER BY Contig_name")
    contigs = [r[0] for r in cur.fetchall()]
    contig_select = Select(value=contigs[0], options=contigs, sizing_mode="stretch_width")

    # Widget Selector for Samples
    cur.execute("SELECT Sample_name FROM Sample ORDER BY Sample_name")
    samples = [r[0] for r in cur.fetchall()]
    sample_select = Select(value=samples[0], options=samples, sizing_mode="stretch_width")

    # Modules and variables
    cur.execute("SELECT DISTINCT Module FROM Variable")
    modules = [r[0] for r in cur.fetchall()]

    # For each module get variables
    module_widgets = []
    variables_widgets = []
    helps_widgets = []
    module_names = []
    for module in modules:
        module_names.append(module)
        cur.execute("SELECT DISTINCT Subplot, Help FROM Variable WHERE Module=?", (module,))
        records = cur.fetchall()
        variables_checkbox = [r[0] for r in records]
        helps_checkbox = [r[1] for r in records]

        if len(variables_checkbox) > 1:
            module_checkbox = CheckboxGroup(labels=[module], active=[])
            module_widgets.append(module_checkbox)
        else:
            module_widgets.append(None)

        # Consolidate help texts for the module into a single HelpButton attached to module title
        combined_help = "\n".join([f"{lbl}: {ht}" for lbl, ht in zip(variables_checkbox, helps_checkbox) if ht and ht.strip() != ""])
        if combined_help:
            tooltip = Tooltip(content=combined_help, position="right")
            help_button = HelpButton(tooltip=tooltip, width=30, height=30, align="center")
        else:
            help_button = None
        helps_widgets.append(help_button)

        # use a CheckboxButtonGroup for selecting individual variables in this module
        cbg = CheckboxButtonGroup(labels=variables_checkbox, active=[], sizing_mode="stretch_width", orientation="vertical")
        variables_widgets.append(cbg)

    widgets = {
        'sample_select': sample_select,
        'contig_select': contig_select,
        'module_names': module_names,
        'module_widgets': module_widgets,
        'helps_widgets': helps_widgets,
        'variables_widgets': variables_widgets
    }
    return widgets

def modify_doc_factory(db_path):
    """Return a modify_doc(doc) function to be used by Bokeh server application."""
    # Load the CSS
    css_path = os.path.join(os.path.dirname(__file__), "static", "bokeh_styles.css")
    with open(css_path) as f:
        css_text = f.read()
    stylesheet = InlineStyleSheet(css=css_text)

    instructions = Div(text="<b>Select elements to plot and click Apply:</b>")

    views_title = Div(text="<b>View</b>")
    views = RadioButtonGroup(labels=["One sample", "All samples"], active=0, sizing_mode="stretch_width")

    conn = sqlite3.connect(db_path)
    widgets = build_controls(conn)

    # Header placeholders (populated below) so view-change callback can toggle them
    header_with_checkbox = []
    header_with_title = []

    # Build presence mappings: sample -> contigs and contig -> samples
    presence_cur = conn.cursor()
    presence_cur.execute("""
    SELECT Contig.Contig_name, Sample.Sample_name FROM Presences
      JOIN Contig ON Presences.Contig_id = Contig.Contig_id
      JOIN Sample ON Presences.Sample_id = Sample.Sample_id
    """)
    sample_to_contigs = {}
    contig_to_samples = {}
    for contig_name, sample_name in presence_cur.fetchall():
        sample_to_contigs.setdefault(sample_name, set()).add(contig_name)
        contig_to_samples.setdefault(contig_name, set()).add(sample_name)

    # Keep original full lists so we can restore when filters are off
    orig_contigs = list(widgets['contig_select'].options)
    orig_samples = list(widgets['sample_select'].options)

    contigs_title = Div(text="<b>Contig</b>")
    filter_contigs = CheckboxGroup(labels=["Only show contigs present with selected sample"], active=[])

    samples_title = Div(text="<b>Sample</b>")
    filter_samples = CheckboxGroup(labels=["Only show samples present with selected contig"], active=[])

    # Helper functions to refresh options based on filters
    def refresh_contig_options():
        if 0 in filter_contigs.active and views.active == 0:
            sel_sample = widgets['sample_select'].value
            allowed = sample_to_contigs.get(sel_sample, set())
            options = [c for c in orig_contigs if c in allowed]
        else:
            options = list(orig_contigs)

        widgets['contig_select'].options = options
        if widgets['contig_select'].value not in options:
            widgets['contig_select'].value = options[0] if options else ""

    def refresh_sample_options():
        if 0 in filter_samples.active and views.active == 0:
            sel_contig = widgets['contig_select'].value
            allowed = contig_to_samples.get(sel_contig, set())
            options = [s for s in orig_samples if s in allowed]
        else:
            options = list(orig_samples)

        widgets['sample_select'].options = options
        if widgets['sample_select'].value not in options:
            widgets['sample_select'].value = options[0] if options else ""

    # Global lock for toggles when enforcing single-variable mode
    global_toggle_lock = {'locked': False}

    # Enforce single-variable selection when in "All samples" view
    def make_global_variable_callback(cbg):
        def callback(attr, old, new):
            if global_toggle_lock['locked']:
                return
            # Only enforce in All samples mode and when some index was selected
            if views.active == 1 and new:
                # pick the last selected index
                sel_index = new[-1]
                global_toggle_lock['locked'] = True
                for other in widgets['variables_widgets']:
                    if other is cbg:
                        # ensure current cbg only has sel_index
                        other.active = [sel_index]
                    else:
                        other.active = []
                global_toggle_lock['locked'] = False
        return callback

    # Attach global variable callbacks to all CheckboxButtonGroups
    for cbg in widgets['variables_widgets']:
        cbg.on_change('active', make_global_variable_callback(cbg))

    # Views (One sample / All samples) callback: show/hide sample-related controls
    def on_view_change(attr, old, new):
        # new == 1 means All samples
        is_all = (new == 1)
        widgets['sample_select'].visible = not is_all
        samples_title.visible = not is_all
        # hide/deactivate presence filters when All samples
        if is_all:
            # deactivate and hide the "only show contigs present with selected sample" checkbox
            filter_contigs.active = []
            filter_contigs.visible = False
            # hide sample filter as well
            filter_samples.active = []
            filter_samples.visible = False
            # hide module checkboxes but keep module titles
            for mw in widgets['module_widgets']:
                if mw is not None:
                    mw.visible = False
                    mw.active = []
            # Clear all variable toggles when switching to All-samples view
            for cbg in widgets['variables_widgets']:
                cbg.active = []
        else:
            filter_contigs.visible = True
            filter_samples.visible = True
            # show module checkboxes again
            for mw in widgets['module_widgets']:
                if mw is not None:
                    mw.visible = True
                    mw.active = mw.active
        # When switching mode, refresh options to ensure consistency
        # Toggle header visibility: show checkbox-header in One-sample, title-header in All-samples
        for i, mw in enumerate(widgets['module_widgets']):
            hdr_cb = header_with_checkbox[i] if i < len(header_with_checkbox) else None
            hdr_title = header_with_title[i] if i < len(header_with_title) else None
            if mw is None:
                # single-variable module: ensure title header visible
                if hdr_cb is not None:
                    hdr_cb.visible = False
                if hdr_title is not None:
                    hdr_title.visible = True
            else:
                # module with checkbox: toggle which header is visible
                if is_all:
                    if hdr_cb is not None:
                        hdr_cb.visible = False
                    if hdr_title is not None:
                        hdr_title.visible = True
                else:
                    if hdr_cb is not None:
                        hdr_cb.visible = True
                    if hdr_title is not None:
                        hdr_title.visible = False
        # When switching mode, refresh options to ensure consistency
        refresh_contig_options()
        refresh_sample_options()

    views.on_change('active', on_view_change)

    # Wire up select/filter interactions
    widgets['sample_select'].on_change('value', lambda attr, old, new: refresh_contig_options())
    widgets['contig_select'].on_change('value', lambda attr, old, new: refresh_sample_options())
    filter_contigs.on_change('active', lambda attr, old, new: refresh_contig_options())
    filter_samples.on_change('active', lambda attr, old, new: refresh_sample_options())

    variables_title = Div(text="<b>Variables</b>")
    controls_children = [instructions, views_title, views, contigs_title, widgets['contig_select'], filter_contigs, samples_title, widgets['sample_select'], filter_samples, variables_title]
    
    # Append variable selectors. For modules that have a module-checkbox widget we
    # show either the checkbox (One-sample) or the plain module title (All-samples).
    header_with_checkbox = []
    header_with_title = []
    for i, module_widget in enumerate(widgets['module_widgets']):
        module_name = widgets['module_names'][i]
        module_title_div = Div(text=f"{module_name}")

        help_btn = widgets['helps_widgets'][i]

        # Build two header variants: one with the checkbox (shows module name as label)
        # and one with the plain title (used when checkbox is hidden). Both may include the help button.
        if module_widget is not None:
            if help_btn is not None:
                hdr_cb = row(module_widget, help_btn, sizing_mode="stretch_width")
                hdr_title = row(module_title_div, help_btn, sizing_mode="stretch_width")
            else:
                hdr_cb = module_widget
                hdr_title = module_title_div

            # Default: show the checkbox header, hide the plain title header
            hdr_cb.visible = True
            hdr_title.visible = False
            controls_children.append(hdr_cb)
            controls_children.append(hdr_title)
        else:
            # No module checkbox exists: show plain title (with help if available)
            if help_btn is not None:
                hdr = row(module_title_div, help_btn, sizing_mode="stretch_width"
                          )
            else:
                hdr = module_title_div
            hdr.visible = True
            controls_children.append(hdr)

        header_with_checkbox.append(hdr_cb if module_widget is not None else None)
        header_with_title.append(hdr_title if module_widget is not None else hdr)

        # Add the module's CheckboxButtonGroup for variables
        cbg = widgets['variables_widgets'][i]
        controls_children.append(cbg)

    apply_button = Button(label="Apply", button_type="primary", align="center")
    controls_children.append(apply_button)

    controls_column = column(*controls_children, width=350, sizing_mode="stretch_height", spacing=0)
    controls_column.css_classes = ["left-col"]

    main_placeholder = column(Div(text="<i>No plot yet. Select options and click Apply.</i>"), sizing_mode="stretch_both")

    # Wrap everything in a Flex container
    layout = row(controls_column, main_placeholder, sizing_mode="stretch_both", spacing = 0)
    layout.stylesheets = [stylesheet]

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

    def apply_clicked():
        try:
            sample = widgets['sample_select'].value
            contig = widgets['contig_select'].value

            # Build requested_features list
            requested_features = []
            for cbg in widgets['variables_widgets']:
                # cbg.active is a list of indices
                for idx in cbg.active:
                    requested_features.append(cbg.labels[idx])

            print(f"[start_bokeh_server] Generating plot for sample={sample}, contig={contig}, features={requested_features}")
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
                grid = generate_bokeh_plot_all_samples(conn, selected_var, contig)
            else:
                # One-sample view: collect possibly-many requested features and call per-sample plot
                requested_features = []
                for cbg in widgets['variables_widgets']:
                    for idx in cbg.active:
                        requested_features.append(cbg.labels[idx])

                print(f"[start_bokeh_server] Generating plot for sample={sample}, contig={contig}, features={requested_features}")
                grid = generate_bokeh_plot_per_sample(conn, requested_features, contig, sample)

            main_placeholder.children = [grid]

        except Exception as e:
            tb = traceback.format_exc()
            main_placeholder.children = [Div(text=f"<pre>Error building plot:\n{tb}</pre>")]

    apply_button.on_click(lambda: apply_clicked())

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

if __name__ == "__main__":
    main()