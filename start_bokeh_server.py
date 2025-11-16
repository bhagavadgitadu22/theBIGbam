import argparse
import os
import sqlite3
import traceback

from bokeh.layouts import column, row
from bokeh.models import Div, InlineStyleSheet, Tooltip, Toggle
from bokeh.models.widgets import Select, CheckboxGroup, HelpButton, Button, RadioButtonGroup

# Import the plotting function from the repo
from plotting_data import generate_bokeh_plot

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
    for module in modules:
        cur.execute("SELECT DISTINCT Subplot, Help FROM Variable WHERE Module=?", (module,))
        records = cur.fetchall()
        variables_checkbox = [r[0] for r in records]
        helps_checkbox = [r[1] for r in records]

        if len(variables_checkbox) > 1:
            module_checkbox = CheckboxGroup(labels=[module], active=[])
            module_widgets.append(module_checkbox)
        else:
            module_widgets.append(None)

        # use CheckboxButtonGroup for selecting individual variables
        checkboxes = []
        helps = []
        for label, help_text in zip(variables_checkbox, helps_checkbox):
            checkbox = Toggle(label=label, active=False, sizing_mode="stretch_width", height=30)
            checkboxes.append(checkbox)

            if help_text != "":
                tooltip = Tooltip(
                    content=f"{help_text}",
                    position="right"
                )
                help_button = HelpButton(tooltip=tooltip, width=30, height=30)
                helps.append(help_button)
            else:
                helps.append(None)

        # Combine into a column
        variables_widgets.append(checkboxes)
        helps_widgets.append(helps)

    apply_button = Button(label="Apply", button_type="primary", align="center")

    widgets = {
        'sample_select': sample_select,
        'contig_select': contig_select,
        'module_widgets': module_widgets,
        'variables_widgets': variables_widgets,
        'helps_widgets': helps_widgets,
        'apply_button': apply_button
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

    contigs_title = Div(text="<b>Contig</b>")
    filter_contigs = CheckboxGroup(labels=["Only show contigs present with selected sample"], active=[])

    samples_title = Div(text="<b>Sample</b>")
    filter_samples = CheckboxGroup(labels=["Only show samples present with selected contig"], active=[])

    variables_title = Div(text="<b>Variables</b>")
    controls_children = [instructions, views_title, views, contigs_title, widgets['contig_select'], filter_contigs, samples_title, widgets['sample_select'], filter_samples, variables_title]
    
    # Append variable selectors
    for i, module_widget in enumerate(widgets['module_widgets']):
        if module_widget is not None:
            controls_children.append(module_widget)

        # Get the variable widget block (no layout inside)
        checkboxes = widgets['variables_widgets'][i]
        helps = widgets['helps_widgets'][i]

        for cb, hb in zip(checkboxes, helps):
            if hb is not None:
                controls_children.append(row(cb, hb, sizing_mode="stretch_width"))
            else:
                controls_children.append(cb)

    controls_children.append(widgets['apply_button'])
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

        # Module → toggles
        def make_module_callback(mc, toggles, lock):
            def callback(attr, old, new):
                # Only act if this change comes from user
                if lock.get("locked", False):
                    return
                lock["locked"] = True
                module_on = 0 in mc.active
                for t in toggles:
                    t.active = module_on
                lock["locked"] = False
            return callback

        mc.on_change("active", make_module_callback(mc, toggles, lock))

        # Variable → module (update module checkbox only)
        def make_variable_callback(mc, toggles, lock):
            def callback(attr, old, new):
                if lock.get("locked", False):
                    return
                # Count toggles that are ON
                total = len(toggles)
                active_count = sum(1 for t in toggles if t.active)

                lock["locked"] = True
                if active_count == total:
                    mc.active = [0]  # all selected → module ON
                else:
                    mc.active = []   # not all selected → module OFF
                lock["locked"] = False
            return callback

        for t in toggles:
            t.on_change("active", make_variable_callback(mc, toggles, lock))

    def apply_clicked():
        try:
            sample = widgets['sample_select'].value
            contig = widgets['contig_select'].value

            # Build requested_features list
            requested_features = []
            for cb_list in widgets['variables_widgets']:
                for cb in cb_list:
                    if cb.active:  # means checkbox is checked
                        requested_features.append(cb.label)

            print(f"[start_bokeh_server] Generating plot for sample={sample}, contig={contig}, features={requested_features}")
            grid = generate_bokeh_plot(conn, requested_features, contig, sample)

            main_placeholder.children = [grid]

        except Exception as e:
            tb = traceback.format_exc()
            main_placeholder.children = [Div(text=f"<pre>Error building plot:\n{tb}</pre>")]

    widgets['apply_button'].on_click(lambda: apply_clicked())

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