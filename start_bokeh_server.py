import argparse
import sqlite3
import traceback

from bokeh.layouts import column, row
from bokeh.models import Div
from bokeh.models.widgets import Select, CheckboxGroup, CheckboxButtonGroup, Button

# Import the plotting function from the repo
from plotting_data import generate_bokeh_plot

def build_controls(conn):
    """Query DB and return widgets and helper mappings."""
    cur = conn.cursor()

    # Widget Selector for Samples
    cur.execute("SELECT Sample_name FROM Sample ORDER BY Sample_name")
    samples = [r[0] for r in cur.fetchall()]
    sample_select = Select(title="Sample", value=samples[0], options=samples)

    # Widget Selector for Contigs
    cur.execute("SELECT Contig_name FROM Contig ORDER BY Contig_name")
    contigs = [r[0] for r in cur.fetchall()]
    contig_select = Select(title="Contig", value=contigs[0], options=contigs)

    # Modules and variables
    cur.execute("SELECT DISTINCT Module FROM Variable")
    modules = [r[0] for r in cur.fetchall()]

    # For each module get variables
    module_widgets = []
    variables_widgets = []
    for module in modules:
        cur.execute("SELECT Variable_name FROM Variable WHERE Module=?", (module,))
        variables_checkbox = [r[0] for r in cur.fetchall()]

        if len(variables_checkbox) > 1:
            module_checkbox = CheckboxGroup(labels=[module], active=[])
            module_widgets.append(module_checkbox)
        else:
            module_widgets.append(None)

        # use CheckboxButtonGroup for selecting individual variables
        vars_buttons = CheckboxButtonGroup(labels=variables_checkbox, active=[])
        variables_widgets.append(vars_buttons)

    apply_button = Button(label="Apply", button_type="primary")

    widgets = {
        'sample_select': sample_select,
        'contig_select': contig_select,
        'module_widgets': module_widgets,
        'variables_widgets': variables_widgets,
        'apply_button': apply_button
    }
    return widgets

def modify_doc_factory(db_path, max_visible_width=1800, subplot_size=130):
    """Return a modify_doc(doc) function to be used by Bokeh server application."""
    instructions = Div(text="<b>Select elements to plot then click Apply:</b>")

    conn = sqlite3.connect(db_path)
    widgets = build_controls(conn)

    variables_title = Div(text="Variables")
    controls_children = [instructions, widgets['sample_select'], widgets['contig_select'], variables_title]
    
    # Append variable selectors
    for i, module_widget in enumerate(widgets['module_widgets']):
        if module_widget is not None:
            controls_children.append(module_widget)
        
        var_widget = widgets['variables_widgets'][i]
        controls_children.append(var_widget)

    controls_children.append(widgets['apply_button'])
    controls_column = column(*controls_children, width=350)

    main_placeholder = column(Div(text="<i>No plot yet. Select options and click Apply.</i>"))
    layout = row(controls_column, main_placeholder)

    ### Attach callbacks
    def make_module_callback(i, mc, vb):
        """Module checkbox clicked → toggle its variable buttons."""
        def callback(attr, old, new):
            if 0 in mc.active:   # module checkbox active
                vb.active = list(range(len(vb.labels)))
            else:
                vb.active = []
        return callback

    def make_variable_callback(i, mc, vb):
        """Variable toggled → update module checkbox."""
        def callback(attr, old, new):
            total = len(vb.labels)

            if len(vb.active) == total:
                # all active → module checkbox should be active
                if mc.active != [0]:
                    mc.active = [0]
            else:
                # some inactive → module checkbox should be inactive
                if mc.active != []:
                    mc.active = []
        return callback
    
    # module checkbox callback
    for i, mc in enumerate(widgets['module_widgets']):
        vb = widgets['variables_widgets'][i]
        if mc is None:
            continue

        # module → variables
        mc.on_change("active", make_module_callback(i, mc, vb))
        # variables → module
        vb.on_change("active", make_variable_callback(i, mc, vb))

    def apply_clicked():
        try:
            sample = widgets['sample_select'].value
            contig = widgets['contig_select'].value

            # Build requested_features list
            requested_features = []
            for variables_widget in widgets['variables_widgets']:
                for vidx in variables_widget.active:
                    requested_features.append(variables_widget.labels[vidx])

            # Call plotting_data.generate_bokeh_plot to build the grid
            print(f"[start_bokeh_server] Generating plot for sample={sample}, contig={contig}, features={requested_features}")
            grid = generate_bokeh_plot(db_path, requested_features, contig, sample, max_visible_width, subplot_size)

            # Replace main panel
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
