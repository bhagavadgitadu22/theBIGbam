from bokeh.models import Range1d,ColumnDataSource, HoverTool
from bokeh.layouts import column
from bokeh.plotting import output_file, save, figure
from dna_features_viewer import BiopythonTranslator
import constants

### Custom translator for coloring and labeling features (with DNAFeaturesViewer python library)
class CustomTranslator(BiopythonTranslator):
    def compute_feature_color(self, feature):
        if ANNOTATION_TOOL == "pharokka":
            function = feature.qualifiers.get("function", [""])[0].lower()
            color_scheme = constants.PHAROKKA_COLORS
            for key, color in color_scheme.items():
                if key in function:
                    return color
        return "#AAAAAA"

    def compute_feature_label(self, feature):
        return None  # fallback to None if missing or invalid
    
    def compute_feature_html(self, feature):
        return feature.qualifiers.get("product", [])
    
### Plotting functions
def make_bokeh_subplot(feature, xx, yy, width, height, x_range):
    type_picked = constants.FEATURE_SUBPLOTS[feature]["type_picked"]
    color_picked = constants.FEATURE_SUBPLOTS[feature]["color_picked"]
    alpha_picked = constants.FEATURE_SUBPLOTS[feature]["alpha_picked"]
    size_picked = constants.FEATURE_SUBPLOTS[feature]["size_picked"]
    title_picked = constants.FEATURE_SUBPLOTS[feature]["title_picked"]

    p = figure(
        width=width,
        height=height,
        title=title_picked,
        x_range=x_range,
        tools="xpan,xwheel_zoom,reset,save"
    )
    source = ColumnDataSource(data=dict(x=xx, y=yy))

    # Part specific to the type of subplot
    if type_picked == "curve":
        p.line(
            x='x',
            y='y',
            source=source,
            line_color=color_picked,
            line_alpha=alpha_picked,
            line_width=size_picked,
        )
    elif type_picked == "bars":
        source = ColumnDataSource(data=dict(x=xx, y=yy))
        p.vbar(
            x='x',
            bottom=0,
            top='y',
            source=source,
            color=color_picked,
            alpha=alpha_picked,
            width=size_picked
        )

    # Add hover
    hover = HoverTool(tooltips=[("Position", "@x"), ("Number", "@y")], mode='vline')
    p.add_tools(hover)

    # A clean style like your matplotlib setup
    p.toolbar.logo = None
    p.xgrid.visible = False

    p.y_range.start = 0
    p.yaxis.axis_label = title_picked
    p.yaxis.axis_label_text_font_size = "10pt"
    p.yaxis.axis_label_standoff = 0
    p.ygrid.grid_line_alpha = 0.2
    p.yaxis.axis_label = None
    
    p.outline_line_color = None  # hides top/right borders
    p.min_border_left = 40
    p.min_border_right = 10

    return p

def prepare_subplot(feature, features_x, features_y, shared_xrange, max_visible_width, subplot_size):
    # Skip if nothing left
    if len(features_x) == 0:
        return None
    feature_subplot = make_bokeh_subplot(feature, features_x, features_y, max_visible_width, subplot_size, shared_xrange)
    return feature_subplot

### One function to rule them all
def prepare_all_subplots(data, max_visible_width, subplot_size, shared_xrange, window_size):
    subplots = []
    for feature in data:
        feature_positions = data[feature]["x"]
        feature_values = data[feature]["y"]

        subplot_feature = prepare_subplot(feature, feature_positions, feature_values, shared_xrange, max_visible_width, subplot_size)
        if subplot_feature is not None:
            subplots.append(subplot_feature)

    return subplots

def prepare_main_plot(data_dictionary, genbank_record, protein_annotation_tool, locus_size, window_size, max_visible_width, subplot_size, output_name):
    global ANNOTATION_TOOL
    ANNOTATION_TOOL = protein_annotation_tool

    # Plotting gene map
    print("Plotting gene map...", flush=True)
    graphic_record = CustomTranslator().translate_record(genbank_record)
    # figure_width and figure_height for the arrow size
    annotation_fig = graphic_record.plot_with_bokeh(figure_width=30, figure_height=40)
    annotation_fig.width = max_visible_width
    annotation_fig.height = subplot_size

    shared_xrange = Range1d(0, locus_size)
    annotation_fig.x_range = shared_xrange

    subplots = prepare_all_subplots(data_dictionary, max_visible_width, subplot_size, shared_xrange, window_size)

    layout = column(annotation_fig, *subplots)
    output_file(output_name)
    save(layout)
    print(f"Saved interactive plot to {output_name}")
