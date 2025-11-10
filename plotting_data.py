import argparse, sys, os, csv
from Bio import SeqIO
from bokeh.models import Range1d,ColumnDataSource, HoverTool
from bokeh.layouts import column, gridplot
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
def prepare_all_subplots(data, max_visible_width, subplot_size, shared_xrange):
    subplots = []
    for feature in data:
        feature_positions = data[feature]["x"]
        feature_values = data[feature]["y"]

        subplot_feature = prepare_subplot(feature, feature_positions, feature_values, shared_xrange, max_visible_width, subplot_size)
        if subplot_feature is not None:
            subplots.append(subplot_feature)

    return subplots

def prepare_main_plot(data_dictionary, genbank_record, protein_annotation_tool, locus_size, max_visible_width, subplot_size, output_name):
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

    subplots = prepare_all_subplots(data_dictionary, max_visible_width, subplot_size, shared_xrange)

    layout = column(annotation_fig, *subplots)
    output_file(output_name)
    save(layout)
    print(f"Saved interactive plot to {output_name}")

def prepare_main_plot(data_dictionary, genbank_record, protein_annotation_tool, locus_size, max_visible_width, subplot_size, output_name):
    global ANNOTATION_TOOL
    ANNOTATION_TOOL = protein_annotation_tool

    # --- Main gene annotation plot
    print("Plotting gene map...", flush=True)
    graphic_record = CustomTranslator().translate_record(genbank_record)
    # figure_width and figure_height for the arrow size
    annotation_fig = graphic_record.plot_with_bokeh(figure_width=30, figure_height=40)
    annotation_fig.width = max_visible_width
    annotation_fig.height = subplot_size

    shared_xrange = Range1d(0, locus_size)
    annotation_fig.x_range = shared_xrange

    # --- Prepare subplots
    subplots = prepare_all_subplots(data_dictionary, max_visible_width,subplot_size, shared_xrange)

    # --- Combine all figures in a single grid with one shared toolbar
    all_plots = [annotation_fig] + subplots
    grid = gridplot([[p] for p in all_plots], merge_tools=True)

    # --- Save interactive HTML
    output_file(output_name)
    save(grid)
    print(f"Saved interactive plot with shared toolbar to {output_name}")

### Parsing features
def parse_requested_features(requested_modules):
    features = []
    if "coverage" in requested_modules:
        features = ["coverage"]
    if "assemblycheck" in requested_modules:
        features.extend(["read_lengths", "insert_sizes", "bad_orientations", "left_clippings", "right_clippings", "insertions", "deletions", "mismatches"])
    if "phagetermini" in requested_modules:
        features.extend(["coverage_reduced", "reads_starts", "reads_ends", "tau"])
    return(features)

### Function to generate one HTML plot per locus
def generate_html_plots(genbank_path, annotation_tool, bam_files, input_csv_dir, requested_features, output_prefix, max_visible_width, subplot_size):
    print("Generating HTML plots per locus...", flush=True)

    allowed_types = ["CDS", "tRNA/tmRNA", "rRNA", "ncRNA", "ncRNA-region", "CRISPR", "Gap", "Misc"]
    # Parse GenBank file
    for record in SeqIO.parse(genbank_path, "genbank"):
        locus_name = record.id
        locus_size = len(record.seq)

        print(f"Processing locus: {locus_name} ({locus_size} bp)", flush=True)

        # Filter allowed feature types
        filtered_features_gbk = [f for f in record.features if f.type in allowed_types]
        record.features = filtered_features_gbk

        for sample_name in bam_files:
            # Reading data from CSV files for all requested features
            data_for_locus_sample = {}
            for feature in requested_features:
                data_for_locus_sample.setdefault(feature, {"x": [], "y": []})

                csv_file_path = os.path.join(input_csv_dir, f"{feature}_values_for_{locus_name}_in_{sample_name}.csv")
                if not os.path.exists(csv_file_path):
                    print(f"WARNING: CSV file {csv_file_path} not found")
                else:
                    with open(csv_file_path, 'r') as csvfile:
                                        reader = csv.reader(csvfile)
                                        next(reader)
                                        for row in reader:
                                            position = int(row[2])
                                            value = float(row[3])
                                            data_for_locus_sample[feature]["x"].append(position)
                                            data_for_locus_sample[feature]["y"].append(value)

            # Generate the HTML plot
            output_html = f"{output_prefix}_{locus_name}_in_{sample_name}.html"
            prepare_main_plot(
                data_for_locus_sample,
                record,
                annotation_tool,
                locus_size,
                max_visible_width,
                subplot_size,
                output_html,
            )
            print(f"Saved: {output_html}", flush=True)

### Main function
def main():
    # Parse command line arguments
    print("Parsing arguments...", flush=True)
    parser = argparse.ArgumentParser(description="Parse input files.")
    parser.add_argument("-t", "--threads", required=True, help="Number of threads available")
    parser.add_argument("-g", "--genbank", required=True, help="Path to genbank file of all investigated contigs")
    parser.add_argument("-a", "--annotation", required=True, help="Annotation tool used (options allowed 'pharokka' or 'other')")
    parser.add_argument("-b", "--bam_files", required=True, help="Path to bam file or directory containing mapping files (BAM format)")
    parser.add_argument("-i", "--input_dir", required=True, help="Input directory where input csv files are stored")
    parser.add_argument("-m", "--modules", required=True, help="List of modules to compute (comma-separated) (options allowed: coverage, phagetermini, assemblycheck)")
    parser.add_argument("-o", "--output_prefix", required=False, default="MGFeaturesViewer", help="Prefix for output files, including complete path if you want to save them in a specific folder")
    parser.add_argument("-pw", "--plot_width", required=False, default=1800, help="Width of the plot (in pixels)")
    parser.add_argument("-sh", "--subplot_height", required=False, default=130, help="Height of each subplot (in pixels)")
    args = parser.parse_args()

    genbank_file = args.genbank
    annotation_method = args.annotation
    bam_files = args.bam_files
    input_csv_dir = args.input_dir
    requested_modules = args.modules.split(",")
    requested_features = parse_requested_features(requested_modules)
    output_prefix = args.output_prefix
    max_visible_width = int(args.plot_width)
    subplot_size = int(args.subplot_height)

    # Get list of BAM files
    if os.path.isdir(args.bam_files):
        bam_names = [os.path.basename(bam_file).replace(".bam", "") for bam_file in os.listdir(bam_files) if bam_file.endswith(".bam")]
    else:
        bam_names = [os.path.basename(bam_files).replace(".bam", "")]

    # Reading values for all requested modules from mapping files
    print("Reading csv files and plotting (one output file per locus and sample)...", flush=True)
    generate_html_plots(genbank_file, annotation_method, bam_names, input_csv_dir, requested_features, output_prefix, max_visible_width, subplot_size)    
    
if __name__ == "__main__":
    main()