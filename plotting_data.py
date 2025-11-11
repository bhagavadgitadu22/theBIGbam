import argparse, sys, os, csv
from Bio import SeqIO
from bokeh.models import Range1d,ColumnDataSource, HoverTool
from bokeh.layouts import column, gridplot
from bokeh.plotting import output_file, save, figure
from dna_features_viewer import BiopythonTranslator
import colors_for_genbank

### Custom translator for coloring and labeling features (with DNAFeaturesViewer python library)
class CustomTranslator(BiopythonTranslator):
    def compute_feature_color(self, feature):
        if ANNOTATION_TOOL == "pharokka":
            function = feature.qualifiers.get("function", [""])[0].lower()
            color_scheme = colors_for_genbank.PHAROKKA_COLORS
            for key, color in color_scheme.items():
                if key in function:
                    return color
        return "#AAAAAA"

    def compute_feature_label(self, feature):
        return None  # fallback to None if missing or invalid
    
    def compute_feature_html(self, feature):
        return feature.qualifiers.get("product", [])
    
### Plotting functions
def make_bokeh_subplot(xx, yy, width, height, x_range, type_picked, color, alpha, size, title):
    p = figure(
        width=width,
        height=height,
        title=title,
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
            line_color=color,
            line_alpha=alpha,
            line_width=size,
        )
    elif type_picked == "bars":
        source = ColumnDataSource(data=dict(x=xx, y=yy))
        p.vbar(
            x='x',
            bottom=0,
            top='y',
            source=source,
            color=color,
            alpha=alpha,
            width=size
        )

    # Add hover
    hover = HoverTool(tooltips=[("Position", "@x"), ("Number", "@y")], mode='vline')
    p.add_tools(hover)

    # A clean style like your matplotlib setup
    p.toolbar.logo = None
    p.xgrid.visible = False

    p.y_range.start = 0
    p.yaxis.axis_label = title
    p.yaxis.axis_label_text_font_size = "10pt"
    p.yaxis.axis_label_standoff = 0
    p.ygrid.grid_line_alpha = 0.2
    p.yaxis.axis_label = None
    
    p.outline_line_color = None  # hides top/right borders
    p.min_border_left = 40
    p.min_border_right = 10

    return p

### One function to rule them all
def prepare_all_subplots(data_all_features, max_visible_width, subplot_size, shared_xrange):
    subplots = []
    for feature_dict in data_all_features:
        feature_positions = feature_dict["x"]
        feature_values = feature_dict["y"]

        # Skip if nothing left
        if len(feature_positions) == 0:
            continue
        
        # Get subplot characteristics
        type_picked = feature_dict["type"]
        color = feature_dict["color"]
        alpha = feature_dict["alpha"]
        size = feature_dict["size"]
        title = feature_dict["title"]

        # Create subplot
        subplot_feature = make_bokeh_subplot(feature_positions, feature_values, 
                                             max_visible_width, subplot_size, shared_xrange,
                                             type_picked, color, alpha, size, title)
        if subplot_feature is not None:
            subplots.append(subplot_feature)

    return subplots

def prepare_main_plot(data_dictionary, genbank_record, locus_size, max_visible_width, subplot_size, output_name):
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
def generate_html_plots(genbank_path, bam_files, input_csv_dir, requested_features, 
                        custom_characs, output_prefix, max_visible_width, subplot_size):
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
            all_features_for_one_locus_one_sample = []
            for feature in requested_features:
                feature_dict = {}
                feature_dict["filename"] = os.path.join(input_csv_dir, f"{feature}_values_for_{locus_name}_in_{sample_name}.csv")
                
                # Subplot characteristics
                feature_dict["type"] = colors_for_genbank.FEATURE_SUBPLOTS[feature]["type"]
                feature_dict["color"] = colors_for_genbank.FEATURE_SUBPLOTS[feature]["color"]
                feature_dict["alpha"] = colors_for_genbank.FEATURE_SUBPLOTS[feature]["alpha"]
                feature_dict["size"] = colors_for_genbank.FEATURE_SUBPLOTS[feature]["size"]
                feature_dict["title"] = colors_for_genbank.FEATURE_SUBPLOTS[feature]["title"]
                
                # Initialize subplot values
                feature_dict["x"] = []
                feature_dict["y"] = []
                all_features_for_one_locus_one_sample.append(feature_dict)

            # Add custom characs to the list
            all_features_for_one_locus_one_sample.extend(custom_characs)

            # For all subplots
            for feature_dict in all_features_for_one_locus_one_sample:
                if not os.path.exists(feature_dict["filename"]):
                    print(f"WARNING: CSV file {feature_dict['filename']} not found")
                else:
                    with open(feature_dict["filename"], 'r') as csvfile:
                        reader = csv.reader(csvfile)
                        for row in reader:
                            # For custom files need to keep only relevant rows
                            if row[0] == sample_name and row[1] == locus_name:
                                feature_dict["x"].append(int(row[2]))
                                feature_dict["y"].append(float(row[3]))

            # Generate the HTML plot
            output_html = f"{output_prefix}_{locus_name}_in_{sample_name}.html"
            prepare_main_plot(
                all_features_for_one_locus_one_sample,
                record,
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
    parser.add_argument("-c", "--custom", required=False, help="List of custom variables to plot. Each variable should be in the format filename:type:color:title. Type can be 'curve' or 'bars'. Use comma between variables")
    parser.add_argument("-o", "--output_prefix", required=False, default="MGFeaturesViewer", help="Prefix for output files, including complete path if you want to save them in a specific folder")
    parser.add_argument("-pw", "--plot_width", required=False, default=1800, help="Width of the plot (in pixels)")
    parser.add_argument("-sh", "--subplot_height", required=False, default=130, help="Height of each subplot (in pixels)")
    args = parser.parse_args()

    genbank_file = args.genbank
    global ANNOTATION_TOOL
    ANNOTATION_TOOL = args.annotation
    
    bam_files = args.bam_files
    input_csv_dir = args.input_dir

    # Get requested features
    requested_modules = args.modules.split(",")
    requested_features = parse_requested_features(requested_modules)

    custom_variables = args.custom
    custom_characs = []
    if custom_variables:
        for custom_var in custom_variables.split(","):
            filename, type_picked, color, title = custom_var.split(":")
            custom_charac = {"filename": filename, "x": [], "y": []}
            custom_charac.update(colors_for_genbank.config_feature_subplot(type_picked, color, title))
            custom_characs.append(custom_charac)
    
    # Optional plotting parameters
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
    generate_html_plots(genbank_file, bam_names, input_csv_dir, requested_features, 
                        custom_characs, output_prefix, max_visible_width, subplot_size)    
    
if __name__ == "__main__":
    main()