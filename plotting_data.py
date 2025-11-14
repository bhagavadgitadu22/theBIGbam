import argparse, sys, os, csv, sqlite3
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from bokeh.models import Range1d,ColumnDataSource, HoverTool
from bokeh.layouts import column, gridplot
from bokeh.plotting import output_file, save, figure
from dna_features_viewer import BiopythonTranslator
import colors_for_genbank

### Custom translator for coloring and labeling features (with DNAFeaturesViewer python library)
class CustomTranslator(BiopythonTranslator):
    def compute_feature_color(self, feature):
        if ANNOTATION_TOOL == "pharokka":
            function = feature.qualifiers.get("function", [""]).lower()
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
def make_bokeh_subplot(feature_dict, width, height, x_range):
    p = figure(
        width=width,
        height=height,
        x_range=x_range,
        tools="xpan,xwheel_zoom,reset,save"
    )
            
    for feature in feature_dict:
        data_feature = feature_dict[feature]
        xx = data_feature["x"]
        yy = data_feature["y"]
        type_picked = data_feature["type"]
        color = data_feature["color"]
        alpha = data_feature["alpha"]
        fill_alpha = data_feature["fill_alpha"]
        size = data_feature["size"]
        title = data_feature["title"]

        source = ColumnDataSource(data=dict(x=xx, y=yy))

        # Part specific to the type of subplot
        if type_picked == "curve":
            p.varea(
                x='x',
                y1=0,
                y2='y',
                source=source,
                fill_color=color,
                fill_alpha=fill_alpha,
                legend_label = title
            )
            p.line(
                x='x',
                y='y',
                source=source,
                line_color=color,
                line_alpha=alpha,
                line_width=size,
                legend_label = title
            )
        elif type_picked == "bars":
            p.vbar(
                x='x',
                bottom=0,
                top='y',
                source=source,
                color=color,
                alpha=alpha,
                width=size,
                legend_label = title
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

    p.legend.location = "top_left"
    if len(feature_dict) > 1:
        p.legend.click_policy="hide"

    return p

### Function to get features of one variable
def get_feature_data(cur, feature, contig_id, sample_id):
    feature_dict = {} 

    # Query Variable table to get rendering info and feature table name
    cur.execute("SELECT Type, Color, Alpha, Fill_alpha, Size, Title, Feature_table_name FROM Variable WHERE Variable_name=?", (feature,))
    row = cur.fetchone()

    type_picked, color, alpha, fill_alpha, size, title, feature_table = row
    feature_dict["type"] = type_picked
    feature_dict["color"] = color
    feature_dict["alpha"] = alpha
    feature_dict["fill_alpha"] = fill_alpha
    feature_dict["size"] = size
    feature_dict["title"] = title
    feature_dict["x"] = []
    feature_dict["y"] = []

    # Query feature table for this sample and contig
    try:
        cur.execute(f"SELECT Position, Value FROM {feature_table} WHERE Sample_id=? AND Contig_id=? ORDER BY Position", (sample_id, contig_id))
        rows = cur.fetchall()
        feature_dict["x"] = [r[0] for r in rows]
        feature_dict["y"] = [r[1] for r in rows]
    except Exception as e:
        print(f"WARNING: Could not read feature table {feature_table}: {e}", flush=True)

    return feature_dict

### Function to generate one HTML plot per locus
def generate_bokeh_plot(db_path, requested_features, contig_name, sample_name, max_visible_width, subplot_size):
    # Connect to DB and gather contigs and samples
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    # Get contig characteristics
    cur.execute("SELECT Contig_id, Contig_name, Contig_length FROM Contig WHERE Contig_name=?", (contig_name,))
    row = cur.fetchone()
    if row is None:
        raise ValueError(f"Contig not found: {contig_name}")
    contig_id, locus_name, locus_size = row
    print(f"Locus {locus_name} validated ({locus_size} bp)", flush=True)

    # Get sample characteristics
    cur.execute("SELECT Sample_id, Sample_name FROM Sample WHERE Sample_name=?", (sample_name,))
    row = cur.fetchone()
    if row is None:
        raise ValueError(f"Sample not found: {sample_name}")
    sample_id, sample_name = row
    print(f"Sample {sample_name} validated.", flush=True)

    # --- Main gene annotation plot ---
    # Build a SeqRecord from Sequence_annotation entries for this contig
    cur.execute("SELECT Start, End, Strand, Type, Product, Function, Phrog FROM Sequence_annotation WHERE Contig_id=?", (contig_id,))
    seq_ann_rows = cur.fetchall()

    sequence_annotations = []
    for start, end, strand, ftype, product, function, phrog in seq_ann_rows:
        # Biopython FeatureLocation is 0-based half-open
        try:
            floc = FeatureLocation(start - 1, end, strand=strand)
        except Exception:
            continue
        qualifiers = {}
        if product:
            qualifiers['product'] = product
        if function:
            qualifiers['function'] = function
        if phrog:
            qualifiers['phrog'] = phrog
        feat = SeqFeature(location=floc, type=ftype, qualifiers=qualifiers)
        sequence_annotations.append(feat)
        
    sequence_records = SeqRecord(Seq('N' * locus_size), id=locus_name, features=sequence_annotations)
    graphic_record = CustomTranslator().translate_record(sequence_records)
    # figure_width and figure_height for the arrow size
    annotation_fig = graphic_record.plot_with_bokeh(figure_width=30, figure_height=40)
    annotation_fig.width = max_visible_width
    annotation_fig.height = subplot_size

    shared_xrange = Range1d(0, locus_size)
    annotation_fig.x_range = shared_xrange

    # --- Add one subplot per feature requested ---
    # Requested features are variables like 'coverage', 'reads_starts', etc.
    subplots = []
    for feature in requested_features:
        feature_dict = {}
        
        if feature in ["clippings", "indels", "reads_ends"]:
            if feature == "clippings":
                feature_dict["right_clippings"] = get_feature_data(cur, "right_clippings", contig_id, sample_id)
                feature_dict["left_clippings"] = get_feature_data(cur, "left_clippings", contig_id, sample_id)
            elif feature == "indels":
                feature_dict["insertions"] = get_feature_data(cur, "insertions", contig_id, sample_id)
                feature_dict["deletions"] = get_feature_data(cur, "deletions", contig_id, sample_id)
            elif feature == "reads_ends":
                feature_dict["reads_starts"] = get_feature_data(cur, "reads_starts", contig_id, sample_id)
                feature_dict["reads_ends"] = get_feature_data(cur, "reads_ends", contig_id, sample_id)
            
            if all(len(vals["x"]) == 0 for vals in feature_dict.values()):
                continue
            subplot_feature = make_bokeh_subplot(feature_dict, max_visible_width, subplot_size, shared_xrange)
            if subplot_feature is not None:
                subplots.append(subplot_feature)
 
        else:            
            feature_dict[feature] = get_feature_data(cur, feature, contig_id, sample_id)
            
            if len(feature_dict[feature]["x"]) == 0:
                continue
            subplot_feature = make_bokeh_subplot(feature_dict, max_visible_width, subplot_size, shared_xrange)
            if subplot_feature is not None:
                subplots.append(subplot_feature)

    conn.close()

    # --- Combine all figures in a single grid with one shared toolbar ---
    all_plots = [annotation_fig] + subplots
    grid = gridplot([[p] for p in all_plots], merge_tools=True)

    return grid

def save_html_plot(db_path, requested_features, contig_name, sample_name, output_prefix, max_visible_width, subplot_size): 
    # --- Save interactive HTML plot ---
    # If output_prefix looks like a directory path (contains a separator or is an existing dir),
    # create that directory and write the HTML files inside it. Otherwise use it as a filename prefix.
    norm_prefix = os.path.normpath(output_prefix)
    looks_like_dir = os.path.isdir(norm_prefix) or (os.path.sep in output_prefix) or ("/" in output_prefix) or ("\\" in output_prefix)
    if looks_like_dir:
        outdir = norm_prefix
        os.makedirs(outdir, exist_ok=True)
        # use the directory basename as a short prefix for files
        base = os.path.basename(outdir) or "output"
        output_html = os.path.join(outdir, f"{base}_{contig_name}_in_{sample_name}.html")
    else:
        output_html = f"{output_prefix}_{contig_name}_in_{sample_name}.html"

    output_file(output_html)
    grid = generate_bokeh_plot(db_path, requested_features, contig_name, sample_name, max_visible_width, subplot_size)
    save(grid)

### Parsing features
def parse_requested_features(list_features):
    features = []
    for feature in list_features:
        if feature == "coverage":
            features.append("coverage")
        elif feature == "assemblycheck":
            features.extend(["read_lengths", "insert_sizes", "bad_orientations", "clippings", "indels", "mismatches"])
        elif feature == "phagetermini":
            features.extend(["coverage_reduced", "reads_ends", "tau"])
        else:
            features.append(feature)
    return(features)

### Main function
def main():
    # Parse command line arguments
    print("Parsing arguments...", flush=True)
    parser = argparse.ArgumentParser(description="Parse input files.")
    parser.add_argument("-t", "--threads", required=True, help="Number of threads available")
    parser.add_argument("-d", "--db", required=True, help="Path to sqlite database file to store results")
    parser.add_argument("-m", "--modules", required=True, help="List of modules to compute (comma-separated) (options allowed: coverage, phagetermini, assemblycheck)")
    parser.add_argument("--contig", required=True, help="Name of the contig to plot")
    parser.add_argument("--sample", required=True, help="Name of the sample to plot")
    parser.add_argument("-o", "--output_prefix", required=False, default="MGFeaturesViewer", help="Prefix for output files, including complete path if you want to save them in a specific folder")
    parser.add_argument("--color", required=False, help="Color system for the sequence annotations (options allowed 'pharokka' or 'other')")
    parser.add_argument("--plot_width", required=False, default=1800, help="Width of the plot (in pixels)")
    parser.add_argument("--subplot_height", required=False, default=130, help="Height of each subplot (in pixels)")
    args = parser.parse_args()

    # Get requested features
    list_features = args.modules.split(",")
    requested_features = parse_requested_features(list_features)

    # Path parameters
    db_path = args.db
    contig_name = args.contig
    sample_name = args.sample
    output_prefix = args.output_prefix

    # Optional plotting parameters
    global ANNOTATION_TOOL
    ANNOTATION_TOOL = args.color if args.color else "other"
    max_visible_width = int(args.plot_width)
    subplot_size = int(args.subplot_height)

    # Reading values from database and plotting (one output file per locus and sample)
    print("Reading database and plotting (one output file per locus and sample)...", flush=True)
    save_html_plot(db_path, requested_features, contig_name, sample_name, output_prefix, max_visible_width, subplot_size)
    
if __name__ == "__main__":
    main()