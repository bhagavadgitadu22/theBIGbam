import argparse, sqlite3
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from bokeh.models import Range1d, ColumnDataSource, HoverTool, WheelZoomTool
from bokeh.layouts import gridplot
from bokeh.plotting import output_file, save, figure
from dna_features_viewer import BiopythonTranslator

### Custom translator for coloring and labeling features (with DNAFeaturesViewer python library)
# Define function-to-color mapping
# Use the color scheme from pharokka
PHAROKKA_CDS_COLORS = {
    "vfdb_card": "#FF0000",
    "unknown function": "#AAAAAA",
    "other": "#4deeea",
    "tail": "#74ee15",
    "transcription regulation": "#ffe700",
    "dna, rna and nucleotide metabolism": "#f000ff",
    "lysis": "#001eff",
    "moron, auxiliary metabolic gene and host takeover": "#8900ff",
    "integration and excision": "#E0B0FF",
    "head and packaging": "#ff008d",
    "connector": "#5A5A5A",
}

# From https://github.com/oschwengers/bakta/blob/d6443639958750c3bece5822e84978271d1a4dc7/bakta/plot.py#L40
BAKTA_TYPE_COLORS = {
    'CDS': '#cccccc',
    'tRNA': '#b2df8a',
    'tmRNA': '#b2df8a',
    'rRNA': '#fb8072',
    'ncRNA': '#fdb462'
}

class CustomTranslator(BiopythonTranslator):
    def compute_feature_color(self, feature):
        type_feature = feature.type

        if type_feature == "CDS":
            if feature.qualifiers.get("annotation_tool", []) == "pharokka":
                # Get the function field safely
                function = feature.qualifiers.get("function")
                if isinstance(function, list):  # Biopython often stores qualifiers as lists
                    function = function[0] if function else None

                if not isinstance(function, str):  # Missing or wrong type
                    return "#cccccc"

                function = function.lower()

                color_scheme = PHAROKKA_CDS_COLORS
                for key, color in color_scheme.items():
                    if key in function:
                        return color
            
            return "#cccccc"
        
        elif type_feature != "CDS":
            if type_feature not in BAKTA_TYPE_COLORS:
                print("Unknown type of feature:" , type_feature, flush=True)
                return "#000000"
            return BAKTA_TYPE_COLORS.get(type_feature, "#cccccc")

    def compute_feature_label(self, feature):
        return None  # fallback to None if missing or invalid
    
    def compute_feature_html(self, feature):
        return feature.qualifiers.get("product", [])
    
### Plotting functions
def get_contig_info(cur, contig_name):
    cur.execute("SELECT Contig_id, Contig_name, Contig_length, Annotation_tool FROM Contig WHERE Contig_name=?", (contig_name,))
    row = cur.fetchone()
    if row is None:
        raise ValueError(f"Contig not found: {contig_name}")
    return row

def make_bokeh_genemap(conn, contig_id, locus_name, locus_size, annotation_tool, subplot_size, shared_xrange):
    cur = conn.cursor()

    cur.execute("SELECT Start, End, Strand, Type, Product, Function, Phrog FROM Contig_annotation WHERE Contig_id=?", (contig_id,))
    seq_ann_rows = cur.fetchall()

    sequence_annotations = []
    for start, end, strand, ftype, product, function, phrog in seq_ann_rows:
        # Biopython FeatureLocation is 0-based half-open
        try:
            floc = FeatureLocation(start-1, end, strand=strand)
        except Exception:
            continue
        qualifiers = {}
        if product:
            qualifiers['product'] = product
        if function:
            qualifiers['function'] = function
        if phrog:
            qualifiers['phrog'] = phrog
        qualifiers['annotation_tool'] = annotation_tool
        feat = SeqFeature(location=floc, type=ftype, qualifiers=qualifiers)
        sequence_annotations.append(feat)
        
    sequence_records = SeqRecord(Seq('N' * locus_size), id=locus_name, features=sequence_annotations)
    graphic_record = CustomTranslator().translate_record(sequence_records)
    # figure_width and figure_height for the arrow size
    annotation_fig = graphic_record.plot_with_bokeh(figure_width=30, figure_height=40)
    annotation_fig.height = subplot_size

    wheel = WheelZoomTool(dimensions='width')  # only x-axis
    annotation_fig.add_tools(wheel)
    # Don't set active_scroll here - will be set when merging with subplots

    annotation_fig.x_range = shared_xrange

    return annotation_fig

def make_bokeh_subplot(feature_dict, height, x_range, sample_title=None, feature_name=None):
    # Create the figure first (even if empty)
    p = figure(
        height=height,
        x_range=x_range,
        tools="xpan,reset,save"
    )
    
    # Check if we have data to plot
    has_data = feature_dict and any(len(vals["x"]) > 0 for vals in feature_dict)
    
    if not has_data:
        # No data: create empty subplot with message
        if feature_name:
            p.yaxis.axis_label = f"{feature_name} (no data)"
        else:
            p.yaxis.axis_label = "(no data)"
        p.yaxis.axis_label_text_font_size = "10pt"
        p.yaxis.axis_label_standoff = 0
        p.y_range.start = 0
        p.y_range.end = 1
        p.toolbar.logo = None
        p.xgrid.visible = False
        p.ygrid.grid_line_alpha = 0.2
        p.outline_line_color = None
        p.min_border_left = 40
        p.min_border_right = 10
        return p
            
    for data_feature in feature_dict:
        xx = data_feature["x"]
        yy = data_feature["y"]
        
        # Warn if dataset is very large (can cause browser issues)
        if len(xx) > 100000:
            print(f"Warning: Large dataset ({len(xx)} points) may cause browser rendering issues. Consider using a smaller region or higher compression ratio.", flush=True)
        
        type_picked = data_feature["type"]
        color = data_feature["color"]
        alpha = data_feature["alpha"]
        fill_alpha = data_feature["fill_alpha"]
        size = data_feature["size"]
        
        title = data_feature["title"]
        if sample_title:
            title = f"{sample_title} {title}"

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

    wheel = WheelZoomTool(dimensions='width')  # or dimensions='both' for full zoom
    p.add_tools(wheel)
    p.toolbar.active_scroll = wheel

    return p

### Function to get features of one variable
def get_feature_data(cur, feature, contig_id, sample_id, contig_name=None, sample_name=None):
    """Get feature data for plotting.

    Args:
        cur: SQLite cursor
        feature: Feature name to query
        contig_id: Contig ID
        sample_id: Sample ID
        contig_name: Contig name (unused, kept for compatibility)
        sample_name: Sample name (unused, kept for compatibility)
    """
    # Get rendering info from Variable table
    cur.execute(
        "SELECT Type, Color, Alpha, Fill_alpha, Size, Title, Feature_table_name "
        "FROM Variable WHERE Subplot=?",
        (feature,)
    )
    rows = cur.fetchall()

    list_feature_dict = []
    for row in rows:
        type_picked, color, alpha, fill_alpha, size, title, feature_table = row
        feature_dict = {
            "type": type_picked,
            "color": color,
            "alpha": alpha,
            "fill_alpha": fill_alpha,
            "size": size,
            "title": title,
            "x": [],
            "y": []
        }

        # Query Feature_* table (RLE format: First_position, Last_position, Value)
        cur.execute(
            f"SELECT First_position, Last_position, Value FROM {feature_table} "
            "WHERE Sample_id=? AND Contig_id=? ORDER BY First_position",
            (sample_id, contig_id)
        )
        data_rows = cur.fetchall()
        
        # Expand RLE runs into individual points for plotting
        x_coords = []
        y_coords = []
        for first_pos, last_pos, value in data_rows:
            if type_picked == "bars":
                # For bars: expand to all positions in the run
                for pos in range(first_pos, last_pos + 1):
                    x_coords.append(pos)
                    y_coords.append(value)
            else:
                # For curves: only need start and end points
                if first_pos == last_pos:
                    x_coords.append(first_pos)
                    y_coords.append(value)
                else:
                    x_coords.extend([first_pos, last_pos])
                    y_coords.extend([value, value])
        
        feature_dict["x"] = x_coords
        feature_dict["y"] = y_coords

        # Only append if we have actual data points
        if x_coords:
            list_feature_dict.append(feature_dict)

    return list_feature_dict

### Function to generate the bokeh plot
def generate_bokeh_plot_per_sample(conn, list_features, contig_name, sample_name, xstart=None, xend=None, subplot_size=130):
    """Generate a Bokeh plot for a single sample.

    Args:
        conn: SQLite connection
        list_features: List of features/modules to plot (can be mix of modules and individual features)
        contig_name: Name of the contig to plot
        sample_name: Name of the sample to plot
        xstart: Optional x-axis start position
        xend: Optional x-axis end position
        subplot_size: Height of each subplot in pixels
    """
    cur = conn.cursor()

    # Get contig characteristics
    contig_id, locus_name, locus_size, annotation_tool = get_contig_info(cur, contig_name)
    print(f"Locus {locus_name} validated ({locus_size} bp)", flush=True)

    # --- Main gene annotation plot ---
    # Build a SeqRecord from Contig_annotation entries for this contig
    shared_xrange = Range1d(-1000, locus_size+1000)
    if xstart is not None and xend is not None:
        shared_xrange.start = xstart
        shared_xrange.end = xend

    annotation_fig = make_bokeh_genemap(conn, contig_id, locus_name, locus_size, annotation_tool, subplot_size, shared_xrange)

    # Get sample characteristics
    cur.execute("SELECT Sample_id, Sample_name FROM Sample WHERE Sample_name=?", (sample_name,))
    row = cur.fetchone()
    if row is None:
        raise ValueError(f"Sample not found: {sample_name}")
    sample_id, sample_name = row
    print(f"Sample {sample_name} validated.", flush=True)

    # --- Add one subplot per feature requested ---
    # Requested features are variables like 'coverage', 'reads_starts', etc.
    subplots = []
    requested_features = parse_requested_features(list_features)

    for feature in requested_features:
        try:
            list_feature_dict = get_feature_data(cur, feature, contig_id, sample_id,
                                                 contig_name=contig_name, sample_name=sample_name)
            # Always create subplot, even if empty
            subplot_feature = make_bokeh_subplot(list_feature_dict, subplot_size, shared_xrange, feature_name=feature)
            subplots.append(subplot_feature)
        except Exception as e:
            print(f"Error processing feature '{feature}': {e}", flush=True)
            # Create empty subplot with error message
            subplot_feature = make_bokeh_subplot([], subplot_size, shared_xrange, feature_name=f"{feature} (error)")
            subplots.append(subplot_feature)

    # --- Combine all figures in a single grid with one shared toolbar ---
    if not subplots:
        grid = gridplot([[annotation_fig]], merge_tools=True, sizing_mode='stretch_width')
    else:
        all_plots = [annotation_fig] + subplots
        grid = gridplot([[p] for p in all_plots], merge_tools=True, sizing_mode='stretch_width')

    return grid

def save_html_plot_per_sample(db_path, list_features, contig_name, sample_name, subplot_size, output_filename):
    # --- Save interactive HTML plot ---
    output_file(filename = output_filename)
    conn = sqlite3.connect(db_path)
    try:
        grid = generate_bokeh_plot_per_sample(conn, list_features, contig_name, sample_name, subplot_size=subplot_size)
        save(grid)
    finally:
        conn.close()

### Parsing features
def parse_requested_features(list_features):
    """Parse requested features, expanding modules to individual features.
    
    Accepts a mix of module names and individual feature names.
    Module names are case-insensitive and can include:
    - "coverage" or "Coverage" -> primary_reads, secondary_reads, supplementary_reads
    - "phagetermini" or "Phage termini" -> coverage_reduced, reads_starts, reads_ends, tau
    - "assemblycheck" or "Assembly check" -> all assembly check features
    
    Returns deduplicated list of individual feature names.
    """
    features = []
    for item in list_features:
        item_lower = item.lower().strip()
        
        # Module: Coverage
        if item_lower in ["coverage"]:
            features.extend(["Primary alignments", "Other alignments", "Other alignments"])
        # Module: Phage termini / phagetermini
        elif item_lower in ["phage termini", "phagetermini", "phage_termini"]:
            features.extend(["Coverage reduced", "Reads termini", "Tau"])
        # Module: Assembly check / assemblycheck
        elif item_lower in ["assembly check", "assemblycheck", "assembly_check"]:
            features.extend(["Clippings", "Indels", "Mismatches", "Read lengths", "Insert sizes", "Bad orientations"])
        # Individual feature
        else:
            features.append(item)
    
    # Deduplicate while preserving order
    seen = set()
    deduped_features = [f for f in features if not (f in seen or seen.add(f))]
    return deduped_features

### Main function
def add_plot_per_sample_args(parser):
    parser.add_argument("--db", required=True, help="Path to sqlite database file to store results")
    parser.add_argument("--variables", required=True, help="List of variables or full modules to compute (comma-separated) (options allowed for modules: coverage, phagetermini, assemblycheck)")
    parser.add_argument("--contig", required=True, help="Name of the contig to plot")
    parser.add_argument("--sample", required=True, help="Name of the sample to plot")
    parser.add_argument("--html", required=False, default="MGFeatureViewer_per_sample.html", help="Name for output html files. A bokeh server will be started if not provided")
    parser.add_argument("--subplot_height", required=False, default=130, help="Height of each subplot (in pixels)")

def run_plot_per_sample(args):
    db_path = args.db
    list_features = args.variables.split(",")
    contig_name = args.contig
    sample_name = args.sample
    output_filename = args.html
    subplot_size = int(args.subplot_height)

    print(f"Saving static HTML to {output_filename}...", flush=True)
    save_html_plot_per_sample(db_path, list_features, contig_name, sample_name, subplot_size, output_filename)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    add_plot_per_sample_args(parser)
    args = parser.parse_args()
    run_plot_per_sample(args)
