import argparse
import duckdb
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from bokeh.models import Range1d, ColumnDataSource, HoverTool, WheelZoomTool, NumeralTickFormatter, Label, TapTool
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

def make_bokeh_genemap(conn, contig_id, locus_name, locus_size, annotation_tool, subplot_size, shared_xrange, xstart=None, xend=None):
    cur = conn.cursor()

    # Build position filter clause for annotations
    position_filter = ""
    params = [contig_id]
    if xstart is not None and xend is not None:
        position_filter = " AND \"End\" >= ? AND \"Start\" <= ?"
        params.extend([xstart, xend])

    cur.execute(
        f"SELECT \"Start\", \"End\", Strand, \"Type\", Product, \"Function\", Phrog FROM Contig_annotation WHERE Contig_id=?{position_filter}",
        tuple(params)
    )
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

    # Remove tap tool added by DNAFeaturesViewer
    annotation_fig.tools = [t for t in annotation_fig.tools if not isinstance(t, TapTool)]

    # Disable scientific notation on x-axis
    annotation_fig.xaxis.formatter = NumeralTickFormatter(format="0,0")
    
    wheel = WheelZoomTool(dimensions='width')  # only x-axis
    annotation_fig.add_tools(wheel)
    annotation_fig.toolbar.active_scroll = wheel

    annotation_fig.x_range = shared_xrange

    return annotation_fig

def make_bokeh_subplot(feature_dict, height, x_range, sample_title=None, feature_name=None):
    # Create the figure first (even if empty)
    p = figure(
        height=height,
        x_range=x_range,
        tools="xpan,reset,save"
    )
    
    # Disable scientific notation on x-axis
    p.xaxis.formatter = NumeralTickFormatter(format="0,0")
    
    # Check if we have data to plot
    # You need one dataset of the subplots to have at least one non-zero points
    has_data = bool(feature_dict) and any(any(y > 0 for y in d["y"]) for d in feature_dict)
    title = ""
    if not has_data:
        return None
    else:      
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

            # Prepare data for ColumnDataSource
            data_dict = dict(x=xx, y=yy)
            
            # Add width for bars if available and any width is different from 1
            if type_picked == "bars" and "width" in data_feature:
                data_dict["first_pos"] = data_feature["first_pos"]
                data_dict["last_pos"] = data_feature["last_pos"]
                # Only use variable width for rendering if this specific feature has width != 1
                if any(w != 1 for w in data_feature["width"]):
                    data_dict["width"] = data_feature["width"]

            # Add duplication-specific fields if available
            if data_feature.get("is_duplication", False):
                data_dict["linked_start"] = data_feature["linked_start"]
                data_dict["linked_end"] = data_feature["linked_end"]
                data_dict["length"] = data_feature["length"]

            # Add statistics if available
            has_stats = data_feature.get("has_stats", False)
            if has_stats:
                data_dict["mean"] = data_feature["mean"]
                data_dict["median"] = data_feature["median"]
                data_dict["std"] = data_feature["std"]
            
            source = ColumnDataSource(data=data_dict)

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
                # Use width from data if available (for RLE spans), otherwise use size parameter
                # Pass column name (no '@') for variable width per bar, or scalar for uniform width
                bar_width = 'width' if "width" in data_dict else size
                p.vbar(
                    x='x',
                    bottom=0,
                    top='y',
                    source=source,
                    color=color,
                    alpha=alpha,
                    width=bar_width,
                    legend_label = title
                )

    # Add hover with conditional tooltips based on whether statistics are available
    # Check if any feature in feature_dict has statistics, variable width, or is duplication
    has_any_stats = any(d.get("has_stats", False) for d in feature_dict)
    has_variable_width = any("width" in d and any(w != 1 for w in d["width"]) for d in feature_dict)
    is_duplication = any(d.get("is_duplication", False) for d in feature_dict)

    if is_duplication:
        # Duplication-specific tooltips
        tooltips = [
            ("First position", "@first_pos{0,0}"),
            ("Last position", "@last_pos{0,0}"),
            ("Linked start", "@linked_start{0,0}"),
            ("Linked end", "@linked_end{0,0}"),
            ("Length", "@length{0,0}"),
            ("Identity", "@y{0.1}%")
        ]
    elif has_variable_width:
        # For bars with spans, show first and last position
        tooltips = [
            ("First position", "@first_pos{0,0}"),
            ("Last position", "@last_pos{0,0}"),
            ("Value", "@y{0.0}")
        ]
    elif has_any_stats:
        tooltips = [
            ("Position", "@x{0,0}"),
            ("Value", "@y{0.0}"),
            ("Mean", "@mean{0.0}"),
            ("Median", "@median{0.0}"),
            ("Std", "@std{0.0}")
        ]
    else:
        tooltips = [("Position", "@x{0,0}"), ("Value", "@y{0.0}")]
    
    hover = HoverTool(tooltips=tooltips, mode='vline')
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

    wheel = WheelZoomTool(dimensions='width')
    p.add_tools(wheel)
    p.toolbar.active_scroll = wheel

    return p

### Function to get repeats (contig-level, sample-independent)
def get_repeats_data(cur, contig_id, variable_name="direct_repeats", xstart=None, xend=None):
    """Get repeats data for plotting, formatted for make_bokeh_subplot().

    Args:
        cur: DuckDB cursor
        contig_id: Contig ID
        variable_name: Either 'direct_repeats' or 'inverted_repeats'
        xstart: Optional start position for filtering (only fetch data intersecting this range)
        xend: Optional end position for filtering (only fetch data intersecting this range)

    Returns:
        List with one feature dict formatted for make_bokeh_subplot()
    """
    # Map variable name to table name
    table_map = {
        "direct_repeats": "Contig_DirectRepeats",
        "inverted_repeats": "Contig_InvertedRepeats",
    }
    table_name = table_map.get(variable_name)
    if not table_name:
        return []

    # Get variable info from database
    cur.execute(
        "SELECT \"Type\", Color, Alpha, Fill_alpha, \"Size\", Title "
        f"FROM Variable WHERE Variable_name='{variable_name}'"
    )
    var_row = cur.fetchone()
    if not var_row:
        return []  # Variable not found

    type_picked, color, alpha, fill_alpha, size, title = var_row

    # Check if table exists and has data
    try:
        # Build position filter clause
        position_filter = ""
        params = [contig_id]
        if xstart is not None and xend is not None:
            position_filter = " AND Position2 >= ? AND Position1 <= ?"
            params.extend([xstart, xend])
        
        cur.execute(
            f"SELECT Position1, Position2, Position1prime, Position2prime, Pident "
            f"FROM {table_name} WHERE Contig_id=?{position_filter} ORDER BY Position1",
            tuple(params)
        )
        rows = cur.fetchall()
    except Exception:
        return []  # Table doesn't exist or no data

    if not rows:
        return []

    # Clip repeat positions to requested range
    if xstart is not None and xend is not None:
        clipped_rows = []
        for row in rows:
            pos1, pos2, pos1p, pos2p, pident_int = row
            clipped_pos1 = max(pos1, xstart)
            clipped_pos2 = min(pos2, xend)
            clipped_rows.append((clipped_pos1, clipped_pos2, pos1p, pos2p, pident_int))
        rows = clipped_rows

    x_coords = []
    y_coords = []
    width_coords = []
    first_pos_coords = []
    last_pos_coords = []
    linked_start_coords = []
    linked_end_coords = []
    length_coords = []

    for row in rows:
        pos1, pos2, pos1p, pos2p, pident_int = row
        # Pident stored as INTEGER (×100), convert back to percentage
        pident = pident_int / 100.0 if pident_int is not None else 0.0
        length = abs(pos2 - pos1) + 1
        midpoint = (pos1 + pos2) / 2.0

        x_coords.append(midpoint)
        y_coords.append(pident)
        width_coords.append(length)
        first_pos_coords.append(min(pos1, pos2))
        last_pos_coords.append(max(pos1, pos2))
        linked_start_coords.append(pos1p)
        linked_end_coords.append(pos2p)
        length_coords.append(length)

    return [{
        "type": type_picked,
        "color": color,
        "alpha": alpha,
        "fill_alpha": fill_alpha,
        "size": size,
        "title": title,
        "x": x_coords,
        "y": y_coords,
        "width": width_coords,
        "first_pos": first_pos_coords,
        "last_pos": last_pos_coords,
        "linked_start": linked_start_coords,
        "linked_end": linked_end_coords,
        "length": length_coords,
        "has_stats": False,
        "is_duplication": True,  # Flag for special tooltip handling
    }]


### Function to get GC content (contig-level, sample-independent)
def get_gc_content_data(cur, contig_id, xstart=None, xend=None):
    """Get GC content data for plotting, formatted for make_bokeh_subplot().

    Args:
        cur: DuckDB cursor
        contig_id: Contig ID
        xstart: Optional start position for filtering (only fetch data intersecting this range)
        xend: Optional end position for filtering (only fetch data intersecting this range)

    Returns:
        List with one feature dict formatted for make_bokeh_subplot()
    """
    # Get variable info from database
    cur.execute(
        "SELECT \"Type\", Color, Alpha, Fill_alpha, \"Size\", Title "
        "FROM Variable WHERE Variable_name='gc_content'"
    )
    var_row = cur.fetchone()
    if not var_row:
        return []  # Variable not found

    type_picked, color, alpha, fill_alpha, size, title = var_row

    # Check if table exists and has data
    try:
        # Build position filter clause
        position_filter = ""
        params = [contig_id]
        if xstart is not None and xend is not None:
            position_filter = " AND Last_position >= ? AND First_position <= ?"
            params.extend([xstart, xend])
        
        cur.execute(
            f"SELECT First_position, Last_position, GC_percentage "
            f"FROM Contig_GCContent WHERE Contig_id=?{position_filter} ORDER BY First_position",
            tuple(params)
        )
        rows = cur.fetchall()
    except Exception:
        return []  # Table doesn't exist or no data

    if not rows:
        return []

    # Clip GC content positions to requested range
    if xstart is not None and xend is not None:
        clipped_rows = []
        for row in rows:
            first_pos, last_pos, gc_pct = row
            clipped_first = max(first_pos, xstart)
            clipped_last = min(last_pos, xend)
            clipped_rows.append((clipped_first, clipped_last, gc_pct))
        rows = clipped_rows

    x_coords = []
    y_coords = []
    first_pos_coords = []
    last_pos_coords = []

    # Create step plot by adding points at both start and end of each window
    # This maintains the GC percentage value across the entire sliding window range
    for first_pos, last_pos, gc_pct in rows:
        # Add point at start of window
        x_coords.append(first_pos)
        y_coords.append(gc_pct)
        first_pos_coords.append(first_pos)
        last_pos_coords.append(last_pos)
        
        # Add point at end of window (creates horizontal line across window)
        x_coords.append(last_pos)
        y_coords.append(gc_pct)
        first_pos_coords.append(first_pos)
        last_pos_coords.append(last_pos)

    return [{
        "type": type_picked,
        "color": color,
        "alpha": alpha,
        "fill_alpha": fill_alpha,
        "size": size,
        "title": title,
        "x": x_coords,
        "y": y_coords,
        "first_pos": first_pos_coords,
        "last_pos": last_pos_coords,
        "has_stats": False,
    }]


def merge_rle_segments(plus_rows, minus_rows):
    """Merge two RLE-encoded arrays by summing overlapping values.

    Args:
        plus_rows: List of (first, last, value) tuples from plus strand
        minus_rows: List of (first, last, value) tuples from minus strand

    Returns:
        List of (first, last, combined_value) tuples with merged segments
    """
    # Collect all boundaries from both sources
    boundaries = set()
    for first, last, _ in plus_rows:
        boundaries.add(first)
        boundaries.add(last + 1)
    for first, last, _ in minus_rows:
        boundaries.add(first)
        boundaries.add(last + 1)

    if not boundaries:
        return []

    # Create segments from sorted boundaries
    sorted_bounds = sorted(boundaries)
    result = []

    for i in range(len(sorted_bounds) - 1):
        seg_start = sorted_bounds[i]
        seg_end = sorted_bounds[i + 1] - 1

        # Sum values from both strands that cover this segment
        value = 0
        for first, last, val in plus_rows:
            if first <= seg_start and last >= seg_end:
                value += val
        for first, last, val in minus_rows:
            if first <= seg_start and last >= seg_end:
                value += val

        result.append((seg_start, seg_end, value))

    return result


### Function to get features of one variable
def get_feature_data(cur, feature, contig_id, sample_id, xstart=None, xend=None):
    """Get feature data for plotting.

    Args:
        cur: DuckDB cursor
        feature: Feature name to query
        contig_id: Contig ID
        sample_id: Sample ID
        xstart: Optional start position for filtering (only fetch data intersecting this range)
        xend: Optional end position for filtering (only fetch data intersecting this range)
    """
    # Get rendering info from Variable table (Type and Size are quoted - reserved words in DuckDB)
    cur.execute(
        "SELECT \"Type\", Color, Alpha, Fill_alpha, \"Size\", Title, Feature_table_name "
        "FROM Variable WHERE Subplot=?",
        (feature,)
    )
    rows = cur.fetchall()

    # list_feature_dict has several elements if multiple variables share the same subplot
    # example the clippings (right vs left)
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

        # Check if this feature has statistics columns
        features_with_stats = ["left_clippings", "right_clippings", "insertions"]
        has_stats = feature_table in [f"Feature_{f}" for f in features_with_stats]

        # Check if this feature stores scaled values (stored as INTEGER ×100)
        scaled_features = ["Feature_tau", "Feature_mapq"]
        is_scaled = feature_table in scaled_features

        # Special handling for primary_reads: combine strand tables in Python
        # This avoids the OOM issues from the complex VIEW that computes on-the-fly
        if feature_table == "Feature_primary_reads":
            # Build position filter clause
            position_filter = ""
            params = [sample_id, contig_id]
            if xstart is not None and xend is not None:
                position_filter = " AND Last_position >= ? AND First_position <= ?"
                params.extend([xstart, xend])
            
            cur.execute(
                f"SELECT First_position, Last_position, Value FROM Feature_primary_reads_plus_only "
                f"WHERE Sample_id=? AND Contig_id=?{position_filter} ORDER BY First_position",
                tuple(params)
            )
            plus_rows = cur.fetchall()

            cur.execute(
                f"SELECT First_position, Last_position, Value FROM Feature_primary_reads_minus_only "
                f"WHERE Sample_id=? AND Contig_id=?{position_filter} ORDER BY First_position",
                tuple(params)
            )
            minus_rows = cur.fetchall()

            # Clip plus and minus rows before merging to ensure boundaries are respected
            if xstart is not None and xend is not None:
                clipped_plus = []
                for first_pos, last_pos, value in plus_rows:
                    clipped_first = max(first_pos, xstart)
                    clipped_last = min(last_pos, xend)
                    clipped_plus.append((clipped_first, clipped_last, value))
                plus_rows = clipped_plus

                clipped_minus = []
                for first_pos, last_pos, value in minus_rows:
                    clipped_first = max(first_pos, xstart)
                    clipped_last = min(last_pos, xend)
                    clipped_minus.append((clipped_first, clipped_last, value))
                minus_rows = clipped_minus

            data_rows = merge_rle_segments(plus_rows, minus_rows)
        # Query Feature_* table (RLE format: First_position, Last_position, Value, and optionally Mean, Median, Std)
        elif has_stats:
            # Build position filter clause
            position_filter = ""
            params = [sample_id, contig_id]
            if xstart is not None and xend is not None:
                position_filter = " AND Last_position >= ? AND First_position <= ?"
                params.extend([xstart, xend])
            
            cur.execute(
                f"SELECT First_position, Last_position, Value, Mean, Median, Std FROM {feature_table} "
                f"WHERE Sample_id=? AND Contig_id=?{position_filter} ORDER BY First_position",
                tuple(params)
            )
            data_rows = cur.fetchall()
        else:
            # Build position filter clause
            position_filter = ""
            params = [sample_id, contig_id]
            if xstart is not None and xend is not None:
                position_filter = " AND Last_position >= ? AND First_position <= ?"
                params.extend([xstart, xend])
            
            cur.execute(
                f"SELECT First_position, Last_position, Value FROM {feature_table} "
                f"WHERE Sample_id=? AND Contig_id=?{position_filter} ORDER BY First_position",
                tuple(params)
            )
            data_rows = cur.fetchall()

        # Clip feature positions to requested range
        if xstart is not None and xend is not None:
            clipped_rows = []
            for row in data_rows:
                if has_stats:
                    first_pos, last_pos, value, mean, median, std = row
                    clipped_first = max(first_pos, xstart)
                    clipped_last = min(last_pos, xend)
                    clipped_rows.append((clipped_first, clipped_last, value, mean, median, std))
                else:
                    first_pos, last_pos, value = row
                    clipped_first = max(first_pos, xstart)
                    clipped_last = min(last_pos, xend)
                    clipped_rows.append((clipped_first, clipped_last, value))
            data_rows = clipped_rows

        # Expand RLE runs into individual points for plotting
        x_coords = []
        y_coords = []
        mean_coords = []
        median_coords = []
        std_coords = []

        # Store widths for bars (needed for plotting spans)
        width_coords = []
        first_pos_coords = []
        last_pos_coords = []

        for row in data_rows:
            if has_stats:
                first_pos, last_pos, value, mean, median, std = row
            else:
                first_pos, last_pos, value = row
                mean = median = std = None

            # Scale back if this is a scaled feature (tau, mapq stored as ×100)
            if is_scaled:
                value = value / 100.0 if value is not None else None

            if type_picked == "bars":
                # For bars: use midpoint as x position and calculate width
                midpoint = (first_pos + last_pos) / 2.0
                width = last_pos - first_pos + 1
                x_coords.append(midpoint)
                y_coords.append(value)
                width_coords.append(width)
                first_pos_coords.append(first_pos)
                last_pos_coords.append(last_pos)
                if has_stats:
                    mean_coords.append(mean)
                    median_coords.append(median)
                    std_coords.append(std)
            else:
                # For curves: only need start and end points
                if first_pos == last_pos:
                    x_coords.append(first_pos)
                    y_coords.append(value)
                    if has_stats:
                        mean_coords.append(mean)
                        median_coords.append(median)
                        std_coords.append(std)
                else:
                    x_coords.extend([first_pos, last_pos])
                    y_coords.extend([value, value])
                    if has_stats:
                        mean_coords.extend([mean, mean])
                        median_coords.extend([median, median])
                        std_coords.extend([std, std])

        feature_dict["x"] = x_coords
        feature_dict["y"] = y_coords
        feature_dict["has_stats"] = has_stats
        if type_picked == "bars":
            feature_dict["width"] = width_coords
            feature_dict["first_pos"] = first_pos_coords
            feature_dict["last_pos"] = last_pos_coords
        if has_stats:
            feature_dict["mean"] = mean_coords
            feature_dict["median"] = median_coords
            feature_dict["std"] = std_coords

        # Only append if we have actual data points
        if x_coords:
            list_feature_dict.append(feature_dict)

    return list_feature_dict

### Function to generate the bokeh plot
def generate_bokeh_plot_per_sample(conn, list_features, contig_name, sample_name, xstart=None, xend=None, subplot_size=100, genbank_path=None):
    """Generate a Bokeh plot for a single sample.

    Args:
        conn: DuckDB connection
        list_features: List of features/modules to plot (can be mix of modules and individual features)
        contig_name: Name of the contig to plot
        sample_name: Name of the sample to plot
        xstart: Optional x-axis start position
        xend: Optional x-axis end position
        subplot_size: Height of each subplot in pixels
        genbank_path: Optional genbank file path (if provided, gene map will be shown)
    """
    cur = conn.cursor()

    # Get contig characteristics
    contig_id, locus_name, locus_size, annotation_tool = get_contig_info(cur, contig_name)
    print(f"Locus {locus_name} validated ({locus_size} bp)", flush=True)

    # --- Main gene annotation plot (only if genbank provided) ---
    shared_xrange = Range1d(-1000, locus_size+1000)
    if xstart is not None and xend is not None:
        shared_xrange.start = xstart
        shared_xrange.end = xend

    annotation_fig = make_bokeh_genemap(conn, contig_id, locus_name, locus_size, annotation_tool, subplot_size, shared_xrange, xstart, xend) if genbank_path else None

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
    requested_features, include_repeats, include_gc_content = parse_requested_features(list_features)

    # Add Repeats subplots if requested (contig-level, sample-independent)
    if include_repeats:
        # Direct repeats
        try:
            direct_feature_dict = get_repeats_data(cur, contig_id, "direct_repeats", xstart, xend)
            if direct_feature_dict:
                direct_subplot = make_bokeh_subplot(direct_feature_dict, subplot_size, shared_xrange)
                if direct_subplot is not None:
                    subplots.append(direct_subplot)
        except Exception as e:
            print(f"Error processing Direct Repeats: {e}", flush=True)
        # Inverted repeats
        try:
            inverted_feature_dict = get_repeats_data(cur, contig_id, "inverted_repeats", xstart, xend)
            if inverted_feature_dict:
                inverted_subplot = make_bokeh_subplot(inverted_feature_dict, subplot_size, shared_xrange)
                if inverted_subplot is not None:
                    subplots.append(inverted_subplot)
        except Exception as e:
            print(f"Error processing Inverted Repeats: {e}", flush=True)

    # Add GC content subplot if requested (contig-level, sample-independent)
    if include_gc_content:
        try:
            gc_feature_dict = get_gc_content_data(cur, contig_id, xstart, xend)
            if gc_feature_dict:
                gc_subplot = make_bokeh_subplot(gc_feature_dict, subplot_size, shared_xrange)
                if gc_subplot is not None:
                    subplots.append(gc_subplot)
        except Exception as e:
            print(f"Error processing GC Content: {e}", flush=True)

    # Add other requested features
    for feature in requested_features:
        try:
            list_feature_dict = get_feature_data(cur, feature, contig_id, sample_id, xstart, xend)
            subplot_feature = make_bokeh_subplot(list_feature_dict, subplot_size, shared_xrange, feature_name=feature)
            if subplot_feature is not None:
                subplots.append(subplot_feature)
        except Exception as e:
            print(f"Error processing feature '{feature}': {e}", flush=True)

    # --- Combine all figures in a single grid with one shared toolbar ---
    if annotation_fig:
        if not subplots:
            grid = gridplot([[annotation_fig]], merge_tools=True, sizing_mode='stretch_width')
        else:
            all_plots = [annotation_fig] + subplots
            grid = gridplot([[p] for p in all_plots], merge_tools=True, sizing_mode='stretch_width')
    else:
        # No gene map - just show subplots
        if not subplots:
            raise ValueError("No plots to display")
        grid = gridplot([[p] for p in subplots], merge_tools=True, sizing_mode='stretch_width')

    return grid

def save_html_plot_per_sample(db_path, list_features, contig_name, sample_name, subplot_size, output_filename, genbank_path=None):
    # --- Save interactive HTML plot ---
    output_file(filename = output_filename)
    conn = duckdb.connect(db_path, read_only=True)
    try:
        grid = generate_bokeh_plot_per_sample(conn, list_features, contig_name, sample_name, subplot_size=subplot_size, genbank_path=genbank_path)
        save(grid)
    finally:
        conn.close()

### Parsing features
def parse_requested_features(list_features):
    """Parse requested features, expanding modules to individual features.

    Accepts a mix of module names and individual feature names.
    Module names are case-insensitive and can include:
    - "coverage" or "Coverage" -> primary_reads, secondary_reads, supplementary_reads
    - "phagetermini" or "Phage termini" -> coverage_reduced, reads_starts, reads_ends, tau + Repeats
    - "assemblycheck" or "Assembly check" -> all assembly check features
    - "genome" or "Genome" -> Repeats + GC content

    Returns tuple of (deduplicated list of individual feature names, include_repeats bool, include_gc_content bool).
    """
    features = []
    include_repeats = False
    include_gc_content = False

    for item in list_features:
        item_lower = item.lower().strip()

        # Module: Genome
        if item_lower in ["genome"]:
            include_repeats = True
            include_gc_content = True
        # Module: Coverage
        elif item_lower in ["coverage"]:
            features.extend(["Primary alignments", "Other alignments", "Other alignments"])
        # Module: Phage termini / phagetermini
        elif item_lower in ["phage termini", "phagetermini", "phage_termini"]:
            features.extend(["Coverage reduced", "Reads termini", "Tau"])
            include_repeats = True
        # Module: Assembly check / assemblycheck
        elif item_lower in ["assembly check", "assemblycheck", "assembly_check"]:
            features.extend(["Clippings", "Indels", "Mismatches", "Read lengths", "Insert sizes", "Bad orientations"])
        # Handle "Repeats" specifically (also accept legacy "duplications")
        elif item_lower in ["repeats", "repeat", "duplications", "duplication"]:
            include_repeats = True
        # Handle "GC content" specifically
        elif item_lower in ["gc_content", "gc content", "gccontent", "gc"]:
            include_gc_content = True
        # Individual feature
        else:
            features.append(item)

    # Deduplicate while preserving order
    seen = set()
    deduped_features = [f for f in features if not (f in seen or seen.add(f))]
    return deduped_features, include_repeats, include_gc_content

### Main function
def add_plot_per_sample_args(parser):
    parser.add_argument("--db", required=True, help="Path to DuckDB database file")
    parser.add_argument("--variables", required=True, help="List of variables or full modules to compute (comma-separated) (options allowed for modules: Coverage, Misalignment, Long-reads, Paired-reads, Phage termini)")
    parser.add_argument("--contig", required=True, help="Name of the contig to plot")
    parser.add_argument("--sample", required=True, help="Name of the sample to plot")
    parser.add_argument("-g", "--genbank", help="Path to genbank file (optional; if provided, gene map will be plotted)")
    parser.add_argument("--html", required=False, default="thebigbam_per_sample.html", help="Name for output html files. A bokeh server will be started if not provided")
    parser.add_argument("--subplot_height", required=False, default=130, help="Height of each subplot (in pixels)")

def run_plot_per_sample(args):
    db_path = args.db
    list_features = args.variables.split(",")
    contig_name = args.contig
    sample_name = args.sample
    output_filename = args.html
    subplot_size = int(args.subplot_height)
    genbank_path = getattr(args, 'genbank', None)

    print(f"Saving static HTML to {output_filename}...", flush=True)
    save_html_plot_per_sample(db_path, list_features, contig_name, sample_name, subplot_size, output_filename, genbank_path)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    add_plot_per_sample_args(parser)
    args = parser.parse_args()
    run_plot_per_sample(args)
