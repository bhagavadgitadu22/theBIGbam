import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from bokeh.models import Range1d, ColumnDataSource, HoverTool, WheelZoomTool, NumeralTickFormatter, Label, TapTool
from bokeh.layouts import gridplot
from bokeh.plotting import figure
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
TYPE_COLORS = {
    # grey for protein-coding genes
    'CDS': '#cccccc',
    'mRNA': '#777777',
    # green for RNA genes
    'tRNA': '#66c2a5',
    'tmRNA': '#99d8c9',
    'rRNA': '#238b45',
    'ncRNA': '#33a02c',
    'precursor_RNA': '#a1d99b',
    'misc_RNA': '#74c476',
    # orange for regulatory / gene structure
    'exon': '#fdae61',
    "5'UTR": '#fee08b',
    "3'UTR": '#f46d43',
    # purple for genome architecture & mobility
    'repeat_region': '#6a3d9a',
    'mobile_element': '#cab2d6',
    # other features
    'misc_feature': '#3c5bfe',
    'gap': '#e5049c',
    'pseudogene': "#e31a1c"
    # features not listed here will get black color by default
}

class CustomTranslator(BiopythonTranslator):

    # Track seen unknown feature types at the class level
    _seen_unknown_types = set()

    def compute_feature_color(self, feature):
        type_feature = feature.type

        if type_feature == "CDS":
            use_phage_colors = feature.qualifiers.get("use_phage_colors", False)

            # Use phage colors if checkbox is checked
            if use_phage_colors:
                # Get the function field safely
                function = feature.qualifiers.get("function")
                if isinstance(function, list):  # Biopython often stores qualifiers as lists
                    function = function[0] if function else None

                if not isinstance(function, str):  # Missing or wrong type
                    return "#cccccc"

                function = function.lower()

                for key, color in PHAROKKA_CDS_COLORS.items():
                    if key in function:
                        return color

            return "#cccccc"

        else:
            if type_feature not in TYPE_COLORS:
                if type_feature not in CustomTranslator._seen_unknown_types:
                    print("Unknown type of feature:", type_feature, flush=True)
                    CustomTranslator._seen_unknown_types.add(type_feature)
                return "#000000"
            return TYPE_COLORS.get(type_feature, "#cccccc")

    def compute_feature_label(self, feature):
        return None  # fallback to None if missing or invalid
    
    def compute_feature_html(self, feature):
        type_feature = feature.type
        if type_feature == "CDS":
            return feature.qualifiers.get("product", [])
        return feature.type
        
    
### Plotting functions
def get_contig_info(cur, contig_name):
    cur.execute("SELECT Contig_id, Contig_name, Contig_length FROM Contig WHERE Contig_name=?", (contig_name,))
    row = cur.fetchone()
    if row is None:
        raise ValueError(f"Contig not found: {contig_name}")
    return row

def make_bokeh_genemap(conn, contig_id, locus_name, locus_size, subplot_size, shared_xrange, xstart=None, xend=None, feature_types=None, use_phage_colors=False, plot_isoforms=True):
    cur = conn.cursor()

    # Build position filter clause for annotations
    position_filter = ""
    params = [contig_id]
    if xstart is not None and xend is not None:
        position_filter = " AND \"End\" >= ? AND \"Start\" <= ?"
        params.extend([xstart, xend])

    # Build feature type filter
    type_filter = ""
    if feature_types:
        placeholders = ','.join('?' * len(feature_types))
        type_filter = f' AND "Type" IN ({placeholders})'
        params.extend(feature_types)

    # When plot_isoforms is False, filter to show only longest isoform per (locus_tag, Type) pair
    # Features without locus_tag always display (Longest_isoform is NULL for them)
    if not plot_isoforms:
        # Use pre-computed Longest_isoform boolean column for efficient filtering
        isoform_filter = " AND (Locus_tag IS NULL OR Longest_isoform = true)"
        query = f'SELECT "Start", "End", Strand, "Type", Product, "Function", Phrog, Locus_tag FROM Contig_annotation WHERE Contig_id=?{position_filter}{type_filter}{isoform_filter}'
    else:
        query = f'SELECT "Start", "End", Strand, "Type", Product, "Function", Phrog, Locus_tag FROM Contig_annotation WHERE Contig_id=?{position_filter}{type_filter}'
    
    cur.execute(query, tuple(params))
    seq_ann_rows = cur.fetchall()

    sequence_annotations = []
    for start, end, strand, ftype, product, function, phrog, locus_tag in seq_ann_rows:
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
        if locus_tag:
            qualifiers['locus_tag'] = locus_tag
        qualifiers['use_phage_colors'] = use_phage_colors
        feat = SeqFeature(location=floc, type=ftype, qualifiers=qualifiers)
        sequence_annotations.append(feat)

    # Return None if no features to plot (avoids empty sequence error in dna_features_viewer)
    if not sequence_annotations:
        return None

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

def make_bokeh_subplot(feature_dict, height, x_range, sample_title=None):
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
    has_data = bool(feature_dict) and any(any(y != 0 for y in d["y"]) for d in feature_dict)
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
                # Set y1 to the minimum of zero and the minimum value in y to allow negative values
                p.varea(
                    x='x',
                    y1=0,
                    y2='y',
                    source=source,
                    fill_color=color,
                    fill_alpha=fill_alpha,
                    legend_label=title
                )
                p.line(
                    x='x',
                    y='y',
                    source=source,
                    line_color=color,
                    line_alpha=alpha,
                    line_width=size,
                    legend_label=title
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
            ("Identity", "@y{0.01}%")
        ]
    elif has_variable_width:
        # For bars with spans, show first and last position
        tooltips = [
            ("First position", "@first_pos{0,0}"),
            ("Last position", "@last_pos{0,0}"),
            ("Value", "@y{0.00}")
        ]
    elif has_any_stats:
        tooltips = [
            ("Position", "@x{0,00}"),
            ("Value", "@y{0.00}"),
            ("Mean", "@mean{0.00}"),
            ("Median", "@median{0.00}"),
            ("Std", "@std{0.00}")
        ]
    else:
        tooltips = [("Position", "@x{0,0}"), ("Value", "@y{0.00}")]
    
    hover = HoverTool(tooltips=tooltips, mode='vline')
    p.add_tools(hover)

    # A clean style like your matplotlib setup
    p.toolbar.logo = None
    p.xgrid.visible = False

    p.y_range.start = min(0, *(y for d in feature_dict for y in d["y"]))
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

### Function to render DNA sequence as colored rectangles
def make_bokeh_sequence_subplot(conn, contig_name, xstart, xend, height, x_range):
    """Create a subplot showing colored nucleotide rectangles for a genomic region.

    Returns None if no sequence data is available
    or the Contig_sequence table doesn't exist.

    Args:
        conn: DuckDB connection
        contig_name: Name of the contig
        xstart: Start position (0 = full genome, or 1-based genome position)
        xend: End position
        height: Height of the subplot in pixels
        x_range: Shared x_range from other subplots
    """
    try:
        cur = conn.cursor()

        # Clamp start to 1 (no nucleotide at position 0)
        seq_start = max(xstart, 1)

        # Query only the needed substring (SUBSTR is 1-based)
        cur.execute(
            "SELECT SUBSTR(cs.Sequence, ?, ? - ? + 1) "
            "FROM Contig_sequence cs "
            "JOIN Contig c ON cs.Contig_id = c.Contig_id "
            "WHERE c.Contig_name = ?",
            (seq_start, xend, seq_start, contig_name)
        )
        row = cur.fetchone()
        if row is None or row[0] is None:
            return None

        seq = row[0]
        if not seq:
            return None

        # Map nucleotides to colors
        color_map = {
            'A': '#d62728', 'a': '#d62728',
            'T': '#2ca02c', 't': '#2ca02c',
            'G': '#ff7f0e', 'g': '#ff7f0e',
            'C': '#1f77b4', 'c': '#1f77b4',
        }

        positions = []
        colors = []
        nucleotides = []
        for i, nt in enumerate(seq):
            positions.append(seq_start + i)          # 1-based position
            colors.append(color_map.get(nt, '#999999'))
            nucleotides.append(nt.upper())

        source = ColumnDataSource(data=dict(
            left=[p - 0.5 for p in positions],       # Center quad on position
            right=[p + 0.5 for p in positions],
            bottom=[0] * len(positions),
            top=[1] * len(positions),
            color=colors,
            nucleotide=nucleotides,
            position=positions,
        ))

        p = figure(
            height=height,
            x_range=x_range,
            y_range=Range1d(0, 1),
            tools="xpan,reset,save"
        )

        p.quad(
            left='left', right='right', bottom='bottom', top='top',
            color='color', source=source, line_color=None
        )

        hover = HoverTool(tooltips=[
            ("Position", "@position{0,0}"),
            ("Nucleotide", "@nucleotide"),
        ])
        p.add_tools(hover)

        # Match styling from make_bokeh_subplot
        p.toolbar.logo = None
        p.xaxis.formatter = NumeralTickFormatter(format="0,0")
        p.yaxis.visible = False
        p.xgrid.visible = False
        p.ygrid.visible = False
        p.outline_line_color = None
        p.min_border_left = 40
        p.min_border_right = 10

        wheel = WheelZoomTool(dimensions='width')
        p.add_tools(wheel)
        p.toolbar.active_scroll = wheel

        return p

    except Exception:
        return None


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
        "direct_repeats": "Contig_directRepeats",
        "inverted_repeats": "Contig_invertedRepeats",
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
        # Both arrays are sorted by first position, so we can break early
        value = 0
        for first, last, val in plus_rows:
            if first > seg_end:
                break
            if first <= seg_start and last >= seg_end:
                value += val
        for first, last, val in minus_rows:
            if first > seg_end:
                break
            if first <= seg_start and last >= seg_end:
                value += val

        result.append((seg_start, seg_end, value))

    return result


### Function to get variable metadata (rendering info from Variable table)
def get_variable_metadata(cur, feature):
    """Get rendering metadata from the Variable table for a given subplot/feature.

    Args:
        cur: DuckDB cursor
        feature: Subplot name to query

    Returns:
        List of tuples: (Type, Color, Alpha, Fill_alpha, Size, Title, Feature_table_name)
    """
    cur.execute(
        "SELECT \"Type\", Color, Alpha, Fill_alpha, \"Size\", Title, Feature_table_name "
        "FROM Variable WHERE Subplot=?",
        (feature,)
    )
    return cur.fetchall()


def get_variable_metadata_batch(cur, subplot_list):
    """Batch fetch variable metadata for multiple subplots in one query.

    Returns:
        Dict mapping subplot_name -> list of metadata tuples
        (same tuple format as get_variable_metadata returns)
    """
    if not subplot_list:
        return {}
    placeholders = ', '.join(['?'] * len(subplot_list))
    cur.execute(
        f'SELECT Subplot, "Type", Color, Alpha, Fill_alpha, "Size", Title, Feature_table_name '
        f'FROM Variable WHERE Subplot IN ({placeholders}) '
        f'ORDER BY Module_order',
        tuple(subplot_list)
    )
    result = {}
    for row in cur.fetchall():
        subplot_name = row[0]
        metadata_tuple = row[1:]  # matches get_variable_metadata format
        result.setdefault(subplot_name, []).append(metadata_tuple)
    return result


### Function to get features of one variable
def get_feature_data(cur, feature, contig_id, sample_id, xstart=None, xend=None, variable_metadata=None):
    """Get feature data for plotting.

    Args:
        cur: DuckDB cursor
        feature: Feature name to query
        contig_id: Contig ID
        sample_id: Sample ID
        xstart: Optional start position for filtering (only fetch data intersecting this range)
        xend: Optional end position for filtering (only fetch data intersecting this range)
        variable_metadata: Optional cached result from get_variable_metadata(); avoids re-querying Variable table
    """
    # Get rendering info from Variable table (Type and Size are quoted - reserved words in DuckDB)
    if variable_metadata is not None:
        rows = variable_metadata
    else:
        rows = get_variable_metadata(cur, feature)

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
        # - tau and mapq are stored as value × 100
        # - gc_skew is stored as value × 100 (range: -100 to +100 representing -1.0 to +1.0)
        scaled_features = ["Feature_tau", "Feature_mapq", "Contig_GCSkew"]
        is_scaled = feature_table in scaled_features

        # Detect contig-level table (no Sample_id column)
        is_contig_table = feature_table.startswith("Contig_")

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
            if is_contig_table:
                params = [contig_id]
            else:
                params = [sample_id, contig_id]
            if xstart is not None and xend is not None:
                position_filter = " AND Last_position >= ? AND First_position <= ?"
                params.extend([xstart, xend])

            if is_contig_table:
                cur.execute(
                    f"SELECT First_position, Last_position, Value FROM {feature_table} "
                    f"WHERE Contig_id=?{position_filter} ORDER BY First_position",
                    tuple(params)
                )
            else:
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
                    if clipped_first <= clipped_last:
                        clipped_rows.append((clipped_first, clipped_last, value, mean, median, std))
                else:
                    first_pos, last_pos, value = row
                    clipped_first = max(first_pos, xstart)
                    clipped_last = min(last_pos, xend)
                    if clipped_first <= clipped_last:
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


### Function to get features for multiple samples in a single batch
def _expand_rle_rows(data_rows, type_picked, has_stats, is_scaled, xstart, xend):
    """Expand RLE rows into plot coordinates (shared logic for single and batch).

    Args:
        data_rows: List of tuples (First_position, Last_position, Value[, Mean, Median, Std])
        type_picked: Plot type ('bars' or 'curve')
        has_stats: Whether rows include statistics columns
        is_scaled: Whether values need to be divided by 100
        xstart: Optional start position for clipping
        xend: Optional end position for clipping

    Returns:
        Dict with x, y, and optional width/stats coordinate lists, or None if no data
    """
    # Clip feature positions to requested range
    if xstart is not None and xend is not None:
        clipped_rows = []
        for row in data_rows:
            if has_stats:
                first_pos, last_pos, value, mean, median, std = row
                clipped_first = max(first_pos, xstart)
                clipped_last = min(last_pos, xend)
                if clipped_first <= clipped_last:
                    clipped_rows.append((clipped_first, clipped_last, value, mean, median, std))
            else:
                first_pos, last_pos, value = row
                clipped_first = max(first_pos, xstart)
                clipped_last = min(last_pos, xend)
                if clipped_first <= clipped_last:
                    clipped_rows.append((clipped_first, clipped_last, value))
        data_rows = clipped_rows

    x_coords = []
    y_coords = []
    mean_coords = []
    median_coords = []
    std_coords = []
    width_coords = []
    first_pos_coords = []
    last_pos_coords = []

    for row in data_rows:
        if has_stats:
            first_pos, last_pos, value, mean, median, std = row
        else:
            first_pos, last_pos, value = row
            mean = median = std = None

        if is_scaled:
            value = value / 100.0 if value is not None else None

        if type_picked == "bars":
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

    if not x_coords:
        return None

    result = {"x": x_coords, "y": y_coords, "has_stats": has_stats}
    if type_picked == "bars":
        result["width"] = width_coords
        result["first_pos"] = first_pos_coords
        result["last_pos"] = last_pos_coords
    if has_stats:
        result["mean"] = mean_coords
        result["median"] = median_coords
        result["std"] = std_coords
    return result


def get_feature_data_batch(cur, feature, contig_id, sample_ids, xstart=None, xend=None, variable_metadata=None):
    """Get feature data for multiple samples in a single batch query.

    Instead of running N separate queries (one per sample), runs one query with
    WHERE Sample_id IN (...) and partitions results in Python.

    Args:
        cur: DuckDB cursor
        feature: Subplot name to query
        contig_id: Contig ID
        sample_ids: List of sample IDs to fetch data for
        xstart: Optional start position for filtering
        xend: Optional end position for filtering
        variable_metadata: Optional cached result from get_variable_metadata()

    Returns:
        Dict mapping sample_id to list_feature_dict (same format as get_feature_data returns)
    """
    if variable_metadata is not None:
        rows = variable_metadata
    else:
        rows = get_variable_metadata(cur, feature)

    # Result: {sample_id: list_feature_dict}
    result = {sid: [] for sid in sample_ids}

    if not sample_ids:
        return result

    # Build the IN clause placeholders
    placeholders = ", ".join(["?"] * len(sample_ids))

    for row in rows:
        type_picked, color, alpha, fill_alpha, size, title, feature_table = row

        features_with_stats = ["left_clippings", "right_clippings", "insertions"]
        has_stats = feature_table in [f"Feature_{f}" for f in features_with_stats]
        scaled_features = ["Feature_tau", "Feature_mapq"]
        is_scaled = feature_table in scaled_features

        # Detect contig-level table (no Sample_id column)
        is_contig_table = feature_table.startswith("Contig_")

        # Build position filter
        position_filter = ""
        extra_params = []
        if xstart is not None and xend is not None:
            position_filter = " AND Last_position >= ? AND First_position <= ?"
            extra_params = [xstart, xend]

        if feature_table == "Feature_primary_reads":
            # Batch query for plus strand
            params = list(sample_ids) + [contig_id] + extra_params
            cur.execute(
                f"SELECT Sample_id, First_position, Last_position, Value "
                f"FROM Feature_primary_reads_plus_only "
                f"WHERE Sample_id IN ({placeholders}) AND Contig_id=?{position_filter} "
                f"ORDER BY Sample_id, First_position",
                tuple(params)
            )
            all_plus = cur.fetchall()

            # Batch query for minus strand
            cur.execute(
                f"SELECT Sample_id, First_position, Last_position, Value "
                f"FROM Feature_primary_reads_minus_only "
                f"WHERE Sample_id IN ({placeholders}) AND Contig_id=?{position_filter} "
                f"ORDER BY Sample_id, First_position",
                tuple(params)
            )
            all_minus = cur.fetchall()

            # Group by sample_id
            plus_by_sample = {sid: [] for sid in sample_ids}
            for sid, first, last, val in all_plus:
                if sid in plus_by_sample:
                    plus_by_sample[sid].append((first, last, val))

            minus_by_sample = {sid: [] for sid in sample_ids}
            for sid, first, last, val in all_minus:
                if sid in minus_by_sample:
                    minus_by_sample[sid].append((first, last, val))

            # Merge and expand per sample
            for sid in sample_ids:
                data_rows = merge_rle_segments(plus_by_sample[sid], minus_by_sample[sid])
                expanded = _expand_rle_rows(data_rows, type_picked, has_stats, is_scaled, xstart, xend)
                if expanded is not None:
                    feature_dict = {
                        "type": type_picked, "color": color, "alpha": alpha,
                        "fill_alpha": fill_alpha, "size": size, "title": title,
                    }
                    feature_dict.update(expanded)
                    result[sid].append(feature_dict)

        elif has_stats:
            params = list(sample_ids) + [contig_id] + extra_params
            cur.execute(
                f"SELECT Sample_id, First_position, Last_position, Value, Mean, Median, Std "
                f"FROM {feature_table} "
                f"WHERE Sample_id IN ({placeholders}) AND Contig_id=?{position_filter} "
                f"ORDER BY Sample_id, First_position",
                tuple(params)
            )
            all_rows = cur.fetchall()

            # Group by sample_id
            rows_by_sample = {sid: [] for sid in sample_ids}
            for sid, first, last, val, mean, median, std in all_rows:
                if sid in rows_by_sample:
                    rows_by_sample[sid].append((first, last, val, mean, median, std))

            for sid in sample_ids:
                expanded = _expand_rle_rows(rows_by_sample[sid], type_picked, has_stats, is_scaled, xstart, xend)
                if expanded is not None:
                    feature_dict = {
                        "type": type_picked, "color": color, "alpha": alpha,
                        "fill_alpha": fill_alpha, "size": size, "title": title,
                    }
                    feature_dict.update(expanded)
                    result[sid].append(feature_dict)

        else:
            if is_contig_table:
                # Contig-level table: query once (no Sample_id), duplicate for all samples
                params = [contig_id] + extra_params
                cur.execute(
                    f"SELECT First_position, Last_position, Value FROM {feature_table} "
                    f"WHERE Contig_id=?{position_filter} ORDER BY First_position",
                    tuple(params)
                )
                contig_rows = cur.fetchall()
                expanded = _expand_rle_rows(contig_rows, type_picked, has_stats, is_scaled, xstart, xend)
                if expanded is not None:
                    for sid in sample_ids:
                        feature_dict = {
                            "type": type_picked, "color": color, "alpha": alpha,
                            "fill_alpha": fill_alpha, "size": size, "title": title,
                        }
                        feature_dict.update(expanded)
                        result[sid].append(feature_dict)
            else:
                params = list(sample_ids) + [contig_id] + extra_params
                cur.execute(
                    f"SELECT Sample_id, First_position, Last_position, Value "
                    f"FROM {feature_table} "
                    f"WHERE Sample_id IN ({placeholders}) AND Contig_id=?{position_filter} "
                    f"ORDER BY Sample_id, First_position",
                    tuple(params)
                )
                all_rows = cur.fetchall()

                # Group by sample_id
                rows_by_sample = {sid: [] for sid in sample_ids}
                for sid, first, last, val in all_rows:
                    if sid in rows_by_sample:
                        rows_by_sample[sid].append((first, last, val))

                for sid in sample_ids:
                    expanded = _expand_rle_rows(rows_by_sample[sid], type_picked, has_stats, is_scaled, xstart, xend)
                    if expanded is not None:
                        feature_dict = {
                            "type": type_picked, "color": color, "alpha": alpha,
                            "fill_alpha": fill_alpha, "size": size, "title": title,
                        }
                        feature_dict.update(expanded)
                        result[sid].append(feature_dict)

    return result


### Function to generate the bokeh plot
def generate_bokeh_plot_per_sample(conn, list_features, contig_name, sample_name, xstart=None, xend=None, subplot_size=100, genbank_path=None, feature_types=None, use_phage_colors=False, plot_isoforms=True, plot_sequence=False, same_y_scale=False):
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
        feature_types: Optional list of feature types to include in gene map (None = all)
        use_phage_colors: Whether to use phage color scheme for CDS features
        plot_isoforms: Whether to show all isoforms (True) or only longest per locus_tag (False)
    """
    cur = conn.cursor()

    # Get contig characteristics
    contig_id, locus_name, locus_size = get_contig_info(cur, contig_name)
    print(f"Locus {locus_name} validated ({locus_size} bp)", flush=True)

    # --- Main gene annotation plot (only if genbank provided) ---
    shared_xrange = Range1d(-1000, locus_size+1000)
    if xstart is not None and xend is not None:
        shared_xrange.start = xstart
        shared_xrange.end = xend

    annotation_fig = make_bokeh_genemap(conn, contig_id, locus_name, locus_size, subplot_size, shared_xrange, xstart, xend, feature_types=feature_types, use_phage_colors=use_phage_colors, plot_isoforms=plot_isoforms) if genbank_path else None

    # Get sample characteristics (optional – contig-level features work without a sample)
    sample_id = None
    if sample_name:
        cur.execute("SELECT Sample_id, Sample_name FROM Sample WHERE Sample_name=?", (sample_name,))
        row = cur.fetchone()
        if row is None:
            raise ValueError(f"Sample not found: {sample_name}")
        sample_id, sample_name = row
        print(f"Sample {sample_name} validated.", flush=True)

    # --- Add one subplot per feature requested ---
    # Requested features are variables like 'coverage', 'reads_starts', etc.
    subplots = []

    # --- Add sequence subplot right after annotation (top of data tracks) ---
    if plot_sequence:
        seq_subplot = make_bokeh_sequence_subplot(conn, contig_name, xstart, xend, subplot_size // 2, shared_xrange)
        if seq_subplot:
            subplots.append(seq_subplot)

    requested_features, include_repeats = parse_requested_features(list_features)

    # Add Repeats subplots if requested (contig-level, sample-independent)
    if include_repeats:
        repeats_feature_dicts = []
        try:
            direct_feature_dict = get_repeats_data(cur, contig_id, "direct_repeats", xstart, xend)
            if direct_feature_dict:
                repeats_feature_dicts.extend(direct_feature_dict)
        except Exception as e:
            print(f"Error processing Direct Repeats: {e}", flush=True)
        try:
            inverted_feature_dict = get_repeats_data(cur, contig_id, "inverted_repeats", xstart, xend)
            if inverted_feature_dict:
                repeats_feature_dicts.extend(inverted_feature_dict)
        except Exception as e:
            print(f"Error processing Inverted Repeats: {e}", flush=True)
        if repeats_feature_dicts:
            repeats_subplot = make_bokeh_subplot(repeats_feature_dicts, subplot_size, shared_xrange)
            if repeats_subplot is not None:
                subplots.append(repeats_subplot)

    # Separate contig-level features (GC content, GC skew) from sample-dependent features
    contig_level_features = ["GC content", "GC skew"]
    contig_features = [f for f in requested_features if f in contig_level_features]
    sample_features = [f for f in requested_features if f not in contig_level_features]

    # Add contig-level features (don't require sample_id)
    if contig_features:
        metadata_cache = get_variable_metadata_batch(cur, contig_features)
        for feature in contig_features:
            try:
                list_feature_dict = get_feature_data(
                    cur, feature, contig_id, None, xstart, xend,
                    variable_metadata=metadata_cache.get(feature)
                )
                subplot_feature = make_bokeh_subplot(list_feature_dict, subplot_size, shared_xrange)
                if subplot_feature is not None:
                    subplots.append(subplot_feature)
            except Exception as e:
                print(f"Error processing feature '{feature}': {e}", flush=True)

    # Subplot names whose y-axis should be relative to Primary alignments max
    PRIMARY_RELATIVE_SUBPLOTS = {
        "Primary alignments", "Other alignments", "Clippings", "Indels",
        "Mismatches", "Bad orientations", "Coverage reduced", "Reads termini",
        "Non-inward pairs", "Missing mates"
    }

    # Add sample-dependent features only when a sample is selected
    feature_subplots = []  # track (feature_name, figure, max_y) for same_y_scale
    if sample_id is not None and sample_features:
        # Pre-fetch metadata for all features in one query
        metadata_cache = get_variable_metadata_batch(cur, sample_features)

        # Add other requested features
        for feature in sample_features:
            try:
                list_feature_dict = get_feature_data(
                    cur, feature, contig_id, sample_id, xstart, xend,
                    variable_metadata=metadata_cache.get(feature)
                )
                subplot_feature = make_bokeh_subplot(list_feature_dict, subplot_size, shared_xrange)
                if subplot_feature is not None:
                    if same_y_scale:
                        max_y = max((y for d in list_feature_dict for y in d["y"]), default=0)
                        feature_subplots.append((feature, subplot_feature, max_y))
                    subplots.append(subplot_feature)
            except Exception as e:
                print(f"Error processing feature '{feature}': {e}", flush=True)

    # Post-process y-ranges for same_y_scale (per-sample view)
    if same_y_scale and sample_id is not None and feature_subplots:
        # Find primary_reads max y
        primary_max = 0
        for fname, fig, my in feature_subplots:
            if fname == "Primary alignments":
                primary_max = max(primary_max, my)

        # If Primary alignments not plotted, fetch its data to find the max
        if primary_max == 0:
            try:
                primary_data = get_feature_data(cur, "Primary alignments", contig_id, sample_id, xstart, xend)
                for d in primary_data:
                    if d["y"]:
                        primary_max = max(primary_max, max(d["y"]))
            except Exception:
                pass

        # Apply shared y-range to all primary-relative subplots
        if primary_max > 0:
            for fname, fig, my in feature_subplots:
                if fname in PRIMARY_RELATIVE_SUBPLOTS:
                    fig.y_range = Range1d(0, primary_max)

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

### Parsing features
def parse_requested_features(list_features):
    """Parse requested features, expanding modules to individual features.

    Accepts a mix of module names and individual feature names.
    Module names are case-insensitive and can include:
    - "coverage" or "Coverage" -> primary_reads, secondary_reads, supplementary_reads
    - "phagetermini" or "Phage termini" -> coverage_reduced, reads_starts, reads_ends, tau + Repeats
    - "assemblycheck" or "Assembly check" -> all assembly check features
    - "genome" or "Genome" -> Repeats + GC content + GC skew

    Returns tuple of (deduplicated list of individual feature names, include_repeats bool).
    Note: GC content and GC skew are returned as regular features in the list.
    """
    features = []
    include_repeats = False

    for item in list_features:
        item_lower = item.lower().strip()

        # Module: Genome
        if item_lower in ["genome"]:
            include_repeats = True
            features.extend(["GC content", "GC skew"])
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
        # Handle "GC content" specifically - add as regular feature
        elif item_lower in ["gc_content", "gc content", "gccontent", "gc"]:
            features.append("GC content")
        # Handle "GC skew" specifically - add as regular feature
        elif item_lower in ["gc_skew", "gc skew", "gcskew", "skew"]:
            features.append("GC skew")
        # Individual feature
        else:
            features.append(item)

    # Deduplicate while preserving order
    seen = set()
    deduped_features = [f for f in features if not (f in seen or seen.add(f))]
    return deduped_features, include_repeats
