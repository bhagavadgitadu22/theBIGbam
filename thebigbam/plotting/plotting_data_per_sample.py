import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from bokeh.models import Range1d, ColumnDataSource, HoverTool, WheelZoomTool, NumeralTickFormatter, TapTool
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

def make_bokeh_subplot(feature_dict, height, x_range, sample_title=None, show_tooltips=True):
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

            # Add repeat positions if available
            if "repeat_positions" in data_feature:
                data_dict["repeat_positions"] = data_feature["repeat_positions"]

            # Add statistics if available
            has_stats = data_feature.get("has_stats", False)
            if has_stats:
                data_dict["mean"] = data_feature["mean"]
                data_dict["median"] = data_feature["median"]
                data_dict["std"] = data_feature["std"]

            # Add sequence/prevalence if available
            has_sequences = data_feature.get("has_sequences", False)
            if has_sequences:
                data_dict["sequence"] = data_feature["sequence"]
                data_dict["sequence_prevalence"] = data_feature["sequence_prevalence"]
                # Pre-format prevalence as 0-1 float string for tooltip;
                # empty string when value is missing so Bokeh shows "" not "???"
                data_dict["sequence_prevalence_str"] = [
                    f"{v / 100.0:.2f}" if v is not None else ""
                    for v in data_feature["sequence_prevalence"]
                ]

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

    # Add hover tooltips only in full-resolution mode (window <= 10kb)
    if show_tooltips:
        has_any_stats = any(d.get("has_stats", False) for d in feature_dict)
        has_any_sequences = any(d.get("has_sequences", False) for d in feature_dict)
        has_variable_width = any("width" in d and any(w != 1 for w in d["width"]) for d in feature_dict)
        is_duplication = any(d.get("is_duplication", False) for d in feature_dict)
        has_repeat_positions = any("repeat_positions" in d for d in feature_dict)

        if is_duplication:
            tooltips = [
                ("First position", "@first_pos{0,0}"),
                ("Last position", "@last_pos{0,0}"),
                ("Linked start", "@linked_start{0,0}"),
                ("Linked end", "@linked_end{0,0}"),
                ("Length", "@length{0,0}"),
                ("Identity", "@y{0.01}%")
            ]
        elif has_repeat_positions:
            tooltips = [
                ("Position", "@x{0,0}"),
                ("Value", "@y{0.00}"),
                ("Repeat positions", "@repeat_positions"),
            ]
        elif has_variable_width:
            tooltips = [
                ("First position", "@first_pos{0,0}"),
                ("Last position", "@last_pos{0,0}"),
                ("Value", "@y{0.00}")
            ]
        elif has_any_stats:
            has_mean = any(
                any(m is not None for m in d.get("mean", []))
                for d in feature_dict if d.get("has_stats", False)
            )
            tooltips = [
                ("Position", "@x{0,00}"),
                ("Value", "@y{0.00}"),
            ]
            if has_mean:
                tooltips.append(("Mean", "@mean{0.00}"))
            tooltips.append(("Median clipping", "@median{0.00}") if not has_mean else ("Median", "@median{0.00}"))
            if has_mean:
                tooltips.append(("Std", "@std{0.00}"))
            if has_any_sequences:
                tooltips.append(("Sequence", "@sequence"))
                tooltips.append(("Prevalence", "@sequence_prevalence_str"))
        elif has_any_sequences:
            tooltips = [
                ("Position", "@x{0,0}"),
                ("Value", "@y{0.00}"),
                ("Sequence", "@sequence"),
                ("Prevalence", "@sequence_prevalence_str"),
            ]
        else:
            tooltips = [("Position", "@x{0,0}"), ("Value", "@y{0.00}")]

        hover = HoverTool(tooltips=tooltips, mode='vline')
        p.add_tools(hover)

    # A clean style like your matplotlib setup
    p.toolbar.logo = None
    p.xgrid.visible = False

    p.y_range.start = min(0, *(y for d in feature_dict for y in d["y"]))
    # Cap y-axis at 1 for relative-to-coverage features (values are ratios between 0 and 1)
    if all(d.get("is_relative_scaled", False) for d in feature_dict):
        p.y_range = Range1d(p.y_range.start, 1)
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


def aggregate_repeats_sweep(rows, original_positions=None):
    """Aggregate repeat intervals using a boundary-based sweep-line.

    Each repeat covers [min(pos1, pos2), max(pos1, pos2)] (first copy only).
    For each segment between consecutive boundaries, compute:
      - count: number of overlapping repeat intervals
      - max_identity: max pident among overlapping intervals

    Args:
        rows: List of (pos1, pos2, pos1prime, pos2prime, pident_int) tuples
        original_positions: Optional list of (orig_start, orig_end) per row,
            used to track which repeat regions contribute to each segment

    Returns:
        (count_segments, identity_segments) - each a list of
        (first, last, value, meta) tuples with adjacent same-value+meta
        segments merged (RLE compression). meta is a tuple of positions for
        count, or a single (start, end) for identity.
    """
    if not rows:
        return [], []

    # Build intervals: [start, end, pident, idx]
    intervals = []
    for idx, (pos1, pos2, _pos1p, _pos2p, pident_int) in enumerate(rows):
        start = min(pos1, pos2)
        end = max(pos1, pos2)
        pident = pident_int / 100.0 if pident_int is not None else 0.0
        intervals.append((start, end, pident, idx))

    # Collect all boundaries
    boundaries = set()
    for start, end, _, _ in intervals:
        boundaries.add(start)
        boundaries.add(end + 1)

    sorted_bounds = sorted(boundaries)
    if len(sorted_bounds) < 2:
        return [], []

    # Sweep: for each segment between consecutive boundaries, find overlapping intervals
    raw_count = []
    raw_identity = []
    for i in range(len(sorted_bounds) - 1):
        seg_start = sorted_bounds[i]
        seg_end = sorted_bounds[i + 1] - 1

        count = 0
        max_ident = 0.0
        max_ident_idx = -1
        overlapping_indices = []
        for iv_start, iv_end, pident, iv_idx in intervals:
            if iv_start > seg_end:
                continue
            if iv_end < seg_start:
                continue
            # Interval overlaps this segment
            count += 1
            overlapping_indices.append(iv_idx)
            if pident > max_ident:
                max_ident = pident
                max_ident_idx = iv_idx

        if count > 0:
            if original_positions is not None:
                count_meta = tuple(original_positions[j] for j in overlapping_indices)
                identity_meta = original_positions[max_ident_idx]
            else:
                count_meta = None
                identity_meta = None
            raw_count.append((seg_start, seg_end, count, count_meta))
            raw_identity.append((seg_start, seg_end, max_ident, identity_meta))

    # RLE merge: collapse adjacent segments with same value AND meta
    def _rle_merge(segments):
        if not segments:
            return []
        merged = [segments[0]]
        for first, last, value, meta in segments[1:]:
            prev_first, prev_last, prev_value, prev_meta = merged[-1]
            if value == prev_value and meta == prev_meta and first == prev_last + 1:
                merged[-1] = (prev_first, last, value, meta)
            else:
                merged.append((first, last, value, meta))
        return merged

    return _rle_merge(raw_count), _rle_merge(raw_identity)


def _format_count_positions(all_positions):
    """Format all overlapping repeat positions for count tooltip."""
    if not all_positions:
        return ""
    MAX_DISPLAY = 5
    parts = [f"({s:,}, {e:,})" for s, e in all_positions[:MAX_DISPLAY]]
    result = ", ".join(parts)
    if len(all_positions) > MAX_DISPLAY:
        result += f", ... (+{len(all_positions) - MAX_DISPLAY} more)"
    return result

def _format_identity_position(best_position):
    """Format the best-identity repeat position for identity tooltip."""
    if best_position is None:
        return ""
    return f"({best_position[0]:,}, {best_position[1]:,})"


def _segments_to_curve_coords(segments, position_formatter=None):
    """Convert (first, last, value[, meta]) segments to x/y arrays for varea plotting.

    Inserts zero-padding between non-contiguous segments so varea doesn't
    fill across gaps.

    Args:
        segments: List of (first, last, value) or (first, last, value, meta) tuples
        position_formatter: Optional callable that formats meta into a string.
            When provided, returns a third list of position strings.

    Returns:
        (x_list, y_list) when position_formatter is None
        (x_list, y_list, pos_strings) when position_formatter is provided
    """
    if not segments:
        if position_formatter is not None:
            return [], [], []
        return [], []

    x_coords = []
    y_coords = []
    pos_strings = [] if position_formatter is not None else None

    for i, seg in enumerate(segments):
        first, last, value = seg[0], seg[1], seg[2]
        meta = seg[3] if len(seg) > 3 else None

        # Insert zero gap before this segment if it's not contiguous with previous
        if i > 0:
            prev_last = segments[i - 1][1]
            if first > prev_last + 1:
                # Drop to zero at end of previous segment
                x_coords.append(prev_last + 1)
                y_coords.append(0)
                if pos_strings is not None:
                    pos_strings.append("")
                # Stay at zero until start of this segment
                x_coords.append(first - 1)
                y_coords.append(0)
                if pos_strings is not None:
                    pos_strings.append("")

        # Segment start and end
        formatted = position_formatter(meta) if position_formatter is not None else None
        x_coords.append(first)
        y_coords.append(value)
        if pos_strings is not None:
            pos_strings.append(formatted)
        if first != last:
            x_coords.append(last)
            y_coords.append(value)
            if pos_strings is not None:
                pos_strings.append(formatted)

    if pos_strings is not None:
        return x_coords, y_coords, pos_strings
    return x_coords, y_coords


def get_repeats_aggregated_data(cur, contig_id, xstart=None, xend=None):
    """Get aggregated repeat data as two sets of curve feature dicts.

    For each repeat type (direct, inverted):
      - Queries raw repeats from Contig_directRepeats / Contig_invertedRepeats
      - Aggregates into count and max-identity curves via sweep-line
      - Returns curve feature dicts ready for make_bokeh_subplot()

    Args:
        cur: DuckDB cursor
        contig_id: Contig ID
        xstart: Optional start position for filtering
        xend: Optional end position for filtering

    Returns:
        (count_feature_dicts, identity_feature_dicts) - each is a list of 0-2
        feature dicts (one per repeat type that has data)
    """
    # Each entry: (count_variable, identity_variable, table_name)
    repeat_types = [
        ("direct_repeat_count", "direct_repeat_identity", "Contig_directRepeats"),
        ("inverted_repeat_count", "inverted_repeat_identity", "Contig_invertedRepeats"),
    ]

    count_feature_dicts = []
    identity_feature_dicts = []

    for count_var, identity_var, table_name in repeat_types:
        # Get variable info for count and identity from Variable table
        cur.execute(
            "SELECT Color, Alpha, Fill_alpha, Title "
            "FROM Variable WHERE Variable_name=?",
            (count_var,)
        )
        count_meta = cur.fetchone()

        cur.execute(
            "SELECT Color, Alpha, Fill_alpha, Title "
            "FROM Variable WHERE Variable_name=?",
            (identity_var,)
        )
        identity_meta = cur.fetchone()

        if not count_meta and not identity_meta:
            continue

        # Query raw repeats
        try:
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
            continue

        if not rows:
            continue

        # Capture original (unclipped) repeat positions before clipping
        original_positions = [(min(p1, p2), max(p1, p2)) for p1, p2, _, _, _ in rows]

        # Clip to requested range
        if xstart is not None and xend is not None:
            clipped = []
            for pos1, pos2, pos1p, pos2p, pident_int in rows:
                clipped.append((max(pos1, xstart), min(pos2, xend), pos1p, pos2p, pident_int))
            rows = clipped

        # Aggregate
        count_segments, identity_segments = aggregate_repeats_sweep(rows, original_positions=original_positions)

        # Convert to curve coordinates (with position formatters for tooltips)
        count_x, count_y, count_pos = _segments_to_curve_coords(count_segments, position_formatter=_format_count_positions)
        ident_x, ident_y, ident_pos = _segments_to_curve_coords(identity_segments, position_formatter=_format_identity_position)

        # Build feature dicts (type="curve" overrides whatever Variable table says)
        if count_x and count_meta:
            color, alpha, fill_alpha, title = count_meta
            count_feature_dicts.append({
                "type": "curve",
                "color": color,
                "alpha": alpha,
                "fill_alpha": fill_alpha,
                "size": 1,
                "title": title,
                "x": count_x,
                "y": count_y,
                "has_stats": False,
                "repeat_positions": count_pos,
            })
        if ident_x and identity_meta:
            color, alpha, fill_alpha, title = identity_meta
            identity_feature_dicts.append({
                "type": "curve",
                "color": color,
                "alpha": alpha,
                "fill_alpha": fill_alpha,
                "size": 1,
                "title": title,
                "x": ident_x,
                "y": ident_y,
                "has_stats": False,
                "repeat_positions": ident_pos,
            })

    return count_feature_dicts, identity_feature_dicts


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


### Downsampling threshold: windows larger than this use SQL-side binning
_DOWNSAMPLE_THRESHOLD = 100_000  # 100 kb
_NUM_BINS = 1000

# Features stored as relative values (value / coverage * 1000) — need ÷1000 scaling
RELATIVE_SCALED_FEATURES = {
    "Feature_left_clippings", "Feature_right_clippings", "Feature_insertions", "Feature_deletions",
    "Feature_mismatches", "Feature_reads_starts", "Feature_reads_ends", "Feature_non_inward_pairs",
    "Feature_mate_not_mapped", "Feature_mate_on_another_contig"
}

# Features stored as INTEGER × 100 — need ÷100 scaling
_SCALED_FEATURES = ["Feature_mapq", "Contig_GCSkew"]

# Features that have statistics columns (min/mean/max per position)
_FEATURES_WITH_STATS = ["left_clippings", "right_clippings", "insertions", "reads_starts", "reads_ends"]

# Features that have Sequence and Sequence_prevalence columns
_FEATURES_WITH_SEQUENCES = ["mismatches", "insertions", "left_clippings", "right_clippings", "reads_starts", "reads_ends"]

# SQL expression for Last_position that handles NULL.
# Single-position bar features (insertions, mismatches, clippings, reads_starts/ends) store
# NULL for Last_position to save space; COALESCE recovers First_position at read time.
_LP = "COALESCE(Last_position, First_position)"


def _rle_weighted_bin_sql(feature_table, is_contig_table, xstart, xend, sample_id, contig_id, cur, num_bins=_NUM_BINS, value_col="Value", has_sequences=False, type_picked="curve"):
    """SQL MAX binning: maximum Value per bin to preserve spikes.

    Each RLE row is expanded to all bins it overlaps via generate_series,
    ensuring long RLE runs contribute to every bin they span.

    For bar features (type_picked="bars"), the full run (First_position, Last_position) of the
    max-value row is preserved. If the same run is the max in multiple adjacent bins, it is
    emitted only once (deduplicated). No zero-fill is performed for bars.

    For curve features, behavior is unchanged: position is clamped to bin boundaries and gaps
    between occupied bins are zero-filled to prevent varea interpolation artifacts.

    Args:
        value_col: Column name to bin (default: "Value")
        has_sequences: If True, also fetch Sequence and Sequence_prevalence of the max-value row.
        type_picked: "bars" or "curve" — controls post-processing (dedup/zero-fill).

    Returns for curves:
        (x_coords, y_coords) if has_sequences=False,
        (x_coords, y_coords, seq_coords, prev_coords) if has_sequences=True.
    Returns for bars:
        (x_coords, y_coords, width_coords, first_pos_coords, last_pos_coords) if has_sequences=False,
        (x_coords, y_coords, width_coords, first_pos_coords, last_pos_coords, seq_coords, prev_coords) if has_sequences=True.
    """
    is_bars = type_picked == "bars"
    window_size = xend - xstart + 1
    bin_width = max(1, window_size // num_bins)

    # Each RLE row is expanded to all bins it overlaps via generate_series.
    # ARG_MAX returns attributes of the max-value row per bin.
    seq_inner = ", f.Sequence AS seq, f.Sequence_prevalence AS seq_prev" if has_sequences else ""
    seq_outer = (
        f", ARG_MAX(seq, CAST(val AS DOUBLE)) AS max_seq"
        f", ARG_MAX(seq_prev, CAST(val AS DOUBLE)) AS max_seq_prev"
    ) if has_sequences else ""

    # For bars: retrieve full run endpoints; mid_pos not needed (fp/lp give exact position).
    # For curves: mid_pos (clamped to bin) prevents non-monotonic x from long RLE runs.
    bars_inner = f", f.First_position AS fp, {_LP} AS lp" if is_bars else ""
    bars_outer = (
        f", ARG_MAX(fp, CAST(val AS DOUBLE)) AS max_fp"
        f", ARG_MAX(lp, CAST(val AS DOUBLE)) AS max_lp"
    ) if is_bars else ""
    curve_mid_inner = "" if is_bars else f", (f.First_position + {_LP}) / 2 AS mid_pos"
    curve_mid_outer = "" if is_bars else f", ARG_MAX(mid_pos, CAST(val AS DOUBLE)) AS max_pos"

    if is_contig_table:
        sql = (
            f"SELECT bin_idx AS bin, "
            f"MAX(CAST(val AS DOUBLE)) AS max_value"
            f"{curve_mid_outer}{bars_outer}{seq_outer} "
            f"FROM ("
            f"  SELECT f.{value_col} AS val{curve_mid_inner}{bars_inner}{seq_inner}, "
            f"  UNNEST(generate_series("
            f"    GREATEST(0, CAST((GREATEST(f.First_position, ?) - ?) / ? AS INTEGER)), "
            f"    LEAST(? - 1, CAST((LEAST({_LP}, ?) - ?) / ? AS INTEGER))"
            f"  )) AS bin_idx "
            f"  FROM {feature_table} f "
            f"  WHERE f.Contig_id = ? AND {_LP} >= ? AND f.First_position <= ?"
            f") sub "
            f"GROUP BY bin_idx ORDER BY bin_idx"
        )
        params = (xstart, xstart, bin_width,
                  num_bins, xend, xstart, bin_width,
                  contig_id, xstart, xend)
    else:
        sql = (
            f"SELECT bin_idx AS bin, "
            f"MAX(CAST(val AS DOUBLE)) AS max_value"
            f"{curve_mid_outer}{bars_outer}{seq_outer} "
            f"FROM ("
            f"  SELECT f.{value_col} AS val{curve_mid_inner}{bars_inner}{seq_inner}, "
            f"  UNNEST(generate_series("
            f"    GREATEST(0, CAST((GREATEST(f.First_position, ?) - ?) / ? AS INTEGER)), "
            f"    LEAST(? - 1, CAST((LEAST({_LP}, ?) - ?) / ? AS INTEGER))"
            f"  )) AS bin_idx "
            f"  FROM {feature_table} f "
            f"  WHERE f.Sample_id = ? AND f.Contig_id = ? AND {_LP} >= ? AND f.First_position <= ?"
            f") sub "
            f"GROUP BY bin_idx ORDER BY bin_idx"
        )
        params = (xstart, xstart, bin_width,
                  num_bins, xend, xstart, bin_width,
                  sample_id, contig_id, xstart, xend)

    # Optionally save EXPLAIN plan
    if os.environ.get("BIGBAMB_PROFILE"):
        try:
            cur.execute(f"EXPLAIN {sql}", params)
            explain_rows = cur.fetchall()
            explain_file = f"explain_{feature_table}_{xstart}_{xend}.txt"
            with open(explain_file, "w") as f:
                for erow in explain_rows:
                    f.write(str(erow[0]) + "\n")
        except Exception as e:
            print(f"[get_feature_data] Could not save EXPLAIN for {feature_table}: {e}", flush=True)

    cur.execute(sql, params)
    binned_rows = cur.fetchall()

    if is_bars:
        # --- BARS: keep full run, dedup across bins ---
        # Row layout: (bin, max_value, max_fp, max_lp[, max_seq, max_seq_prev])
        # Collect unique runs: same (fp, lp) appearing as max in multiple bins → emit once
        seen = set()
        x_coords = []
        y_coords = []
        width_coords = []
        first_pos_coords = []
        last_pos_coords = []
        seq_coords = [] if has_sequences else None
        prev_coords = [] if has_sequences else None
        for row in binned_rows:
            if has_sequences:
                bin_idx, max_value, max_fp, max_lp, max_seq, max_seq_prev = row
            else:
                bin_idx, max_value, max_fp, max_lp = row
            if max_value is None or bin_idx is None or max_fp is None or max_lp is None:
                continue
            fp = int(max_fp)
            lp = int(max_lp)
            key = (fp, lp)
            if key in seen:
                continue
            seen.add(key)
            midpoint = (fp + lp) / 2.0
            width = lp - fp + 1
            x_coords.append(midpoint)
            y_coords.append(float(max_value))
            width_coords.append(width)
            first_pos_coords.append(fp)
            last_pos_coords.append(lp)
            if has_sequences:
                seq_coords.append(max_seq if max_seq is not None else "")
                prev_coords.append(max_seq_prev / 10.0 if max_seq_prev is not None else None)
        if not x_coords:
            if has_sequences:
                return [], [], [], [], [], [], []
            return [], [], [], [], []
        if has_sequences:
            return x_coords, y_coords, width_coords, first_pos_coords, last_pos_coords, seq_coords, prev_coords
        return x_coords, y_coords, width_coords, first_pos_coords, last_pos_coords
    else:
        # --- CURVES: clamped position + zero-fill (unchanged) ---
        # Build lookup of bin_idx -> (position, max_value[, seq, seq_prev])
        bin_data = {}
        for row in binned_rows:
            if has_sequences:
                bin_idx, max_value, max_pos, max_seq, max_seq_prev = row
            else:
                bin_idx, max_value, max_pos = row
            if max_value is not None and bin_idx is not None:
                bi = int(bin_idx)
                bin_start = xstart + bi * bin_width
                position = max(bin_start, min(bin_start + bin_width - 1, int(max_pos))) if max_pos is not None else int(bin_start + bin_width // 2)
                if has_sequences:
                    bin_data[bi] = (position, float(max_value), max_seq, max_seq_prev)
                else:
                    bin_data[bi] = (position, float(max_value))

        if not bin_data:
            if has_sequences:
                return [], [], [], []
            return [], []

        # Zero-fill all bins between min and max occupied bins
        # This prevents linear interpolation artifacts (triangles) in varea/line rendering
        # For zero-filled bins, use the bin center as position
        min_bin = min(bin_data)
        max_bin = max(bin_data)
        x_coords = []
        y_coords = []
        seq_coords = [] if has_sequences else None
        prev_coords = [] if has_sequences else None
        for bi in range(min_bin, max_bin + 1):
            if bi in bin_data:
                if has_sequences:
                    position, value, seq, seq_prev = bin_data[bi]
                    seq_coords.append(seq if seq is not None else "")
                    prev_coords.append(seq_prev / 10.0 if seq_prev is not None else None)
                else:
                    position, value = bin_data[bi]
                x_coords.append(position)
                y_coords.append(value)
            else:
                # Zero-filled bin: use bin center (bin_width matches SQL integer division)
                bin_center = int(xstart + (bi + 0.5) * bin_width)
                x_coords.append(bin_center)
                y_coords.append(0.0)
                if has_sequences:
                    seq_coords.append("")
                    prev_coords.append(None)
        if has_sequences:
            return x_coords, y_coords, seq_coords, prev_coords
        return x_coords, y_coords


def _rle_weighted_bin_primary_reads(xstart, xend, sample_id, contig_id, cur, num_bins=_NUM_BINS):
    """Bin primary_reads by binning plus and minus strands separately, then summing per bin.
    
    Each RLE row is expanded to all bins it overlaps via generate_series.
    Returns (x_coords, y_coords) where x is the bin center.
    """
    window_size = xend - xstart + 1
    bin_width = max(1, window_size // num_bins)

    # Store per-bin: {bin_idx: [(position, value), ...]} for each strand
    bin_contributions = {}

    for table_name in ["Feature_primary_reads_plus_only", "Feature_primary_reads_minus_only"]:
        sql = (
            f"SELECT bin_idx AS bin, "
            f"MAX(CAST(Value AS DOUBLE)) AS max_value, "
            f"ARG_MAX(mid_pos, CAST(Value AS DOUBLE)) AS max_pos "
            f"FROM ("
            f"  SELECT f.Value, (f.First_position + {_LP}) / 2 AS mid_pos, "
            f"  UNNEST(generate_series("
            f"    GREATEST(0, CAST((GREATEST(f.First_position, ?) - ?) / ? AS INTEGER)), "
            f"    LEAST(? - 1, CAST((LEAST({_LP}, ?) - ?) / ? AS INTEGER))"
            f"  )) AS bin_idx "
            f"  FROM {table_name} f "
            f"  WHERE f.Sample_id = ? AND f.Contig_id = ? AND {_LP} >= ? AND f.First_position <= ?"
            f") sub "
            f"GROUP BY bin_idx ORDER BY bin_idx"
        )
        params = (xstart, xstart, bin_width,
                  num_bins, xend, xstart, bin_width,
                  sample_id, contig_id, xstart, xend)
        cur.execute(sql, params)
        for bin_idx, max_value, max_pos in cur.fetchall():
            if max_value is not None and bin_idx is not None:
                bi = int(bin_idx)
                bin_start = xstart + bi * bin_width
                position = max(bin_start, min(bin_start + bin_width - 1, int(max_pos))) if max_pos is not None else int(bin_start + bin_width // 2)
                if bi not in bin_contributions:
                    bin_contributions[bi] = []
                bin_contributions[bi].append((position, float(max_value)))

    if not bin_contributions:
        return [], []

    # Sum contributions per bin and track position of max contributor
    results = {}  # bin_idx -> (sum_value, position_of_max_contributor)
    for bi, contributions in bin_contributions.items():
        total_value = sum(val for pos, val in contributions)
        # Use position from the strand that contributed more
        max_contributor = max(contributions, key=lambda x: x[1])
        results[bi] = (total_value, max_contributor[0])

    # Zero-fill all bins between min and max occupied bins
    min_bin = min(results)
    max_bin = max(results)
    x_coords = []
    y_coords = []
    for bi in range(min_bin, max_bin + 1):
        if bi in results:
            value, position = results[bi]
            x_coords.append(position)
            y_coords.append(value)
        else:
            # Zero-filled bin: use bin center (bin_width matches SQL integer division)
            bin_center = int(xstart + (bi + 0.5) * bin_width)
            x_coords.append(bin_center)
            y_coords.append(0.0)
    return x_coords, y_coords


### Function to get features of one variable
def get_feature_data(cur, feature, contig_id, sample_id, xstart=None, xend=None, variable_metadata=None, downsample_threshold=None, min_relative_value=0.0):
    """Get feature data for plotting.

    Uses SQL-side binning (midpoint assignment, 1000 bins) for windows larger than the
    downsampling threshold (default: _DOWNSAMPLE_THRESHOLD = 100kb). For small windows or
    undefined ranges, full-resolution RLE expansion is used.

    Args:
        cur: DuckDB cursor
        feature: Feature name to query
        contig_id: Contig ID
        sample_id: Sample ID
        xstart: Optional start position for filtering (only fetch data intersecting this range)
        xend: Optional end position for filtering (only fetch data intersecting this range)
        variable_metadata: Optional cached result from get_variable_metadata(); avoids re-querying Variable table
        downsample_threshold: Override the module-level _DOWNSAMPLE_THRESHOLD for this call
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
            "y": [],
            "is_relative_scaled": False,  # Will be set to True for relative-scaled features
        }

        # Check if this feature stores scaled values (stored as INTEGER ×100)
        is_scaled = feature_table in _SCALED_FEATURES

        # Check if this feature stores relative values (stored as INTEGER ×1000)
        is_relative_scaled = feature_table in RELATIVE_SCALED_FEATURES

        # Detect contig-level table (no Sample_id column)
        is_contig_table = feature_table.startswith("Contig_")

        # --- DOWNSAMPLING PATH: large windows (> threshold) ---
        _threshold = downsample_threshold if downsample_threshold is not None else _DOWNSAMPLE_THRESHOLD
        use_binning = (
            xstart is not None and xend is not None
            and (xend - xstart) > _threshold
        )

        if use_binning:
            # Check if this feature has sequence columns (needed for binning path too)
            has_sequences_bin = feature_table in [f"Feature_{f}" for f in _FEATURES_WITH_SEQUENCES]
            is_bars = type_picked == "bars"

            # Special case: primary_reads combines two strand tables (always curve)
            if feature_table == "Feature_primary_reads":
                x_coords, y_coords = _rle_weighted_bin_primary_reads(
                    xstart, xend, sample_id, contig_id, cur
                )
                seq_coords_bin = []
                prev_coords_bin = []
                width_coords_bin = []
                fps_bin = []
                lps_bin = []
            elif is_bars:
                # Bars: get full run with dedup
                if has_sequences_bin:
                    x_coords, y_coords, width_coords_bin, fps_bin, lps_bin, seq_coords_bin, prev_coords_bin = _rle_weighted_bin_sql(
                        feature_table, is_contig_table, xstart, xend,
                        sample_id, contig_id, cur, has_sequences=True, type_picked="bars"
                    )
                else:
                    x_coords, y_coords, width_coords_bin, fps_bin, lps_bin = _rle_weighted_bin_sql(
                        feature_table, is_contig_table, xstart, xend,
                        sample_id, contig_id, cur, type_picked="bars"
                    )
                    seq_coords_bin = []
                    prev_coords_bin = []
            elif has_sequences_bin:
                x_coords, y_coords, seq_coords_bin, prev_coords_bin = _rle_weighted_bin_sql(
                    feature_table, is_contig_table, xstart, xend,
                    sample_id, contig_id, cur, has_sequences=True
                )
                width_coords_bin = []
                fps_bin = []
                lps_bin = []
            else:
                x_coords, y_coords = _rle_weighted_bin_sql(
                    feature_table, is_contig_table, xstart, xend,
                    sample_id, contig_id, cur
                )
                seq_coords_bin = []
                prev_coords_bin = []
                width_coords_bin = []
                fps_bin = []
                lps_bin = []

            # Apply scaling after binning (MAX was on raw INTEGER values)
            if is_relative_scaled and y_coords:
                # Relative features stored as INTEGER × 1000
                y_coords = [v / 1000.0 for v in y_coords]
            elif is_scaled and y_coords:
                # Scaled features stored as INTEGER × 100
                y_coords = [v / 100.0 for v in y_coords]

            # Filter by min_relative_value — keep all parallel arrays aligned
            if is_relative_scaled and min_relative_value > 0.0 and x_coords:
                filtered = [
                    i for i, y in enumerate(y_coords) if y is None or y >= min_relative_value
                ]
                x_coords = [x_coords[i] for i in filtered]
                y_coords = [y_coords[i] for i in filtered]
                if seq_coords_bin:
                    seq_coords_bin = [seq_coords_bin[i] for i in filtered]
                if prev_coords_bin:
                    prev_coords_bin = [prev_coords_bin[i] for i in filtered]
                if width_coords_bin:
                    width_coords_bin = [width_coords_bin[i] for i in filtered]
                    fps_bin = [fps_bin[i] for i in filtered]
                    lps_bin = [lps_bin[i] for i in filtered]

            feature_dict["x"] = x_coords
            feature_dict["y"] = y_coords
            feature_dict["has_stats"] = False  # stats don't aggregate meaningfully
            feature_dict["has_sequences"] = has_sequences_bin
            feature_dict["is_relative_scaled"] = is_relative_scaled
            if is_bars and width_coords_bin:
                feature_dict["width"] = width_coords_bin
                feature_dict["first_pos"] = fps_bin
                feature_dict["last_pos"] = lps_bin
            if has_sequences_bin and x_coords:
                feature_dict["sequence"] = seq_coords_bin
                feature_dict["sequence_prevalence"] = prev_coords_bin
            if x_coords:
                list_feature_dict.append(feature_dict)

        else:
            # --- FULL RESOLUTION PATH: small windows or undefined range ---
            # Check if this feature has statistics columns
            has_stats = feature_table in [f"Feature_{f}" for f in _FEATURES_WITH_STATS]
            has_sequences = feature_table in [f"Feature_{f}" for f in _FEATURES_WITH_SEQUENCES]

            # Special handling for primary_reads: combine strand tables in Python
            if feature_table == "Feature_primary_reads":
                position_filter = ""
                params = [sample_id, contig_id]
                if xstart is not None and xend is not None:
                    position_filter = f" AND {_LP} >= ? AND First_position <= ?"
                    params.extend([xstart, xend])
                cur.execute(
                    f"SELECT First_position, {_LP}, Value FROM Feature_primary_reads_plus_only "
                    f"WHERE Sample_id=? AND Contig_id=?{position_filter} ORDER BY First_position",
                    tuple(params)
                )
                plus_rows = cur.fetchall()
                cur.execute(
                    f"SELECT First_position, {_LP}, Value FROM Feature_primary_reads_minus_only "
                    f"WHERE Sample_id=? AND Contig_id=?{position_filter} ORDER BY First_position",
                    tuple(params)
                )
                minus_rows = cur.fetchall()
                data_rows = merge_rle_segments(plus_rows, minus_rows)
            elif has_stats and has_sequences:
                position_filter = ""
                params = [sample_id, contig_id]
                if xstart is not None and xend is not None:
                    position_filter = f" AND {_LP} >= ? AND First_position <= ?"
                    params.extend([xstart, xend])
                cur.execute(
                    f"SELECT First_position, {_LP}, Value, Mean, Median, Std, Sequence, Sequence_prevalence FROM {feature_table} "
                    f"WHERE Sample_id=? AND Contig_id=?{position_filter} ORDER BY First_position",
                    tuple(params)
                )
                data_rows = cur.fetchall()
            elif has_stats:
                position_filter = ""
                params = [sample_id, contig_id]
                if xstart is not None and xend is not None:
                    position_filter = f" AND {_LP} >= ? AND First_position <= ?"
                    params.extend([xstart, xend])
                cur.execute(
                    f"SELECT First_position, {_LP}, Value, Mean, Median, Std FROM {feature_table} "
                    f"WHERE Sample_id=? AND Contig_id=?{position_filter} ORDER BY First_position",
                    tuple(params)
                )
                data_rows = cur.fetchall()
            elif has_sequences:
                position_filter = ""
                params = [sample_id, contig_id]
                if xstart is not None and xend is not None:
                    position_filter = f" AND {_LP} >= ? AND First_position <= ?"
                    params.extend([xstart, xend])
                cur.execute(
                    f"SELECT First_position, {_LP}, Value, Sequence, Sequence_prevalence FROM {feature_table} "
                    f"WHERE Sample_id=? AND Contig_id=?{position_filter} ORDER BY First_position",
                    tuple(params)
                )
                data_rows = cur.fetchall()
            else:
                position_filter = ""
                if is_contig_table:
                    params = [contig_id]
                else:
                    params = [sample_id, contig_id]
                if xstart is not None and xend is not None:
                    position_filter = f" AND {_LP} >= ? AND First_position <= ?"
                    params.extend([xstart, xend])
                if is_contig_table:
                    cur.execute(
                        f"SELECT First_position, {_LP}, Value FROM {feature_table} "
                        f"WHERE Contig_id=?{position_filter} ORDER BY First_position",
                        tuple(params)
                    )
                else:
                    cur.execute(
                        f"SELECT First_position, {_LP}, Value FROM {feature_table} "
                        f"WHERE Sample_id=? AND Contig_id=?{position_filter} ORDER BY First_position",
                        tuple(params)
                    )
                data_rows = cur.fetchall()

            # Clip feature positions to requested range
            if xstart is not None and xend is not None:
                clipped_rows = []
                for row in data_rows:
                    if has_stats and has_sequences:
                        first_pos, last_pos, value, mean, median, std, seq, seq_prev = row
                        clipped_first = max(first_pos, xstart)
                        clipped_last = min(last_pos, xend)
                        if clipped_first <= clipped_last:
                            clipped_rows.append((clipped_first, clipped_last, value, mean, median, std, seq, seq_prev))
                    elif has_stats:
                        first_pos, last_pos, value, mean, median, std = row
                        clipped_first = max(first_pos, xstart)
                        clipped_last = min(last_pos, xend)
                        if clipped_first <= clipped_last:
                            clipped_rows.append((clipped_first, clipped_last, value, mean, median, std))
                    elif has_sequences:
                        first_pos, last_pos, value, seq, seq_prev = row
                        clipped_first = max(first_pos, xstart)
                        clipped_last = min(last_pos, xend)
                        if clipped_first <= clipped_last:
                            clipped_rows.append((clipped_first, clipped_last, value, seq, seq_prev))
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
            sequence_coords = []
            prevalence_coords = []
            width_coords = []
            first_pos_coords = []
            last_pos_coords = []
            for row in data_rows:
                if has_stats and has_sequences:
                    first_pos, last_pos, value, mean, median, std, seq, seq_prev = row
                elif has_stats:
                    first_pos, last_pos, value, mean, median, std = row
                    seq = seq_prev = None
                elif has_sequences:
                    first_pos, last_pos, value, seq, seq_prev = row
                    mean = median = std = None
                else:
                    first_pos, last_pos, value = row
                    mean = median = std = seq = seq_prev = None
                # Convert sequence prevalence from ×10 integer to percentage
                if seq_prev is not None:
                    seq_prev = seq_prev / 10.0
                if seq is None:
                    seq = ""
                # Scale values based on storage format
                if is_relative_scaled:
                    # Relative-to-coverage features stored as INTEGER × 1000
                    value = value / 1000.0 if value is not None else None
                elif is_scaled:
                    # Scaled features (mapq, gc_skew) stored as INTEGER × 100
                    value = value / 100.0 if value is not None else None
                if is_relative_scaled and min_relative_value > 0.0 and value is not None and value < min_relative_value:
                    continue
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
                    if has_sequences:
                        sequence_coords.append(seq)
                        prevalence_coords.append(seq_prev)
                else:
                    if first_pos == last_pos:
                        x_coords.append(first_pos)
                        y_coords.append(value)
                        if has_stats:
                            mean_coords.append(mean)
                            median_coords.append(median)
                            std_coords.append(std)
                        if has_sequences:
                            sequence_coords.append(seq)
                            prevalence_coords.append(seq_prev)
                    else:
                        x_coords.extend([first_pos, last_pos])
                        y_coords.extend([value, value])
                        if has_stats:
                            mean_coords.extend([mean, mean])
                            median_coords.extend([median, median])
                            std_coords.extend([std, std])
                        if has_sequences:
                            sequence_coords.extend([seq, seq])
                            prevalence_coords.extend([seq_prev, seq_prev])
            feature_dict["x"] = x_coords
            feature_dict["y"] = y_coords
            feature_dict["has_stats"] = has_stats
            feature_dict["has_sequences"] = has_sequences
            feature_dict["is_relative_scaled"] = is_relative_scaled
            if type_picked == "bars":
                feature_dict["width"] = width_coords
                feature_dict["first_pos"] = first_pos_coords
                feature_dict["last_pos"] = last_pos_coords
            if has_stats:
                feature_dict["mean"] = mean_coords
                feature_dict["median"] = median_coords
                feature_dict["std"] = std_coords
            if has_sequences:
                feature_dict["sequence"] = sequence_coords
                feature_dict["sequence_prevalence"] = prevalence_coords
            if x_coords:
                list_feature_dict.append(feature_dict)

    return list_feature_dict


### Function to get features for multiple samples in a single batch
def _expand_rle_rows(data_rows, type_picked, has_stats, is_scaled, xstart, xend, is_relative_scaled=False, min_relative_value=0.0, has_sequences=False):
    """Expand RLE rows into plot coordinates (shared logic for single and batch).

    Args:
        data_rows: List of tuples (First_position, Last_position, Value[, Mean, Median, Std][, Sequence, Sequence_prevalence])
        type_picked: Plot type ('bars' or 'curve')
        has_stats: Whether rows include statistics columns
        is_scaled: Whether values need to be divided by 100
        xstart: Optional start position for clipping
        xend: Optional end position for clipping
        is_relative_scaled: Whether values are relative-to-coverage (need to be divided by 1000)
        has_sequences: Whether rows include Sequence and Sequence_prevalence columns

    Returns:
        Dict with x, y, and optional width/stats/sequence coordinate lists, or None if no data
    """
    # Clip feature positions to requested range
    if xstart is not None and xend is not None:
        clipped_rows = []
        for row in data_rows:
            # Determine the number of base columns (first, last, value)
            # then stats (mean, median, std), then sequences (seq, seq_prev)
            if has_stats and has_sequences:
                first_pos, last_pos, value, mean, median, std, seq, seq_prev = row
                clipped_first = max(first_pos, xstart)
                clipped_last = min(last_pos, xend)
                if clipped_first <= clipped_last:
                    clipped_rows.append((clipped_first, clipped_last, value, mean, median, std, seq, seq_prev))
            elif has_stats:
                first_pos, last_pos, value, mean, median, std = row
                clipped_first = max(first_pos, xstart)
                clipped_last = min(last_pos, xend)
                if clipped_first <= clipped_last:
                    clipped_rows.append((clipped_first, clipped_last, value, mean, median, std))
            elif has_sequences:
                first_pos, last_pos, value, seq, seq_prev = row
                clipped_first = max(first_pos, xstart)
                clipped_last = min(last_pos, xend)
                if clipped_first <= clipped_last:
                    clipped_rows.append((clipped_first, clipped_last, value, seq, seq_prev))
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
    sequence_coords = []
    prevalence_coords = []
    width_coords = []
    first_pos_coords = []
    last_pos_coords = []

    for row in data_rows:
        if has_stats and has_sequences:
            first_pos, last_pos, value, mean, median, std, seq, seq_prev = row
        elif has_stats:
            first_pos, last_pos, value, mean, median, std = row
            seq = seq_prev = None
        elif has_sequences:
            first_pos, last_pos, value, seq, seq_prev = row
            mean = median = std = None
        else:
            first_pos, last_pos, value = row
            mean = median = std = seq = seq_prev = None

        # Convert sequence prevalence from ×10 integer to percentage
        if seq_prev is not None:
            seq_prev = seq_prev / 10.0
        if seq is None:
            seq = ""

        # Scale values based on storage format
        if is_relative_scaled:
            # Relative-to-coverage features stored as INTEGER × 1000
            value = value / 1000.0 if value is not None else None
        elif is_scaled:
            # Scaled features (mapq, gc_skew) stored as INTEGER × 100
            value = value / 100.0 if value is not None else None

        if is_relative_scaled and min_relative_value > 0.0 and value is not None and value < min_relative_value:
            continue

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
            if has_sequences:
                sequence_coords.append(seq)
                prevalence_coords.append(seq_prev)
        else:
            if first_pos == last_pos:
                x_coords.append(first_pos)
                y_coords.append(value)
                if has_stats:
                    mean_coords.append(mean)
                    median_coords.append(median)
                    std_coords.append(std)
                if has_sequences:
                    sequence_coords.append(seq)
                    prevalence_coords.append(seq_prev)
            else:
                x_coords.extend([first_pos, last_pos])
                y_coords.extend([value, value])
                if has_stats:
                    mean_coords.extend([mean, mean])
                    median_coords.extend([median, median])
                    std_coords.extend([std, std])
                if has_sequences:
                    sequence_coords.extend([seq, seq])
                    prevalence_coords.extend([seq_prev, seq_prev])

    if not x_coords:
        return None

    result = {"x": x_coords, "y": y_coords, "has_stats": has_stats, "has_sequences": has_sequences}
    if type_picked == "bars":
        result["width"] = width_coords
        result["first_pos"] = first_pos_coords
        result["last_pos"] = last_pos_coords
    if has_stats:
        result["mean"] = mean_coords
        result["median"] = median_coords
        result["std"] = std_coords
    if has_sequences:
        result["sequence"] = sequence_coords
        result["sequence_prevalence"] = prevalence_coords
    return result


def _rle_weighted_bin_batch_sql(feature_table, xstart, xend, sample_ids, contig_id, cur, num_bins=_NUM_BINS, value_col="Value", has_sequences=False, type_picked="curve"):
    """Run SQL-side binning for multiple samples at once using MAX aggregation to preserve spikes.
    
    Each RLE row is expanded to all bins it overlaps via generate_series.

    For bar features (type_picked="bars"), the full run (First_position, Last_position) of the
    max-value row is preserved and deduplicated across bins. No zero-fill for bars.

    Args:
        value_col: Column name to bin (default: "Value")
        has_sequences: If True, also fetch Sequence and Sequence_prevalence of the max-value row.
        type_picked: "bars" or "curve" — controls post-processing.
    
    Returns dict per sample:
        curves: {sid: (x, y)} or {sid: (x, y, seq, prev)}
        bars:   {sid: (x, y, width, first_pos, last_pos)} or {sid: (x, y, width, first_pos, last_pos, seq, prev)}
    """
    is_bars = type_picked == "bars"
    window_size = xend - xstart + 1
    bin_width = max(1, window_size // num_bins)
    placeholders = ", ".join(["?"] * len(sample_ids))

    seq_inner = ", f.Sequence AS seq, f.Sequence_prevalence AS seq_prev" if has_sequences else ""
    seq_outer = (
        f", ARG_MAX(seq, CAST(val AS DOUBLE)) AS max_seq"
        f", ARG_MAX(seq_prev, CAST(val AS DOUBLE)) AS max_seq_prev"
    ) if has_sequences else ""

    # For bars: retrieve full run endpoints; mid_pos not needed (fp/lp give exact position).
    # For curves: mid_pos (clamped to bin) prevents non-monotonic x from long RLE runs.
    bars_inner = f", f.First_position AS fp, {_LP} AS lp" if is_bars else ""
    bars_outer = (
        f", ARG_MAX(fp, CAST(val AS DOUBLE)) AS max_fp"
        f", ARG_MAX(lp, CAST(val AS DOUBLE)) AS max_lp"
    ) if is_bars else ""
    curve_mid_inner = "" if is_bars else f", (f.First_position + {_LP}) / 2 AS mid_pos"
    curve_mid_outer = "" if is_bars else f", ARG_MAX(mid_pos, CAST(val AS DOUBLE)) AS max_pos"

    sql = (
        f"SELECT Sample_id, bin_idx AS bin, "
        f"MAX(CAST(val AS DOUBLE)) AS max_value"
        f"{curve_mid_outer}{bars_outer}{seq_outer} "
        f"FROM ("
        f"  SELECT f.Sample_id, f.{value_col} AS val{curve_mid_inner}{bars_inner}{seq_inner}, "
        f"  UNNEST(generate_series("
        f"    GREATEST(0, CAST((GREATEST(f.First_position, ?) - ?) / ? AS INTEGER)), "
        f"    LEAST(? - 1, CAST((LEAST({_LP}, ?) - ?) / ? AS INTEGER))"
        f"  )) AS bin_idx "
        f"  FROM {feature_table} f "
        f"  WHERE f.Sample_id IN ({placeholders}) AND f.Contig_id = ? AND {_LP} >= ? AND f.First_position <= ?"
        f") sub "
        f"GROUP BY Sample_id, bin_idx ORDER BY Sample_id, bin_idx"
    )
    params = (xstart, xstart, bin_width,
              num_bins, xend, xstart, bin_width) + tuple(sample_ids) + (contig_id, xstart, xend)
    cur.execute(sql, params)

    if is_bars:
        # --- BARS: dedup runs per sample ---
        # Row layout: (sid, bin, max_value, max_fp, max_lp[, max_seq, max_seq_prev])
        raw_by_sample = {sid: [] for sid in sample_ids}
        for row in cur.fetchall():
            if has_sequences:
                sid, bin_idx, max_value, max_fp, max_lp, max_seq, max_seq_prev = row
            else:
                sid, bin_idx, max_value, max_fp, max_lp = row
            if sid in raw_by_sample and max_value is not None and max_fp is not None and max_lp is not None:
                if has_sequences:
                    raw_by_sample[sid].append((int(max_fp), int(max_lp), float(max_value), max_seq, max_seq_prev))
                else:
                    raw_by_sample[sid].append((int(max_fp), int(max_lp), float(max_value)))

        result = {}
        for sid in sample_ids:
            seen = set()
            x_coords = []
            y_coords = []
            width_coords = []
            first_pos_coords = []
            last_pos_coords = []
            seq_coords = [] if has_sequences else None
            prev_coords = [] if has_sequences else None
            for entry in raw_by_sample[sid]:
                if has_sequences:
                    fp, lp, val, seq, seq_prev = entry
                else:
                    fp, lp, val = entry
                key = (fp, lp)
                if key in seen:
                    continue
                seen.add(key)
                midpoint = (fp + lp) / 2.0
                width = lp - fp + 1
                x_coords.append(midpoint)
                y_coords.append(val)
                width_coords.append(width)
                first_pos_coords.append(fp)
                last_pos_coords.append(lp)
                if has_sequences:
                    seq_coords.append(seq if seq is not None else "")
                    prev_coords.append(seq_prev / 10.0 if seq_prev is not None else None)
            if has_sequences:
                result[sid] = (x_coords, y_coords, width_coords, first_pos_coords, last_pos_coords, seq_coords, prev_coords)
            else:
                result[sid] = (x_coords, y_coords, width_coords, first_pos_coords, last_pos_coords)
        return result
    else:
        # --- CURVES: clamped position + zero-fill (unchanged) ---
        by_sample = {sid: {} for sid in sample_ids}
        for row in cur.fetchall():
            if has_sequences:
                sid, bin_idx, max_value, max_pos, max_seq, max_seq_prev = row
            else:
                sid, bin_idx, max_value, max_pos = row
            if sid in by_sample and max_value is not None and bin_idx is not None:
                bi = int(bin_idx)
                bin_start = xstart + bi * bin_width
                position = max(bin_start, min(bin_start + bin_width - 1, int(max_pos))) if max_pos is not None else int(bin_start + bin_width // 2)
                if has_sequences:
                    by_sample[sid][bi] = (position, float(max_value), max_seq, max_seq_prev)
                else:
                    by_sample[sid][bi] = (position, float(max_value))

        result = {}
        for sid in sample_ids:
            bin_data = by_sample[sid]
            if not bin_data:
                result[sid] = ([], [], [], []) if has_sequences else ([], [])
                continue
            # Zero-fill all bins between min and max occupied bins
            min_bin = min(bin_data)
            max_bin = max(bin_data)
            x_coords = []
            y_coords = []
            seq_coords = [] if has_sequences else None
            prev_coords = [] if has_sequences else None
            for bi in range(min_bin, max_bin + 1):
                if bi in bin_data:
                    if has_sequences:
                        position, value, seq, seq_prev = bin_data[bi]
                        seq_coords.append(seq if seq is not None else "")
                        prev_coords.append(seq_prev / 10.0 if seq_prev is not None else None)
                    else:
                        position, value = bin_data[bi]
                    x_coords.append(position)
                    y_coords.append(value)
                else:
                    # Zero-filled bin: use bin center (bin_width matches SQL integer division)
                    bin_center = int(xstart + (bi + 0.5) * bin_width)
                    x_coords.append(bin_center)
                    y_coords.append(0.0)
                    if has_sequences:
                        seq_coords.append("")
                        prev_coords.append(None)
            if has_sequences:
                result[sid] = (x_coords, y_coords, seq_coords, prev_coords)
            else:
                result[sid] = (x_coords, y_coords)
        return result


def get_feature_data_batch(cur, feature, contig_id, sample_ids, xstart=None, xend=None, variable_metadata=None, downsample_threshold=None, min_relative_value=0.0):
    """Get feature data for multiple samples in a single batch query.

    Uses SQL-side MAX binning for windows larger than the downsampling threshold
    (default: _DOWNSAMPLE_THRESHOLD = 100kb). For small windows, uses full-resolution RLE expansion.

    Args:
        cur: DuckDB cursor
        feature: Subplot name to query
        contig_id: Contig ID
        sample_ids: List of sample IDs to fetch data for
        xstart: Optional start position for filtering
        xend: Optional end position for filtering
        variable_metadata: Optional cached result from get_variable_metadata()
        downsample_threshold: Override the module-level _DOWNSAMPLE_THRESHOLD for this call

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

    # Determine if we should use binning
    _threshold = downsample_threshold if downsample_threshold is not None else _DOWNSAMPLE_THRESHOLD
    use_binning = (
        xstart is not None and xend is not None
        and (xend - xstart) > _threshold
    )

    for row in rows:
        type_picked, color, alpha, fill_alpha, size, title, feature_table = row

        is_relative_scaled = feature_table in RELATIVE_SCALED_FEATURES

        has_stats = feature_table in [f"Feature_{f}" for f in _FEATURES_WITH_STATS]
        has_sequences = feature_table in [f"Feature_{f}" for f in _FEATURES_WITH_SEQUENCES]
        is_scaled = feature_table in _SCALED_FEATURES

        # Detect contig-level table (no Sample_id column)
        is_contig_table = feature_table.startswith("Contig_")

        if use_binning:
            # --- DOWNSAMPLING PATH ---
            if feature_table == "Feature_primary_reads":
                # Bin each strand table separately, then sum per bin per sample
                # Track position of max contributor per bin
                window_size = xend - xstart + 1
                bin_width = max(1, window_size // _NUM_BINS)

                # Store {sid: {bin: [(position, value), ...]}} for each strand
                bins_by_sample = {sid: {} for sid in sample_ids}
                for table_name in ["Feature_primary_reads_plus_only", "Feature_primary_reads_minus_only"]:
                    sql = (
                        f"SELECT Sample_id, bin_idx AS bin, "
                        f"MAX(CAST(Value AS DOUBLE)) AS max_value, "
                        f"ARG_MAX(mid_pos, CAST(Value AS DOUBLE)) AS max_pos "
                        f"FROM ("
                        f"  SELECT f.Sample_id, f.Value, (f.First_position + {_LP}) / 2 AS mid_pos, "
                        f"  UNNEST(generate_series("
                        f"    GREATEST(0, CAST((GREATEST(f.First_position, ?) - ?) / ? AS INTEGER)), "
                        f"    LEAST(? - 1, CAST((LEAST({_LP}, ?) - ?) / ? AS INTEGER))"
                        f"  )) AS bin_idx "
                        f"  FROM {table_name} f "
                        f"  WHERE f.Sample_id IN ({placeholders}) AND f.Contig_id = ? "
                        f"  AND {_LP} >= ? AND f.First_position <= ?"
                        f") sub "
                        f"GROUP BY Sample_id, bin_idx ORDER BY Sample_id, bin_idx"
                    )
                    params = (xstart, xstart, bin_width,
                              num_bins, xend, xstart, bin_width) + tuple(sample_ids) + (contig_id, xstart, xend)
                    cur.execute(sql, params)
                    for sid, bin_idx, max_value, max_pos in cur.fetchall():
                        if sid in bins_by_sample and max_value is not None and bin_idx is not None:
                            bi = int(bin_idx)
                            bin_start = xstart + bi * bin_width
                            position = max(bin_start, min(bin_start + bin_width - 1, int(max_pos))) if max_pos is not None else int(bin_start + bin_width // 2)
                            if bi not in bins_by_sample[sid]:
                                bins_by_sample[sid][bi] = []
                            bins_by_sample[sid][bi].append((position, float(max_value)))

                # Aggregate per sample: sum values and use position of max contributor
                for sid in sample_ids:
                    x_coords = []
                    y_coords = []
                    for bi in sorted(bins_by_sample[sid].keys()):
                        contributions = bins_by_sample[sid][bi]
                        total_value = sum(val for pos, val in contributions)
                        # Use position from the strand that contributed more
                        max_contributor = max(contributions, key=lambda x: x[1])
                        x_coords.append(max_contributor[0])
                        y_coords.append(total_value)
                    if x_coords:
                        feature_dict = {
                            "type": type_picked, "color": color, "alpha": alpha,
                            "fill_alpha": fill_alpha, "size": size, "title": title,
                            "x": x_coords, "y": y_coords, "has_stats": False,
                            "has_sequences": False,
                            "is_relative_scaled": is_relative_scaled,
                        }
                        result[sid].append(feature_dict)

            elif is_contig_table:
                # Contig-level: bin once, share across all samples
                # Contig tables are always curves (gc_content, gc_skew, repeats)
                x_coords, y_coords = _rle_weighted_bin_sql(
                    feature_table, True, xstart, xend, None, contig_id, cur,
                    type_picked=type_picked
                )
                # Apply scaling
                if is_relative_scaled and y_coords:
                    y_coords = [v / 1000.0 for v in y_coords]
                elif is_scaled and y_coords:
                    y_coords = [v / 100.0 for v in y_coords]
                if is_relative_scaled and min_relative_value > 0.0 and x_coords:
                    filtered = [i for i, y in enumerate(y_coords) if y is None or y >= min_relative_value]
                    x_coords = [x_coords[i] for i in filtered]
                    y_coords = [y_coords[i] for i in filtered]
                if x_coords:
                    for sid in sample_ids:
                        feature_dict = {
                            "type": type_picked, "color": color, "alpha": alpha,
                            "fill_alpha": fill_alpha, "size": size, "title": title,
                            "x": list(x_coords), "y": list(y_coords), "has_stats": False,
                            "has_sequences": False,
                            "is_relative_scaled": is_relative_scaled,
                        }
                        result[sid].append(feature_dict)

            else:
                # Sample-level table: batch bin all samples at once
                is_bars = type_picked == "bars"
                binned = _rle_weighted_bin_batch_sql(
                    feature_table, xstart, xend, sample_ids, contig_id, cur,
                    has_sequences=has_sequences, type_picked=type_picked
                )
                for sid in sample_ids:
                    if is_bars:
                        if has_sequences:
                            x_coords, y_coords, width_coords_b, fps_b, lps_b, seq_coords_b, prev_coords_b = binned[sid]
                        else:
                            x_coords, y_coords, width_coords_b, fps_b, lps_b = binned[sid]
                            seq_coords_b = []
                            prev_coords_b = []
                    elif has_sequences:
                        x_coords, y_coords, seq_coords_b, prev_coords_b = binned[sid]
                    else:
                        x_coords, y_coords = binned[sid]
                        seq_coords_b = []
                        prev_coords_b = []
                    # Apply scaling
                    if is_relative_scaled and y_coords:
                        y_coords = [v / 1000.0 for v in y_coords]
                    elif is_scaled and y_coords:
                        y_coords = [v / 100.0 for v in y_coords]
                    # Filter by min_relative_value — keep all parallel arrays aligned
                    if is_relative_scaled and min_relative_value > 0.0 and x_coords:
                        filtered = [i for i, y in enumerate(y_coords) if y is None or y >= min_relative_value]
                        x_coords = [x_coords[i] for i in filtered]
                        y_coords = [y_coords[i] for i in filtered]
                        if seq_coords_b:
                            seq_coords_b = [seq_coords_b[i] for i in filtered]
                        if prev_coords_b:
                            prev_coords_b = [prev_coords_b[i] for i in filtered]
                        if is_bars:
                            width_coords_b = [width_coords_b[i] for i in filtered]
                            fps_b = [fps_b[i] for i in filtered]
                            lps_b = [lps_b[i] for i in filtered]
                    if x_coords:
                        feature_dict = {
                            "type": type_picked, "color": color, "alpha": alpha,
                            "fill_alpha": fill_alpha, "size": size, "title": title,
                            "x": x_coords, "y": y_coords, "has_stats": False,
                            "has_sequences": has_sequences,
                            "is_relative_scaled": is_relative_scaled,
                        }
                        if is_bars:
                            feature_dict["width"] = width_coords_b
                            feature_dict["first_pos"] = fps_b
                            feature_dict["last_pos"] = lps_b
                        if has_sequences:
                            feature_dict["sequence"] = seq_coords_b
                            feature_dict["sequence_prevalence"] = prev_coords_b
                        result[sid].append(feature_dict)

        else:
            # --- FULL RESOLUTION PATH ---
            # Build position filter
            position_filter = ""
            extra_params = []
            if xstart is not None and xend is not None:
                position_filter = f" AND {_LP} >= ? AND First_position <= ?"
                extra_params = [xstart, xend]

            if feature_table == "Feature_primary_reads":
                # Batch query for plus strand
                params = list(sample_ids) + [contig_id] + extra_params
                cur.execute(
                    f"SELECT Sample_id, First_position, {_LP}, Value "
                    f"FROM Feature_primary_reads_plus_only "
                    f"WHERE Sample_id IN ({placeholders}) AND Contig_id=?{position_filter} "
                    f"ORDER BY Sample_id, First_position",
                    tuple(params)
                )
                all_plus = cur.fetchall()

                # Batch query for minus strand
                cur.execute(
                    f"SELECT Sample_id, First_position, {_LP}, Value "
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
                    expanded = _expand_rle_rows(data_rows, type_picked, has_stats, is_scaled, xstart, xend, False, min_relative_value)
                    if expanded is not None:
                        feature_dict = {
                            "type": type_picked, "color": color, "alpha": alpha,
                            "fill_alpha": fill_alpha, "size": size, "title": title,
                        }
                        feature_dict.update(expanded)
                        feature_dict["is_relative_scaled"] = False  # primary_reads not relative
                        result[sid].append(feature_dict)

            elif has_stats and has_sequences:
                params = list(sample_ids) + [contig_id] + extra_params
                cur.execute(
                    f"SELECT Sample_id, First_position, {_LP}, Value, Mean, Median, Std, Sequence, Sequence_prevalence "
                    f"FROM {feature_table} "
                    f"WHERE Sample_id IN ({placeholders}) AND Contig_id=?{position_filter} "
                    f"ORDER BY Sample_id, First_position",
                    tuple(params)
                )
                all_rows = cur.fetchall()

                # Group by sample_id
                rows_by_sample = {sid: [] for sid in sample_ids}
                for sid, first, last, val, mean, median, std, seq, seq_prev in all_rows:
                    if sid in rows_by_sample:
                        rows_by_sample[sid].append((first, last, val, mean, median, std, seq, seq_prev))

                for sid in sample_ids:
                    expanded = _expand_rle_rows(rows_by_sample[sid], type_picked, has_stats, is_scaled, xstart, xend, is_relative_scaled, min_relative_value, has_sequences)
                    if expanded is not None:
                        feature_dict = {
                            "type": type_picked, "color": color, "alpha": alpha,
                            "fill_alpha": fill_alpha, "size": size, "title": title,
                        }
                        feature_dict.update(expanded)
                        feature_dict["is_relative_scaled"] = is_relative_scaled
                        result[sid].append(feature_dict)

            elif has_stats:
                params = list(sample_ids) + [contig_id] + extra_params
                cur.execute(
                    f"SELECT Sample_id, First_position, {_LP}, Value, Mean, Median, Std "
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
                    expanded = _expand_rle_rows(rows_by_sample[sid], type_picked, has_stats, is_scaled, xstart, xend, is_relative_scaled, min_relative_value)
                    if expanded is not None:
                        feature_dict = {
                            "type": type_picked, "color": color, "alpha": alpha,
                            "fill_alpha": fill_alpha, "size": size, "title": title,
                        }
                        feature_dict.update(expanded)
                        feature_dict["is_relative_scaled"] = is_relative_scaled
                        result[sid].append(feature_dict)

            elif has_sequences:
                params = list(sample_ids) + [contig_id] + extra_params
                cur.execute(
                    f"SELECT Sample_id, First_position, {_LP}, Value, Sequence, Sequence_prevalence "
                    f"FROM {feature_table} "
                    f"WHERE Sample_id IN ({placeholders}) AND Contig_id=?{position_filter} "
                    f"ORDER BY Sample_id, First_position",
                    tuple(params)
                )
                all_rows = cur.fetchall()

                # Group by sample_id
                rows_by_sample = {sid: [] for sid in sample_ids}
                for sid, first, last, val, seq, seq_prev in all_rows:
                    if sid in rows_by_sample:
                        rows_by_sample[sid].append((first, last, val, seq, seq_prev))

                for sid in sample_ids:
                    expanded = _expand_rle_rows(rows_by_sample[sid], type_picked, has_stats, is_scaled, xstart, xend, is_relative_scaled, min_relative_value, has_sequences)
                    if expanded is not None:
                        feature_dict = {
                            "type": type_picked, "color": color, "alpha": alpha,
                            "fill_alpha": fill_alpha, "size": size, "title": title,
                        }
                        feature_dict.update(expanded)
                        feature_dict["is_relative_scaled"] = is_relative_scaled
                        result[sid].append(feature_dict)

            else:
                if is_contig_table:
                    # Contig-level table: query once (no Sample_id), duplicate for all samples
                    params = [contig_id] + extra_params
                    cur.execute(
                        f"SELECT First_position, {_LP}, Value FROM {feature_table} "
                        f"WHERE Contig_id=?{position_filter} ORDER BY First_position",
                        tuple(params)
                    )
                    contig_rows = cur.fetchall()
                    expanded = _expand_rle_rows(contig_rows, type_picked, has_stats, is_scaled, xstart, xend, is_relative_scaled, min_relative_value)
                    if expanded is not None:
                        for sid in sample_ids:
                            feature_dict = {
                                "type": type_picked, "color": color, "alpha": alpha,
                                "fill_alpha": fill_alpha, "size": size, "title": title,
                            }
                            feature_dict.update(expanded)
                            feature_dict["is_relative_scaled"] = is_relative_scaled
                            result[sid].append(feature_dict)
                else:
                    # Sample-level table: batch query
                    params = list(sample_ids) + [contig_id] + extra_params
                    cur.execute(
                        f"SELECT Sample_id, First_position, {_LP}, Value "
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
                        expanded = _expand_rle_rows(rows_by_sample[sid], type_picked, has_stats, is_scaled, xstart, xend, is_relative_scaled, min_relative_value)
                        if expanded is not None:
                            feature_dict = {
                                "type": type_picked, "color": color, "alpha": alpha,
                                "fill_alpha": fill_alpha, "size": size, "title": title,
                            }
                            feature_dict.update(expanded)
                            feature_dict["is_relative_scaled"] = is_relative_scaled
                            result[sid].append(feature_dict)

    return result


### Parsing features
def parse_requested_features(list_features):
    """Parse requested features, expanding modules to individual features.

    Accepts a mix of module names and individual feature names.
    Module names are case-insensitive and can include:
    - "coverage" or "Coverage" -> primary_reads, secondary_reads, supplementary_reads
    - "phagetermini" or "Phage termini" -> coverage_reduced, reads_starts, reads_ends, tau + Repeats
    - "assemblycheck" or "Assembly check" -> all assembly check features
    - "genome" or "Genome" -> Repeat count + Max repeat identity + GC content + GC skew

    Returns tuple of (deduplicated list of individual feature names, include_repeat_count, include_repeat_identity).
    Note: GC content and GC skew are returned as regular features in the list.
    """
    features = []
    include_repeat_count = False
    include_repeat_identity = False

    for item in list_features:
        item_lower = item.lower().strip()

        # Module: Genome
        if item_lower in ["genome"]:
            include_repeat_count = True
            include_repeat_identity = True
            features.extend(["GC content", "GC skew"])
        # Module: Coverage
        elif item_lower in ["coverage"]:
            features.extend(["Primary alignments", "Other alignments"])
        # Module: Phage termini / phagetermini
        elif item_lower in ["phage termini", "phagetermini", "phage_termini"]:
            features.extend(["Coverage reduced", "Reads termini", "Read termini transformation"])
            include_repeat_count = True
            include_repeat_identity = True
        # Module: Assembly check / assemblycheck
        elif item_lower in ["assembly check", "assemblycheck", "assembly_check"]:
            features.extend(["Clippings", "Indels", "Mismatches", "Read lengths", "Insert sizes", "Bad orientations"])
        # Handle individual repeat subplot buttons
        elif item_lower in ["repeat count"]:
            include_repeat_count = True
        elif item_lower in ["max repeat identity"]:
            include_repeat_identity = True
        # Handle legacy "Repeats" (also accept "duplications")
        elif item_lower in ["repeats", "repeat", "duplications", "duplication"]:
            include_repeat_count = True
            include_repeat_identity = True
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
    return deduped_features, include_repeat_count, include_repeat_identity


### Function to generate the bokeh plot
def generate_bokeh_plot_per_sample(conn, list_features, contig_name, sample_name, xstart=None, xend=None, subplot_size=100, genbank_path=None, feature_types=None, use_phage_colors=False, plot_isoforms=True, plot_sequence=False, same_y_scale=False, genemap_size=None, downsample_threshold=None, max_genemap_window=None, max_sequence_window=None, min_relative_value=0.0):
    """Generate a Bokeh plot for a single sample."""
    cur = conn.cursor()

    # Get contig characteristics
    contig_id, locus_name, locus_size = get_contig_info(cur, contig_name)
    print(f"Locus {locus_name} validated ({locus_size} bp)", flush=True)

    # --- Main gene annotation plot (only if genbank provided and window <= 100kb) ---
    shared_xrange = Range1d(0, locus_size)
    if xstart is not None and xend is not None:
        shared_xrange.start = xstart
        shared_xrange.end = xend

    _genemap_threshold = max_genemap_window if max_genemap_window is not None else 100_000
    annotation_fig = None
    if genbank_path and xstart is not None and xend is not None and (xend - xstart) <= _genemap_threshold:
        annotation_fig = make_bokeh_genemap(
            conn, contig_id, locus_name, locus_size,
            genemap_size if genemap_size is not None else subplot_size,
            shared_xrange, xstart, xend,
            feature_types=feature_types, use_phage_colors=use_phage_colors, plot_isoforms=plot_isoforms
        )
    elif genbank_path and xstart is not None and xend is not None and (xend - xstart) > _genemap_threshold:
        print(f"Gene map not plotted: window > {_genemap_threshold} bp", flush=True)

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
    _seq_threshold = max_sequence_window if max_sequence_window is not None else 1_000
    if plot_sequence and xstart is not None and xend is not None and (xend - xstart) <= _seq_threshold:
        seq_subplot = make_bokeh_sequence_subplot(conn, contig_name, xstart, xend, subplot_size // 2, shared_xrange)
        if seq_subplot:
            subplots.append(seq_subplot)
    elif plot_sequence and xstart is not None and xend is not None and (xend - xstart) > _seq_threshold:
        print(f"Sequence not plotted: window > {_seq_threshold} bp", flush=True)

    requested_features, include_repeat_count, include_repeat_identity = parse_requested_features(list_features)

    # Add Repeats subplots if requested (contig-level, sample-independent)
    if include_repeat_count or include_repeat_identity:
        try:
            count_dicts, identity_dicts = get_repeats_aggregated_data(cur, contig_id, xstart, xend)
            # Track 1: Repeat count (direct + inverted as separate legend entries)
            if include_repeat_count and count_dicts:
                subplot = make_bokeh_subplot(count_dicts, subplot_size, shared_xrange, show_tooltips=True)
                if subplot is not None:
                    subplots.append(subplot)
            # Track 2: Repeat max identity (direct + inverted as separate legend entries)
            if include_repeat_identity and identity_dicts:
                subplot = make_bokeh_subplot(identity_dicts, subplot_size, shared_xrange, show_tooltips=True)
                if subplot is not None:
                    subplots.append(subplot)
        except Exception as e:
            print(f"Error processing Repeats: {e}", flush=True)

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
                    variable_metadata=metadata_cache.get(feature),
                    downsample_threshold=downsample_threshold,
                    min_relative_value=min_relative_value
                )
                subplot_feature = make_bokeh_subplot(list_feature_dict, subplot_size, shared_xrange, show_tooltips=True)
                if subplot_feature is not None:
                    subplots.append(subplot_feature)
            except Exception as e:
                print(f"Error processing feature '{feature}': {e}", flush=True)

    # Subplot names whose y-axis should be relative to Primary alignments max
    PRIMARY_RELATIVE_SUBPLOTS = {
        "Primary alignments", "Alignments by strand", "Other alignments", "Clippings", "Indels",
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
                    variable_metadata=metadata_cache.get(feature),
                    downsample_threshold=downsample_threshold,
                    min_relative_value=min_relative_value
                )
                subplot_feature = make_bokeh_subplot(list_feature_dict, subplot_size, shared_xrange, show_tooltips=True)
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
