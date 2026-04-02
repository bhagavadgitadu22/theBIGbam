import os
import duckdb
import numpy as np
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
        tooltip_key = feature.qualifiers.get("_tooltip_key", "product")
        value = feature.qualifiers.get(tooltip_key)
        if value:
            return value
        return feature.type
        
    
### Plotting functions
def get_contig_info(cur, contig_name):
    cur.execute("SELECT Contig_id, Contig_name, Contig_length FROM Contig WHERE Contig_name=?", (contig_name,))
    row = cur.fetchone()
    if row is None:
        raise ValueError(f"Contig not found: {contig_name}")
    return row

def make_bokeh_genemap(conn, contig_id, locus_name, locus_size, subplot_size, shared_xrange, xstart=None, xend=None, feature_types=None, use_phage_colors=False, plot_isoforms=True, feature_label_key=None):
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
        isoform_filter = " AND (Locus_tag IS NULL OR Longest_isoform = true)"
        query = f'SELECT Annotation_id, "Start", "End", Strand, "Type", Product, "Function", Phrog, Locus_tag FROM Contig_annotation WHERE Contig_id=?{position_filter}{type_filter}{isoform_filter}'
    else:
        query = f'SELECT Annotation_id, "Start", "End", Strand, "Type", Product, "Function", Phrog, Locus_tag FROM Contig_annotation WHERE Contig_id=?{position_filter}{type_filter}'

    cur.execute(query, tuple(params))
    seq_ann_rows = cur.fetchall()

    # Fetch tooltip qualifier values from the KV table if a label key is selected
    label_map = {}
    if feature_label_key and seq_ann_rows:
        ann_ids = [row[0] for row in seq_ann_rows]
        placeholders = ','.join('?' * len(ann_ids))
        rows = cur.execute(
            f'SELECT Annotation_id, "Value" FROM Annotation_qualifier '
            f'WHERE "Key" = ? AND Annotation_id IN ({placeholders})',
            [feature_label_key] + ann_ids
        ).fetchall()
        label_map = {aid: val for aid, val in rows}

    sequence_annotations = []
    for ann_id, start, end, strand, ftype, product, function, phrog, locus_tag in seq_ann_rows:
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
        # Set tooltip key and value
        if feature_label_key:
            qualifiers['_tooltip_key'] = feature_label_key
            if ann_id in label_map:
                qualifiers[feature_label_key] = label_map[ann_id]
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

            # Add codon annotation if available
            if "codon_category" in data_feature:
                data_dict["codon_category"] = data_feature["codon_category"]
                data_dict["codon_change"] = data_feature["codon_change"]
                data_dict["aa_change"] = data_feature["aa_change"]

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
            has_codon_data = any("codon_category" in d for d in feature_dict)
            tooltips = [
                ("Position", "@x{0,0}"),
                ("Value", "@y{0.00}"),
                ("Sequence", "@sequence"),
                ("Prevalence", "@sequence_prevalence_str"),
            ]
            if has_codon_data:
                tooltips.append(("Category", "@codon_category"))
                tooltips.append(("Codon", "@codon_change"))
                tooltips.append(("Amino acid", "@aa_change"))
        else:
            tooltips = [("Position", "@x{0,0}"), ("Value", "@y{0.00}")]

        hover = HoverTool(tooltips=tooltips, mode='vline')
        p.add_tools(hover)

    # A clean style like your matplotlib setup
    p.toolbar.logo = None
    p.xgrid.visible = False

    p.y_range.start = min(0, *(y for d in feature_dict for y in d["y"] if y is not None))
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


### Function to render translated amino acid sequence as colored rectangles
def make_bokeh_translated_sequence_subplot(conn, contig_name, xstart, xend, height, x_range):
    """Create a subplot showing color-coded amino acid rectangles for CDS annotations.

    Forward-strand CDS are drawn in the top half (y: 0.5–1.0),
    reverse-strand CDS in the bottom half (y: 0.0–0.5).

    Returns None if no translated annotation data is available.

    Args:
        conn: DuckDB connection
        contig_name: Name of the contig
        xstart: Start position (1-based genome position)
        xend: End position
        height: Height of the subplot in pixels
        x_range: Shared x_range from other subplots
    """
    try:
        cur = conn.cursor()

        # Check Contig_annotation has Protein_sequence column
        cur.execute("SELECT 1 FROM information_schema.columns WHERE table_name = 'Contig_annotation' AND column_name = 'Protein_sequence'")
        if cur.fetchone() is None:
            return None

        # Load codon color/name lookup
        cur.execute("SELECT Codon, AminoAcid, AminoAcid_name, Color FROM Codon_table")
        codon_info = {}
        for codon, aa, aa_name, color in cur.fetchall():
            codon_info[codon.upper()] = (aa, aa_name, color)

        # Get contig_id
        cur.execute("SELECT Contig_id FROM Contig WHERE Contig_name = ?", (contig_name,))
        row = cur.fetchone()
        if row is None:
            return None
        contig_id = row[0]

        # Query CDS rows that overlap the visible window (only longest isoforms)
        cur.execute("""
            SELECT Start, "End", Strand, Nucleotide_sequence, Protein_sequence, Product
            FROM Contig_annotation
            WHERE Contig_id = ? AND Type = 'CDS'
              AND Protein_sequence IS NOT NULL
              AND "End" >= ? AND Start <= ?
              AND (Locus_tag IS NULL OR Longest_isoform = true)
        """, (contig_id, xstart, xend))
        cds_rows = cur.fetchall()

        if not cds_rows:
            return None

        # --- Pass 1: assign each CDS a lane using greedy interval graph coloring ---
        # Within each strand, CDS that overlap share the half-band; non-overlapping
        # CDS reuse the same lane.  Greedy coloring (sorted by start) gives an
        # optimal 2-coloring whenever possible (A→0, B→1, C→0, …).
        def assign_lanes(strand_cds_list):
            """Return {cds_idx: lane} using greedy interval-graph coloring.

            strand_cds_list: [(cds_idx, start, end), …] sorted by start.
            """
            lane_ends = []          # lane_ends[lane] = end of last CDS in that lane
            assignment = {}         # cds_idx -> lane number
            for cds_idx, cds_start, cds_end in strand_cds_list:
                assigned_lane = None
                for lane, end in enumerate(lane_ends):
                    if end < cds_start:          # no overlap — reuse lane
                        assigned_lane = lane
                        lane_ends[lane] = cds_end
                        break
                if assigned_lane is None:
                    assigned_lane = len(lane_ends)
                    lane_ends.append(cds_end)
                assignment[cds_idx] = assigned_lane
            return assignment

        forward_cds = []
        reverse_cds = []
        for cds_idx, (cds_start, cds_end, strand, nuc_seq, prot_seq, product) in enumerate(cds_rows):
            strand = int(strand) if strand is not None else 1
            if strand >= 0:
                forward_cds.append((cds_idx, cds_start, cds_end))
            else:
                reverse_cds.append((cds_idx, cds_start, cds_end))
        forward_cds.sort(key=lambda x: x[1])
        reverse_cds.sort(key=lambda x: x[1])

        fwd_lanes = assign_lanes(forward_cds)
        rev_lanes = assign_lanes(reverse_cds)
        cds_lane = {**fwd_lanes, **rev_lanes}

        fwd_num_lanes = max(fwd_lanes.values(), default=-1) + 1  # 0 if empty
        rev_num_lanes = max(rev_lanes.values(), default=-1) + 1
        total_lanes = fwd_num_lanes + rev_num_lanes

        # --- Pass 2: build rectangles with lane-aware y-coords ---
        # All CDS share the full 0.0–1.0 band: forward lanes first (top), reverse lanes below.
        lefts, rights, bottoms, tops = [], [], [], []
        colors, amino_acids, aa_names, codons_list = [], [], [], []
        pos_starts, pos_ends = [], []

        for cds_idx, (cds_start, cds_end, strand, nuc_seq, prot_seq, product) in enumerate(cds_rows):
            strand = int(strand) if strand is not None else 1
            lane = cds_lane[cds_idx]
            # Global lane index: forward lanes 0..fwd-1, then reverse lanes fwd..fwd+rev-1
            global_lane = lane if strand >= 0 else fwd_num_lanes + lane

            for i, aa in enumerate(prot_seq):
                # Compute genomic coordinates of this codon
                if strand >= 0:  # forward
                    left = cds_start + i * 3 - 0.5
                    right = cds_start + i * 3 + 2.5
                else:  # reverse
                    left = cds_end - (i + 1) * 3 + 1 - 0.5
                    right = cds_end - i * 3 + 0.5

                # Clip to visible window
                if right < xstart or left > xend:
                    continue

                # Extract the codon triplet from nucleotide sequence
                codon_str = nuc_seq[i * 3:i * 3 + 3].upper() if nuc_seq and i * 3 + 3 <= len(nuc_seq) else "???"

                info = codon_info.get(codon_str, (aa, 'Unknown', '#999999'))

                # Compute y-coords: full 0.0–1.0 band divided among all lanes
                lane_height = 1.0 / total_lanes
                top_y = 1.0 - global_lane * lane_height
                bottom_y = top_y - lane_height

                lefts.append(left)
                rights.append(right)
                bottoms.append(bottom_y)
                tops.append(top_y)
                colors.append(info[2])
                amino_acids.append(aa)
                aa_names.append(info[1])
                codons_list.append(codon_str)
                pos_starts.append(int(left + 0.5))
                pos_ends.append(int(right - 0.5))

        if not lefts:
            return None

        source = ColumnDataSource(data=dict(
            left=lefts, right=rights, bottom=bottoms, top=tops,
            color=colors, amino_acid=amino_acids, amino_acid_name=aa_names,
            codon=codons_list, position_start=pos_starts, position_end=pos_ends,
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
            ("Position", "@position_start{0,0}–@position_end{0,0}"),
            ("Codon", "@codon"),
            ("Amino acid", "@amino_acid (@amino_acid_name)"),
        ])
        p.add_tools(hover)

        # Match styling from make_bokeh_sequence_subplot
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

    except Exception as e:
        print(f"  WARNING: Could not create translated sequence subplot: {e}", flush=True)
        import traceback
        traceback.print_exc()
        return None



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


### Default threshold: windows larger than this use zoom-level binning instead of base resolution
_DEFAULT_MAX_BASE_RESOLUTION = 10_000  # 10 kb


# ============================================================================
# BLOB-based feature data retrieval
# ============================================================================

def _get_feature_blob(cur, contig_id, sample_id, feature_name):
    """Fetch a single feature BLOB from the Feature_blob table.

    Returns raw bytes or None if not found.
    """
    from thebigbam.database.blob_decoder import feature_name_to_id
    fid = feature_name_to_id(feature_name)
    if fid is None:
        return None
    cur.execute(
        "SELECT Data FROM Feature_blob WHERE Contig_id=? AND Sample_id=? AND Feature_id=?",
        (contig_id, sample_id, fid)
    )
    row = cur.fetchone()
    return row[0] if row else None


def _get_feature_blobs_batch(cur, contig_id, sample_ids, feature_name):
    """Fetch feature BLOBs for multiple samples in one query.

    Returns dict mapping sample_id -> raw bytes. Missing samples are omitted.
    """
    from thebigbam.database.blob_decoder import feature_name_to_id
    fid = feature_name_to_id(feature_name)
    if fid is None:
        return {}
    placeholders = ", ".join(["?"] * len(sample_ids))
    cur.execute(
        f"SELECT Sample_id, Data FROM Feature_blob WHERE Contig_id=? AND Feature_id=? AND Sample_id IN ({placeholders})",
        [contig_id, fid] + list(sample_ids)
    )
    return {row[0]: row[1] for row in cur.fetchall()}


def _get_contig_blob(cur, contig_id, feature_name):
    """Fetch a single Contig_blob entry for a contig-level feature.

    Contig_blob stores contig-level features (gc_content, gc_skew, repeat stats).
    Unlike Feature_blob, Contig_blob has no sample dimension.

    Returns raw bytes or None if not found.
    """
    from thebigbam.database.blob_decoder import contig_blob_name_to_id
    fid = contig_blob_name_to_id(feature_name)
    if fid is None:
        return None
    cur.execute(
        "SELECT Data FROM Contig_blob WHERE Contig_id=? AND Feature_id=?",
        (contig_id, fid)
    )
    row = cur.fetchone()
    return row[0] if row else None


def _blob_to_feature_dict(blob_bytes, type_picked, xstart=None, xend=None, max_base_resolution=None):
    """Decode a BLOB and format for Bokeh plotting.

    Args:
        max_base_resolution: Window size (bp) below which base resolution is used.
            Above this, zoom levels are tried. Default: 10_000.

    Returns dict with x, y, and optional metadata arrays, or None if empty.
    """
    from thebigbam.database.blob_decoder import decode_blob, decode_zoom_by_bin_size

    threshold = max_base_resolution if max_base_resolution is not None else _DEFAULT_MAX_BASE_RESOLUTION
    window = (xend - xstart) if (xstart is not None and xend is not None) else 0

    # Try zoom levels for windows larger than threshold
    # Pick smallest bin size that keeps points under ~10,000
    if window > threshold:
        for bin_size in [100, 1000, 10000]:
            if window // bin_size <= 10_000:
                data = decode_zoom_by_bin_size(blob_bytes, bin_size)
                if data is not None:
                    return _format_zoom_for_bokeh(data, type_picked, xstart, xend)
        # Fallback to coarsest zoom
        data = decode_zoom_by_bin_size(blob_bytes, 10000)
        if data is not None:
            return _format_zoom_for_bokeh(data, type_picked, xstart, xend)

    # Base resolution — window is small enough for per-position data
    data = decode_blob(blob_bytes)
    if len(data.get("x", [])) == 0:
        return None

    x = data["x"]
    y = data["y"]

    # Slice to window if specified
    if xstart is not None and xend is not None:
        # x is 0-indexed in BLOB, positions are 1-indexed in display
        mask = (x >= max(0, xstart - 1)) & (x <= (xend - 1))
        x = x[mask]
        y = y[mask]
        # Also slice metadata arrays
        for key in list(data.keys()):
            if key not in ("x", "y") and hasattr(data[key], '__len__') and len(data[key]) == len(mask):
                if isinstance(data[key], np.ndarray):
                    data[key] = data[key][mask]
                elif isinstance(data[key], list):
                    indices = [i for i, m in enumerate(mask) if m]
                    data[key] = [data[key][i] for i in indices]

    if len(x) == 0:
        return None

    # Convert 0-indexed to 1-indexed for display
    result = {
        "x": (x + 1).tolist(),
        "y": y.tolist(),
    }

    # Width for bars
    if type_picked == "bars":
        result["width"] = [1] * len(x)
        result["first_pos"] = result["x"]
        result["last_pos"] = result["x"]

    # Stats
    if "mean" in data:
        result["mean"] = data["mean"].tolist()
        result["median"] = data["median"].tolist()
        result["std"] = data["std"].tolist()

    # Sequences
    if "sequence" in data:
        result["sequence"] = data["sequence"]
        result["sequence_prevalence"] = data["sequence_prevalence"].tolist() if hasattr(data["sequence_prevalence"], 'tolist') else data["sequence_prevalence"]

    # Codons
    if "codon_category" in data:
        result["codon_category"] = data["codon_category"]
        result["codon_change"] = data["codon_change"]
        result["aa_change"] = data["aa_change"]

    return result


def _format_zoom_for_bokeh(zoom_data, type_picked, xstart=None, xend=None):
    """Format zoom level data for Bokeh plotting."""
    bin_starts = zoom_data["bin_start"]
    bin_ends = zoom_data["bin_end"]

    # Filter to window
    if xstart is not None and xend is not None:
        mask = (bin_ends >= max(0, xstart - 1)) & (bin_starts <= (xend - 1))
        bin_starts = bin_starts[mask]
        bin_ends = bin_ends[mask]
        for key in zoom_data:
            if key not in ("bin_start", "bin_end"):
                if isinstance(zoom_data[key], np.ndarray):
                    zoom_data[key] = zoom_data[key][mask]

    if len(bin_starts) == 0:
        return None

    # Convert to 1-indexed midpoints
    midpoints = ((bin_starts + bin_ends) / 2.0 + 1).tolist()

    if type_picked == "bars":
        # For bars: use max to preserve spikes
        y_values = zoom_data.get("max", zoom_data.get("mean", [])).tolist()
        widths = (bin_ends - bin_starts + 1).tolist()
        return {
            "x": midpoints,
            "y": y_values,
            "width": widths,
            "first_pos": (bin_starts + 1).tolist(),
            "last_pos": (bin_ends + 1).tolist(),
        }
    else:
        # For curves: use mean
        y_values = zoom_data.get("mean", []).tolist()
        return {
            "x": midpoints,
            "y": y_values,
        }


### Function to get features of one variable
def get_feature_data(cur, feature, contig_id, sample_id, xstart=None, xend=None, variable_metadata=None, max_base_resolution=None, min_relative_value=0.0):
    """Get feature data for plotting from BLOB storage.

    Decodes Feature_blob (sample-level) or Contig_blob (contig-level) entries.
    Uses zoom levels for large windows, base resolution for small windows.

    Args:
        cur: DuckDB cursor
        feature: Feature name to query
        contig_id: Contig ID
        sample_id: Sample ID
        xstart: Optional start position for filtering
        xend: Optional end position for filtering
        variable_metadata: Optional cached result from get_variable_metadata()
        max_base_resolution: Window size (bp) below which base resolution is used
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
            "is_relative_scaled": False,
        }

        # Detect contig-level table (no Sample_id column)
        is_contig_table = feature_table.startswith("Contig_")

        # --- BLOB PATH: try Feature_blob table first ---
        if feature_table == "Feature_blob" and not is_contig_table:
            # Extract feature variable name from title→Variable_name mapping
            cur.execute(
                "SELECT Variable_name FROM Variable WHERE Title=? AND Feature_table_name='Feature_blob'",
                (title,)
            )
            vname_row = cur.fetchone()
            feature_var_name = vname_row[0] if vname_row else None

            if feature_var_name:
                blob_bytes = _get_feature_blob(cur, contig_id, sample_id, feature_var_name)
                if blob_bytes is not None:
                    blob_dict = _blob_to_feature_dict(blob_bytes, type_picked, xstart, xend, max_base_resolution)

                    if blob_dict is not None:
                        feature_dict.update(blob_dict)
                        feature_dict["is_relative_scaled"] = False  # BLOB values already descaled
                        feature_dict["has_stats"] = "mean" in blob_dict
                        feature_dict["has_sequences"] = "sequence" in blob_dict
                        list_feature_dict.append(feature_dict)
                    continue  # Skip legacy path for this feature

            # BLOB feature not found or failed — skip (no legacy fallback)
            continue

        # --- CONTIG BLOB PATH: contig-level features (GC content, GC skew, repeats) ---
        if feature_table == "Contig_blob":
            cur.execute(
                "SELECT Variable_name FROM Variable WHERE Title=? AND Feature_table_name='Contig_blob'",
                (title,)
            )
            vname_row = cur.fetchone()
            contig_feature_name = vname_row[0] if vname_row else None

            if contig_feature_name:
                blob_bytes = _get_contig_blob(cur, contig_id, contig_feature_name)
                if blob_bytes is not None:
                    # GC content/skew base resolution is already windowed (500/1000 bp),
                    # so use base resolution up to 10 Mbp window
                    res_threshold = 10_000_000 if contig_feature_name in ("gc_content", "gc_skew") else max_base_resolution
                    blob_dict = _blob_to_feature_dict(blob_bytes, type_picked, xstart, xend, res_threshold)
                    if blob_dict is not None:
                        feature_dict.update(blob_dict)
                        feature_dict["is_relative_scaled"] = False
                        feature_dict["has_stats"] = False
                        feature_dict["has_sequences"] = False
                        list_feature_dict.append(feature_dict)
            continue

        # All features should be handled by BLOB paths above.
        # If we reach here, the feature is not supported or missing.

    return list_feature_dict




def get_feature_data_batch(cur, feature, contig_id, sample_ids, xstart=None, xend=None, variable_metadata=None, max_base_resolution=None, min_relative_value=0.0):
    """Get feature data for multiple samples in a single batch query.

    Decodes Feature_blob (sample-level) or Contig_blob (contig-level) entries.
    Uses zoom levels for large windows, base resolution for small windows.

    Args:
        cur: DuckDB cursor
        feature: Subplot name to query
        contig_id: Contig ID
        sample_ids: List of sample IDs to fetch data for
        xstart: Optional start position for filtering
        xend: Optional end position for filtering
        variable_metadata: Optional cached result from get_variable_metadata()
        max_base_resolution: Window size (bp) below which base resolution is used

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

    for row in rows:
        type_picked, color, alpha, fill_alpha, size, title, feature_table = row

        # --- BLOB path: fetch all samples in one query ---
        if feature_table == "Feature_blob":
            cur.execute(
                "SELECT Variable_name FROM Variable WHERE Title=? AND Feature_table_name='Feature_blob'",
                (title,)
            )
            vname_row = cur.fetchone()
            feature_var_name = vname_row[0] if vname_row else None

            if feature_var_name:
                # Single query for all samples instead of N individual queries
                blobs_by_sample = _get_feature_blobs_batch(cur, contig_id, sample_ids, feature_var_name)
                for sid, blob_bytes in blobs_by_sample.items():
                    blob_dict = _blob_to_feature_dict(blob_bytes, type_picked, xstart, xend, max_base_resolution)
                    if blob_dict is not None:
                        feature_dict = {
                            "type": type_picked, "color": color, "alpha": alpha,
                            "fill_alpha": fill_alpha, "size": size, "title": title,
                            "is_relative_scaled": False,
                            "has_stats": "mean" in blob_dict,
                            "has_sequences": "sequence" in blob_dict,
                        }
                        feature_dict.update(blob_dict)
                        result[sid].append(feature_dict)
            continue

        # --- CONTIG BLOB PATH: contig-level features (GC content, GC skew, repeats) ---
        if feature_table == "Contig_blob":
            cur.execute(
                "SELECT Variable_name FROM Variable WHERE Title=? AND Feature_table_name='Contig_blob'",
                (title,)
            )
            vname_row = cur.fetchone()
            contig_feature_name = vname_row[0] if vname_row else None

            if contig_feature_name:
                blob_bytes = _get_contig_blob(cur, contig_id, contig_feature_name)
                if blob_bytes is not None:
                    # GC content/skew base resolution is already windowed (500/1000 bp),
                    # so use base resolution up to 10 Mbp window
                    res_threshold = 10_000_000 if contig_feature_name in ("gc_content", "gc_skew") else max_base_resolution
                    blob_dict = _blob_to_feature_dict(blob_bytes, type_picked, xstart, xend, res_threshold)
                    if blob_dict is not None:
                        for sid in sample_ids:
                            feature_dict = {
                                "type": type_picked, "color": color, "alpha": alpha,
                                "fill_alpha": fill_alpha, "size": size, "title": title,
                                "is_relative_scaled": False,
                                "has_stats": False,
                                "has_sequences": False,
                            }
                            feature_dict.update(blob_dict)
                            result[sid].append(feature_dict)
            continue

        # All features should be handled by BLOB paths above.
        # If we reach here, the feature is not supported or missing.

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

    Returns deduplicated list of individual feature names (subplot names).
    All features including repeats are returned as regular features.
    """
    features = []

    for item in list_features:
        item_lower = item.lower().strip()

        # Module: Genome
        if item_lower in ["genome"]:
            features.extend(["Repeat count", "Max repeat identity", "GC content", "GC skew"])
        # Module: Coverage
        elif item_lower in ["coverage"]:
            features.extend(["Primary alignments", "Other alignments"])
        # Module: Phage termini / phagetermini
        elif item_lower in ["phage termini", "phagetermini", "phage_termini"]:
            features.extend(["Coverage reduced", "Reads termini", "Read termini transformation", "Repeat count", "Max repeat identity"])
        # Module: Assembly check / assemblycheck
        elif item_lower in ["assembly check", "assemblycheck", "assembly_check"]:
            features.extend(["Clippings", "Indels", "Mismatches", "Read lengths", "Insert sizes", "Bad orientations"])
        # Handle individual repeat subplot buttons
        elif item_lower in ["repeat count"]:
            features.append("Repeat count")
        elif item_lower in ["max repeat identity"]:
            features.append("Max repeat identity")
        # Handle legacy "Repeats" (also accept "duplications")
        elif item_lower in ["repeats", "repeat", "duplications", "duplication"]:
            features.extend(["Repeat count", "Max repeat identity"])
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
    return deduped_features


### Function to generate the bokeh plot
def generate_bokeh_plot_per_sample(conn, list_features, contig_name, sample_name, xstart=None, xend=None, subplot_size=100, genbank_path=None, feature_types=None, use_phage_colors=False, plot_isoforms=True, plot_sequence=False, plot_translated_sequence=False, same_y_scale=False, genemap_size=None, sequence_size=None, translated_sequence_size=None, max_base_resolution=None, max_genemap_window=None, max_sequence_window=None, min_relative_value=0.0, feature_label_key=None):
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
            feature_types=feature_types, use_phage_colors=use_phage_colors, plot_isoforms=plot_isoforms,
            feature_label_key=feature_label_key
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
        _seq_height = sequence_size if sequence_size is not None else subplot_size // 2
        seq_subplot = make_bokeh_sequence_subplot(conn, contig_name, xstart, xend, _seq_height, shared_xrange)
        if seq_subplot:
            subplots.append(seq_subplot)
    elif plot_sequence and xstart is not None and xend is not None and (xend - xstart) > _seq_threshold:
        print(f"Sequence not plotted: window > {_seq_threshold} bp", flush=True)

    # --- Add translated sequence subplot right after DNA sequence ---
    if plot_translated_sequence and xstart is not None and xend is not None and (xend - xstart) <= _seq_threshold:
        _trans_height = translated_sequence_size if translated_sequence_size is not None else (sequence_size if sequence_size is not None else subplot_size // 2)
        trans_subplot = make_bokeh_translated_sequence_subplot(conn, contig_name, xstart, xend, _trans_height, shared_xrange)
        if trans_subplot:
            subplots.append(trans_subplot)
    elif plot_translated_sequence and xstart is not None and xend is not None and (xend - xstart) > _seq_threshold:
        print(f"Translated sequence not plotted: window > {_seq_threshold} bp", flush=True)

    requested_features = parse_requested_features(list_features)

    # Separate contig-level features from sample-dependent features
    # Repeat features now use SQL views with standard binning (like GC content/skew)
    contig_level_features = ["GC content", "GC skew", "Repeat count", "Max repeat identity"]
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
                    max_base_resolution=max_base_resolution,
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
                    max_base_resolution=max_base_resolution,
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
