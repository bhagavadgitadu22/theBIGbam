import os
import colorsys
import hashlib
import duckdb
import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from bokeh.models import Range1d, ColumnDataSource, HoverTool, WheelZoomTool, NumeralTickFormatter, TapTool
from bokeh.layouts import gridplot
from bokeh.palettes import Viridis256
from bokeh.plotting import figure
from dna_features_viewer import BiopythonTranslator


_NULL_STRINGS = {'', 'none', 'null', 'nan'}


def make_bokeh_mag_track(conn, mag_name, height=30, shared_xrange=None):
    """Build a MAG overview track: horizontal grey line spanning the MAG with
    black vertical tick marks at each contig boundary. Contigs are arranged
    longest-first (matches MAG_contigs_association write order).

    Only boundaries and segments that overlap the visible window are rendered,
    so that large MAGs don't send unnecessary data to the browser.
    """
    from ..database.database_getters import get_mag_contigs
    members = get_mag_contigs(conn, mag_name)
    if not members:
        return None
    total_len = sum(length for _n, length, _o in members)

    # Determine the visible window (clip rendering to this range)
    if shared_xrange is not None:
        vis_start = shared_xrange.start
        vis_end = shared_xrange.end
    else:
        vis_start = 0
        vis_end = total_len

    # In MAG view (shared_xrange passed) the track shares x_range with the data
    # subplots, so giving it xpan + xwheel_zoom lets the user click/drag on the
    # track itself to slide the entire linked view. In contig view the track is
    # informational only (different coordinate space from the data plots), so
    # we omit the drag tools there.
    tools = "xpan,reset" if shared_xrange is not None else ""
    fig = figure(
        height=int(height), sizing_mode='stretch_width',
        x_range=shared_xrange if shared_xrange is not None else Range1d(0, total_len),
        y_range=Range1d(-1, 1), tools=tools,
        toolbar_location=None, outline_line_color=None,
    )
    if shared_xrange is not None:
        wheel = WheelZoomTool(dimensions='width')
        fig.add_tools(wheel)
        fig.toolbar.active_scroll = wheel
    fig.yaxis.visible = False
    fig.xgrid.grid_line_color = None
    fig.ygrid.grid_line_color = None

    # Grey backbone: filled rectangle as tall as the contig border ticks
    fig.quad(left=vis_start, right=vis_end, bottom=-0.8, top=0.8,
             color="#DCDCDC", line_color=None)

    # Collect all boundary positions
    boundaries = []
    running = 0
    for _name, length, _offset in members:
        boundaries.append(running)
        running += length
    boundaries.append(running)

    # Black tick marks: only those inside the visible window
    for x in boundaries:
        if vis_start <= x <= vis_end:
            fig.line([x, x], [-0.8, 0.8], color="black", line_width=3)

    # Hover segments: only those overlapping the visible window
    visible_segments = [
        (name, length, offset, offset + length)
        for name, length, offset in members
        if offset < vis_end and offset + length > vis_start
    ]
    if visible_segments:
        src = ColumnDataSource(dict(
            name=[s[0] for s in visible_segments],
            length=[s[1] for s in visible_segments],
            x0=[s[2] for s in visible_segments],
            x1=[s[3] for s in visible_segments],
        ))
        quads = fig.quad(left='x0', right='x1', bottom=-0.8, top=0.8,
                         source=src, fill_alpha=0.0, line_color=None)
        fig.add_tools(HoverTool(renderers=[quads], tooltips=[
            ("Contig", "@name"), ("Length", "@length{0,0}"), ("Start", "@x0{0,0}"), ("End", "@x1{0,0}"),
        ]))
    return fig


def _is_nullish(val):
    if val is None:
        return True
    if isinstance(val, str) and val.strip().lower() in _NULL_STRINGS:
        return True
    return False


def _hash_color(s):
    """Deterministic hex color — same input always yields the same hue."""
    h = int(hashlib.md5(str(s).encode()).hexdigest()[:8], 16)
    hue = (h % 360) / 360.0
    r, g, b = colorsys.hls_to_rgb(hue, 0.55, 0.65)
    return f'#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}'


def _try_numeric_map(values):
    """Return {str(raw_val): float(raw_val)} iff every non-null value parses as number."""
    out = {}
    for v in values:
        if _is_nullish(v):
            continue
        try:
            out[str(v)] = float(v)
        except (TypeError, ValueError):
            return None
    return out or None


def _gradient_palette(numeric_map):
    """Assign each distinct value a Viridis256 color based on its float position."""
    floats = list(numeric_map.values())
    vmin, vmax = min(floats), max(floats)
    span = (vmax - vmin) if vmax != vmin else 1.0
    last = len(Viridis256) - 1
    return {
        k: Viridis256[int(((f - vmin) / span) * last)]
        for k, f in numeric_map.items()
    }

### Custom translator for coloring and labeling features (with DNAFeaturesViewer python library)
# Define function-to-color mapping
def _normalize_segments(segments, fallback_start, fallback_end):
    """Return list of (start, end) 1-based inclusive tuples in genomic order.

    DuckDB returns Segments as None (unspliced) or a list of
    {'start': int, 'end': int} dicts. Unspliced features get a single-element
    list built from the bounding box so callers use one loop for both cases.
    """
    if not segments:
        return [(int(fallback_start), int(fallback_end))]
    out = [(int(s['start']), int(s['end'])) for s in segments]
    out.sort(key=lambda se: se[0])
    return out


def _codon_quads_for_aa(segments, strand, aa_index):
    """Return [(left, right), ...] quads for the aa_index-th amino acid.

    segments: genomic-ordered [(start, end), ...] 1-based inclusive (from
    _normalize_segments). Walks 3 nt starting at aa_index*3 along the
    concatenated mRNA coordinate system. Forward strand walks left-to-right;
    reverse strand walks right-to-left through reversed segments.

    Codons straddling an exon boundary yield multiple quads. Each quad's
    (left, right) is inset by 0.5 to match the existing codon-rectangle
    convention.
    """
    ordered = list(segments) if strand >= 0 else list(reversed(segments))

    codon_start_mrna = aa_index * 3
    codon_end_mrna = codon_start_mrna + 3  # exclusive

    quads = []
    cursor = 0
    for (g_start, g_end) in ordered:
        seg_len = g_end - g_start + 1
        seg_mrna_end = cursor + seg_len

        lo = max(codon_start_mrna, cursor)
        hi = min(codon_end_mrna, seg_mrna_end)
        if lo >= hi:
            cursor = seg_mrna_end
            if cursor >= codon_end_mrna:
                break
            continue

        offset_lo = lo - cursor
        offset_hi = hi - cursor - 1  # inclusive
        if strand >= 0:
            g_lo = g_start + offset_lo
            g_hi = g_start + offset_hi
        else:
            g_lo = g_end - offset_hi
            g_hi = g_end - offset_lo
        quads.append((g_lo - 0.5, g_hi + 0.5))

        cursor = seg_mrna_end
        if cursor >= codon_end_mrna:
            break

    return quads


class CustomTranslator(BiopythonTranslator):

    # Track seen unknown feature types at the class level
    _seen_unknown_types = set()

    def compute_feature_color(self, feature):
        # Custom user-defined color (from coloring rules / templates) takes priority
        custom_color = feature.qualifiers.get("_custom_color")
        if custom_color:
            return custom_color
        # Default fallback: grey for CDS, black for other types
        return "#cccccc" if feature.type == "CDS" else "#000000"

    def compute_feature_label(self, feature):
        return None  # fallback to None if missing or invalid

    def compute_feature_html(self, feature):
        tooltip_key = feature.qualifiers.get("_tooltip_key", "product")
        value = feature.qualifiers.get(tooltip_key)
        # Compare against None explicitly so numeric 0 and other falsy-but-
        # meaningful values (e.g. Main_isoform=False) still render instead of
        # falling back to the feature type.
        if value is not None and value != "":
            return str(value)
        return feature.type
        
    
### Plotting functions
def get_contig_info(cur, contig_name):
    cur.execute("SELECT Contig_id, Contig_name, Contig_length FROM Contig WHERE Contig_name=?", (contig_name,))
    row = cur.fetchone()
    if row is None:
        raise ValueError(f"Contig not found: {contig_name}")
    return row

def make_bokeh_genemap(conn, contig_id, locus_name, locus_size, subplot_size, shared_xrange, xstart=None, xend=None, feature_types=None, plot_isoforms=True, feature_label_key=None, custom_colors=None, figure_width=30):
    cur = conn.cursor()

    # Build position filter clause for annotations
    position_filter = ""
    params = [contig_id]
    if xstart is not None and xend is not None:
        position_filter = ' AND ca."End" >= ? AND ca."Start" <= ?'
        params.extend([xstart, xend])

    # Build feature type filter
    type_filter = ""
    if feature_types:
        placeholders = ','.join('?' * len(feature_types))
        type_filter = f' AND ca."Type" IN ({placeholders})'
        params.extend(feature_types)

    # Product/Function/Locus_tag live in the Annotation_qualifier KV table.
    # LEFT JOIN once per key to expose them alongside the structural columns.
    base_select = (
        'SELECT ca.Annotation_id, ca."Start", ca."End", ca.Strand, ca."Type", '
        'pq.Value AS Product, fq.Value AS Function, lq.Value AS Locus_tag '
        'FROM Contig_annotation_core ca '
        'LEFT JOIN Annotation_qualifier pq ON pq.Annotation_id = ca.Annotation_id AND pq."Key" = \'product\' '
        'LEFT JOIN Annotation_qualifier fq ON fq.Annotation_id = ca.Annotation_id AND fq."Key" = \'function\' '
        'LEFT JOIN Annotation_qualifier lq ON lq.Annotation_id = ca.Annotation_id AND lq."Key" = \'locus_tag\' '
        'WHERE ca.Contig_id=?'
    )

    # When plot_isoforms is False, filter to show only the main isoform per
    # (locus_tag, Type) pair. For spliced eukaryotic annotations this returns
    # the first mRNA seen in file order plus every CDS/exon linked to it via
    # Parent_annotation_id. Features without locus_tag always display
    # (Main_isoform is NULL for them).
    if not plot_isoforms:
        isoform_filter = ' AND (lq.Value IS NULL OR ca.Main_isoform = true)'
        query = f'{base_select}{position_filter}{type_filter}{isoform_filter}'
    else:
        query = f'{base_select}{position_filter}{type_filter}'

    cur.execute(query, tuple(params))
    seq_ann_rows = cur.fetchall()

    # Direct columns in seq_ann_rows: Strand (col 3), Type (col 4).
    direct_col_idx = {'Type': 4, 'type': 4, 'Strand': 3, 'strand': 3}
    # Whitelist of Contig_annotation_core columns that the label/color paths
    # may query directly. Prevents SQL injection via qualifier names.
    ALLOWED_CORE_COLUMNS = {'Main_isoform'}

    ann_ids = [r[0] for r in seq_ann_rows]
    placeholders = ','.join('?' * len(ann_ids)) if ann_ids else ''

    def fetch_values_for_qkey(qkey):
        """Return {annotation_id: value}: KV first, then direct row col, then core table.

        Shared by the 'Label features with:' tooltip map and the custom color
        rules, so every qualifier exposed in filtering_metadata['Annotations']
        resolves through the same fallback chain.
        """
        if not ann_ids:
            return {}
        try:
            kv_rows = cur.execute(
                f'SELECT Annotation_id, "Value" FROM Annotation_qualifier '
                f'WHERE "Key" = ? AND Annotation_id IN ({placeholders})',
                [qkey] + ann_ids
            ).fetchall()
            kv = {aid: val for aid, val in kv_rows if val is not None}
            if kv:
                return kv
        except Exception:
            pass
        col_idx = direct_col_idx.get(qkey)
        if col_idx is not None:
            return {r[0]: r[col_idx] for r in seq_ann_rows if r[col_idx] is not None}
        if qkey in ALLOWED_CORE_COLUMNS:
            try:
                rows = cur.execute(
                    f'SELECT Annotation_id, "{qkey}" FROM Contig_annotation_core '
                    f'WHERE Annotation_id IN ({placeholders})',
                    ann_ids
                ).fetchall()
                return {aid: val for aid, val in rows if val is not None}
            except Exception:
                pass
        return {}

    # Fetch tooltip values for the label key; uses the shared fallback so
    # direct columns (Type, Strand, Main_isoform, …) work alongside KV keys.
    label_map = {}
    if feature_label_key and seq_ann_rows:
        label_map = fetch_values_for_qkey(feature_label_key)

    # Build custom color map: annotation_id -> hex color (first matching rule wins)
    custom_color_map = {}
    if custom_colors and seq_ann_rows:
        from collections import defaultdict
        rules_by_key = defaultdict(list)
        for rule in custom_colors:
            rules_by_key[rule['qualifier_key']].append(rule)

        for qkey, rules in rules_by_key.items():
            values = fetch_values_for_qkey(qkey)
            for rule in rules:
                mode = rule.get('match_mode', 'exact')

                if mode == 'random':
                    # Skip nullish values; they fall through to the default color.
                    filtered = {aid: val for aid, val in values.items() if not _is_nullish(val)}
                    numeric_map = _try_numeric_map(filtered.values())
                    if numeric_map is not None:
                        value_to_color = _gradient_palette(numeric_map)
                    else:
                        distinct = {str(v) for v in filtered.values()}
                        value_to_color = {v: _hash_color(v) for v in distinct}
                    for aid, val in filtered.items():
                        if aid in custom_color_map:
                            continue
                        color = value_to_color.get(str(val))
                        if color:
                            custom_color_map[aid] = color
                    continue

                rule_val = rule['value']

                if mode in ('lt', 'gt'):
                    # Numeric comparison: both sides must parse as floats.
                    try:
                        rule_num = float(rule_val)
                    except (TypeError, ValueError):
                        continue
                    for aid, val in values.items():
                        if aid in custom_color_map or _is_nullish(val):
                            continue
                        try:
                            v_num = float(val)
                        except (TypeError, ValueError):
                            continue
                        if mode == 'lt' and v_num < rule_num:
                            custom_color_map[aid] = rule['color']
                        elif mode == 'gt' and v_num > rule_num:
                            custom_color_map[aid] = rule['color']
                    continue

                rule_val_str = str(rule_val)
                rule_val_lower = rule_val_str.lower()
                # For numeric rule values (Spinner), attempt float equality so
                # a rule value of 1.0 matches a stored int of 1. Fall back to
                # string equality when either side isn't numeric.
                rule_val_num = float(rule_val) if isinstance(rule_val, (int, float)) else None

                def _values_equal(val):
                    if rule_val_num is not None:
                        try:
                            return float(val) == rule_val_num
                        except (TypeError, ValueError):
                            pass
                    return str(val) == rule_val_str

                is_negation = mode in ('has_not', 'not_equal')
                for aid, val in values.items():
                    if aid in custom_color_map or _is_nullish(val):
                        continue
                    val_str = str(val)
                    val_lower = val_str.lower()
                    if mode == 'has':
                        if rule_val_lower in val_lower:
                            custom_color_map[aid] = rule['color']
                    elif mode == 'has_not':
                        if rule_val_lower not in val_lower:
                            custom_color_map[aid] = rule['color']
                    elif mode == 'not_equal':
                        if not _values_equal(val):
                            custom_color_map[aid] = rule['color']
                    else:  # 'exact'
                        if _values_equal(val):
                            custom_color_map[aid] = rule['color']

                if is_negation:
                    for aid in ann_ids:
                        if aid not in custom_color_map and aid not in values:
                            custom_color_map[aid] = rule['color']

    sequence_annotations = []
    for ann_id, start, end, strand, ftype, product, function, locus_tag in seq_ann_rows:
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
        if locus_tag:
            qualifiers['locus_tag'] = locus_tag
        if custom_color_map and ann_id in custom_color_map:
            qualifiers['_custom_color'] = custom_color_map[ann_id]
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
    annotation_fig = graphic_record.plot_with_bokeh(figure_width=figure_width, figure_height=40)
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

            # Add partner contig name if available (hit_identity_within_mag)
            if "partner_contig" in data_feature:
                data_dict["partner_contig"] = data_feature["partner_contig"]

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
        has_partner_contig = any("partner_contig" in d for d in feature_dict)

        if is_duplication:
            tooltips = [
                ("First position", "@first_pos{0,0}"),
                ("Last position", "@last_pos{0,0}"),
                ("Linked start", "@linked_start{0,0}"),
                ("Linked end", "@linked_end{0,0}"),
                ("Length", "@length{0,0}"),
                ("Identity", "@y{0.01}%")
            ]
        elif has_partner_contig:
            tooltips = [
                ("Position", "@x{0,0}"),
                ("Identity (%)", "@y{0.01}"),
                ("Closest contig:", "@partner_contig"),
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
_NUCLEOTIDE_COLOR_MAP = {
    'A': '#d62728', 'a': '#d62728',
    'T': '#2ca02c', 't': '#2ca02c',
    'G': '#ff7f0e', 'g': '#ff7f0e',
    'C': '#1f77b4', 'c': '#1f77b4',
}


def _build_nucleotide_figure(positions, colors, nucleotides, height, x_range):
    """Render a sequence track from aligned lists of 1-based positions, colors
    and nucleotide letters. Returns None on empty input.
    """
    if not positions:
        return None
    source = ColumnDataSource(data=dict(
        left=[p - 0.5 for p in positions],
        right=[p + 0.5 for p in positions],
        bottom=[0] * len(positions),
        top=[1] * len(positions),
        color=colors,
        nucleotide=nucleotides,
        position=positions,
        x_center=[float(p) for p in positions],
        y_center=[0.5] * len(positions),
    ))
    p = figure(height=height, x_range=x_range, y_range=Range1d(0, 1),
               tools="xpan,reset,save")
    p.quad(left='left', right='right', bottom='bottom', top='top',
           color='color', source=source, line_color=None)
    p.text(x='x_center', y='y_center', text='nucleotide', source=source,
           text_align='center', text_baseline='middle',
           text_font_size='8pt', text_color='white', text_font_style='bold')
    p.add_tools(HoverTool(tooltips=[
        ("Position", "@position{0,0}"),
        ("Nucleotide", "@nucleotide"),
    ]))
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


def make_bokeh_sequence_subplot(conn, contig_name, xstart, xend, height, x_range):
    """Colored nucleotide rectangles for a contig region. Returns None if no
    sequence data is available.
    """
    try:
        cur = conn.cursor()
        seq_start = max(xstart, 1)
        cur.execute(
            "SELECT SUBSTR(cs.Sequence, ?, ? - ? + 1) "
            "FROM Contig_sequence cs "
            "JOIN Contig c ON cs.Contig_id = c.Contig_id "
            "WHERE c.Contig_name = ?",
            (seq_start, xend, seq_start, contig_name)
        )
        row = cur.fetchone()
        if row is None or not row[0]:
            return None
        seq = row[0]
        positions = [seq_start + i for i in range(len(seq))]
        colors = [_NUCLEOTIDE_COLOR_MAP.get(nt, '#999999') for nt in seq]
        nucleotides = [nt.upper() for nt in seq]
        return _build_nucleotide_figure(positions, colors, nucleotides, height, x_range)
    except Exception:
        return None


def make_bokeh_sequence_subplot_mag(conn, mag_name, xstart, xend, height, x_range):
    """MAG-wide sequence track. For each member contig that overlaps the visible
    MAG window, fetches the overlapping substring and places nucleotide quads
    at their MAG-space positions.
    """
    from ..database.database_getters import get_mag_contigs
    try:
        members = get_mag_contigs(conn, mag_name)
        if not members:
            return None
        cur = conn.cursor()
        win_start = max(xstart, 1)
        positions, colors, nucleotides = [], [], []
        for contig_name, contig_length, offset in members:
            contig_mag_start = offset + 1
            contig_mag_end = offset + contig_length
            lo = max(win_start, contig_mag_start)
            hi = min(xend, contig_mag_end)
            if lo > hi:
                continue
            local_start = lo - offset
            local_len = hi - lo + 1
            cur.execute(
                "SELECT SUBSTR(cs.Sequence, ?, ?) "
                "FROM Contig_sequence cs "
                "JOIN Contig c ON cs.Contig_id = c.Contig_id "
                "WHERE c.Contig_name = ?",
                (local_start, local_len, contig_name),
            )
            row = cur.fetchone()
            if row is None or not row[0]:
                continue
            for i, nt in enumerate(row[0]):
                positions.append(lo + i)
                colors.append(_NUCLEOTIDE_COLOR_MAP.get(nt, '#999999'))
                nucleotides.append(nt.upper())
        return _build_nucleotide_figure(positions, colors, nucleotides, height, x_range)
    except Exception:
        return None


def _load_codon_info(cur):
    """Return {codon_upper: (aa, aa_name, aa_label, color)} from Codon_table."""
    cur.execute("SELECT Codon, AminoAcid, AminoAcid_name, AminoAcid_label, Color FROM Codon_table")
    return {codon.upper(): (aa, aa_name, aa_label, color)
            for codon, aa, aa_name, aa_label, color in cur.fetchall()}


def _build_translated_figure(cds_entries, codon_info, xstart, xend, height, x_range):
    """Render CDS quads from a list of entries (start, end, strand, nuc_seq,
    prot_seq, segments_normalized) all in the target coordinate space. Returns
    None if no quads end up in the visible window.

    Segments are 1-based inclusive (start, end) tuples, ordered by start.
    """
    if not cds_entries:
        return None

    def _assign_lanes(all_cds):
        lane_ends, assignment = [], {}
        for cds_idx, cds_start, cds_end in all_cds:
            chosen = None
            for lane, end in enumerate(lane_ends):
                if end < cds_start:
                    chosen = lane
                    lane_ends[lane] = cds_end
                    break
            if chosen is None:
                chosen = len(lane_ends)
                lane_ends.append(cds_end)
            assignment[cds_idx] = chosen
        return assignment

    all_cds = [(i, e[0], e[1]) for i, e in enumerate(cds_entries)]
    all_cds.sort(key=lambda x: x[1])
    cds_lane = _assign_lanes(all_cds)
    total_lanes = max(cds_lane.values(), default=0) + 1

    lefts, rights, bottoms, tops = [], [], [], []
    colors, amino_acids, aa_names, aa_labels, codons_list = [], [], [], [], []
    pos_starts, pos_ends = [], []

    for cds_idx, (cds_start, cds_end, strand, nuc_seq, prot_seq, segments) in enumerate(cds_entries):
        strand = int(strand) if strand is not None else 1
        lane = cds_lane[cds_idx]
        for i, aa in enumerate(prot_seq):
            quads = _codon_quads_for_aa(segments, strand, i)
            if not quads:
                continue
            codon_str = nuc_seq[i*3:i*3+3].upper() if nuc_seq and i*3+3 <= len(nuc_seq) else "???"
            info = codon_info.get(codon_str, (aa, 'Unknown', '???', '#999999'))
            lane_height = 1.0 / total_lanes
            top_y = 1.0 - lane * lane_height
            bottom_y = top_y - lane_height
            for (left, right) in quads:
                if right < xstart or left > xend:
                    continue
                lefts.append(left)
                rights.append(right)
                bottoms.append(bottom_y)
                tops.append(top_y)
                colors.append(info[3])
                amino_acids.append(aa)
                aa_names.append(info[1])
                aa_labels.append(info[2])
                codons_list.append(codon_str)
                pos_starts.append(int(left + 0.5))
                pos_ends.append(int(right - 0.5))

    if not lefts:
        return None

    x_centers = [(l + r) / 2 for l, r in zip(lefts, rights)]
    y_centers = [(b + t) / 2 for b, t in zip(bottoms, tops)]
    source = ColumnDataSource(data=dict(
        left=lefts, right=rights, bottom=bottoms, top=tops,
        color=colors, amino_acid=amino_acids, amino_acid_name=aa_names,
        aa_label=aa_labels, codon=codons_list,
        position_start=pos_starts, position_end=pos_ends,
        x_center=x_centers, y_center=y_centers,
    ))
    p = figure(height=height, x_range=x_range, y_range=Range1d(0, 1),
               tools="xpan,reset,save")
    p.quad(left='left', right='right', bottom='bottom', top='top',
           color='color', source=source, line_color=None)
    p.text(x='x_center', y='y_center', text='aa_label', source=source,
           text_align='center', text_baseline='middle',
           text_font_size='7pt', text_color='white', text_font_style='bold')
    p.add_tools(HoverTool(tooltips=[
        ("Position", "@position_start{0,0}–@position_end{0,0}"),
        ("Codon", "@codon"),
        ("Amino acid", "@amino_acid (@amino_acid_name)"),
    ]))
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


def _fetch_cds_rows(cur, contig_name, xstart, xend):
    """Return raw CDS rows (start, end, strand, nuc_seq, prot_seq, segments_raw)
    for one contig, filtered to the main isoform and the visible window.
    Returns [] if Contig_annotation has no Protein_sequence column.
    """
    cur.execute(
        "SELECT 1 FROM information_schema.columns "
        "WHERE table_name = 'Contig_annotation' AND column_name = 'Protein_sequence'"
    )
    if cur.fetchone() is None:
        return []
    cur.execute("""
        SELECT ca.Start, ca."End", ca.Strand, ca.Nucleotide_sequence, ca.Protein_sequence, ca.Segments
        FROM Contig_annotation ca
        JOIN Contig c ON c.Contig_id = ca.Contig_id
        LEFT JOIN Annotation_qualifier lq ON lq.Annotation_id = ca.Annotation_id AND lq."Key" = 'locus_tag'
        WHERE c.Contig_name = ? AND ca.Type = 'CDS'
          AND ca.Protein_sequence IS NOT NULL
          AND ca."End" >= ? AND ca.Start <= ?
          AND (lq.Value IS NULL OR ca.Main_isoform = true)
    """, (contig_name, xstart, xend))
    return cur.fetchall()


def make_bokeh_translated_sequence_subplot(conn, contig_name, xstart, xend, height, x_range):
    """Color-coded amino-acid rectangles for CDS annotations on a single contig.
    Returns None if no translated annotation data is available.
    """
    try:
        cur = conn.cursor()
        rows = _fetch_cds_rows(cur, contig_name, xstart, xend)
        if not rows:
            return None
        codon_info = _load_codon_info(cur)
        entries = [
            (start, end, strand, nuc_seq, prot_seq,
             _normalize_segments(segments_raw, start, end))
            for start, end, strand, nuc_seq, prot_seq, segments_raw in rows
        ]
        return _build_translated_figure(entries, codon_info, xstart, xend, height, x_range)
    except Exception as e:
        print(f"  WARNING: Could not create translated sequence subplot: {e}", flush=True)
        import traceback
        traceback.print_exc()
        return None


def make_bokeh_translated_sequence_subplot_mag(conn, mag_name, xstart, xend, height, x_range):
    """MAG-wide translated sequence track. Collects CDS rows across every
    member contig that overlaps the visible MAG window, shifts positions +
    segments by each contig's MAG offset, then renders with a global lane
    assignment across all CDS in the window.
    """
    from ..database.database_getters import get_mag_contigs
    try:
        cur = conn.cursor()
        members = get_mag_contigs(conn, mag_name)
        if not members:
            return None
        entries = []
        for contig_name, contig_length, offset in members:
            contig_mag_start = offset + 1
            contig_mag_end = offset + contig_length
            if contig_mag_start > xend or contig_mag_end < xstart:
                continue
            local_xstart = max(xstart - offset, 1)
            local_xend = min(xend - offset, contig_length)
            for start, end, strand, nuc_seq, prot_seq, segments_raw in _fetch_cds_rows(
                cur, contig_name, local_xstart, local_xend
            ):
                seg_local = _normalize_segments(segments_raw, start, end)
                seg_mag = [(s + offset, e + offset) for s, e in seg_local]
                entries.append((start + offset, end + offset, strand, nuc_seq, prot_seq, seg_mag))
        if not entries:
            return None
        codon_info = _load_codon_info(cur)
        return _build_translated_figure(entries, codon_info, xstart, xend, height, x_range)
    except Exception as e:
        print(f"  WARNING: Could not create translated sequence subplot (MAG): {e}", flush=True)
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


def split_contig_vs_sample_features(metadata_cache, requested_features):
    """Split requested subplots into contig-level vs sample-dependent.

    A subplot is contig-level when every one of its Variable rows is backed
    by Contig_blob (no per-sample data). Falls back to sample-dependent if
    the subplot has no metadata rows.
    """
    contig_features = []
    sample_features = []
    for feature in requested_features:
        rows = metadata_cache.get(feature)
        if rows and all(r[6] == "Contig_blob" for r in rows):
            contig_features.append(feature)
        else:
            sample_features.append(feature)
    return contig_features, sample_features


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

### Default window thresholds for optional plot components
DEFAULT_GENEMAP_WINDOW = 100_000   # gene-map shown when window <= this (bp)
DEFAULT_SEQUENCE_WINDOW = 1_000    # sequence track shown when window <= this (bp)


# ============================================================================
# BLOB-based feature data retrieval
# ============================================================================

def _get_feature_blob(cur, contig_id, sample_id, feature_name):
    """Fetch Zoom_data for a single feature from Feature_blob.

    Returns raw bytes or None if not found.
    """
    from thebigbam.database.blob_decoder import feature_name_to_id
    fid = feature_name_to_id(feature_name, cur)
    if fid is None:
        return None
    cur.execute(
        "SELECT Zoom_data FROM Feature_blob WHERE Contig_id=? AND Sample_id=? AND Feature_id=?",
        (contig_id, sample_id, fid)
    )
    row = cur.fetchone()
    return row[0] if row else None


def _get_feature_blobs_batch(cur, contig_id, sample_ids, feature_name):
    """Fetch Zoom_data for multiple samples in one query.

    Returns dict mapping sample_id -> raw bytes. Missing samples are omitted.
    """
    from thebigbam.database.blob_decoder import feature_name_to_id
    fid = feature_name_to_id(feature_name, cur)
    if fid is None:
        return {}
    placeholders = ", ".join(["?"] * len(sample_ids))
    cur.execute(
        f"SELECT Sample_id, Zoom_data FROM Feature_blob WHERE Contig_id=? AND Feature_id=? AND Sample_id IN ({placeholders})",
        [contig_id, fid] + list(sample_ids)
    )
    return {row[0]: row[1] for row in cur.fetchall()}


def _get_contig_blob(cur, contig_id, feature_name):
    """Fetch Zoom_data for a contig-level feature from Contig_blob.

    Returns raw bytes or None if not found.
    """
    from thebigbam.database.blob_decoder import feature_name_to_id
    fid = feature_name_to_id(feature_name, cur)
    if fid is None:
        return None
    cur.execute(
        "SELECT Zoom_data FROM Contig_blob WHERE Contig_id=? AND Feature_id=?",
        (contig_id, fid)
    )
    row = cur.fetchone()
    return row[0] if row else None


def _get_feature_chunks(cur, contig_id, sample_id, feature_name, start_pos, end_pos):
    """Fetch only the base-resolution chunks overlapping [start_pos, end_pos) from Feature_blob_chunk.

    Returns list of (chunk_idx, raw_bytes) sorted by chunk_idx, or None if table doesn't exist.
    """
    from thebigbam.database.blob_decoder import feature_name_to_id, CHUNK_SIZE
    fid = feature_name_to_id(feature_name, cur)
    if fid is None:
        return None
    first_chunk = max(0, start_pos // CHUNK_SIZE)
    last_chunk = end_pos // CHUNK_SIZE
    try:
        cur.execute(
            "SELECT Chunk_idx, Data FROM Feature_blob_chunk "
            "WHERE Contig_id=? AND Sample_id=? AND Feature_id=? AND Chunk_idx BETWEEN ? AND ? "
            "ORDER BY Chunk_idx",
            (contig_id, sample_id, fid, first_chunk, last_chunk)
        )
        rows = cur.fetchall()
        return [(r[0], r[1]) for r in rows] if rows else None
    except Exception:
        return None


def _get_contig_chunks(cur, contig_id, feature_name, start_pos, end_pos, window_size=1):
    """Fetch only the base-resolution chunks overlapping [start_pos, end_pos) from Contig_blob_chunk.

    For windowed data (GC content/skew), pass window_size to convert genomic positions
    to window indices before computing chunk indices.

    Returns list of (chunk_idx, raw_bytes) sorted by chunk_idx, or None if table doesn't exist.
    """
    from thebigbam.database.blob_decoder import feature_name_to_id, CHUNK_SIZE
    fid = feature_name_to_id(feature_name, cur)
    if fid is None:
        return None
    # Convert genomic positions to array indices (for windowed data like GC)
    start_idx = start_pos // window_size
    end_idx = end_pos // window_size
    first_chunk = max(0, start_idx // CHUNK_SIZE)
    last_chunk = end_idx // CHUNK_SIZE
    try:
        cur.execute(
            "SELECT Chunk_idx, Data FROM Contig_blob_chunk "
            "WHERE Contig_id=? AND Feature_id=? AND Chunk_idx BETWEEN ? AND ? "
            "ORDER BY Chunk_idx",
            (contig_id, fid, first_chunk, last_chunk)
        )
        rows = cur.fetchall()
        return [(r[0], r[1]) for r in rows] if rows else None
    except Exception:
        return None


_REPEAT_IDENTITY_FEATURES = {"direct_repeat_identity", "inverted_repeat_identity"}


def _resolve_partner_lane(cur, blob_dict, feature_name):
    if "partner_contig_id" not in blob_dict:
        return
    partner_ids = list(blob_dict.pop("partner_contig_id"))
    if feature_name in _REPEAT_IDENTITY_FEATURES:
        # Events arrive as (seg_start, seg_end) pairs per repeat segment. The seg_start
        # event carries Partner_start and seg_end carries Partner_end (see db.rs
        # write_repeat_features_to_blob_static). Pair them into a single "start-end"
        # label and assign it to both events so the tooltip is uniform across the plateau.
        labels = [None] * len(partner_ids)
        for i in range(0, len(partner_ids) - 1, 2):
            s, e = partner_ids[i], partner_ids[i + 1]
            if s is None or e is None or int(s) < 0 or int(e) < 0:
                continue
            s1, e1 = int(s) + 1, int(e) + 1
            label = f"{s1}-{e1}" if s1 != e1 else f"{s1}"
            labels[i] = labels[i + 1] = label
        blob_dict["repeat_positions"] = labels
        return
    unique_ids = {int(i) for i in partner_ids if i is not None and int(i) >= 0}
    if unique_ids:
        placeholders = ",".join("?" * len(unique_ids))
        cur.execute(
            f"SELECT Contig_id, Contig_name FROM Contig WHERE Contig_id IN ({placeholders})",
            list(unique_ids),
        )
        id_to_name = {row[0]: row[1] for row in cur.fetchall()}
    else:
        id_to_name = {}
    blob_dict["partner_contig"] = [
        id_to_name.get(int(i), "unknown") if i is not None and int(i) >= 0 else None
        for i in partner_ids
    ]


def _format_chunks_for_bokeh(chunk_dict, xstart, xend, type_picked="line"):
    """Slice decoded chunk data to view window and format for Bokeh.

    chunk_dict has 0-indexed numpy arrays {"x": ..., "y": ..., "sparse": bool, ...metadata}.
    Returns dict with 1-indexed lists, or None if empty after slicing.
    """
    x = chunk_dict["x"]
    y = chunk_dict["y"]
    if len(x) == 0:
        return None
    n_orig = len(x)
    # Slice to view window (0-indexed)
    lo = xstart - 1
    hi = xend - 1
    left_idx = int(np.searchsorted(x, lo, side="right")) - 1
    right_idx = int(np.searchsorted(x, hi, side="left"))
    start = max(0, left_idx)
    stop = min(n_orig, right_idx + 1)
    x = x[start:stop]
    y = y[start:stop]
    if len(x) == 0:
        return None

    # Slice metadata arrays before zero-anchor insertion (lengths match pre-anchor x/y)
    sliced_meta = {}
    for key in list(chunk_dict.keys()):
        if key in ("x", "y", "sparse"):
            continue
        val = chunk_dict[key]
        if hasattr(val, "__len__") and len(val) == n_orig:
            if isinstance(val, np.ndarray):
                sliced_meta[key] = val[start:stop].tolist()
            elif isinstance(val, list):
                sliced_meta[key] = val[start:stop]

    # For sparse features, insert zero-anchor points between distant events
    # so Bokeh lines/areas drop to zero instead of interpolating across gaps.
    # Metadata lanes are padded with None at anchors so per-event tooltip
    # fields (partner_contig_id, repeat_positions, stats, sequence, codons…)
    # stay aligned with x/y.
    if chunk_dict.get("sparse") and type_picked != "bars" and len(x) > 0:
        new_x = []
        new_y = []
        anchor_meta = {k: [] for k in sliced_meta}

        def _emit_anchor(px):
            new_x.append(px)
            new_y.append(0.0)
            for k in anchor_meta:
                anchor_meta[k].append(None)

        for i in range(len(x)):
            if i == 0 and x[i] > 0:
                _emit_anchor(x[i] - 1)
            elif i > 0 and x[i] > x[i - 1] + 1:
                if float(y[i]) != float(y[i - 1]):
                    _emit_anchor(x[i - 1] + 1)
                    _emit_anchor(x[i] - 1)
                else:
                    # Same y across a gap. Preserve plateau only for segment
                    # pairs (repeat start/end): both sides of the gap are
                    # isolated events, not part of dense per-position runs.
                    prev_is_dense = (i >= 2 and x[i - 1] == x[i - 2] + 1)
                    next_is_dense = (i + 1 < len(x) and x[i + 1] == x[i] + 1)
                    if prev_is_dense or next_is_dense:
                        _emit_anchor(x[i - 1] + 1)
                        _emit_anchor(x[i] - 1)
            new_x.append(x[i])
            new_y.append(float(y[i]))
            for k, vals in sliced_meta.items():
                anchor_meta[k].append(vals[i])
        _emit_anchor(x[-1] + 1)
        x = np.array(new_x, dtype=np.uint32)
        y = np.array(new_y, dtype=np.float64)
        sliced_meta = anchor_meta

    result = {
        "x": (x + 1).tolist(),
        "y": y.tolist(),
    }
    result.update(sliced_meta)
    return result


def _blob_to_feature_dict(blob_bytes, type_picked, xstart=None, xend=None, max_base_resolution=None, zoom_blob_bytes=None):
    """Decode a BLOB and format for Bokeh plotting.

    Args:
        max_base_resolution: Window size (bp) below which base resolution is used.
            Above this, zoom levels are tried. Default: 10_000.
        zoom_blob_bytes: Optional standalone zoom BLOB (TBZ format). When provided and
            the window is large enough for zoom, this is decoded instead of the full BLOB.

    Returns dict with x, y, and optional metadata arrays, or None if empty.
    """
    from thebigbam.database.blob_decoder import decode_blob, decode_zoom_by_bin_size, decode_zoom_standalone

    threshold = max_base_resolution if max_base_resolution is not None else _DEFAULT_MAX_BASE_RESOLUTION
    window = (xend - xstart) if (xstart is not None and xend is not None) else 0

    # Try zoom levels for windows larger than threshold
    # Pick smallest bin size that keeps points under ~10,000
    if window > threshold:
        # Use standalone zoom blob if available (avoids fetching full Data BLOB)
        zoom_decoder = decode_zoom_standalone if zoom_blob_bytes is not None else decode_zoom_by_bin_size
        zoom_src = zoom_blob_bytes if zoom_blob_bytes is not None else blob_bytes
        for bin_size in [100, 1000, 10000]:
            if window // bin_size <= 10_000:
                data = zoom_decoder(zoom_src, bin_size)
                if data is not None:
                    return _format_zoom_for_bokeh(data, type_picked, xstart, xend)
        # Fallback to coarsest zoom
        data = zoom_decoder(zoom_src, 10000)
        if data is not None:
            return _format_zoom_for_bokeh(data, type_picked, xstart, xend)

    # Base resolution — full decode (chunk optimization handled at caller level).
    # This path is only reached for GC content/skew (small windowed BLOBs).
    data = decode_blob(blob_bytes)
    if len(data.get("x", [])) == 0:
        return None

    x = data["x"]
    y = data["y"]

    # Detect windowed data (GC content/skew) from spacing between positions
    if len(x) >= 2:
        spacing = int(x[1] - x[0])
    else:
        spacing = 1

    # Slice to window if specified.
    # Include one anchor point on each side (last x <= xstart, first x >= xend)
    # so Bokeh always has a segment crossing the viewport — critical for windowed
    # data (GC content 500bp, GC skew 1000bp) when zoomed below the window spacing.
    if xstart is not None and xend is not None:
        # x is 0-indexed in BLOB; xstart/xend are 1-indexed display positions.
        lo = xstart - 1
        hi = xend - 1
        n_orig = len(x)
        # x is sorted ascending (BLOB decode order) — use searchsorted for O(log n).
        left_idx = int(np.searchsorted(x, lo, side="right")) - 1   # last x <= lo (may be -1)
        right_idx = int(np.searchsorted(x, hi, side="left"))        # first x >= hi
        start = max(0, left_idx)
        stop = min(n_orig, right_idx + 1)
        sl = slice(start, stop)
        x = x[sl]
        y = y[sl]
        # Also slice metadata arrays with the same slice
        for key in list(data.keys()):
            if key in ("x", "y"):
                continue
            val = data[key]
            if hasattr(val, "__len__") and len(val) == n_orig:
                if isinstance(val, np.ndarray):
                    data[key] = val[sl]
                elif isinstance(val, list):
                    data[key] = val[start:stop]

    if len(x) == 0:
        return None

    # For sparse features, insert zero-anchor points between distant events
    # so Bokeh lines/areas drop to zero instead of interpolating across gaps.
    if data.get("sparse") and type_picked != "bars" and len(x) > 0:
        new_x = []
        new_y = []
        for i in range(len(x)):
            if i == 0 and x[i] > 0:
                new_x.append(x[i] - 1)
                new_y.append(0.0)
            elif i > 0 and x[i] > x[i - 1] + 1:
                new_x.append(x[i - 1] + 1)
                new_y.append(0.0)
                new_x.append(x[i] - 1)
                new_y.append(0.0)
            new_x.append(x[i])
            new_y.append(float(y[i]))
        new_x.append(x[-1] + 1)
        new_y.append(0.0)
        x = np.array(new_x, dtype=np.uint32)
        y = np.array(new_y, dtype=np.float64)

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

        # For sparse zoom data, insert zero-anchor points between non-contiguous bins
        if zoom_data.get("sparse") and len(midpoints) > 0:
            bin_size = int(bin_ends[0] - bin_starts[0] + 1) if len(bin_starts) > 0 else 1
            new_x = []
            new_y = []
            for i in range(len(midpoints)):
                if i == 0:
                    # Zero just before first bin
                    new_x.append(midpoints[i] - bin_size)
                    new_y.append(0.0)
                elif bin_starts[i] > bin_ends[i - 1] + 1:
                    # Gap between bins — insert zeros on both sides
                    new_x.append(midpoints[i - 1] + bin_size)
                    new_y.append(0.0)
                    new_x.append(midpoints[i] - bin_size)
                    new_y.append(0.0)
                new_x.append(midpoints[i])
                new_y.append(y_values[i])
            # Zero just after last bin
            new_x.append(midpoints[-1] + bin_size)
            new_y.append(0.0)
            midpoints = new_x
            y_values = new_y

        return {
            "x": midpoints,
            "y": y_values,
        }


### Function to get features of one variable
def get_feature_data(cur, feature, contig_id, sample_id, xstart=None, xend=None, variable_metadata=None, max_base_resolution=None, min_relative_value=0.0, window_for_zoom=None):
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
                # Decide whether to fetch lightweight zoom-only blob or full data
                _threshold = max_base_resolution if max_base_resolution is not None else _DEFAULT_MAX_BASE_RESOLUTION
                _window = window_for_zoom if window_for_zoom is not None else (
                    (xend - xstart) if (xstart is not None and xend is not None) else 0
                )
                _use_zoom = _window > _threshold
                if _use_zoom:
                    zoom_blob = _get_feature_blob(cur, contig_id, sample_id, feature_var_name)
                    blob_dict = _blob_to_feature_dict(None, type_picked, xstart, xend, max_base_resolution, zoom_blob_bytes=zoom_blob) if zoom_blob else None
                else:
                    # Zoomed in: fetch only 1-2 x 65kbp chunks from chunk table
                    from thebigbam.database.blob_decoder import decode_raw_chunks, decode_raw_sparse_chunks, get_scale_from_zoom_blob, is_sparse_zoom_blob
                    chunk_rows = _get_feature_chunks(cur, contig_id, sample_id, feature_var_name, xstart - 1, xend - 1)
                    zoom_blob = _get_feature_blob(cur, contig_id, sample_id, feature_var_name)
                    scale_div = get_scale_from_zoom_blob(zoom_blob) if zoom_blob else 1
                    if chunk_rows:
                        if zoom_blob and is_sparse_zoom_blob(zoom_blob):
                            blob_dict = _format_chunks_for_bokeh(decode_raw_sparse_chunks(chunk_rows, scale_div), xstart, xend, type_picked)
                        else:
                            blob_dict = _format_chunks_for_bokeh(decode_raw_chunks(chunk_rows, scale_div), xstart, xend, type_picked)
                    else:
                        blob_dict = None

                if blob_dict is not None:
                    # Filter sparse positions: zero out y values below min_relative_value * max(y)
                    if min_relative_value > 0.0 and "y" in blob_dict and blob_dict["y"]:
                        max_y = max(blob_dict["y"])
                        if max_y > 0:
                            threshold = min_relative_value * max_y
                            blob_dict["y"] = [v if v >= threshold else 0 for v in blob_dict["y"]]
                    feature_dict.update(blob_dict)
                    feature_dict["is_relative_scaled"] = False  # BLOB values already descaled
                    feature_dict["has_stats"] = "mean" in blob_dict
                    feature_dict["has_sequences"] = "sequence" in blob_dict
                    list_feature_dict.append(feature_dict)
                    continue

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
                # GC content/skew base resolution is already windowed (500/1000 bp),
                # so use base resolution up to 10 Mbp window
                is_gc = contig_feature_name in ("gc_content", "gc_skew")
                res_threshold = 10_000_000 if is_gc else max_base_resolution
                _window = window_for_zoom if window_for_zoom is not None else (
                    (xend - xstart) if (xstart is not None and xend is not None) else 0
                )
                _use_zoom = _window > (res_threshold if res_threshold is not None else _DEFAULT_MAX_BASE_RESOLUTION)
                if _use_zoom:
                    zoom_blob = _get_contig_blob(cur, contig_id, contig_feature_name)
                    blob_dict = _blob_to_feature_dict(None, type_picked, xstart, xend, res_threshold, zoom_blob_bytes=zoom_blob) if zoom_blob else None
                elif xstart is not None and xend is not None:
                    # Zoomed in: fetch chunks
                    from thebigbam.database.blob_decoder import decode_raw_chunks, decode_raw_sparse_chunks, get_scale_from_zoom_blob, is_sparse_zoom_blob, gc_window_size
                    gc_window = gc_window_size(contig_feature_name)
                    chunk_rows = _get_contig_chunks(cur, contig_id, contig_feature_name, xstart - 1, xend - 1, window_size=gc_window)
                    zoom_blob = _get_contig_blob(cur, contig_id, contig_feature_name)
                    scale_div = get_scale_from_zoom_blob(zoom_blob) if zoom_blob else 1
                    is_sparse = bool(zoom_blob and is_sparse_zoom_blob(zoom_blob))
                    if chunk_rows:
                        if is_sparse:
                            blob_dict = _format_chunks_for_bokeh(decode_raw_sparse_chunks(chunk_rows, scale_div), xstart, xend, type_picked)
                        else:
                            chunk_dict = decode_raw_chunks(chunk_rows, scale_div)
                            # For windowed data, convert array indices to 0-indexed window starts
                            # (kept as-is for _format_chunks_for_bokeh slicing), then shift to
                            # midpoints after slicing (add half - 1 because _format_chunks_for_bokeh
                            # adds 1 when converting to 1-indexed genomic coordinates).
                            if gc_window > 1:
                                chunk_dict["x"] = chunk_dict["x"] * gc_window
                            blob_dict = _format_chunks_for_bokeh(chunk_dict, xstart, xend, type_picked)
                            if gc_window > 1 and blob_dict is not None:
                                half = gc_window // 2 - 1
                                blob_dict["x"] = [v + half for v in blob_dict["x"]]
                    else:
                        blob_dict = None
                else:
                    # No window specified — fetch all chunks
                    from thebigbam.database.blob_decoder import decode_raw_chunks, decode_raw_sparse_chunks, get_scale_from_zoom_blob, is_sparse_zoom_blob, gc_window_size
                    gc_window = gc_window_size(contig_feature_name)
                    chunk_rows = _get_contig_chunks(cur, contig_id, contig_feature_name, 0, 2**31, window_size=gc_window)
                    zoom_blob = _get_contig_blob(cur, contig_id, contig_feature_name)
                    scale_div = get_scale_from_zoom_blob(zoom_blob) if zoom_blob else 1
                    is_sparse = bool(zoom_blob and is_sparse_zoom_blob(zoom_blob))
                    if chunk_rows:
                        if is_sparse:
                            raw = decode_raw_sparse_chunks(chunk_rows, scale_div)
                            blob_dict = {"x": (raw["x"] + 1).tolist(), "y": raw["y"].tolist()}
                            for k in raw:
                                if k not in ("x", "y", "sparse"):
                                    blob_dict[k] = raw[k].tolist() if hasattr(raw[k], "tolist") else raw[k]
                        else:
                            chunk_dict = decode_raw_chunks(chunk_rows, scale_div)
                            if gc_window > 1:
                                # Convert indices to 1-indexed midpoints directly
                                chunk_dict["x"] = chunk_dict["x"] * gc_window + gc_window // 2
                                blob_dict = {"x": chunk_dict["x"].tolist(), "y": chunk_dict["y"].tolist()}
                            else:
                                blob_dict = {"x": (chunk_dict["x"] + 1).tolist(), "y": chunk_dict["y"].tolist()}
                    else:
                        blob_dict = None
                if blob_dict is not None:
                    _resolve_partner_lane(cur, blob_dict, contig_feature_name)
                    feature_dict.update(blob_dict)
                    feature_dict["is_relative_scaled"] = False
                    feature_dict["has_stats"] = False
                    feature_dict["has_sequences"] = False
                    list_feature_dict.append(feature_dict)
            continue

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
                # Decide whether to fetch lightweight zoom-only blob or full data
                _threshold = max_base_resolution if max_base_resolution is not None else _DEFAULT_MAX_BASE_RESOLUTION
                _window = (xend - xstart) if (xstart is not None and xend is not None) else 0
                _use_zoom = _window > _threshold
                _use_chunks = not _use_zoom and xstart is not None and xend is not None

                if _use_chunks:
                    # Chunk path per sample (fast: only 1-2 small rows each)
                    from thebigbam.database.blob_decoder import decode_raw_chunks, decode_raw_sparse_chunks, get_scale_from_zoom_blob, is_sparse_zoom_blob
                    zoom_blob = _get_feature_blob(cur, contig_id, sample_ids[0], feature_var_name)
                    scale_div = get_scale_from_zoom_blob(zoom_blob) if zoom_blob else 1
                    _sparse = zoom_blob and is_sparse_zoom_blob(zoom_blob)
                    for sid in sample_ids:
                        chunk_rows = _get_feature_chunks(cur, contig_id, sid, feature_var_name, xstart - 1, xend - 1)
                        if chunk_rows:
                            if _sparse:
                                blob_dict = _format_chunks_for_bokeh(decode_raw_sparse_chunks(chunk_rows, scale_div), xstart, xend, type_picked)
                            else:
                                blob_dict = _format_chunks_for_bokeh(decode_raw_chunks(chunk_rows, scale_div), xstart, xend, type_picked)
                        else:
                            blob_dict = None
                        if blob_dict is not None:
                            if min_relative_value > 0.0 and "y" in blob_dict and blob_dict["y"]:
                                max_y = max(blob_dict["y"])
                                if max_y > 0:
                                    threshold = min_relative_value * max_y
                                    blob_dict["y"] = [v if v >= threshold else 0 for v in blob_dict["y"]]
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

                # Single query for all samples instead of N individual queries
                blobs_by_sample = _get_feature_blobs_batch(cur, contig_id, sample_ids, feature_var_name)
                for sid, blob_bytes in blobs_by_sample.items():
                    if _use_zoom:
                        blob_dict = _blob_to_feature_dict(None, type_picked, xstart, xend, max_base_resolution, zoom_blob_bytes=blob_bytes)
                    else:
                        blob_dict = _blob_to_feature_dict(blob_bytes, type_picked, xstart, xend, max_base_resolution)
                    if blob_dict is not None:
                        # Filter sparse positions: zero out y values below min_relative_value * max(y)
                        if min_relative_value > 0.0 and "y" in blob_dict and blob_dict["y"]:
                            max_y = max(blob_dict["y"])
                            if max_y > 0:
                                threshold = min_relative_value * max_y
                                blob_dict["y"] = [v if v >= threshold else 0 for v in blob_dict["y"]]
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
                # GC content/skew base resolution is already windowed (500/1000 bp),
                # so use base resolution up to 10 Mbp window
                is_gc = contig_feature_name in ("gc_content", "gc_skew")
                res_threshold = 10_000_000 if is_gc else max_base_resolution
                _window = (xend - xstart) if (xstart is not None and xend is not None) else 0
                _use_zoom = _window > (res_threshold if res_threshold is not None else _DEFAULT_MAX_BASE_RESOLUTION)
                if _use_zoom:
                    zoom_blob = _get_contig_blob(cur, contig_id, contig_feature_name)
                    blob_dict = _blob_to_feature_dict(None, type_picked, xstart, xend, res_threshold, zoom_blob_bytes=zoom_blob) if zoom_blob else None
                elif xstart is not None and xend is not None:
                    # Zoomed in: fetch chunks
                    from thebigbam.database.blob_decoder import decode_raw_chunks, decode_raw_sparse_chunks, get_scale_from_zoom_blob, is_sparse_zoom_blob, gc_window_size
                    gc_window = gc_window_size(contig_feature_name)
                    chunk_rows = _get_contig_chunks(cur, contig_id, contig_feature_name, xstart - 1, xend - 1, window_size=gc_window)
                    zoom_blob = _get_contig_blob(cur, contig_id, contig_feature_name)
                    scale_div = get_scale_from_zoom_blob(zoom_blob) if zoom_blob else 1
                    is_sparse = bool(zoom_blob and is_sparse_zoom_blob(zoom_blob))
                    if chunk_rows:
                        if is_sparse:
                            blob_dict = _format_chunks_for_bokeh(decode_raw_sparse_chunks(chunk_rows, scale_div), xstart, xend, type_picked)
                        else:
                            chunk_dict = decode_raw_chunks(chunk_rows, scale_div)
                            if gc_window > 1:
                                chunk_dict["x"] = chunk_dict["x"] * gc_window
                            if contig_feature_name == "gc_skew":
                                chunk_dict["y"] = chunk_dict["y"] / 100.0
                            blob_dict = _format_chunks_for_bokeh(chunk_dict, xstart, xend, type_picked)
                    else:
                        blob_dict = None
                else:
                    # No window specified — fetch all chunks
                    from thebigbam.database.blob_decoder import decode_raw_chunks, decode_raw_sparse_chunks, get_scale_from_zoom_blob, is_sparse_zoom_blob, gc_window_size
                    gc_window = gc_window_size(contig_feature_name)
                    chunk_rows = _get_contig_chunks(cur, contig_id, contig_feature_name, 0, 2**31, window_size=gc_window)
                    zoom_blob = _get_contig_blob(cur, contig_id, contig_feature_name)
                    scale_div = get_scale_from_zoom_blob(zoom_blob) if zoom_blob else 1
                    is_sparse = bool(zoom_blob and is_sparse_zoom_blob(zoom_blob))
                    if chunk_rows:
                        if is_sparse:
                            raw = decode_raw_sparse_chunks(chunk_rows, scale_div)
                            blob_dict = {"x": (raw["x"] + 1).tolist(), "y": raw["y"].tolist()}
                            for k in raw:
                                if k not in ("x", "y", "sparse"):
                                    blob_dict[k] = raw[k].tolist() if hasattr(raw[k], "tolist") else raw[k]
                        else:
                            chunk_dict = decode_raw_chunks(chunk_rows, scale_div)
                            if gc_window > 1:
                                chunk_dict["x"] = chunk_dict["x"] * gc_window
                            blob_dict = {"x": (chunk_dict["x"] + 1).tolist(), "y": chunk_dict["y"].tolist()}
                    else:
                        blob_dict = None
                if blob_dict is not None:
                    _resolve_partner_lane(cur, blob_dict, contig_feature_name)
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
def generate_bokeh_plot_per_sample(conn, list_features, contig_name, sample_name, xstart=None, xend=None, subplot_size=100, genbank_path=None, feature_types=None, plot_isoforms=True, plot_sequence=False, plot_translated_sequence=False, same_y_scale=False, genemap_size=None, sequence_size=None, translated_sequence_size=None, max_base_resolution=None, max_genemap_window=None, max_sequence_window=None, min_relative_value=0.0, feature_label_key=None, custom_colors=None, mag_name=None):
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

    _genemap_threshold = max_genemap_window if max_genemap_window is not None else DEFAULT_GENEMAP_WINDOW
    annotation_fig = None
    if genbank_path and xstart is not None and xend is not None and (xend - xstart) <= _genemap_threshold:
        annotation_fig = make_bokeh_genemap(
            conn, contig_id, locus_name, locus_size,
            genemap_size if genemap_size is not None else subplot_size,
            shared_xrange, xstart, xend,
            feature_types=feature_types, plot_isoforms=plot_isoforms,
            feature_label_key=feature_label_key, custom_colors=custom_colors
        )

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
    _seq_threshold = max_sequence_window if max_sequence_window is not None else DEFAULT_SEQUENCE_WINDOW
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

    # Separate contig-level features from sample-dependent features by
    # consulting the Variable table: any subplot whose rows all live in
    # Contig_blob is plotted once per contig, irrespective of samples.
    metadata_cache = get_variable_metadata_batch(cur, requested_features)
    contig_features, sample_features = split_contig_vs_sample_features(metadata_cache, requested_features)

    # Add contig-level features (don't require sample_id)
    if contig_features:
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

    # --- Optional MAG track (only when a MAG is selected in MAG-mode DBs) ---
    mag_fig = None
    if mag_name:
        try:
            mag_fig = make_bokeh_mag_track(conn, mag_name, height=30)
        except Exception as e:
            print(f"Error building MAG track for '{mag_name}': {e}", flush=True)
            mag_fig = None

    # --- Combine all figures in a single grid with one shared toolbar ---
    top_plots = [mag_fig] if mag_fig is not None else []
    if annotation_fig:
        if not subplots:
            grid = gridplot([[p] for p in (top_plots + [annotation_fig])], merge_tools=True, sizing_mode='stretch_width')
        else:
            all_plots = top_plots + [annotation_fig] + subplots
            grid = gridplot([[p] for p in all_plots], merge_tools=True, sizing_mode='stretch_width')
    else:
        if not subplots and not top_plots:
            raise ValueError("No plots to display")
        grid = gridplot([[p] for p in (top_plots + subplots)], merge_tools=True, sizing_mode='stretch_width')

    return grid


def _merge_decoded_chunks(decoded_list):
    """Merge decoded chunk dicts from multiple contigs into one sorted result."""
    import numpy as np
    all_x = np.concatenate([d["x"] for d in decoded_list])
    all_y = np.concatenate([d["y"] for d in decoded_list])
    is_sparse = any(d.get("sparse", False) for d in decoded_list)

    order = np.argsort(all_x)
    result = {"x": all_x[order], "y": all_y[order], "sparse": is_sparse}

    meta_keys = set()
    for d in decoded_list:
        meta_keys.update(k for k in d if k not in ("x", "y", "sparse"))
    for key in meta_keys:
        arrays = []
        for d in decoded_list:
            if key in d:
                val = d[key]
                if isinstance(val, np.ndarray):
                    arrays.append(val)
                elif isinstance(val, list):
                    arrays.append(np.array(val))
        if arrays:
            merged = np.concatenate(arrays)
            result[key] = merged[order]

    return result


def _assemble_mag_chunks_from_contigs(cur, mag_id, sample_id, variable_name,
                                       zoom_blob, xs, xe, type_picked,
                                       is_contig_blob=False):
    """Build MAG base-resolution data on the fly from per-contig chunks.

    Instead of reading from (now-removed) MAG_blob_chunk / MAG_contig_blob_chunk
    tables, fetch the relevant per-contig chunks and shift positions by the
    contig's offset in the MAG.
    """
    import numpy as np
    from thebigbam.database.database_getters import get_mag_members
    from thebigbam.database.blob_decoder import (
        decode_raw_chunks, decode_raw_sparse_chunks,
        get_scale_from_zoom_blob, is_sparse_zoom_blob,
    )

    members = get_mag_members(cur, mag_id)
    scale_div = get_scale_from_zoom_blob(zoom_blob)
    sparse = is_sparse_zoom_blob(zoom_blob)

    all_decoded = []
    for cid, clen, offset in members:
        # Contig covers [offset+1, offset+clen] in 1-indexed MAG space
        if offset + clen < xs or offset + 1 > xe:
            continue
        # Convert MAG window to contig-local 0-indexed coords
        local_start = max(0, xs - 1 - offset)
        local_end = min(clen, xe - 1 - offset)

        if is_contig_blob:
            chunk_rows = _get_contig_chunks(cur, cid, variable_name, local_start, local_end)
        else:
            chunk_rows = _get_feature_chunks(cur, cid, sample_id, variable_name, local_start, local_end)

        if not chunk_rows:
            continue

        if sparse:
            decoded = decode_raw_sparse_chunks(chunk_rows, scale_div)
        else:
            decoded = decode_raw_chunks(chunk_rows, scale_div)

        if len(decoded["x"]) > 0:
            decoded["x"] = decoded["x"] + offset
            all_decoded.append(decoded)

    if all_decoded:
        merged = _merge_decoded_chunks(all_decoded)
        return _format_chunks_for_bokeh(merged, xs, xe, type_picked)
    return None


def get_mag_feature_data(cur, feature, mag_id, sample_id, mag_length,
                         xstart=None, xend=None, variable_metadata=None,
                         max_base_resolution=None, min_relative_value=0.0):
    """Read one MAG-wide blob per subplot from MAG_blob / MAG_contig_blob.

    Shape of the returned list_feature_dict is identical to
    get_feature_data(): `make_bokeh_subplot` consumes both the same way.
    X-coordinates are already in MAG space — no offset arithmetic required.
    """
    from thebigbam.database.database_getters import (
        get_mag_feature_zoom, get_mag_contig_zoom, get_mag_members,
    )
    from thebigbam.database.blob_decoder import (
        decode_raw_chunks, decode_raw_sparse_chunks,
        get_scale_from_zoom_blob, is_sparse_zoom_blob, CHUNK_SIZE,
    )

    rows = variable_metadata if variable_metadata is not None else get_variable_metadata(cur, feature)
    xs = xstart if xstart is not None else 1
    xe = xend   if xend   is not None else mag_length
    _threshold = max_base_resolution if max_base_resolution is not None else _DEFAULT_MAX_BASE_RESOLUTION
    _window = xe - xs
    _use_zoom = _window > _threshold

    list_feature_dict = []
    for row in rows:
        type_picked, color, alpha, fill_alpha, size, title, feature_table = row
        feature_dict = {
            "type": type_picked, "color": color, "alpha": alpha,
            "fill_alpha": fill_alpha, "size": size, "title": title,
            "x": [], "y": [], "is_relative_scaled": False,
        }

        cur.execute(
            "SELECT Variable_name FROM Variable WHERE Title=? AND Feature_table_name=?",
            (title, feature_table),
        )
        vname_row = cur.fetchone()
        variable_name = vname_row[0] if vname_row else None
        if not variable_name:
            continue

        if feature_table == "Feature_blob":
            zoom_blob = get_mag_feature_zoom(cur, mag_id, sample_id, variable_name)
            if zoom_blob is None:
                continue
            if _use_zoom:
                blob_dict = _blob_to_feature_dict(None, type_picked, xs, xe, max_base_resolution, zoom_blob_bytes=zoom_blob)
            else:
                blob_dict = _assemble_mag_chunks_from_contigs(
                    cur, mag_id, sample_id, variable_name, zoom_blob,
                    xs, xe, type_picked, is_contig_blob=False,
                )
        elif feature_table == "Contig_blob":
            zoom_blob = get_mag_contig_zoom(cur, mag_id, variable_name)
            if zoom_blob is None:
                continue
            if _use_zoom:
                blob_dict = _blob_to_feature_dict(None, type_picked, xs, xe, max_base_resolution, zoom_blob_bytes=zoom_blob)
            else:
                blob_dict = _assemble_mag_chunks_from_contigs(
                    cur, mag_id, None, variable_name, zoom_blob,
                    xs, xe, type_picked, is_contig_blob=True,
                )
        else:
            continue

        if blob_dict is None:
            continue

        if min_relative_value > 0.0 and blob_dict.get("y"):
            max_y = max(blob_dict["y"])
            if max_y > 0:
                cutoff = min_relative_value * max_y
                blob_dict["y"] = [v if v >= cutoff else 0 for v in blob_dict["y"]]

        _resolve_partner_lane(cur, blob_dict, variable_name)

        feature_dict.update(blob_dict)
        feature_dict["is_relative_scaled"] = False
        feature_dict["has_stats"] = "mean" in blob_dict
        feature_dict["has_sequences"] = "sequence" in blob_dict
        list_feature_dict.append(feature_dict)

    return list_feature_dict


def make_bokeh_genemap_mag(conn, mag_id, mag_name, mag_length, subplot_size,
                           shared_xrange, xstart=None, xend=None,
                           feature_types=None, plot_isoforms=True,
                           feature_label_key=None, custom_colors=None,
                           figure_width=30):
    """MAG-wide gene map. Builds a single SeqRecord spanning the entire MAG
    using MAG_annotation_core (coordinates already offset into MAG space),
    then runs it through the same DNAFeaturesViewer path as
    make_bokeh_genemap. This preserves the y_range DNAFeaturesViewer sets,
    so features render at identical thickness to the Contig view.
    """
    cur = conn.cursor()

    position_filter = ""
    params = [mag_id]
    if xstart is not None and xend is not None:
        position_filter = ' AND mac."End" >= ? AND mac."Start" <= ?'
        params.extend([xstart, xend])

    type_filter = ""
    if feature_types:
        type_filter = f' AND mac."Type" IN ({",".join("?" * len(feature_types))})'
        params.extend(feature_types)

    base_select = (
        'SELECT mac.Annotation_id, mac."Start", mac."End", mac.Strand, mac."Type", '
        'pq.Value AS Product, fq.Value AS Function, lq.Value AS Locus_tag, '
        'mac.Main_isoform '
        'FROM MAG_annotation_core mac '
        'LEFT JOIN Annotation_qualifier pq ON pq.Annotation_id = mac.Annotation_id AND pq."Key" = \'product\' '
        'LEFT JOIN Annotation_qualifier fq ON fq.Annotation_id = mac.Annotation_id AND fq."Key" = \'function\' '
        'LEFT JOIN Annotation_qualifier lq ON lq.Annotation_id = mac.Annotation_id AND lq."Key" = \'locus_tag\' '
        'WHERE mac.MAG_id = ?'
    )
    if not plot_isoforms:
        query = f'{base_select}{position_filter}{type_filter} AND (lq.Value IS NULL OR mac.Main_isoform = true)'
    else:
        query = f'{base_select}{position_filter}{type_filter}'

    cur.execute(query, tuple(params))
    seq_ann_rows = cur.fetchall()

    ann_ids = [r[0] for r in seq_ann_rows]
    placeholders = ','.join('?' * len(ann_ids)) if ann_ids else ''
    direct_col_idx = {'Type': 4, 'type': 4, 'Strand': 3, 'strand': 3,
                      'Main_isoform': 8, 'main_isoform': 8}

    def fetch_values_for_qkey(qkey):
        if not ann_ids:
            return {}
        try:
            kv_rows = cur.execute(
                f'SELECT Annotation_id, "Value" FROM Annotation_qualifier '
                f'WHERE "Key" = ? AND Annotation_id IN ({placeholders})',
                [qkey] + ann_ids,
            ).fetchall()
            kv = {aid: val for aid, val in kv_rows if val is not None}
            if kv:
                return kv
        except Exception:
            pass
        col_idx = direct_col_idx.get(qkey)
        if col_idx is not None:
            return {r[0]: r[col_idx] for r in seq_ann_rows if r[col_idx] is not None}
        return {}

    label_map = {}
    if feature_label_key and seq_ann_rows:
        label_map = fetch_values_for_qkey(feature_label_key)

    custom_color_map = {}
    if custom_colors and seq_ann_rows:
        from collections import defaultdict
        rules_by_key = defaultdict(list)
        for rule in custom_colors:
            rules_by_key[rule['qualifier_key']].append(rule)
        for qkey, rules in rules_by_key.items():
            values = fetch_values_for_qkey(qkey)
            for rule in rules:
                mode = rule.get('match_mode', 'exact')
                if mode == 'random':
                    filtered = {aid: val for aid, val in values.items() if not _is_nullish(val)}
                    numeric_map = _try_numeric_map(filtered.values())
                    if numeric_map is not None:
                        value_to_color = _gradient_palette(numeric_map)
                    else:
                        distinct = {str(v) for v in filtered.values()}
                        value_to_color = {v: _hash_color(v) for v in distinct}
                    for aid, val in filtered.items():
                        if aid in custom_color_map:
                            continue
                        color = value_to_color.get(str(val))
                        if color:
                            custom_color_map[aid] = color
                    continue
                rule_val = rule['value']
                if mode in ('lt', 'gt'):
                    try:
                        rule_num = float(rule_val)
                    except (TypeError, ValueError):
                        continue
                    for aid, val in values.items():
                        if aid in custom_color_map or _is_nullish(val):
                            continue
                        try:
                            v_num = float(val)
                        except (TypeError, ValueError):
                            continue
                        if mode == 'lt' and v_num < rule_num:
                            custom_color_map[aid] = rule['color']
                        elif mode == 'gt' and v_num > rule_num:
                            custom_color_map[aid] = rule['color']
                    continue
                rule_val_str = str(rule_val)
                rule_val_lower = rule_val_str.lower()
                rule_val_num = float(rule_val) if isinstance(rule_val, (int, float)) else None

                def _values_equal(val):
                    if rule_val_num is not None:
                        try:
                            return float(val) == rule_val_num
                        except (TypeError, ValueError):
                            pass
                    return str(val) == rule_val_str

                is_negation = mode in ('has_not', 'not_equal')
                for aid, val in values.items():
                    if aid in custom_color_map or _is_nullish(val):
                        continue
                    val_str = str(val)
                    val_lower = val_str.lower()
                    if mode == 'has':
                        if rule_val_lower in val_lower:
                            custom_color_map[aid] = rule['color']
                    elif mode == 'has_not':
                        if rule_val_lower not in val_lower:
                            custom_color_map[aid] = rule['color']
                    elif mode == 'not_equal':
                        if not _values_equal(val):
                            custom_color_map[aid] = rule['color']
                    else:
                        if _values_equal(val):
                            custom_color_map[aid] = rule['color']

                if is_negation:
                    for aid in ann_ids:
                        if aid not in custom_color_map and aid not in values:
                            custom_color_map[aid] = rule['color']

    sequence_annotations = []
    for ann_id, start, end, strand, ftype, product, function, locus_tag, _main in seq_ann_rows:
        try:
            floc = FeatureLocation(start - 1, end, strand=strand)
        except Exception:
            continue
        qualifiers = {}
        if product:
            qualifiers['product'] = product
        if function:
            qualifiers['function'] = function
        if locus_tag:
            qualifiers['locus_tag'] = locus_tag
        if custom_color_map and ann_id in custom_color_map:
            qualifiers['_custom_color'] = custom_color_map[ann_id]
        if feature_label_key:
            qualifiers['_tooltip_key'] = feature_label_key
            if ann_id in label_map:
                qualifiers[feature_label_key] = label_map[ann_id]
        sequence_annotations.append(SeqFeature(location=floc, type=ftype, qualifiers=qualifiers))

    if not sequence_annotations:
        return None

    sequence_records = SeqRecord(Seq('N' * mag_length), id=mag_name, features=sequence_annotations)
    graphic_record = CustomTranslator().translate_record(sequence_records)
    annotation_fig = graphic_record.plot_with_bokeh(figure_width=figure_width, figure_height=40)
    annotation_fig.height = subplot_size
    annotation_fig.tools = [t for t in annotation_fig.tools if not isinstance(t, TapTool)]
    annotation_fig.xaxis.formatter = NumeralTickFormatter(format="0,0")
    wheel = WheelZoomTool(dimensions='width')
    annotation_fig.add_tools(wheel)
    annotation_fig.toolbar.active_scroll = wheel
    annotation_fig.x_range = shared_xrange
    return annotation_fig


def generate_bokeh_plot_mag_view(conn, list_features, mag_name, sample_name, xstart=None, xend=None, subplot_size=100, genbank_path=None, feature_types=None, plot_isoforms=True, plot_sequence=False, plot_translated_sequence=False, same_y_scale=False, genemap_size=None, sequence_size=None, translated_sequence_size=None, max_base_resolution=None, max_genemap_window=None, max_sequence_window=None, min_relative_value=0.0, feature_label_key=None, custom_colors=None, is_all=False, allowed_samples=None):
    """Generate a concatenated Bokeh plot for a MAG (all contigs, longest-first).

    All contigs are placed consecutively on a shared x-axis (longest first), with
    a MAG overview track at the top showing contig boundaries.
    """
    from ..database.database_getters import get_mag_contigs, get_mag_id

    cur = conn.cursor()
    members = get_mag_contigs(conn, mag_name)
    if not members:
        raise ValueError(f"MAG not found or empty: {mag_name}")

    mag_id = get_mag_id(conn, mag_name)
    if mag_id is None:
        raise ValueError(f"MAG_id lookup failed for: {mag_name}")

    total_len = sum(length for _n, length, _o in members)
    print(f"MAG {mag_name}: {len(members)} contigs, {total_len} bp total", flush=True)

    shared_xrange = Range1d(
        xstart if xstart is not None else 1,
        xend if xend is not None else total_len,
    )

    # --- MAG track at top (same height as nucleotide sequence track) ---
    _seq_height = sequence_size if sequence_size is not None else subplot_size // 2
    try:
        mag_fig = make_bokeh_mag_track(conn, mag_name, height=_seq_height, shared_xrange=shared_xrange)
    except Exception as e:
        print(f"Error building MAG track for '{mag_name}': {e}", flush=True)
        mag_fig = None

    # --- Combined gene map (only when window is small enough) ---
    _genemap_threshold = max_genemap_window if max_genemap_window is not None else DEFAULT_GENEMAP_WINDOW
    annotation_fig = None
    if genbank_path and xstart is not None and xend is not None and (xend - xstart) <= _genemap_threshold:
        _gm_size = genemap_size if genemap_size is not None else subplot_size
        try:
            annotation_fig = make_bokeh_genemap_mag(
                conn, mag_id, mag_name, total_len, _gm_size, shared_xrange, xstart, xend,
                feature_types=feature_types, plot_isoforms=plot_isoforms,
                feature_label_key=feature_label_key, custom_colors=custom_colors,
            )
        except Exception as e:
            print(f"Error building combined gene map for MAG '{mag_name}': {e}", flush=True)

    # --- Resolve sample ---
    sample_id = None
    if sample_name:
        cur.execute("SELECT Sample_id, Sample_name FROM Sample WHERE Sample_name=?", (sample_name,))
        row = cur.fetchone()
        if row is None:
            raise ValueError(f"Sample not found: {sample_name}")
        sample_id, sample_name = row
        print(f"Sample {sample_name} validated.", flush=True)

    requested_features = parse_requested_features(list_features)
    metadata_cache = get_variable_metadata_batch(cur, requested_features)
    contig_features, sample_features = split_contig_vs_sample_features(metadata_cache, requested_features)

    subplots = []

    # --- Sequence subplot (stitched across member contigs) ---
    _seq_threshold = max_sequence_window if max_sequence_window is not None else DEFAULT_SEQUENCE_WINDOW
    if plot_sequence and xstart is not None and xend is not None and (xend - xstart) <= _seq_threshold:
        seq_subplot = make_bokeh_sequence_subplot_mag(conn, mag_name, xstart, xend, _seq_height, shared_xrange)
        if seq_subplot:
            subplots.append(seq_subplot)
    elif plot_sequence and xstart is not None and xend is not None and (xend - xstart) > _seq_threshold:
        print(f"Sequence not plotted: window > {_seq_threshold} bp", flush=True)

    # --- Translated-sequence subplot (CDS across member contigs) ---
    if plot_translated_sequence and xstart is not None and xend is not None and (xend - xstart) <= _seq_threshold:
        _trans_height = translated_sequence_size if translated_sequence_size is not None else (sequence_size if sequence_size is not None else subplot_size // 2)
        trans_subplot = make_bokeh_translated_sequence_subplot_mag(conn, mag_name, xstart, xend, _trans_height, shared_xrange)
        if trans_subplot:
            subplots.append(trans_subplot)
    elif plot_translated_sequence and xstart is not None and xend is not None and (xend - xstart) > _seq_threshold:
        print(f"Translated sequence not plotted: window > {_seq_threshold} bp", flush=True)

    # --- Contig-level features (GC content, GC skew, repeats, hits, …) ---
    if contig_features:
        for feature in contig_features:
            try:
                combined_dicts = get_mag_feature_data(
                    cur, feature, mag_id, None, total_len,
                    xstart=xstart, xend=xend,
                    variable_metadata=metadata_cache.get(feature),
                    max_base_resolution=max_base_resolution,
                    min_relative_value=min_relative_value,
                )
                subplot = make_bokeh_subplot(combined_dicts, subplot_size, shared_xrange, show_tooltips=True)
                if subplot is not None:
                    subplots.append(subplot)
            except Exception as e:
                print(f"Error processing contig feature '{feature}' for MAG: {e}", flush=True)

    # --- Sample-dependent features ---
    if is_all and sample_features:
        # ALL SAMPLES mode: one subplot per sample, each showing the concatenated MAG
        all_rows = cur.execute(
            "SELECT DISTINCT s.Sample_id, s.Sample_name FROM Sample s "
            "JOIN MAG_blob mb ON mb.Sample_id = s.Sample_id "
            "WHERE mb.MAG_id = ?",
            (mag_id,),
        ).fetchall()
        if allowed_samples is not None:
            all_rows = [(sid, sname) for sid, sname in all_rows if sname in allowed_samples]
        all_rows.sort(key=lambda r: r[1])
        for feature in sample_features:
            for sid, sname in all_rows:
                try:
                    combined_dicts = get_mag_feature_data(
                        cur, feature, mag_id, sid, total_len,
                        xstart=xstart, xend=xend,
                        variable_metadata=metadata_cache.get(feature),
                        max_base_resolution=max_base_resolution,
                        min_relative_value=min_relative_value,
                    )
                    for d in combined_dicts:
                        d['title'] = sname
                    subplot = make_bokeh_subplot(combined_dicts, subplot_size, shared_xrange, show_tooltips=True)
                    if subplot is not None:
                        subplots.append(subplot)
                except Exception as e:
                    print(f"Error processing feature '{feature}' for sample '{sname}' in MAG: {e}", flush=True)
    elif sample_id is not None and sample_features:
        for feature in sample_features:
            try:
                combined_dicts = get_mag_feature_data(
                    cur, feature, mag_id, sample_id, total_len,
                    xstart=xstart, xend=xend,
                    variable_metadata=metadata_cache.get(feature),
                    max_base_resolution=max_base_resolution,
                    min_relative_value=min_relative_value,
                )
                subplot = make_bokeh_subplot(combined_dicts, subplot_size, shared_xrange, show_tooltips=True)
                if subplot is not None:
                    subplots.append(subplot)
            except Exception as e:
                print(f"Error processing sample feature '{feature}' for MAG: {e}", flush=True)

    # --- Assemble grid ---
    top_plots = [mag_fig] if mag_fig is not None else []
    all_plots = top_plots + ([annotation_fig] if annotation_fig is not None else []) + subplots

    if not all_plots:
        raise ValueError("No plots to display for MAG view")

    return gridplot([[p] for p in all_plots], merge_tools=True, sizing_mode='stretch_width')
