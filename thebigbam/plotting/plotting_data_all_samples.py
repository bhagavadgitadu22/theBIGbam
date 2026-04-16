import time

from bokeh.models import Range1d
from bokeh.layouts import gridplot

from .plotting_data_per_sample import get_contig_info, get_feature_data, get_feature_data_batch, get_variable_metadata, make_bokeh_subplot, make_bokeh_genemap, make_bokeh_sequence_subplot, make_bokeh_translated_sequence_subplot

### Function to generate the bokeh plot
def generate_bokeh_plot_all_samples(conn, variable, contig_name, xstart=None, xend=None, subplot_size=130, genbank_path=None, genome_features=None, allowed_samples=None, feature_types=None, plot_isoforms=True, plot_sequence=False, plot_translated_sequence=False, same_y_scale=False, genemap_size=None, sequence_size=None, translated_sequence_size=None, order_by_column=None, max_base_resolution=None, max_genemap_window=None, min_relative_value=0.0, feature_label_key=None, custom_colors=None):
    """Generate a Bokeh plot showing all samples for a single variable.

    Args:
        conn: DuckDB connection
        variable: Variable/feature to plot (from non-Genome modules)
        contig_name: Name of the contig to plot
        xstart: Optional x-axis start position
        xend: Optional x-axis end position
        subplot_size: Height of each subplot in pixels
        genbank_path: Path to genbank file (optional; if provided, gene map will be plotted)
        genome_features: List of additional Genome module features to plot (optional)
        allowed_samples: Set of sample names to include (optional; if None, all samples are included)
        feature_types: Optional list of feature types to include in gene map (None = all)
        order_by_column: Optional Sample table column name to order samples by (alphabetically)
        max_base_resolution: Window size (bp) below which base resolution is used
    """
    cur = conn.cursor()

    # Get contig characteristics
    contig_id, locus_name, locus_size = get_contig_info(cur, contig_name)
    print(f"Locus {locus_name} validated ({locus_size} bp)", flush=True)

    # --- Main gene annotation plot (only if genbank provided) ---
    shared_xrange = Range1d(0, locus_size)
    if xstart is not None and xend is not None:
        shared_xrange.start = xstart
        shared_xrange.end = xend

    _genemap_threshold = max_genemap_window if max_genemap_window is not None else 100_000
    annotation_fig = None
    if genbank_path and xstart is not None and xend is not None and (xend - xstart) <= _genemap_threshold:
        annotation_fig = make_bokeh_genemap(conn, contig_id, locus_name, locus_size, genemap_size if genemap_size is not None else subplot_size, shared_xrange, xstart, xend, feature_types=feature_types, plot_isoforms=plot_isoforms, feature_label_key=feature_label_key, custom_colors=custom_colors)

    # Get list of samples
    cur.execute("SELECT Coverage.Sample_id, Sample_name FROM Coverage JOIN Sample ON Coverage.Sample_id = Sample.Sample_id WHERE Contig_id=?", (contig_id,))
    rows = cur.fetchall()
    if rows is None:
        raise ValueError(f"No sample comprised this contig in the database: {contig_name}")

    # Filter to allowed samples if specified (respects Filtering section criteria)
    if allowed_samples is not None:
        rows = [(sid, sname) for sid, sname in rows if sname in allowed_samples]
        if not rows:
            raise ValueError("No samples match the current filters")

    sample_ids, sample_names = [list(t) for t in zip(*rows)]

    # Order samples by the specified column if requested
    if order_by_column:
        try:
            # Query Sample table for the ordering column values
            sample_id_placeholders = ",".join(["?"] * len(sample_ids))
            cur.execute(
                f'SELECT Sample_id, "{order_by_column}" FROM Sample WHERE Sample_id IN ({sample_id_placeholders}) ORDER BY "{order_by_column}" ASC NULLS LAST',
                sample_ids
            )
            ordered_rows = cur.fetchall()
            # Rebuild sample_ids and sample_names in the sorted order
            ordered_ids = [r[0] for r in ordered_rows]
            # Map sample_id to sample_name
            id_to_name = {sid: sname for sid, sname in zip(sample_ids, sample_names)}
            sample_ids = ordered_ids
            sample_names = [id_to_name[sid] for sid in ordered_ids]
            print(f"Samples ordered by {order_by_column}", flush=True)
        except Exception as e:
            print(f"Warning: Could not order samples by '{order_by_column}': {e}", flush=True)

    # --- Add subplots for additional Genome features (contig-level, not per-sample) ---
    genome_subplots = []
    if genome_features:
        for genome_feature in genome_features:
            try:
                feature_lower = genome_feature.lower().strip()

                # Handle Repeat count - now uses SQL views with standard binning
                if feature_lower in ["repeat count"]:
                    list_feature_dict = get_feature_data(cur, "Repeat count", contig_id, sample_id=None, xstart=xstart, xend=xend, max_base_resolution=max_base_resolution)
                    if list_feature_dict:
                        subplot = make_bokeh_subplot(list_feature_dict, subplot_size, shared_xrange, show_tooltips=True)
                        if subplot is not None:
                            genome_subplots.append(subplot)
                # Handle Max repeat identity - now uses SQL views with standard binning
                elif feature_lower in ["max repeat identity"]:
                    list_feature_dict = get_feature_data(cur, "Max repeat identity", contig_id, sample_id=None, xstart=xstart, xend=xend, max_base_resolution=max_base_resolution)
                    if list_feature_dict:
                        subplot = make_bokeh_subplot(list_feature_dict, subplot_size, shared_xrange, show_tooltips=True)
                        if subplot is not None:
                            genome_subplots.append(subplot)
                # Handle legacy "Repeats" - add both subplots
                elif feature_lower in ["repeats", "repeat", "direct repeats", "inverted repeats"]:
                    # Track 1: Repeat count
                    count_dicts = get_feature_data(cur, "Repeat count", contig_id, sample_id=None, xstart=xstart, xend=xend, max_base_resolution=max_base_resolution)
                    if count_dicts:
                        subplot = make_bokeh_subplot(count_dicts, subplot_size, shared_xrange, show_tooltips=True)
                        if subplot is not None:
                            genome_subplots.append(subplot)
                    # Track 2: Repeat max identity
                    identity_dicts = get_feature_data(cur, "Max repeat identity", contig_id, sample_id=None, xstart=xstart, xend=xend, max_base_resolution=max_base_resolution)
                    if identity_dicts:
                        subplot = make_bokeh_subplot(identity_dicts, subplot_size, shared_xrange, show_tooltips=True)
                        if subplot is not None:
                            genome_subplots.append(subplot)
                # Handle GC content and GC skew - contig-level tables, use get_feature_data
                elif feature_lower in ["gc_content", "gc content", "gccontent", "gc"]:
                    list_feature_dict = get_feature_data(cur, "GC content", contig_id, sample_id=None, xstart=xstart, xend=xend, max_base_resolution=max_base_resolution)
                    if list_feature_dict:
                        gc_subplot = make_bokeh_subplot(list_feature_dict, subplot_size, shared_xrange, show_tooltips=True)
                        if gc_subplot is not None:
                            genome_subplots.append(gc_subplot)
                elif feature_lower in ["gc_skew", "gc skew", "gcskew", "skew"]:
                    list_feature_dict = get_feature_data(cur, "GC skew", contig_id, sample_id=None, xstart=xstart, xend=xend, max_base_resolution=max_base_resolution)
                    if list_feature_dict:
                        gc_skew_subplot = make_bokeh_subplot(list_feature_dict, subplot_size, shared_xrange, show_tooltips=True)
                        if gc_skew_subplot is not None:
                            genome_subplots.append(gc_skew_subplot)
                else:
                    # Other Genome features - try to get data (may fail if sample-dependent)
                    list_feature_dict = get_feature_data(cur, genome_feature, contig_id, sample_id=None, xstart=xstart, xend=xend, max_base_resolution=max_base_resolution)
                    if not list_feature_dict:
                        continue

                    subplot_feature = make_bokeh_subplot(list_feature_dict, subplot_size, shared_xrange, sample_title=genome_feature, show_tooltips=True)
                    if subplot_feature is not None:
                        genome_subplots.append(subplot_feature)
            except Exception as e:
                print(f"Error processing genome feature '{genome_feature}': {e}", flush=True)
                continue

    # --- Add sequence subplot right after annotation (top of genome tracks) ---
    seq_subplot_added = False
    _seq_height = sequence_size if sequence_size is not None else subplot_size // 2
    if plot_sequence:
        seq_subplot = make_bokeh_sequence_subplot(conn, contig_name, xstart, xend, _seq_height, shared_xrange)
        if seq_subplot:
            genome_subplots.insert(0, seq_subplot)
            seq_subplot_added = True

    # --- Add translated sequence subplot after DNA sequence ---
    if plot_translated_sequence:
        _trans_height = translated_sequence_size if translated_sequence_size is not None else _seq_height
        trans_subplot = make_bokeh_translated_sequence_subplot(conn, contig_name, xstart, xend, _trans_height, shared_xrange)
        if trans_subplot:
            insert_pos = 1 if seq_subplot_added else 0
            genome_subplots.insert(insert_pos, trans_subplot)

    # --- Add one subplot per sample for the main variable ---
    # Requested features are variables like 'coverage', 'reads_starts', etc.
    # Batch-fetch data for all samples in one query, then create subplots
    subplots = []
    all_max_y = 0
    all_min_y = 0
    try:
        var_metadata = get_variable_metadata(cur, variable)
        t_batch = time.perf_counter()
        batch_results = get_feature_data_batch(cur, variable, contig_id, sample_ids, xstart, xend, variable_metadata=var_metadata, max_base_resolution=max_base_resolution, min_relative_value=min_relative_value)
        print(f"[timing]   get_feature_data_batch ({len(sample_ids)} samples): {time.perf_counter() - t_batch:.3f}s", flush=True)
        t_subplots = time.perf_counter()
        for sample_id, sample_name in zip(sample_ids, sample_names):
            list_feature_dict = batch_results.get(sample_id, [])
            if not list_feature_dict:
                continue
            if same_y_scale:
                for d in list_feature_dict:
                    if d["y"]:
                        all_max_y = max(all_max_y, max(d["y"]))
                        all_min_y = min(all_min_y, min(d["y"]))
            subplot_feature = make_bokeh_subplot(list_feature_dict, subplot_size, shared_xrange, sample_title=sample_name, show_tooltips=True)
            if subplot_feature is not None:
                subplots.append(subplot_feature)
        print(f"[timing]   make_bokeh_subplot x{len(subplots)}: {time.perf_counter() - t_subplots:.3f}s", flush=True)
    except Exception as e:
        print(f"Error batch-processing variable '{variable}': {e}", flush=True)

    # Post-process y-ranges for same_y_scale (all-samples view)
    if same_y_scale and all_max_y > 0 and subplots:
        y_start = min(0, all_min_y)
        for p in subplots:
            p.y_range = Range1d(y_start, all_max_y)

    # --- Combine all figures in a single grid with one shared toolbar ---
    all_plots = []
    if annotation_fig:
        all_plots.append(annotation_fig)
    all_plots.extend(genome_subplots)
    all_plots.extend(subplots)

    if not all_plots:
        raise ValueError("No plots to display")

    t_grid = time.perf_counter()
    grid = gridplot([[p] for p in all_plots], merge_tools=True, sizing_mode='stretch_width')
    print(f"[timing]   gridplot ({len(all_plots)} figures): {time.perf_counter() - t_grid:.3f}s", flush=True)

    return grid
