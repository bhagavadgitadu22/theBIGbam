import duckdb
from bokeh.models import Range1d
from bokeh.layouts import gridplot

from .plotting_data_per_sample import get_contig_info, get_feature_data, get_feature_data_batch, get_variable_metadata, get_repeats_data, make_bokeh_subplot, make_bokeh_genemap

### Function to generate the bokeh plot
def generate_bokeh_plot_all_samples(conn, variable, contig_name, xstart=None, xend=None, subplot_size=130, genbank_path=None, genome_features=None, allowed_samples=None, feature_types=None, use_phage_colors=False):
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
        use_phage_colors: Whether to use phage color scheme for CDS features
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

    annotation_fig = make_bokeh_genemap(conn, contig_id, locus_name, locus_size, subplot_size, shared_xrange, xstart, xend, feature_types=feature_types, use_phage_colors=use_phage_colors) if genbank_path else None

    # Get list of samples
    cur.execute("SELECT Presences.Sample_id, Sample_name FROM Presences JOIN Sample ON Presences.Sample_id = Sample.Sample_id WHERE Contig_id=?", (contig_id,))
    rows = cur.fetchall()
    if rows is None:
        raise ValueError(f"No sample comprised this contig in the database: {contig_name}")

    # Filter to allowed samples if specified (respects Filtering section criteria)
    if allowed_samples is not None:
        rows = [(sid, sname) for sid, sname in rows if sname in allowed_samples]
        if not rows:
            raise ValueError("No samples match the current filters")

    sample_ids, sample_names = [list(t) for t in zip(*rows)]

    # --- Add subplots for additional Genome features (contig-level, not per-sample) ---
    genome_subplots = []
    if genome_features:
        for genome_feature in genome_features:
            try:
                feature_lower = genome_feature.lower().strip()

                # Handle Repeats specially - they use contig-only tables with linked positions
                if feature_lower in ["repeats", "repeat", "direct repeats", "inverted repeats"]:
                    repeats_feature_dicts = []
                    if feature_lower in ["repeats", "repeat", "direct repeats"]:
                        direct_feature_dict = get_repeats_data(cur, contig_id, "direct_repeats", xstart, xend)
                        if direct_feature_dict:
                            repeats_feature_dicts.extend(direct_feature_dict)
                    if feature_lower in ["repeats", "repeat", "inverted repeats"]:
                        inverted_feature_dict = get_repeats_data(cur, contig_id, "inverted_repeats", xstart, xend)
                        if inverted_feature_dict:
                            repeats_feature_dicts.extend(inverted_feature_dict)
                    if repeats_feature_dicts:
                        repeats_subplot = make_bokeh_subplot(repeats_feature_dicts, subplot_size, shared_xrange)
                        if repeats_subplot is not None:
                            genome_subplots.append(repeats_subplot)
                # Handle GC content and GC skew - contig-level tables, use get_feature_data
                elif feature_lower in ["gc_content", "gc content", "gccontent", "gc"]:
                    list_feature_dict = get_feature_data(cur, "gc_content", contig_id, sample_id=None, xstart=xstart, xend=xend)
                    if list_feature_dict:
                        gc_subplot = make_bokeh_subplot(list_feature_dict, subplot_size, shared_xrange)
                        if gc_subplot is not None:
                            genome_subplots.append(gc_subplot)
                elif feature_lower in ["gc_skew", "gc skew", "gcskew", "skew"]:
                    list_feature_dict = get_feature_data(cur, "gc_skew", contig_id, sample_id=None, xstart=xstart, xend=xend)
                    if list_feature_dict:
                        gc_skew_subplot = make_bokeh_subplot(list_feature_dict, subplot_size, shared_xrange)
                        if gc_skew_subplot is not None:
                            genome_subplots.append(gc_skew_subplot)
                else:
                    # Other Genome features - try to get data (may fail if sample-dependent)
                    list_feature_dict = get_feature_data(cur, genome_feature, contig_id, sample_id=None, xstart=xstart, xend=xend)
                    if not list_feature_dict:
                        continue

                    subplot_feature = make_bokeh_subplot(list_feature_dict, subplot_size, shared_xrange, sample_title=genome_feature)
                    if subplot_feature is not None:
                        genome_subplots.append(subplot_feature)
            except Exception as e:
                print(f"Error processing genome feature '{genome_feature}': {e}", flush=True)
                continue

    # --- Add one subplot per sample for the main variable ---
    # Requested features are variables like 'coverage', 'reads_starts', etc.
    # Batch-fetch data for all samples in one query, then create subplots
    subplots = []
    try:
        var_metadata = get_variable_metadata(cur, variable)
        batch_results = get_feature_data_batch(cur, variable, contig_id, sample_ids, xstart, xend, variable_metadata=var_metadata)
        for sample_id, sample_name in zip(sample_ids, sample_names):
            list_feature_dict = batch_results.get(sample_id, [])
            if not list_feature_dict:
                continue
            subplot_feature = make_bokeh_subplot(list_feature_dict, subplot_size, shared_xrange, sample_title=sample_name)
            if subplot_feature is not None:
                subplots.append(subplot_feature)
    except Exception as e:
        print(f"Error batch-processing variable '{variable}': {e}", flush=True)

    # --- Combine all figures in a single grid with one shared toolbar ---
    all_plots = []
    if annotation_fig:
        all_plots.append(annotation_fig)
    all_plots.extend(genome_subplots)
    all_plots.extend(subplots)

    if not all_plots:
        raise ValueError("No plots to display")

    grid = gridplot([[p] for p in all_plots], merge_tools=True, sizing_mode='stretch_width')

    return grid