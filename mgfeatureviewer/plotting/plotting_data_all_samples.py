import argparse
import duckdb
from bokeh.models import Range1d
from bokeh.layouts import gridplot
from bokeh.plotting import output_file, save

from .plotting_data_per_sample import get_contig_info, get_feature_data, get_repeats_data, make_bokeh_subplot, make_bokeh_genemap

### Function to generate the bokeh plot
def generate_bokeh_plot_all_samples(conn, variable, contig_name, xstart=None, xend=None, subplot_size=130, genbank_path=None, genome_features=None, allowed_samples=None):
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

    annotation_fig = make_bokeh_genemap(conn, contig_id, locus_name, locus_size, annotation_tool, subplot_size, shared_xrange) if genbank_path else None

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

                # Handle Repeats specially - they use contig-only tables
                if feature_lower in ["repeats", "repeat", "direct repeats", "inverted repeats"]:
                    # Direct repeats
                    if feature_lower in ["repeats", "repeat", "direct repeats"]:
                        direct_feature_dict = get_repeats_data(cur, contig_id, "direct_repeats")
                        if direct_feature_dict:
                            direct_subplot = make_bokeh_subplot(direct_feature_dict, subplot_size, shared_xrange)
                            if direct_subplot is not None:
                                genome_subplots.append(direct_subplot)
                    # Inverted repeats
                    if feature_lower in ["repeats", "repeat", "inverted repeats"]:
                        inverted_feature_dict = get_repeats_data(cur, contig_id, "inverted_repeats")
                        if inverted_feature_dict:
                            inverted_subplot = make_bokeh_subplot(inverted_feature_dict, subplot_size, shared_xrange)
                            if inverted_subplot is not None:
                                genome_subplots.append(inverted_subplot)
                else:
                    # Other Genome features - try to get data (may fail if sample-dependent)
                    list_feature_dict = get_feature_data(cur, genome_feature, contig_id, sample_id=None)
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
    subplots = []
    for sample_id, sample_name in zip(sample_ids, sample_names):
        try:
            list_feature_dict = get_feature_data(cur, variable, contig_id, sample_id)
            if not list_feature_dict:
                continue

            subplot_feature = make_bokeh_subplot(list_feature_dict, subplot_size, shared_xrange, sample_title=sample_name)
            if subplot_feature is not None:
                subplots.append(subplot_feature)
        except Exception as e:
            print(f"Error processing variable '{variable}' for sample '{sample_name}': {e}", flush=True)
            continue

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

def save_html_plot_all_samples(db_path, variable, contig_name, subplot_size, output_filename, genbank_path=None):
    # --- Save interactive HTML plot ---
    output_file(filename = output_filename)
    conn = duckdb.connect(db_path, read_only=True)
    try:
        grid = generate_bokeh_plot_all_samples(conn, variable, contig_name, subplot_size=subplot_size, genbank_path=genbank_path)
        save(grid)
    finally:
        conn.close()

### Main function
def main():
    # Parse command line arguments
    print("Parsing arguments...", flush=True)
    parser = argparse.ArgumentParser(description="Parse input files.")
    parser.add_argument("-d", "--db", required=True, help="Path to DuckDB database file")
    parser.add_argument("-v", "--variable", required=True, help="Variable to compute (only one variable allowed)")
    parser.add_argument("--contig", required=True, help="Name of the contig to plot")
    parser.add_argument("--html", required=False, default="MGFeatureViewer_all_samples.html", help="Name for output html files. A bokeh server will be started if not provided")
    parser.add_argument("--subplot_height", required=False, default=130, help="Height of each subplot (in pixels)")
    args = parser.parse_args()

    # Path parameters
    db_path = args.db
    variable = args.variable
    contig_name = args.contig
    output_filename = args.html

    # Optional plotting parameters
    subplot_size = int(args.subplot_height)

    # Reading values from database and plotting
    print(f"Saving static HTML to {output_filename}...", flush=True)
    save_html_plot_all_samples(db_path, variable, contig_name, subplot_size, output_filename)

def add_plot_all_args(parser):
    parser.add_argument("-d", "--db", required=True, help="Path to DuckDB database file")
    parser.add_argument("-v", "--variable", required=True, help="Variable to compute (only one variable allowed)")
    parser.add_argument("--contig", required=True, help="Name of the contig to plot")
    parser.add_argument("-g", "--genbank", help="Path to genbank file (optional; if provided, gene map will be plotted)")
    parser.add_argument("--html", required=False, default="MGFeatureViewer_all_samples.html", help="Name for output html files. A bokeh server will be started if not provided")
    parser.add_argument("--subplot_height", required=False, default=130, help="Height of each subplot (in pixels)")

def run_plot_all(args):
    db_path = args.db
    variable = args.variable
    contig_name = args.contig
    output_filename = args.html
    subplot_size = int(args.subplot_height)
    genbank_path = getattr(args, 'genbank', None)

    print(f"Saving static HTML to {output_filename}...", flush=True)
    save_html_plot_all_samples(db_path, variable, contig_name, subplot_size, output_filename, genbank_path)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    add_plot_all_args(parser)
    args = parser.parse_args()
    run_plot_all(args)