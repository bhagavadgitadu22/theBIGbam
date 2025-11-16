import argparse, sqlite3
from bokeh.models import Range1d
from bokeh.layouts import gridplot
from bokeh.plotting import output_file, save

from plotting_data_per_sample import get_contig_info, get_feature_data, make_bokeh_subplot, make_bokeh_genemap

### Function to generate the bokeh plot
def generate_bokeh_plot_all_samples(conn, variable, contig_name, subplot_size=130):
    cur = conn.cursor()

    # Get contig characteristics
    contig_id, locus_name, locus_size, annotation_tool = get_contig_info(cur, contig_name)
    print(f"Locus {locus_name} validated ({locus_size} bp)", flush=True)

    # --- Main gene annotation plot ---
    # Build a SeqRecord from Contig_annotation entries for this contig
    shared_xrange = Range1d(-1000, locus_size+1000)
    annotation_fig = make_bokeh_genemap(conn, contig_id, locus_name, locus_size, annotation_tool, subplot_size, shared_xrange)

    # Get list of samples
    cur.execute("SELECT Presences.Sample_id, Sample_name FROM Presences JOIN Sample ON Presences.Sample_id = Sample.Sample_id WHERE Contig_id=?", (contig_id,))
    rows = cur.fetchall()
    if rows is None:
        raise ValueError(f"No sample comprised this contig in the database: {contig_name}")
    sample_ids, sample_names = [list(t) for t in zip(*rows)]

    # --- Add one subplot per feature requested ---
    # Requested features are variables like 'coverage', 'reads_starts', etc.
    subplots = []
    for sample_id, sample_name in zip(sample_ids, sample_names):
        list_feature_dict = get_feature_data(cur, variable, contig_id, sample_id)
        if all(len(vals["x"]) == 0 for vals in list_feature_dict):
            continue

        subplot_feature = make_bokeh_subplot(list_feature_dict, subplot_size, shared_xrange, sample_title=sample_name)
        if subplot_feature is not None:
            subplots.append(subplot_feature)

    # --- Combine all figures in a single grid with one shared toolbar ---
    all_plots = [annotation_fig] + subplots
    grid = gridplot([[p] for p in all_plots], merge_tools=True, sizing_mode='stretch_width')

    return grid

def save_html_plot_all_samples(db_path, variable, contig_name, subplot_size, output_filename):
    # --- Save interactive HTML plot ---
    output_file(filename = output_filename)
    conn = sqlite3.connect(db_path)
    grid = generate_bokeh_plot_all_samples(conn, variable, contig_name, subplot_size)
    save(grid)
    conn.close()

### Main function
def main():
    # Parse command line arguments
    print("Parsing arguments...", flush=True)
    parser = argparse.ArgumentParser(description="Parse input files.")
    parser.add_argument("-d", "--db", required=True, help="Path to sqlite database file to store results")
    parser.add_argument("-v", "--variable", required=True, help="Variable to compute (only one variable allowed)")
    parser.add_argument("--contig", required=True, help="Name of the contig to plot")
    parser.add_argument("--html", required=False, default="MGFeaturesViewer_all_samples.html", help="Name for output html files. A bokeh server will be started if not provided")
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
    
if __name__ == "__main__":
    main()