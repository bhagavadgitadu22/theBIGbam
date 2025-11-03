import argparse
from pyexpat import features
import constants, plotting_features
from Bio import SeqIO
from dna_features_viewer import BiopythonTranslator
from bokeh.models import Range1d
from bokeh.layouts import column
from bokeh.plotting import output_file, save

def check_single_locus(gbk_file):
    try:
        record = SeqIO.read(gbk_file, "genbank")  # expects exactly 1 record
        return record
    except ValueError:
        raise SystemExit(f"Error: {gbk_file} must contain exactly 1 locus.")

class CustomTranslator(BiopythonTranslator):
    def compute_feature_color(self, feature):
        if ANNOTATION_TOOL == "pharokka":
            function = feature.qualifiers.get("function", [""])[0].lower()
            color_scheme = constants.PHAROKKA_COLORS
            for key, color in color_scheme.items():
                if key in function:
                    return color
        return "#AAAAAA"

    def compute_feature_label(self, feature):
        return None  # fallback to None if missing or invalid
    
    def compute_feature_html(self, feature):
        return feature.qualifiers.get("product", [])
    
def main():
    # Parse command line arguments
    print("Parsing arguments...", flush=True)
    parser = argparse.ArgumentParser(description="Parse input files.")
    parser.add_argument("-g", "--genbank", required=True, help="Path to genbank file")
    parser.add_argument("-a", "--annotation", required=True, help="Annotation tool used (options allowed 'pharokka' or 'other')")
    parser.add_argument("-m", "--mapping", required=True, help="Path to mapping file")
    parser.add_argument("-s", "--sequencing", required=True, choices=["short-paired", "short-single", "long"], help="Type of sequencing (options allowed 'short-paired' and 'short-single' for short-read sequencing and 'long' for long-read sequencing)")
    parser.add_argument("-f", "--features", required=False, help="List of feature subplots to include (comma-separated) (options allowed: coverage, starts, misassembly)")
    parser.add_argument("-p", "--precision", required=False, default=100, help="Size of the average window for the features' calculation (in bp)")
    parser.add_argument("-w", "--width", required=False, default=1800, help="Width of the plot (in pixels)")
    parser.add_argument("-sh", "--subplot_height", required=False, default=130, help="Height of each subplot (in pixels)")
    args = parser.parse_args()

    # Read annotation file
    print("Reading genbank and mapping file...", flush=True)
    genbank_file = args.genbank
    record = check_single_locus(genbank_file)

    locus_name = record.name
    locus_size = len(record.seq)
    print(f"Locus: {record.name} ({locus_size} bp)", flush=True)

    allowed_types = ["CDS", "tRNA/tmRNA", "rRNA", "ncRNA", "ncRNA-region", "CRISPR", "Gap", "Misc"]
    filtered_features = [f for f in record.features if f.type in allowed_types]
    record.features = filtered_features

    mapping_file = args.mapping
    bam_name = args.mapping.split("/")[-1].replace(".bam", "")

    global ANNOTATION_TOOL
    ANNOTATION_TOOL = args.annotation

    # Plotting gene map
    print("Plotting gene map...", flush=True)
    max_visible_width = int(args.width)
    subplot_size = int(args.subplot_height)

    graphic_record = CustomTranslator().translate_record(record)
    # figure_width and figure_height for the arrow size
    annotation_fig = graphic_record.plot_with_bokeh(figure_width=30, figure_height=40)
    annotation_fig.width = max_visible_width
    annotation_fig.height = subplot_size

    shared_xrange = Range1d(0, locus_size)
    annotation_fig.x_range = shared_xrange

    # Adding the feature subplots
    requested_features = args.features.split(",") if args.features else []
    sequencing_type = args.sequencing

    # if starts in feature_list replace it by starts_plus, starts_minus, ends_plus, ends_minus
    feature_list = []
    for feature in requested_features:
        if feature == "starts":
            feature_list.extend(["coverage_reduced", "reads_starts", "reads_ends", "tau"])

        elif feature == "misassembly":
            if sequencing_type == "long":
                feature_list.extend(["insert_sizes"])
            if sequencing_type == "short-paired":
                feature_list.extend(["insert_sizes", "bad_orientations"])
            feature_list.extend(["left_clippings", "right_clippings", "insertions", "deletions", "mismatches"])

        else:
            feature_list.append(feature)

    sequencing_type = "short" if sequencing_type.startswith("short") else "long"
    window_size = int(args.precision)
    
    subplots = []
    if feature_list:
        subplots = plotting_features.adding_subplots(mapping_file, feature_list, locus_name, locus_size, sequencing_type, max_visible_width, subplot_size, shared_xrange, window_size)

    layout = column(annotation_fig, *subplots)
    output_path = f"MGFeaturesViewer_{bam_name}_mapped_on_{locus_name}.html"
    output_file(output_path)
    save(layout)
    print(f"Saved interactive plot to {output_path}")

if __name__ == "__main__":
    main()