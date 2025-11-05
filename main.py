import argparse
import csv
from Bio import SeqIO
from calculating_data import generating_data_from_bam
from plotting_data import prepare_main_plot

def check_single_locus(gbk_file):
    try:
        record = SeqIO.read(gbk_file, "genbank")  # expects exactly 1 record
        return record
    except ValueError:
        raise SystemExit(f"Error: {gbk_file} must contain exactly 1 locus.")

def save_data_dictionary(data_dictionary, bam_name, locus_name, output_file):
    """
    Save a nested data dictionary to a CSV file.
    """
    # Define CSV headers
    fieldnames = ["sample", "contig", "variable", "position", "value"]

    with open(output_file, mode="w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for variable_name, content in data_dictionary.items():
            xs = content.get("x", [])
            ys = content.get("y", [])
            # Ensure lengths match
            for pos, val in zip(xs, ys):
                writer.writerow({
                    "sample": bam_name,
                    "contig": locus_name,
                    "variable": variable_name,
                    "position": pos,
                    "value": val
                })

def main():
    # Parse command line arguments
    print("Parsing arguments...", flush=True)
    parser = argparse.ArgumentParser(description="Parse input files.")
    parser.add_argument("-g", "--genbank", required=True, help="Path to genbank file")
    parser.add_argument("-a", "--annotation", required=True, help="Annotation tool used (options allowed 'pharokka' or 'other')")
    parser.add_argument("-m", "--mapping", required=True, help="Path to mapping file")
    parser.add_argument("-s", "--sequencing", required=True, choices=["short-paired", "short-single", "long"], help="Type of sequencing (options allowed 'short-paired' and 'short-single' for short-read sequencing and 'long' for long-read sequencing)")
    parser.add_argument("-f", "--features", required=False, help="List of feature subplots to include (comma-separated) (options allowed: starts, misassembly)")
    parser.add_argument("-p", "--precision", required=False, default=100, help="Size of the average window for the features' calculation (in bp)")
    parser.add_argument("-pw", "--plot_width", required=False, default=1800, help="Width of the plot (in pixels)")
    parser.add_argument("-sh", "--subplot_height", required=False, default=130, help="Height of each subplot (in pixels)")
    args = parser.parse_args()

    # Read files
    print("Reading genbank and mapping file...", flush=True)
    genbank_file = args.genbank
    record = check_single_locus(genbank_file)

    locus_name = record.name
    locus_size = len(record.seq)
    print(f"Locus: {record.name} ({locus_size} bp)", flush=True)

    allowed_types = ["CDS", "tRNA/tmRNA", "rRNA", "ncRNA", "ncRNA-region", "CRISPR", "Gap", "Misc"]
    filtered_genbank_features = [f for f in record.features if f.type in allowed_types]
    record.features = filtered_genbank_features

    mapping_file = args.mapping
    bam_name = args.mapping.split("/")[-1].replace(".bam", "")

    requested_features = args.features.split(",") if args.features else []

    protein_annotation_tool = args.annotation
    sequencing_type = args.sequencing
    window_size = int(args.precision)

    max_visible_width = int(args.plot_width)
    subplot_size = int(args.subplot_height)

    # Getting list of features requested
    # if starts in feature_list replace it by starts_plus, starts_minus, ends_plus, ends_minus
    feature_list = ["coverage"]  # always include coverage subplot
    for feature in requested_features:
        if feature == "starts":
            feature_list.extend(["coverage_reduced", "reads_starts", "reads_ends", "tau"])

        elif feature == "misassembly":
            if sequencing_type == "long":
                feature_list.extend(["read_lengths"])
            if sequencing_type == "short-paired":
                feature_list.extend(["insert_sizes", "bad_orientations"])
            feature_list.extend(["left_clippings", "right_clippings", "insertions", "deletions", "mismatches"])

        else:
            feature_list.append(feature)

    # Generating data from bam
    print("Generating data from mapping file...", flush=True)
    data_dictionary = generating_data_from_bam(mapping_file, locus_name, locus_size, feature_list, sequencing_type, window_size)
    
    print("Saving data as csv files...", flush=True)
    output_csv = f"MGFeaturesViewer_{bam_name}_mapped_on_{locus_name}.csv"
    save_data_dictionary(data_dictionary, bam_name, locus_name, output_csv)

    # Prepare main plot
    print("Preparing main plot...", flush=True)
    output_html = f"MGFeaturesViewer_{bam_name}_mapped_on_{locus_name}.html"
    prepare_main_plot(data_dictionary, record, protein_annotation_tool, locus_size, window_size, max_visible_width, subplot_size, output_html)

if __name__ == "__main__":
    main()
