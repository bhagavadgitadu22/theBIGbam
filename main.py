import argparse
import os
import sys
from Bio import SeqIO
import pysam
from calculating_data import calculating_all_features_parallel
from plotting_data import prepare_main_plot

def generate_html_plots(genbank_path, data_dictionary, annotation_tool, window_size, output_prefix, max_visible_width, subplot_size):
    print("Generating HTML plots per locus...", flush=True)

    allowed_types = ["CDS", "tRNA/tmRNA", "rRNA", "ncRNA", "ncRNA-region", "CRISPR", "Gap", "Misc"]
    # Parse GenBank file
    for record in SeqIO.parse(genbank_path, "genbank"):
        locus_name = record.id
        locus_size = len(record.seq)

        print(f"Processing locus: {locus_name} ({locus_size} bp)", flush=True)

        # Filter allowed feature types
        filtered_features = [f for f in record.features if f.type in allowed_types]
        record.features = filtered_features

        for bam_name, contigs_list in data_dictionary.items():
            # Create output filename
            output_html = f"{output_prefix}_{locus_name}_in_{bam_name}.html"

            # Extract contig data for this locus
            if locus_name not in contigs_list:
                print(f"Warning: No mapping data found for locus '{locus_name}'. Skipping.", flush=True)
                continue
            contig_data = contigs_list[locus_name]

            # Generate the HTML plot
            prepare_main_plot(
                contig_data,
                record,
                annotation_tool,
                locus_size,
                window_size,
                max_visible_width,
                subplot_size,
                output_html,
            )
            print(f"Saved: {output_html}", flush=True)

def main():
    # Parse command line arguments
    print("Parsing arguments...", flush=True)
    parser = argparse.ArgumentParser(description="Parse input files.")
    parser.add_argument("-t", "--threads", required=True, help="Number of threads available")
    parser.add_argument("-g", "--genbank", required=True, help="Path to genbank file of all investigated contigs")
    parser.add_argument("-m", "--mapping", required=True, help="Path to bam file or directory containing mapping files (BAM format)")
    parser.add_argument("-a", "--annotation", required=True, help="Annotation tool used (options allowed 'pharokka' or 'other')")
    parser.add_argument("-s", "--sequencing", required=True, choices=["short-paired", "short-single", "long"], help="Type of sequencing (options allowed 'short-paired' and 'short-single' for short-read sequencing and 'long' for long-read sequencing)")
    parser.add_argument("-f", "--features", required=False, help="List of feature subplots to include (comma-separated) (options allowed: starts, misassembly)")
    parser.add_argument("-p", "--precision", required=False, default=100, help="Size of the average window for the features' calculation (in bp)")
    parser.add_argument("-o", "--output_prefix", required=False, default="MGFeaturesViewer", help="Prefix for output files, including complete path if you want to save them in a specific folder")
    parser.add_argument("-pw", "--plot_width", required=False, default=1800, help="Width of the plot (in pixels)")
    parser.add_argument("-sh", "--subplot_height", required=False, default=130, help="Height of each subplot (in pixels)")
    args = parser.parse_args()

    # Checking genbank and mapping files are compatible... (sanity checks)
    print("### Checking that genbank and mapping files are compatible...", flush=True)
    genbank_loci = set(rec.id for rec in SeqIO.parse(args.genbank, "genbank"))
    if not genbank_loci:
        sys.exit("ERROR: No loci found in the provided GenBank file.")

    # Get list of BAM files
    if os.path.isdir(args.mapping):
        bam_files = [os.path.join(args.mapping, f) for f in os.listdir(args.mapping) if f.endswith(".bam")]
    else:
        bam_files = [args.mapping]
    if not bam_files:
        sys.exit("ERROR: No BAM files found in the specified mapping path.")

    # Check that references in BAM headers match loci in GenBank
    for bam_file in bam_files:
        try:
            bam = pysam.AlignmentFile(bam_file, "rb")
        except Exception as e:
            sys.exit(f"ERROR: Could not open BAM file '{bam_file}': {e}")
    bam_refs = set(bam.references)
    bam.close()

    # Check if all BAM references exist in the GenBank loci
    missing_in_genbank = bam_refs - genbank_loci
    missing_in_bam = genbank_loci - bam_refs
    if missing_in_genbank:
        raise ValueError(
            f"ERROR: References in BAM file '{os.path.basename(bam_file)}' not found in GenBank:\n"
            f"{', '.join(sorted(missing_in_genbank))}"
        )
    if missing_in_bam:
        print(f"Warning: Some GenBank loci not present in BAM '{os.path.basename(bam_file)}': "
              f"{', '.join(sorted(missing_in_bam))}")

    print("Sanity checks passed. All BAM reference names match GenBank loci.", flush=True)

    # Getting list of features requested
    # if starts in feature_list replace it by starts_plus, starts_minus, ends_plus, ends_minus
    requested_features = args.features.split(",") if args.features else []
    sequencing_type = args.sequencing
    window_size = int(args.precision)
    n_cores = int(args.threads)

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

    # Calculating values for all requested features from mapping files
    print("### Calculating values for all requested features from mapping files...", flush=True)
    data_dictionary = calculating_all_features_parallel(bam_files, feature_list, sequencing_type, window_size, args.output_prefix, n_cores)

    # Preparing final html files
    print("### Preparing final html files...", flush=True)
    max_visible_width = int(args.plot_width)
    subplot_size = int(args.subplot_height)
    generate_html_plots(args.genbank, data_dictionary, args.annotation, window_size, args.output_prefix, max_visible_width, subplot_size)
    
if __name__ == "__main__":
    main()