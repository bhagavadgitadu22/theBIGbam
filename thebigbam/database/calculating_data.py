import argparse, sys, os
from pathlib import Path
from multiprocessing import cpu_count

# Import Rust bindings (required)
try:
    import thebigbam_rs as _rust
    HAS_RUST = True
except ImportError:
    HAS_RUST = False
    _rust = None


def calculating_all_features_parallel(list_modules, bam_files, output_db, min_aligned_fraction, min_coverage_depth, curve_ratio, bar_ratio, contig_variation_percentage=0.1, n_sample_cores=None, sequencing_type=None, genbank_path=None, assembly_path=None):
    """Process all BAM files in parallel using Rust bindings."""
    if not HAS_RUST:
        sys.exit("ERROR: Rust bindings (thebigbam_rs) are required but not available. Please install them first.")

    if n_sample_cores is None:
        n_sample_cores = max(1, cpu_count() - 1)

    print(f"Using Rust bindings to process {len(bam_files)} samples with rayon ({n_sample_cores} threads)...", flush=True)

    try:
        result = _rust.process_all_samples(
            genbank_path=genbank_path if genbank_path else "",
            bam_files=bam_files,
            output_db=output_db,
            modules=list_modules,
            threads=n_sample_cores,
            sequencing_type=sequencing_type,
            min_aligned_fraction=float(min_aligned_fraction),
            min_coverage_depth=float(min_coverage_depth),
            curve_ratio=float(curve_ratio),
            bar_ratio=float(bar_ratio),
            contig_variation_percentage=float(contig_variation_percentage),
            create_indexes=True,
            assembly_path=assembly_path if assembly_path else "",
        )
    except Exception as e:
        print(f"ERROR: Rust processing failed: {e}", flush=True)
        import traceback
        traceback.print_exc()
        sys.exit(1)

    samples_failed = result.get("samples_failed", 0)
    if samples_failed > 0:
        print(f"Warning: {samples_failed} samples failed to process", flush=True)

def add_calculate_args(parser):
    parser.add_argument('-t', '--threads', type=int, default=4, help='Number of threads available (default: 4)')
    parser.add_argument("-g", "--genbank", help="Path to annotation file: GenBank (.gbk, .gbff) or GFF3 (.gff) format. Required if no BAM files provided.")
    parser.add_argument("-b", "--bam_files", help="Path to bam file or directory containing mapping files (BAM format). Optional if genbank is provided.")
    parser.add_argument("-o", "--output", required=True, help="Output database file path (.db)")
    parser.add_argument("-m", "--modules", required=False, default=None, help="List of modules to compute (comma-separated). If not provided, all modules are computed. Options: Coverage, Misalignment, Long-reads, Paired-reads, Phage termini")
    parser.add_argument("-a", "--assembly", help="Path to assembly FASTA file (only needed for autoblast when genbank lacks sequence data)")
    parser.add_argument('-s', '--sequencing_type', choices=['long', 'paired-short', 'single-short'], help='Sequencing type (long or short allowed)')
    parser.add_argument("--min_aligned_fraction", type=int, default=50, help="Minimum alignment-length coverage proportion for contig inclusion (default: 50%%)")
    parser.add_argument("--min_coverage_depth", type=int, default=0, help="Minimum mean coverage depth for contig inclusion (disabled by default, e.g. 5 to filter low-depth contigs)")
    parser.add_argument('--variation_percentage', type=float, default=50, help='Run-length encoding ratio for independent features like coverage (default: 50%%)')
    parser.add_argument('--coverage_percentage', type=float, default=10, help='Compressing ratio for features depending on coverage: only values above this %% of the local coverage are kept (default: 10%%)')
    parser.add_argument('--contig_variation_percentage', type=float, default=10, help='Run-length encoding ratio for contig-level features like GC content (default: 10%%)')

VALID_MODULES = ["Coverage", "Misalignment", "Long-reads", "Paired-reads", "Phage termini"]

def run_calculate_args(args):
    genbank_path = getattr(args, 'genbank', None)
    assembly_path = getattr(args, 'assembly', None)

    # Handle optional bam_files
    bam_files = []
    if args.bam_files:
        if os.path.isdir(args.bam_files):
            bam_files = [os.path.join(args.bam_files, f) for f in os.listdir(args.bam_files) if f.endswith(".bam")]
            if not bam_files:
                print(f"WARNING: No .bam files found in directory '{args.bam_files}'", flush=True)
        elif os.path.isfile(args.bam_files):
            bam_files = [args.bam_files]
        else:
            sys.exit(f"ERROR: BAM path not found: {args.bam_files}")
    
    # Validate: need at least genbank OR bam_files
    if not bam_files and not genbank_path:
        sys.exit("ERROR: You must provide either --bam_files or --genbank (or both).")
    
    if not bam_files:
        print("No BAM files provided. Will only populate contig-level data from GenBank.", flush=True)

    # Handle modules - if no BAM files, modules are ignored
    if args.modules is None:
        # Default: all modules (but they require BAM files)
        requested_modules = VALID_MODULES.copy() if bam_files else []
    else:
        requested_modules = [m.strip() for m in args.modules.split(",")]
        # Validate module names
        for module in requested_modules:
            if module not in VALID_MODULES:
                sys.exit(f"ERROR: Unknown module '{module}'. Valid modules are: {', '.join(VALID_MODULES)}")
        # Warn if modules specified but no BAM files
        if not bam_files and requested_modules:
            print(f"WARNING: Modules {requested_modules} require BAM files - they will be skipped.", flush=True)
            requested_modules = []
    min_aligned_fraction = args.min_aligned_fraction
    min_coverage_depth = args.min_coverage_depth
    variation_percentage = args.variation_percentage
    coverage_percentage = args.coverage_percentage
    contig_variation_percentage = args.contig_variation_percentage
    n_cores = int(args.threads)

    output_db = args.output
    if os.path.exists(output_db):
        sys.exit(f"ERROR: Output file '{output_db}' already exists. Please provide a new path to avoid overwriting.")

    if genbank_path and not os.path.exists(genbank_path):
        sys.exit(f"ERROR: Annotation file not found: {genbank_path}")

    # Validate annotation file extension
    if genbank_path:
        valid_extensions = ('.gbk', '.gbff', '.gb', '.genbank', '.gff', '.gff3')
        if not genbank_path.lower().endswith(valid_extensions):
            sys.exit(f"ERROR: Unsupported annotation file format. Supported extensions: {', '.join(valid_extensions)}")

    if assembly_path and not os.path.exists(assembly_path):
        sys.exit(f"ERROR: Assembly file not found: {assembly_path}")

    print("\nCalculating values for all requested features from mapping files...", flush=True)
    calculating_all_features_parallel(
        requested_modules, bam_files, output_db, min_aligned_fraction, min_coverage_depth, variation_percentage, coverage_percentage,
        contig_variation_percentage=contig_variation_percentage, n_sample_cores=n_cores,
        sequencing_type=args.sequencing_type, genbank_path=genbank_path,
        assembly_path=assembly_path,
    )


def main():
    print("Parsing arguments...", flush=True)
    parser = argparse.ArgumentParser(description="Parse input files.")
    add_calculate_args(parser)
    args = parser.parse_args()
    run_calculate_args(args)


if __name__ == "__main__":
    main()