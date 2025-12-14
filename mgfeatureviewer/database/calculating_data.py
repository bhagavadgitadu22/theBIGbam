import argparse, sys, os
from multiprocessing import cpu_count

# Import Rust bindings (required)
try:
    import mgfeatureviewer_rs as _rust
    HAS_RUST = True
except ImportError:
    HAS_RUST = False
    _rust = None


def calculating_all_features_parallel(list_modules, bam_files, output_db, min_coverage, curve_ratio, bar_ratio, circular=False, n_sample_cores=None, genbank_path=None, annotation_tool=""):
    """Process all BAM files in parallel using Rust bindings."""
    if not HAS_RUST:
        sys.exit("ERROR: Rust bindings (mgfeatureviewer_rs) are required but not available. Please install them first.")

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
            annotation_tool=annotation_tool,
            min_coverage=float(min_coverage),
            curve_ratio=float(curve_ratio),
            bar_ratio=float(bar_ratio),
            circular=circular,
            create_indexes=True,
        )
    except Exception as e:
        print(f"ERROR: Rust processing failed: {e}", flush=True)
        import traceback
        traceback.print_exc()
        sys.exit(1)

    total_time = result.get("total_time", 0.0)
    samples_processed = result.get("samples_processed", 0)
    samples_failed = result.get("samples_failed", 0)
    print(f"Finished {samples_processed} samples in {total_time:.2f}s ({total_time/max(1, samples_processed):.2f}s per sample avg)", flush=True)
    if samples_failed > 0:
        print(f"Warning: {samples_failed} samples failed to process", flush=True)

def add_calculate_args(parser):
    parser.add_argument("-t", "--threads", required=True, help="Number of threads available")
    parser.add_argument("-g", "--genbank", help="Path to genbank file (optional; if not provided, no gene annotations will be stored)")
    parser.add_argument("-b", "--bam_files", required=True, help="Path to bam file or directory containing mapping files (BAM format)")
    parser.add_argument("-m", "--modules", required=True, help="List of modules to compute (comma-separated) (options allowed: coverage, phagetermini, assemblycheck)")
    parser.add_argument("-o", "--output", required=True, help="Output database file path (.db)")
    parser.add_argument("--annotation_tool", default="", help="Optional: to color the contigs specify the annotation tool used (options allowed: pharokka)")
    parser.add_argument("--min_coverage", type=int, default=50, help="Minimum alignment-length coverage proportion for contig inclusion (default 50%% change threshold)")
    parser.add_argument('--curve_ratio', type=float, default=10, help='Compression ratio for curve plots (default: 10%%)')
    parser.add_argument('--bar_ratio', type=float, default=10, help='Compression ratio for bar plots (default: 10%%)')
    parser.add_argument("--circular", action="store_true", help="Set if assembly was doubled during mapping (enables modulo logic)")

def run_calculate_args(args):
    annotation_tool = args.annotation_tool
    genbank_path = getattr(args, 'genbank', None)

    if os.path.isdir(args.bam_files):
        bam_files = [os.path.join(args.bam_files, f) for f in os.listdir(args.bam_files) if f.endswith(".bam")]
    else:
        bam_files = [args.bam_files]
    if not bam_files:
        sys.exit("ERROR: No BAM files found in the specified mapping path.")

    requested_modules = args.modules.split(",")
    min_coverage = args.min_coverage
    curve_ratio = args.curve_ratio
    bar_ratio = args.bar_ratio
    circular = args.circular
    n_cores = int(args.threads)

    output_db = args.output
    if os.path.exists(output_db):
        sys.exit(f"ERROR: Output file '{output_db}' already exists. Please provide a new path to avoid overwriting.")

    if genbank_path and not os.path.exists(genbank_path):
        sys.exit(f"ERROR: GenBank file not found: {genbank_path}")

    print("Calculating values for all requested features from mapping files...", flush=True)
    calculating_all_features_parallel(
        requested_modules, bam_files, output_db, min_coverage, curve_ratio, bar_ratio, circular, n_cores,
        genbank_path=genbank_path, annotation_tool=annotation_tool if genbank_path else ""
    )

    print(f"\nOutput written to: {output_db}", flush=True)

def main():
    print("Parsing arguments...", flush=True)
    parser = argparse.ArgumentParser(description="Parse input files.")
    add_calculate_args(parser)
    args = parser.parse_args()
    run_calculate_args(args)


if __name__ == "__main__":
    main()