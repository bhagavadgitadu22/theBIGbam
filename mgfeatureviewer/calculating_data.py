import argparse, sys, os
from multiprocessing import cpu_count

# Import Rust bindings (required)
try:
    import mgfeatureviewer_rs as _rust
    HAS_RUST = True
except ImportError:
    HAS_RUST = False
    _rust = None


def calculating_all_features_parallel(list_modules, bam_files, output_db, min_coverage, step, z_thresh, deriv_thresh, max_points, n_sample_cores=None, genbank_path=None, annotation_tool=""):
    """Process all BAM files in parallel using Rust bindings."""
    if not HAS_RUST:
        sys.exit("ERROR: Rust bindings (mgfeatureviewer_rs) are required but not available. Please install them first.")

    if n_sample_cores is None:
        n_sample_cores = max(1, cpu_count() - 1)

    # Get parent directory of BAM files (assumes all in same directory)
    bam_dir = os.path.dirname(bam_files[0]) if bam_files else ""

    print(f"Using Rust bindings to process {len(bam_files)} samples with rayon ({n_sample_cores} threads)...", flush=True)

    # Call the Rust function that handles everything: GenBank parsing, parallel BAM processing,
    # and writing all data to SQLite
    result = _rust.process_all_samples(
        genbank_path=genbank_path,
        bam_dir=bam_dir,
        output_db=output_db,
        modules=list_modules,
        threads=n_sample_cores,
        annotation_tool=annotation_tool,
        min_coverage=float(min_coverage),
        step=step,
        z_thresh=float(z_thresh),
        deriv_thresh=float(deriv_thresh),
        max_points=max_points,
    )

    # Rust handles SQLite updates internally, just report results
    total_time = result.get("total_time", 0.0)
    samples_processed = result.get("samples_processed", 0)
    samples_failed = result.get("samples_failed", 0)
    print(f"Finished {samples_processed} samples in {total_time:.2f}s ({total_time/max(1, samples_processed):.2f}s per sample avg)", flush=True)
    if samples_failed > 0:
        print(f"Warning: {samples_failed} samples failed to process", flush=True)

### Main function helpers (shared-args)
def add_calculate_args(parser):
    parser.add_argument("-t", "--threads", required=True, help="Number of threads available")
    parser.add_argument("-g", "--genbank", required=True, help="Path to genbank file of all investigated contigs")
    parser.add_argument("-b", "--bam_files", required=True, help="Path to bam file or directory containing mapping files (BAM format)")
    parser.add_argument("-m", "--modules", required=True, help="List of modules to compute (comma-separated) (options allowed: coverage, phagetermini, assemblycheck)")
    parser.add_argument("-o", "--output", required=True, help="Output database file path (.db)")
    parser.add_argument("-a", "--annotation_tool", default="", help="Optional: to color the contigs specify the annotation tool used (options allowed: pharokka)")
    parser.add_argument("--min_coverage", type=int, default=50, help="Minimum alignment-length coverage proportion for contig inclusion")
    parser.add_argument("--step", type=int, default=50, help="Step size for compression (keep every Nth point in addition to the outliers)")
    parser.add_argument("--outlier_threshold", type=int, default=3, help="Points beyond mean+std*N are kept as outliers")
    parser.add_argument("--derivative_threshold", type=int, default=3, help="Points were the derivative is beyond mean+std*N are kept as outliers")
    parser.add_argument("--max_points", type=int, default=10000, help="Maximum number of points kept during compression")

def run_calculate_args(args):
    annotation_tool = args.annotation_tool

    # Get list of BAM files
    if os.path.isdir(args.bam_files):
        bam_files = [os.path.join(args.bam_files, f) for f in os.listdir(args.bam_files) if f.endswith(".bam")]
    else:
        bam_files = [args.bam_files]
    if not bam_files:
        sys.exit("ERROR: No BAM files found in the specified mapping path.")

    # Requested modules
    requested_modules = args.modules.split(",")

    # Parameters for compression
    min_coverage = args.min_coverage
    step = args.step
    z_thresh = args.outlier_threshold
    deriv_thresh = args.derivative_threshold
    max_points = args.max_points

    n_cores = int(args.threads)

    # Setup output database path
    output_db = args.output
    if os.path.exists(output_db):
        sys.exit(f"ERROR: Output file '{output_db}' already exists. Please provide a new path to avoid overwriting.")

    # Rust handles everything: GenBank parsing, database creation, parallel BAM processing, and output writing
    print("Calculating values for all requested features from mapping files...", flush=True)
    calculating_all_features_parallel(
        requested_modules, bam_files, output_db, min_coverage, step, z_thresh, deriv_thresh, max_points, n_cores,
        genbank_path=args.genbank, annotation_tool=annotation_tool
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