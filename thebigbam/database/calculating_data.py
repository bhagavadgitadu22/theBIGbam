import argparse, sys, os
import tempfile
from pathlib import Path
from multiprocessing import cpu_count

# Import Rust bindings (required)
try:
    import thebigbam_rs as _rust
    HAS_RUST = True
except ImportError:
    HAS_RUST = False
    _rust = None


def calculating_all_features_parallel(list_modules, bam_files, output_db, min_coverage, curve_ratio, bar_ratio, circular=False, n_sample_cores=None, sequencing_type=None, genbank_path=None, annotation_tool="", max_samples_in_memory=10, autoblast_file=None):
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
            annotation_tool=annotation_tool,
            min_coverage=float(min_coverage),
            curve_ratio=float(curve_ratio),
            bar_ratio=float(bar_ratio),
            circular=circular,
            create_indexes=True,
            max_samples_in_memory=max_samples_in_memory,
            autoblast_file=autoblast_file if autoblast_file else "",
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
    parser.add_argument("-t", "--threads", required=True, help="Number of threads available")
    parser.add_argument("-g", "--genbank", help="Path to annotation file: GenBank (.gbk, .gbff) or GFF3 (.gff) format (optional; if not provided, no gene annotations will be stored)")
    parser.add_argument("-b", "--bam_files", required=True, help="Path to bam file or directory containing mapping files (BAM format)")
    parser.add_argument("--circular", action="store_true", help="Set if assembly was doubled during mapping (enables modulo logic)")
    parser.add_argument("-o", "--output", required=True, help="Output database file path (.db)")
    parser.add_argument("-m", "--modules", required=False, default=None, help="List of modules to compute (comma-separated). If not provided, all modules are computed. Options: Coverage, Mapping metrics per position, Long-read metrics, Paired-read metrics, Phage termini")
    parser.add_argument("-a", "--assembly", help="Path to assembly FASTA file (only needed for autoblast when genbank lacks sequence data)")
    parser.add_argument("--annotation_tool", default="", help="Optional: to color the contigs specify the annotation tool used (options allowed: pharokka)")
    parser.add_argument('-s', '--sequencing_type', choices=['long', 'paired-short', 'single-short'], help='Sequencing type (long or short allowed)')
    parser.add_argument("--min_coverage", type=int, default=50, help="Minimum alignment-length coverage proportion for contig inclusion (default 50%% change threshold)")
    parser.add_argument('--variation_percentage', type=float, default=50, help='Run-length encoding ratio for independent features like coverage (default: 50%%)')
    parser.add_argument('--coverage_percentage', type=float, default=10, help='Compressing ratio for features depending on coverage: only values above this %% of the local coverage are kept (default: 10%%)')
    parser.add_argument("--max-samples-in-memory", type=int, default=10, help="Max samples to hold in memory during processing (default: 10). Lower for large datasets to reduce memory usage.")

VALID_MODULES = ["Coverage", "Mapping metrics per position", "Long-read metrics", "Paired-read metrics", "Phage termini"]

def run_calculate_args(args):
    annotation_tool = args.annotation_tool
    genbank_path = getattr(args, 'genbank', None)
    assembly_path = getattr(args, 'assembly', None)

    if os.path.isdir(args.bam_files):
        bam_files = [os.path.join(args.bam_files, f) for f in os.listdir(args.bam_files) if f.endswith(".bam")]
    else:
        bam_files = [args.bam_files]
    if not bam_files:
        sys.exit("ERROR: No BAM files found in the specified mapping path.")

    if args.modules is None:
        # Default: all modules
        requested_modules = VALID_MODULES.copy()
    else:
        requested_modules = [m.strip() for m in args.modules.split(",")]
        # Validate module names
        for module in requested_modules:
            if module not in VALID_MODULES:
                sys.exit(f"ERROR: Unknown module '{module}'. Valid modules are: {', '.join(VALID_MODULES)}")
    min_coverage = args.min_coverage
    variation_percentage = args.variation_percentage
    coverage_percentage = args.coverage_percentage
    circular = args.circular
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

    max_samples_in_memory = getattr(args, 'max_samples_in_memory', 10)

    # Warn if "Phage termini" is requested without genbank
    if "Phage termini" in requested_modules and not genbank_path:
        print("WARNING: No annotation file was provided: phage packaging will not be determined properly for contigs harboring a terminal repeat at both ends. Rerun with an annotation file or at least an assembly file to perform the \"Phage termini\" module properly.", flush=True)

    # Run autoblast if Phage termini module is requested
    autoblast_file = None
    if "Phage termini" in requested_modules and genbank_path:
        from thebigbam.utils.autoblast import perform_autoblast, extract_fasta_from_genbank
        print("Running autoblast for Phage termini analysis...", flush=True)

        # Determine FASTA source: use assembly file if provided, otherwise extract from genbank
        fasta_path = None
        temp_fasta = None

        if assembly_path:
            fasta_path = Path(assembly_path)
            print(f"  Using assembly file: {fasta_path}", flush=True)
        elif genbank_path:
            # Try to extract FASTA from genbank
            temp_fasta = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
            temp_fasta.close()
            temp_fasta_path = Path(temp_fasta.name)

            if extract_fasta_from_genbank(Path(genbank_path), temp_fasta_path):
                fasta_path = temp_fasta_path
                print(f"  Extracted sequences from GenBank file", flush=True)
            else:
                os.unlink(temp_fasta_path)
                print("  WARNING: Could not extract sequences from GenBank file (no sequence data found)", flush=True)
                print("  Provide --assembly parameter for autoblast functionality", flush=True)

        if fasta_path:
            # Create autoblast output file next to the output database
            output_dir = Path(output_db).parent
            output_dir.mkdir(parents=True, exist_ok=True)
            autoblast_output = output_dir / f"{Path(output_db).stem}_autoblast.tsv"

            try:
                perform_autoblast(n_cores, fasta_path, autoblast_output)
                autoblast_file = str(autoblast_output)
                print(f"  Autoblast results written to: {autoblast_file}", flush=True)
            except Exception as e:
                print(f"  WARNING: Autoblast failed: {e}", flush=True)
                print("  Continuing without autoblast results...", flush=True)

            # Clean up temp file if we created one
            if temp_fasta and os.path.exists(temp_fasta.name):
                os.unlink(temp_fasta.name)

    print("\nCalculating values for all requested features from mapping files...", flush=True)
    calculating_all_features_parallel(
        requested_modules, bam_files, output_db, min_coverage, variation_percentage, coverage_percentage, circular, n_cores,
        sequencing_type=args.sequencing_type, genbank_path=genbank_path, annotation_tool=annotation_tool if genbank_path else "",
        max_samples_in_memory=max_samples_in_memory, autoblast_file=autoblast_file
    )

def main():
    print("Parsing arguments...", flush=True)
    parser = argparse.ArgumentParser(description="Parse input files.")
    add_calculate_args(parser)
    args = parser.parse_args()
    run_calculate_args(args)


if __name__ == "__main__":
    main()