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


def calculating_all_features_parallel(list_modules, bam_files, output_db, min_aligned_fraction, min_coverage_depth, curve_ratio, bar_ratio, n_sample_cores=None, sequencing_type=None, genbank_path=None, assembly_path=None, extend_db=None, min_occurrences=2, enable_timing=False):
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
            create_indexes=True,
            assembly_path=assembly_path if assembly_path else "",
            extend_db=extend_db if extend_db else "",
            min_occurrences=int(min_occurrences),
            enable_timing=enable_timing,
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
    parser.add_argument("-m", "--modules", required=False, default=None, help="List of modules to compute (comma-separated). If not provided, all modules are computed. Options: coverage, misalignment, rna, longreads, pairedreads, termini")
    parser.add_argument("-a", "--assembly", help="Path to assembly FASTA file (only needed for autoblast when genbank lacks sequence data)")
    parser.add_argument('-s', '--sequencing_type', choices=['long', 'paired-short', 'single-short'], help='Sequencing type (long or short allowed)')
    parser.add_argument("--min_aligned_fraction", type=int, default=50, help="Minimum alignment-length coverage proportion for contig inclusion (default: 50%%)")
    parser.add_argument("--min_coverage_depth", type=int, default=0, help="Minimum mean coverage depth for contig inclusion (disabled by default, e.g. 5 to filter low-depth contigs)")
    parser.add_argument('--coverage_percentage', type=float, default=10, help='Compressing ratio for features depending on coverage: only values above this %% of the local coverage are kept (default: 10%%)')
    parser.add_argument('--min_occurrences', type=int, default=2, help='Minimum absolute event count for sparse features (default: 2). Position kept only if value > coverage × coverage_percentage AND value > min_occurrences.')
    parser.add_argument('--extend', action='store_true', help='Extend an existing database with new samples (and optionally new contigs)')
    parser.add_argument('--time', action='store_true', default=False, help=argparse.SUPPRESS)

# CLI names → internal module names (stored in DB/Rust)
MODULE_ALIASES = {
    "coverage": "Coverage",
    "misalignment": "Misalignment",
    "rna": "RNA",
    "longreads": "Long-reads",
    "pairedreads": "Paired-reads",
    "phagetermini": "Phage termini",
}
VALID_MODULES = list(MODULE_ALIASES.keys())

def run_calculate_args(args):
    genbank_path = getattr(args, 'genbank', None)
    assembly_path = getattr(args, 'assembly', None)
    is_extending = getattr(args, 'extend', False)

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

    output_db = args.output

    # --- Extend mode validation ---
    if is_extending:
        if not os.path.exists(output_db):
            sys.exit(f"ERROR: Database '{output_db}' not found. --extend requires an existing database.")

        if not bam_files:
            sys.exit("ERROR: --extend requires BAM files to add new samples.")

        import duckdb
        conn = duckdb.connect(output_db, read_only=True)

        # Get existing sample names
        existing_samples = set()
        if conn.execute("SELECT 1 FROM information_schema.tables WHERE table_name = 'Sample'").fetchone():
            existing_samples = {r[0] for r in conn.execute("SELECT Sample_name FROM Sample").fetchall()}

        # Detect modules from existing Variable table
        existing_modules = sorted({
            r[0] for r in conn.execute("SELECT DISTINCT Module FROM Variable").fetchall()
        } & set(MODULE_ALIASES.values()))

        conn.close()

        # Check for sample name collisions
        new_sample_names = [Path(bam).stem for bam in bam_files]
        for name in new_sample_names:
            if name in existing_samples:
                sys.exit(
                    f"ERROR: Sample '{name}' already exists in database.\n"
                    f"  If this is a different sample, rename your BAM file.\n"
                    f"  If it is the same sample with updated mappings, first run:\n"
                    f"    thebigbam remove-sample --db {output_db} --name {name}\n"
                    f"  Then re-run calculate --extend."
                )

        # Override modules with existing DB modules
        if args.modules is not None:
            print(f"WARNING: --modules is ignored in extend mode. Using modules from existing database: {', '.join(existing_modules)}", flush=True)
        requested_modules = existing_modules

        # Print warning
        print(
            "\nWARNING: Extending an existing database. Please ensure:\n"
            "  - Contig names are consistent between your former and new BAM mappings\n"
            "  - If you added custom variables or metadata for existing samples/contigs,\n"
            "    you will need to add them for the new data as well\n"
            "  - If new contigs are added, existing samples will NOT have mapping data\n"
            "    for those contigs. To add that data, remove those samples with\n"
            "    `thebigbam remove-sample` and re-add them with updated BAM mappings\n"
            "    that include the new contigs.\n",
            flush=True,
        )
    else:
        # Normal mode
        if os.path.exists(output_db):
            sys.exit(f"ERROR: Output file '{output_db}' already exists. Please provide a new path to avoid overwriting.")

        # Handle modules - if no BAM files, modules are ignored
        if args.modules is None:
            requested_modules = list(MODULE_ALIASES.values()) if bam_files else []
        else:
            requested_modules = [m.strip().lower() for m in args.modules.split(",")]
            for module in requested_modules:
                if module not in MODULE_ALIASES:
                    sys.exit(f"ERROR: Unknown module '{module}'. Valid modules are: {', '.join(VALID_MODULES)}")
            requested_modules = [MODULE_ALIASES[m] for m in requested_modules]
            if not bam_files and requested_modules:
                print(f"WARNING: Modules {requested_modules} require BAM files - they will be skipped.", flush=True)
                requested_modules = []

    min_aligned_fraction = args.min_aligned_fraction
    min_coverage_depth = args.min_coverage_depth
    coverage_percentage = args.coverage_percentage
    n_cores = int(args.threads)

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
        requested_modules, bam_files, output_db, min_aligned_fraction, min_coverage_depth, coverage_percentage, coverage_percentage,
        n_sample_cores=n_cores,
        sequencing_type=args.sequencing_type, genbank_path=genbank_path,
        assembly_path=assembly_path,
        extend_db=output_db if is_extending else None,
        min_occurrences=args.min_occurrences,
        enable_timing=getattr(args, 'time', False),
    )


def main():
    print("Parsing arguments...", flush=True)
    parser = argparse.ArgumentParser(description="Parse input files.")
    add_calculate_args(parser)
    args = parser.parse_args()
    run_calculate_args(args)


if __name__ == "__main__":
    main()
