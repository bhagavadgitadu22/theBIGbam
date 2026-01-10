import argparse
import os
import sys

# Import command modules so we can share arg definitions and run functions
from mgfeatureviewer.utils import (
    read_mapping, assembly_annotation
)
from mgfeatureviewer.database import add_variable, calculating_data
from mgfeatureviewer.plotting import plotting_data_all_samples, plotting_data_per_sample, start_bokeh_server

# Path helpers
BASE_DIR = os.path.dirname(__file__)

SCRIPTS = {
    'calculate': "Run feature calculations over BAMs",
    'add-variable': "Add an external variable from CSV to DB",

    'plot-per-sample': "Produce per-sample static HTML plot",
    'plot-all-samples': "Produce all-samples static HTML plot",
    'serve': "Start interactive Bokeh server",

    'list-variables': 'List variables and metadata from DB',
    'list-samples': 'List samples from DB',
    'list-contigs': 'List contigs from DB',

    'mapping-per-sample': 'Map reads for a single sample (one CSV row)',
    'mapping-all-samples': 'Map reads for multiple samples listed in a CSV',
    'annotate-assemblies': 'Combine assemblies and run an annotator (pharokka/bakta)'
}

def build_argparser():
    p = argparse.ArgumentParser(prog="mgfeatureviewer", description="MGFeatureViewer command-line front-end")
    sub = p.add_subparsers(dest="cmd", required=True)

    # modify database commands
    sp = sub.add_parser('calculate', help=SCRIPTS['calculate'])
    calculating_data.add_calculate_args(sp)

    sp = sub.add_parser('add-variable', help=SCRIPTS['add-variable'])
    add_variable.add_add_variable_args(sp)

    # plotting commands
    sp = sub.add_parser('plot-per-sample', help=SCRIPTS['plot-per-sample'])
    plotting_data_per_sample.add_plot_per_sample_args(sp)

    sp = sub.add_parser('plot-all-samples', help=SCRIPTS['plot-all-samples'])
    plotting_data_all_samples.add_plot_all_args(sp)

    sp = sub.add_parser('serve', help=SCRIPTS['serve'])
    start_bokeh_server.add_serve_args(sp)

    # database inspection
    sp = sub.add_parser('list-variables', help=SCRIPTS['list-variables'])
    sp.add_argument('-d', '--db', required=True)
    sp.add_argument('--detailed', action='store_true', help='Enable detailed output')

    sp = sub.add_parser('list-samples', help=SCRIPTS['list-samples'])
    sp.add_argument('-d', '--db', required=True)

    sp = sub.add_parser('list-contigs', help=SCRIPTS['list-contigs'])
    sp.add_argument('-d', '--db', required=True)

    # mapping commands (use shared add_*_args functions from mapping modules)
    sp = sub.add_parser('mapping-per-sample', help=SCRIPTS['mapping-per-sample'])
    read_mapping.add_mapping_per_sample_args(sp)

    sp = sub.add_parser('mapping-all-samples', help=SCRIPTS['mapping-all-samples'])
    read_mapping.add_mapping_all_args(sp)

    sp = sub.add_parser('annotate-assemblies', help=SCRIPTS['annotate-assemblies'])
    assembly_annotation.add_annotation_args(sp)

    # pipeline macro: mapping -> annotation -> calculate
    # Reuses argument definitions from subcommands to avoid duplication
    sp = sub.add_parser('run-pipeline', help='Run mapping, annotation, and calculation in sequence')
    
    # Mapping inputs - add both per-sample and all-samples args, make csv/read1 mutually exclusive
    mapping_group = sp.add_mutually_exclusive_group(required=True)
    mapping_group.add_argument('--csv', help='CSV file for mapping multiple samples (uses mapping-all-samples logic)')
    mapping_group.add_argument('-r1', '--read1', help='Read1 file for single-sample mapping (uses mapping-per-sample logic)')
    
    # Common mapping args (shared between both mapping modes)
    sp.add_argument('-r2', '--read2', help='Read2 file (optional, for paired-end reads)')
    sp.add_argument('-a', '--assembly', required=True, help='Reference assembly (fasta file)')
    sp.add_argument('-s', '--sequencing_type', required=True, choices=['long', 'paired-short', 'single-short'], help='Sequencing type: use "paired-short", "single-short" or "long". Use optionally if all your samples were sequenced the same way. If not specified, the tool will infer the sequencing type per bam file.')
    sp.add_argument('--circular', action='store_true', help='Treat assemblies as circular')
    sp.add_argument('-t', '--threads', type=int, default=4, help='Number of threads (default: 4)')
    
    # Annotation inputs (optional - skip annotation step if not provided)
    sp.add_argument('--annotation_tool', choices=['pharokka', 'bakta'], help='Annotation tool (optional - skip annotation if not provided)')
    sp.add_argument('--annotation_db', help='Annotation database path (required if --annotation_tool is set)')
    sp.add_argument('--meta', action='store_true', help='Pass --meta to annotator for multi-organism assemblies')
    
    # Calculation inputs
    sp.add_argument('-m', '--modules', required=True, help='Comma-separated modules (coverage,phagetermini,assemblycheck)')
    sp.add_argument('--min_coverage', type=int, default=50, help='Minimum coverage for contig inclusion (default: 50%%)')
    sp.add_argument('--variation_percentage', type=float, default=50, help='Run-length encoding ratio for independent features like coverage (default: 50%%)')
    sp.add_argument('--coverage_percentage', type=float, default=10, help='Compressing ratio for features depending on coverage: only values above this %% of the local coverage are kept (default: 10%%)')

    # Output
    sp.add_argument('-o', '--output', required=True, help='Output directory (must NOT exist)')

    return p

def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    parser = build_argparser()
    args, extras = parser.parse_known_args(argv)

    # Warn about unused arguments
    if extras:
        print(f"Warning: Unknown/unused arguments provided: {' '.join(extras)}", file=sys.stderr)

    # Dispatch to module run functions (shared-args approach)
    if args.cmd == 'calculate':
        return calculating_data.run_calculate_args(args)

    if args.cmd == 'add-variable':
        return add_variable.run_add_variable(args)

    if args.cmd == 'plot-per-sample':
        return plotting_data_per_sample.run_plot_per_sample(args)

    if args.cmd == 'plot-all-samples':
        return plotting_data_all_samples.run_plot_all(args)

    if args.cmd == 'serve':
        return start_bokeh_server.run_serve(args)

    # mapping commands
    if args.cmd == 'mapping-per-sample':
        try:
            return read_mapping.run_mapping_per_sample(args)
        except Exception as e:
            print(f"Error running mapping-per-sample: {e}")
            return 2

    if args.cmd == 'mapping-all-samples':
        try:
            return read_mapping.run_mapping_all(args)
        except Exception as e:
            print(f"Error running mapping-all-samples: {e}")
            return 2

    if args.cmd == 'annotate-assemblies':
        try:
            return assembly_annotation.run_annotation(args)
        except Exception as e:
            print(f"Error running annotate-assemblies: {e}")
            return 2

    if args.cmd == 'run-pipeline':
        # Run mapping -> [annotation] -> calculate sequentially
        # Annotation is optional - skip if --annotation_tool not provided
        try:
            output_dir = os.path.abspath(args.output)
            os.makedirs(output_dir, exist_ok=True)

            # Mapping outputs go into a bams subdirectory
            map_outdir = os.path.join(output_dir, 'bams')

            # Determine step count based on whether annotation is enabled
            run_annotation = args.annotation_tool is not None
            total_steps = 3 if run_annotation else 2
            step = 1

            # Determine if single-sample or multi-sample mapping
            if args.csv:
                # Multi-sample mapping using CSV
                map_ns = argparse.Namespace(
                    csv=args.csv,
                    assembly=args.assembly,
                    sequencing_type=args.sequencing_type,
                    circular=args.circular,
                    output_dir=map_outdir,
                    threads=args.threads,
                )
                print(f"[{step}/{total_steps}] Mapping: multiple samples from CSV -> {map_outdir}")
                read_mapping.run_mapping_all(map_ns)
            else:
                # Single-sample mapping using read1/read2
                if not args.read1:
                    raise ValueError("Either --csv or --read1 must be provided")

                os.makedirs(map_outdir, exist_ok=True)
                output_bam = os.path.join(map_outdir, f"{os.path.basename(args.assembly).split('.')[0]}.bam")

                map_ns = argparse.Namespace(
                    read1=args.read1,
                    read2=args.read2,
                    assembly=args.assembly,
                    sequencing_type=args.sequencing_type,
                    circular=args.circular,
                    output=output_bam,
                    threads=args.threads,
                )
                print(f"[{step}/{total_steps}] Mapping: single sample -> {output_bam}")
                read_mapping.run_mapping_per_sample(map_ns)
            step += 1

            # Annotation step (optional)
            anno_target = ""
            if run_annotation:
                if not args.annotation_db:
                    raise ValueError("--annotation_db is required when --annotation_tool is specified")

                anno_target = os.path.join(output_dir, 'annotation.gbk')
                anno_ns = argparse.Namespace(
                    csv=args.csv if args.csv else None,
                    assembly=args.assembly if not args.csv else None,
                    annotation_tool=args.annotation_tool,
                    annotation_db=args.annotation_db,
                    meta=args.meta,
                    threads=args.threads,
                    genbank=anno_target,
                )
                print(f"[{step}/{total_steps}] Annotation: {args.annotation_tool} -> {anno_target}")
                assembly_annotation.run_annotation(anno_ns)

                if not os.path.exists(anno_target):
                    raise FileNotFoundError(f"Annotation failed: {anno_target} not created")
                step += 1

            # Calculate step - use actual database path as output
            final_db = os.path.join(output_dir, 'features.db')
            calc_ns = argparse.Namespace(
                threads=args.threads,
                genbank=anno_target,  # Empty string if no annotation
                bam_files=map_outdir,
                modules=args.modules,
                output=final_db,
                annotation_tool=args.annotation_tool if args.annotation_tool else "",
                sequencing_type=args.sequencing_type,
                min_coverage=args.min_coverage,
                variation_percentage=args.variation_percentage,
                coverage_percentage=args.coverage_percentage,
                circular=args.circular,
            )

            print(f"[{step}/{total_steps}] Calculate: modules={args.modules} -> {final_db}")
            calculating_data.run_calculate_args(calc_ns)

            print(f"\nPipeline complete! Output: {final_db}")
            return 0
        except Exception as e:
            print(f"Pipeline error: {e}")
            import traceback
            traceback.print_exc()
            return 2

    # DB inspection commands (call into package functions)
    if args.cmd == 'list-variables':
        try:
            from mgfeatureviewer.database import database_getters
            database_getters.list_variables(args.db, args.detailed)
            return 0
        except Exception as e:
            print(f"Error listing variables: {e}")
            return 2

    if args.cmd == 'list-samples':
        try:
            from mgfeatureviewer.database import database_getters
            database_getters.list_samples(args.db)
            return 0
        except Exception as e:
            print(f"Error listing samples: {e}")
            return 2

    if args.cmd == 'list-contigs':
        try:
            from mgfeatureviewer.database import database_getters
            database_getters.list_contigs(args.db)
            return 0
        except Exception as e:
            print(f"Error listing contigs: {e}")
            return 2

    # fallback
    print("Unknown command", args.cmd)
    return 2

if __name__ == '__main__':
    raise SystemExit(main())
