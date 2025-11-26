import argparse
import os
import sys

# Import command modules so we can share arg definitions and run functions
from mgfeatureviewer import (
    assembly_annotation,
    calculating_data,
    add_variable,
    plotting_data_per_sample,
    plotting_data_all_samples,
    read_mapping,
    start_bokeh_server,
)

# Path helpers
BASE_DIR = os.path.dirname(__file__)

SCRIPTS = {
    'mapping-per-sample': 'Map reads for a single sample (one CSV row)',
    'mapping-all-samples': 'Map reads for multiple samples listed in a CSV',
    'annotate-assemblies': 'Combine assemblies and run an annotator (pharokka/bakta)',

    'calculate': "Run feature calculations over BAMs",
    'add-variable': "Add an external variable from CSV to DB",

    'plot-per-sample': "Produce per-sample static HTML plot",
    'plot-all-samples': "Produce all-samples static HTML plot",
    'serve': "Start interactive Bokeh server",

    'list-variables': 'List variables and metadata from DB',
    'list-samples': 'List samples from DB',
    'list-contigs': 'List contigs from DB'
}

def build_argparser():
    p = argparse.ArgumentParser(prog="mgfeatureviewer", description="MGFeatureViewer command-line front-end")
    sub = p.add_subparsers(dest="cmd", required=True)

    # mapping commands (use shared add_*_args functions from mapping modules)
    sp = sub.add_parser('mapping-per-sample', help=SCRIPTS['mapping-per-sample'])
    read_mapping.add_mapping_per_sample_args(sp)

    sp = sub.add_parser('mapping-all-samples', help=SCRIPTS['mapping-all-samples'])
    read_mapping.add_mapping_all_args(sp)

    sp = sub.add_parser('annotate-assemblies', help=SCRIPTS['annotate-assemblies'])
    assembly_annotation.add_annotation_args(sp)

    # Pipeline macro: mapping -> annotation -> calculate
    sp = sub.add_parser('run-pipeline', help='Run mapping-all-samples, annotate-assemblies, then calculate sequentially')
    # Inputs for mapping step
    sp.add_argument('--csv', required=True, help='CSV for mapping (assembly,read1,read2,sequencing_type)')
    sp.add_argument('--circular', action='store_true', help='Treat assemblies as circular')
    # Inputs for annotation step
    sp.add_argument('--annotation-tool', required=True, choices=['pharokka', 'bakta'], help='Annotation tool to run')
    sp.add_argument('--annotation-db', required=True, help='Annotation tool DB identifier/path')
    sp.add_argument('--meta', action='store_true', help='Pass --meta to annotator if the combined assembly contain contigs from several organisms')
    # Inputs for calculate step
    sp.add_argument('--modules', required=True, help='Comma-separated modules for calculation (coverage,phagetermini,assemblycheck)')
    sp.add_argument('--db', required=True, help='Path to sqlite DB file to create for calculation (must NOT exist)')
    # Parallelization options
    sp.add_argument('--threads', type=int, default=4, help='Available CPUs (or CPUs per Slurm job if using Slurm)')
    sp.add_argument('--use-slurm', action='store_true', help='Dispatch mapping and calculation tasks as Slurm array (per-sample)')
    sp.add_argument('--max-concurrent', type=int, default=20, help='Maximal number of concurrent Slurm tasks')
    sp.add_argument('--max-time', default='02:00:00', help='Maximum time per Slurm job')
    sp.add_argument('--mem-per-cpu', type=int, default=8, help='Memory per CPU for Slurm jobs (in GB)')
    sp.add_argument('--parallelize-contigs', action='store_true', help='When running calculation per-sample, parallelize contigs inside the job')

    # calculate
    sp = sub.add_parser('calculate', help=SCRIPTS['calculate'])
    calculating_data.add_calculate_args(sp)

    # add-variable
    sp = sub.add_parser('add-variable', help=SCRIPTS['add-variable'])
    add_variable.add_add_variable_args(sp)

    # plot-per-sample
    sp = sub.add_parser('plot-per-sample', help=SCRIPTS['plot-per-sample'])
    plotting_data_per_sample.add_plot_per_sample_args(sp)

    # plot-all-samples
    sp = sub.add_parser('plot-all-samples', help=SCRIPTS['plot-all-samples'])
    plotting_data_all_samples.add_plot_all_args(sp)

    # serve
    sp = sub.add_parser('serve', help=SCRIPTS['serve'])
    start_bokeh_server.add_serve_args(sp)

    # Database inspection (kept simple)
    sp = sub.add_parser('list-variables', help=SCRIPTS['list-variables'])
    sp.add_argument('-d', '--db', required=True)
    sp.add_argument('--detailed', action='store_true', help='Enable detailed output')

    sp = sub.add_parser('list-samples', help=SCRIPTS['list-samples'])
    sp.add_argument('-d', '--db', required=True)

    sp = sub.add_parser('list-contigs', help=SCRIPTS['list-contigs'])
    sp.add_argument('-d', '--db', required=True)

    return p

def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    parser = build_argparser()
    args, extras = parser.parse_known_args(argv)

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
        # Use the simple, user-provided run-pipeline args.
        try:
            db_path = args.db
            db_dir = os.path.abspath(os.path.dirname(db_path)) or os.getcwd()

            # Mapping outputs go into mgfeatureviewer_bams inside the DB directory
            map_outdir = os.path.join(db_dir, 'mgfeatureviewer_bams')

            # Build mapping namespace expected by read_mapping.run_mapping_all
            map_ns = argparse.Namespace(
                csv=args.csv,
                assembly=None,
                circular=bool(getattr(args, 'circular', False)),
                output_dir=map_outdir,
                threads=int(getattr(args, 'threads_per_job', 4)),
                threads=int(getattr(args, 'threads_per_job', 4)),
                use_slurm=bool(getattr(args, 'use_slurm', False)),
                max_concurrent=int(getattr(args, 'max_concurrent', 20)),
                max_time=getattr(args, 'max_time', '02:00:00'),
                mem_per_cpu=f"{getattr(args, 'mem_per_cpu', 8)}G",
            )
            print(f"Starting mapping step (CSV={map_ns.csv}) -> outputs: {map_ns.output_dir}")
            print("Command args for mapping step:", dict(vars(map_ns)))
            #read_mapping.run_mapping_all(map_ns)

            # Annotation output target: place file in DB directory. Use .gbk target name;
            # assembly_annotation will copy whatever bakta/pharokka produced into this path.
            anno_target = os.path.join(db_dir, 'annotation_output.gbk')
            anno_ns = argparse.Namespace(
                csv=args.csv,
                assembly=None,
                annotation_tool=args.annotation_tool,
                annotation_db=args.annotation_db,
                meta=bool(getattr(args, 'meta', False)),
                threads=int(getattr(args, 'threads_per_job', 4)),
                genbank=anno_target,
            )
            print(f"Starting annotation step (annotation_tool={anno_ns.annotation_tool}) -> output: {anno_ns.genbank}")
            print("Command args for annotation step:", dict(vars(anno_ns)))
            #assembly_annotation.run_annotation(anno_ns)

            # Ensure annotation output exists (assembly_annotation should have created it)
            if not os.path.exists(anno_target):
                raise FileNotFoundError(f"Annotation output not found at expected location: {anno_target}")

            # Calculation namespace: point to the annotation file and the mapping BAM directory
            calc_ns = argparse.Namespace(
                threads=int(getattr(args, 'threads_per_job', 4)),
                genbank=anno_target,
                bam_files=map_outdir,
                modules=args.modules,
                db=db_path,
                annotation_tool=args.annotation_tool,
                min_coverage=50,
                step=50,
                outlier_threshold=3,
                derivative_threshold=3,
                max_points=10000,
                threads=int(getattr(args, 'threads_per_job', 4)),
                use_slurm=bool(getattr(args, 'use_slurm', False)),
                max_concurrent=int(getattr(args, 'max_concurrent', 20)),
                max_time=getattr(args, 'max_time', '02:00:00'),
                mem_per_cpu=f"{getattr(args, 'mem_per_cpu', 8)}G",
                parallelize_contigs=bool(getattr(args, 'parallelize_contigs', False)),
            )
            print(f"Starting calculation step (DB={calc_ns.db}) using genbank {calc_ns.genbank} and bams in {calc_ns.bam_files}")
            print("Command args for calculation step:", dict(vars(calc_ns)))
            calculating_data.run_calculate_args(calc_ns)
            return 0
        except Exception as e:
            print(f"Pipeline error: {e}")
            return 2

    # DB inspection commands (call into package functions)
    if args.cmd == 'list-variables':
        try:
            from mgfeatureviewer import database_getters
            database_getters.list_variables(args.db, args.detailed)
            return 0
        except Exception as e:
            print(f"Error listing variables: {e}")
            return 2

    if args.cmd == 'list-samples':
        try:
            from mgfeatureviewer import database_getters
            database_getters.list_samples(args.db)
            return 0
        except Exception as e:
            print(f"Error listing samples: {e}")
            return 2

    if args.cmd == 'list-contigs':
        try:
            from mgfeatureviewer import database_getters
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