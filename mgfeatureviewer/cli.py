import argparse
import os
import subprocess
import sys

# Import command modules so we can share arg definitions and run functions
from mgfeatureviewer import (
    calculating_data,
    add_variable,
    plotting_data_per_sample,
    plotting_data_all_samples,
    start_bokeh_server,
)

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
    'list-contigs': 'List contigs from DB'
}

def build_argparser():
    p = argparse.ArgumentParser(prog="mgfeatureviewer", description="MGFeatureViewer command-line front-end")
    sub = p.add_subparsers(dest="cmd", required=True)

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