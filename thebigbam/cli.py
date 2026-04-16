import argparse
import sys

# Import command modules so we can share arg definitions and run functions
from thebigbam.utils import (
    read_mapping, add_sample_metadata, add_contig_metadata, add_mag_metadata,
)
from thebigbam.database import add_variable, calculating_data, export_data, inspect_blob
from thebigbam.plotting import start_bokeh_server
from thebigbam.analysis import (
    mapping_patterns_per_CDS_per_contig_per_sample,
    annotations_per_CDS_per_contig,
)

ANALYSIS_SCRIPTS = {
    'cds-mapping-patterns': (
        'Export per-CDS mapping signals (coverage, mismatches, indels, clippings)',
        mapping_patterns_per_CDS_per_contig_per_sample,
    ),
    'cds-annotations': (
        'Export all annotation info per CDS (qualifiers, GC, coordinates)',
        annotations_per_CDS_per_contig,
    ),
}

SCRIPTS = {
    'mapping-per-sample': 'Map reads for a single sample (one CSV row)',
    'calculate': "Run feature calculations over alignment files",
    'serve': "Start interactive Bokeh server",

    'export': 'Export a metric as contig x sample matrix (TSV)',

    'add-sample-metadata': "Add sample metadata from CSV as new columns in Sample table",
    'add-contig-metadata': "Add contig metadata from CSV as new columns in Contig table",
    'add-mag-metadata': "Add MAG metadata from CSV as new columns in MAG table",
    'add-variable': "Add an external variable from CSV to DB",

    'remove-sample': 'Remove a sample and all its data from DB',
    'remove-contig': 'Remove a contig and all its data from DB',
    'remove-mag': 'Remove a MAG and all its data from DB',
    'remove-sample-metadata': 'Remove a user-added metadata column from Sample table',
    'remove-contig-metadata': 'Remove a user-added metadata column from Contig table',
    'remove-mag-metadata': 'Remove a user-added metadata column from MAG table',
    'remove-variable': "Remove a Custom variable from DB",

    'list-samples': 'List samples from DB',
    'list-contigs': 'List contigs from DB',
    'list-mags': 'List MAGs from DB',
    'list-sample-metadata': 'List user-added metadata columns on Sample table',
    'list-contig-metadata': 'List user-added metadata columns on Contig table',
    'list-mag-metadata': 'List user-added metadata columns on MAG table',
    'list-variables': 'List variables and metadata from DB',

    'inspect': 'List values taken by a feature in a given contig-sample pair (decode and display BLOB contents)',
}

def build_argparser():
    p = argparse.ArgumentParser(prog="thebigbam", description="theBIGbam command-line front-end")
    sub = p.add_subparsers(dest="cmd", required=True)

    # main commands
    sp = sub.add_parser('mapping-per-sample', help=SCRIPTS['mapping-per-sample'])
    read_mapping.add_mapping_per_sample_args(sp)

    sp = sub.add_parser('calculate', help=SCRIPTS['calculate'])
    calculating_data.add_calculate_args(sp)

    sp = sub.add_parser('serve', help=SCRIPTS['serve'])
    start_bokeh_server.add_serve_args(sp)

    # export command
    sp = sub.add_parser('export', help=SCRIPTS['export'])
    export_data.add_export_args(sp)

    # adding data to the db
    sp = sub.add_parser('add-sample-metadata', help=SCRIPTS['add-sample-metadata'])
    add_sample_metadata.add_add_sample_metadata_args(sp)

    sp = sub.add_parser('add-contig-metadata', help=SCRIPTS['add-contig-metadata'])
    add_contig_metadata.add_add_contig_metadata_args(sp)

    sp = sub.add_parser('add-mag-metadata', help=SCRIPTS['add-mag-metadata'])
    add_mag_metadata.add_add_mag_metadata_args(sp)

    sp = sub.add_parser('add-variable', help=SCRIPTS['add-variable'])
    add_variable.add_add_variable_args(sp)

    # removing data from the db
    sp = sub.add_parser('remove-sample', help=SCRIPTS['remove-sample'])
    sp.add_argument('-d', '--db', required=True)
    sp.add_argument('--name', required=True, help='Name of the sample to remove')

    sp = sub.add_parser('remove-contig', help=SCRIPTS['remove-contig'])
    sp.add_argument('-d', '--db', required=True)
    sp.add_argument('--name', required=True, help='Name of the contig to remove')

    sp = sub.add_parser('remove-mag', help=SCRIPTS['remove-mag'])
    sp.add_argument('-d', '--db', required=True)
    sp.add_argument('--name', required=True, help='Name of the MAG to remove')

    sp = sub.add_parser('remove-sample-metadata', help=SCRIPTS['remove-sample-metadata'])
    sp.add_argument('-d', '--db', required=True)
    sp.add_argument('--colname', required=True, help='Name of the column to remove')

    sp = sub.add_parser('remove-contig-metadata', help=SCRIPTS['remove-contig-metadata'])
    sp.add_argument('-d', '--db', required=True)
    sp.add_argument('--colname', required=True, help='Name of the column to remove')

    sp = sub.add_parser('remove-mag-metadata', help=SCRIPTS['remove-mag-metadata'])
    sp.add_argument('-d', '--db', required=True)
    sp.add_argument('--colname', required=True, help='Name of the column to remove')

    sp = sub.add_parser('remove-variable', help=SCRIPTS['remove-variable'])
    add_variable.add_remove_variable_args(sp)

    # database inspection
    sp = sub.add_parser('list-samples', help=SCRIPTS['list-samples'])
    sp.add_argument('-d', '--db', required=True)

    sp = sub.add_parser('list-contigs', help=SCRIPTS['list-contigs'])
    sp.add_argument('-d', '--db', required=True)

    sp = sub.add_parser('list-mags', help=SCRIPTS['list-mags'])
    sp.add_argument('-d', '--db', required=True)

    sp = sub.add_parser('list-sample-metadata', help=SCRIPTS['list-sample-metadata'])
    sp.add_argument('-d', '--db', required=True)

    sp = sub.add_parser('list-contig-metadata', help=SCRIPTS['list-contig-metadata'])
    sp.add_argument('-d', '--db', required=True)

    sp = sub.add_parser('list-mag-metadata', help=SCRIPTS['list-mag-metadata'])
    sp.add_argument('-d', '--db', required=True)

    sp = sub.add_parser('list-variables', help=SCRIPTS['list-variables'])
    sp.add_argument('-d', '--db', required=True)
    sp.add_argument('--detailed', action='store_true', help='Enable detailed output')

    # inspect command
    sp = sub.add_parser('inspect', help=SCRIPTS['inspect'])
    inspect_blob.add_inspect_args(sp)

    # analysis subcommand group
    analysis_parser = sub.add_parser('analysis', help='Run analysis scripts on a theBIGbam database')
    analysis_sub = analysis_parser.add_subparsers(dest="analysis_cmd", required=True)
    for name, (desc, module) in ANALYSIS_SCRIPTS.items():
        asp = analysis_sub.add_parser(name, help=desc,
                                      description=getattr(module, 'DESCRIPTION', desc),
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
        module.add_args(asp)

    return p

def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    parser = build_argparser()
    args, extras = parser.parse_known_args(argv)

    # Warn about unused arguments
    if extras:
        print(f"Warning: Unknown/unused arguments provided: {' '.join(extras)}", file=sys.stderr)

    # Dispatch to module run functions (shared-args approach)
    if args.cmd == 'mapping-per-sample':
        try:
            return read_mapping.run_mapping_per_sample(args)
        except Exception as e:
            print(f"Error running mapping-per-sample: {e}")
            return 2

    if args.cmd == 'calculate':
        return calculating_data.run_calculate_args(args)

    if args.cmd == 'serve':
        return start_bokeh_server.run_serve(args)

    if args.cmd == 'export':
        return export_data.run_export(args)
    
    # adding data to the db (call into package functions)
    if args.cmd == 'add-sample-metadata':
        return add_sample_metadata.run_add_sample_metadata(args)

    if args.cmd == 'add-contig-metadata':
        return add_contig_metadata.run_add_contig_metadata(args)

    if args.cmd == 'add-mag-metadata':
        return add_mag_metadata.run_add_mag_metadata(args)

    if args.cmd == 'add-variable':
        return add_variable.run_add_variable(args)

    # removing data from the db (call into package functions)
    if args.cmd == 'remove-sample-metadata':
            try:
                from thebigbam.database import database_getters
                database_getters.remove_sample_metadata(args.db, args.colname)
                return 0
            except Exception as e:
                print(f"Error removing sample metadata: {e}")
                return 2

    if args.cmd == 'remove-contig-metadata':
        try:
            from thebigbam.database import database_getters
            database_getters.remove_contig_metadata(args.db, args.colname)
            return 0
        except Exception as e:
            print(f"Error removing contig metadata: {e}")
            return 2

    if args.cmd == 'remove-mag-metadata':
        try:
            from thebigbam.database import database_getters
            database_getters.remove_mag_metadata(args.db, args.colname)
            return 0
        except Exception as e:
            print(f"Error removing MAG metadata: {e}")
            return 2

    if args.cmd == 'remove-sample':
        try:
            from thebigbam.database import database_getters
            database_getters.remove_sample(args.db, args.name)
            return 0
        except Exception as e:
            print(f"Error removing sample: {e}")
            return 2

    if args.cmd == 'remove-contig':
        try:
            from thebigbam.database import database_getters
            database_getters.remove_contig(args.db, args.name)
            return 0
        except Exception as e:
            print(f"Error removing contig: {e}")
            return 2

    if args.cmd == 'remove-mag':
        try:
            from thebigbam.database import database_getters
            database_getters.remove_mag(args.db, args.name)
            return 0
        except Exception as e:
            print(f"Error removing MAG: {e}")
            return 2
        
    if args.cmd == 'remove-variable':
        return add_variable.run_remove_variable(args)

    # DB inspection commands (call into package functions)
    if args.cmd == 'list-samples':
        try:
            from thebigbam.database import database_getters
            database_getters.list_samples(args.db)
            return 0
        except Exception as e:
            print(f"Error listing samples: {e}")
            return 2

    if args.cmd == 'list-contigs':
        try:
            from thebigbam.database import database_getters
            database_getters.list_contigs(args.db)
            return 0
        except Exception as e:
            print(f"Error listing contigs: {e}")
            return 2

    if args.cmd == 'list-mags':
        try:
            from thebigbam.database import database_getters
            database_getters.list_mags_cli(args.db)
            return 0
        except Exception as e:
            print(f"Error listing MAGs: {e}")
            return 2

    if args.cmd == 'list-sample-metadata':
        try:
            from thebigbam.database import database_getters
            database_getters.list_sample_metadata(args.db)
            return 0
        except Exception as e:
            print(f"Error listing sample metadata: {e}")
            return 2

    if args.cmd == 'list-contig-metadata':
        try:
            from thebigbam.database import database_getters
            database_getters.list_contig_metadata(args.db)
            return 0
        except Exception as e:
            print(f"Error listing contig metadata: {e}")
            return 2

    if args.cmd == 'list-mag-metadata':
        try:
            from thebigbam.database import database_getters
            database_getters.list_mag_metadata(args.db)
            return 0
        except Exception as e:
            print(f"Error listing MAG metadata: {e}")
            return 2

    if args.cmd == 'inspect':
        return inspect_blob.run_inspect(args)

    if args.cmd == 'analysis':
        _, module = ANALYSIS_SCRIPTS[args.analysis_cmd]
        return module.run(args)

    if args.cmd == 'list-variables':
            try:
                from thebigbam.database import database_getters
                database_getters.list_variables(args.db, args.detailed)
                return 0
            except Exception as e:
                print(f"Error listing variables: {e}")
                return 2

    # fallback
    print("Unknown command", args.cmd)
    return 2

if __name__ == '__main__':
    raise SystemExit(main())
