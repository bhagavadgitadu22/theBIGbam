import duckdb

# Columns to exclude from annotation filtering UI (internal/metadata columns)
ANNOTATION_EXCLUDED_COLUMNS = {'Contig_id', 'Start', 'End', 'Strand', 'Longest_isoform', 'Locus_tag'}


def update_database_metadata(conn):
    """Update Date_of_last_modification and Tool_version_used_for_last_modification."""
    from importlib.metadata import version
    import subprocess
    from datetime import datetime

    tool_version = version('thebigbam')
    try:
        h = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'],
                          capture_output=True, text=True).stdout.strip()
        if h:
            tool_version = f"{tool_version}+{h}"
    except Exception:
        pass

    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    try:
        conn.execute("UPDATE Database_metadata SET Value = ? WHERE Key = 'Date_of_last_modification'", [now])
        conn.execute("UPDATE Database_metadata SET Value = ? WHERE Key = 'Tool_version_used_for_last_modification'", [tool_version])
    except Exception:
        pass  # Older databases without Database_metadata table


def get_filtering_metadata(db_path: str) -> dict:
    """
    Get column metadata for Filtering2 UI.

    Returns cached metadata for Filtering2.
    Structure: {
        category: {
            'source': 'table_or_view_name',
            'columns': {
                col_name: {
                    'type': 'text' | 'numeric',
                    'distinct_values': [...] | None  # Only for text columns
                }
            }
        }
    }
    """
    conn = duckdb.connect(db_path, read_only=True)

    # Define category mappings
    category_config = {
        'Contig': {
            'source': 'Contig',
            'exclude': ['Contig_id']
        },
        'Sample': {
            'source': 'Sample',
            'exclude': ['Sample_id']
        },
        'Coverage': {
            'source': 'Explicit_coverage',
            'exclude': ['Contig_name', 'Sample_name']
        },
        'Misassembly': {
            'source': 'Explicit_misassembly',
            'exclude': ['Contig_name', 'Sample_name']
        },
        'Microdiversity': {
            'source': 'Explicit_microdiversity',
            'exclude': ['Contig_name', 'Sample_name']
        },
        'Side misassembly': {
            'source': 'Explicit_side_misassembly',
            'exclude': ['Contig_name', 'Sample_name']
        },
        'Topology': {
            'source': 'Explicit_topology',
            'exclude': ['Contig_name', 'Sample_name']
        },
        'Termini': {
            'source': 'Explicit_phage_mechanisms',
            'exclude': ['Contig_name', 'Sample_name']
        }
    }

    # Text type names in DuckDB
    text_types = {'VARCHAR', 'TEXT', 'STRING'}
    bool_types = {'BOOLEAN', 'BOOL'}

    result = {}
    for category, config in category_config.items():
        source = config['source']
        exclude = set(config['exclude'])

        # Check if table/view exists
        try:
            cols_info = conn.execute(f"DESCRIBE {source}").fetchall()
        except Exception:
            # Table/view doesn't exist, skip this category
            continue

        columns = {}
        for col_name, col_type, *_ in cols_info:
            if col_name in exclude:
                continue

            is_bool = col_type.upper() in bool_types
            is_text = is_bool or any(t in col_type.upper() for t in text_types)
            col_data = {'type': 'text' if is_text else 'numeric'}

            # For boolean columns, expose as text with yes/no values
            if is_bool:
                col_data['distinct_values'] = ['yes', 'no']
                col_data['is_bool'] = True

            # For text columns, get distinct values
            elif is_text:
                try:
                    distinct = conn.execute(
                        f"SELECT DISTINCT \"{col_name}\" FROM {source} WHERE \"{col_name}\" IS NOT NULL ORDER BY \"{col_name}\""
                    ).fetchall()
                    col_data['distinct_values'] = [row[0] for row in distinct]
                except Exception:
                    col_data['distinct_values'] = []
                
                # Skip columns with only NULL values (no distinct non-NULL values)
                if not col_data['distinct_values']:
                    continue
            else:
                # For numeric columns, check if there are any non-NULL values
                try:
                    has_values = conn.execute(
                        f"SELECT 1 FROM {source} WHERE \"{col_name}\" IS NOT NULL LIMIT 1"
                    ).fetchone()
                    if not has_values:
                        continue  # Skip columns with only NULL values
                except Exception:
                    continue

            columns[col_name] = col_data

        if columns:  # Only add category if it has columns
            result[category] = {
                'source': source,
                'columns': columns
            }

    # Add annotation columns to Contig category from Contig_annotation table
    if 'Contig' in result:
        try:
            ann_cols_info = conn.execute("DESCRIBE Contig_annotation").fetchall()
            for col_name, col_type, *_ in ann_cols_info:
                if col_name in ANNOTATION_EXCLUDED_COLUMNS:
                    continue

                is_text = any(t in col_type.upper() for t in text_types)
                col_data = {
                    'type': 'text' if is_text else 'numeric',
                    'source': 'Contig_annotation'  # Mark as annotation column
                }

                # For text columns, get distinct values
                if is_text:
                    try:
                        distinct = conn.execute(
                            f'SELECT DISTINCT "{col_name}" FROM Contig_annotation WHERE "{col_name}" IS NOT NULL ORDER BY "{col_name}"'
                        ).fetchall()
                        col_data['distinct_values'] = [row[0] for row in distinct]
                    except Exception:
                        col_data['distinct_values'] = []
                    
                    # Skip columns with only NULL values
                    if not col_data['distinct_values']:
                        continue
                else:
                    # For numeric columns, check if there are any non-NULL values
                    try:
                        has_values = conn.execute(
                            f'SELECT 1 FROM Contig_annotation WHERE "{col_name}" IS NOT NULL LIMIT 1'
                        ).fetchone()
                        if not has_values:
                            continue  # Skip columns with only NULL values
                    except Exception:
                        continue

                result['Contig']['columns'][col_name] = col_data
        except Exception:
            pass  # Contig_annotation table doesn't exist

    conn.close()
    return result


def list_variables(db_path, detailed=False):
    """Print variables and detailed metadata from Variable table (excluding Feature_table_name)."""
    conn = duckdb.connect(db_path, read_only=True)
    cur = conn.cursor()

    # DuckDB: get column names from DESCRIBE statement
    cur.execute("DESCRIBE Variable")
    cols = [r[0] for r in cur.fetchall()]

    # Fields we will display (exclude Feature_table_name)
    display_fields = [c for c in cols if not(c in ['Variable_id', 'Feature_table_name'])]

    # Query the table for the display fields
    sel = ", ".join(display_fields)
    cur.execute(f"SELECT {sel} FROM Variable ORDER BY Variable_name")
    rows = cur.fetchall()

    if not rows:
        print("No variables found in the database.")
        conn.close()
        return

    # Print header
    for row in rows:
        # Pair field name and value and print nicely
        print(row[0])  # Variable_name as header
        if detailed:
            # Print other fields minus variable name
            for fname, val in zip(display_fields[1:], row[1:]):
                print(f"- {fname}: {val}")
            print("")

    conn.close()

def list_samples(db_path):
    """Print Sample_name values from Sample table."""
    conn = duckdb.connect(db_path, read_only=True)
    cur = conn.cursor()
    cur.execute("SELECT 1 FROM information_schema.tables WHERE table_name = 'Sample'")
    if cur.fetchone() is None:
        print("No samples in the database (genbank-only mode).")
        conn.close()
        return
    cur.execute("SELECT Sample_name FROM Sample ORDER BY Sample_name")
    rows = [r[0] for r in cur.fetchall()]
    if not rows:
        print("No samples found in the database.")
    else:
        for s in rows:
            print(f"{s}")
    conn.close()

def list_contigs(db_path):
    """Print Contig_name values from Contig table."""
    conn = duckdb.connect(db_path, read_only=True)
    cur = conn.cursor()
    cur.execute("SELECT Contig_name FROM Contig ORDER BY Contig_name")
    rows = [r[0] for r in cur.fetchall()]
    if not rows:
        print("No contigs found in the database.")
    else:
        for c in rows:
            print(f"{c}")
    conn.close()

SAMPLE_INTERNAL_COLUMNS = {
    'Sample_id', 'Sample_name', 'Sequencing_type',
    'Number_of_reads', 'Number_of_mapped_reads',
}

CONTIG_INTERNAL_COLUMNS = {
    'Contig_id', 'Contig_name', 'Contig_length',
    'Duplication_percentage', 'GC_mean', 'GC_sd', 'GC_skew_amplitude', 'Positive_GC_skew_windows_percentage',
}


def list_sample_metadata(db_path):
    """Print user-added column names on the Sample table."""
    conn = duckdb.connect(db_path, read_only=True)
    if conn.execute("SELECT 1 FROM information_schema.tables WHERE table_name = 'Sample'").fetchone() is None:
        conn.close()
        print("No Sample table in database (genbank-only mode).")
        return
    cols = [r[0] for r in conn.execute("DESCRIBE Sample").fetchall()]
    conn.close()
    user_cols = [c for c in cols if c not in SAMPLE_INTERNAL_COLUMNS]
    if not user_cols:
        print("No user-added metadata columns on Sample table.")
    else:
        for c in user_cols:
            print(c)


def list_contig_metadata(db_path):
    """Print user-added column names on the Contig table."""
    conn = duckdb.connect(db_path, read_only=True)
    cols = [r[0] for r in conn.execute("DESCRIBE Contig").fetchall()]
    conn.close()
    user_cols = [c for c in cols if c not in CONTIG_INTERNAL_COLUMNS]
    if not user_cols:
        print("No user-added metadata columns on Contig table.")
    else:
        for c in user_cols:
            print(c)


def remove_sample_metadata(db_path, colname):
    """Remove a user-added column from the Sample table."""
    if colname in SAMPLE_INTERNAL_COLUMNS:
        print(f"Error: '{colname}' is a built-in column and cannot be removed.")
        return
    conn = duckdb.connect(db_path)
    if conn.execute("SELECT 1 FROM information_schema.tables WHERE table_name = 'Sample'").fetchone() is None:
        conn.close()
        print("No Sample table in database (genbank-only mode).")
        return
    cols = [r[0] for r in conn.execute("DESCRIBE Sample").fetchall()]
    if colname not in cols:
        conn.close()
        print(f"Error: column '{colname}' does not exist on Sample table.")
        return
    conn.execute(f'ALTER TABLE Sample DROP COLUMN "{colname}"')
    update_database_metadata(conn)
    conn.close()
    print(f"Removed column '{colname}' from Sample table.")


def remove_contig_metadata(db_path, colname):
    """Remove a user-added column from the Contig table."""
    if colname in CONTIG_INTERNAL_COLUMNS:
        print(f"Error: '{colname}' is a built-in column and cannot be removed.")
        return
    conn = duckdb.connect(db_path)
    cols = [r[0] for r in conn.execute("DESCRIBE Contig").fetchall()]
    if colname not in cols:
        conn.close()
        print(f"Error: column '{colname}' does not exist on Contig table.")
        return
    conn.execute(f'ALTER TABLE Contig DROP COLUMN "{colname}"')
    update_database_metadata(conn)
    conn.close()
    print(f"Removed column '{colname}' from Contig table.")


def _table_exists(conn, table_name):
    """Check if a table exists in the database."""
    return conn.execute(
        "SELECT 1 FROM information_schema.tables WHERE table_name = ?", [table_name]
    ).fetchone() is not None


def _delete_from(conn, table_name, column, value):
    """Delete rows from a table if it exists. Returns number of deleted rows."""
    if not _table_exists(conn, table_name):
        return 0
    count = conn.execute(
        f'SELECT COUNT(*) FROM "{table_name}" WHERE "{column}" = ?', [value]
    ).fetchone()[0]
    if count > 0:
        conn.execute(f'DELETE FROM "{table_name}" WHERE "{column}" = ?', [value])
    return count


def _get_feature_tables(conn):
    """Get all feature table names from the Variable table."""
    if not _table_exists(conn, 'Variable'):
        return []
    return [
        r[0] for r in conn.execute(
            "SELECT Feature_table_name FROM Variable WHERE Feature_table_name IS NOT NULL"
        ).fetchall()
    ]


def _table_has_column(conn, table_name, column_name):
    """Check if a table has a specific column."""
    cols = [r[0] for r in conn.execute(f'DESCRIBE "{table_name}"').fetchall()]
    return column_name in cols


def remove_sample(db_path, sample_name):
    """Remove a sample and all its associated data from the database."""
    conn = duckdb.connect(db_path)

    if not _table_exists(conn, 'Sample'):
        conn.close()
        print("No Sample table in database (genbank-only mode).")
        return

    row = conn.execute(
        "SELECT Sample_id FROM Sample WHERE Sample_name = ?", [sample_name]
    ).fetchone()
    if row is None:
        conn.close()
        print(f"Error: sample '{sample_name}' not found in database.")
        return
    sample_id = row[0]

    # Delete PhageTermini via PhageMechanisms packaging IDs
    if _table_exists(conn, 'PhageMechanisms') and _table_exists(conn, 'PhageTermini'):
        conn.execute(
            "DELETE FROM PhageTermini WHERE Packaging_id IN "
            "(SELECT Packaging_id FROM PhageMechanisms WHERE Sample_id = ?)",
            [sample_id],
        )

    # Delete from all fixed tables that reference Sample_id
    for table in [
        'PhageMechanisms', 'Coverage', 'Misassembly', 'Microdiversity',
        'Side_misassembly', 'Topology',
    ]:
        _delete_from(conn, table, 'Sample_id', sample_id)

    # Delete from dynamic feature tables
    for ft in _get_feature_tables(conn):
        if _table_exists(conn, ft) and _table_has_column(conn, ft, 'Sample_id'):
            _delete_from(conn, ft, 'Sample_id', sample_id)

    # Delete the sample row itself
    conn.execute("DELETE FROM Sample WHERE Sample_id = ?", [sample_id])
    update_database_metadata(conn)
    conn.close()
    print(f"Removed sample '{sample_name}' and all associated data.")


def remove_contig(db_path, contig_name):
    """Remove a contig and all its associated data from the database."""
    conn = duckdb.connect(db_path)

    row = conn.execute(
        "SELECT Contig_id FROM Contig WHERE Contig_name = ?", [contig_name]
    ).fetchone()
    if row is None:
        conn.close()
        print(f"Error: contig '{contig_name}' not found in database.")
        return
    contig_id = row[0]

    # Delete PhageTermini via PhageMechanisms packaging IDs
    if _table_exists(conn, 'PhageMechanisms') and _table_exists(conn, 'PhageTermini'):
        conn.execute(
            "DELETE FROM PhageTermini WHERE Packaging_id IN "
            "(SELECT Packaging_id FROM PhageMechanisms WHERE Contig_id = ?)",
            [contig_id],
        )

    # Delete from all fixed tables that reference Contig_id
    for table in [
        'PhageMechanisms', 'Coverage', 'Misassembly', 'Microdiversity',
        'Side_misassembly', 'Topology',
        'Contig_sequence', 'Contig_annotation',
        'Contig_directRepeats', 'Contig_invertedRepeats',
        'Contig_GCContent', 'Contig_GCSkew',
        'Contig_direct_repeat_count', 'Contig_inverted_repeat_count',
        'Contig_direct_repeat_identity', 'Contig_inverted_repeat_identity',
    ]:
        _delete_from(conn, table, 'Contig_id', contig_id)

    # Delete from dynamic feature tables
    for ft in _get_feature_tables(conn):
        if _table_exists(conn, ft) and _table_has_column(conn, ft, 'Contig_id'):
            _delete_from(conn, ft, 'Contig_id', contig_id)

    # Delete the contig row itself
    conn.execute("DELETE FROM Contig WHERE Contig_id = ?", [contig_id])
    update_database_metadata(conn)
    conn.close()
    print(f"Removed contig '{contig_name}' and all associated data.")


def main(argv=None):
    import argparse

    parser = argparse.ArgumentParser(prog="database_getters", description="Inspect database contents")
    sub = parser.add_subparsers(dest="cmd", required=True)

    sp = sub.add_parser('list-variables', help='List variables and metadata')
    sp.add_argument('-d', '--db', required=True)

    sp = sub.add_parser('list-samples', help='List samples')
    sp.add_argument('-d', '--db', required=True)

    sp = sub.add_parser('list-contigs', help='List contigs')
    sp.add_argument('-d', '--db', required=True)

    args = parser.parse_args(argv)
    if args.cmd == 'list-variables':
        list_variables(args.db)
    elif args.cmd == 'list-samples':
        list_samples(args.db)
    elif args.cmd == 'list-contigs':
        list_contigs(args.db)

if __name__ == '__main__':
    main()
