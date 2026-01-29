import duckdb


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
            'exclude': ['Contig_id', 'Contig_name', 'Annotation_tool']
        },
        'Sample': {
            'source': 'Sample',
            'exclude': ['Sample_id', 'Sample_name']
        },
        'Presences': {
            'source': 'Explicit_presences',
            'exclude': ['Contig_name', 'Sample_name', 'Read_number_correction_ratio', 'Read_mapped_correction_ratio']
        },
        'Completeness': {
            'source': 'Explicit_completeness',
            'exclude': ['Contig_name', 'Sample_name', 'Score_completeness', 'Score_contamination']
        },
        'Termini': {
            'source': 'Explicit_phage_mechanisms',
            'exclude': ['Contig_name', 'Sample_name']
        }
    }

    # Text type names in DuckDB
    text_types = {'VARCHAR', 'TEXT', 'STRING'}

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

            is_text = any(t in col_type.upper() for t in text_types)
            col_data = {'type': 'text' if is_text else 'numeric'}

            # For text columns, get distinct values
            if is_text:
                try:
                    distinct = conn.execute(
                        f"SELECT DISTINCT \"{col_name}\" FROM {source} WHERE \"{col_name}\" IS NOT NULL ORDER BY \"{col_name}\""
                    ).fetchall()
                    col_data['distinct_values'] = [row[0] for row in distinct]
                except Exception:
                    col_data['distinct_values'] = []

            columns[col_name] = col_data

        if columns:  # Only add category if it has columns
            result[category] = {
                'source': source,
                'columns': columns
            }

    # Add annotation columns to Contig category from Contig_annotation table
    if 'Contig' in result:
        annotation_exclude = {'Contig_id', 'Start', 'End', 'Strand', 'Type'}
        try:
            ann_cols_info = conn.execute("DESCRIBE Contig_annotation").fetchall()
            for col_name, col_type, *_ in ann_cols_info:
                if col_name in annotation_exclude:
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
