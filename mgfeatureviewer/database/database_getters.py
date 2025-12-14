import sqlite3

def list_variables(db_path, detailed=False):
    """Print variables and detailed metadata from Variable table (excluding Feature_table_name)."""
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    cur.execute("PRAGMA table_info(Variable)")
    cols = [r[1] for r in cur.fetchall()]

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
    conn = sqlite3.connect(db_path)
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
    conn = sqlite3.connect(db_path)
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
