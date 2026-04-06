import duckdb
import csv
import re
import traceback

from thebigbam.database.database_getters import update_database_metadata


def config_feature_subplot(subplot, module, plot_type, color, title, alpha=0.8, fill_alpha=0.4, size=1, help=""):
    """Configure feature subplot parameters."""
    if plot_type == "bars":
        alpha = 0.6
        size = 1
    return {
        "subplot": subplot,
        "module": module,
        "type": plot_type,
        "color": color,
        "alpha": alpha,
        "fill_alpha": fill_alpha,
        "size": size,
        "title": title,
        "help": help
    }


def add_add_variable_args(parser):
    parser.add_argument('--db', required=True, help='Path to DuckDB database')
    parser.add_argument('--name', dest='variable_name', required=True, help='Name for the new variable')
    parser.add_argument('--type', dest='type', required=True, help='Plot type (bars or curve)')
    parser.add_argument('--color', required=True, help='Color for the variable (hex code: example #FF0000 for red)')
    parser.add_argument('--title', required=True, help='Title for the variable')
    parser.add_argument('--csv', dest='csv_file', required=True, help='CSV file with variable data. Required columns Contig,Sample,First_position,Last_position,Value. If all values for Sample are empty a Contig feature will be added instead of Contig/Sample pair feature.')
    parser.add_argument('--force', action='store_true', default=False, help='Remove existing variable with the same name before adding')


def run_add_variable(args):
    DB = args.db
    var_name = args.variable_name
    vtype = args.type
    color = args.color
    title = args.title
    csv_file = args.csv_file

    # --- Validate plot type ---
    if vtype not in ("bars", "curve"):
        print("ERROR: Invalid plot type. Must be 'bars' or 'curve'.")
        return 1

    # --- Validate variable name ---
    if not re.match(r"^[A-Za-z_][A-Za-z0-9_]*$", var_name):
        print("ERROR: Invalid variable name. Use only letters, digits, and underscores, not starting with a digit.")
        return 1

    cfg = config_feature_subplot(var_name, "Custom", vtype, color, title)
    subplot = cfg.get("subplot")
    alpha = cfg.get("alpha")
    fill_alpha = cfg.get("fill_alpha")
    size = cfg.get("size")

    # --- Connect ---
    conn = duckdb.connect(DB)
    cur = conn.cursor()

    try:
        # --- Check variable existence ---
        cur.execute("SELECT Variable_id, Feature_table_name FROM Variable WHERE Variable_name=?", (var_name,))
        existing = cur.fetchone()
        if existing:
            if args.force:
                existing_id, existing_table = existing
                cur.execute(f"DROP TABLE IF EXISTS {existing_table}")
                cur.execute("DELETE FROM Variable WHERE Variable_id=?", (existing_id,))
                print(f"Removed existing variable '{var_name}' (--force)")
            else:
                raise ValueError(f"ERROR: Variable '{var_name}' already exists in the database. Use --force to replace it.")

        # --- Read known samples and contigs ---
        cur.execute("SELECT 1 FROM information_schema.tables WHERE table_name = 'Sample'")
        if cur.fetchone() is not None:
            cur.execute("SELECT Sample_name, Sample_id FROM Sample")
            samples = dict(cur.fetchall())
        else:
            samples = {}
        cur.execute("SELECT Contig_name, Contig_id, Contig_length FROM Contig")
        contigs_info = {row[0]: (row[1], row[2]) for row in cur.fetchall()}

        # --- Pre-read CSV and detect contig-only mode ---
        with open(csv_file, newline="") as f:
            reader = csv.DictReader(f)
            raw_rows = list(reader)

        sample_values = [r.get("Sample", "").strip() for r in raw_rows]
        all_samples_empty = all(s == "" for s in sample_values)
        any_samples_empty = any(s == "" for s in sample_values)
        if any_samples_empty and not all_samples_empty:
            raise ValueError("ERROR: Mixed Sample values — either fill Sample for all rows or leave empty for all.")
        contig_only = all_samples_empty
        if not contig_only and not samples:
            raise ValueError("ERROR: CSV has Sample values but database has no samples (genbank-only mode). Leave the Sample column empty for contig-level variables.")

        # Set feature table name based on mode
        feature_table = f"Contig_{var_name}" if contig_only else f"Feature_{var_name}"

        # --- Validate rows ---
        rows_to_insert = []
        presences_validated = []
        absences_validated = []

        for line_no, row in enumerate(raw_rows, start=1):
            try:
                sample = row["Sample"].strip()
                contig = row["Contig"].strip()
                first_pos = int(row["First_position"])
                last_pos = int(row["Last_position"])
                value = float(row["Value"])
            except KeyError:
                raise ValueError("ERROR: CSV missing required columns (Contig,Sample,First_position,Last_position,Value)")
            except ValueError:
                raise ValueError(f"ERROR: Invalid number at line {line_no}: {row}")

            # --- Check contig existence ---
            if contig not in contigs_info:
                raise ValueError(f"ERROR: Contig '{contig}' not found in database (line {line_no})")

            contig_id, contig_length = contigs_info[contig]
            if first_pos > last_pos:
                raise ValueError(f"ERROR: first_position > last_position at line {line_no}")
            if first_pos < 0 or last_pos > contig_length:
                raise ValueError(f"ERROR: Position out of bounds for contig '{contig}' (length {contig_length}) at line {line_no}")

            if contig_only:
                rows_to_insert.append((contig_id, first_pos, last_pos, value))
            else:
                # --- Check sample existence ---
                if sample not in samples:
                    raise ValueError(f"ERROR: Sample '{sample}' not found in database (line {line_no})")

                # --- Check sample/contig pair is present ---
                sample_id = samples[sample]
                if (contig_id, sample_id) not in presences_validated and (contig_id, sample_id) not in absences_validated:
                    cur.execute("SELECT 1 FROM Coverage WHERE Contig_id=? AND Sample_id=? LIMIT 1", (contig_id, sample_id))
                    exists = cur.fetchone() is not None

                    if not exists:
                        absences_validated.append((contig_id, sample_id))
                        print(f"WARNING: No presence record for contig '{contig}' and sample '{sample}'. Associated data was not written into the database")
                    else:
                        presences_validated.append((contig_id, sample_id))

                if (contig_id, sample_id) in presences_validated:
                    rows_to_insert.append((contig_id, sample_id, first_pos, last_pos, value))

        # --- Check for overlapping ranges ---
        if contig_only:
            # Group rows by (contig_id,)
            groups = {}
            for contig_id, first_pos, last_pos, value in rows_to_insert:
                key = (contig_id,)
                if key not in groups:
                    groups[key] = []
                groups[key].append((first_pos, last_pos))

            for (contig_id,), ranges in groups.items():
                sorted_ranges = sorted(ranges)
                for i in range(len(sorted_ranges) - 1):
                    curr_start, curr_end = sorted_ranges[i]
                    next_start, next_end = sorted_ranges[i + 1]
                    if curr_end >= next_start:
                        contig_name = [name for name, (cid, _) in contigs_info.items() if cid == contig_id][0]
                        raise ValueError(
                            f"ERROR: Overlapping ranges detected for contig '{contig_name}': "
                            f"[{curr_start}, {curr_end}] overlaps with [{next_start}, {next_end}]"
                        )
        else:
            # Group rows by (contig_id, sample_id)
            groups = {}
            for contig_id, sample_id, first_pos, last_pos, value in rows_to_insert:
                key = (contig_id, sample_id)
                if key not in groups:
                    groups[key] = []
                groups[key].append((first_pos, last_pos))

            for (contig_id, sample_id), ranges in groups.items():
                sorted_ranges = sorted(ranges)
                for i in range(len(sorted_ranges) - 1):
                    curr_start, curr_end = sorted_ranges[i]
                    next_start, next_end = sorted_ranges[i + 1]
                    if curr_end >= next_start:
                        contig_name = [name for name, (cid, _) in contigs_info.items() if cid == contig_id][0]
                        sample_name = [name for name, sid in samples.items() if sid == sample_id][0]
                        raise ValueError(
                            f"ERROR: Overlapping ranges detected for contig '{contig_name}' and sample '{sample_name}': "
                            f"[{curr_start}, {curr_end}] overlaps with [{next_start}, {next_end}]"
                        )

        # --- Create new variable ---
        # Get the next Variable_id (table uses explicit integer PK, not auto-increment)
        cur.execute("SELECT COALESCE(MAX(Variable_id), 0) + 1 FROM Variable")
        next_var_id = cur.fetchone()[0]
        # Get the next Module_order for the Custom module
        cur.execute("SELECT COALESCE(MAX(Module_order), 0) + 1 FROM Variable WHERE Module='Custom'")
        next_module_order = cur.fetchone()[0]
        cur.execute("""
        INSERT INTO Variable (Variable_id, Variable_name, Subplot, Module, Module_order, "Type", Color, Alpha, Fill_alpha, "Size", Title, Help, Feature_table_name)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, (next_var_id, var_name, subplot, "Custom", next_module_order, vtype, color, float(alpha), float(fill_alpha), float(size), title, "", feature_table))

        # --- Create associated feature table ---
        if contig_only:
            cur.execute(f"""
            CREATE TABLE {feature_table} (
                Contig_id INTEGER,
                First_position INTEGER,
                Last_position INTEGER,
                Value REAL
            );
            """)
        else:
            cur.execute(f"""
            CREATE TABLE {feature_table} (
                Contig_id INTEGER,
                Sample_id INTEGER,
                First_position INTEGER,
                Last_position INTEGER,
                Value REAL
            );
            """)

        # --- Insert data ---
        if contig_only:
            if rows_to_insert:
                cur.executemany(f"INSERT INTO {feature_table} (Contig_id, First_position, Last_position, Value) VALUES (?, ?, ?, ?)", rows_to_insert)
        else:
            if rows_to_insert:
                cur.executemany(f"INSERT INTO {feature_table} (Contig_id, Sample_id, First_position, Last_position, Value) VALUES (?, ?, ?, ?, ?)", rows_to_insert)

        update_database_metadata(conn)
        conn.commit()
        print(f"Variable '{var_name}' added and {len(rows_to_insert)} records inserted into '{feature_table}'")

    except Exception as e:
        traceback.print_exc()
        return 1

    finally:
        conn.close()

    return 0


def add_remove_variable_args(parser):
    parser.add_argument('--db', required=True, help='Path to DuckDB database')
    parser.add_argument('--name', dest='variable_name', required=True, help='Name of the variable to remove')


def run_remove_variable(args):
    DB = args.db
    var_name = args.variable_name

    conn = duckdb.connect(DB)
    cur = conn.cursor()

    try:
        # Check that the variable exists
        cur.execute(
            "SELECT Variable_id, Feature_table_name, Module FROM Variable WHERE Variable_name=?",
            (var_name,)
        )
        row = cur.fetchone()
        if not row:
            print(f"Error: variable '{var_name}' does not exist in the database.")
            return 1

        variable_id, feature_table_name, module = row

        if module != 'Custom':
            print(f"Error: variable '{var_name}' is a built-in variable and cannot be removed.")
            return 1

        # Drop the associated feature table
        cur.execute(f"DROP TABLE IF EXISTS {feature_table_name}")

        # Remove the Variable metadata row
        cur.execute("DELETE FROM Variable WHERE Variable_id=?", (variable_id,))

        update_database_metadata(conn)
        conn.commit()
        print(f"Variable '{var_name}' removed from the database.")

    except Exception as e:
        traceback.print_exc()
        return 1

    finally:
        conn.close()

    return 0


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    add_add_variable_args(parser)
    args = parser.parse_args()
    raise SystemExit(run_add_variable(args))