import sqlite3
import csv
import re


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
    parser.add_argument('db')
    parser.add_argument('variable_name')
    parser.add_argument('type')
    parser.add_argument('color')
    parser.add_argument('title')
    parser.add_argument('csv_file')


def run_add_variable(args):
    DB = args.db
    var_name = args.variable_name
    vtype = args.type
    color = args.color
    title = args.title
    csv_file = args.csv_file

    # --- Validate variable name ---
    if not re.match(r"^[A-Za-z_][A-Za-z0-9_]*$", var_name):
        print("ERROR: Invalid variable name. Use only letters, digits, and underscores, not starting with a digit.")
        return 1

    cfg = config_feature_subplot(var_name, "External", vtype, color, title)
    feature_table = f"Feature_{var_name}"
    subplot = cfg.get("subplot")
    alpha = cfg.get("alpha")
    fill_alpha = cfg.get("fill_alpha")
    size = cfg.get("size")

    # --- Connect ---
    conn = sqlite3.connect(DB)
    cur = conn.cursor()

    try:
        # --- Check variable existence ---
        cur.execute("SELECT Variable_id FROM Variable WHERE Variable_name=?", (var_name,))
        if cur.fetchone():
            raise ValueError(f"ERROR: Variable '{var_name}' already exists in the database.")

        # --- Read known samples and contigs ---
        cur.execute("SELECT Sample_name, Sample_id FROM Sample")
        samples = dict(cur.fetchall())
        cur.execute("SELECT Contig_name, Contig_id, Contig_length FROM Contig")
        contigs_info = {row[0]: (row[1], row[2]) for row in cur.fetchall()}

        # --- Read CSV and validate ---
        rows_to_insert = []
        presences_validated = []
        absences_validated = []
        line_no = 0

        with open(csv_file, newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                line_no += 1
                try:
                    sample = row["sample"].strip()
                    contig = row["contig"].strip()
                    first_pos = int(row["first_position"])
                    last_pos = int(row["last_position"])
                    value = float(row["value"])
                except KeyError:
                    raise ValueError(f"ERROR: CSV missing required columns (sample,contig,first_position,last_position,value)")
                except ValueError:
                    raise ValueError(f"ERROR: Invalid number at line {line_no}: {row}")

                # --- Check sample/contig existence ---
                if sample not in samples:
                    raise ValueError(f"ERROR: Sample '{sample}' not found in database (line {line_no})")
                if contig not in contigs_info:
                    raise ValueError(f"ERROR: Contig '{contig}' not found in database (line {line_no})")

                contig_id, contig_length = contigs_info[contig]
                if first_pos > last_pos:
                    raise ValueError(f"ERROR: first_position > last_position at line {line_no}")
                if first_pos < 0 or last_pos > contig_length:
                    raise ValueError(f"ERROR: Position out of bounds for contig '{contig}' (length {contig_length}) at line {line_no}")

                # --- Check sample/contig pair is present ---
                sample_id = samples[sample]
                if (contig_id, sample_id) not in presences_validated and (contig_id, sample_id) not in absences_validated:
                    cur.execute("SELECT 1 FROM Presences WHERE Contig_id=? AND Sample_id=? LIMIT 1", (contig_id, sample_id))
                    exists = cur.fetchone() is not None

                    if not exists:
                        absences_validated.append((contig_id, sample_id))
                        print(f"WARNING: No presence record for contig '{contig}' and sample '{sample}'. Associated data was not written into the database")
                    else:
                        presences_validated.append((contig_id, sample_id))

                if (contig_id, sample_id) in presences_validated:
                    rows_to_insert.append((contig_id, sample_id, first_pos, last_pos, value))

        # --- Check for overlapping ranges within each contig/sample pair ---
        # Group rows by (contig_id, sample_id)
        groups = {}
        for contig_id, sample_id, first_pos, last_pos, value in rows_to_insert:
            key = (contig_id, sample_id)
            if key not in groups:
                groups[key] = []
            groups[key].append((first_pos, last_pos))
        
        # Check for overlaps in each group
        for (contig_id, sample_id), ranges in groups.items():
            # Sort ranges by start position
            sorted_ranges = sorted(ranges)
            for i in range(len(sorted_ranges) - 1):
                curr_start, curr_end = sorted_ranges[i]
                next_start, next_end = sorted_ranges[i + 1]
                # Check if ranges overlap: current range ends after next range starts
                if curr_end >= next_start:
                    # Find the contig and sample names for better error message
                    contig_name = [name for name, (cid, _) in contigs_info.items() if cid == contig_id][0]
                    sample_name = [name for name, sid in samples.items() if sid == sample_id][0]
                    raise ValueError(
                        f"ERROR: Overlapping ranges detected for contig '{contig_name}' and sample '{sample_name}': "
                        f"[{curr_start}, {curr_end}] overlaps with [{next_start}, {next_end}]"
                    )

        # --- Create new variable ---
        cur.execute("""
        INSERT INTO Variable (Variable_name, Subplot, Module, Type, Color, Alpha, Fill_alpha, Size, Title, Help, Feature_table_name)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, (var_name, subplot, "External", vtype, color, float(alpha), float(fill_alpha), float(size), title, "", feature_table))

        # --- Create associated feature table ---
        cur.execute(f"""
        CREATE TABLE {feature_table} (
            Feature_id INTEGER PRIMARY KEY AUTOINCREMENT,
            Contig_id INTEGER,
            Sample_id INTEGER,
            First_position INTEGER,
            Last_position INTEGER,
            Value REAL,
            FOREIGN KEY(Contig_id) REFERENCES Contig(Contig_id),
            FOREIGN KEY(Sample_id) REFERENCES Sample(Sample_id)
        );
        """)

        # --- Insert data ---
        if rows_to_insert:
            cur.executemany(f"INSERT INTO {feature_table} (Contig_id, Sample_id, First_position, Last_position, Value) VALUES (?, ?, ?, ?, ?)", rows_to_insert)

        # --- Get Variable_id for the newly created variable ---
        cur.execute("SELECT Variable_id FROM Variable WHERE Variable_name=?", (var_name,))
        variable_id = cur.fetchone()[0]

        # --- Update Summary table with row counts per contig/sample ---
        if rows_to_insert:
            # Count rows per (contig_id, sample_id) pair
            summary_counts = {}
            for contig_id, sample_id, first_pos, last_pos, value in rows_to_insert:
                key = (contig_id, sample_id)
                summary_counts[key] = summary_counts.get(key, 0) + 1
            
            # Insert into Summary table
            summary_rows = [(contig_id, sample_id, variable_id, count) 
                          for (contig_id, sample_id), count in summary_counts.items()]
            cur.executemany(
                "INSERT INTO Summary (Contig_id, Sample_id, Variable_id, Row_count) VALUES (?, ?, ?, ?)",
                summary_rows
            )

        conn.commit()
        print(f"Variable '{var_name}' added and {len(rows_to_insert)} records inserted into '{feature_table}'")

    except Exception as e:
        print(str(e))
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