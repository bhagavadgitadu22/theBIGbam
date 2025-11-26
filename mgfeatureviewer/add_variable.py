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
            next(reader)
            for row in reader:
                line_no += 1
                try:
                    sample = row["sample"].strip()
                    contig = row["contig"].strip()
                    pos = int(row["position"])
                    value = float(row["value"])
                except KeyError:
                    raise ValueError(f"ERROR: CSV missing required columns (sample, contig, position, value)")
                except ValueError:
                    raise ValueError(f"ERROR: Invalid number at line {line_no}: {row}")

                # --- Check sample/contig existence ---
                if sample not in samples:
                    raise ValueError(f"ERROR: Sample '{sample}' not found in database (line {line_no})")
                if contig not in contigs_info:
                    raise ValueError(f"ERROR: Contig '{contig}' not found in database (line {line_no})")

                contig_id, contig_length = contigs_info[contig]
                if pos < 0 or pos > contig_length:
                    raise ValueError(f"ERROR: Position {pos} out of range for contig '{contig}' (length={contig_length}) at line {line_no}")

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
                    rows_to_insert.append((contig_id, sample_id, pos, value))

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
            Position INTEGER,
            Value REAL,
            FOREIGN KEY(Contig_id) REFERENCES Contig(Contig_id),
            FOREIGN KEY(Sample_id) REFERENCES Sample(Sample_id)
        );
        """)

        # --- Insert data ---
        if rows_to_insert:
            cur.executemany(f"INSERT INTO {feature_table} (Contig_id, Sample_id, Position, Value) VALUES (?, ?, ?, ?)", rows_to_insert)

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