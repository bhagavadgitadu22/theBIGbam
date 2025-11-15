import sqlite3
import sys
import csv
import re

from generate_database import config_feature_subplot 

# --- Args ---
if len(sys.argv) != 7:
    print("Usage: python add_variable.py <database.db> <variable_name> <type> <color> <title> <data.csv>")
    sys.exit(1)

DB, var_name, vtype, color, title, csv_file = sys.argv[1:9]

# --- Validate variable name ---
if not re.match(r"^[A-Za-z_][A-Za-z0-9_]*$", var_name):
    print("ERROR: Invalid variable name. Use only letters, digits, and underscores, not starting with a digit.")
    sys.exit(1)

cfg = config_feature_subplot("External", vtype, color, title)
feature_table = f"Feature_{var_name}"
module = cfg.get("module")
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
        print(f"WARNING: Variable '{var_name}' already exists in the database.")
        conn.close()
        sys.exit(0)

    # --- Read known samples and contigs ---
    cur.execute("SELECT Sample_name, Sample_id FROM Sample")
    samples = dict(cur.fetchall())
    cur.execute("SELECT Contig_name, Contig_id, Contig_length FROM Contig")
    contigs_info = {row[0]: (row[1], row[2]) for row in cur.fetchall()}

    # --- Create new variable ---
    cur.execute("""
    INSERT INTO Variable (Variable_name, Module, Type, Color, Alpha, Fill_alpha, Size, Title, Feature_table_name)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, (var_name, "External", vtype, color, float(alpha), float(fill_alpha), float(size), title, feature_table))

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

    # --- Read CSV and validate ---
    rows_to_insert = []
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
                conn.rollback()
                raise ValueError(f"ERROR: CSV missing required columns (sample, contig, position, value)")
            except ValueError:
                conn.rollback()
                raise ValueError(f"ERROR: Invalid number at line {line_no}: {row}")

            # --- Check sample/contig existence ---
            if sample not in samples:
                conn.rollback()
                raise ValueError(f"ERROR: Sample '{sample}' not found in database (line {line_no})")
            if contig not in contigs_info:
                conn.rollback()
                raise ValueError(f"ERROR: Contig '{contig}' not found in database (line {line_no})")

            contig_id, contig_length = contigs_info[contig]
            if pos < 0 or pos > contig_length:
                conn.rollback()
                raise ValueError(
                    f"ERROR: Position {pos} out of range for contig '{contig}' (length={contig_length}) at line {line_no}"
                )

            rows_to_insert.append((contig_id, samples[sample], pos, value))

    # --- Insert data ---
    cur.executemany(
        f"INSERT INTO {feature_table} (Contig_id, Sample_id, Position, Value) VALUES (?, ?, ?, ?)",
        rows_to_insert,
    )

    conn.commit()
    print(f"Variable '{var_name}' added and {len(rows_to_insert)} records inserted into '{feature_table}'")

except Exception as e:
    print(str(e))
    print("Error encountered — reverting changes...")
    conn.rollback()

    # Try to drop newly created feature table if it exists
    try:
        cur.execute(f"DROP TABLE IF EXISTS {feature_table};")
        cur.execute("DELETE FROM Variable WHERE Variable_name=?", (var_name,))
        conn.commit()
        print(f"Rolled back and removed partially created variable '{var_name}' and table '{feature_table}'")
    except Exception as cleanup_error:
        print(f"WARNING: Cleanup error: {cleanup_error}")

finally:
    conn.close()