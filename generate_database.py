import sqlite3
import sys

if len(sys.argv) != 2:
    print("Usage: python generate_database.py <database_name.db>")
    sys.exit(1)

DB = sys.argv[1]
conn = sqlite3.connect(DB)
cur = conn.cursor()

# --- Drop old tables if they exist ---
cur.executescript("""
DROP TABLE IF EXISTS Contig;
DROP TABLE IF EXISTS Sample;
DROP TABLE IF EXISTS Variable;
""")

# --- Create Contig and Sample tables ---
cur.execute("""
CREATE TABLE Contig (
    Contig_id INTEGER PRIMARY KEY AUTOINCREMENT,
    Contig_name TEXT UNIQUE,
    Contig_length INTEGER
);
""")

cur.execute("""
CREATE TABLE Sample (
    Sample_id INTEGER PRIMARY KEY AUTOINCREMENT,
    Sample_name TEXT UNIQUE
);
""")

# --- Create Variable table ---
cur.execute("""
CREATE TABLE Variable (
    Variable_id INTEGER PRIMARY KEY AUTOINCREMENT,
    Variable_name TEXT UNIQUE,
    Status TEXT,
    Type TEXT,
    Color TEXT,
    Alpha REAL,
    Size REAL,
    Title TEXT,
    Feature_table_name TEXT
);
""")

# Populate Variable table from constants.FEATURE_SUBPLOTS
# Define subplots characteristics
def config_feature_subplot(plot_type, color, title, alpha=0.7, size=1):
    if plot_type == "bars":
        alpha = 0.5
        size = 3
    return {
        "type": plot_type,
        "color": color,
        "alpha": alpha,
        "size": size,
        "title": title
    }

FEATURE_SUBPLOTS = {
    "coverage": config_feature_subplot("curve", "black", "Coverage Depth"),

    # Starts subplots
    "coverage_reduced": config_feature_subplot("curve", "black", "Coverage Depth (only reads starting and ending with a match)"),
    "reads_starts": config_feature_subplot("bars", "blue", "Reads' Starts"),
    "reads_ends": config_feature_subplot("bars", "blue", "Reads' Ends"),
    "tau": config_feature_subplot("bars", "blue", "Tau"),

    # Misassembly subplots
    "read_lengths": config_feature_subplot("curve", "green", "Read Lengths"),
    "insert_sizes": config_feature_subplot("curve", "green", "Insert Sizes"),
    "bad_orientations": config_feature_subplot("bars", "green", "Bad Orientations"),
    "left_clippings": config_feature_subplot("bars", "purple", "Left Clippings"),
    "right_clippings": config_feature_subplot("bars", "purple", "Right Clippings"),
    "insertions": config_feature_subplot("bars", "red", "Insertions"),
    "deletions": config_feature_subplot("bars", "red", "Deletions"),
    "mismatches": config_feature_subplot("bars", "red", "Mismatches"),
}

for var_name, cfg in FEATURE_SUBPLOTS.items():
    feature_table = f"Feature_{var_name}"
    status = "internal"
    vtype = cfg.get("type")
    color = cfg.get("color")
    alpha = cfg.get("alpha")
    size = cfg.get("size")
    title = cfg.get("title")
    cur.execute("INSERT OR IGNORE INTO Variable (Variable_name, Status, Type, Color, Alpha, Size, Title, Feature_table_name) VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                (var_name, status, vtype, color, alpha, size, title, feature_table))
    
    # Create the feature table for this variable
    cur.execute(f"""
    CREATE TABLE IF NOT EXISTS {feature_table} (
        Feature_id INTEGER PRIMARY KEY AUTOINCREMENT,
        Contig_id INTEGER,
        Sample_id INTEGER,
        Position INTEGER,
        Value REAL,
        FOREIGN KEY(Contig_id) REFERENCES Contig(Contig_id),
        FOREIGN KEY(Sample_id) REFERENCES Sample(Sample_id)
    );
    """)

conn.commit()
conn.close()
print(f"Empty database created and populated with variables: {DB}")