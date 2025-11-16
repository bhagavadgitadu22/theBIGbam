import sqlite3
import sys

# Define subplots characteristics
def config_feature_subplot(subplot, module, plot_type, color, title, alpha=0.8, fill_alpha=0.4, size=1, help=""):    
    if plot_type == "bars":
        alpha = 0.5
        size = 5
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

def main():
    if len(sys.argv) != 2:
        print("Usage: python generate_database.py <database_name.db>")
        sys.exit(1)

    DB = sys.argv[1]
    conn = sqlite3.connect(DB)
    cur = conn.cursor()

    # --- Create Contig and Sample tables ---
    cur.execute("""
    CREATE TABLE Contig (
        Contig_id INTEGER PRIMARY KEY AUTOINCREMENT,
        Contig_name TEXT UNIQUE,
        Contig_length INTEGER,
        Annotation_tool TEXT
    );
    """)

    cur.execute("""
    CREATE TABLE Sample (
        Sample_id INTEGER PRIMARY KEY AUTOINCREMENT,
        Sample_name TEXT UNIQUE
    );
    """)

    cur.execute("""
    CREATE TABLE Presences (
        Presence_id INTEGER PRIMARY KEY AUTOINCREMENT,
        Contig_id INTEGER,
        Sample_id INTEGER,
        Coverage_percentage REAL,
        FOREIGN KEY(Contig_id) REFERENCES Contig(Contig_id),
        FOREIGN KEY(Sample_id) REFERENCES Sample(Sample_id)
    );
    """)

    # --- Create Contig tables ---
    cur.execute("""
    CREATE TABLE Contig_annotation (
        Contig_annotation_id INTEGER PRIMARY KEY AUTOINCREMENT,
        Contig_id INTEGER,
        Start INTEGER,
        End INTEGER,
        Strand INTEGER,
        Type TEXT,
        Product TEXT,
        Function TEXT,
        Phrog INTEGER,
        FOREIGN KEY(Contig_id) REFERENCES Contig(Contig_id)
    );
    """)

    # --- Create Variable table ---
    cur.execute("""
    CREATE TABLE Variable (
        Variable_id INTEGER PRIMARY KEY AUTOINCREMENT,
        Variable_name TEXT UNIQUE,
        Subplot TEXT,
        Module TEXT,
        Type TEXT,
        Color TEXT,
        Alpha REAL,
        Fill_alpha REAL,
        Size REAL,
        Title TEXT,
        Help TEXT,
        Feature_table_name TEXT
    );
    """)

    # Populate Variable table from constants.FEATURE_SUBPLOTS
    FEATURE_SUBPLOTS = {
        "coverage": config_feature_subplot("Coverage", "Coverage", "curve", "#333333", "Coverage depth"),

        # Starts subplots
        "coverage_reduced": config_feature_subplot("Coverage reduced", "Phage termini", "curve", "#00c53b", "Coverage reduced", 
                                                   help = "Coverage depth calculated only from reads that start with a match to the reference (and end with a match for long reads)"),
        "reads_starts": config_feature_subplot("Reads termini", "Phage termini", "bars", "#215732", "Read Starts"),
        "reads_ends": config_feature_subplot("Reads termini", "Phage termini", "bars", "#6cc24a", "Read Ends"),
        "tau": config_feature_subplot("Tau", "Phage termini", "bars", "#44883e", "Tau"),

        # Misassembly subplots
        "read_lengths": config_feature_subplot("Read lengths", "Assembly check", "curve", "#ed8b00", "Read Lengths"),
        "insert_sizes": config_feature_subplot("Insert sizes", "Assembly check", "curve", "#ed8b00", "Insert Sizes"),
        "bad_orientations": config_feature_subplot("Bad orientations", "Assembly check", "bars", "#c94009", "Bad Orientations"),
        "left_clippings": config_feature_subplot("Clippings", "Assembly check", "bars", "#7f0091", "Left Clippings"),
        "right_clippings": config_feature_subplot("Clippings", "Assembly check", "bars", "#8e43e7", "Right Clippings",),
        "insertions": config_feature_subplot("Indels", "Assembly check", "bars", "#e50001", "Insertions"),
        "deletions": config_feature_subplot("Indels", "Assembly check", "bars", "#97011a", "Deletions"),
        "mismatches": config_feature_subplot("Mismatches", "Assembly check", "bars", "#5a0f0b", "Mismatches"),
    }

    for var_name, cfg in FEATURE_SUBPLOTS.items():
        feature_table = f"Feature_{var_name}"
        subplot = cfg.get("subplot")
        module = cfg.get("module")
        vtype = cfg.get("type")
        color = cfg.get("color")
        alpha = cfg.get("alpha")
        fill_alpha = cfg.get("fill_alpha")
        size = cfg.get("size")
        title = cfg.get("title")
        help = cfg.get("help")

        cur.execute("INSERT OR IGNORE INTO Variable (Variable_name, Subplot, Module, Type, Color, Alpha, Fill_alpha, Size, Title, Help, Feature_table_name) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                    (var_name, subplot, module, vtype, color, alpha, fill_alpha, size, title, help, feature_table))
        
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

if __name__ == "__main__":
    main()