import duckdb

# Columns to exclude from annotation filtering UI (internal/metadata columns)
ANNOTATION_EXCLUDED_COLUMNS = {'Contig_id', 'Start', 'End', 'Parent_annotation_id', 'Annotation_id', 'Segments', 'Nucleotide_sequence', 'Protein_sequence', 'S_sites', 'N_sites'}


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
    except (subprocess.CalledProcessError, OSError):
        pass

    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    conn.execute("UPDATE Database_metadata SET Value = ? WHERE Key = 'Date_of_last_modification'", [now])
    conn.execute("UPDATE Database_metadata SET Value = ? WHERE Key = 'Tool_version_used_for_last_modification'", [tool_version])


def get_view_mode(conn):
    """Return 'mag' or 'contig'. Defaults to 'contig' for legacy DBs."""
    try:
        row = conn.execute("SELECT Value FROM Database_metadata WHERE Key = 'View_mode'").fetchone()
        if row and row[0]:
            return row[0]
    except duckdb.Error:
        pass
    return 'contig'


def is_mag_mode(conn):
    return get_view_mode(conn) == 'mag'


def list_mags(conn):
    """Return list of MAG names (alphabetically sorted)."""
    try:
        rows = conn.execute("SELECT MAG_name FROM MAG ORDER BY MAG_name").fetchall()
    except duckdb.Error:
        return []
    return [r[0] for r in rows]


def get_mag_metadata(conn, mag_name):
    """Return a dict of column_name → value for the given MAG, or None if not found."""
    try:
        row = conn.execute(
            "SELECT * FROM MAG WHERE MAG_name = ?", [mag_name]
        ).fetchone()
        if row is None:
            return None
        cols = [d[0] for d in conn.description]
        return dict(zip(cols, row))
    except duckdb.Error:
        return None


def get_mag_contigs(conn, mag_name):
    """Return list of (contig_name, contig_length, cumulative_offset) for contigs
    in the MAG, ordered longest-first. Offset is read straight from
    MAG_contigs_association.Offset_in_MAG (written at ingest time).
    """
    try:
        rows = conn.execute(
            "SELECT c.Contig_name, c.Contig_length, mca.Offset_in_MAG "
            "FROM MAG_contigs_association mca "
            "JOIN MAG mg ON mg.MAG_id = mca.MAG_id "
            "JOIN Contig c ON c.Contig_id = mca.Contig_id "
            "WHERE mg.MAG_name = ? "
            "ORDER BY mca.Offset_in_MAG ASC",
            [mag_name],
        ).fetchall()
    except duckdb.Error:
        return []
    return [(name, int(length), int(offset)) for name, length, offset in rows]


def list_mag_samples(conn, mag_name):
    """Samples that have any Coverage row on at least one contig of the given MAG."""
    try:
        rows = conn.execute(
            "SELECT DISTINCT s.Sample_name "
            "FROM Coverage cov "
            "JOIN Sample s ON s.Sample_id = cov.Sample_id "
            "JOIN MAG_contigs_association mca ON mca.Contig_id = cov.Contig_id "
            "JOIN MAG mg ON mg.MAG_id = mca.MAG_id "
            "WHERE mg.MAG_name = ? "
            "ORDER BY s.Sample_name",
            [mag_name],
        ).fetchall()
    except duckdb.Error:
        return []
    return [r[0] for r in rows]


def get_mag_id(conn, mag_name):
    """Return MAG_id for the given MAG name, or None if not found."""
    try:
        row = conn.execute(
            "SELECT MAG_id FROM MAG WHERE MAG_name = ?", [mag_name]
        ).fetchone()
    except duckdb.Error:
        return None
    return int(row[0]) if row else None


def _feature_id_for(conn, feature_name):
    """Look up Variable_id for a feature name."""
    try:
        row = conn.execute(
            "SELECT Variable_id FROM Variable WHERE Variable_name = ?", [feature_name]
        ).fetchone()
    except duckdb.Error:
        return None
    return int(row[0]) if row else None


def get_mag_feature_zoom(cur, mag_id, sample_id, feature_name):
    """Zoom blob bytes for one (MAG, Sample, Feature), or None."""
    fid = _feature_id_for(cur, feature_name)
    if fid is None:
        return None
    try:
        row = cur.execute(
            "SELECT Zoom_data FROM MAG_blob "
            "WHERE MAG_id = ? AND Sample_id = ? AND Feature_id = ?",
            [mag_id, sample_id, fid],
        ).fetchone()
    except duckdb.Error:
        return None
    return bytes(row[0]) if row else None


def get_mag_feature_chunks(cur, mag_id, sample_id, feature_name, chunk_lo, chunk_hi):
    """Return [(chunk_idx, bytes)] for Chunk_idx in [chunk_lo, chunk_hi]."""
    fid = _feature_id_for(cur, feature_name)
    if fid is None:
        return []
    try:
        rows = cur.execute(
            "SELECT Chunk_idx, Data FROM MAG_blob_chunk "
            "WHERE MAG_id = ? AND Sample_id = ? AND Feature_id = ? "
            "AND Chunk_idx BETWEEN ? AND ? "
            "ORDER BY Chunk_idx",
            [mag_id, sample_id, fid, chunk_lo, chunk_hi],
        ).fetchall()
    except duckdb.Error:
        return []
    return [(int(idx), bytes(data)) for idx, data in rows]


def get_mag_contig_zoom(cur, mag_id, feature_name):
    """Zoom blob bytes for one (MAG, Feature) in MAG_contig_blob, or None."""
    fid = _feature_id_for(cur, feature_name)
    if fid is None:
        return None
    try:
        row = cur.execute(
            "SELECT Zoom_data FROM MAG_contig_blob "
            "WHERE MAG_id = ? AND Feature_id = ?",
            [mag_id, fid],
        ).fetchone()
    except duckdb.Error:
        return None
    return bytes(row[0]) if row else None


def get_mag_contig_chunks(cur, mag_id, feature_name, chunk_lo, chunk_hi):
    """Return [(chunk_idx, bytes)] for Chunk_idx in [chunk_lo, chunk_hi]."""
    fid = _feature_id_for(cur, feature_name)
    if fid is None:
        return []
    try:
        rows = cur.execute(
            "SELECT Chunk_idx, Data FROM MAG_contig_blob_chunk "
            "WHERE MAG_id = ? AND Feature_id = ? "
            "AND Chunk_idx BETWEEN ? AND ? "
            "ORDER BY Chunk_idx",
            [mag_id, fid, chunk_lo, chunk_hi],
        ).fetchall()
    except duckdb.Error:
        return []
    return [(int(idx), bytes(data)) for idx, data in rows]


def get_mag_contig_map(conn):
    """Return (mag_to_contigs, contig_to_mag). Empty dicts when not MAG-mode."""
    if not is_mag_mode(conn):
        return {}, {}
    try:
        rows = conn.execute(
            "SELECT mg.MAG_name, c.Contig_name, c.Contig_length "
            "FROM MAG_contigs_association mca "
            "JOIN MAG mg ON mg.MAG_id = mca.MAG_id "
            "JOIN Contig c ON c.Contig_id = mca.Contig_id "
            "ORDER BY mg.MAG_name, c.Contig_length DESC, c.Contig_id ASC"
        ).fetchall()
    except duckdb.Error:
        return {}, {}
    mag_to_contigs = {}
    contig_to_mag = {}
    for mag_name, contig_name, _length in rows:
        mag_to_contigs.setdefault(mag_name, []).append(contig_name)
        contig_to_mag[contig_name] = mag_name
    return mag_to_contigs, contig_to_mag


def get_filtering_metadata(db_path: str) -> dict:
    """
    Get column metadata for Filtering2 UI.

    Returns metadata structure with lazy-loaded distinct values.
    On startup, only column names and types are collected (fast DESCRIBE queries).
    Distinct values for text columns are fetched on-demand via resolve_distinct_values().

    Structure: {
        category: {
            'source': 'table_or_view_name',
            'columns': {
                col_name: {
                    'type': 'text' | 'numeric',
                    'distinct_values': None  # Populated lazily
                    'is_bool': True  # Only for boolean columns
                    'source': ...  # For qualifier-based columns
                    'qualifier_key': ...  # For qualifier-based columns
                }
            }
        }
    }
    """
    conn = duckdb.connect(db_path, read_only=True)
    has_mags = is_mag_mode(conn)

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

    if has_mags:
        category_config.update({
            'MAG': {
                'source': 'MAG',
                'exclude': ['MAG_id']
            },
            'MAG coverage': {
                'source': 'Explicit_coverage_per_MAG',
                'exclude': ['MAG_name', 'Sample_name']
            },
            'MAG misassembly': {
                'source': 'Explicit_misassembly_per_MAG',
                'exclude': ['MAG_name', 'Sample_name']
            },
            'MAG microdiversity': {
                'source': 'Explicit_microdiversity_per_MAG',
                'exclude': ['MAG_name', 'Sample_name']
            },
        })

    # Text type names in DuckDB
    text_types = {'VARCHAR', 'TEXT', 'STRING'}
    bool_types = {'BOOLEAN', 'BOOL'}

    result = {}
    for category, config in category_config.items():
        source = config['source']
        exclude = set(config['exclude'])

        # Check if table/view exists and get column info (DESCRIBE is instant)
        try:
            cols_info = conn.execute(f"DESCRIBE {source}").fetchall()
        except duckdb.Error:
            continue

        # Classify columns by type — no data queries needed
        columns = {}
        for col_name, col_type, *_ in cols_info:
            if col_name in exclude:
                continue
            is_bool = col_type.upper() in bool_types
            is_text = is_bool or any(t in col_type.upper() for t in text_types)

            col_data = {'type': 'text' if is_text else 'numeric'}
            if is_bool:
                col_data['distinct_values'] = ['yes', 'no']
                col_data['is_bool'] = True
            else:
                col_data['distinct_values'] = None  # Loaded lazily
            columns[col_name] = col_data

        if columns:
            result[category] = {
                'source': source,
                'columns': columns
            }

    # Add Contig_qualifier keys (cheap: just get key names, not values)
    if 'Contig' in result:
        try:
            cq_keys = [
                row[0] for row in conn.execute(
                    'SELECT DISTINCT "Key" FROM Contig_qualifier WHERE "Key" IS NOT NULL ORDER BY "Key"'
                ).fetchall()
            ]
        except duckdb.Error:
            cq_keys = []

        for key in cq_keys:
            if key in result['Contig']['columns']:
                continue
            result['Contig']['columns'][key] = {
                'type': 'text',
                'source': 'Contig_qualifier',
                'qualifier_key': key,
                'distinct_values': None,  # Loaded lazily
            }

    # Build Annotations category — schema only (DESCRIBE + key listing)
    annotations_columns = {}

    try:
        ann_cols_info = conn.execute("DESCRIBE Contig_annotation").fetchall()
    except duckdb.Error:
        ann_cols_info = []

    for col_name, col_type, *_ in ann_cols_info:
        if col_name in ANNOTATION_EXCLUDED_COLUMNS:
            continue
        is_text = any(t in col_type.upper() for t in text_types)
        annotations_columns[col_name] = {
            'type': 'text' if is_text else 'numeric',
            'source': 'Contig_annotation',
            'distinct_values': None,  # Loaded lazily
        }

    # Annotation_qualifier keys (just key names)
    pivoted_lower = {c.lower() for c in annotations_columns.keys()}
    try:
        aq_keys = [
            row[0] for row in conn.execute(
                'SELECT DISTINCT "Key" FROM Annotation_qualifier WHERE "Key" IS NOT NULL ORDER BY "Key"'
            ).fetchall()
        ]
    except duckdb.Error:
        aq_keys = []

    for key in aq_keys:
        if key.lower() in pivoted_lower or key in annotations_columns:
            continue
        annotations_columns[key] = {
            'type': 'text',
            'source': 'Annotation_qualifier',
            'qualifier_key': key,
            'distinct_values': None,  # Loaded lazily
        }

    if annotations_columns:
        # Insert Annotations right after Contig in the ordered dict
        new_result = {}
        for key, val in result.items():
            new_result[key] = val
            if key == 'Contig':
                new_result['Annotations'] = {
                    'source': 'Contig_annotation',
                    'columns': annotations_columns,
                }
        if 'Annotations' not in new_result:
            new_result['Annotations'] = {
                'source': 'Contig_annotation',
                'columns': annotations_columns,
            }
        result = new_result

    conn.close()
    return result


def resolve_distinct_values(db_path: str, filtering_metadata: dict,
                            category: str, col_name: str) -> list:
    """Fetch and cache distinct values for a text column on demand.

    First call for a given (category, col_name) runs the query; subsequent
    calls return the cached list.  Returns [] for numeric columns or on error.
    """
    cat_meta = filtering_metadata.get(category, {})
    col_info = cat_meta.get('columns', {}).get(col_name, {})
    if not col_info:
        return []

    # Already resolved
    if col_info.get('distinct_values') is not None:
        return col_info['distinct_values']

    # Numeric columns don't have distinct values
    if col_info.get('type') != 'text':
        return []

    conn = duckdb.connect(db_path, read_only=True)
    distinct_values = []

    try:
        source_override = col_info.get('source')
        qualifier_key = col_info.get('qualifier_key')

        if qualifier_key and source_override == 'Contig_qualifier':
            distinct = conn.execute(
                'SELECT DISTINCT "Value" FROM Contig_qualifier '
                'WHERE "Key" = ? AND "Value" IS NOT NULL ORDER BY "Value"',
                [qualifier_key]
            ).fetchall()
            distinct_values = [row[0] for row in distinct]
        elif qualifier_key and source_override == 'Annotation_qualifier':
            distinct = conn.execute(
                'SELECT DISTINCT "Value" FROM Annotation_qualifier '
                'WHERE "Key" = ? AND "Value" IS NOT NULL ORDER BY "Value"',
                [qualifier_key]
            ).fetchall()
            distinct_values = [row[0] for row in distinct]
        else:
            # Regular column from a table or view
            source = source_override or cat_meta.get('source', '')
            if source:
                distinct = conn.execute(
                    f'SELECT DISTINCT "{col_name}" FROM {source} '
                    f'WHERE "{col_name}" IS NOT NULL ORDER BY "{col_name}"'
                ).fetchall()
                distinct_values = [row[0] for row in distinct]
    except duckdb.Error:
        distinct_values = []
    finally:
        conn.close()

    # Cache result so we don't query again
    col_info['distinct_values'] = distinct_values
    return distinct_values


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

def list_mags_cli(db_path):
    """Print MAG_name values from MAG table."""
    conn = duckdb.connect(db_path, read_only=True)
    if conn.execute("SELECT 1 FROM information_schema.tables WHERE table_name = 'MAG'").fetchone() is None:
        print("No MAG table in database (not in MAG mode).")
        conn.close()
        return
    rows = [r[0] for r in conn.execute("SELECT MAG_name FROM MAG ORDER BY MAG_name").fetchall()]
    if not rows:
        print("No MAGs found in the database.")
    else:
        for m in rows:
            print(f"{m}")
    conn.close()

SAMPLE_INTERNAL_COLUMNS = {
    'Sample_id', 'Sample_name', 'Sequencing_type',
    'Number_of_reads', 'Number_of_mapped_reads',
}

CONTIG_INTERNAL_COLUMNS = {
    'Contig_id', 'Contig_name', 'Contig_length', 'Duplication_percentage', 
    'GC_mean', 'GC_sd', 'GC_skew_amplitude', 'Positive_GC_skew_windows_percentage',
    'Number_of_samples',
}

MAG_INTERNAL_COLUMNS = {
    'MAG_id', 'MAG_name', 'MAG_length', 'Number_of_contigs', 'N50', 'Duplication_percentage', 
    'GC_mean', 'GC_sd', 'GC_skew_amplitude', 'Positive_GC_skew_windows_percentage',
    'Number_of_samples',
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


def list_mag_metadata(db_path):
    """Print user-added column names on the MAG table."""
    conn = duckdb.connect(db_path, read_only=True)
    if conn.execute("SELECT 1 FROM information_schema.tables WHERE table_name = 'MAG'").fetchone() is None:
        conn.close()
        print("No MAG table in database (not in MAG mode).")
        return
    cols = [r[0] for r in conn.execute("DESCRIBE MAG").fetchall()]
    conn.close()
    user_cols = [c for c in cols if c not in MAG_INTERNAL_COLUMNS]
    if not user_cols:
        print("No user-added metadata columns on MAG table.")
    else:
        for c in user_cols:
            print(c)


def remove_mag_metadata(db_path, colname):
    """Remove a user-added column from the MAG table."""
    if colname in MAG_INTERNAL_COLUMNS:
        print(f"Error: '{colname}' is a built-in column and cannot be removed.")
        return
    conn = duckdb.connect(db_path)
    if conn.execute("SELECT 1 FROM information_schema.tables WHERE table_name = 'MAG'").fetchone() is None:
        conn.close()
        print("No MAG table in database (not in MAG mode).")
        return
    cols = [r[0] for r in conn.execute("DESCRIBE MAG").fetchall()]
    if colname not in cols:
        conn.close()
        print(f"Error: column '{colname}' not found in MAG table.")
        return
    conn.execute(f'ALTER TABLE MAG DROP COLUMN "{colname}"')
    update_database_metadata(conn)
    conn.close()
    print(f"Removed column '{colname}' from MAG table.")


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

    # Delete Phage_termini via Phage_mechanisms packaging IDs
    if _table_exists(conn, 'Phage_mechanisms') and _table_exists(conn, 'Phage_termini'):
        conn.execute(
            "DELETE FROM Phage_termini WHERE Packaging_id IN "
            "(SELECT Packaging_id FROM Phage_mechanisms WHERE Sample_id = ?)",
            [sample_id],
        )

    # Delete from all fixed tables that reference Sample_id
    for table in [
        'Phage_mechanisms', 'Coverage', 'Misassembly', 'Microdiversity',
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

    # Delete Phage_termini via Phage_mechanisms packaging IDs
    if _table_exists(conn, 'Phage_mechanisms') and _table_exists(conn, 'Phage_termini'):
        conn.execute(
            "DELETE FROM Phage_termini WHERE Packaging_id IN "
            "(SELECT Packaging_id FROM Phage_mechanisms WHERE Contig_id = ?)",
            [contig_id],
        )

    # Delete from all fixed tables that reference Contig_id
    for table in [
        'Phage_mechanisms', 'Coverage', 'Misassembly', 'Microdiversity',
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




def remove_mag(db_path, mag_name):
    """Remove a MAG and all its associated data from the database."""
    conn = duckdb.connect(db_path)

    if not _table_exists(conn, 'MAG'):
        conn.close()
        print("No MAG table in database (not in MAG mode).")
        return

    row = conn.execute(
        "SELECT MAG_id FROM MAG WHERE MAG_name = ?", [mag_name]
    ).fetchone()
    if row is None:
        conn.close()
        print(f"Error: MAG '{mag_name}' not found in database.")
        return
    mag_id = row[0]

    # Delete MAG-level blob data
    for table in [
        'MAG_blob', 'MAG_blob_chunk',
        'MAG_contig_blob', 'MAG_contig_blob_chunk',
        'MAG_coverage',
    ]:
        _delete_from(conn, table, 'MAG_id', mag_id)

    # Get contigs belonging to this MAG
    contig_ids = [
        r[0] for r in conn.execute(
            "SELECT Contig_id FROM MAG_contigs_association WHERE MAG_id = ?",
            [mag_id],
        ).fetchall()
    ]

    # Delete per-contig data for each member contig
    for contig_id in contig_ids:
        # Phage termini via Phage_mechanisms
        if _table_exists(conn, 'Phage_mechanisms') and _table_exists(conn, 'Phage_termini'):
            conn.execute(
                "DELETE FROM Phage_termini WHERE Packaging_id IN "
                "(SELECT Packaging_id FROM Phage_mechanisms WHERE Contig_id = ?)",
                [contig_id],
            )

        for table in [
            'Phage_mechanisms', 'Coverage', 'Misassembly', 'Microdiversity',
            'Side_misassembly', 'Topology',
            'Contig_sequence', 'Contig_annotation',
            'Contig_directRepeats', 'Contig_invertedRepeats',
            'Contig_GCContent', 'Contig_GCSkew',
            'Contig_direct_repeat_count', 'Contig_inverted_repeat_count',
            'Contig_direct_repeat_identity', 'Contig_inverted_repeat_identity',
        ]:
            _delete_from(conn, table, 'Contig_id', contig_id)

        # Dynamic feature tables
        for ft in _get_feature_tables(conn):
            if _table_exists(conn, ft) and _table_has_column(conn, ft, 'Contig_id'):
                _delete_from(conn, ft, 'Contig_id', contig_id)

        # Delete the contig row itself
        conn.execute("DELETE FROM Contig WHERE Contig_id = ?", [contig_id])

    # Delete association rows and the MAG row itself
    _delete_from(conn, 'MAG_contigs_association', 'MAG_id', mag_id)
    conn.execute("DELETE FROM MAG WHERE MAG_id = ?", [mag_id])
    update_database_metadata(conn)
    conn.close()
    print(f"Removed MAG '{mag_name}' and all associated data ({len(contig_ids)} contigs).")

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
