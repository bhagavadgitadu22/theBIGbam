import time
from dataclasses import dataclass

import duckdb


@dataclass(frozen=True)
class HistogramResult:
    """Plain histogram data, including whether it represents a bounded sample."""

    edges: object
    counts: object
    sampled_rows: int
    approximate: bool

    @property
    def bin_count(self) -> int:
        return len(self.counts)

# Columns to exclude from annotation filtering UI (internal/metadata columns)
ANNOTATION_EXCLUDED_COLUMNS = {'Contig_id', 'Start', 'End', 'Parent_annotation_id', 'Annotation_id', 'Segments', 'Nucleotide_sequence', 'Protein_sequence', 'S_sites', 'N_sites'}


def update_database_metadata(conn):
    """Update Date_of_last_modification and Tool_version_used_for_last_modification."""
    import subprocess
    from datetime import datetime
    from importlib.metadata import version

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



def get_mag_members(conn, mag_id):
    """Return [(contig_id, contig_length, offset)] for a MAG, ordered by offset."""
    return [(cid, clen, off) for cid, _name, clen, off in get_mag_members_full(conn, mag_id)]


def get_mag_members_full(conn, mag_id):
    """Return [(contig_id, contig_name, contig_length, offset)] for a MAG, ordered by offset."""
    try:
        rows = conn.execute(
            "SELECT mca.Contig_id, c.Contig_name, c.Contig_length, mca.Offset_in_MAG "
            "FROM MAG_contigs_association mca "
            "JOIN Contig c ON c.Contig_id = mca.Contig_id "
            "WHERE mca.MAG_id = ? ORDER BY mca.Offset_in_MAG",
            [mag_id],
        ).fetchall()
    except duckdb.Error:
        return []
    return [(int(cid), name, int(clen), int(off)) for cid, name, clen, off in rows]


def _validate_column(conn, table, column):
    """Check that column exists in table to prevent SQL injection."""
    try:
        cols = [r[0] for r in conn.execute(f"DESCRIBE {table}").fetchall()]
    except duckdb.Error:
        return False
    return column in cols


def _sorted_members_query(conn, mag_id, source_table, metric, ascending, sample_name=None):
    """Shared implementation for sorted member queries.

    Returns [(contig_id, contig_name, contig_length)] in the requested order.
    """
    if not _validate_column(conn, source_table, metric):
        raise ValueError(f"Column {metric!r} not found in {source_table}")
    direction = "ASC" if ascending else "DESC"
    if source_table == "Contig":
        rows = conn.execute(
            f'SELECT mca.Contig_id, c.Contig_name, c.Contig_length '
            f'FROM MAG_contigs_association mca '
            f'JOIN Contig c ON c.Contig_id = mca.Contig_id '
            f'WHERE mca.MAG_id = ? '
            f'ORDER BY c."{metric}" {direction}',
            [mag_id],
        ).fetchall()
    else:
        if sample_name is None:
            raise ValueError("sample_name required for sample-dependent sort")
        rows = conn.execute(
            f'SELECT mca.Contig_id, c.Contig_name, c.Contig_length '
            f'FROM MAG_contigs_association mca '
            f'JOIN Contig c ON c.Contig_id = mca.Contig_id '
            f'LEFT JOIN {source_table} sv ON sv.Contig_name = c.Contig_name '
            f'  AND sv.Sample_name = ? '
            f'WHERE mca.MAG_id = ? '
            f'ORDER BY COALESCE(sv."{metric}", 0) {direction}',
            [sample_name, mag_id],
        ).fetchall()
    return rows


def get_mag_members_sorted(conn, mag_id, source_table, metric, ascending, sample_name=None):
    """Return [(contig_id, contig_length, offset)] sorted by metric with recomputed offsets."""
    return [
        (contig_id, contig_length, offset)
        for contig_id, _contig_name, contig_length, offset in get_mag_members_full_sorted(
            conn, mag_id, source_table, metric, ascending, sample_name
        )
    ]


def get_mag_contigs_sorted(conn, mag_id, source_table, metric, ascending, sample_name=None):
    """Return [(contig_name, contig_length, offset)] sorted by metric with recomputed offsets."""
    return [
        (contig_name, contig_length, offset)
        for _contig_id, contig_name, contig_length, offset in get_mag_members_full_sorted(
            conn, mag_id, source_table, metric, ascending, sample_name
        )
    ]


def get_mag_members_full_sorted(conn, mag_id, source_table, metric, ascending, sample_name=None):
    """Return sorted members with names and IDs from one database query."""
    rows = _sorted_members_query(conn, mag_id, source_table, metric, ascending, sample_name)
    offset = 0
    result = []
    for contig_id, contig_name, contig_length in rows:
        length = int(contig_length)
        result.append((int(contig_id), str(contig_name), length, offset))
        offset += length
    return result





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


def get_filtering_metadata(db_path: str, enable_timing: bool = False) -> dict:
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
    t0 = time.perf_counter()
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
            'exclude': ['Contig_id', 'Contig_name', 'Sample_id', 'Sample_name']
        },
        'Misassembly': {
            'source': 'Explicit_misassembly',
            'exclude': ['Contig_id', 'Contig_name', 'Sample_id', 'Sample_name']
        },
        'Microdiversity': {
            'source': 'Explicit_microdiversity',
            'exclude': ['Contig_id', 'Contig_name', 'Sample_id', 'Sample_name']
        },
        'Side misassembly': {
            'source': 'Explicit_side_misassembly',
            'exclude': ['Contig_id', 'Contig_name', 'Sample_id', 'Sample_name']
        },
        'Topology': {
            'source': 'Explicit_topology',
            'exclude': ['Contig_id', 'Contig_name', 'Sample_id', 'Sample_name']
        },
        'Termini': {
            'source': 'Explicit_phage_mechanisms',
            'exclude': ['Contig_id', 'Contig_name', 'Sample_id', 'Sample_name']
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
                'exclude': ['MAG_id', 'MAG_name', 'Sample_id', 'Sample_name']
            },
            'MAG misassembly': {
                'source': 'Explicit_misassembly_per_MAG',
                'exclude': ['MAG_id', 'MAG_name', 'Sample_id', 'Sample_name']
            },
            'MAG microdiversity': {
                'source': 'Explicit_microdiversity_per_MAG',
                'exclude': ['MAG_id', 'MAG_name', 'Sample_id', 'Sample_name']
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
    if enable_timing:
        print(f"[timing] get_filtering_metadata: {time.perf_counter() - t0:.3f}s", flush=True)
    return result


def resolve_distinct_values(db_path: str, filtering_metadata: dict,
                            category: str, col_name: str, enable_timing: bool = False) -> list:
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

    t0 = time.perf_counter()
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
    if enable_timing:
        print(f"[timing] resolve_distinct_values({category}.{col_name}): {time.perf_counter() - t0:.3f}s ({len(distinct_values)} values)", flush=True)
    return distinct_values


def search_distinct_values(
    db_path: str,
    filtering_metadata: dict,
    category: str,
    col_name: str,
    search_term: str = "",
    limit: int = 100,
) -> list:
    """Return a bounded text-value search without materializing the full cardinality."""
    cat_meta = filtering_metadata.get(category, {})
    col_info = cat_meta.get('columns', {}).get(col_name, {})
    if col_info.get('type') != 'text':
        return []
    limit = max(1, min(int(limit), 500))
    term = (search_term or '').strip()
    source_override = col_info.get('source')
    qualifier_key = col_info.get('qualifier_key')
    conn = duckdb.connect(db_path, read_only=True)
    try:
        if qualifier_key and source_override in {'Contig_qualifier', 'Annotation_qualifier'}:
            source = source_override
            value_sql = 'CAST("Value" AS VARCHAR)'
            base_sql = f'SELECT DISTINCT "Value" FROM {source} WHERE "Key" = ? AND "Value" IS NOT NULL'
            base_params = [qualifier_key]
        else:
            source = source_override or cat_meta.get('source', '')
            if not source:
                return []
            value_sql = f'CAST("{col_name}" AS VARCHAR)'
            base_sql = f'SELECT DISTINCT "{col_name}" FROM {source} WHERE "{col_name}" IS NOT NULL'
            base_params = []

        if not term:
            rows = conn.execute(f'{base_sql} LIMIT ?', [*base_params, limit]).fetchall()
            return sorted((row[0] for row in rows), key=lambda value: str(value).casefold())

        # Run increasingly broad bounded searches. Exact and prefix matches
        # normally fill a useful result window without scanning all contains
        # matches or globally sorting a high-cardinality column.
        results = []
        seen = set()
        searches = (
            (f'{value_sql} = ?', term),
            (f'{value_sql} ILIKE ? || \'%\'', term),
            (f'{value_sql} ILIKE \'%\' || ? || \'%\'', term),
        )
        for predicate, parameter in searches:
            if len(results) >= limit:
                break
            rows = conn.execute(
                f'{base_sql} AND {predicate} LIMIT ?',
                [*base_params, parameter, limit],
            ).fetchall()
            for row in rows:
                value = row[0]
                identity = str(value)
                if identity not in seen:
                    seen.add(identity)
                    results.append(value)
                    if len(results) == limit:
                        break
        return sorted(results, key=lambda value: str(value).casefold())
    except duckdb.Error:
        return []
    finally:
        conn.close()


def list_variables(db_path, detailed=False):
    """Print variables and detailed metadata from Variable table (excluding Feature_table_name)."""
    conn = duckdb.connect(db_path, read_only=True)
    cur = conn.cursor()

    # DuckDB: get column names from DESCRIBE statement
    cur.execute("DESCRIBE Variable")
    cols = [r[0] for r in cur.fetchall()]

    # Fields we will display (exclude Feature_table_name)
    display_fields = [c for c in cols if c not in ['Variable_id', 'Feature_table_name']]

    # Query the table for the display fields
    sel = ", ".join(display_fields)
    cur.execute(f"SELECT {sel} FROM Variable ORDER BY Variable_name")
    rows = cur.fetchall()

    if not rows:
        print("No variables found in the database.", flush=True)
        conn.close()
        return

    # Print header
    for row in rows:
        # Pair field name and value and print nicely
        print(row[0], flush=True)  # Variable_name as header
        if detailed:
            # Print other fields minus variable name
            for fname, val in zip(display_fields[1:], row[1:]):
                print(f"- {fname}: {val}", flush=True)
            print("", flush=True)

    conn.close()

def list_samples(db_path):
    """Print Sample_name values from Sample table."""
    conn = duckdb.connect(db_path, read_only=True)
    cur = conn.cursor()
    cur.execute("SELECT 1 FROM information_schema.tables WHERE table_name = 'Sample'")
    if cur.fetchone() is None:
        print("No samples in the database (genbank-only mode).", flush=True)
        conn.close()
        return
    cur.execute("SELECT Sample_name FROM Sample ORDER BY Sample_name")
    rows = [r[0] for r in cur.fetchall()]
    if not rows:
        print("No samples found in the database.", flush=True)
    else:
        for s in rows:
            print(f"{s}", flush=True)
    conn.close()

def list_contigs(db_path):
    """Print Contig_name values from Contig table."""
    conn = duckdb.connect(db_path, read_only=True)
    cur = conn.cursor()
    cur.execute("SELECT Contig_name FROM Contig ORDER BY Contig_name")
    rows = [r[0] for r in cur.fetchall()]
    if not rows:
        print("No contigs found in the database.", flush=True)
    else:
        for c in rows:
            print(f"{c}", flush=True)
    conn.close()

def list_mags_cli(db_path):
    """Print MAG_name values from MAG table."""
    conn = duckdb.connect(db_path, read_only=True)
    if conn.execute("SELECT 1 FROM information_schema.tables WHERE table_name = 'MAG'").fetchone() is None:
        print("No MAG table in database (not in MAG mode).", flush=True)
        conn.close()
        return
    rows = [r[0] for r in conn.execute("SELECT MAG_name FROM MAG ORDER BY MAG_name").fetchall()]
    if not rows:
        print("No MAGs found in the database.", flush=True)
    else:
        for m in rows:
            print(f"{m}", flush=True)
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
        print("No Sample table in database (genbank-only mode).", flush=True)
        return
    cols = [r[0] for r in conn.execute("DESCRIBE Sample").fetchall()]
    conn.close()
    user_cols = [c for c in cols if c not in SAMPLE_INTERNAL_COLUMNS]
    if not user_cols:
        print("No user-added metadata columns on Sample table.", flush=True)
    else:
        for c in user_cols:
            print(c, flush=True)


def list_contig_metadata(db_path):
    """Print user-added column names on the Contig table."""
    conn = duckdb.connect(db_path, read_only=True)
    cols = [r[0] for r in conn.execute("DESCRIBE Contig").fetchall()]
    conn.close()
    user_cols = [c for c in cols if c not in CONTIG_INTERNAL_COLUMNS]
    if not user_cols:
        print("No user-added metadata columns on Contig table.", flush=True)
    else:
        for c in user_cols:
            print(c, flush=True)


def remove_sample_metadata(db_path, colname):
    """Remove a user-added column from the Sample table."""
    if colname in SAMPLE_INTERNAL_COLUMNS:
        print(f"Error: '{colname}' is a built-in column and cannot be removed.", flush=True)
        return
    conn = duckdb.connect(db_path)
    if conn.execute("SELECT 1 FROM information_schema.tables WHERE table_name = 'Sample'").fetchone() is None:
        conn.close()
        print("No Sample table in database (genbank-only mode).", flush=True)
        return
    cols = [r[0] for r in conn.execute("DESCRIBE Sample").fetchall()]
    if colname not in cols:
        conn.close()
        print(f"Error: column '{colname}' does not exist on Sample table.", flush=True)
        return
    conn.execute(f'ALTER TABLE Sample DROP COLUMN "{colname}"')
    update_database_metadata(conn)
    conn.close()
    print(f"Removed column '{colname}' from Sample table.", flush=True)


def remove_contig_metadata(db_path, colname):
    """Remove a user-added column from the Contig table."""
    if colname in CONTIG_INTERNAL_COLUMNS:
        print(f"Error: '{colname}' is a built-in column and cannot be removed.", flush=True)
        return
    conn = duckdb.connect(db_path)
    cols = [r[0] for r in conn.execute("DESCRIBE Contig").fetchall()]
    if colname not in cols:
        conn.close()
        print(f"Error: column '{colname}' does not exist on Contig table.", flush=True)
        return
    conn.execute(f'ALTER TABLE Contig DROP COLUMN "{colname}"')
    update_database_metadata(conn)
    conn.close()
    print(f"Removed column '{colname}' from Contig table.", flush=True)


def list_mag_metadata(db_path):
    """Print user-added column names on the MAG table."""
    conn = duckdb.connect(db_path, read_only=True)
    if conn.execute("SELECT 1 FROM information_schema.tables WHERE table_name = 'MAG'").fetchone() is None:
        conn.close()
        print("No MAG table in database (not in MAG mode).", flush=True)
        return
    cols = [r[0] for r in conn.execute("DESCRIBE MAG").fetchall()]
    conn.close()
    user_cols = [c for c in cols if c not in MAG_INTERNAL_COLUMNS]
    if not user_cols:
        print("No user-added metadata columns on MAG table.", flush=True)
    else:
        for c in user_cols:
            print(c, flush=True)


def remove_mag_metadata(db_path, colname):
    """Remove a user-added column from the MAG table."""
    if colname in MAG_INTERNAL_COLUMNS:
        print(f"Error: '{colname}' is a built-in column and cannot be removed.", flush=True)
        return
    conn = duckdb.connect(db_path)
    if conn.execute("SELECT 1 FROM information_schema.tables WHERE table_name = 'MAG'").fetchone() is None:
        conn.close()
        print("No MAG table in database (not in MAG mode).", flush=True)
        return
    cols = [r[0] for r in conn.execute("DESCRIBE MAG").fetchall()]
    if colname not in cols:
        conn.close()
        print(f"Error: column '{colname}' not found in MAG table.", flush=True)
        return
    conn.execute(f'ALTER TABLE MAG DROP COLUMN "{colname}"')
    update_database_metadata(conn)
    conn.close()
    print(f"Removed column '{colname}' from MAG table.", flush=True)


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


def _recompute_sample_counts(conn):
    """Recompute Number_of_samples on Contig (and MAG if in MAG mode)."""
    if _table_exists(conn, 'Contig') and _table_has_column(conn, 'Contig', 'Number_of_samples'):
        conn.execute("""
            UPDATE Contig SET Number_of_samples = (
                SELECT COUNT(DISTINCT Sample_id) FROM Coverage
                WHERE Coverage.Contig_id = Contig.Contig_id)
        """)
    if _table_exists(conn, 'MAG') and _table_has_column(conn, 'MAG', 'Number_of_samples'):
        conn.execute("""
            UPDATE MAG SET Number_of_samples = (
                SELECT COUNT(DISTINCT c.Sample_id)
                FROM Coverage c
                JOIN MAG_contigs_association mca ON mca.Contig_id = c.Contig_id
                WHERE mca.MAG_id = MAG.MAG_id)
        """)


def _recompute_mapped_reads(conn):
    """Recompute Number_of_mapped_reads from remaining Coverage entries."""
    if _table_exists(conn, 'Sample') and _table_exists(conn, 'Coverage'):
        conn.execute("""
            UPDATE Sample SET Number_of_mapped_reads = (
                SELECT COALESCE(SUM(Read_count), 0) FROM Coverage
                WHERE Coverage.Sample_id = Sample.Sample_id)
        """)


def _delete_contig_annotations(conn, contig_id):
    """Delete annotation data for a contig (core table + sub-tables)."""
    if not _table_exists(conn, 'Contig_annotation_core'):
        return
    anno_ids = [r[0] for r in conn.execute(
        "SELECT Annotation_id FROM Contig_annotation_core WHERE Contig_id = ?", [contig_id]
    ).fetchall()]
    if anno_ids:
        placeholders = ','.join('?' * len(anno_ids))
        for table in ['Annotation_segments', 'Annotation_qualifier', 'Annotation_sequence']:
            if _table_exists(conn, table):
                conn.execute(f'DELETE FROM "{table}" WHERE Annotation_id IN ({placeholders})', anno_ids)
    _delete_from(conn, 'Contig_qualifier', 'Contig_id', contig_id)
    _delete_from(conn, 'Contig_annotation_core', 'Contig_id', contig_id)


def remove_sample(db_path, sample_name):
    """Remove a sample and all its associated data from the database."""
    conn = duckdb.connect(db_path)

    if not _table_exists(conn, 'Sample'):
        conn.close()
        print("No Sample table in database (genbank-only mode).", flush=True)
        return

    row = conn.execute(
        "SELECT Sample_id FROM Sample WHERE Sample_name = ?", [sample_name]
    ).fetchone()
    if row is None:
        conn.close()
        print(f"Error: sample '{sample_name}' not found in database.", flush=True)
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
        'Feature_blob', 'Feature_blob_chunk',
        'MAG_coverage',
    ]:
        _delete_from(conn, table, 'Sample_id', sample_id)

    # Delete from dynamic feature tables
    for ft in _get_feature_tables(conn):
        if _table_exists(conn, ft) and _table_has_column(conn, ft, 'Sample_id'):
            _delete_from(conn, ft, 'Sample_id', sample_id)

    # Delete the sample row itself
    conn.execute("DELETE FROM Sample WHERE Sample_id = ?", [sample_id])

    _recompute_sample_counts(conn)
    update_database_metadata(conn)
    conn.close()
    print(f"Removed sample '{sample_name}' and all associated data.", flush=True)


def remove_contig(db_path, contig_name):
    """Remove a contig and all its associated data from the database (contig-mode only)."""
    conn = duckdb.connect(db_path)

    if _table_exists(conn, 'MAG'):
        conn.close()
        print("Error: cannot remove individual contigs in a MAG-mode database. "
              "Use remove-mag instead.", flush=True)
        return

    row = conn.execute(
        "SELECT Contig_id FROM Contig WHERE Contig_name = ?", [contig_name]
    ).fetchone()
    if row is None:
        conn.close()
        print(f"Error: contig '{contig_name}' not found in database.", flush=True)
        return
    contig_id = row[0]

    # Delete Phage_termini via Phage_mechanisms packaging IDs
    if _table_exists(conn, 'Phage_mechanisms') and _table_exists(conn, 'Phage_termini'):
        conn.execute(
            "DELETE FROM Phage_termini WHERE Packaging_id IN "
            "(SELECT Packaging_id FROM Phage_mechanisms WHERE Contig_id = ?)",
            [contig_id],
        )

    # Delete annotation sub-tables (Contig_annotation is a VIEW, not a table)
    _delete_contig_annotations(conn, contig_id)

    # Delete from all fixed tables that reference Contig_id
    for table in [
        'Phage_mechanisms', 'Coverage', 'Misassembly', 'Microdiversity',
        'Side_misassembly', 'Topology',
        'Contig_sequence',
        'Contig_directRepeats', 'Contig_invertedRepeats',
        'Contig_GCContent', 'Contig_GCSkew',
        'Contig_direct_repeat_count', 'Contig_inverted_repeat_count',
        'Contig_direct_repeat_identity', 'Contig_inverted_repeat_identity',
        'Contig_blob', 'Contig_blob_chunk',
        'Feature_blob', 'Feature_blob_chunk',
    ]:
        _delete_from(conn, table, 'Contig_id', contig_id)

    # Delete from dynamic feature tables
    for ft in _get_feature_tables(conn):
        if _table_exists(conn, ft) and _table_has_column(conn, ft, 'Contig_id'):
            _delete_from(conn, ft, 'Contig_id', contig_id)

    # Delete the contig row itself
    conn.execute("DELETE FROM Contig WHERE Contig_id = ?", [contig_id])

    _recompute_mapped_reads(conn)
    update_database_metadata(conn)
    conn.close()
    print(f"Removed contig '{contig_name}' and all associated data.", flush=True)




def remove_mag(db_path, mag_name):
    """Remove a MAG and all its associated data from the database."""
    conn = duckdb.connect(db_path)

    if not _table_exists(conn, 'MAG'):
        conn.close()
        print("No MAG table in database (not in MAG mode).", flush=True)
        return

    row = conn.execute(
        "SELECT MAG_id FROM MAG WHERE MAG_name = ?", [mag_name]
    ).fetchone()
    if row is None:
        conn.close()
        print(f"Error: MAG '{mag_name}' not found in database.", flush=True)
        return
    mag_id = row[0]

    # Delete MAG-level data
    for table in [
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

        # Annotation sub-tables (Contig_annotation is a VIEW, not a table)
        _delete_contig_annotations(conn, contig_id)

        for table in [
            'Phage_mechanisms', 'Coverage', 'Misassembly', 'Microdiversity',
            'Side_misassembly', 'Topology',
            'Contig_sequence',
            'Contig_directRepeats', 'Contig_invertedRepeats',
            'Contig_GCContent', 'Contig_GCSkew',
            'Contig_direct_repeat_count', 'Contig_inverted_repeat_count',
            'Contig_direct_repeat_identity', 'Contig_inverted_repeat_identity',
            'Contig_blob', 'Contig_blob_chunk',
            'Feature_blob', 'Feature_blob_chunk',
        ]:
            _delete_from(conn, table, 'Contig_id', contig_id)

        # Inter-contig BLAST hits (dual-column foreign key)
        if _table_exists(conn, 'Contig_blast_hits'):
            conn.execute(
                "DELETE FROM Contig_blast_hits WHERE Contig_id_1 = ? OR Contig_id_2 = ?",
                [contig_id, contig_id],
            )

        # Dynamic feature tables
        for ft in _get_feature_tables(conn):
            if _table_exists(conn, ft) and _table_has_column(conn, ft, 'Contig_id'):
                _delete_from(conn, ft, 'Contig_id', contig_id)

        # Delete the contig row itself
        conn.execute("DELETE FROM Contig WHERE Contig_id = ?", [contig_id])

    # Delete association rows and the MAG row itself
    _delete_from(conn, 'MAG_contigs_association', 'MAG_id', mag_id)
    conn.execute("DELETE FROM MAG WHERE MAG_id = ?", [mag_id])

    _recompute_mapped_reads(conn)
    _recompute_sample_counts(conn)
    update_database_metadata(conn)
    conn.close()
    print(f"Removed MAG '{mag_name}' and all associated data ({len(contig_ids)} contigs).", flush=True)

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

def resolve_histogram_bins(db_path: str, filtering_metadata: dict,
                           category: str, col_name: str,
                           n_bins: int = 50, log_mode: bool = False,
                           scale: float | None = None, enable_timing: bool = False) -> HistogramResult | None:
    """Compute a bounded histogram while evaluating its source relation once."""
    import numpy as np
    cat_meta = filtering_metadata.get(category, {})
    col_info = cat_meta.get('columns', {}).get(col_name, {})
    if not col_info or col_info.get('type') != 'numeric' or col_info.get('is_bool'):
        return None

    source = col_info.get('source') or cat_meta.get('source', '')
    if not source:
        return None
    direct_sources = {
        'Explicit_coverage': 'Coverage',
        'Explicit_misassembly': 'Misassembly',
        'Explicit_microdiversity': 'Microdiversity',
        'Explicit_side_misassembly': 'Side_misassembly',
        'Explicit_topology': 'Topology',
        'Explicit_phage_mechanisms': 'Phage_mechanisms',
        'Explicit_coverage_per_MAG': 'MAG_coverage',
    }
    direct_source = direct_sources.get(source)

    if scale:
        val_expr = f'("{col_name}" / {scale})'
    else:
        val_expr = f'"{col_name}"'

    null_filter = f'"{col_name}" IS NOT NULL'
    if log_mode:
        if scale:
            val_expr = f'LOG10("{col_name}" / {scale})'
        else:
            val_expr = f'LOG10("{col_name}")'
        null_filter += f' AND "{col_name}" > 0'

    t0 = time.perf_counter()
    conn = duckdb.connect(db_path, read_only=True)
    try:
        if direct_source and _validate_column(conn, direct_source, col_name):
            source = direct_source

        plan = conn.execute(f'EXPLAIN SELECT "{col_name}" FROM {source}').fetchall()
        import re
        estimates = [int(v.replace(',', '')) for row in plan for v in re.findall(r'~([\d,]+) rows', str(row))]
        estimated_rows = max(estimates, default=0)
        approximate = estimated_rows > 250_000
        sample_clause = ''
        if approximate:
            # Oversample blocks before the hard row cap. This avoids empty
            # SYSTEM samples near the threshold while retaining bounded data.
            sample_percent = min(100.0, max(0.01, 250_000 / estimated_rows * 100))
            sample_clause = (
                f' USING SAMPLE system({sample_percent:.8f} PERCENT) REPEATABLE(42) LIMIT 100000'
            )
        rows = conn.execute(
            f'''WITH vals AS MATERIALIZED (
                    SELECT {val_expr} AS value FROM {source} WHERE {null_filter}{sample_clause}
                ), stats AS (
                    SELECT MIN(value) mn, MAX(value) mx, COUNT(*) total FROM vals
                ), buckets AS (
                    SELECT CASE WHEN s.mn = s.mx THEN 0 ELSE LEAST(
                        CAST(FLOOR((v.value - s.mn) / ((s.mx - s.mn) / LEAST({int(n_bins)}, s.total))) AS INTEGER),
                        CAST(LEAST({int(n_bins)}, s.total) - 1 AS INTEGER)
                    ) END bucket, s.mn, s.mx, s.total
                    FROM vals v CROSS JOIN stats s
                )
                SELECT bucket, COUNT(*), MIN(mn), MIN(mx), MIN(total)
                FROM buckets GROUP BY bucket ORDER BY bucket'''
        ).fetchall()
        if not rows or not rows[0][4]:
            return None
        mn, mx, total = float(rows[0][2]), float(rows[0][3]), int(rows[0][4])
        actual_bins = 1 if mn == mx else min(n_bins, total)
        edges = np.array([mn - 0.5, mx + 0.5]) if mn == mx else np.linspace(mn, mx, actual_bins + 1)
        counts = np.zeros(actual_bins, dtype=int)
        for bucket_num, count, *_ in rows:
            if 0 <= int(bucket_num) < actual_bins:
                counts[int(bucket_num)] = int(count)
        result = HistogramResult(edges, counts, total, approximate)
    except duckdb.Error as e:
        print(f"[resolve_histogram_bins] {e}", flush=True)
        result = None
    finally:
        conn.close()

    if enable_timing:
        mode = 'log' if log_mode else 'linear'
        strategy = 'sampled' if result and result.approximate else 'exact'
        bins = result.bin_count if result else 0
        sampled = result.sampled_rows if result else 0
        print(f"[timing] resolve_histogram_bins({category}.{col_name}, {mode}): "
              f"{time.perf_counter() - t0:.3f}s ({strategy}, {bins} bins, {sampled} rows)", flush=True)
    return result


def resolve_value_counts(db_path: str, filtering_metadata: dict,
                          category: str, col_name: str, enable_timing: bool = False) -> list | None:
    """Fetch and cache (value, count) pairs for a text column (for treemap display)."""
    cat_meta = filtering_metadata.get(category, {})
    col_info = cat_meta.get('columns', {}).get(col_name, {})
    if not col_info or col_info.get('type') != 'text':
        return None
    if 'value_counts' in col_info:
        return col_info['value_counts']

    source_override = col_info.get('source')
    qualifier_key = col_info.get('qualifier_key')
    t0 = time.perf_counter()
    conn = duckdb.connect(db_path, read_only=True)
    try:
        if qualifier_key and source_override in ('Contig_qualifier', 'Annotation_qualifier'):
            rows = conn.execute(
                f'SELECT "Value", COUNT(*) FROM {source_override} '
                'WHERE "Key" = ? AND "Value" IS NOT NULL '
                'GROUP BY "Value" ORDER BY COUNT(*) DESC',
                [qualifier_key]
            ).fetchall()
        else:
            source = source_override or cat_meta.get('source', '')
            if not source:
                col_info['value_counts'] = None
                return None
            rows = conn.execute(
                f'SELECT "{col_name}", COUNT(*) FROM {source} '
                f'WHERE "{col_name}" IS NOT NULL '
                f'GROUP BY "{col_name}" ORDER BY COUNT(*) DESC'
            ).fetchall()
        counts = [(r[0], r[1]) for r in rows]
    except duckdb.Error:
        counts = None
    finally:
        conn.close()

    col_info['value_counts'] = counts
    if enable_timing:
        n_vals = len(counts) if counts else 0
        print(f"[timing] resolve_value_counts({category}.{col_name}): {time.perf_counter() - t0:.3f}s ({n_vals} values)", flush=True)
    return counts

def resolve_column_null_stats(
    db_path: str,
    filtering_metadata: dict,
    category: str,
    col_name: str,
) -> tuple[int, int] | None:
    """Return (non_null_count, total_possible) for a column.

    ``total_possible`` is the total number of entities in the category
    (contigs, samples, MAGs, contig/sample pairs, etc.).
    For qualifier columns this is the row count of the parent table,
    not just rows that have that qualifier key.

    Queries the underlying table once and caches the result in
    ``filtering_metadata[category]['columns'][col_name]['null_stats']`` so
    repeated calls (e.g. log-toggle rebuilds) are free.

    Returns ``None`` when the column or its source cannot be resolved.
    """
    cat_meta = filtering_metadata.get(category, {})
    col_info = cat_meta.get("columns", {}).get(col_name, {})
    if not col_info:
        return None
    if "null_stats" in col_info:
        return col_info["null_stats"]

    source_override = col_info.get("source")
    qualifier_key = col_info.get("qualifier_key")
    source = source_override or cat_meta.get("source", "")

    conn = duckdb.connect(db_path, read_only=True)
    try:
        if qualifier_key and source_override in ("Contig_qualifier", "Annotation_qualifier"):
            row = conn.execute(
                f'SELECT COUNT(*) FILTER (WHERE "Value" IS NOT NULL) '
                f'FROM {source_override} WHERE "Key" = ?',
                [qualifier_key],
            ).fetchone()
            non_null = int(row[0]) if row else 0
            parent = "Contig" if source_override == "Contig_qualifier" else "Contig_annotation"
            total_row = conn.execute(f'SELECT COUNT(*) FROM {parent}').fetchone()
            total = int(total_row[0]) if total_row else 0
        elif source:
            row = conn.execute(
                f'SELECT '
                f'  COUNT(*) FILTER (WHERE "{col_name}" IS NOT NULL), '
                f'  COUNT(*) '
                f'FROM {source}'
            ).fetchone()
            non_null = int(row[0]) if row else 0
            total = int(row[1]) if row else 0
        else:
            col_info["null_stats"] = None
            return None
        result: tuple[int, int] | None = (non_null, total)
    except duckdb.Error as exc:
        print(f"[resolve_column_null_stats] {exc}", flush=True)
        result = None
    finally:
        conn.close()

    col_info["null_stats"] = result
    return result
