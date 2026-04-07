"""Functions for downloading data as CSV files.

This module contains functions for exporting data from the theBIGbam database
to CSV files. Returns CSV content as strings for browser download.
"""


def make_safe_filename(name):
    """Sanitize a string for use in filenames.
    
    Args:
        name: String to sanitize
        
    Returns:
        Safe string with only alphanumeric, dash, and underscore characters
    """
    return "".join(c if c.isalnum() or c in "-_" else "_" for c in name)


def download_contig_summary_csv(db_path, contig_name):
    """Generate CSV content with contig summary.
    
    Uses DuckDB's native COPY TO for maximum performance.
    
    Args:
        db_path: Path to DuckDB database file
        contig_name: Name of the contig
        
    Returns:
        CSV content as string, or None if error
    """
    import duckdb
    import tempfile
    import os
    
    try:
        # Open connection (read_only=True works with COPY TO since it writes to external file)
        conn = duckdb.connect(db_path, read_only=True)
        
        # Check if contig exists
        row = conn.execute("SELECT 1 FROM Contig WHERE Contig_name = ?", [contig_name]).fetchone()
        if not row:
            print(f"[downloading_data] Contig not found: {contig_name}", flush=True)
            conn.close()
            return None
        
        # DuckDB COPY TO requires the full query as a string literal — parameterized
        # queries cannot be used here. Values come from the database, not user input.
        safe_contig = contig_name.replace("'", "''")
        
        # Use DuckDB COPY TO directly with embedded query
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            temp_path = f.name.replace('\\', '/')  # Use forward slashes for DuckDB on Windows
        
        try:
            # GC_sd and GC_skew_amplitude are stored as int × 100, decode them in query
            query = f"""
                COPY (
                    SELECT
                        Contig_name as "Contig",
                        Contig_length as "Contig length",
                        ROUND(Duplication_percentage / 10.0, 1) as "Duplication (%)",
                        GC_mean as "GC mean",
                        ROUND(GC_sd / 100.0, 2) as "GC sd",
                        GC_skew_amplitude / 100.0 as "GC skew amplitude",
                        ROUND(Positive_GC_skew_windows_percentage / 10.0, 1) as "Positive GC skew windows (%)"
                    FROM Contig
                    WHERE Contig_name = '{safe_contig}'
                ) TO '{temp_path}' (HEADER, DELIMITER ',')
            """
            conn.execute(query)
            conn.close()
            
            with open(temp_path, 'r', encoding='utf-8') as f:
                csv_content = f.read()
            
            print(f"[downloading_data] Contig CSV generated", flush=True)
            return csv_content
            
        finally:
            if os.path.exists(temp_path):
                os.unlink(temp_path)
        
    except Exception as e:
        import traceback
        tb = traceback.format_exc()
        print(f"[downloading_data] Download exception: {tb}", flush=True)
        return None


def download_metrics_summary_csv(db_path, contig_name, sample_names):
    """Generate CSV content with metrics summary for all samples.
    
    Uses DuckDB's native COPY TO for maximum performance.
    
    Args:
        db_path: Path to DuckDB database file
        contig_name: Name of the contig
        sample_names: List of sample names to include
        
    Returns:
        CSV content as string, or None if error
    """
    import duckdb
    import tempfile
    import os
    
    try:
        # Open connection (read_only=True works with COPY TO since it writes to external file)
        conn = duckdb.connect(db_path, read_only=True)

        # DuckDB COPY TO requires the full query as a string literal — parameterized
        # queries cannot be used here. Values come from the database, not user input.
        safe_contig = contig_name.replace("'", "''")
        safe_samples = [s.replace("'", "''") for s in sample_names]
        samples_list = ", ".join(f"'{s}'" for s in safe_samples)

        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            temp_path = f.name.replace('\\', '/')  # Use forward slashes for DuckDB on Windows

        try:
            # Check which views exist — build query dynamically
            existing_views = set()
            for view in ['Explicit_coverage', 'Explicit_misassembly', 'Explicit_microdiversity',
                         'Explicit_side_misassembly', 'Explicit_topology', 'Explicit_phage_mechanisms']:
                try:
                    conn.execute(f"SELECT 1 FROM {view} LIMIT 0")
                    existing_views.add(view)
                except Exception:
                    pass

            select_parts = ['s.Sample_name as "Sample"']
            join_parts = []

            if 'Explicit_coverage' in existing_views:
                select_parts.extend([
                    'cov.Aligned_fraction_percentage', 'cov.Above_expected_aligned_fraction',
                    'cov.Read_count', 'cov.Coverage_mean', 'cov.Coverage_median', 'cov.Coverage_trimmed_mean',
                    'cov.RPKM', 'cov.TPM', 'cov.Coverage_coefficient_of_variation', 'cov.Relative_coverage_roughness',
                ])
                join_parts.append(f"LEFT JOIN Explicit_coverage cov ON s.Sample_name = cov.Sample_name AND cov.Contig_name = '{safe_contig}'")

            if 'Explicit_misassembly' in existing_views:
                select_parts.extend([
                    'mis.Mismatches_per_100kbp as "Misassembly_mismatches_per_100kbp"',
                    'mis.Deletions_per_100kbp as "Misassembly_deletions_per_100kbp"',
                    'mis.Insertions_per_100kbp as "Misassembly_insertions_per_100kbp"',
                    'mis.Clippings_per_100kbp as "Misassembly_clippings_per_100kbp"',
                    'mis.Collapse_bp', 'mis.Collapse_per_100kbp',
                    'mis.Expansion_bp', 'mis.Expansion_per_100kbp',
                ])
                join_parts.append(f"LEFT JOIN Explicit_misassembly mis ON s.Sample_name = mis.Sample_name AND mis.Contig_name = '{safe_contig}'")

            if 'Explicit_microdiversity' in existing_views:
                select_parts.extend([
                    'mic.Mismatches_per_100kbp as "Microdiversity_mismatches_per_100kbp"',
                    'mic.Deletions_per_100kbp as "Microdiversity_deletions_per_100kbp"',
                    'mic.Insertions_per_100kbp as "Microdiversity_insertions_per_100kbp"',
                    'mic.Clippings_per_100kbp as "Microdiversity_clippings_per_100kbp"',
                    'mic.Microdiverse_bp_on_reads', 'mic.Microdiverse_bp_per_100kbp_on_reads',
                    'mic.Microdiverse_bp_on_reference', 'mic.Microdiverse_bp_per_100kbp_on_reference',
                ])
                join_parts.append(f"LEFT JOIN Explicit_microdiversity mic ON s.Sample_name = mic.Sample_name AND mic.Contig_name = '{safe_contig}'")

            if 'Explicit_side_misassembly' in existing_views:
                select_parts.extend([
                    'sm.Coverage_first_position',
                    'sm.Contig_start_collapse_prevalence', 'sm.Contig_start_collapse_bp', 'sm.Contig_start_expansion_bp',
                    'sm.Coverage_last_position',
                    'sm.Contig_end_collapse_prevalence', 'sm.Contig_end_collapse_bp', 'sm.Contig_end_expansion_bp',
                    'sm.Contig_end_misjoint_mates', 'sm.Normalized_contig_end_misjoint_mates',
                ])
                join_parts.append(f"LEFT JOIN Explicit_side_misassembly sm ON s.Sample_name = sm.Sample_name AND sm.Contig_name = '{safe_contig}'")

            if 'Explicit_topology' in existing_views:
                select_parts.extend([
                    't.Circularising_reads', 't.Circularising_reads_prevalence',
                    't.Circularising_inserts', 't.Circularising_insert_size_deviation',
                    't.Normalized_circularising_inserts',
                ])
                join_parts.append(f"LEFT JOIN Explicit_topology t ON s.Sample_name = t.Sample_name AND t.Contig_name = '{safe_contig}'")

            if 'Explicit_phage_mechanisms' in existing_views:
                select_parts.extend([
                    'ph.Packaging_mechanism', 'ph.Left_termini', 'ph.Right_termini',
                ])
                join_parts.append(f"LEFT JOIN Explicit_phage_mechanisms ph ON s.Sample_name = ph.Sample_name AND ph.Contig_name = '{safe_contig}'")

            select_clause = ",\n                        ".join(select_parts)
            join_clause = "\n                    ".join(join_parts)

            query = f"""
                COPY (
                    WITH samples AS (
                        SELECT Sample_name FROM Sample WHERE Sample_name IN ({samples_list})
                    )
                    SELECT
                        {select_clause}
                    FROM samples s
                    {join_clause}
                    ORDER BY s.Sample_name
                ) TO '{temp_path}' (HEADER, DELIMITER ',')
            """
            conn.execute(query)
            conn.close()
            
            with open(temp_path, 'r', encoding='utf-8') as f:
                csv_content = f.read()
            
            print(f"[downloading_data] Metrics CSV generated ({len(sample_names)} samples)", flush=True)
            return csv_content
            
        finally:
            if os.path.exists(temp_path):
                os.unlink(temp_path)
        
    except Exception as e:
        import traceback
        tb = traceback.format_exc()
        print(f"[downloading_data] Download metrics exception: {tb}", flush=True)
        return None


def download_feature_data_csv(db_path, contig_name, sample_names, xstart=None, xend=None, is_all_samples=False):
    """Generate CSV content with feature data from DuckDB.

    Queries feature blobs and legacy per-position tables for the given
    contig/samples within the specified position range and returns the data as CSV.

    Args:
        db_path: Path to DuckDB database file
        contig_name: Name of the contig
        sample_names: List of sample names to include
        xstart: Optional start position for filtering (inclusive)
        xend: Optional end position for filtering (inclusive)
        is_all_samples: If True, include sample_name column for each row

    Returns:
        CSV content as string, or None if error
    """
    import duckdb
    import tempfile
    import os

    try:
        conn = duckdb.connect(db_path, read_only=True)
        cur = conn.cursor()

        # Resolve contig_id
        cur.execute("SELECT Contig_id FROM Contig WHERE Contig_name = ?", [contig_name])
        row = cur.fetchone()
        if not row:
            print(f"[downloading_data] Download data: Contig '{contig_name}' not found", flush=True)
            conn.close()
            return None
        contig_id = row[0]

        # Resolve sample_ids (may be empty in no-bam mode)
        sample_id_to_name = {}
        sample_ids = []
        sample_ids_sql = ""
        if sample_names:
            placeholders = ",".join(["?"] * len(sample_names))
            cur.execute(f"SELECT Sample_id, Sample_name FROM Sample WHERE Sample_name IN ({placeholders})", sample_names)
            sample_rows = cur.fetchall()
            sample_id_to_name = {r[0]: r[1] for r in sample_rows}
            sample_ids = list(sample_id_to_name.keys())
            sample_ids_sql = ", ".join(str(sid) for sid in sample_ids)

        # Get all feature tables and their metadata from Variable table
        cur.execute("SELECT DISTINCT Feature_table_name, Variable_name, Subplot FROM Variable WHERE Feature_table_name IS NOT NULL")
        variable_rows = cur.fetchall()
        if not variable_rows:
            print(f"[downloading_data] Download data: No feature tables found", flush=True)
            conn.close()
            return None

        # Position filter SQL fragment
        _lp = "COALESCE(Last_position, First_position)"
        pos_filter = ""
        if xstart is not None and xend is not None:
            pos_filter = f" AND {_lp} >= {int(xstart)} AND First_position <= {int(xend)}"

        # Scaled features (stored as INTEGER ×100)
        scaled_tables = {"Contig_GCSkew"}

        # Tables with a different column schema (Position1/Position2 instead of First_position/Last_position)
        non_standard_schema_tables = {"Contig_directRepeats", "Contig_invertedRepeats"}

        # Build UNION ALL query for any remaining Contig_* per-position tables
        # (Feature_blob and Contig_blob are handled below via Python decoding)
        union_parts = []
        for table_name, var_name, subplot_name in variable_rows:
            if table_name in ("Feature_blob", "Contig_blob"):
                continue
            if not table_name.startswith("Contig_"):
                continue
            if table_name in non_standard_schema_tables:
                continue

            scale_expr = " / 100.0" if table_name in scaled_tables else ""
            safe_var_name = var_name.replace("'", "''")

            if is_all_samples:
                for sid in sample_ids:
                    safe_sample = sample_id_to_name[sid].replace("'", "''")
                    union_parts.append(
                        f"SELECT '{safe_sample}' as sample_name, '{safe_var_name}' as feature_name, "
                        f"f.First_position as start_position, COALESCE(f.Last_position, f.First_position) as last_position, "
                        f"f.Value{scale_expr} as value "
                        f"FROM {table_name} f "
                        f"WHERE f.Contig_id = {contig_id}{pos_filter}"
                    )
            else:
                union_parts.append(
                    f"SELECT '{safe_var_name}' as feature_name, "
                    f"f.First_position as start_position, COALESCE(f.Last_position, f.First_position) as last_position, "
                    f"f.Value{scale_expr} as value "
                    f"FROM {table_name} f "
                    f"WHERE f.Contig_id = {contig_id}{pos_filter}"
                )

        # === BLOB export: decode Feature_blob_chunk entries and add as CSV rows ===
        blob_csv_rows = []
        blob_features = [(tn, vn, sn) for tn, vn, sn in variable_rows if tn == "Feature_blob"]
        if blob_features:
            from thebigbam.database.blob_decoder import (
                feature_name_to_id, decode_raw_chunks, decode_raw_sparse_chunks,
                get_scale_from_zoom_blob, is_sparse_zoom_blob,
            )
            import numpy as np
            for _, var_name_b, _ in blob_features:
                fid = feature_name_to_id(var_name_b, cur)
                if fid is None:
                    continue
                for sid in sample_ids:
                    try:
                        # Fetch zoom blob for scale/sparse detection
                        cur.execute(
                            "SELECT Zoom_data FROM Feature_blob WHERE Contig_id=? AND Sample_id=? AND Feature_id=?",
                            (contig_id, sid, fid)
                        )
                        zoom_row = cur.fetchone()
                        if zoom_row is None:
                            continue
                        zoom_blob = zoom_row[0]
                        scale_div = get_scale_from_zoom_blob(zoom_blob)
                        sparse = is_sparse_zoom_blob(zoom_blob)

                        # Fetch all chunks
                        cur.execute(
                            "SELECT Chunk_idx, Data FROM Feature_blob_chunk "
                            "WHERE Contig_id=? AND Sample_id=? AND Feature_id=? ORDER BY Chunk_idx",
                            (contig_id, sid, fid)
                        )
                        chunk_rows = [(r[0], r[1]) for r in cur.fetchall()]
                        if not chunk_rows:
                            continue

                        if sparse:
                            data = decode_raw_sparse_chunks(chunk_rows, scale_div)
                        else:
                            data = decode_raw_chunks(chunk_rows, scale_div)
                        x_arr = data["x"]
                        y_arr = data["y"]

                        # Apply position filter (x is 0-indexed)
                        if xstart is not None and xend is not None:
                            mask = (x_arr >= (xstart - 1)) & (x_arr <= (xend - 1))
                            x_arr = x_arr[mask]
                            y_arr = y_arr[mask]
                        sample_name_str = sample_id_to_name[sid]
                        for i in range(len(x_arr)):
                            pos = int(x_arr[i]) + 1  # 0-indexed → 1-indexed
                            val = float(y_arr[i])
                            if is_all_samples:
                                blob_csv_rows.append(f"{sample_name_str},{var_name_b},{pos},{pos},{val}")
                            else:
                                blob_csv_rows.append(f"{var_name_b},{pos},{pos},{val}")
                    except Exception:
                        pass

        # === Contig_blob export: decode contig-level features from chunks ===
        contig_blob_features = [(tn, vn, sn) for tn, vn, sn in variable_rows if tn == "Contig_blob"]
        if contig_blob_features:
            from thebigbam.database.blob_decoder import (
                contig_blob_name_to_id, decode_raw_chunks, get_scale_from_zoom_blob,
            )
            import numpy as np
            for _, var_name_b, _ in contig_blob_features:
                fid = contig_blob_name_to_id(var_name_b)
                if fid is None:
                    continue
                try:
                    # Fetch zoom blob for scale detection
                    cur.execute(
                        "SELECT Zoom_data FROM Contig_blob WHERE Contig_id=? AND Feature_id=?",
                        (contig_id, fid)
                    )
                    zoom_row = cur.fetchone()
                    if zoom_row is None:
                        continue
                    scale_div = get_scale_from_zoom_blob(zoom_row[0])

                    # Fetch all chunks
                    cur.execute(
                        "SELECT Chunk_idx, Data FROM Contig_blob_chunk "
                        "WHERE Contig_id=? AND Feature_id=? ORDER BY Chunk_idx",
                        (contig_id, fid)
                    )
                    chunk_rows = [(r[0], r[1]) for r in cur.fetchall()]
                    if not chunk_rows:
                        continue

                    data = decode_raw_chunks(chunk_rows, scale_div)
                    x_arr = data["x"]
                    y_arr = data["y"]

                    # For windowed data (GC content/skew), convert array indices to genomic positions
                    gc_window = 500 if var_name_b == "gc_content" else (1000 if var_name_b == "gc_skew" else 1)
                    if gc_window > 1:
                        x_arr = x_arr * gc_window

                    # Apply position filter (x is 0-indexed)
                    if xstart is not None and xend is not None:
                        mask = (x_arr >= (xstart - 1)) & (x_arr <= (xend - 1))
                        x_arr = x_arr[mask]
                        y_arr = y_arr[mask]
                    for i in range(len(x_arr)):
                        pos = int(x_arr[i]) + 1  # 0-indexed → 1-indexed
                        val = float(y_arr[i])
                        if is_all_samples:
                            # Contig_blob is sample-independent, but we output for all samples
                            for sid in sample_ids:
                                sample_name_str = sample_id_to_name[sid]
                                blob_csv_rows.append(f"{sample_name_str},{var_name_b},{pos},{pos},{val}")
                        else:
                            blob_csv_rows.append(f"{var_name_b},{pos},{pos},{val}")
                except Exception:
                    pass

        if not union_parts and not blob_csv_rows:
            print(f"[downloading_data] Download data: No exportable feature tables", flush=True)
            conn.close()
            return None

        if not union_parts:
            # Only BLOB data — build CSV manually
            if is_all_samples:
                header = "sample_name,feature_name,start_position,last_position,value"
            else:
                header = "feature_name,start_position,last_position,value"
            csv_content = header + "\n" + "\n".join(blob_csv_rows) + "\n"
            conn.close()
            return csv_content

        full_query = " UNION ALL ".join(union_parts)

        # Wrap with ORDER BY for deterministic output
        if is_all_samples:
            wrapper = f"SELECT sample_name, feature_name, start_position, last_position, value FROM ({full_query}) ORDER BY sample_name, feature_name, start_position"
        else:
            wrapper = f"SELECT feature_name, start_position, last_position, value FROM ({full_query}) ORDER BY feature_name, start_position"

        # Use DuckDB's native COPY TO for fast CSV export
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            temp_path = f.name.replace('\\', '/')

        try:
            copy_query = f"COPY ({wrapper}) TO '{temp_path}' (HEADER, DELIMITER ',')"
            conn.execute(copy_query)
            conn.close()

            with open(temp_path, 'r', encoding='utf-8') as f:
                csv_content = f.read()

            # Append BLOB-decoded rows if any
            if blob_csv_rows:
                csv_content = csv_content.rstrip('\n') + '\n' + '\n'.join(blob_csv_rows) + '\n'

            row_count = csv_content.count('\n') - 1
            print(f"[downloading_data] Data CSV generated ({row_count} rows)", flush=True)
            return csv_content

        finally:
            if os.path.exists(temp_path):
                os.unlink(temp_path)

    except Exception as e:
        import traceback
        tb = traceback.format_exc()
        print(f"[downloading_data] Download data exception: {tb}", flush=True)
        return None
