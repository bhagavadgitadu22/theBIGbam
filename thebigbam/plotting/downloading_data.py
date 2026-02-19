"""Functions for downloading data as CSV files.

This module contains functions for exporting data from the theBIGbam database
to CSV files. Returns CSV content as strings for browser download.
"""
import io
import csv


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
        
        # Escape single quotes in contig_name for SQL
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
        
        # Escape single quotes for SQL embedding
        safe_contig = contig_name.replace("'", "''")
        safe_samples = [s.replace("'", "''") for s in sample_names]
        samples_list = ", ".join(f"'{s}'" for s in safe_samples)
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            temp_path = f.name.replace('\\', '/')  # Use forward slashes for DuckDB on Windows
        
        try:
            query = f"""
                COPY (
                    WITH samples AS (
                        SELECT Sample_name FROM Sample WHERE Sample_name IN ({samples_list})
                    )
                    SELECT
                        s.Sample_name as "Sample",
                        cov.Aligned_fraction_percentage, cov.Above_expected_aligned_fraction,
                        cov.Read_count, cov.Coverage_mean, cov.Coverage_median, cov.Coverage_trimmed_mean,
                        cov.RPKM, cov.TPM, cov.Coverage_sd, cov.Coverage_variation,
                        mis.Mismatches_per_100kbp as "Misassembly_mismatches_per_100kbp",
                        mis.Deletions_per_100kbp as "Misassembly_deletions_per_100kbp",
                        mis.Insertions_per_100kbp as "Misassembly_insertions_per_100kbp",
                        mis.Clippings_per_100kbp as "Misassembly_clippings_per_100kbp",
                        mis.Collapse_bp, mis.Collapse_per_100kbp,
                        mis.Expansion_bp, mis.Expansion_per_100kbp,
                        mic.Mismatches_per_100kbp as "Microdiversity_mismatches_per_100kbp",
                        mic.Deletions_per_100kbp as "Microdiversity_deletions_per_100kbp",
                        mic.Insertions_per_100kbp as "Microdiversity_insertions_per_100kbp",
                        mic.Clippings_per_100kbp as "Microdiversity_clippings_per_100kbp",
                        mic.Microdiverse_bp_on_reads, mic.Microdiverse_bp_per_100kbp_on_reads,
                        mic.Microdiverse_bp_on_reference, mic.Microdiverse_bp_per_100kbp_on_reference,
                        sm.Coverage_first_position,
                        sm.Contig_start_collapse_prevalence, sm.Contig_start_collapse_bp, sm.Contig_start_expansion_bp,
                        sm.Coverage_last_position,
                        sm.Contig_end_collapse_prevalence, sm.Contig_end_collapse_bp, sm.Contig_end_expansion_bp,
                        sm.Contig_end_misjoint_mates, sm.Normalized_contig_end_misjoint_mates,
                        t.Circularising_reads, t.Circularising_reads_prevalence,
                        t.Circularising_inserts, t.Circularising_insert_size_deviation,
                        t.Normalized_circularising_inserts,
                        ph.Packaging_mechanism, ph.Left_termini, ph.Right_termini
                    FROM samples s
                    LEFT JOIN Explicit_coverage cov ON s.Sample_name = cov.Sample_name AND cov.Contig_name = '{safe_contig}'
                    LEFT JOIN Explicit_misassembly mis ON s.Sample_name = mis.Sample_name AND mis.Contig_name = '{safe_contig}'
                    LEFT JOIN Explicit_microdiversity mic ON s.Sample_name = mic.Sample_name AND mic.Contig_name = '{safe_contig}'
                    LEFT JOIN Explicit_side_misassembly sm ON s.Sample_name = sm.Sample_name AND sm.Contig_name = '{safe_contig}'
                    LEFT JOIN Explicit_topology t ON s.Sample_name = t.Sample_name AND t.Contig_name = '{safe_contig}'
                    LEFT JOIN Explicit_phage_mechanisms ph ON s.Sample_name = ph.Sample_name AND ph.Contig_name = '{safe_contig}'
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
    """Generate CSV content with raw RLE feature data from DuckDB.

    Queries all feature tables for the given contig/samples within the
    specified position range and returns the data as CSV.

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

        # Resolve sample_ids
        placeholders = ",".join(["?"] * len(sample_names))
        cur.execute(f"SELECT Sample_id, Sample_name FROM Sample WHERE Sample_name IN ({placeholders})", sample_names)
        sample_rows = cur.fetchall()
        if not sample_rows:
            print(f"[downloading_data] Download data: No valid samples found", flush=True)
            conn.close()
            return None
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

        # Check which tables have stats columns
        all_tables = set(r[0] for r in variable_rows)
        table_has_stats = {}
        for table_name in all_tables:
            try:
                cur.execute(f"PRAGMA table_info({table_name})")
                cols = {r[1] for r in cur.fetchall()}
                table_has_stats[table_name] = all(c in cols for c in ['Mean', 'Median', 'Std'])
            except Exception:
                table_has_stats[table_name] = False

        any_has_stats = any(table_has_stats.values())

        # Position filter SQL fragment
        # Use COALESCE to handle NULL Last_position (single-position bar features)
        _lp = "COALESCE(Last_position, First_position)"
        pos_filter = ""
        if xstart is not None and xend is not None:
            pos_filter = f" AND {_lp} >= {int(xstart)} AND First_position <= {int(xend)}"

        # Scaled features (stored as INTEGER ×100)
        scaled_tables = {"Feature_mapq", "Contig_GCSkew"}

        # Relative-scaled features (stored as INTEGER ×1000)
        RELATIVE_SCALED_FEATURES = {
            "Feature_left_clippings", "Feature_right_clippings", "Feature_insertions", "Feature_deletions",
            "Feature_mismatches", "Feature_reads_starts", "Feature_reads_ends", "Feature_non_inward_pairs",
            "Feature_mate_not_mapped", "Feature_mate_on_another_contig"
        }

        # Tables with a different column schema (Position1/Position2 instead of First_position/Last_position)
        non_standard_schema_tables = {"Contig_directRepeats", "Contig_invertedRepeats"}

        # Build UNION ALL query across all feature tables
        union_parts = []
        for table_name, var_name, subplot_name in variable_rows:
            is_contig_table = table_name.startswith("Contig_")

            # Skip tables with non-standard column schema
            if table_name == "Feature_primary_reads":
                continue
            if table_name in non_standard_schema_tables:
                continue

            # Determine scaling factor
            if table_name in RELATIVE_SCALED_FEATURES:
                scale_expr = " / 1000.0"
            elif table_name in scaled_tables:
                scale_expr = " / 100.0"
            else:
                scale_expr = ""

            safe_var_name = var_name.replace("'", "''")
            has_stats = table_has_stats.get(table_name, False)

            if is_contig_table:
                # Contig-level tables have no Sample_id — emit one row per sample with same data
                if is_all_samples:
                    for sid in sample_ids:
                        safe_sample = sample_id_to_name[sid].replace("'", "''")
                        stats_cols = f", f.Mean{scale_expr} as mean, f.Median{scale_expr} as median, f.Std{scale_expr} as std" if has_stats else (", NULL as mean, NULL as median, NULL as std" if any_has_stats else "")
                        union_parts.append(
                            f"SELECT '{safe_sample}' as sample_name, '{safe_var_name}' as feature_name, "
                            f"f.First_position as start_position, COALESCE(f.Last_position, f.First_position) as last_position, "
                            f"f.Value{scale_expr} as value{stats_cols} "
                            f"FROM {table_name} f "
                            f"WHERE f.Contig_id = {contig_id}{pos_filter}"
                        )
                else:
                    stats_cols = f", f.Mean{scale_expr} as mean, f.Median{scale_expr} as median, f.Std{scale_expr} as std" if has_stats else (", NULL as mean, NULL as median, NULL as std" if any_has_stats else "")
                    union_parts.append(
                        f"SELECT '{safe_var_name}' as feature_name, "
                        f"f.First_position as start_position, COALESCE(f.Last_position, f.First_position) as last_position, "
                        f"f.Value{scale_expr} as value{stats_cols} "
                        f"FROM {table_name} f "
                        f"WHERE f.Contig_id = {contig_id}{pos_filter}"
                    )
            else:
                # Sample-level tables
                stats_cols = f", f.Mean{scale_expr} as mean, f.Median{scale_expr} as median, f.Std{scale_expr} as std" if has_stats else (", NULL as mean, NULL as median, NULL as std" if any_has_stats else "")
                if is_all_samples:
                    union_parts.append(
                        f"SELECT s.Sample_name as sample_name, '{safe_var_name}' as feature_name, "
                        f"f.First_position as start_position, COALESCE(f.Last_position, f.First_position) as last_position, "
                        f"f.Value{scale_expr} as value{stats_cols} "
                        f"FROM {table_name} f "
                        f"JOIN Sample s ON f.Sample_id = s.Sample_id "
                        f"WHERE f.Contig_id = {contig_id} AND f.Sample_id IN ({sample_ids_sql}){pos_filter}"
                    )
                else:
                    union_parts.append(
                        f"SELECT '{safe_var_name}' as feature_name, "
                        f"f.First_position as start_position, COALESCE(f.Last_position, f.First_position) as last_position, "
                        f"f.Value{scale_expr} as value{stats_cols} "
                        f"FROM {table_name} f "
                        f"WHERE f.Contig_id = {contig_id} AND f.Sample_id IN ({sample_ids_sql}){pos_filter}"
                    )

        if not union_parts:
            print(f"[downloading_data] Download data: No exportable feature tables", flush=True)
            conn.close()
            return None

        full_query = " UNION ALL ".join(union_parts)

        # Wrap with ORDER BY for deterministic output
        if is_all_samples:
            if any_has_stats:
                wrapper = f"SELECT sample_name, feature_name, start_position, last_position, value, mean, median, std FROM ({full_query}) ORDER BY sample_name, feature_name, start_position"
            else:
                wrapper = f"SELECT sample_name, feature_name, start_position, last_position, value FROM ({full_query}) ORDER BY sample_name, feature_name, start_position"
        else:
            if any_has_stats:
                wrapper = f"SELECT feature_name, start_position, last_position, value, mean, median, std FROM ({full_query}) ORDER BY feature_name, start_position"
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
