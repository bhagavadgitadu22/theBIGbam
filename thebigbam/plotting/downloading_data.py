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
                        Duplication_percentage as "Duplication (%)",
                        GC_mean as "GC mean",
                        GC_sd / 100.0 as "GC sd",
                        GC_skew_amplitude / 100.0 as "GC skew amplitude",
                        Positive_GC_skew_windows_percentage / 100.0 as "Positive GC skew windows (%)"
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
                        p.Aligned_fraction_percentage as "Aligned_fraction_percentage",
                        p.Coverage_mean, p.Coverage_median, p.Coverage_sd, p.Coverage_variation,
                        p.Coverage_mean_corrected_by_number_of_reads, p.Coverage_median_corrected_by_number_of_reads,
                        p.Coverage_mean_corrected_by_number_of_mapped_reads, p.Coverage_median_corrected_by_number_of_mapped_reads,
                        c.Completeness_percentage, c.Contamination_percentage,
                        c.Mismatch_frequency, c.Insertion_frequency, c.Deletion_frequency,
                        c.Read_based_clipping_frequency, c.Reference_based_clippings_frequency,
                        c.Left_completeness_percentage, c.Left_contamination_length, c.Left_missing_length,
                        c.Right_completeness_percentage, c.Right_contamination_length, c.Right_missing_length,
                        c.Circularising_reads, c.Circularising_reads_percentage,
                        ph.Packaging_mechanism, ph.Left_termini, ph.Right_termini
                    FROM samples s
                    LEFT JOIN Explicit_presences p ON s.Sample_name = p.Sample_name AND p.Contig_name = '{safe_contig}'
                    LEFT JOIN Explicit_completeness c ON s.Sample_name = c.Sample_name AND c.Contig_name = '{safe_contig}'
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


def download_feature_data_csv(db_path, contig_name, sample_names, feature_names, is_all_samples=False):
    """Generate CSV content with raw feature data (positions, values, stats).
    
    Uses DuckDB's native COPY TO for maximum performance.
    
    Args:
        db_path: Path to DuckDB database file
        contig_name: Name of the contig
        sample_names: List of sample names to include
        feature_names: List of feature/variable names to export
        is_all_samples: If True, first column is sample_name
                       If False, first column is feature_name
                       
    Returns:
        CSV content as string, or None if error
    """
    import duckdb
    import tempfile
    import os
    
    try:
        # Open connection (read_only=True works with COPY TO since it writes to external file)
        conn = duckdb.connect(db_path, read_only=True)
        cur = conn.cursor()
        
        # Get contig_id
        cur.execute("SELECT Contig_id FROM Contig WHERE Contig_name = ?", [contig_name])
        row = cur.fetchone()
        if not row:
            print(f"[downloading_data] Download data: Contig '{contig_name}' not found", flush=True)
            conn.close()
            return None
        contig_id = row[0]
        
        # Get sample_ids
        placeholders = ",".join(["?"] * len(sample_names))
        cur.execute(f"SELECT Sample_id, Sample_name FROM Sample WHERE Sample_name IN ({placeholders})", sample_names)
        sample_rows = cur.fetchall()
        if not sample_rows:
            print(f"[downloading_data] Download data: No valid samples found", flush=True)
            conn.close()
            return None
        sample_id_to_name = {r[0]: r[1] for r in sample_rows}
        sample_ids = list(sample_id_to_name.keys())
        
        # Get feature table names from Variable table (query by Subplot)
        # Also get Variable_name for each table to use in CSV output
        feature_to_tables = {}
        table_to_variable_name = {}  # Map table name to Variable_name for CSV output
        for feat in feature_names:
            cur.execute("SELECT DISTINCT Feature_table_name, Variable_name FROM Variable WHERE Subplot = ?", [feat])
            rows = cur.fetchall()
            tables = [r[0] for r in rows if r[0]]
            if tables:
                feature_to_tables[feat] = tables
                for r in rows:
                    if r[0]:
                        table_to_variable_name[r[0]] = r[1]  # table_name -> Variable_name
        
        if not feature_to_tables:
            print(f"[downloading_data] Download data: No feature tables found for {feature_names}", flush=True)
            conn.close()
            return None
        
        # Check which tables have stats columns (do once per unique table)
        all_tables = set()
        for tables in feature_to_tables.values():
            all_tables.update(tables)
        
        table_has_stats = {}
        for table_name in all_tables:
            cur.execute(f"PRAGMA table_info({table_name})")
            cols = {r[1] for r in cur.fetchall()}
            table_has_stats[table_name] = all(c in cols for c in ['Mean', 'Median', 'Std'])
        
        # Build sample_ids list as SQL string for embedding in query
        sample_ids_sql = ", ".join(str(sid) for sid in sample_ids)
        
        # Build UNION ALL query with all values embedded (no parameters)
        union_parts = []
        
        # Determine if any table has stats
        any_has_stats = any(table_has_stats.values())
        
        for feat, tables in feature_to_tables.items():
            needs_scaling = feat.lower() in ('tau', 'mapq')
            scale_factor = "/ 100.0" if needs_scaling else ""
            
            for table_name in tables:
                has_stats = table_has_stats[table_name]
                
                if is_all_samples:
                    # First column is sample_name
                    name_col = "s.Sample_name"
                else:
                    # First column is feature_name (actual Variable_name, not Subplot)
                    var_name = table_to_variable_name.get(table_name, feat)
                    safe_var_name = var_name.replace("'", "''")
                    name_col = f"'{safe_var_name}'"
                
                if has_stats:
                    select = f"""
                        SELECT {name_col} as name,
                               f.First_position as start_position,
                               f.Last_position as last_position,
                               f.Value {scale_factor} as value,
                               f.Mean {scale_factor} as mean,
                               f.Median {scale_factor} as median,
                               f.Std {scale_factor} as std
                        FROM {table_name} f
                        JOIN Sample s ON f.Sample_id = s.Sample_id
                        WHERE f.Contig_id = {contig_id} AND f.Sample_id IN ({sample_ids_sql})
                    """
                else:
                    # Pad with NULLs for consistent columns
                    if any_has_stats:
                        select = f"""
                            SELECT {name_col} as name,
                                   f.First_position as start_position,
                                   f.Last_position as last_position,
                                   f.Value {scale_factor} as value,
                                   NULL as mean,
                                   NULL as median,
                                   NULL as std
                            FROM {table_name} f
                            JOIN Sample s ON f.Sample_id = s.Sample_id
                            WHERE f.Contig_id = {contig_id} AND f.Sample_id IN ({sample_ids_sql})
                        """
                    else:
                        select = f"""
                            SELECT {name_col} as name,
                                   f.First_position as start_position,
                                   f.Last_position as last_position,
                                   f.Value {scale_factor} as value
                            FROM {table_name} f
                            JOIN Sample s ON f.Sample_id = s.Sample_id
                            WHERE f.Contig_id = {contig_id} AND f.Sample_id IN ({sample_ids_sql})
                        """
                
                union_parts.append(select)
        
        # Build combined query
        full_query = " UNION ALL ".join(union_parts)
        
        # Rename first column in the query itself
        first_col_name = 'sample_name' if is_all_samples else 'feature_name'
        
        # Wrap query to rename column and optionally drop stats columns
        if any_has_stats:
            wrapper_query = f"""
                SELECT name AS {first_col_name}, start_position, last_position, value, mean, median, std
                FROM ({full_query})
            """
        else:
            wrapper_query = f"""
                SELECT name AS {first_col_name}, start_position, last_position, value
                FROM ({full_query})
            """
        
        # Use DuckDB's native COPY TO directly (no temp table needed)
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            temp_path = f.name.replace('\\', '/')  # Use forward slashes for DuckDB on Windows
        
        try:
            # COPY directly from the query - no temp table, all done in DuckDB
            copy_query = f"COPY ({wrapper_query}) TO '{temp_path}' (HEADER, DELIMITER ',')"
            conn.execute(copy_query)
            conn.close()
            
            # Read the CSV file
            with open(temp_path, 'r', encoding='utf-8') as f:
                csv_content = f.read()
            
            # Count rows (excluding header)
            row_count = csv_content.count('\n') - 1
            print(f"[downloading_data] Data CSV generated ({row_count} rows)", flush=True)
            return csv_content
            
        finally:
            # Clean up temp file
            if os.path.exists(temp_path):
                os.unlink(temp_path)
        
    except Exception as e:
        import traceback
        tb = traceback.format_exc()
        print(f"[downloading_data] Download data exception: {tb}", flush=True)
        return None
