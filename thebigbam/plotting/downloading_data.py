"""Functions for downloading data as CSV files.

This module contains functions for exporting contig/metrics summary data
from the theBIGbam database to CSV files.
"""


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
