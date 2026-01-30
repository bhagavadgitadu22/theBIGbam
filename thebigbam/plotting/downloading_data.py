"""Functions for downloading data as CSV files.

This module contains functions for exporting data from the theBIGbam database
to CSV files via Save As dialogs.
"""


def _save_csv_with_dialog(csv_content, default_filename):
    """Show Save As dialog and save CSV content to user-selected location.
    
    Args:
        csv_content: The CSV content as a string
        default_filename: Suggested filename for the save dialog
        
    Returns:
        True if file was saved, False if cancelled or error
    """
    import os
    try:
        import tkinter as tk
        from tkinter import filedialog
        
        # Create hidden root window
        root = tk.Tk()
        root.withdraw()
        root.attributes('-topmost', True)  # Bring dialog to front
        
        # Get user's home directory as initial location
        initial_dir = os.path.expanduser("~")
        # Try Documents folder if it exists
        docs_dir = os.path.join(initial_dir, "Documents")
        if os.path.isdir(docs_dir):
            initial_dir = docs_dir
        
        # Show Save As dialog
        file_path = filedialog.asksaveasfilename(
            parent=root,
            initialdir=initial_dir,
            initialfile=default_filename,
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
            title="Save CSV file"
        )
        
        root.destroy()
        
        if not file_path:
            print(f"[downloading_data] Save cancelled by user", flush=True)
            return False
        
        # Write the file
        with open(file_path, 'w', encoding='utf-8') as f:
            f.write(csv_content)
        
        print(f"[downloading_data] CSV saved to {file_path}", flush=True)
        return True
        
    except Exception as e:
        print(f"[downloading_data] Error saving file: {e}", flush=True)
        # Fallback: save to temp directory
        import tempfile
        temp_dir = tempfile.gettempdir()
        temp_path = os.path.join(temp_dir, default_filename)
        try:
            with open(temp_path, 'w', encoding='utf-8') as f:
                f.write(csv_content)
            print(f"[downloading_data] Fallback: CSV saved to {temp_path}", flush=True)
            return True
        except Exception as e2:
            print(f"[downloading_data] Fallback save also failed: {e2}", flush=True)
            return False


def make_safe_filename(name):
    """Sanitize a string for use in filenames.
    
    Args:
        name: String to sanitize
        
    Returns:
        Safe string with only alphanumeric, dash, and underscore characters
    """
    return "".join(c if c.isalnum() or c in "-_" else "_" for c in name)


def download_contig_summary_csv(conn, contig_name):
    """Generate CSV with contig summary and download it.
    
    Args:
        conn: DuckDB connection
        contig_name: Name of the contig
        
    Returns:
        True if successful, False otherwise
    """
    try:
        # Get contig info
        cur = conn.cursor()
        cur.execute("SELECT Contig_length, Duplication_percentage, GC_mean, GC_sd, GC_median FROM Contig WHERE Contig_name = ?", (contig_name,))
        result = cur.fetchone()
        
        if not result:
            print(f"[downloading_data] Contig not found: {contig_name}", flush=True)
            return False
        
        contig_length = result[0] if result[0] is not None else ""
        contig_duplication = result[1] if result[1] is not None else ""
        gc_mean = result[2] if result[2] is not None else ""
        gc_sd = result[3] if result[3] is not None else ""
        gc_median = result[4] if result[4] is not None else ""
        
        # Generate CSV content
        csv_lines = [
            "Contig,Contig length,Duplication (%),GC mean,GC sd,GC median",
            f"{contig_name},{contig_length},{contig_duplication},{gc_mean},{gc_sd},{gc_median}"
        ]
        csv_content = "\n".join(csv_lines)
        
        # Create filename from contig name (sanitize for filesystem)
        filename = f"{make_safe_filename(contig_name)}_summary.csv"
        
        # Show Save As dialog
        return _save_csv_with_dialog(csv_content, filename)
        
    except Exception as e:
        import traceback
        tb = traceback.format_exc()
        print(f"[downloading_data] Download exception: {tb}", flush=True)
        return False


def download_metrics_summary_csv(conn, contig_name, sample_names, filename):
    """Generate CSV with metrics summary for all samples and download it.
    
    Args:
        conn: DuckDB connection
        contig_name: Name of the contig
        sample_names: List of sample names to include
        filename: Name for the output CSV file
        
    Returns:
        True if successful, False otherwise
    """
    try:
        # Define all metric columns in order
        coverage_cols = [
            "Aligned_fraction_percentage", "Coverage_mean", "Coverage_median",
            "Coverage_sd", "Coverage_variation", 
            "Coverage_mean_corrected_by_number_of_reads", "Coverage_median_corrected_by_number_of_reads", 
            "Coverage_mean_corrected_by_number_of_mapped_reads", "Coverage_median_corrected_by_number_of_mapped_reads"
        ]

        completeness_cols = [
            "Completeness_percentage", "Contamination_percentage",
            "Mismatch_frequency", "Insertion_frequency", "Deletion_frequency",
            "Read_based_clipping_frequency", "Reference_based_clippings_frequency"
        ]

        side_completeness_cols = [
            "Left_completeness_percentage", "Left_contamination_length", "Left_missing_length",
            "Right_completeness_percentage", "Right_contamination_length", "Right_missing_length",
            "Circularising_reads", "Circularising_reads_percentage"
        ]

        phage_cols = [
            "Packaging_mechanism", "Left_termini", "Right_termini"
        ]
        
        all_metric_cols = coverage_cols + completeness_cols + side_completeness_cols + phage_cols
        
        cur = conn.cursor()
        
        # Initialize data dict
        n = len(sample_names)
        sample_idx = {name: i for i, name in enumerate(sample_names)}
        data = {col: [None] * n for col in all_metric_cols}
        data["Sample"] = list(sample_names)
        
        # Query Explicit_presences view for coverage data
        try:
            cols_str = ", ".join(coverage_cols)
            query = f"""
                SELECT Sample_name, {cols_str}
                FROM Explicit_presences
                WHERE Contig_name = ? AND Sample_name IN ({','.join(['?'] * len(sample_names))})
            """
            cur.execute(query, [contig_name] + list(sample_names))
            for row in cur.fetchall():
                sample_name, *values = row
                idx = sample_idx.get(sample_name)
                if idx is None:
                    continue
                for col, val in zip(coverage_cols, values):
                    data[col][idx] = val
        except Exception:
            pass
        
        # Query Explicit_completeness view
        try:
            comp_cols = completeness_cols + side_completeness_cols
            cols_str = ", ".join(comp_cols)
            query = f"""
                SELECT Sample_name, {cols_str}
                FROM Explicit_completeness
                WHERE Contig_name = ? AND Sample_name IN ({','.join(['?'] * len(sample_names))})
            """
            cur.execute(query, [contig_name] + list(sample_names))
            for row in cur.fetchall():
                sample_name, *values = row
                idx = sample_idx.get(sample_name)
                if idx is None:
                    continue
                for col, val in zip(comp_cols, values):
                    data[col][idx] = val
        except Exception:
            pass
        
        # Query Explicit_phage_mechanisms view
        try:
            cols_str = ", ".join(phage_cols)
            query = f"""
                SELECT Sample_name, {cols_str}
                FROM Explicit_phage_mechanisms
                WHERE Contig_name = ? AND Sample_name IN ({','.join(['?'] * len(sample_names))})
            """
            cur.execute(query, [contig_name] + list(sample_names))
            for row in cur.fetchall():
                sample_name, *values = row
                idx = sample_idx.get(sample_name)
                if idx is None:
                    continue
                for col, val in zip(phage_cols, values):
                    data[col][idx] = val
        except Exception:
            pass
        
        # Build CSV content
        # Header row: Sample + all metric columns
        header = ["Sample"] + all_metric_cols
        csv_lines = [",".join(header)]
        
        # Data rows
        for i in range(n):
            row_values = [data["Sample"][i]]
            for col in all_metric_cols:
                val = data[col][i]
                if val is None:
                    row_values.append("")
                elif isinstance(val, str):
                    # Escape commas and quotes in string values
                    escaped = val.replace('"', '""')
                    if ',' in escaped or '"' in escaped:
                        row_values.append(f'"{escaped}"')
                    else:
                        row_values.append(escaped)
                else:
                    row_values.append(str(val))
            csv_lines.append(",".join(row_values))
        
        csv_content = "\n".join(csv_lines)
        
        # Show Save As dialog
        return _save_csv_with_dialog(csv_content, filename)
        
    except Exception as e:
        import traceback
        tb = traceback.format_exc()
        print(f"[downloading_data] Download metrics exception: {tb}", flush=True)
        return False


def download_feature_data_csv(conn, contig_name, sample_names, feature_names, filename, is_all_samples=False):
    """Generate CSV with raw feature data (positions, values, stats) and download it.
    
    Args:
        conn: DuckDB connection
        contig_name: Name of the contig
        sample_names: List of sample names to include
        feature_names: List of feature/variable names to export
        filename: Output filename for CSV
        is_all_samples: If True, columns are sample_name,start,end,value,...
                       If False, columns are feature_name,start,end,value,...
    """
    try:
        cur = conn.cursor()
        
        # Get contig_id
        cur.execute("SELECT Contig_id FROM Contig WHERE Contig_name = ?", [contig_name])
        row = cur.fetchone()
        if not row:
            print(f"[downloading_data] Download data: Contig '{contig_name}' not found", flush=True)
            return False
        contig_id = row[0]
        
        # Build sample_id mapping
        placeholders = ",".join(["?"] * len(sample_names))
        cur.execute(f"SELECT Sample_id, Sample_name FROM Sample WHERE Sample_name IN ({placeholders})", sample_names)
        sample_id_to_name = {r[0]: r[1] for r in cur.fetchall()}
        sample_name_to_id = {v: k for k, v in sample_id_to_name.items()}
        sample_ids = [sample_name_to_id[s] for s in sample_names if s in sample_name_to_id]
        
        if not sample_ids:
            print(f"[downloading_data] Download data: No valid samples found", flush=True)
            return False
        
        # Get feature table names from Variable table (query by Subplot, which is what UI labels use)
        feature_to_table = {}
        for feat in feature_names:
            cur.execute("SELECT DISTINCT Feature_table_name FROM Variable WHERE Subplot = ?", [feat])
            rows = cur.fetchall()
            for row in rows:
                if row[0]:
                    # Multiple variables may share the same subplot (e.g., left/right clippings)
                    # Store all unique table names for this subplot
                    if feat not in feature_to_table:
                        feature_to_table[feat] = []
                    if row[0] not in feature_to_table[feat]:
                        feature_to_table[feat].append(row[0])
        
        if not feature_to_table:
            print(f"[downloading_data] Download data: No feature tables found for {feature_names}", flush=True)
            return False
        
        # Count total rows to warn user if large
        total_rows = 0
        for feat, table_names in feature_to_table.items():
            for table_name in table_names:
                sample_placeholders = ",".join(["?"] * len(sample_ids))
                cur.execute(
                    f"SELECT COUNT(*) FROM {table_name} WHERE Contig_id = ? AND Sample_id IN ({sample_placeholders})",
                    [contig_id] + sample_ids
                )
                total_rows += cur.fetchone()[0]
        
        # If large dataset, ask user for confirmation
        ROW_THRESHOLD = 100000
        if total_rows > ROW_THRESHOLD:
            estimated_size_mb = total_rows * 50 / (1024 * 1024)  # ~50 bytes per row estimate
            try:
                import tkinter as tk
                from tkinter import messagebox
                root = tk.Tk()
                root.withdraw()  # Hide the main window
                root.attributes('-topmost', True)  # Bring dialog to front
                result = messagebox.askyesno(
                    "Large Data Export",
                    f"The export contains {total_rows:,} rows (~{estimated_size_mb:.1f} MB).\n\n"
                    f"Do you want to continue?",
                    parent=root
                )
                root.destroy()
                if not result:
                    print(f"[downloading_data] Download data: User cancelled large export ({total_rows} rows)", flush=True)
                    return False
            except Exception as tk_err:
                # If tkinter fails (e.g., no display), proceed anyway
                print(f"[downloading_data] Warning: Could not show confirmation dialog: {tk_err}", flush=True)
        
        # Collect all data
        all_rows = []
        stats_columns_present = False
        
        for feat, table_names in feature_to_table.items():
            for table_name in table_names:
                # Check if table has stats columns (Mean, Median, Std)
                cur.execute(f"PRAGMA table_info({table_name})")
                columns_info = {r[1]: r for r in cur.fetchall()}
                has_stats = all(col in columns_info for col in ['Mean', 'Median', 'Std'])
                if has_stats:
                    stats_columns_present = True
                
                # Build query
                sample_placeholders = ",".join(["?"] * len(sample_ids))
                if has_stats:
                    query = f"SELECT Sample_id, First_position, Last_position, Value, Mean, Median, Std FROM {table_name} WHERE Contig_id = ? AND Sample_id IN ({sample_placeholders}) ORDER BY Sample_id, First_position"
                else:
                    query = f"SELECT Sample_id, First_position, Last_position, Value FROM {table_name} WHERE Contig_id = ? AND Sample_id IN ({sample_placeholders}) ORDER BY Sample_id, First_position"
                
                cur.execute(query, [contig_id] + sample_ids)
                rows = cur.fetchall()
                
                # Scale tau and mapq values (stored as int * 100)
                needs_scaling = feat.lower() in ('tau', 'mapq')
                
                for r in rows:
                    sample_id = r[0]
                    sample_name = sample_id_to_name.get(sample_id, str(sample_id))
                    start_pos = r[1]
                    end_pos = r[2]
                    value = r[3] / 100.0 if needs_scaling else r[3]
                    
                    if has_stats:
                        mean_val = r[4] / 100.0 if needs_scaling and r[4] is not None else r[4]
                        median_val = r[5] / 100.0 if needs_scaling and r[5] is not None else r[5]
                        std_val = r[6] / 100.0 if needs_scaling and r[6] is not None else r[6]
                        if is_all_samples:
                            all_rows.append((sample_name, start_pos, end_pos, value, mean_val, median_val, std_val))
                        else:
                            all_rows.append((feat, start_pos, end_pos, value, mean_val, median_val, std_val))
                    else:
                        if is_all_samples:
                            all_rows.append((sample_name, start_pos, end_pos, value, None, None, None))
                        else:
                            all_rows.append((feat, start_pos, end_pos, value, None, None, None))
        
        # Build CSV content
        if is_all_samples:
            if stats_columns_present:
                header = "sample_name,start_position,last_position,value,mean,median,std"
            else:
                header = "sample_name,start_position,last_position,value"
        else:
            if stats_columns_present:
                header = "feature_name,start_position,last_position,value,mean,median,std"
            else:
                header = "feature_name,start_position,last_position,value"
        
        csv_lines = [header]
        for r in all_rows:
            if stats_columns_present:
                # Format values, handling None for missing stats
                vals = [
                    str(r[0]),  # name
                    str(r[1]),  # start
                    str(r[2]),  # end
                    str(r[3]) if r[3] is not None else "",  # value
                    str(r[4]) if r[4] is not None else "",  # mean
                    str(r[5]) if r[5] is not None else "",  # median
                    str(r[6]) if r[6] is not None else ""   # std
                ]
            else:
                vals = [
                    str(r[0]),  # name
                    str(r[1]),  # start
                    str(r[2]),  # end
                    str(r[3]) if r[3] is not None else ""   # value
                ]
            csv_lines.append(",".join(vals))
        
        csv_content = "\n".join(csv_lines)
        
        print(f"[downloading_data] Data CSV generated ({len(all_rows)} rows)", flush=True)
        
        # Show Save As dialog
        return _save_csv_with_dialog(csv_content, filename)
        
    except Exception as e:
        import traceback
        tb = traceback.format_exc()
        print(f"[downloading_data] Download data exception: {tb}", flush=True)
        return False
