from bokeh.models import Div, ColumnDataSource, DataTable, TableColumn


## Function to generate HTML and open in new browser window
def generate_and_open_peruse_html(conn, contig_name, sample_names):
    """Generate HTML summary and open in new browser window.
    
    Args:
        conn: DuckDB connection
        contig_name: Name of the contig
        sample_names: List of sample names to include
        
    Returns:
        True if successful, False otherwise
    """
    import tempfile
    import webbrowser
    
    try:
        # Get contig info
        cur = conn.cursor()
        cur.execute("SELECT Contig_length, Duplication_percentage, GC_mean, GC_sd, GC_skew_amplitude, Positive_GC_skew_windows_percentage FROM Contig WHERE Contig_name = ?", (contig_name,))
        result = cur.fetchone()
        contig_length = result[0] if result else "unknown"
        contig_duplication = round(result[1] / 10.0, 1) if result and result[1] is not None else "unknown"
        gc_mean = result[2] if result and result[2] is not None else "N/A"
        # GC_sd and GC_skew_amplitude are stored as int × 100, decode them
        gc_sd = round(result[3] / 100.0, 2) if result and result[3] is not None else "N/A"
        gc_skew_amplitude = result[4] / 100.0 if result and result[4] is not None else "N/A"
        gc_skew_percent_positive = round(result[5] / 10.0, 1) if result and result[5] is not None else "N/A"

        # Build data
        content = build_summary_data(conn, contig_name, sample_names)

        if content is None or len(content) == 0:
            print("[perusing_data] No data available", flush=True)
            return False

        # Build contig features table HTML
        contig_table_html = f"""
        <table style="border-collapse: collapse; margin-bottom: 20px; background-color: white; border-radius: 5px;">
            <thead>
                <tr style="background-color: #f0f0f0;">
                    <th style="border: 1px solid #ddd; padding: 8px; text-align: left;">Contig</th>
                    <th style="border: 1px solid #ddd; padding: 8px; text-align: left;">Contig length</th>
                    <th style="border: 1px solid #ddd; padding: 8px; text-align: left;">Duplication (%)</th>
                    <th style="border: 1px solid #ddd; padding: 8px; text-align: left;">GC mean</th>
                    <th style="border: 1px solid #ddd; padding: 8px; text-align: left;">GC sd</th>
                    <th style="border: 1px solid #ddd; padding: 8px; text-align: left;">GC skew amplitude</th>
                    <th style="border: 1px solid #ddd; padding: 8px; text-align: left;">% positive GC skew</th>
                </tr>
            </thead>
            <tbody>
                <tr>
                    <td style="border: 1px solid #ddd; padding: 8px;">{contig_name}</td>
                    <td style="border: 1px solid #ddd; padding: 8px;">{contig_length}</td>
                    <td style="border: 1px solid #ddd; padding: 8px;">{contig_duplication}</td>
                    <td style="border: 1px solid #ddd; padding: 8px;">{gc_mean}</td>
                    <td style="border: 1px solid #ddd; padding: 8px;">{gc_sd}</td>
                    <td style="border: 1px solid #ddd; padding: 8px;">{gc_skew_amplitude}</td>
                    <td style="border: 1px solid #ddd; padding: 8px;">{gc_skew_percent_positive}</td>
                </tr>
            </tbody>
        </table>
        """

        # Generate HTML for new window
        html_parts = [
            "<!DOCTYPE html>",
            "<html>",
            "<head>",
            "<meta charset='utf-8'>",
            f"<title>theBIGbam - {contig_name} Summary</title>",
            "<style>",
            "body { font-family: Arial, sans-serif; margin: 20px; background-color: #f5f5f5; }",
            "h2, h3 { color: #333; }",
            "b { display: block; margin-top: 15px; margin-bottom: 5px; }",
            "i { display: block; margin-bottom: 15px; font-size: 0.9em; color: #666; }",
            "</style>",
            "</head>",
            "<body>",
            f"<h2>{contig_name}</h2>",
            contig_table_html,
            "<h3>Metrics summary:</h3>",
        ]
        
        # Add each content item (all are HTML strings now)
        for item in content:
            html_parts.append(item)
        
        html_parts.extend([
            "</body>",
            "</html>"
        ])
        
        html_content = "\n".join(html_parts)
        
        # Write to temporary file and open in new window
        with tempfile.NamedTemporaryFile(mode='w', suffix='.html', delete=False, encoding='utf-8') as f:
            f.write(html_content)
            temp_path = f.name
        
        # Open in new browser window/tab
        from pathlib import Path
        import subprocess
        import sys
        import os
        
        opened = False
        
        # Check if running in WSL (Windows Subsystem for Linux)
        is_wsl = 'microsoft' in os.uname().release.lower() if hasattr(os, 'uname') else False
        
        if sys.platform == 'win32':
            # Windows: use os.startfile
            os.startfile(temp_path)
            opened = True
        elif is_wsl:
            # WSL: convert path and use Windows browser via cmd.exe
            # Convert /mnt/c/... to C:\...
            if temp_path.startswith('/mnt/'):
                win_path = temp_path.replace('/mnt/', '').replace('/', '\\')
                win_path = win_path[0].upper() + ':' + win_path[1:]
            else:
                # For /tmp paths, use wslpath or keep as-is for wslview
                try:
                    result = subprocess.run(['wslpath', '-w', temp_path], capture_output=True, text=True)
                    win_path = result.stdout.strip() if result.returncode == 0 else temp_path
                except Exception:
                    win_path = temp_path
            
            # Try wslview first (from wslu package), then cmd.exe
            try:
                subprocess.Popen(['wslview', temp_path], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                opened = True
            except (FileNotFoundError, PermissionError):
                try:
                    subprocess.Popen(['cmd.exe', '/c', 'start', '', win_path], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                    opened = True
                except Exception:
                    pass
        elif sys.platform == 'darwin':
            # macOS: use 'open' command
            subprocess.Popen(['open', temp_path])
            opened = True
        else:
            # Native Linux: try browsers directly
            file_url = Path(temp_path).as_uri()
            for cmd in ['firefox', 'google-chrome', 'chromium', 'chromium-browser', 'xdg-open']:
                try:
                    subprocess.Popen([cmd, file_url], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                    opened = True
                    break
                except (FileNotFoundError, PermissionError):
                    continue
        
        if not opened:
            # Fallback to webbrowser module
            file_url = Path(temp_path).as_uri()
            webbrowser.open(file_url, new=2)
        
        print(f"[perusing_data] Opened {temp_path} in browser", flush=True)
        return True
        
    except Exception as e:
        import traceback
        tb = traceback.format_exc()
        print(f"[perusing_data] Exception: {tb}", flush=True)
        return False


## Helper function to build summary data for peruse button
def build_summary_data(conn, contig_name, sample_names):
    """Query database views and return dict suitable for ColumnDataSource.

    Args:
        contig_name: Name of the contig to query
        sample_names: List of sample names to include

    Returns:
        Dict with column names as keys and lists of values
    """
    # Define column groups for each subsection
    coverage_cols = {
        "Aligned_fraction_percentage": "Aligned fraction (%)",
        "Above_expected_aligned_fraction": "Above expected AF",
        "Read_count": "Read count",
        "Coverage_mean": "Coverage mean",
        "Coverage_median": "Coverage median",
        "Coverage_trimmed_mean": "Coverage trimmed mean",
        "RPKM": "RPKM",
        "TPM": "TPM",
        "Coverage_sd": "Coverage sd",
        "Coverage_variation": "Coverage variation",
    }

    misassembly_cols = {
        "Mismatches_per_100kbp": "Mismatches (per 100kbp)",
        "Deletions_per_100kbp": "Deletions (per 100kbp)",
        "Insertions_per_100kbp": "Insertions (per 100kbp)",
        "Clippings_per_100kbp": "Clippings (per 100kbp)",
        "Collapse_bp": "Collapse (bp)",
        "Collapse_per_100kbp": "Collapse (per 100kbp)",
        "Expansion_bp": "Expansion (bp)",
        "Expansion_per_100kbp": "Expansion (per 100kbp)",
    }

    microdiversity_cols = {
        "Micro_Mismatches_per_100kbp": "Mismatches (per 100kbp)",
        "Micro_Deletions_per_100kbp": "Deletions (per 100kbp)",
        "Micro_Insertions_per_100kbp": "Insertions (per 100kbp)",
        "Micro_Clippings_per_100kbp": "Clippings (per 100kbp)",
        "Microdiverse_bp_on_reads": "Microdiverse bp (reads)",
        "Microdiverse_bp_per_100kbp_on_reads": "Microdiverse bp/100kbp (reads)",
        "Microdiverse_bp_on_reference": "Microdiverse bp (reference)",
        "Microdiverse_bp_per_100kbp_on_reference": "Microdiverse bp/100kbp (reference)",
    }

    side_misassembly_cols = {
        "Coverage_first_position": "Coverage first position",
        "Contig_start_collapse_prevalence": "Start collapse prevalence",
        "Contig_start_collapse_bp": "Start collapse (bp)",
        "Contig_start_expansion_bp": "Start expansion (bp)",
        "Coverage_last_position": "Coverage last position",
        "Contig_end_collapse_prevalence": "End collapse prevalence",
        "Contig_end_collapse_bp": "End collapse (bp)",
        "Contig_end_expansion_bp": "End expansion (bp)",
        "Contig_end_misjoint_mates": "End misjoint mates",
        "Normalized_contig_end_misjoint_mates": "Normalized end misjoint mates",
    }

    phage_cols = {
        "Packaging_mechanism": "Mechanism", "Left_termini": "Left termini", "Right_termini": "Right termini"
    }

    topology_cols = {
        "Circularising_reads": "Circularising reads",
        "Circularising_reads_prevalence": "Circularising reads prevalence",
        "Circularising_inserts": "Circularising inserts",
        "Circularising_insert_size_deviation": "Circularising insert size deviation",
        "Normalized_circularising_inserts": "Normalized circularising inserts",
    }

    cur = conn.cursor()

    # Initialize data dict with Sample column
    data = {"Sample": list(sample_names)}
    n = len(sample_names)
    # Create sample name to index mapping
    sample_idx = {name: i for i, name in enumerate(sample_names)}

    if not sample_names:
        return  # nothing to query

    # Check if any sample is paired-read
    placeholders = ','.join(['?'] * len(sample_names))
    cur.execute(
        f"SELECT COUNT(*) FROM Sample WHERE Sample_name IN ({placeholders}) AND Sequencing_type = 'paired-short'",
        list(sample_names)
    )
    has_paired = cur.fetchone()[0] > 0

    if not has_paired:
        for key in ("Contig_end_misjoint_mates", "Normalized_contig_end_misjoint_mates"):
            side_misassembly_cols.pop(key, None)
        for key in ("Circularising_inserts", "Circularising_insert_size_deviation", "Normalized_circularising_inserts"):
            topology_cols.pop(key, None)

    # Initialize all columns with None
    all_cols = coverage_cols | misassembly_cols | microdiversity_cols | side_misassembly_cols | phage_cols | topology_cols
    for col in all_cols.keys():
        data[col] = [None] * n

    # SQL column name mapping for microdiversity (data-dict keys differ from SQL columns)
    microdiversity_sql_map = {
        "Micro_Mismatches_per_100kbp": "Mismatches_per_100kbp",
        "Micro_Deletions_per_100kbp": "Deletions_per_100kbp",
        "Micro_Insertions_per_100kbp": "Insertions_per_100kbp",
        "Micro_Clippings_per_100kbp": "Clippings_per_100kbp",
    }

    # Helper to query a view and fill data
    def query_view(view_name, col_dict, sql_col_map=None):
        try:
            if sql_col_map:
                cols_str = ", ".join(sql_col_map.get(k, k) for k in col_dict.keys())
            else:
                cols_str = ", ".join(col_dict.keys())
            query = f"""
                SELECT Sample_name, {cols_str}
                FROM {view_name}
                WHERE Contig_name = ? AND Sample_name IN ({','.join(['?'] * len(sample_names))})
            """
            cur.execute(query, [contig_name] + list(sample_names))
            for row in cur.fetchall():
                sample_name, *values = row
                idx = sample_idx.get(sample_name)
                if idx is None:
                    continue
                for value_col, cell in zip(col_dict.keys(), values):
                    data[value_col][idx] = cell
        except Exception:
            pass  # View might not exist or have no data

    query_view("Explicit_coverage", coverage_cols)
    query_view("Explicit_misassembly", misassembly_cols)
    query_view("Explicit_microdiversity", microdiversity_cols, sql_col_map=microdiversity_sql_map)
    query_view("Explicit_side_misassembly", side_misassembly_cols)
    query_view("Explicit_phage_mechanisms", phage_cols)
    query_view("Explicit_topology", topology_cols)

    # Create content (as HTML strings)
    content = []

    # Coverage subsection
    coverage_table = generate_summary_table_html(data, coverage_cols)
    if coverage_table:
        content.append("<b>Coverage:</b>")
        content.append(coverage_table)

    # Misassembly subsection
    misassembly_table = generate_summary_table_html(data, misassembly_cols)
    if misassembly_table:
        content.append("<b>Misassembly:</b>")
        content.append(misassembly_table)
        content.append(
            "<i>Positions with ≥50% prevalence. "
            "Collapse: extra bases in reads (insertions + clippings + mismatches). "
            "Expansion: missing bases in reads (deletions + mismatches + paired clip distances).</i>"
        )

    # Microdiversity subsection
    microdiversity_table = generate_summary_table_html(data, microdiversity_cols)
    if microdiversity_table:
        content.append("<b>Microdiversity:</b>")
        content.append(microdiversity_table)
        content.append(
            "<i>Positions with ≥10% prevalence.</i>"
        )

    # Side misassembly subsection
    side_misassembly_table = generate_summary_table_html(data, side_misassembly_cols)
    if side_misassembly_table:
        content.append("<b>Side misassembly:</b>")
        content.append(side_misassembly_table)
        content.append(
            "<i>Clipping events at contig extremities with ≥50% prevalence. "
            "Collapse: clipped bases extending beyond contig boundary. "
            "Expansion: distance from clipping site to contig end.</i>"
        )

    # Topology subsection
    topology_table = generate_summary_table_html(data, topology_cols)
    if topology_table:
        content.append("<b>Topology:</b>")
        content.append(topology_table)

    # Phage mechanism subsection
    phage_table = generate_summary_table_html(data, phage_cols)
    if phage_table:
        content.append("<b>Phage mechanism:</b>")
        content.append(phage_table)

    return content

## Helper function to round to N significant figures
def round_to_n_sigfigs(value, n=2):
    """Round a number to n significant figures."""
    if value is None or not isinstance(value, (float)):
        return value
    if value == 0:
        return 0
    from math import log10, floor
    magnitude = floor(log10(abs(value)))
    return round(value, -int(magnitude) + (n - 1))

## Helper function to generate summary table as HTML string
def generate_summary_table_html(data, column_list):
    """Create HTML table from data dict with specified columns.

    Args:
        data: Dict with column names as keys and lists of values
        column_list: Dict mapping column names to display titles

    Returns:
        HTML string for the table, or None if no data columns exist
    """
    # Filter to columns that exist in data
    available_cols = ["Sample"] + [col for col in column_list if col in data]
    if len(available_cols) <= 1:  # Only Sample column
        return None

    # Round float values to 2 significant figures
    formatted_data = {}
    for col in available_cols:
        formatted_data[col] = [round_to_n_sigfigs(v, 2) for v in data[col]]

    # Build HTML table
    html_parts = [
        '<table style="border-collapse: collapse; margin-bottom: 10px; background-color: white; border-radius: 5px; width: 100%;">',
        '<thead>',
        '<tr style="background-color: #f0f0f0;">'
    ]
    
    # Header row
    for col in available_cols:
        title = column_list.get(col, col) if col != "Sample" else "Sample"
        html_parts.append(f'<th style="border: 1px solid #ddd; padding: 8px; text-align: left;">{title}</th>')
    
    html_parts.append('</tr>')
    html_parts.append('</thead>')
    html_parts.append('<tbody>')
    
    # Data rows
    num_rows = len(data["Sample"])
    for i in range(num_rows):
        html_parts.append('<tr>')
        for col in available_cols:
            val = formatted_data[col][i]
            display_val = "" if val is None else str(val)
            html_parts.append(f'<td style="border: 1px solid #ddd; padding: 8px;">{display_val}</td>')
        html_parts.append('</tr>')
    
    html_parts.append('</tbody>')
    html_parts.append('</table>')
    
    return "\n".join(html_parts)

## Helper function to generate summary table (Bokeh DataTable - kept for compatibility)
def generate_summary_table(data, column_list):
    """Create Bokeh DataTable from data dict with specified columns.

    Args:
        data: Dict with column names as keys and lists of values
        column_list: List of column names to include (Sample is always prepended)

    Returns:
        Bokeh DataTable widget, or None if no data columns exist
    """
    # Filter to columns that exist in data
    available_cols = ["Sample"] + [col for col in column_list if col in data]
    if len(available_cols) <= 1:  # Only Sample column
        return None

    # Round float values to 2 significant figures
    formatted_data = {}
    for col in available_cols:
        formatted_data[col] = [round_to_n_sigfigs(v, 2) for v in data[col]]

    source = ColumnDataSource(formatted_data)

    # Create TableColumn for each column with width based on max content length
    columns = []
    for col in available_cols:
        title = column_list.get(col, col) if col != "Sample" else "Sample"
        # Find max length across title and all values
        max_len = len(title)
        for val in formatted_data[col]:
            if val is not None:
                max_len = max(max_len, len(str(val)))
        # Calculate width (~8px per character + padding)
        width = max_len * 8
        columns.append(TableColumn(field=col, title=title, width=width))

    # Calculate height based on number of rows (header ~30px + ~25px per row)
    num_rows = len(data["Sample"])
    table_height = 30 + (num_rows * 25)

    table = DataTable(
        source=source,
        columns=columns,
        sortable=True,
        reorderable=True,
        index_position=None,
        autosize_mode="none",
        sizing_mode="stretch_width",
        height=table_height
    )

    return table