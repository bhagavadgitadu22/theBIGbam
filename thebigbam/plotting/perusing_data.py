from bokeh.models import Div, ColumnDataSource, DataTable, TableColumn

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
        "Coverage_percentage": "Aligned fraction (%)", "Coverage_mean": "Coverage mean", "Coverage_median": "Coverage median",
        "Coverage_sd": "Coverage sd", "Coverage_variation": "Coverage variation", 
        "Coverage_mean_corrected_by_read_number": "Coverage mean (2)", "Coverage_median_corrected_by_read_number": "Coverage median (2)", 
        "Coverage_mean_corrected_by_read_mapped": "Coverage mean (3)", "Coverage_median_corrected_by_read_mapped": "Coverage median (3)"
    }

    completeness_cols = {
        "Percentage_completeness": "Completeness (%)", "Percentage_contamination": "Contamination (%)",
        "Total_mismatches": "Mismatches (bp)", "Total_insertions": "Insertions* (bp)", "Total_deletions": "Deletions** (bp)", 
        "Total_reads_clipped": "Read clippings (bp)", "Total_reference_clipped": "Reference clippings (bp)"
    }

    side_completeness_cols = {
        "Prevalence_completeness_left": "Left completeness (%)", "Distance_contaminated_left": "Left expansion* (bp)", "Min_missing_left": "Left collapse** (bp)",
        "Prevalence_completeness_right": "Right completeness (%)", "Distance_contaminated_right": "Right expansion* (bp)", "Min_missing_right": "Right collapse** (bp)",
    }

    phage_cols = {
        "Phage_packaging_mechanism": "Mechanism", "Phage_left_terminus": "Left termini", "Phage_right_terminus": "Right termini"
    }

    cur = conn.cursor()

    # Initialize data dict with Sample column
    data = {"Sample": list(sample_names)}
    n = len(sample_names)
    # Create sample name to index mapping
    sample_idx = {name: i for i, name in enumerate(sample_names)}

    if not sample_names:
        return  # nothing to query

    # Initialize all columns with None
    all_cols = coverage_cols | completeness_cols | side_completeness_cols | phage_cols
    for col in all_cols.keys():
        data[col] = [None] * n

    # Query Explicit_presences view
    try:
        cols_str = ", ".join(coverage_cols.keys())
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

            for value_col, cell in zip(coverage_cols.keys(), values):
                data[value_col][idx] = cell
    except Exception:
        pass  # View might not exist or have no data

    # Query Explicit_completeness view
    try:
        comp_cols = completeness_cols | side_completeness_cols
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

            for value_col, cell in zip(comp_cols.keys(), values):
                data[value_col][idx] = cell
    except Exception:
        pass  # View might not exist or have no data

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

            for value_col, cell in zip(phage_cols.keys(), values):
                data[value_col][idx] = cell
    except Exception:
        pass  # View might not exist or have no data

    # Create content
    content = []

    # Coverage subsection
    coverage_table = generate_summary_table(data, coverage_cols)
    if coverage_table:
        coverage_header = Div(text="<b>Coverage:</b>")
        content.append(coverage_header)
        content.append(coverage_table)
        # Add explanation for corrected coverage columns
        correction_explanation = Div(
            text="<i>(2) Coverage mean/median corrected by number of reads or number of reads mapped<br>"
                 "(3) Coverage mean/median corrected by number of reads or number of reads mapped</i>",
            margin=(0, 0, 10, 0)
        )
        content.append(correction_explanation)

    # Completeness subsection
    completeness_table = generate_summary_table(data, completeness_cols)
    side_completeness_table = generate_summary_table(data, side_completeness_cols)
    if completeness_table or side_completeness_table:
        completeness_header = Div(text="<b>Completeness:</b>")
        content.append(completeness_header)
        content.append(completeness_table)
        # Add explanation for expansion and collapse
        explanations = Div(
            text="<i>* Insertion: Extra bases are present in the read but not in the contig<br>"
                 "** Deletion: A stretch of the contig has no corresponding bases in the read</i>",
            margin=(0, 0, 10, 0)
        )
        content.append(explanations)

        content.append(side_completeness_table)
        # Add explanation for expansion and collapse
        explanations = Div(
            text="<i>* Expansion: alignment gaps in the aligned reads (ie extra sequences in the contigs)<br>"
                 "** Collapse: extra sequences in the aligned reads (ie alignment gaps in the contigs)</i>",
            margin=(0, 0, 10, 0)
        )
        content.append(explanations)

    # Phage mechanism subsection
    phage_table = generate_summary_table(data, phage_cols)
    if phage_table:
        phage_header = Div(text="<b>Phage mechanism:</b>")
        content.append(phage_header)
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

## Helper function to generate summary table
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