import argparse
import base64
import os
import duckdb
import traceback

import panel as pn

from bokeh.layouts import column, row
from bokeh.models import Div, InlineStyleSheet, Tooltip
from bokeh.models.widgets import CheckboxGroup, HelpButton, Button, RadioButtonGroup, CheckboxButtonGroup, Select, TextInput, Spinner
from bokeh.models.plots import GridPlot

# Import the plotting function from the repo
from .plotting_data_per_sample import generate_bokeh_plot_per_sample
from .plotting_data_all_samples import generate_bokeh_plot_all_samples
from .perusing_data import build_summary_data, generate_summary_table
from ..database.database_getters import get_filtering_metadata

def build_controls(conn):

    """Query DB and return widgets and helper mappings."""
    cur = conn.cursor()

    # Check if PhageMechanisms table exists and has data
    phage_mechanisms_list = []
    try:
        cur.execute("SELECT DISTINCT Phage_packaging_mechanism FROM PhageMechanisms WHERE Phage_packaging_mechanism IS NOT NULL")
        phage_mechanisms_list = [r[0] for r in cur.fetchall()]
    except Exception:
        pass  # Table doesn't exist or has no data

    # Check if Completeness table exists and has data
    has_completeness = False
    try:
        cur.execute("SELECT 1 FROM Completeness LIMIT 1")
        has_completeness = cur.fetchone() is not None
    except Exception:
        pass

    # Get Coverage_mean min/max from Presences table (stored as INTEGER)
    coverage_mean_max = 0
    try:
        cur.execute("SELECT MIN(Coverage_mean), MAX(Coverage_mean) FROM Presences WHERE Coverage_mean IS NOT NULL")
        result = cur.fetchone()
        if result and result[0] is not None and result[1] is not None:
            coverage_mean_min = result[0]
            coverage_mean_max = result[1]
    except Exception:
        pass

    # Get Coverage_median min/max from Presences table (stored as INTEGER)
    coverage_median_max = 0
    try:
        cur.execute("SELECT MIN(Coverage_median), MAX(Coverage_median) FROM Presences WHERE Coverage_median IS NOT NULL")
        result = cur.fetchone()
        if result and result[0] is not None and result[1] is not None:
            coverage_median_min = result[0]
            coverage_median_max = result[1]
    except Exception:
        pass

    # Get Coverage_variation min/max from Explicit_presences view (already divided by 1000000)
    coverage_variation_max = 1.0
    try:
        cur.execute("SELECT MIN(Coverage_variation), MAX(Coverage_variation) FROM Explicit_presences WHERE Coverage_variation IS NOT NULL")
        result = cur.fetchone()
        if result and result[0] is not None and result[1] is not None:
            coverage_variation_min = result[0]
            coverage_variation_max = result[1]
    except Exception:
        pass

    # Get Coverage_sd min/max from Explicit_presences view (already divided by 1000000)
    coverage_sd_max = 1.0
    try:
        cur.execute("SELECT MIN(Coverage_sd), MAX(Coverage_sd) FROM Explicit_presences WHERE Coverage_sd IS NOT NULL")
        result = cur.fetchone()
        if result and result[0] is not None and result[1] is not None:
            coverage_sd_min = result[0]
            coverage_sd_max = result[1]
    except Exception:
        pass

    whole_contamination_max = 100.0
    try:
        cur.execute("SELECT MAX(Percentage_contamination) FROM Explicit_completeness WHERE Percentage_contamination IS NOT NULL")
        result = cur.fetchone()
        if result and result[0] is not None:
            whole_contamination_max = result[0]
    except Exception:
        pass

    # Get Duplication_percentage min/max from Contig table
    duplication_percentage_min = 0
    duplication_percentage_max = 100
    has_duplication_data = False
    try:
        cur.execute("SELECT MIN(Duplication_percentage), MAX(Duplication_percentage) FROM Contig WHERE Duplication_percentage IS NOT NULL")
        result = cur.fetchone()
        if result and result[0] is not None and result[1] is not None:
            duplication_percentage_min = result[0]
            duplication_percentage_max = result[1]
            has_duplication_data = True
    except Exception:
        pass

    # Get annotation columns and their distinct values from Contig_annotation table
    annotation_filters = {}
    excluded_columns = {'Contig_id', 'Start', 'End', 'Strand', 'Type'}
    try:
        # Get column names from Contig_annotation table
        cur.execute("PRAGMA table_info(Contig_annotation)")
        columns = [row[1] for row in cur.fetchall() if row[1] not in excluded_columns]
        
        # Get distinct non-null values for each column
        for column in columns:
            cur.execute(f'SELECT DISTINCT "{column}" FROM Contig_annotation WHERE "{column}" IS NOT NULL ORDER BY "{column}"')
            values = [str(r[0]) for r in cur.fetchall()]  # Convert to strings for AutocompleteInput
            if values:  # Only add if there are values
                annotation_filters[column] = values
    except Exception:
        pass  # Table doesn't exist or has no data

    # Widget Selector for Contigs (autocomplete with max 20 suggestions)
    cur.execute("SELECT Contig_name, Contig_length, Duplication_percentage FROM Contig ORDER BY Contig_name")
    rows = cur.fetchall()
    contigs = [r[0] for r in rows]
    contig_lengths = {r[0]: r[1] for r in rows}  # Dictionary mapping contig_name -> length
    contig_duplications = {r[0]: r[2] for r in rows}  # Dictionary mapping contig_name -> duplication percentage (can be None)
    
    # If only one contig in database, pre-fill the field
    contig_select = pn.widgets.AutocompleteInput(
        value=contigs[0] if len(contigs) == 1 else "",
        options=contigs,
        min_characters=0,
        case_sensitive=False,
        placeholder="Type to search contigs...",
        search_strategy="includes",
        sizing_mode="stretch_width"
    )

    # Widget Selector for Samples (autocomplete with max 20 suggestions)
    cur.execute("SELECT Sample_name FROM Sample ORDER BY Sample_name")
    samples = [r[0] for r in cur.fetchall()]
    
    # If only one sample in database, pre-fill the field
    sample_select = pn.widgets.AutocompleteInput(
        value=samples[0] if len(samples) == 1 else "",
        options=samples,
        min_characters=0,
        case_sensitive=False,
        placeholder="Type to search samples...",
        search_strategy="includes",
        sizing_mode="stretch_width"
    )

    # Build presence mappings: sample -> contigs and contig -> samples
    cur.execute("""
    SELECT Contig.Contig_name, Sample.Sample_name FROM Presences
      JOIN Contig ON Presences.Contig_id = Contig.Contig_id
      JOIN Sample ON Presences.Sample_id = Sample.Sample_id
    """)
    sample_to_contigs = {}
    contig_to_samples = {}
    for contig_name, sample_name in cur.fetchall():
        sample_to_contigs.setdefault(sample_name, set()).add(contig_name)
        contig_to_samples.setdefault(contig_name, set()).add(sample_name)

    # Modules and variables - only show those with data in database
    # First, identify which feature tables have at least one row with non-zero value
    cur.execute("SELECT DISTINCT Feature_table_name FROM Variable")
    feature_tables = [r[0] for r in cur.fetchall()]

    tables_with_data = set()
    for table_name in feature_tables:
        tables_with_data.add(table_name)

    # Check if repeat tables have data
    # Add table names to tables_with_data so repeat variables are included
    try:
        cur.execute("SELECT 1 FROM Contig_DirectRepeats LIMIT 1")
        if cur.fetchone() is not None:
            tables_with_data.add("Contig_DirectRepeats")
    except Exception:
        pass

    try:
        cur.execute("SELECT 1 FROM Contig_InvertedRepeats LIMIT 1")
        if cur.fetchone() is not None:
            tables_with_data.add("Contig_InvertedRepeats")
    except Exception:
        pass

    # Get variables that have data (their feature table has rows)
    cur.execute("SELECT DISTINCT Variable_name, Feature_table_name FROM Variable")
    variables = [r[0] for r in cur.fetchall() if r[1] in tables_with_data]

    # Get modules that have at least one variable with data
    # Define the display order for modules
    MODULE_ORDER = ["Genome", "Coverage", "Misalignment", "Long-reads", "Paired-reads", "Phage termini"]

    cur.execute("SELECT DISTINCT Module FROM Variable WHERE Feature_table_name IN ({})".format(
        ','.join('?' * len(tables_with_data))
    ), tuple(tables_with_data))
    modules_from_db = [r[0] for r in cur.fetchall()]

    # Sort modules according to MODULE_ORDER, keeping any unknown modules at the end
    modules = sorted(modules_from_db, key=lambda m: MODULE_ORDER.index(m) if m in MODULE_ORDER else len(MODULE_ORDER))

    # For each module get variables (only those with data)
    # Create TWO sets of widgets: one for "One Sample" view, one for "All Samples" view
    module_names = []
    module_widgets_one = []  # Module checkboxes for One Sample view
    variables_widgets_one = []  # Variable button groups for One Sample view
    variables_widgets_all = []  # Variable button groups for All Samples view
    helps_widgets = []
    variables_labels = []  # Store labels for each module
    for module in modules:
        # Get distinct subplots (deduplicate by subplot name only)
        cur.execute(
            "SELECT DISTINCT Subplot FROM Variable WHERE Module=? AND Feature_table_name IN ({}) ORDER BY Module_order".format(
                ','.join('?' * len(tables_with_data))
            ),
            (module,) + tuple(tables_with_data)
        )
        variables_checkbox = [r[0] for r in cur.fetchall()]

        # For Genome module, add "Gene map" as the first option
        if module == "Genome":
            variables_checkbox = ["Gene map"] + variables_checkbox

        # Skip this module if no variables with data
        if not variables_checkbox:
            continue

        module_names.append(module)
        variables_labels.append(variables_checkbox)

        # Module checkbox (for One Sample view only)
        module_checkbox = CheckboxGroup(labels=[module], active=[])
        module_widgets_one.append(module_checkbox)

        # CheckboxButtonGroup for One Sample view
        cbg_one = CheckboxButtonGroup(labels=variables_checkbox, active=[], sizing_mode="stretch_width", orientation="vertical")
        variables_widgets_one.append(cbg_one)

        # CheckboxButtonGroup for All Samples view (separate instance)
        cbg_all = CheckboxButtonGroup(labels=variables_checkbox, active=[], sizing_mode="stretch_width", orientation="vertical")
        variables_widgets_all.append(cbg_all)

        # Consolidate help texts for the module into a single HelpButton attached to module title
        combined_help = ""
        cur.execute(
            "SELECT DISTINCT Subplot, Title, Help FROM Variable WHERE Module=? AND Feature_table_name IN ({}) ORDER BY Subplot".format(
                ','.join('?' * len(tables_with_data))
            ),
            (module,) + tuple(tables_with_data)
        )
        records = cur.fetchall()
        for subplot, title, help_text in records:
            if help_text is None or not help_text.strip():
                continue  # Skip empty or None help texts
            combined_help += f"{title} ({subplot} subplot): {help_text}\n"

        if combined_help:
            tooltip = Tooltip(content=combined_help, position="right")
        else:
            tooltip = None
        helps_widgets.append(tooltip)

    widgets = {
        'sample_select': sample_select,
        'contig_select': contig_select,
        'sample_to_contigs': sample_to_contigs,
        'contig_to_samples': contig_to_samples,
        'module_names': module_names,
        'module_widgets_one': module_widgets_one,
        'helps_widgets': helps_widgets,
        'variables_widgets_one': variables_widgets_one,
        'variables_widgets_all': variables_widgets_all,
        'variables_labels': variables_labels,
        'contigs': contigs,
        'contig_lengths': contig_lengths,
        'contig_duplications': contig_duplications,
        'samples': samples,
        'variables': variables,
        'phage_mechanisms': phage_mechanisms_list,
        'has_completeness': has_completeness,
        'coverage_mean_max': coverage_mean_max,
        'coverage_median_max': coverage_median_max,
        'coverage_variation_max': coverage_variation_max,
        'coverage_sd_max': coverage_sd_max,
        'whole_contamination_max': whole_contamination_max,
        'duplication_percentage_min': duplication_percentage_min,
        'duplication_percentage_max': duplication_percentage_max,
        'has_duplication_data': has_duplication_data,
        'annotation_filters': annotation_filters,
        'full_contigs': contigs,  # Store full list for substring matching
        'full_samples': samples   # Store full list for substring matching
    }
    return widgets

def create_layout(db_path):
    """Create and return the application layout for Panel serve."""

    ### Event functions
    ## Helper function to create Panel EditableRangeSlider
    def create_editable_range_slider(title, start, end, value, step, on_change=None):
        """Create a Panel EditableRangeSlider.

        Args:
            title: Slider title/name
            start: Minimum value
            end: Maximum value
            value: Tuple of (min_value, max_value)
            step: Step size
            on_change: Optional callback to call when value changes

        Returns:
            Panel EditableRangeSlider widget
        """
        slider = pn.widgets.EditableRangeSlider(
            name=title,
            start=start,
            end=end,
            value=value,
            step=step,
            sizing_mode="stretch_width",
            margin=(0, 10, 0, 0)
        )
        if on_change:
            slider.param.watch(on_change, 'value_throttled')
        return slider

    ## Helper function for substring matching in autocomplete
    def substring_filter(items, query, max_results=20):
        """Filter items by substring match (case-insensitive)."""
        if not query:
            return [""] + items[:max_results]
        query_lower = query.lower()
        matches = [item for item in items if query_lower in item.lower()]
        return [""] + matches[:max_results]

    ## Helper function to create collapsible section toggle callbacks
    def make_toggle_callback(btn, content):
        def callback():
            content.visible = not content.visible
            if content.visible:
                btn.label = "▼"
            else:
                btn.label = "▶"
        return callback
    
    # Initialize annotation_inputs early so it's accessible in filter functions
    annotation_inputs = {}
    
    ## Helper functions to refresh completions based on filters
    def query_summary_filtered_items(entity_type, var_name, comparison, threshold, filter_type="#Points"):
        """Query Summary table for filtered contig or sample names.
        
        Args:
            entity_type: "Contig" or "Sample"
            var_name: Variable name to filter by
            comparison: ">" or "<"
            threshold: Row count or max value threshold
            filter_type: "#Points" for row count or "Max" for maximum value
        
        Returns:
            Set of matching entity names
        """
        cur = conn.cursor()
        operator = ">" if comparison == ">" else "<"
        
        if filter_type == "#Points":
            # Original behavior: filter by number of data points (Row_count)
            query = f"""
                SELECT DISTINCT {entity_type}.{entity_type}_name 
                FROM Summary
                JOIN {entity_type} ON Summary.{entity_type}_id = {entity_type}.{entity_type}_id
                JOIN Variable ON Summary.Variable_id = Variable.Variable_id
                WHERE Variable.Variable_name = ? AND Summary.Row_count {operator} ?
            """
            cur.execute(query, (var_name, threshold))
        else:  # filter_type == "Max"
            # New behavior: filter by maximum value in the feature table
            # Get the feature table name for this variable
            cur.execute("SELECT Feature_table_name FROM Variable WHERE Variable_name = ? LIMIT 1", (var_name,))
            result = cur.fetchone()
            if not result:
                return set()
            
            feature_table = result[0]
            query = f"""
                SELECT DISTINCT {entity_type}.{entity_type}_name
                FROM {feature_table} ft
                JOIN {entity_type} ON ft.{entity_type}_id = {entity_type}.{entity_type}_id
                WHERE ft.Value {operator} ?
            """
            cur.execute(query, (threshold,))
        
        return {row[0] for row in cur.fetchall()}
    
    def get_variable_filtered_contigs():
        """Apply variable-based filters to get allowed contigs."""
        if not variable_filter_rows:
            return set(orig_contigs)
        
        allowed_contigs = set(orig_contigs)
        
        for filter_row in variable_filter_rows:
            filter_type = filter_row.children[0].value  # "#Points" or "Max"
            var_name = filter_row.children[1].value
            comparison = filter_row.children[2].value
            threshold_str = filter_row.children[3].value
            
            # Skip invalid filters - only apply if variable name is valid and exists
            if not var_name or not threshold_str or var_name not in widgets['variables']:
                continue
            
            try:
                threshold = int(threshold_str) if filter_type == "#Points" else float(threshold_str)
                if threshold < 0:
                    continue
            except ValueError:
                continue
            
            matching = query_summary_filtered_items("Contig", var_name, comparison, threshold, filter_type)
            allowed_contigs &= matching  # AND logic: all filters must match
        
        return allowed_contigs
    
    def get_variable_filtered_samples():
        """Apply variable-based filters to get allowed samples."""
        if not variable_filter_rows:
            return set(orig_samples)
        
        allowed_samples = set(orig_samples)
        
        for filter_row in variable_filter_rows:
            filter_type = filter_row.children[0].value  # "#Points" or "Max"
            var_name = filter_row.children[1].value
            comparison = filter_row.children[2].value
            threshold_str = filter_row.children[3].value
            
            # Skip invalid filters - only apply if variable name is valid and exists
            if not var_name or not threshold_str or var_name not in widgets['variables']:
                continue
            
            try:
                threshold = int(threshold_str) if filter_type == "#Points" else float(threshold_str)
                if threshold < 0:
                    continue
            except ValueError:
                continue
            
            matching = query_summary_filtered_items("Sample", var_name, comparison, threshold, filter_type)
            allowed_samples &= matching  # AND logic: all filters must match
        
        return allowed_samples
    
    def get_variable_filtered_contigs_any_sample():
        """Apply variable-based filters to get contigs where at least one sample meets criteria.
        
        Used in All samples view to suggest contigs that have sufficient data in at least one sample.
        """
        if not variable_filter_rows:
            return set(orig_contigs)
        
        allowed_contigs = set(orig_contigs)
        
        for filter_row in variable_filter_rows:
            filter_type = filter_row.children[0].value  # "#Points" or "Max"
            var_name = filter_row.children[1].value
            comparison = filter_row.children[2].value
            threshold_str = filter_row.children[3].value
            
            # Skip invalid filters
            if not var_name or not threshold_str or var_name not in widgets['variables']:
                continue
            
            try:
                threshold = int(threshold_str) if filter_type == "#Points" else float(threshold_str)
                if threshold < 0:
                    continue
            except ValueError:
                continue
            
            cur = conn.cursor()
            operator = ">" if comparison == ">" else "<"
            
            if filter_type == "#Points":
                # Query Summary for contigs where at least one sample meets the threshold
                query = f"""
                    SELECT DISTINCT Contig.Contig_name 
                    FROM Summary
                    JOIN Contig ON Summary.Contig_id = Contig.Contig_id
                    JOIN Variable ON Summary.Variable_id = Variable.Variable_id
                    WHERE Variable.Variable_name = ? AND Summary.Row_count {operator} ?
                """
                cur.execute(query, (var_name, threshold))
            else:  # filter_type == "Max"
                # Query feature table for contigs where value exceeds threshold
                cur.execute("SELECT Feature_table_name FROM Variable WHERE Variable_name = ? LIMIT 1", (var_name,))
                result = cur.fetchone()
                if not result:
                    continue
                
                feature_table = result[0]
                query = f"""
                    SELECT DISTINCT Contig.Contig_name
                    FROM {feature_table} ft
                    JOIN Contig ON ft.Contig_id = Contig.Contig_id
                    WHERE ft.Value {operator} ?
                """
                cur.execute(query, (threshold,))
            
            matching = {row[0] for row in cur.fetchall()}
            allowed_contigs &= matching  # AND logic: all filters must match
        
        return allowed_contigs

    ## Helper function to get coverage mean column based on correction selection
    def get_coverage_mean_column():
        """Return appropriate Coverage_mean column based on correction_filter selection."""
        try:
            if correction_filter is not None and hasattr(correction_filter, 'value'):
                val = correction_filter.value
                if val == "Correct by number of reads":
                    return "Coverage_mean_corrected_by_read_number"
                elif val == "Correct by number of mapped reads":
                    return "Coverage_mean_corrected_by_read_mapped"
        except NameError:
            pass
        return "Coverage_mean"

    ## Helper function to get coverage median column based on correction selection
    def get_coverage_median_column():
        """Return appropriate Coverage_median column based on correction_filter selection."""
        try:
            if correction_filter is not None and hasattr(correction_filter, 'value'):
                val = correction_filter.value
                if val == "Correct by number of reads":
                    return "Coverage_median_corrected_by_read_number"
                elif val == "Correct by number of mapped reads":
                    return "Coverage_median_corrected_by_read_mapped"
        except NameError:
            pass
        return "Coverage_median"

    ## Helper functions for module-based filtering (phage mechanisms, completeness, coverage)
    def get_module_filtered_contigs(sample_name=None):
        """Apply module-based filters (coverage, completeness, phage mechanism) to get allowed contigs.

        Args:
            sample_name: Optional sample name to filter by. If provided, only considers
                         filter values for that specific sample.
        """
        allowed_contigs = set(orig_contigs)
        cur = conn.cursor()

        # Apply coverage filters using Explicit_presences view (single combined query)
        conditions = []
        params = []

        # Add sample filter if specified
        if sample_name:
            conditions.append("Sample_name = ?")
            params.append(sample_name)

        if coverage_percentage_slider is not None:
            min_val, max_val = coverage_percentage_slider.value
            if min_val > 0 or max_val < 100:
                conditions.append("Coverage_percentage >= ? AND Coverage_percentage <= ?")
                params.extend([min_val, max_val])

        if coverage_mean_slider is not None:
            min_val, max_val = coverage_mean_slider.value
            cov_mean_max = coverage_mean_slider.end
            if min_val > 0 or max_val < cov_mean_max:
                cov_col = get_coverage_mean_column()
                conditions.append(f"{cov_col} >= ? AND {cov_col} <= ?")
                params.extend([min_val, max_val])

        if coverage_median_slider is not None:
            min_val, max_val = coverage_median_slider.value
            cov_median_max = coverage_median_slider.end
            if min_val > 0 or max_val < cov_median_max:
                cov_col = get_coverage_median_column()
                conditions.append(f"{cov_col} >= ? AND {cov_col} <= ?")
                params.extend([min_val, max_val])

        if coverage_variation_slider is not None:
            min_val, max_val = coverage_variation_slider.value
            cov_var_max = widgets['coverage_variation_max']
            if min_val > 0 or max_val < cov_var_max:
                conditions.append("Coverage_variation >= ? AND Coverage_variation <= ?")
                params.extend([min_val, max_val])

        if coverage_sd_slider is not None:
            min_val, max_val = coverage_sd_slider.value
            cov_sd_max = widgets['coverage_sd_max']
            if min_val > 0 or max_val < cov_sd_max:
                conditions.append("Coverage_sd >= ? AND Coverage_sd <= ?")
                params.extend([min_val, max_val])

        if conditions:
            query = f"SELECT DISTINCT Contig_name FROM Explicit_presences WHERE {' AND '.join(conditions)}"
            cur.execute(query, params)
            matching = {row[0] for row in cur.fetchall()}
            allowed_contigs &= matching

        # Apply completeness filters using Explicit_completeness view (single combined query)
        conditions = []
        params = []

        # Add sample filter if specified
        if sample_name:
            conditions.append("Sample_name = ?")
            params.append(sample_name)

        if prevalence_left_slider is not None:
            min_val, max_val = prevalence_left_slider.value
            if min_val > 0 or max_val < 100:
                conditions.append("Prevalence_completeness_left >= ? AND Prevalence_completeness_left <= ?")
                params.extend([min_val, max_val])

        if prevalence_right_slider is not None:
            min_val, max_val = prevalence_right_slider.value
            if min_val > 0 or max_val < 100:
                conditions.append("Prevalence_completeness_right >= ? AND Prevalence_completeness_right <= ?")
                params.extend([min_val, max_val])

        if pct_completeness_slider is not None:
            min_val, max_val = pct_completeness_slider.value
            if min_val > 0 or max_val < 100:
                conditions.append("Percentage_completeness >= ? AND Percentage_completeness <= ?")
                params.extend([min_val, max_val])

        if pct_contamination_slider is not None:
            min_val, max_val = pct_contamination_slider.value
            conta_max = widgets['whole_contamination_max']
            if min_val > 0 or max_val < conta_max:
                conditions.append("Percentage_contamination >= ? AND Percentage_contamination <= ?")
                params.extend([min_val, max_val])

        if circularising_slider is not None:
            min_val, max_val = circularising_slider.value
            if min_val > 0 or max_val < 100:
                conditions.append("Circularising_reads_percentage >= ? AND Circularising_reads_percentage <= ?")
                params.extend([min_val, max_val])

        if conditions:
            query = f"SELECT DISTINCT Contig_name FROM Explicit_completeness WHERE {' AND '.join(conditions)}"
            cur.execute(query, params)
            matching = {row[0] for row in cur.fetchall()}
            allowed_contigs &= matching

        # Apply phage mechanism filter using Explicit_phage_mechanisms view
        if phage_mechanism_filter is not None and phage_mechanism_filter.value:
            if sample_name:
                query = """
                    SELECT DISTINCT Contig_name
                    FROM Explicit_phage_mechanisms
                    WHERE Phage_packaging_mechanism = ? AND Sample_name = ?
                """
                cur.execute(query, (phage_mechanism_filter.value, sample_name))
            else:
                query = """
                    SELECT DISTINCT Contig_name
                    FROM Explicit_phage_mechanisms
                    WHERE Phage_packaging_mechanism = ?
                """
                cur.execute(query, (phage_mechanism_filter.value,))
            matching = {row[0] for row in cur.fetchall()}
            allowed_contigs &= matching

        # Apply annotation filters
        for column_name, input_widget in annotation_inputs.items():
            if input_widget.value:
                query = f'''
                    SELECT DISTINCT Contig.Contig_name 
                    FROM Contig_annotation
                    JOIN Contig ON Contig_annotation.Contig_id = Contig.Contig_id
                    WHERE "{column_name}" = ?
                '''
                cur.execute(query, (input_widget.value,))
                matching = {row[0] for row in cur.fetchall()}
                allowed_contigs &= matching

        return allowed_contigs

    def get_module_filtered_samples(contig_name=None):
        """Apply module-based filters (coverage, completeness, phage mechanism) to get allowed samples.

        Args:
            contig_name: Optional contig name to filter by. If provided, only considers
                         filter values for that specific contig.
        """
        allowed_samples = set(orig_samples)
        cur = conn.cursor()

        # Apply coverage filters using Explicit_presences view (single combined query)
        conditions = []
        params = []

        # Add contig filter if specified
        if contig_name:
            conditions.append("Contig_name = ?")
            params.append(contig_name)

        if coverage_percentage_slider is not None:
            min_val, max_val = coverage_percentage_slider.value
            if min_val > 0 or max_val < 100:
                conditions.append("Coverage_percentage >= ? AND Coverage_percentage <= ?")
                params.extend([min_val, max_val])

        if coverage_mean_slider is not None:
            min_val, max_val = coverage_mean_slider.value
            cov_mean_max = coverage_mean_slider.end
            if min_val > 0 or max_val < cov_mean_max:
                cov_col = get_coverage_mean_column()
                conditions.append(f"{cov_col} >= ? AND {cov_col} <= ?")
                params.extend([min_val, max_val])

        if coverage_median_slider is not None:
            min_val, max_val = coverage_median_slider.value
            cov_median_max = coverage_median_slider.end
            if min_val > 0 or max_val < cov_median_max:
                cov_col = get_coverage_median_column()
                conditions.append(f"{cov_col} >= ? AND {cov_col} <= ?")
                params.extend([min_val, max_val])

        if coverage_variation_slider is not None:
            min_val, max_val = coverage_variation_slider.value
            cov_var_max = widgets['coverage_variation_max']
            if min_val > 0 or max_val < cov_var_max:
                conditions.append("Coverage_variation >= ? AND Coverage_variation <= ?")
                params.extend([min_val, max_val])

        if coverage_sd_slider is not None:
            min_val, max_val = coverage_sd_slider.value
            cov_sd_max = widgets['coverage_sd_max']
            if min_val > 0 or max_val < cov_sd_max:
                conditions.append("Coverage_sd >= ? AND Coverage_sd <= ?")
                params.extend([min_val, max_val])

        if conditions:
            query = f"SELECT DISTINCT Sample_name FROM Explicit_presences WHERE {' AND '.join(conditions)}"
            cur.execute(query, params)
            matching = {row[0] for row in cur.fetchall()}
            allowed_samples &= matching

        # Apply completeness filters using Explicit_completeness view (single combined query)
        conditions = []
        params = []

        # Add contig filter if specified
        if contig_name:
            conditions.append("Contig_name = ?")
            params.append(contig_name)

        if prevalence_left_slider is not None:
            min_val, max_val = prevalence_left_slider.value
            if min_val > 0 or max_val < 100:
                conditions.append("Prevalence_completeness_left >= ? AND Prevalence_completeness_left <= ?")
                params.extend([min_val, max_val])

        if prevalence_right_slider is not None:
            min_val, max_val = prevalence_right_slider.value
            if min_val > 0 or max_val < 100:
                conditions.append("Prevalence_completeness_right >= ? AND Prevalence_completeness_right <= ?")
                params.extend([min_val, max_val])

        if pct_completeness_slider is not None:
            min_val, max_val = pct_completeness_slider.value
            if min_val > 0 or max_val < 100:
                conditions.append("Percentage_completeness >= ? AND Percentage_completeness <= ?")
                params.extend([min_val, max_val])

        if pct_contamination_slider is not None:
            min_val, max_val = pct_contamination_slider.value
            conta_max = widgets['whole_contamination_max']
            if min_val > 0 or max_val < conta_max:
                conditions.append("Percentage_contamination >= ? AND Percentage_contamination <= ?")
                params.extend([min_val, max_val])

        if circularising_slider is not None:
            min_val, max_val = circularising_slider.value
            if min_val > 0 or max_val < 100:
                conditions.append("Circularising_reads_percentage >= ? AND Circularising_reads_percentage <= ?")
                params.extend([min_val, max_val])

        if conditions:
            query = f"SELECT DISTINCT Sample_name FROM Explicit_completeness WHERE {' AND '.join(conditions)}"
            cur.execute(query, params)
            matching = {row[0] for row in cur.fetchall()}
            allowed_samples &= matching

        # Apply phage mechanism filter using Explicit_phage_mechanisms view
        if phage_mechanism_filter is not None and phage_mechanism_filter.value:
            if contig_name:
                query = """
                    SELECT DISTINCT Sample_name
                    FROM Explicit_phage_mechanisms
                    WHERE Phage_packaging_mechanism = ? AND Contig_name = ?
                """
                cur.execute(query, (phage_mechanism_filter.value, contig_name))
            else:
                query = """
                    SELECT DISTINCT Sample_name
                    FROM Explicit_phage_mechanisms
                    WHERE Phage_packaging_mechanism = ?
                """
                cur.execute(query, (phage_mechanism_filter.value,))
            matching = {row[0] for row in cur.fetchall()}
            allowed_samples &= matching

        # Apply annotation filters (samples that have the contig with matching annotations)
        for column_name, input_widget in annotation_inputs.items():
            if input_widget.value:
                if contig_name:
                    # Filter samples that have this specific contig with the annotation
                    query = f'''
                        SELECT DISTINCT Sample.Sample_name
                        FROM Contig_annotation
                        JOIN Contig ON Contig_annotation.Contig_id = Contig.Contig_id
                        JOIN Presences ON Contig.Contig_id = Presences.Contig_id
                        JOIN Sample ON Presences.Sample_id = Sample.Sample_id
                        WHERE "{column_name}" = ? AND Contig.Contig_name = ?
                    '''
                    cur.execute(query, (input_widget.value, contig_name))
                else:
                    # Filter samples that have any contig with the annotation
                    query = f'''
                        SELECT DISTINCT Sample.Sample_name
                        FROM Contig_annotation
                        JOIN Contig ON Contig_annotation.Contig_id = Contig.Contig_id
                        JOIN Presences ON Contig.Contig_id = Presences.Contig_id
                        JOIN Sample ON Presences.Sample_id = Sample.Sample_id
                        WHERE "{column_name}" = ?
                    '''
                    cur.execute(query, (input_widget.value,))
                matching = {row[0] for row in cur.fetchall()}
                allowed_samples &= matching

        return allowed_samples

    def get_filtering_filtered_pairs():
        """Apply Filtering query rows to get allowed contig/sample pairs.

        Returns set of (contig_name, sample_name) tuples that match all conditions.
        Returns None if no filters are active (meaning all pairs are allowed).
        """
        if not or_sections:
            return None

        # Check if any filter has a meaningful value
        has_active_filter = False
        for section_data in or_sections:
            for row_data in section_data['rows']:
                input_ref = row_data['input_ref']
                value = input_ref['widget'].value
                # Skip empty values
                if value is None or value == "":
                    continue
                if input_ref['is_panel'] and isinstance(value, str) and value.strip() == "":
                    continue
                if not input_ref['is_panel'] and value == 0:
                    continue
                has_active_filter = True
                break
            if has_active_filter:
                break

        if not has_active_filter:
            return None

        cur = conn.cursor()

        # Source table mapping from filtering_metadata
        source_table_map = {
            'Contig': 'Contig',
            'Sample': 'Sample',
            'Presences': 'Explicit_presences',
            'Completeness': 'Explicit_completeness',
            'Termini': 'Explicit_phage_mechanisms'
        }

        def get_pairs_for_condition(category, column_name, operator, value):
            """Query database for contig/sample pairs matching a single condition."""
            source = source_table_map.get(category)
            if not source:
                return set()

            # Check if this column is from Contig_annotation table
            col_info = filtering_metadata.get(category, {}).get('columns', {}).get(column_name, {})
            col_source = col_info.get('source')

            if col_source == 'Contig_annotation':
                # Annotation column - join with Contig_annotation table
                query = f'''
                    SELECT DISTINCT c.Contig_name, s.Sample_name
                    FROM Contig_annotation ca
                    JOIN Contig c ON ca.Contig_id = c.Contig_id
                    CROSS JOIN Sample s
                    JOIN Presences p ON c.Contig_id = p.Contig_id AND s.Sample_id = p.Sample_id
                    WHERE ca."{column_name}" {operator} ?
                '''
            elif category == 'Contig':
                # Contig table has no Sample_name, cross-join with all samples
                query = f'''
                    SELECT DISTINCT c.Contig_name, s.Sample_name
                    FROM Contig c
                    CROSS JOIN Sample s
                    JOIN Presences p ON c.Contig_id = p.Contig_id AND s.Sample_id = p.Sample_id
                    WHERE c."{column_name}" {operator} ?
                '''
            elif category == 'Sample':
                # Sample table has no Contig_name, cross-join with all contigs
                query = f'''
                    SELECT DISTINCT c.Contig_name, s.Sample_name
                    FROM Sample s
                    CROSS JOIN Contig c
                    JOIN Presences p ON c.Contig_id = p.Contig_id AND s.Sample_id = p.Sample_id
                    WHERE s."{column_name}" {operator} ?
                '''
            else:
                # Tables with both Contig_name and Sample_name (views)
                query = f'''
                    SELECT DISTINCT Contig_name, Sample_name
                    FROM {source}
                    WHERE "{column_name}" {operator} ?
                '''

            try:
                cur.execute(query, [value])
                return {(row[0], row[1]) for row in cur.fetchall()}
            except Exception as e:
                print(f"[get_filtering_filtered_pairs] Query error: {e}")
                return set()

        def evaluate_section(section_data):
            """Evaluate all rows in a section using AND/OR logic based on Select widgets."""
            if not section_data['rows']:
                return None

            result_pairs = None

            for i, row_data in enumerate(section_data['rows']):
                category = row_data['category_select'].value
                column_name = row_data['subcategory_select'].value
                operator = row_data['comparison_select'].value
                input_ref = row_data['input_ref']

                value = input_ref['widget'].value

                # Skip rows with no value
                if value is None or value == "":
                    continue
                if input_ref['is_panel'] and isinstance(value, str) and value.strip() == "":
                    continue
                if not input_ref['is_panel'] and value == 0:
                    continue

                pairs = get_pairs_for_condition(category, column_name, operator, value)

                if result_pairs is None:
                    result_pairs = pairs
                else:
                    # Check the AND/OR select widget for this row
                    and_div = row_data.get('and_div')
                    if and_div is not None and and_div.value == "OR":
                        result_pairs = result_pairs | pairs  # OR logic
                    else:
                        result_pairs = result_pairs & pairs  # AND logic (default)

            return result_pairs

        # Evaluate all sections using AND/OR logic between sections
        final_pairs = None

        for i, section_data in enumerate(or_sections):
            section_pairs = evaluate_section(section_data)

            if section_pairs is None:
                continue

            if final_pairs is None:
                final_pairs = section_pairs
            else:
                # Check inter-section AND/OR - look at the Select widget between sections
                # The inter-section Select widgets are in filtering_content.children
                # They appear before each section except the first
                # For now, default to AND between sections
                final_pairs = final_pairs & section_pairs

        return final_pairs

    def update_widget_completions(widget, completions):
        """Update widget completions. Clear value if not in completions."""
        # Always add empty option at top to work around Bokeh AutocompleteInput click bug
        widget.options = [""] + completions
        # Clear value if it's not in the new completions (empty is always valid)
        if widget.value and widget.value not in completions:
            widget.value = ""
    
    def refresh_contig_options():
        # Skip if locked (during view transitions to avoid cascading updates)
        if global_toggle_lock.get('locked', False):
            return
        # Apply presence filter only if a sample is selected (One Sample view)
        if views.active == 0 and widgets['sample_select'].value:
            sel_sample = widgets['sample_select'].value
            allowed = widgets['sample_to_contigs'].get(sel_sample, set())
            completions = [c for c in orig_contigs if c in allowed]
        else:
            completions = list(orig_contigs)

        # Apply length filter
        if length_slider is not None:
            min_length, max_length = length_slider.value
            completions = [c for c in completions if min_length <= widgets['contig_lengths'].get(c, 0) <= max_length]

        # Apply duplication percentage filter
        if duplication_slider is not None:
            min_dup, max_dup = duplication_slider.value
            dup_min_default = widgets['duplication_percentage_min']
            dup_max_default = widgets['duplication_percentage_max']
            # Only filter if slider is not at default (full range)
            if min_dup > dup_min_default or max_dup < dup_max_default:
                def passes_dup_filter(contig):
                    dup = widgets['contig_duplications'].get(contig)
                    if dup is None:
                        return False  # Exclude contigs without duplication data when filter is active
                    return min_dup <= dup <= max_dup
                completions = [c for c in completions if passes_dup_filter(c)]

        # Apply module-based filters (phage mechanisms, completeness) for the selected sample
        sel_sample = widgets['sample_select'].value if views.active == 0 else None
        module_allowed = get_module_filtered_contigs(sample_name=sel_sample)
        completions = [c for c in completions if c in module_allowed]

        # Apply variable-based filters
        if views.active == 0:
            # "One sample" view: filter contigs based on selected sample's data
            var_allowed = get_variable_filtered_contigs()
            completions = [c for c in completions if c in var_allowed]
        else:
            # "All samples" view: filter contigs where at least one sample has sufficient data
            var_allowed = get_variable_filtered_contigs_any_sample()
            completions = [c for c in completions if c in var_allowed]

        # Apply Filtering2 query builder filters
        filtered_pairs = get_filtering_filtered_pairs()
        if filtered_pairs is not None:
            allowed_contigs = {pair[0] for pair in filtered_pairs}
            completions = [c for c in completions if c in allowed_contigs]

        update_widget_completions(widgets['contig_select'], completions)
        update_section_titles()

    def refresh_sample_options():
        # Skip if locked (during view transitions to avoid cascading updates)
        if global_toggle_lock.get('locked', False):
            return
        # Apply presence filter only if a contig is selected (One Sample view)
        if views.active == 0 and widgets['contig_select'].value:
            sel_contig = widgets['contig_select'].value
            allowed = widgets['contig_to_samples'].get(sel_contig, set())
            completions = [s for s in orig_samples if s in allowed]
        else:
            completions = list(orig_samples)

        # Apply module-based filters (phage mechanisms, completeness) for the selected contig
        sel_contig = widgets['contig_select'].value if views.active == 0 else None
        module_allowed = get_module_filtered_samples(contig_name=sel_contig)
        completions = [s for s in completions if s in module_allowed]

        # Apply variable-based filters (only in "One sample" view)
        if views.active == 0:
            var_allowed = get_variable_filtered_samples()
            completions = [s for s in completions if s in var_allowed]

        # Apply Filtering2 query builder filters
        filtered_pairs = get_filtering_filtered_pairs()
        if filtered_pairs is not None:
            allowed_samples = {pair[1] for pair in filtered_pairs}
            completions = [s for s in completions if s in allowed_samples]

        update_widget_completions(widgets['sample_select'], completions)
        update_section_titles()

    def update_section_titles():
        """Update Filtering, Contigs, and Samples section titles with current counts."""
        filtered_contigs = set(widgets['contig_select'].options) - {""}
        filtered_samples = set(widgets['sample_select'].options) - {""}

        # If a contig is selected, only count pairs for that contig
        selected_contig = widgets['contig_select'].value
        if selected_contig:
            filtered_contigs = {selected_contig}

        # If a sample is selected, only count pairs for that sample
        selected_sample = widgets['sample_select'].value
        if selected_sample:
            filtered_samples = {selected_sample}

        # Count presences (valid contig/sample pairs within filtered sets)
        presences_count = sum(
            1 for contig in filtered_contigs
            for sample in widgets['contig_to_samples'].get(contig, set())
            if sample in filtered_samples
        )

        contigs_count = len(filtered_contigs)
        samples_count = len(filtered_samples)

        filtering_title.text = f"<b>Filtering</b> ({presences_count} contig/sample pairs)"
        contig_title.text = f"<b>Contigs</b> ({contigs_count} available)"
        sample_title.text = f"<b>Samples</b> ({samples_count} available)"

    ## Views function
    # Enforce single-variable selection when in "All samples" view
    # Genome module is in Contigs section and can be selected freely
    def make_global_variable_callback_all(cbg, _unused=None):
        """Callback for All Samples view - enforces single variable selection.

        Only one variable can be selected at a time across all modules in Variables section.
        Genome module is in the Contigs section and is handled separately.
        """
        def callback(attr, old, new):
            if global_toggle_lock['locked']:
                return
            # Only enforce in All samples mode
            if views.active != 1:
                return

            # Determine which index was most-recently changed
            sel_index = None
            old_set = set(old) if old else set()
            new_set = set(new) if new else set()
            added = new_set - old_set

            if added:
                # pick the (one) newly added index
                sel_index = next(iter(added))
            elif new:
                # no clear addition, fall back to last element
                sel_index = new[-1]

            # Enforce single selection across all modules in Variables section
            global_toggle_lock['locked'] = True
            for other in widgets['variables_widgets_all']:
                if other is cbg:
                    if sel_index is None:
                        other.active = []
                    else:
                        other.active = [sel_index]
                else:
                    other.active = []
            global_toggle_lock['locked'] = False
        return callback

    # Views (One sample / All samples) callback: show/hide sample-related controls
    def on_view_change(attr, old, new):
        is_all = (new == 1)  # True means All samples

        # Lock callbacks during view change to prevent cascading updates
        global_toggle_lock['locked'] = True

        # Toggle Sample section - hide entirely in All Samples view
        separator_samples.visible = not is_all
        sample_header.visible = not is_all
        above_sample_content.visible = not is_all
        widgets['sample_select'].visible = not is_all

        # Toggle visibility between the two variables sections
        # Each section maintains its own state independently
        variables_section_one.visible = not is_all
        variables_section_all.visible = is_all

        # Unlock callbacks
        global_toggle_lock['locked'] = False

        # Refresh options after view change
        refresh_contig_options()
        if not is_all:
            refresh_sample_options()

    def create_variable_filter_row():
        """Create a new row of variable filter widgets."""
        type_select = Select(
            options=["#Points", "Max"],
            value="#Points",
            width=70
        )

        var_input = pn.widgets.AutocompleteInput(
            options=widgets['variables'],
            placeholder="Select variable...",
            min_characters=0,
            case_sensitive=False,
            search_strategy="includes",
            sizing_mode="stretch_width"
        )
        
        comparison_select = Select(
            options=[">", "<"],
            value=">",
            width=45
        )
        
        threshold_input = TextInput(
            value="0",
            placeholder="Threshold",
            width=50
        )
        
        plus_btn = Button(label="+", width=30, height=30)
        # len(variable_filter_rows)>0 to make - button invisible initially
        minus_btn = Button(label="−", width=30, height=30, visible=len(variable_filter_rows)>0)
        
        filter_row = row(type_select, var_input.get_root(), comparison_select, threshold_input, \
                         plus_btn, minus_btn, sizing_mode="stretch_width")
        
        # Adjust margins for better spacing
        # (top, right, bottom, left)
        children = filter_row.children
        for i, w in enumerate(children):
            if i == 0:
                # leftmost: keep left margin
                w.margin = (0, 0, 0, 5)
            else:
                # middle widgets: no horizontal margin
                w.margin = (0, 0, 0, 0)
        
        def add_row_callback():
            new_row = create_variable_filter_row()
            variable_filter_rows.append(new_row)
            variable_filters_column.children = list(variable_filters_column.children) + [new_row]
            # Make minus buttons visible when there's more than one row
            for row_widget in variable_filter_rows:
                row_widget.children[-1].visible = True  # Last child is minus button

            # Refresh options after adding a new line
            refresh_on_filter_change()
        
        def remove_row_callback():
            if filter_row in variable_filter_rows:
                variable_filter_rows.remove(filter_row)
                variable_filters_column.children = [r for r in variable_filters_column.children if r != filter_row]
                # Hide minus buttons if only one row remains
                if len(variable_filter_rows) == 1:
                    variable_filter_rows[0].children[-1].visible = False
                    
                # Refresh options after removing filter
                refresh_on_filter_change()
        
        plus_btn.on_click(add_row_callback)
        minus_btn.on_click(remove_row_callback)
        
        # Create a shared callback that refreshes both contig and sample options
        def refresh_on_filter_change(attr, old, new):
            refresh_contig_options()
            refresh_sample_options()
        
        # Attach to all inputs (var_input is Panel widget, others are Bokeh)
        type_select.on_change('value', refresh_on_filter_change)
        var_input.param.watch(lambda event: refresh_on_filter_change(None, None, event.new), 'value')
        comparison_select.on_change('value', refresh_on_filter_change)
        threshold_input.on_change('value', refresh_on_filter_change)

        return filter_row
    
    ## Apply button function
    def apply_clicked():
        try:
            sample = widgets['sample_select'].value
            contig = widgets['contig_select'].value
            is_all = (views.active == 1)

            # Select the correct widget set based on current view
            active_variables_widgets = widgets['variables_widgets_all'] if is_all else widgets['variables_widgets_one']

            # Genome module is shared between views (in Contigs section)
            # Check if gene map should be shown (Gene map is first label in Genome module's cbg)
            genbank_path = db_path if (genome_cbg_one is not None and 0 in genome_cbg_one.active) else None

            # Parse and validate position inputs
            xstart = None
            xend = None
            
            # Validate contig is selected
            if not contig:
                main_placeholder.children = [Div(text="<pre>Error: Please select a contig.</pre>")]
                return
            
            # Get contig length for validation
            contig_length = widgets['contig_lengths'].get(contig, 0)
            
            # Parse position inputs
            try:
                xstart = int(from_position_input.value) if from_position_input.value.strip() else 0
                xend = int(to_position_input.value) if to_position_input.value.strip() else contig_length
            except ValueError:
                main_placeholder.children = [Div(text="<pre>Error: Invalid position range - positions must be integers.</pre>")]
                return
            
            # Validate position range (1-indexed, positions must be within contig bounds)
            if xstart < 0 or xend > contig_length or xstart >= xend:
                main_placeholder.children = [Div(text=f"<pre>Error: Invalid position range - positions must satisfy 0 ≤ start &lt; end ≤ {contig_length}.</pre>")]
                return

            if is_all:
                # All-samples view: require exactly one variable selected from non-Genome modules
                # Genome module is shared (in Contigs section):
                #   - Gene map (index 0) is handled via genbank_path
                #   - Other Genome features (Repeats, etc.) are passed via genome_features parameter
                selected_var = None
                genome_features = []

                # Collect Genome features from shared genome_cbg_one (except Gene map)
                if genome_cbg_one is not None:
                    for idx in genome_cbg_one.active:
                        if idx == 0:  # Skip Gene map (handled via genbank_path)
                            continue
                        genome_features.append(genome_cbg_one.labels[idx])

                # Get the selected variable from non-Genome modules
                for cbg in active_variables_widgets:
                    if cbg.active and selected_var is None:
                        selected_var = cbg.labels[cbg.active[-1]]

                if not selected_var and not genome_features:
                    raise ValueError("When in 'All samples' view you must select at least one variable to plot.")

                # Compute filtered samples (same logic as refresh_sample_options)
                # Start with samples that have the selected contig
                filtered_samples = [s for s in orig_samples if s in widgets['contig_to_samples'].get(contig, set())]
                # Apply module-based filters (coverage, completeness, phage mechanism) for this specific contig
                module_allowed = get_module_filtered_samples(contig_name=contig)
                filtered_samples = [s for s in filtered_samples if s in module_allowed]
                # Apply variable-based filters
                var_allowed = get_variable_filtered_samples()
                filtered_samples = [s for s in filtered_samples if s in var_allowed]

                print(f"[start_bokeh_server] Generating plot for all samples with variable={selected_var}, contig={contig}, genome_features={genome_features}, filtered_samples={len(filtered_samples)}")
                grid = generate_bokeh_plot_all_samples(conn, selected_var, contig, xstart=xstart, xend=xend, genbank_path=genbank_path, genome_features=genome_features if genome_features else None, allowed_samples=set(filtered_samples))
            else:
                # One-sample view: collect possibly-many requested features and call per-sample plot
                requested_features = []

                # Collect features from Variables section
                for cbg in active_variables_widgets:
                    for idx in cbg.active:
                        requested_features.append(cbg.labels[idx])

                # Collect Genome features from shared genome_cbg_one (in Contigs section)
                if genome_cbg_one is not None:
                    for idx in genome_cbg_one.active:
                        requested_features.append(genome_cbg_one.labels[idx])

                print(f"[start_bokeh_server] Generating plot for sample={sample}, contig={contig}, features={requested_features}")
                grid = generate_bokeh_plot_per_sample(conn, requested_features, contig, sample, xstart=xstart, xend=xend, genbank_path=genbank_path)

            main_placeholder.children = [grid]

        except Exception as e:
            tb = traceback.format_exc()
            print(f"[start_bokeh_server] Exception: {tb}", flush=True)
            main_placeholder.children = [Div(text=f"<pre>Error building plot:\n{tb}</pre>")]

    ## Peruse button callback function
    def peruse_clicked():
        """Display summary tables showing coverage, completeness, and phage mechanism data."""
        try:
            contig = widgets['contig_select'].value

            # Check if contig is selected
            if not contig:
                main_placeholder.children = [Div(text="<i>Please select a contig first.</i>")]
                return

            is_all = (views.active == 1)

            if is_all:
                # All Samples view: get filtered samples
                # Start with samples that have the selected contig
                filtered_samples = [s for s in orig_samples if s in widgets['contig_to_samples'].get(contig, set())]
                # Apply module-based filters for this specific contig
                module_allowed = get_module_filtered_samples(contig_name=contig)
                filtered_samples = [s for s in filtered_samples if s in module_allowed]
                # Apply variable-based filters
                var_allowed = get_variable_filtered_samples()
                filtered_samples = [s for s in filtered_samples if s in var_allowed]

                if not filtered_samples:
                    main_placeholder.children = [Div(text="<i>No samples match current filters.</i>")]
                    return

                sample_names = filtered_samples
            else:
                # One Sample view: use selected sample
                sample = widgets['sample_select'].value
                if not sample:
                    main_placeholder.children = [Div(text="<i>Please select a sample first.</i>")]
                    return
                sample_names = [sample]

            # Get contig info
            cur = conn.cursor()
            cur.execute("SELECT Contig_length, Duplication_percentage FROM Contig WHERE Contig_name = ?", (contig,))
            result = cur.fetchone()
            contig_length = result[0] if result else "unknown"
            contig_duplication = result[1] if result else "unknown"

            # Build data
            content = build_summary_data(conn, contig, sample_names)

            if content is None or len(content) == 0:
                main_placeholder.children = [Div(text="<i>No data available for this contig/sample combination.</i>")]
                return

            # Create main header
            header_div = Div(text=f"<h2>{contig} ({contig_length} bp, {contig_duplication} % duplication)</h2><h3>Summary:</h3>")

            # Build content list
            content = [header_div] + content

            # Display
            main_placeholder.children = [column(*content, sizing_mode="stretch_width")]

        except Exception as e:
            tb = traceback.format_exc()
            print(f"[start_bokeh_server] Peruse exception: {tb}", flush=True)
            main_placeholder.children = [Div(text=f"<pre>Error displaying data:\n{tb}</pre>")]

    ### Creating all DOM elements
    # Open DuckDB database connection to build widgets depending on data
    conn = duckdb.connect(db_path, read_only=True)
    widgets = build_controls(conn)

    # Load the CSS and logo
    static_path = os.path.join(os.path.dirname(__file__), "..", "..", "static")
    css_path = os.path.join(static_path, "bokeh_styles.css")
    with open(css_path) as f:
        css_text = f.read()
    stylesheet = InlineStyleSheet(css=css_text)

    # Pink buttons css
    pink_buttons_css_path = os.path.join(static_path, "pink_buttons.css")
    with open(pink_buttons_css_path) as f:
        pink_buttons_css_text = f.read()
    pink_buttons_stylesheet = InlineStyleSheet(css=pink_buttons_css_text)

    # Separate stylesheet for toggle buttons (minimal styling)
    toggle_css_path = os.path.join(static_path, "toggle_styles.css")
    with open(toggle_css_path) as f:
        toggle_css_text = f.read()
    toggle_stylesheet = InlineStyleSheet(css=toggle_css_text)

    # Load logo as base64 to avoid static file serving issues
    logo_path = os.path.join(static_path, "logo.png")
    with open(logo_path, "rb") as f:
        logo_b64 = base64.b64encode(f.read()).decode("utf-8")

    # Create main elements
    ## Views section
    logo = Div(text=f"""<img src="data:image/png;base64,{logo_b64}" style="width:100%; max-width:800px; padding: 0 15%;">""")
    views = RadioButtonGroup(labels=["ONE SAMPLE", "ALL SAMPLES"], active=0, sizing_mode="stretch_width", stylesheets=[stylesheet])

    # Global lock for toggles when enforcing "All samples" view (single-variable mode)
    global_toggle_lock = {'locked': False}

    # Attach global variable callbacks to All Samples CheckboxButtonGroups only
    # (One Sample view allows multiple selections, so no callback needed there)
    # Genome module is in Contigs section now, so pass None as genome_cbg_ref
    for cbg in widgets['variables_widgets_all']:
        cbg.on_change('active', make_global_variable_callback_all(cbg, None))
    views.on_change('active', on_view_change)


    ## Build Filtering section (dynamic query builder with AND/OR logic)
    filtering_toggle_btn = Button(label="▼", width=20, height=20, button_type="primary", align="center", margin=0, stylesheets=[toggle_stylesheet])
    filtering_title = Div(text="<b>Filtering</b>", align="center")
    filtering_header = row(filtering_toggle_btn, filtering_title, sizing_mode="stretch_width", align="center")

    # Cache filtering metadata once when document loads
    filtering_metadata = get_filtering_metadata(db_path)

    # Store all OR sections in a list for dynamic management
    or_sections = []

    def count_total_query_rows():
        """Count total number of query rows across all OR sections."""
        total = 0
        for section_data in or_sections:
            total += len(section_data['rows'])
        return total

    def create_query_row(section_data):
        """Create a single query row with cascading selects, comparison, dynamic input and remove button."""
        # Get categories from metadata
        categories = list(filtering_metadata.keys())
        if not categories:
            categories = ["No data"]

        initial_category = categories[0]
        initial_columns = list(filtering_metadata.get(initial_category, {}).get('columns', {}).keys())
        if not initial_columns:
            initial_columns = ["No columns"]
        initial_column = initial_columns[0]

        # Determine initial column type
        initial_col_info = filtering_metadata.get(initial_category, {}).get('columns', {}).get(initial_column, {})
        initial_is_text = initial_col_info.get('type') == 'text'

        # First level select (categories)
        category_select = Select(
            options=categories,
            value=initial_category,
            width=110,
            margin=(0, 2, 0, 0)
        )

        # Second level select (columns)
        subcategory_select = Select(
            options=initial_columns,
            value=initial_column,
            sizing_mode="stretch_width",
            margin=(0, 2, 0, 0)
        )

        # Comparison operator select - "=" and "!=" for text, all operators for numeric
        comparison_select = Select(
            options=["=", "!="] if initial_is_text else ["=", ">", "<", "!="],
            value="=",
            margin=(0, 2, 0, 0)
        )

        # Container for the dynamic input widget
        input_container = column(sizing_mode="stretch_width", margin=(0, 2, 0, 0))

        def refresh_on_filter_change():
            """Refresh contig and sample options when Filtering2 values change."""
            refresh_contig_options()
            refresh_sample_options()

        # Create initial input widget based on column type
        if initial_is_text:
            distinct_values = initial_col_info.get('distinct_values', [])
            initial_input = pn.widgets.AutocompleteInput(
                value="",
                options=distinct_values,
                min_characters=0,
                case_sensitive=False,
                placeholder="Value...",
                search_strategy="includes",
                sizing_mode="stretch_width",
                margin=(0, 2, 0, 0)
            )
            input_container.children = [initial_input.get_root()]
            # Add callback for Panel AutocompleteInput
            initial_input.param.watch(lambda event: refresh_on_filter_change(), 'value')
        else:
            initial_input = Spinner(value=0, placeholder="Value...", width=80, margin=(0, 2, 0, 0))
            input_container.children = [initial_input]
            # Add callback for Bokeh Spinner
            initial_input.on_change('value', lambda attr, old, new: refresh_on_filter_change())

        # Remove button
        minus_btn = Button(label="−", width=30, height=30, stylesheets=[stylesheet], margin=(0, 10, 0, 0))

        # Store reference to current input widget (for later retrieval)
        current_input_ref = {'widget': initial_input, 'is_panel': initial_is_text}

        def update_input_widget(col_name):
            """Update the input widget based on column type."""
            category = category_select.value
            col_info = filtering_metadata.get(category, {}).get('columns', {}).get(col_name, {})
            is_text = col_info.get('type') == 'text'

            # Update comparison options based on type
            if is_text:
                comparison_select.options = ["=", "!="]
                if comparison_select.value not in comparison_select.options:
                    comparison_select.value = "="
            else:
                comparison_select.options = ["=", ">", "<", "!="]
                if comparison_select.value not in comparison_select.options:
                    comparison_select.value = "="

            # Create new input widget
            if is_text:
                distinct_values = col_info.get('distinct_values', [])
                new_input = pn.widgets.AutocompleteInput(
                    value="",
                    options=distinct_values,
                    min_characters=0,
                    case_sensitive=False,
                    placeholder="Value...",
                    search_strategy="includes",
                    sizing_mode="stretch_width",
                    margin=(0, 2, 0, 0)
                )
                input_container.children = [new_input.get_root()]
                current_input_ref['widget'] = new_input
                current_input_ref['is_panel'] = True
                # Add callback for Panel AutocompleteInput
                new_input.param.watch(lambda event: refresh_on_filter_change(), 'value')
            else:
                new_input = Spinner(value=0, placeholder="Value...", width=80, margin=(0, 2, 0, 0))
                input_container.children = [new_input]
                current_input_ref['widget'] = new_input
                current_input_ref['is_panel'] = False
                # Add callback for Bokeh Spinner
                new_input.on_change('value', lambda attr, old, new: refresh_on_filter_change())

        def update_subcategories(attr, old, new):
            """Update column options when category changes."""
            columns = list(filtering_metadata.get(new, {}).get('columns', {}).keys())
            if not columns:
                columns = ["No columns"]
            subcategory_select.options = columns
            subcategory_select.value = columns[0]
            # Update input widget for new column
            update_input_widget(columns[0])

        def update_input_on_column_change(attr, old, new):
            """Update input widget when column changes."""
            update_input_widget(new)

        category_select.on_change('value', update_subcategories)
        subcategory_select.on_change('value', update_input_on_column_change)
        # Add callback to comparison_select to refresh when operator changes
        comparison_select.on_change('value', lambda attr, old, new: refresh_on_filter_change())

        query_row = row(category_select, subcategory_select, comparison_select, input_container, minus_btn,
                       sizing_mode="stretch_width", margin=(5, 0, 5, 0))

        # Store reference to this row
        row_data = {
            'query_row': query_row,
            'category_select': category_select,
            'subcategory_select': subcategory_select,
            'comparison_select': comparison_select,
            'input_ref': current_input_ref,
            'and_div': None  # Will be set when AND is added above this row
        }

        def remove_row_callback():
            # Don't allow removal if this is the only query row across all sections
            if count_total_query_rows() <= 1:
                return

            # Find which section this row belongs to
            for sec_data in or_sections:
                if row_data in sec_data['rows']:
                    # Remove the row
                    sec_data['rows'].remove(row_data)

                    # Rebuild the section
                    rebuild_section(sec_data)

                    # If section is now empty, remove it
                    if len(sec_data['rows']) == 0:
                        or_sections.remove(sec_data)
                        rebuild_filtering_content()

                    # Refresh filter options after row removal
                    refresh_on_filter_change()
                    break

        minus_btn.on_click(remove_row_callback)
        return row_data

    def rebuild_section(section_data):
        """Rebuild a section's content with all its query rows and the Add AND button."""
        section_children = []

        for i, row_data in enumerate(section_data['rows']):
            # Add AND div before each row except the first
            if i > 0:
                select_widget = Select(
                    options=["AND", "OR"],
                    value="AND",
                    margin=(5, 0, 5, 0)
                )
                # Add callback to refresh when AND/OR changes
                select_widget.on_change('value', lambda attr, old, new: (refresh_contig_options(), refresh_sample_options()))
                section_children.append(select_widget)
                row_data['and_div'] = select_widget
            else:
                row_data['and_div'] = None

            section_children.append(row_data['query_row'])

        # Add the "+ Add AND" button
        section_children.append(section_data['add_and_btn'])

        # Update the section column's children
        section_data['column'].children = section_children

    def create_or_section():
        """Create a new OR section with one query row and Add AND/OR button."""
        section_column = column(
            sizing_mode="stretch_width",
            styles={'border-left': '3px solid #00b17c', 'padding-left': '10px', 'margin-left': '5px'}
        )

        add_and_btn = Button(
            label="+ Add AND/OR",
            margin=(5, 0, 5, 0),
            stylesheets=[stylesheet]
        )

        section_data = {
            'column': section_column,
            'rows': [],
            'add_and_btn': add_and_btn
        }

        def add_and_or_callback():
            # Create a new query row
            new_row = create_query_row(section_data)
            section_data['rows'].append(new_row)
            rebuild_section(section_data)

        add_and_btn.on_click(add_and_or_callback)

        # Create initial query row
        initial_row = create_query_row(section_data)
        section_data['rows'].append(initial_row)

        rebuild_section(section_data)

        return section_data

    # Create the global "+ Add AND/OR" button
    global_add_btn = Button(
        label="+ Add AND/OR",
        margin=(10, 0, 5, 0),
        stylesheets=[pink_buttons_stylesheet]
    )

    # Store reference to the button widget for replacement
    global_widget_state = {'widget': global_add_btn}

    def rebuild_filtering_content():
        """Rebuild the entire Filtering content with all OR sections."""
        content_children = []

        for i, section_data in enumerate(or_sections):
            # Add OR div before each section except the first
            if i > 0:
                select_widget = Select(
                    options=["AND", "OR"],
                    value="AND",
                    margin=(5, 0, 5, 0)
                )
                # Add callback to refresh when AND/OR changes
                select_widget.on_change('value', lambda attr, old, new: (refresh_contig_options(), refresh_sample_options()))
                content_children.append(select_widget)

            content_children.append(section_data['column'])

        # Add the current global widget (button or select) at the end
        content_children.append(global_widget_state['widget'])

        filtering_content.children = content_children

    def global_add_and_or_callback():
        # Create a new section
        new_section = create_or_section()
        or_sections.append(new_section)
        rebuild_filtering_content()

    global_add_btn.on_click(global_add_and_or_callback)

    # Create initial OR section
    initial_section = create_or_section()
    or_sections.append(initial_section)

    # Create the main content container (using Bokeh column, not Panel)
    filtering_content = column(
        sizing_mode="stretch_width"
    )

    # Initial build of content
    rebuild_filtering_content()

    # Add toggle callback for collapsible Filtering section
    filtering_toggle_btn.on_click(make_toggle_callback(filtering_toggle_btn, filtering_content))


    # Initialize slider variables to None (these were used by the old Filtering section
    # but are still referenced by get_module_filtered_contigs/samples functions)
    phage_mechanism_filter = None
    prevalence_left_slider = None
    prevalence_right_slider = None
    circularising_slider = None
    pct_completeness_slider = None
    pct_contamination_slider = None
    coverage_percentage_slider = None
    coverage_mean_slider = None
    coverage_median_slider = None
    coverage_variation_slider = None
    coverage_sd_slider = None
    length_slider = None
    duplication_slider = None
    correction_filter = None
    variable_filter_rows = []


    ## Build Sample section
    sample_toggle_btn = Button(label="▼", width=20, height=20, button_type="primary", align="center", margin=0, stylesheets=[toggle_stylesheet])
    sample_title = Div(text="<b>Samples</b>", align="center")
    sample_header = row(sample_toggle_btn, sample_title, sizing_mode="stretch_width", align="center")

    above_sample_children = []
    above_sample_content = column(
        *above_sample_children,
        visible=True, sizing_mode="stretch_width"
    )

    sample_toggle_btn.on_click(make_toggle_callback(sample_toggle_btn, above_sample_content))
    widgets['sample_select'].param.watch(lambda event: refresh_contig_options(), 'value')


    ## Build Contig section
    # Keep original full lists so we can restore when filters are off
    orig_contigs = list(widgets['contigs'])
    orig_samples = list(widgets['samples'])

    # Create Contigs section header
    contig_toggle_btn = Button(label="▼", width=20, height=20, button_type="primary", align="center", margin=0, stylesheets=[toggle_stylesheet])
    contig_toggle_btn.styles = {'padding': '0px', 'line-height': '20px'}
    contig_title = Div(text="<b>Contigs</b>", align="center")
    contig_header = row(contig_toggle_btn, contig_title, sizing_mode="stretch_width", align="center")
    above_contig_children = []
    
    above_contig_content = column(
        *above_contig_children,
        visible=True, sizing_mode="stretch_width"
    )
    contig_toggle_btn.on_click(make_toggle_callback(contig_toggle_btn, above_contig_content))

    def on_contig_change(event):
        new = event.new
        refresh_sample_options()
        # Update position inputs when contig changes
        if new and new in widgets['contig_lengths']:
            from_position_input.value = "0"
            to_position_input.value = str(widgets['contig_lengths'][new])
        else:
            from_position_input.value = "0"
            to_position_input.value = ""
    
    widgets['contig_select'].param.watch(on_contig_change, 'value')


    ## Build Variables section - TWO SEPARATE SECTIONS for each view
    variables_title_one = Div(text="<b>Variables</b>")
    variables_title_all = Div(text="<b>Variables</b>")
    genome_cbg_one = None  # Will store reference to Genome module's CheckboxButtonGroup (shared between views)

    # Build "One Sample" view variables section
    # Has module checkboxes + collapsible variable groups
    # NOTE: Genome module is built separately and placed in Contigs section
    controls_variables_one = []
    module_toggles_one = []
    module_contents_one = []
    genome_index_one = None  # Track Genome module index for separate handling

    for i, module_widget in enumerate(widgets['module_widgets_one']):
        module_name = widgets['module_names'][i]
        help_tooltip = widgets['helps_widgets'][i]

        # Skip Genome module - will be added to Contigs section
        if module_name == "Genome":
            genome_index_one = i
            genome_cbg_one = widgets['variables_widgets_one'][i]
            continue

        # Create toggle button for collapsible section
        toggle_btn = Button(label="▶", width=20, height=20, button_type="primary", align="center", margin=0, stylesheets=[toggle_stylesheet])
        toggle_btn.styles = {'padding': '0px', 'line-height': '20px'}
        module_toggles_one.append(toggle_btn)

        # Build header with checkbox (module_widget has module name as label)
        if help_tooltip is not None:
            help_btn = HelpButton(tooltip=help_tooltip, width=20, height=20, align="center", button_type="light", stylesheets=[toggle_stylesheet])
            hdr = row(toggle_btn, module_widget, help_btn, sizing_mode="stretch_width", align="center")
        else:
            hdr = row(toggle_btn, module_widget, sizing_mode="stretch_width", align="center")

        controls_variables_one.append(hdr)

        # Add the module's CheckboxButtonGroup for variables (collapsible, starts folded)
        cbg = widgets['variables_widgets_one'][i]
        cbg.visible = False
        module_contents_one.append(cbg)
        controls_variables_one.append(cbg)

    # Add toggle callbacks for One Sample view
    for i, toggle_btn in enumerate(module_toggles_one):
        content = module_contents_one[i]
        toggle_btn.on_click(make_toggle_callback(toggle_btn, content))

    # Build "All Samples" view variables section
    # NOTE: Genome module is built separately and placed in Contigs section
    # Other modules have title-only headers (no checkbox)
    controls_variables_all = []
    module_toggles_all = []
    module_contents_all = []

    for i in range(len(widgets['module_names'])):
        module_name = widgets['module_names'][i]
        help_tooltip = widgets['helps_widgets'][i]

        # Skip Genome module - it's in the Contigs section (shared between views)
        if module_name == "Genome":
            continue

        # Create toggle button for collapsible section
        toggle_btn = Button(label="▶", width=20, height=20, button_type="primary", align="center", margin=0, stylesheets=[toggle_stylesheet])
        toggle_btn.styles = {'padding': '0px', 'line-height': '20px'}
        module_toggles_all.append(toggle_btn)

        # Build header with title only (no checkbox) for non-Genome modules
        module_title_div = Div(text=f"{module_name}", align="center")
        if help_tooltip is not None:
            # Create new tooltip instance to avoid "already in doc" error
            help_text = help_tooltip.content
            tooltip_all = Tooltip(content=help_text, position="right")
            help_btn = HelpButton(tooltip=tooltip_all, width=20, height=20, align="center", button_type="light", stylesheets=[toggle_stylesheet])
            hdr = row(toggle_btn, module_title_div, help_btn, sizing_mode="stretch_width", align="center")
        else:
            hdr = row(toggle_btn, module_title_div, sizing_mode="stretch_width", align="center")

        controls_variables_all.append(hdr)

        # Add the module's CheckboxButtonGroup for variables (collapsible, starts folded)
        cbg = widgets['variables_widgets_all'][i]
        cbg.visible = False
        module_contents_all.append(cbg)
        controls_variables_all.append(cbg)

    # Add toggle callbacks for All Samples view
    for i, toggle_btn in enumerate(module_toggles_all):
        content = module_contents_all[i]
        toggle_btn.on_click(make_toggle_callback(toggle_btn, content))

    # Create the two variables section containers
    variables_section_one = column(variables_title_one, *controls_variables_one, visible=True, sizing_mode="stretch_width")
    variables_section_all = column(variables_title_all, *controls_variables_all, visible=False, sizing_mode="stretch_width")


    ## Build Genome module controls (placed in Contigs section, shared between views)
    genome_section = None

    if genome_index_one is not None:
        # Get Genome module info
        genome_help_tooltip = widgets['helps_widgets'][genome_index_one]

        # Build simple Genome section header (no collapse needed - Contigs section handles that)
        genome_title = Div(text="<b>Plot genomic features:</b>", align="center")
        
        if genome_help_tooltip is not None:
            help_btn = HelpButton(tooltip=genome_help_tooltip, width=20, height=20, align="center", button_type="light", stylesheets=[toggle_stylesheet])
            genome_hdr = row(genome_title, help_btn, sizing_mode="stretch_width", align="center")
        else:
            genome_hdr = genome_title

        genome_cbg_one.visible = True
        genome_section = column(genome_hdr, genome_cbg_one, visible=True, sizing_mode="stretch_width")

    # Add Genome section to contig_content
    below_contig_children = []
    
    # Create position range inputs
    from_position_input = TextInput(value="0", placeholder="Start position", sizing_mode="stretch_width", margin=(0, 0, 0, 0))
    to_position_input = TextInput(value="", placeholder="End position", sizing_mode="stretch_width", margin=(0, 10, 0, 0))
    
    position_label_from = Div(text="From", width=40, margin=(5, 5, 5, 0))
    position_label_to = Div(text="to", width=25, margin=(5, 5, 5, 5))
    
    position_row = row(
        position_label_from, from_position_input, 
        position_label_to, to_position_input,
        sizing_mode="stretch_width"
    )
    
    below_contig_children.append(position_row)
    
    if genome_section is not None:
        below_contig_children = list(below_contig_children) + [genome_section]

    below_contig_content = column(
        *below_contig_children,
        visible=True, sizing_mode="stretch_width"
    )
    contig_toggle_btn.on_click(make_toggle_callback(contig_toggle_btn, below_contig_content))

    # Initialize position inputs if contig is pre-filled
    if widgets['contig_select'].value and widgets['contig_select'].value in widgets['contig_lengths']:
        to_position_input.value = str(widgets['contig_lengths'][widgets['contig_select'].value])


    ### Attach callbacks for One Sample view (module checkbox ↔ variable bidirectional sync)
    for i, mc in enumerate(widgets['module_widgets_one']):
        toggles = widgets['variables_widgets_one'][i]
        lock = {"locked": False}  # per-module lock

        # Module → toggles (CheckboxButtonGroup)
        def make_module_callback(mc, toggles, lock):
            def callback(attr, old, new):
                if lock.get("locked", False) or global_toggle_lock.get("locked", False):
                    return
                lock["locked"] = True
                module_on = 0 in mc.active
                if module_on:
                    toggles.active = list(range(len(toggles.labels)))
                else:
                    toggles.active = []
                lock["locked"] = False
            return callback

        mc.on_change("active", make_module_callback(mc, toggles, lock))

        # Variable → module (update module checkbox only)
        def make_variable_callback(mc, toggles, lock):
            def callback(attr, old, new):
                if lock.get("locked", False) or global_toggle_lock.get("locked", False):
                    return
                total = len(toggles.labels)
                active_count = len(toggles.active)

                lock["locked"] = True
                if active_count == total and total > 0:
                    mc.active = [0]
                else:
                    mc.active = []
                lock["locked"] = False
            return callback

        toggles.on_change("active", make_variable_callback(mc, toggles, lock))

    ## Create final Apply and Peruse data buttons
    apply_button = Button(label="APPLY", align="center", sizing_mode="stretch_width", stylesheets=[stylesheet], css_classes=["apply-btn"])
    apply_button.on_click(lambda: apply_clicked())

    peruse_button = Button(label="PERUSE DATA", align="center", sizing_mode="stretch_width", stylesheets=[stylesheet], css_classes=["apply-btn"])
    peruse_button.on_click(lambda: peruse_clicked())

    buttons_row = row(peruse_button, apply_button, sizing_mode="stretch_width", spacing=10)


    ## Initialize section titles with counts
    update_section_titles()

    ## Put together all DOM elements
    # Create visual separators (horizontal lines) using background color
    separator_filtering = Div(text="", height=1, sizing_mode="stretch_width", styles={'background-color': '#333', 'margin': '10px 0'})
    separator_samples = Div(text="", height=1, sizing_mode="stretch_width", styles={'background-color': '#333', 'margin': '10px 0'})
    separator_contigs = Div(text="", height=1, sizing_mode="stretch_width", styles={'background-color': '#333', 'margin': '10px 0'})
    separator_variables = Div(text="", height=1, sizing_mode="stretch_width", styles={'background-color': '#333', 'margin': '10px 0'})

    # Gene map is now part of the Genome module's CheckboxButtonGroup
    # Include both variables sections - visibility is toggled by on_view_change
    controls_children = [logo, views, separator_filtering, filtering_header, filtering_content,
                         separator_contigs, contig_header, above_contig_content, widgets['contig_select'], below_contig_content,
                         separator_samples, sample_header, above_sample_content, widgets['sample_select'],
                         separator_variables,
                         variables_section_one,  # One Sample view (with module checkboxes)
                         variables_section_all,  # All Samples view (title headers only)
                         buttons_row]

    controls_column = pn.Column(*controls_children, sizing_mode="stretch_height", css_classes=["left-col"])

    main_placeholder = column(Div(text="<i>No plot yet. Select one sample, one contig and at least one variable in \"One sample\" mode or one contig and one variable in \"All samples\" mode and click Apply.</i>"), sizing_mode="stretch_both")

    # Wrap everything in a Flex container
    layout = pn.Row(controls_column, main_placeholder, sizing_mode="stretch_both")
    layout.stylesheets = [stylesheet]

    return layout

def add_serve_args(parser):
    parser.add_argument("--db", required=True, help="Path to DuckDB database")
    parser.add_argument("--port", type=int, default=5006, help="Port to serve Panel app")

def run_serve(args):
    # Create a factory function that Panel will call for each session
    def create_app():
        return create_layout(args.db)

    static_path = os.path.join(os.path.dirname(__file__), "..", "..", "static")
    pn.serve(
        create_app,
        port=args.port,
        show=True,
        title="theBIGbam",
        static_dirs={'assets': static_path}
    )
    return 0

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    add_serve_args(parser)
    args = parser.parse_args()
    raise SystemExit(run_serve(args))