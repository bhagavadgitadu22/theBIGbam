import argparse
import base64
import os
import duckdb
import traceback

import panel as pn

from bokeh.layouts import column, row
from bokeh.models import Div, InlineStyleSheet, Tooltip
from bokeh.models.widgets import AutocompleteInput, CheckboxGroup, HelpButton, Button, RadioButtonGroup, CheckboxButtonGroup, Select, TextInput
from bokeh.models.plots import GridPlot

# Import the plotting function from the repo
from .plotting_data_per_sample import generate_bokeh_plot_per_sample
from .plotting_data_all_samples import generate_bokeh_plot_all_samples
from .perusing_data import build_summary_data, generate_summary_table

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

    # Widget Selector for Contigs (autocomplete with max 20 suggestions)
    cur.execute("SELECT Contig_name, Contig_length, Duplication_percentage FROM Contig ORDER BY Contig_name")
    rows = cur.fetchall()
    contigs = [r[0] for r in rows]
    contig_lengths = {r[0]: r[1] for r in rows}  # Dictionary mapping contig_name -> length
    contig_duplications = {r[0]: r[2] for r in rows}  # Dictionary mapping contig_name -> duplication percentage (can be None)
    
    # If only one contig in database, pre-fill the field
    contig_select = AutocompleteInput(value=contigs[0] if len(contigs) == 1 else "",
                                      completions=[""] + contigs,
                                      min_characters=0,
                                      case_sensitive=False,
                                      restrict=False,
                                      max_completions=20,
                                      placeholder="Type to search contigs...",
                                      sizing_mode="stretch_width")

    # Widget Selector for Samples (autocomplete with max 20 suggestions)
    cur.execute("SELECT Sample_name FROM Sample ORDER BY Sample_name")
    samples = [r[0] for r in cur.fetchall()]
    
    # If only one sample in database, pre-fill the field
    sample_select = AutocompleteInput(value=samples[0] if len(samples) == 1 else "",
                                      completions=[""] + samples,
                                      min_characters=0,
                                      case_sensitive=False,
                                      restrict=False,
                                      max_completions=20,
                                      placeholder="Type to search samples...",
                                      sizing_mode="stretch_width")

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
    MODULE_ORDER = ["Genome", "Coverage", "Mapping metrics per position", "Long-read metrics", "Paired-read metrics", "Phage termini"]

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
        'has_duplication_data': has_duplication_data
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

    ## Helper function to create collapsible section toggle callbacks
    def make_toggle_callback(btn, content):
        def callback():
            content.visible = not content.visible
            if content.visible:
                btn.label = "▼"
            else:
                btn.label = "▶"
        return callback
    
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
    def get_module_filtered_contigs():
        """Apply module-based filters (coverage, completeness, phage mechanism) to get allowed contigs."""
        allowed_contigs = set(orig_contigs)
        cur = conn.cursor()

        # Apply coverage filters using Explicit_presences view (single combined query)
        conditions = []
        params = []

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

        if conditions:
            query = f"SELECT DISTINCT Contig_name FROM Explicit_completeness WHERE {' AND '.join(conditions)}"
            cur.execute(query, params)
            matching = {row[0] for row in cur.fetchall()}
            allowed_contigs &= matching

        # Apply phage mechanism filter using Explicit_phage_mechanisms view
        if phage_mechanism_filter is not None and phage_mechanism_filter.value:
            query = """
                SELECT DISTINCT Contig_name
                FROM Explicit_phage_mechanisms
                WHERE Phage_packaging_mechanism = ?
            """
            cur.execute(query, (phage_mechanism_filter.value,))
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

        return allowed_samples

    def update_widget_completions(widget, completions):
        """Update widget completions. Clear value if not in completions."""
        # Always add empty option at top to work around Bokeh AutocompleteInput click bug
        widget.completions = [""] + completions
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

        # Apply module-based filters (phage mechanisms, completeness)
        module_allowed = get_module_filtered_contigs()
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

        # Apply module-based filters (phage mechanisms, completeness)
        module_allowed = get_module_filtered_samples()
        completions = [s for s in completions if s in module_allowed]

        # Apply variable-based filters (only in "One sample" view)
        if views.active == 0:
            var_allowed = get_variable_filtered_samples()
            completions = [s for s in completions if s in var_allowed]

        update_widget_completions(widgets['sample_select'], completions)
        update_section_titles()

    def update_section_titles():
        """Update Filtering, Contigs, and Samples section titles with current counts."""
        filtered_contigs = set(widgets['contig_select'].completions) - {""}
        filtered_samples = set(widgets['sample_select'].completions) - {""}

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

        var_input = AutocompleteInput(
            completions=widgets['variables'],
            placeholder="Select variable...",
            min_characters=0,
            case_sensitive=False,
            restrict=False,
            max_completions=10,
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
        
        filter_row = row(type_select, var_input, comparison_select, threshold_input, \
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
        
        # Attach to all three inputs
        type_select.on_change('value', refresh_on_filter_change)
        var_input.on_change('value', refresh_on_filter_change)
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

            # Save current positions of the plot for restoration after re-plot
            fig = main_placeholder.children[0]

            xstart = None
            xend = None
            if isinstance(fig, GridPlot):
                subplot = fig.children[0][0]
                xstart = subplot.x_range.start
                xend = subplot.x_range.end

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
    logo = Div(text=f"""<img src="data:image/png;base64,{logo_b64}" style="width:100%; max-width:600px; padding: 0 20%;">""")
    views = RadioButtonGroup(labels=["ONE SAMPLE", "ALL SAMPLES"], active=0, sizing_mode="stretch_width", stylesheets=[stylesheet])

    # Global lock for toggles when enforcing "All samples" view (single-variable mode)
    global_toggle_lock = {'locked': False}

    # Attach global variable callbacks to All Samples CheckboxButtonGroups only
    # (One Sample view allows multiple selections, so no callback needed there)
    # Genome module is in Contigs section now, so pass None as genome_cbg_ref
    for cbg in widgets['variables_widgets_all']:
        cbg.on_change('active', make_global_variable_callback_all(cbg, None))
    views.on_change('active', on_view_change)


    ## Build filtering section
    filtering_title = Div(text="<b>Filtering</b>")
    filtering_header = filtering_title  # No toggle button - always visible

    filtering_children = []

    # Initialize module filter variables (will be set below if applicable)
    phage_mechanism_filter = None
    prevalence_left_slider = None
    prevalence_right_slider = None
    pct_completeness_slider = None
    pct_contamination_slider = None
    coverage_percentage_slider = None
    coverage_mean_slider = None
    coverage_median_slider = None
    coverage_variation_slider = None
    coverage_sd_slider = None

    # Add "Contig filters" collapsible subsection (if multiple contigs)
    min_len = min(widgets['contig_lengths'].values())
    max_len = max(widgets['contig_lengths'].values())
    length_slider = None
    duplication_slider = None
    if len(widgets["contigs"]) > 1:
        contig_toggle_btn = Button(label="▶", width=20, height=20, button_type="primary", align="center", margin=0, stylesheets=[toggle_stylesheet])
        contig_filters_title = Div(text="Contig filters", align="center")
        contig_filters_header = row(contig_toggle_btn, contig_filters_title, sizing_mode="stretch_width", align="center")

        contig_filters_children = []

        length_slider = create_editable_range_slider(
            "Contig length", min_len, max_len, (min_len, max_len), step=1,
            on_change=lambda event: refresh_contig_options()
        )
        contig_filters_children.append(length_slider)

        # Add duplication percentage slider if data exists
        if widgets['has_duplication_data']:
            dup_min = widgets['duplication_percentage_min']
            dup_max = widgets['duplication_percentage_max']
            duplication_slider = create_editable_range_slider(
                "Duplication (%)", dup_min, dup_max, (dup_min, dup_max), step=1,
                on_change=lambda event: refresh_contig_options()
            )
            contig_filters_children.append(duplication_slider)

        contig_filters_content = pn.Column(*contig_filters_children, visible=False, sizing_mode="stretch_width")
        contig_toggle_btn.on_click(make_toggle_callback(contig_toggle_btn, contig_filters_content))

        filtering_children.append(contig_filters_header)
        filtering_children.append(contig_filters_content)

    # Add "Coverage filters" collapsible subsection (always available if any module filters exist)
    has_module_filters = widgets['phage_mechanisms'] or widgets['has_completeness']
    if has_module_filters:
        # Create collapsible Coverage filters section
        coverage_toggle_btn = Button(label="▶", width=20, height=20, button_type="primary", align="center", margin=0, stylesheets=[toggle_stylesheet])
        coverage_filters_title = Div(text="Coverage filters", align="center")
        coverage_filters_header = row(coverage_toggle_btn, coverage_filters_title, sizing_mode="stretch_width", align="center")

        coverage_filters_children = []

        correction_label = Div(text="Harmonise coverage mean/median between samples:", margin=(5, 0, 0, 10))
        correction_filter = Select(
                options=["No correction", "Correct by number of mapped reads", "Correct by number of reads"],
                value="No correction",   # default = no filter (empty = all)
                sizing_mode="stretch_width",
            )
        coverage_filters_children.append(correction_label)
        coverage_filters_children.append(correction_filter)

        # Coverage patterns sliders (always available from Presences table)
        cov_mean_min = 0
        cov_mean_max = widgets['coverage_mean_max']
        coverage_mean_slider = create_editable_range_slider(
            "Coverage mean", cov_mean_min, cov_mean_max, (cov_mean_min, cov_mean_max), step=1,
            on_change=lambda event: (refresh_contig_options(), refresh_sample_options())
        )
        coverage_filters_children.append(coverage_mean_slider)

        # Function to update coverage_mean_slider range when correction method changes
        def update_coverage_mean_slider():
            cov_col = get_coverage_mean_column()
            cur = conn.cursor()
            cur.execute(f"SELECT MIN({cov_col}), MAX({cov_col}) FROM Explicit_presences WHERE {cov_col} IS NOT NULL")
            result = cur.fetchone()
            if result and result[0] is not None and result[1] is not None:
                new_max = result[1]
                coverage_mean_slider.start = 0
                coverage_mean_slider.end = new_max
                coverage_mean_slider.value = (0, new_max)
            refresh_contig_options()
            refresh_sample_options()

        # Coverage median slider
        cov_median_min = 0
        cov_median_max = widgets['coverage_median_max']
        coverage_median_slider = create_editable_range_slider(
            "Coverage median", cov_median_min, cov_median_max, (cov_median_min, cov_median_max), step=1,
            on_change=lambda event: (refresh_contig_options(), refresh_sample_options())
        )
        coverage_filters_children.append(coverage_median_slider)

        # Function to update coverage_median_slider range when correction method changes
        def update_coverage_median_slider():
            cov_col = get_coverage_median_column()
            cur = conn.cursor()
            cur.execute(f"SELECT MIN({cov_col}), MAX({cov_col}) FROM Explicit_presences WHERE {cov_col} IS NOT NULL")
            result = cur.fetchone()
            if result and result[0] is not None and result[1] is not None:
                new_max = result[1]
                coverage_median_slider.start = 0
                coverage_median_slider.end = new_max
                coverage_median_slider.value = (0, new_max)
            refresh_contig_options()
            refresh_sample_options()

        # Wire up correction_filter to update sliders when changed
        if hasattr(correction_filter, 'on_change'):
            correction_filter.on_change('value', lambda attr, old, new: (update_coverage_mean_slider(), update_coverage_median_slider()))

        coverage_percentage_slider = create_editable_range_slider(
            "Aligned fraction (%)", 0, 100, (0, 100), step=1,
            on_change=lambda event: (refresh_contig_options(), refresh_sample_options())
        )
        coverage_filters_children.append(coverage_percentage_slider)

        cov_var_min = 0.0
        cov_var_max = widgets['coverage_variation_max']
        # Slider shows normalized variance (stored value / 1000000)
        coverage_variation_slider = create_editable_range_slider(
            "Coverage variation", cov_var_min, cov_var_max, (cov_var_min, cov_var_max), step=0.01,
            on_change=lambda event: (refresh_contig_options(), refresh_sample_options())
        )
        coverage_filters_children.append(coverage_variation_slider)

        cov_sd_min = 0.0
        cov_sd_max = widgets['coverage_sd_max']
        # Slider shows coefficient of variation (stored value / 1000000)
        coverage_sd_slider = create_editable_range_slider(
            "Coverage SD", cov_sd_min, cov_sd_max, (cov_sd_min, cov_sd_max), step=0.01,
            on_change=lambda event: (refresh_contig_options(), refresh_sample_options())
        )
        coverage_filters_children.append(coverage_sd_slider)

        # Create coverage filters content container (collapsible)
        coverage_filters_content = pn.Column(*coverage_filters_children, visible=False, sizing_mode="stretch_width")
        coverage_toggle_btn.on_click(make_toggle_callback(coverage_toggle_btn, coverage_filters_content))

        filtering_children.append(coverage_filters_header)
        filtering_children.append(coverage_filters_content)

    # Add "Completeness filters" collapsible subsection (if Completeness table has data)
    if widgets['has_completeness']:
        # Create collapsible Completeness filters section
        completeness_toggle_btn = Button(label="▶", width=20, height=20, button_type="primary", align="center", margin=0, stylesheets=[toggle_stylesheet])
        completeness_filters_title = Div(text="Completeness filters", align="center")
        completeness_filters_header = row(completeness_toggle_btn, completeness_filters_title, sizing_mode="stretch_width", align="center")

        completeness_filters_children = []

        pct_completeness_slider = create_editable_range_slider(
            "Whole completeness (%)", 0, 100, (0, 100), step=1,
            on_change=lambda event: (refresh_contig_options(), refresh_sample_options())
        )

        conta_max = max(100, widgets['whole_contamination_max'])
        pct_contamination_slider = create_editable_range_slider(
            "Whole contamination (%)", 0, conta_max, (0, conta_max), step=1,
            on_change=lambda event: (refresh_contig_options(), refresh_sample_options())
        )

        prevalence_left_slider = create_editable_range_slider(
            "Left completeness (%)", 0, 100, (0, 100), step=1,
            on_change=lambda event: (refresh_contig_options(), refresh_sample_options())
        )

        prevalence_right_slider = create_editable_range_slider(
            "Right completeness (%)", 0, 100, (0, 100), step=1,
            on_change=lambda event: (refresh_contig_options(), refresh_sample_options())
        )

        completeness_filters_children.extend([
            pct_completeness_slider,
            pct_contamination_slider,
            prevalence_left_slider,
            prevalence_right_slider
        ])

        # Create completeness filters content container (collapsible)
        completeness_filters_content = pn.Column(*completeness_filters_children, visible=False, sizing_mode="stretch_width")
        completeness_toggle_btn.on_click(make_toggle_callback(completeness_toggle_btn, completeness_filters_content))

        filtering_children.append(completeness_filters_header)
        filtering_children.append(completeness_filters_content)

    # Add "Phage mechanisms" collapsible subsection (if PhageMechanisms table has data)
    if widgets['phage_mechanisms']:
        phage_toggle_btn = Button(label="▶", width=20, height=20, button_type="primary", align="center", margin=0, stylesheets=[toggle_stylesheet])
        phage_title = Div(text="Phage mechanisms", align="center")
        phage_header = row(phage_toggle_btn, phage_title, sizing_mode="stretch_width", align="center")

        phage_mechanism_filter = Select(
            options=[""] + widgets['phage_mechanisms'],
            value="",   # default = no filter (empty = all)
            sizing_mode="stretch_width",
        )
        phage_mechanism_filter.on_change('value', lambda attr, old, new: (refresh_contig_options(), refresh_sample_options()))

        phage_content = pn.Column(phage_mechanism_filter, visible=False, sizing_mode="stretch_width")
        phage_toggle_btn.on_click(make_toggle_callback(phage_toggle_btn, phage_content))

        filtering_children.append(phage_header)
        filtering_children.append(phage_content)

    # Add "Per variable" collapsible filtering subsection
    combined_help = "Filter by variable characteristics.\n\n#Points: Filter by number of data points.\n  - One sample view: Show contigs where more/less than threshold points exist for the variable in selected sample\n  - All samples view: Show contigs where at least one sample has more/less than threshold points\n\nMax: Filter by maximum value.\n  - One sample view: Show contigs where the variable reaches above/below threshold value in selected sample\n  - All samples view: Show contigs where the variable reaches above/below threshold value in at least one sample"
    tooltip = Tooltip(content=combined_help, position="right")
    help_per_variable = HelpButton(tooltip=tooltip, width=20, height=20, align="center", button_type="light", stylesheets=[toggle_stylesheet])
    per_variable_toggle_btn = Button(label="▶", width=20, height=20, button_type="primary", align="center", margin=0, stylesheets=[toggle_stylesheet])
    per_variable_title = Div(text="Per variable", align="center")
    per_variable_header = row(per_variable_toggle_btn, per_variable_title, help_per_variable, sizing_mode="stretch_width", align="center")

    # Store all variable filter rows for dynamic management (already initialized above)
    variable_filters_column = column(sizing_mode="stretch_width")

    # Create initial filter row
    variable_filter_rows = []
    initial_row = create_variable_filter_row()
    variable_filter_rows.append(initial_row)
    variable_filters_column.children = [initial_row]

    per_variable_content = pn.Column(variable_filters_column, visible=False, sizing_mode="stretch_width")
    per_variable_toggle_btn.on_click(make_toggle_callback(per_variable_toggle_btn, per_variable_content))

    filtering_children.append(per_variable_header)
    filtering_children.append(per_variable_content)

    filtering_content = pn.Column(
        *filtering_children,
        visible=True, sizing_mode="stretch_width"  # Always visible, no toggle
    )


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
    widgets['sample_select'].on_change('value', lambda attr, old, new: refresh_contig_options())


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

    widgets['contig_select'].on_change('value', lambda attr, old, new: refresh_sample_options())


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
    if genome_section is not None:
        below_contig_children = list(below_contig_children) + [genome_section]

    below_contig_content = column(
        *below_contig_children,
        visible=True, sizing_mode="stretch_width"
    )
    contig_toggle_btn.on_click(make_toggle_callback(contig_toggle_btn, below_contig_content))


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
    separator_filtering = Div(text="", height=1, width=350, styles={'background-color': '#333', 'margin': '10px 0'})
    separator_samples = Div(text="", height=1, width=350, styles={'background-color': '#333', 'margin': '10px 0'})
    separator_contigs = Div(text="", height=1, width=350, styles={'background-color': '#333', 'margin': '10px 0'})
    separator_variables = Div(text="", height=1, width=350, styles={'background-color': '#333', 'margin': '10px 0'})

    # Gene map is now part of the Genome module's CheckboxButtonGroup
    # Include both variables sections - visibility is toggled by on_view_change
    controls_children = [logo, views, separator_filtering, filtering_header, filtering_content,
                         separator_contigs, contig_header, above_contig_content, widgets['contig_select'], below_contig_content,
                         separator_samples, sample_header, above_sample_content, widgets['sample_select'],
                         separator_variables,
                         variables_section_one,  # One Sample view (with module checkboxes)
                         variables_section_all,  # All Samples view (title headers only)
                         buttons_row]

    controls_column = pn.Column(*controls_children, width=350, sizing_mode="stretch_height", css_classes=["left-col"])

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
        title="MGFeatureViewer",
        static_dirs={'assets': static_path}
    )
    return 0

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    add_serve_args(parser)
    args = parser.parse_args()
    raise SystemExit(run_serve(args))