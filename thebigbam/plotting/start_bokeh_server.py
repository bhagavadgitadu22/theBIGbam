import argparse
import base64
import os
import duckdb
import traceback

import panel as pn

from bokeh.layouts import column, row
from bokeh.models import Div, InlineStyleSheet, Tooltip
from bokeh.models.widgets import CheckboxGroup, HelpButton, Button, RadioButtonGroup, CheckboxButtonGroup, Select, TextInput, Spinner, MultiChoice

# Import the plotting function from the repo
from .plotting_data_per_sample import generate_bokeh_plot_per_sample
from .plotting_data_all_samples import generate_bokeh_plot_all_samples
from ..database.database_getters import get_filtering_metadata, ANNOTATION_EXCLUDED_COLUMNS
from .searchable_select import SearchableSelect

def build_controls(conn):
    """Query DB and return widgets and helper mappings."""
    cur = conn.cursor()

    # Check if PhageMechanisms table exists and has data
    phage_mechanisms_list = []
    try:
        cur.execute("SELECT DISTINCT Packaging_mechanism FROM PhageMechanisms WHERE Packaging_mechanism IS NOT NULL")
        phage_mechanisms_list = [r[0] for r in cur.fetchall()]
    except Exception:
        pass  # Table doesn't exist or has no data

    # Check if Misassembly table exists and has data
    has_misassembly = False
    try:
        cur = conn.cursor()
        cur.execute("SELECT 1 FROM information_schema.tables WHERE table_name = 'Misassembly'")
        has_misassembly = cur.fetchone() is not None
    except Exception:
        pass

    # Get Coverage_mean min/max from Coverage table (stored as INTEGER)
    coverage_mean_max = 0
    try:
        cur.execute("SELECT MIN(Coverage_mean), MAX(Coverage_mean) FROM Coverage WHERE Coverage_mean IS NOT NULL")
        result = cur.fetchone()
        if result and result[0] is not None and result[1] is not None:
            coverage_mean_min = result[0]
            coverage_mean_max = result[1]
    except Exception:
        pass

    # Get Coverage_median min/max from Coverage table (stored as INTEGER)
    coverage_median_max = 0
    try:
        cur.execute("SELECT MIN(Coverage_median), MAX(Coverage_median) FROM Coverage WHERE Coverage_median IS NOT NULL")
        result = cur.fetchone()
        if result and result[0] is not None and result[1] is not None:
            coverage_median_min = result[0]
            coverage_median_max = result[1]
    except Exception:
        pass

    # Get Coverage_variation min/max from Explicit_coverage view (already divided by 1000000)
    coverage_variation_max = 1.0
    try:
        cur.execute("SELECT MIN(Coverage_variation), MAX(Coverage_variation) FROM Explicit_coverage WHERE Coverage_variation IS NOT NULL")
        result = cur.fetchone()
        if result and result[0] is not None and result[1] is not None:
            coverage_variation_min = result[0]
            coverage_variation_max = result[1]
    except Exception:
        pass

    # Get Coverage_sd min/max from Explicit_coverage view (already divided by 1000000)
    coverage_sd_max = 1.0
    try:
        cur.execute("SELECT MIN(Coverage_sd), MAX(Coverage_sd) FROM Explicit_coverage WHERE Coverage_sd IS NOT NULL")
        result = cur.fetchone()
        if result and result[0] is not None and result[1] is not None:
            coverage_sd_min = result[0]
            coverage_sd_max = result[1]
    except Exception:
        pass

    whole_expansion_max = 100.0
    try:
        cur.execute("SELECT MAX(Expansion_per_100kbp) FROM Explicit_misassembly WHERE Expansion_per_100kbp IS NOT NULL")
        result = cur.fetchone()
        if result and result[0] is not None:
            whole_expansion_max = result[0]
    except Exception:
        pass

    # Get Duplication_percentage min/max from Contig table
    duplication_percentage_min = 0
    duplication_percentage_max = 1000
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
    try:
        # Get column names from Contig_annotation table
        cur.execute("PRAGMA table_info(Contig_annotation)")
        columns = [row[1] for row in cur.fetchall() if row[1] not in ANNOTATION_EXCLUDED_COLUMNS]
        
        # Get distinct non-null values for each column
        for column in columns:
            cur.execute(f'SELECT DISTINCT "{column}" FROM Contig_annotation WHERE "{column}" IS NOT NULL ORDER BY "{column}"')
            values = [str(r[0]) for r in cur.fetchall()]
            if values:  # Only add if there are values
                annotation_filters[column] = values
    except Exception:
        pass  # Table doesn't exist or has no data

    # Get annotation feature types from Annotated_types table
    annotation_types = []
    try:
        cur.execute("SELECT Type_name FROM Annotated_types ORDER BY Frequency DESC")
        annotation_types = [r[0] for r in cur.fetchall()]
    except Exception:
        pass

    # Widget Selector for Contigs (autocomplete with max 20 suggestions)
    cur.execute("SELECT Contig_name, Contig_length, Duplication_percentage FROM Contig ORDER BY Contig_name")
    rows = cur.fetchall()
    contigs = [r[0] for r in rows]
    contig_lengths = {r[0]: r[1] for r in rows}  # Dictionary mapping contig_name -> length
    contig_duplications = {r[0]: r[2] for r in rows}  # Dictionary mapping contig_name -> duplication percentage (can be None)
    
    # If only one contig in database, pre-fill the field
    contig_select = SearchableSelect(
        value=contigs[0] if len(contigs) == 1 else "",
        options=contigs,
        placeholder="Type to search contigs...",
        sizing_mode="stretch_width",
        margin=(0, 5, 0, 5)
    )

    # Widget Selector for Samples (autocomplete with max 20 suggestions)
    cur.execute("SELECT Sample_name FROM Sample ORDER BY Sample_name")
    samples = [r[0] for r in cur.fetchall()]
    has_samples = len(samples) > 0
    
    # If only one sample in database, pre-fill the field
    sample_select = SearchableSelect(
        value=samples[0] if len(samples) == 1 else "",
        options=samples,
        placeholder="Type to search samples...",
        sizing_mode="stretch_width",
        margin=(0, 5, 0, 5)
    )

    # Build presence mappings: sample -> contigs and contig -> samples
    cur.execute("""
    SELECT Contig.Contig_name, Sample.Sample_name FROM Coverage
      JOIN Contig ON Coverage.Contig_id = Contig.Contig_id
      JOIN Sample ON Coverage.Sample_id = Sample.Sample_id
    """)
    sample_to_contigs = {}
    contig_to_samples = {}
    for contig_name, sample_name in cur.fetchall():
        sample_to_contigs.setdefault(sample_name, set()).add(contig_name)
        contig_to_samples.setdefault(contig_name, set()).add(sample_name)

    # Get variables that have data (their feature table exists)
    cur.execute("SELECT DISTINCT Variable_name, Feature_table_name FROM Variable")
    rows = cur.fetchall()
    variables = [r[0] for r in rows]
    tables_with_data = [r[1] for r in rows]

    # Get modules that have at least one variable with data
    # Define the display order for modules
    MODULE_ORDER = ["Genome", "Coverage", "Misalignment", "Long-reads", "Paired-reads", "Phage termini"]

    cur.execute("SELECT DISTINCT Module FROM Variable WHERE Feature_table_name IN ({})".format(
        ','.join('?' * len(tables_with_data))
    ), tuple(tables_with_data))  # Pass only variable names to match with Feature_table_name
    modules_from_db = [r[0] for r in cur.fetchall()]

    # Sort modules according to MODULE_ORDER, keeping any unknown modules at the end
    modules = sorted(modules_from_db, key=lambda m: MODULE_ORDER.index(m) if m in MODULE_ORDER else len(MODULE_ORDER))

    # Identify custom contig-level subplots (Contig_* tables in Custom module)
    cur.execute("SELECT Subplot FROM Variable WHERE Module='Custom' AND Feature_table_name LIKE 'Contig_%'")
    custom_contig_subplots = [r[0] for r in cur.fetchall()]

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

        # For Genome module, "Gene map" is handled separately (not in the checkbox group)
        # So we don't add it here - it will get its own dedicated button

        # For Custom module, exclude contig-level subplots (they go in genome section)
        if module == "Custom" and custom_contig_subplots:
            variables_checkbox = [v for v in variables_checkbox if v not in custom_contig_subplots]

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
        'custom_contig_subplots': custom_contig_subplots,
        'phage_mechanisms': phage_mechanisms_list,
        'has_misassembly': has_misassembly,
        'coverage_mean_max': coverage_mean_max,
        'coverage_median_max': coverage_median_max,
        'coverage_variation_max': coverage_variation_max,
        'coverage_sd_max': coverage_sd_max,
        'whole_expansion_max': whole_expansion_max,
        'duplication_percentage_min': duplication_percentage_min,
        'duplication_percentage_max': duplication_percentage_max,
        'has_duplication_data': has_duplication_data,
        'annotation_filters': annotation_filters,
        'annotation_types': annotation_types,
        'full_contigs': contigs,  # Store full list for substring matching
        'full_samples': samples,  # Store full list for substring matching
        'has_samples': has_samples  # True if database has any samples
    }
    return widgets

def create_layout(db_path):
    """Create and return the application layout for Panel serve."""

    ### Event functions
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

    _filtering_cache = {'result': None, 'valid': False}

    def get_filtering_filtered_pairs():
        """Apply Filtering query rows to get allowed contig/sample pairs.

        Returns set of (contig_name, sample_name) tuples that match all conditions.
        Returns None if no filters are active (meaning all pairs are allowed).
        """
        if _filtering_cache['valid']:
            return _filtering_cache['result']

        if not or_sections:
            _filtering_cache['result'] = None
            _filtering_cache['valid'] = True
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
                has_active_filter = True
                break
            if has_active_filter:
                break

        if not has_active_filter:
            _filtering_cache['result'] = None
            _filtering_cache['valid'] = True
            return None

        cur = conn.cursor()

        # Source table mapping from filtering_metadata
        source_table_map = {
            'Contig': 'Contig',
            'Sample': 'Sample',
            'Coverage': 'Explicit_coverage',
            'Misassembly': 'Explicit_misassembly',
            'Microdiversity': 'Explicit_microdiversity',
            'Side misassembly': 'Explicit_side_misassembly',
            'Topology': 'Explicit_topology',
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
                    LEFT JOIN Coverage p ON c.Contig_id = p.Contig_id
                    LEFT JOIN Sample s ON p.Sample_id = s.Sample_id
                    WHERE ca."{column_name}" {operator} ?
                '''
            elif category == 'Contig':
                # Contig table has no Sample_name, left-join to preserve contigs with 0 samples
                query = f'''
                    SELECT DISTINCT c.Contig_name, s.Sample_name
                    FROM Contig c
                    LEFT JOIN Coverage p ON c.Contig_id = p.Contig_id
                    LEFT JOIN Sample s ON p.Sample_id = s.Sample_id
                    WHERE c."{column_name}" {operator} ?
                '''
            elif category == 'Sample':
                # Sample table has no Contig_name, left-join to preserve samples with 0 contigs
                query = f'''
                    SELECT DISTINCT c.Contig_name, s.Sample_name
                    FROM Sample s
                    LEFT JOIN Coverage p ON s.Sample_id = p.Sample_id
                    LEFT JOIN Contig c ON p.Contig_id = c.Contig_id
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
                # Check inter-section AND/OR select widget
                if i - 1 < len(inter_section_selects) and inter_section_selects[i - 1].value == "OR":
                    final_pairs = final_pairs | section_pairs
                else:
                    final_pairs = final_pairs & section_pairs

        _filtering_cache['result'] = final_pairs
        _filtering_cache['valid'] = True
        return final_pairs

    def update_widget_completions(widget, completions):
        """Update widget completions. Clear value if not in completions."""
        widget.options = completions
        if widget.value and widget.value not in completions:
            widget.value = ""
    
    def refresh_contig_options_unlocked():
        """Core logic — does NOT check the lock."""
        # Apply presence filter only if a sample is selected (One Sample view)
        if views.active == 0 and widgets['sample_select'].value:
            sel_sample = widgets['sample_select'].value
            allowed = widgets['sample_to_contigs'].get(sel_sample, set())
            completions = [c for c in orig_contigs if c in allowed]
        else:
            completions = list(orig_contigs)

        # Apply Filtering2 query builder filters
        filtered_pairs = get_filtering_filtered_pairs()
        if filtered_pairs is not None:
            allowed_contigs = {pair[0] for pair in filtered_pairs}
            completions = [c for c in completions if c in allowed_contigs]

        update_widget_completions(widgets['contig_select'], completions)

    def refresh_sample_options_unlocked():
        """Core logic — does NOT check the lock."""
        # Apply presence filter only if a contig is selected (One Sample view)
        if views.active == 0 and widgets['contig_select'].value:
            sel_contig = widgets['contig_select'].value
            allowed = widgets['contig_to_samples'].get(sel_contig, set())
            completions = [s for s in orig_samples if s in allowed]
        else:
            completions = list(orig_samples)

        # Apply Filtering2 query builder filters
        filtered_pairs = get_filtering_filtered_pairs()
        if filtered_pairs is not None:
            allowed_samples = {pair[1] for pair in filtered_pairs}
            completions = [s for s in completions if s in allowed_samples]

        update_widget_completions(widgets['sample_select'], completions)

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

        filtering_title.text = f"<span style='font-size: 1.2em;'><b>Filtering</b></span> ({presences_count} contig/sample pairs)"
        contig_title.text = f"<span style='font-size: 1.2em;'><b>Contigs</b></span> ({contigs_count} available)"
        sample_title.text = f"<span style='font-size: 1.2em;'><b>Samples</b></span> ({samples_count} available)"

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
        sample_title.visible = not is_all
        above_sample_content.visible = not is_all
        widgets['sample_select'].visible = not is_all

        # Toggle visibility between the two variables sections
        # Each section maintains its own state independently
        variables_section_one.visible = not is_all
        variables_section_all.visible = is_all

        # Refresh options while still locked (suppresses cascading callbacks)
        _filtering_cache['valid'] = False
        refresh_contig_options_unlocked()
        if not is_all:
            refresh_sample_options_unlocked()

        # Unlock AFTER refreshes complete
        global_toggle_lock['locked'] = False
        update_section_titles()

    ## Apply button function
    def apply_clicked():
        try:
            contig = widgets['contig_select'].value
            has_samples = widgets['has_samples']
            
            # When no samples exist, treat as "One Sample" mode with no sample/variables
            is_all = (views.active == 1) if has_samples else False
            sample = widgets['sample_select'].value if has_samples else None

            # Select the correct widget set based on current view
            active_variables_widgets = widgets['variables_widgets_all'] if is_all else widgets['variables_widgets_one']

            # Genome module is shared between views (in Contigs section)
            # Gene map is shown if at least one feature type is selected in the multichoice
            selected_feature_types = feature_type_multichoice.value if feature_type_multichoice is not None else None
            genbank_path = db_path if (selected_feature_types and len(selected_feature_types) > 0) else None
            use_phage_colors = (0 in phage_colors_cbg.active) if (phage_colors_cbg is not None and genbank_path) else False
            plot_isoforms = (0 in plot_isoforms_cbg.active) if (plot_isoforms_cbg is not None and genbank_path) else True

            # Parse and validate position inputs
            xstart = None
            xend = None
            
            # Validate contig is selected
            if not contig:
                peruse_button.visible = False
                main_placeholder.objects = [pn.pane.HTML("<pre>Error: Please select a contig.</pre>")]
                return
            
            # Get contig length for validation
            contig_length = widgets['contig_lengths'].get(contig, 0)
            
            # Parse position inputs
            try:
                xstart = int(from_position_input.value) if from_position_input.value.strip() else 0
                xend = int(to_position_input.value) if to_position_input.value.strip() else contig_length
            except ValueError:
                peruse_button.visible = False
                main_placeholder.objects = [pn.pane.HTML("<pre>Error: Invalid position range - positions must be integers.</pre>")]
                return
            
            # Validate position range
            if xstart >= xend:
                peruse_button.visible = False
                main_placeholder.objects = [pn.pane.HTML(f"<pre>Error: Invalid position range - start must be less than end.</pre>")]
                return

            # Check whether to plot sequence (checkbox checked)
            plot_sequence = (
                sequence_cbg is not None
                and 0 in sequence_cbg.active
            )

            # Check whether to use same y scale for all subplots
            same_y_scale = (0 in same_y_scale_cbg.active)

            # Read subplot height from spinner
            subplot_size = int(subplot_height_input.value)
            genemap_size = int(genemap_height_input.value)

            # Check whether to preserve x-range from previous plot
            preserve_xrange = (
                current_plot_state['shared_xrange'] is not None
                and current_plot_state['contig'] == contig
                and current_plot_state['is_all'] == is_all
                and (is_all or current_plot_state['sample'] == sample)
                and current_plot_state['data_xstart'] == xstart
                and current_plot_state['data_xend'] == xend
            )
            if preserve_xrange:
                prev_xstart = current_plot_state['shared_xrange'].start
                prev_xend = current_plot_state['shared_xrange'].end

            if is_all:
                # All-samples view: require exactly one variable selected from non-Genome modules
                # Genome module is shared (in Contigs section):
                #   - Gene map (index 0) is handled via genbank_path
                #   - Other Genome features (Repeats, etc.) are passed via genome_features parameter
                selected_var = None
                genome_features = []

                # Collect genomic features from combined_features_cbg (Genome features + Custom contig features)
                if combined_features_cbg is not None:
                    for idx in combined_features_cbg.active:
                        genome_features.append(combined_features_cbg.labels[idx])

                # Get the selected variable from non-Genome modules
                for cbg in active_variables_widgets:
                    if cbg.active and selected_var is None:
                        selected_var = cbg.labels[cbg.active[-1]]

                if not selected_var and not genome_features:
                    raise ValueError("When in 'All samples' view you must select at least one variable to plot.")

                # Compute filtered samples (same logic as refresh_sample_options)
                # Start with samples that have the selected contig
                filtered_samples = [s for s in orig_samples if s in widgets['contig_to_samples'].get(contig, set())]
                # Apply Filtering2 query builder conditions
                filtering_pairs = get_filtering_filtered_pairs()
                if filtering_pairs is not None:
                    allowed_samples = {pair[1] for pair in filtering_pairs}
                    filtered_samples = [s for s in filtered_samples if s in allowed_samples]

                print(f"[start_bokeh_server] Generating plot for all samples with variable={selected_var}, contig={contig}, genome_features={genome_features}, filtered_samples={len(filtered_samples)}")
                grid = generate_bokeh_plot_all_samples(conn, selected_var, contig, xstart=xstart, xend=xend, genbank_path=genbank_path, genome_features=genome_features if genome_features else None, allowed_samples=set(filtered_samples), feature_types=selected_feature_types, use_phage_colors=use_phage_colors, plot_sequence=plot_sequence, same_y_scale=same_y_scale, subplot_size=subplot_size, genemap_size=genemap_size)
            else:
                # One-sample view: collect possibly-many requested features and call per-sample plot
                requested_features = []

                # Collect genomic features from combined_features_cbg (in Contigs section)
                if combined_features_cbg is not None:
                    for idx in combined_features_cbg.active:
                        requested_features.append(combined_features_cbg.labels[idx])

                # Collect features from Variables section
                for cbg in active_variables_widgets:
                    for idx in cbg.active:
                        requested_features.append(cbg.labels[idx])

                print(f"[start_bokeh_server] Generating plot for sample={sample}, contig={contig}, features={requested_features}")
                grid = generate_bokeh_plot_per_sample(conn, requested_features, contig, sample, xstart=xstart, xend=xend, genbank_path=genbank_path, feature_types=selected_feature_types, use_phage_colors=use_phage_colors, plot_isoforms=plot_isoforms, plot_sequence=plot_sequence, same_y_scale=same_y_scale, subplot_size=subplot_size, genemap_size=genemap_size)

            # Extract and cache plot data for instant downloads (no database round-trip)
            try:
                plot_data_cache['data'] = extract_plot_data_from_grid(grid)
                plot_data_cache['is_all_samples'] = is_all
                # Set filename based on context
                safe_contig = "".join(c if c.isalnum() or c in "-_" else "_" for c in contig)
                if is_all:
                    safe_var = "".join(c if c.isalnum() or c in "-_" else "_" for c in (selected_var or 'data'))
                    plot_data_cache['filename'] = f"{safe_contig}_all_samples_{safe_var}.csv"
                else:
                    safe_sample = "".join(c if c.isalnum() or c in "-_" else "_" for c in sample)
                    plot_data_cache['filename'] = f"{safe_contig}_{safe_sample}_data.csv"
                print(f"[start_bokeh_server] Cached {len(plot_data_cache['data'])} data sources for download", flush=True)
            except Exception as cache_err:
                print(f"[start_bokeh_server] Warning: Could not cache plot data: {cache_err}", flush=True)
                plot_data_cache['data'] = None

            # Restore preserved x-range and update state
            new_xrange = _get_shared_xrange(grid)
            if preserve_xrange and new_xrange is not None:
                new_xrange.start = prev_xstart
                new_xrange.end = prev_xend

            current_plot_state['contig'] = contig
            current_plot_state['sample'] = sample
            current_plot_state['is_all'] = is_all
            current_plot_state['shared_xrange'] = new_xrange
            current_plot_state['data_xstart'] = xstart
            current_plot_state['data_xend'] = xend

            # Sync From/To inputs when user zooms/pans
            if new_xrange is not None:
                def _sync_from(attr, old, new):
                    from_position_input.value = str(int(new))

                def _sync_to(attr, old, new):
                    to_position_input.value = str(int(new))

                new_xrange.on_change('start', _sync_from)
                new_xrange.on_change('end', _sync_to)

            # Create toolbar-style row with buttons positioned top-right
            # Use Panel Row to mix Bokeh and Panel widgets
            toolbar_row = pn.Row(
                pn.Spacer(sizing_mode="stretch_width"),  # Push buttons to right
                peruse_button,
                download_contig_button,
                download_metrics_button,
                download_data_button,
                margin=(0, 0, 5, 0)
            )
            # Show buttons when plot exists (but some are hidden in 0-sample mode)
            download_contig_button.visible = True  # Always show download contig button
            if has_samples:
                peruse_button.visible = True
                download_metrics_button.visible = True
                download_data_button.visible = True
            else:
                # No samples: hide sample-related buttons
                peruse_button.visible = False
                download_metrics_button.visible = False
                download_data_button.visible = False
            
            # Stack toolbar row above grid
            main_placeholder.objects = [pn.Column(toolbar_row, grid, sizing_mode="stretch_both")]

        except Exception as e:
            peruse_button.visible = False
            download_contig_button.visible = False
            download_metrics_button.visible = False
            download_data_button.visible = False
            tb = traceback.format_exc()
            print(f"[start_bokeh_server] Exception: {tb}", flush=True)
            main_placeholder.objects = [pn.pane.HTML(f"<pre>Error building plot:\n{tb}</pre>")]

    ## Peruse button callback function
    def peruse_clicked():
        """Generate and open summary tables in a new browser window."""
        from .perusing_data import generate_and_open_peruse_html
        
        contig = widgets['contig_select'].value
        has_samples = widgets['has_samples']

        # Check if contig is selected
        if not contig:
            print("[start_bokeh_server] Peruse: No contig selected", flush=True)
            return

        # If no samples in database, cannot peruse
        if not has_samples:
            print("[start_bokeh_server] Peruse: No samples in database", flush=True)
            return

        is_all = (views.active == 1)

        if is_all:
            # All Samples view: get filtered samples
            filtered_samples = [s for s in orig_samples if s in widgets['contig_to_samples'].get(contig, set())]
            # Apply Filtering2 query builder conditions
            filtering_pairs = get_filtering_filtered_pairs()
            if filtering_pairs is not None:
                allowed_samples = {pair[1] for pair in filtering_pairs}
                filtered_samples = [s for s in filtered_samples if s in allowed_samples]

            if not filtered_samples:
                print("[start_bokeh_server] Peruse: No samples match filters", flush=True)
                return

            sample_names = filtered_samples
        else:
            # One Sample view: use selected sample
            sample = widgets['sample_select'].value
            if not sample:
                print("[start_bokeh_server] Peruse: No sample selected", flush=True)
                return
            sample_names = [sample]

        # Generate and open HTML in new window
        generate_and_open_peruse_html(conn, contig, sample_names)

    ## Download functionality using Panel FileDownload widgets
    import io
    
    # Store references to download widgets (created later, after widgets dict exists)
    download_widgets = {'contig': None, 'metrics': None, 'data': None}
    
    # Cache for plot data - extracted from ColumnDataSources after plotting
    # This avoids re-querying the database for downloads (instant downloads)
    plot_data_cache = {'data': None, 'is_all_samples': False, 'filename': 'data.csv'}

    # Track current plot state for x-range preservation across APPLY clicks
    current_plot_state = {
        'contig': None,
        'sample': None,
        'is_all': None,
        'shared_xrange': None,
        'data_xstart': None,
        'data_xend': None,
    }

    def _get_shared_xrange(grid):
        """Extract the shared Range1d from a gridplot's first figure."""
        if hasattr(grid, 'children'):
            for child_spec in grid.children:
                child = child_spec[0] if isinstance(child_spec, tuple) else child_spec
                if hasattr(child, 'x_range'):
                    return child.x_range
        return None

    def extract_plot_data_from_grid(grid):
        """Extract all data from ColumnDataSources in a GridPlot.
        
        Returns a list of dicts, each containing data from one renderer's source.
        """
        from bokeh.models import GlyphRenderer
        
        all_data = []
        figures = []
        
        # Get all figures from the grid
        if hasattr(grid, 'children'):
            for item in grid.children:
                if hasattr(item, '__iter__'):
                    for subitem in item:
                        if hasattr(subitem, 'renderers'):
                            figures.append(subitem)
                elif hasattr(item, 'renderers'):
                    figures.append(item)
        
        for fig in figures:
            if not hasattr(fig, 'renderers'):
                continue
            
            # Get legend label from figure's legend if available
            legend_map = {}
            if hasattr(fig, 'legend') and fig.legend:
                for legend in fig.legend:
                    if hasattr(legend, 'items'):
                        for item in legend.items:
                            if hasattr(item, 'label') and hasattr(item, 'renderers'):
                                # Extract string from Bokeh Value object
                                label_obj = item.label
                                if hasattr(label_obj, 'value'):
                                    label_val = label_obj.value  # Bokeh Value object
                                elif isinstance(label_obj, dict):
                                    label_val = label_obj.get('value', '')
                                else:
                                    label_val = str(label_obj)
                                for r in item.renderers:
                                    legend_map[id(r)] = label_val
            
            for renderer in fig.renderers:
                if isinstance(renderer, GlyphRenderer) and hasattr(renderer, 'data_source'):
                    source = renderer.data_source
                    if hasattr(source, 'data') and source.data:
                        data = {k: list(v) for k, v in source.data.items()}  # Copy data
                        # Skip empty sources or sources without x/y
                        if 'x' not in data or 'y' not in data:
                            continue
                        if not data['x'] or not data['y']:
                            continue
                        # Get feature name from legend
                        feature_name = legend_map.get(id(renderer), 'unknown')
                        data['_feature_name'] = feature_name
                        all_data.append(data)
        
        return all_data
    
    def convert_plot_data_to_csv(plot_data, is_all_samples):
        """Convert extracted plot data to CSV string."""
        if not plot_data:
            return None
        
        import csv as csv_module
        output = io.StringIO()
        writer = csv_module.writer(output)
        
        first_col = 'sample_name' if is_all_samples else 'feature_name'
        
        # Check what columns exist across all data
        has_first_last = any('first_pos' in d for d in plot_data)
        has_stats = any('mean' in d for d in plot_data)
        
        # Build header
        if has_first_last:
            if has_stats:
                header = [first_col, 'start_position', 'last_position', 'value', 'mean', 'median', 'std']
            else:
                header = [first_col, 'start_position', 'last_position', 'value']
        else:
            if has_stats:
                header = [first_col, 'position', 'value', 'mean', 'median', 'std']
            else:
                header = [first_col, 'position', 'value']
        
        writer.writerow(header)
        
        for data in plot_data:
            name = data.get('_feature_name', 'unknown')
            x_vals = data.get('x', [])
            y_vals = data.get('y', [])
            first_pos = data.get('first_pos', None)
            last_pos = data.get('last_pos', None)
            means = data.get('mean', None)
            medians = data.get('median', None)
            stds = data.get('std', None)
            
            # Step-plot data has pairs of points (start, end) with same y value
            # Consolidate pairs into single rows: start_pos, end_pos, value
            # This is faster (half the rows) and produces cleaner output
            n = len(x_vals)
            i = 0
            while i < n:
                # Check if this is a step-plot pair (consecutive points with same y)
                if i + 1 < n and y_vals[i] == y_vals[i + 1]:
                    # Step-plot pair: use x[i] as start, x[i+1] as end
                    start_x = x_vals[i]
                    end_x = x_vals[i + 1]
                    row = [name, start_x, end_x, y_vals[i]]
                    if has_stats:
                        row.append(means[i] if means and i < len(means) else '')
                        row.append(medians[i] if medians and i < len(medians) else '')
                        row.append(stds[i] if stds and i < len(stds) else '')
                    writer.writerow(row)
                    i += 2  # Skip both points of the pair
                else:
                    # Single point (not a pair)
                    if has_first_last and first_pos and last_pos:
                        row = [name, first_pos[i], last_pos[i], y_vals[i]]
                    else:
                        row = [name, x_vals[i], x_vals[i], y_vals[i]]
                    if has_stats:
                        row.append(means[i] if means and i < len(means) else '')
                        row.append(medians[i] if medians and i < len(medians) else '')
                        row.append(stds[i] if stds and i < len(stds) else '')
                    writer.writerow(row)
                    i += 1
        
        return output.getvalue()

    def make_contig_download_callback():
        """Create callback for contig summary download."""
        from .downloading_data import download_contig_summary_csv
        
        contig = widgets['contig_select'].value
        if not contig:
            print("[start_bokeh_server] Download: No contig selected", flush=True)
            return io.StringIO("")
        
        csv_content = download_contig_summary_csv(db_path, contig)
        if csv_content:
            safe_contig = "".join(c if c.isalnum() or c in "-_" else "_" for c in contig)
            if download_widgets['contig']:
                download_widgets['contig'].filename = f"{safe_contig}_contig_summary.csv"
            return io.StringIO(csv_content)
        return io.StringIO("")

    def make_metrics_download_callback():
        """Create callback for metrics summary download."""
        from .downloading_data import download_metrics_summary_csv
        
        contig = widgets['contig_select'].value
        if not contig:
            print("[start_bokeh_server] Download metrics: No contig selected", flush=True)
            return io.StringIO("")

        is_all = (views.active == 1)

        if is_all:
            filtered_samples = [s for s in orig_samples if s in widgets['contig_to_samples'].get(contig, set())]
            filtering_pairs = get_filtering_filtered_pairs()
            if filtering_pairs is not None:
                allowed_samples = {pair[1] for pair in filtering_pairs}
                filtered_samples = [s for s in filtered_samples if s in allowed_samples]

            if not filtered_samples:
                print("[start_bokeh_server] Download metrics: No samples match filters", flush=True)
                return io.StringIO("")

            sample_names = filtered_samples
            safe_contig = "".join(c if c.isalnum() or c in "-_" else "_" for c in contig)
            if download_widgets['metrics']:
                download_widgets['metrics'].filename = f"{safe_contig}_in_all_samples_metrics.csv"
        else:
            sample = widgets['sample_select'].value
            if not sample:
                print("[start_bokeh_server] Download metrics: No sample selected", flush=True)
                return io.StringIO("")
            sample_names = [sample]
            safe_contig = "".join(c if c.isalnum() or c in "-_" else "_" for c in contig)
            safe_sample = "".join(c if c.isalnum() or c in "-_" else "_" for c in sample)
            if download_widgets['metrics']:
                download_widgets['metrics'].filename = f"{safe_contig}_in_{safe_sample}_metrics.csv"

        csv_content = download_metrics_summary_csv(db_path, contig, sample_names)
        if csv_content:
            return io.StringIO(csv_content)
        return io.StringIO("")

    def make_data_download_callback():
        """Create callback for feature data download.
        
        Uses cached plot data for instant downloads (no database round-trip).
        """
        # Use cached plot data (extracted from ColumnDataSources after plotting)
        if plot_data_cache['data']:
            csv_content = convert_plot_data_to_csv(
                plot_data_cache['data'], 
                plot_data_cache['is_all_samples']
            )
            if csv_content:
                if download_widgets['data']:
                    download_widgets['data'].filename = plot_data_cache['filename']
                print(f"[start_bokeh_server] Data download from cache ({len(plot_data_cache['data'])} sources)", flush=True)
                return io.StringIO(csv_content)
        
        # No cache available - return empty (user needs to click Apply first)
        print("[start_bokeh_server] Download data: No cached data (click Apply first)", flush=True)
        return io.StringIO("")

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
    logo = Div(text=f"""<img src="data:image/png;base64,{logo_b64}" style="width:100%; max-width:800px; padding: 0 25%;">""")
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
    # Store inter-section AND/OR Select widgets
    inter_section_selects = []

    def count_total_query_rows():
        """Count total number of query rows across all OR sections."""
        total = 0
        for section_data in or_sections:
            total += len(section_data['rows'])
        return total

    def refresh_on_filter_change():
        """Refresh contig and sample options when Filtering2 values change."""
        _filtering_cache['valid'] = False
        global_toggle_lock['locked'] = True
        refresh_contig_options_unlocked()
        refresh_sample_options_unlocked()
        global_toggle_lock['locked'] = False
        update_section_titles()

    def create_query_row(section_data):
        """Create a single query row with cascading selects, comparison, dynamic input and remove button."""
        # Get categories from metadata
        categories = list(filtering_metadata.keys())
        if not categories:
            categories = ["No data"]

        initial_category = categories[0]
        initial_columns_raw = list(filtering_metadata.get(initial_category, {}).get('columns', {}).keys())
        if not initial_columns_raw:
            initial_columns_raw = ["No columns"]
        initial_columns = [(c, c.replace("_", " ").replace("percentage", "(%)")) for c in initial_columns_raw]
        initial_column = initial_columns_raw[0]

        # Determine initial column type
        initial_col_info = filtering_metadata.get(initial_category, {}).get('columns', {}).get(initial_column, {})
        initial_is_text = initial_col_info.get('type') == 'text'

        # First level select (categories)
        category_select = Select(
            options=categories,
            value=initial_category,
            width=70,
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
            value="=" if initial_is_text else ">",
            width=50,
            margin=(0, 2, 0, 0)
        )

        # Container for the dynamic input widget
        input_container = pn.Column(width=90, margin=(0, 2, 0, 0))

        # Create initial input widget based on column type
        if initial_is_text:
            distinct_values = initial_col_info.get('distinct_values', [])
            initial_input = SearchableSelect(
                value="", options=distinct_values,
                placeholder="Search...", width=90
            )
            input_container.objects = [initial_input]
            initial_input.param.watch(lambda event: refresh_on_filter_change(), 'value')
            initial_is_panel = True
        else:
            initial_input = Spinner(value=0, placeholder="Value...", width=90, margin=(0, 2, 0, 0))
            input_container.objects = [initial_input]
            # Add callback for Bokeh Spinner
            initial_input.on_change('value', lambda attr, old, new: refresh_on_filter_change())
            initial_is_panel = False

        # Remove button (Panel button for proper dynamic event handling)
        minus_btn = pn.widgets.Button(name="−", width=30, height=30, margin=(0, 10, 0, 0), stylesheets=[stylesheet])

        # Store reference to current input widget (for later retrieval)
        current_input_ref = {'widget': initial_input, 'is_panel': initial_is_panel}

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
                new_input = SearchableSelect(
                    value="", options=distinct_values,
                    placeholder="Search...", width=90
                )
                input_container.objects = [new_input]
                current_input_ref['widget'] = new_input
                current_input_ref['is_panel'] = True
                new_input.param.watch(lambda event: refresh_on_filter_change(), 'value')
            else:
                new_input = Spinner(value=0, placeholder="Value...", width=90, margin=(0, 2, 0, 0))
                input_container.objects = [new_input]
                current_input_ref['widget'] = new_input
                current_input_ref['is_panel'] = False
                # Add callback for Bokeh Spinner
                new_input.on_change('value', lambda attr, old, new: refresh_on_filter_change())

            # Immediately apply the new default value to contig/sample filtering
            refresh_on_filter_change()

        def update_subcategories(attr, old, new):
            """Update column options when category changes."""
            columns = list(filtering_metadata.get(new, {}).get('columns', {}).keys())
            if not columns:
                columns = ["No columns"]
            subcategory_select.options = [(c, c.replace("_", " ").replace("percentage", "(%)")) for c in columns]
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

        query_row = pn.Row(category_select, subcategory_select, comparison_select, input_container, minus_btn,
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

        def remove_row_callback(event):
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
                def _on_and_or_change(attr, old, new):
                    _filtering_cache['valid'] = False
                    global_toggle_lock['locked'] = True
                    refresh_contig_options_unlocked()
                    refresh_sample_options_unlocked()
                    global_toggle_lock['locked'] = False
                    update_section_titles()
                select_widget.on_change('value', _on_and_or_change)
                section_children.append(select_widget)
                row_data['and_div'] = select_widget
            else:
                row_data['and_div'] = None

            section_children.append(row_data['query_row'])

        # Add the "+ Add AND" button
        section_children.append(section_data['add_and_btn'])

        # Update the section column's objects
        section_data['column'].objects = section_children

    def create_or_section():
        """Create a new OR section with one query row and Add AND/OR button."""
        section_column = pn.Column(
            sizing_mode="stretch_width",
            styles={'border-left': '3px solid #00b17c', 'padding-left': '10px', 'margin-left': '5px'}
        )

        # Use Panel button instead of Bokeh button for proper dynamic event handling
        add_and_btn = pn.widgets.Button(
            name="+ Add AND/OR",
            margin=(5, 0, 5, 0),
            button_type="success",
            stylesheets=[stylesheet]
        )

        section_data = {
            'column': section_column,
            'rows': [],
            'add_and_btn': add_and_btn
        }

        def add_and_or_callback(event):
            # Create a new query row
            new_row = create_query_row(section_data)
            section_data['rows'].append(new_row)
            rebuild_section(section_data)
            refresh_on_filter_change()

        add_and_btn.on_click(add_and_or_callback)

        # Create initial query row
        initial_row = create_query_row(section_data)
        section_data['rows'].append(initial_row)

        rebuild_section(section_data)

        return section_data

    # Create the global "+ Add AND/OR" button (Panel button for proper dynamic handling)
    global_add_btn = pn.widgets.Button(
        name="+ Add AND/OR",
        margin=(10, 0, 5, 0),
        button_type="primary",
        stylesheets=[pink_buttons_stylesheet]
    )

    # Store reference to the button widget for replacement
    global_widget_state = {'widget': global_add_btn}

    def rebuild_filtering_content():
        """Rebuild the entire Filtering content with all OR sections."""
        content_children = []
        inter_section_selects.clear()

        for i, section_data in enumerate(or_sections):
            # Add OR div before each section except the first
            if i > 0:
                select_widget = Select(
                    options=["AND", "OR"],
                    value="AND",
                    margin=(5, 0, 5, 0)
                )
                # Add callback to refresh when AND/OR changes
                def _on_and_or_change(attr, old, new):
                    _filtering_cache['valid'] = False
                    global_toggle_lock['locked'] = True
                    refresh_contig_options_unlocked()
                    refresh_sample_options_unlocked()
                    global_toggle_lock['locked'] = False
                    update_section_titles()
                select_widget.on_change('value', _on_and_or_change)
                inter_section_selects.append(select_widget)
                content_children.append(select_widget)

            content_children.append(section_data['column'])

        # Add the current global widget (button or select) at the end
        content_children.append(global_widget_state['widget'])

        filtering_content.objects = content_children

    def global_add_and_or_callback(event):
        # Create a new section
        new_section = create_or_section()
        or_sections.append(new_section)
        rebuild_filtering_content()
        refresh_on_filter_change()

    global_add_btn.on_click(global_add_and_or_callback)

    # Create initial OR section
    initial_section = create_or_section()
    or_sections.append(initial_section)

    # Create the main content container (Panel Column to support Panel widgets inside)
    filtering_content = pn.Column(
        sizing_mode="stretch_width"
    )

    # Initial build of content
    rebuild_filtering_content()

    # Add toggle callback for collapsible Filtering section
    filtering_toggle_btn.on_click(make_toggle_callback(filtering_toggle_btn, filtering_content))



    ## Build Sample section
    #sample_toggle_btn = Button(label="▼", width=20, height=20, button_type="primary", align="center", margin=0, stylesheets=[toggle_stylesheet])
    sample_title = Div(text="<b>Samples</b>")
    #sample_header = row(sample_toggle_btn, sample_title, sizing_mode="stretch_width", align="center", margin=(0, 0, 5, 0))

    above_sample_children = []
    above_sample_content = column(
        *above_sample_children,
        visible=True, sizing_mode="stretch_width"
    )

    #sample_toggle_btn.on_click(make_toggle_callback(sample_toggle_btn, above_sample_content))
    def _on_sample_change(event):
        global_toggle_lock['locked'] = True
        refresh_contig_options_unlocked()
        global_toggle_lock['locked'] = False
        update_section_titles()
    widgets['sample_select'].param.watch(_on_sample_change, 'value')


    ## Build Contig section
    # Keep original full lists so we can restore when filters are off
    orig_contigs = list(widgets['contigs'])
    orig_samples = list(widgets['samples'])

    # Create Contigs section header
    contig_toggle_btn = Button(label="▼", width=20, height=20, button_type="primary", align="center", margin=0, stylesheets=[toggle_stylesheet])
    contig_toggle_btn.styles = {'padding': '0px', 'line-height': '20px'}
    contig_title = Div(text="<b>Contigs</b>", align="center")
    contig_header = row(contig_toggle_btn, contig_title, sizing_mode="stretch_width", align="center", margin=(0, 0, 0, 0))
    above_contig_children = []
    
    above_contig_content = column(
        *above_contig_children,
        visible=True, sizing_mode="stretch_width", margin=(0, 0, 0, 0)
    )
    contig_toggle_btn.on_click(make_toggle_callback(contig_toggle_btn, above_contig_content))

    def on_contig_change(event):
        new = event.new
        global_toggle_lock['locked'] = True
        refresh_sample_options_unlocked()
        global_toggle_lock['locked'] = False
        update_section_titles()
        # Update position inputs when contig changes
        if new and new in widgets['contig_lengths']:
            from_position_input.value = "0"
            to_position_input.value = str(widgets['contig_lengths'][new])
        else:
            from_position_input.value = "0"
            to_position_input.value = ""
    
    widgets['contig_select'].param.watch(on_contig_change, 'value')


    ## Build Variables section - TWO SEPARATE SECTIONS for each view
    variables_title_one = Div(text="<span style='font-size: 1.2em;'><b>Variables</b></span>")
    variables_title_all = Div(text="<span style='font-size: 1.2em;'><b>Variables</b></span>")
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
    combined_features_cbg = None

    # Feature type filter (MultiChoice) - only show if annotation types exist
    # Gene map is plotted if at least one feature type is selected
    feature_type_multichoice = None
    if widgets['annotation_types']:
        # Initially select only CDS if available, otherwise nothing
        initial_value = ["CDS"] if "CDS" in widgets['annotation_types'] else []
        multichoice_stylesheet = InlineStyleSheet(css=":host { background-color: white; }")
        feature_type_multichoice = MultiChoice(
            options=widgets['annotation_types'],
            value=initial_value,
            placeholder="Choose feature types to plot",
            sizing_mode="stretch_width",
            stylesheets=[multichoice_stylesheet]
        )

    # Phage color scheme checkbox - only show if database has CDS features with PHAROKKA functions
    phage_colors_cbg = None
    try:
        cur = conn.cursor()
        result = cur.execute(
            "SELECT Status FROM Constants WHERE Constant = 'pharokka'"
        ).fetchone()
        has_pharokka_functions = result[0] if result else False
        
        if has_pharokka_functions:
            phage_colors_cbg = CheckboxGroup(
                labels=["Use phage color scheme for CDS"],
                active=[]
            )
    except Exception:
        pass  # Constants table might not exist in older databases

    # Plot isoforms checkbox - only show if at least one locus_tag appears more than once
    plot_isoforms_cbg = None
    try:
        cur = conn.cursor()
        result = cur.execute(
            "SELECT Status FROM Constants WHERE Constant = 'isoforms'"
        ).fetchone()
        has_isoforms = result[0] if result else False
        
        if has_isoforms:
            plot_isoforms_cbg = CheckboxGroup(
                labels=["Plot isoforms"],
                active=[]  # Unchecked by default
            )
    except Exception:
        pass  # Constants table or isoforms constant might not exist in older databases

    # Build combined labels: Genome features (without Gene map) + Custom contig features
    combined_labels = []
    if genome_cbg_one is not None:
        combined_labels.extend(genome_cbg_one.labels)  # Already without Gene map
    if widgets['custom_contig_subplots']:
        combined_labels.extend(widgets['custom_contig_subplots'])

    if combined_labels:
        combined_features_cbg = CheckboxButtonGroup(
            labels=combined_labels, active=[],
            sizing_mode="stretch_width", orientation="vertical"
        )

    # Add Genome section to contig_content
    below_contig_children = []
    
    # Create position range inputs
    from_position_input = TextInput(value="0", placeholder="Start position", sizing_mode="stretch_width", margin=(0, 0, 0, 0))
    to_position_input = TextInput(value="", placeholder="End position", sizing_mode="stretch_width", margin=(0, 5, 0, 0))
    
    position_label_from = Div(text="From", width=40, margin=(5, 0, 5, 5))
    position_label_to = Div(text="to", width=25, margin=(5, 0, 5, 5))
    
    position_row = row(
        position_label_from, from_position_input, 
        position_label_to, to_position_input,
        sizing_mode="stretch_width",
        margin=(10, 0, 5, 0)
    )
    
    below_contig_children.append(position_row)

    # Check if sequence data is available in the database
    has_sequence_data = False
    try:
        cur = conn.cursor()
        cur.execute("SELECT 1 FROM information_schema.tables WHERE table_name = 'Contig_sequence'")
        has_sequence_data = cur.fetchone() is not None
    except Exception:
        pass

    sequence_cbg = None
    sequence_row = None
    if has_sequence_data:
        sequence_cbg = CheckboxGroup(labels=["Plot sequence"], active=[])
        sequence_help_tooltip = Tooltip(content="Not recommended for regions larger than 1000 bp as it may slow down rendering", position="right")
        sequence_help_btn = HelpButton(
            tooltip=sequence_help_tooltip, width=20, height=20,
            align="center", button_type="light", stylesheets=[toggle_stylesheet]
        )
        sequence_row = row(sequence_cbg, sequence_help_btn, sizing_mode="stretch_width")

    if feature_type_multichoice is not None or combined_features_cbg is not None:
        # Build simple Genome section header (no collapse needed - Contigs section handles that)
        genome_title = Div(text="<b>Other genomic features to plot:</b>", align="center")

        genome_help_tooltip = widgets['helps_widgets'][genome_index_one] if genome_index_one is not None else None
        if genome_help_tooltip is not None:
            help_btn = HelpButton(tooltip=genome_help_tooltip, width=20, height=20, align="center", button_type="light", stylesheets=[toggle_stylesheet])
            genome_hdr = row(genome_title, help_btn, sizing_mode="stretch_width", align="center")
        else:
            genome_hdr = genome_title

        # Build genome section
        # Layout: feature_type_multichoice, sequence_row (if available), phage_colors_cbg (if available), plot_isoforms_cbg (if available), genome_hdr, combined_features_cbg
        genome_children = []
        if feature_type_multichoice is not None:
            genome_children.append(feature_type_multichoice)
        if sequence_row is not None:
            genome_children.append(sequence_row)
        if phage_colors_cbg is not None:
            genome_children.append(phage_colors_cbg)
        if plot_isoforms_cbg is not None:
            genome_children.append(plot_isoforms_cbg)
        genome_children.append(genome_hdr)
        if combined_features_cbg is not None:
            combined_features_cbg.visible = True
            genome_children.append(combined_features_cbg)
        genome_section = column(*genome_children, visible=True, sizing_mode="stretch_width", margin=(0, 0, 0, 0))

    if genome_section is not None:
        below_contig_children = list(below_contig_children) + [genome_section]

    # Fallback: if no genome section exists, append sequence_row directly
    if sequence_row is not None and genome_section is None:
        below_contig_children.append(sequence_row)

    below_contig_content = column(
        *below_contig_children,
        visible=True, sizing_mode="stretch_width", margin=(0, 0, 0, 0)
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

    ## Plotting parameters section
    separator_plotting_params = Div(text="", height=2, sizing_mode="stretch_width",
        styles={'background-color': '#333', 'margin-top': '10px', 'margin-bottom': '10px'})
    plotting_params_toggle_btn = Button(label="▼", width=20, height=20, button_type="primary", align="center", margin=0, stylesheets=[toggle_stylesheet])
    plotting_params_title = Div(text="<b>Plotting parameters</b>", align="center")
    plotting_params_header = row(plotting_params_toggle_btn, plotting_params_title, sizing_mode="stretch_width", align="center")

    same_y_scale_cbg = CheckboxGroup(labels=["Use same scale for all y axis"], active=[])
    same_y_scale_tooltip = Tooltip(content=(
        "All Samples view: all y axis share the global maximum y value. "
        "One Sample view: all coverage-related y axis share the coverage maximum y value."
    ), position="right")
    same_y_scale_help = HelpButton(tooltip=same_y_scale_tooltip, width=20, height=20,
        align="center", button_type="light", stylesheets=[toggle_stylesheet])
    same_y_scale_row = row(same_y_scale_cbg, same_y_scale_help, sizing_mode="stretch_width")
    
    genemap_height_input = Spinner(value=100, low=10, high=1000, step=10, width=80, margin=(0, 2, 0, 0))
    genemap_height_label = Div(text="Height of gene map (px)", width=160, margin=(5, 0, 5, 5))
    genemap_height_row = row(genemap_height_input, genemap_height_label, sizing_mode="stretch_width", margin=(0, 0, 5, 0))

    subplot_height_input = Spinner(value=100, low=10, high=1000, step=10, width=80, margin=(0, 2, 0, 0))
    subplot_height_label = Div(text="Height per subplot (px)", width=160, margin=(5, 0, 5, 5))
    subplot_height_row = row(subplot_height_input, subplot_height_label, sizing_mode="stretch_width")

    plotting_params_content = pn.Column(
        same_y_scale_row, genemap_height_row, subplot_height_row, 
        sizing_mode="stretch_width"
    )
    plotting_params_toggle_btn.on_click(make_toggle_callback(plotting_params_toggle_btn, plotting_params_content))

    ## Create final Apply and Peruse data buttons
    apply_button = Button(label="APPLY", align="center", stylesheets=[stylesheet], css_classes=["apply-btn"])
    apply_button.on_click(lambda: apply_clicked())

    # Peruse button will be positioned in the plot area, styled to match toolbar
    peruse_button = pn.widgets.Button(
        name="SHOW SUMMARY",
        height=30,
        stylesheets=[stylesheet],
        css_classes=["apply-btn"],
        visible=False  # Hidden until plot is generated
    )
    peruse_button.on_click(lambda event: peruse_clicked())

    # Download contig summary - Panel FileDownload widget
    download_contig_button = pn.widgets.FileDownload(
        callback=make_contig_download_callback,
        filename="contig_summary.csv",
        label="DOWNLOAD CONTIG SUMMARY",
        button_type="primary",
        height=30,
        visible=False
    )
    download_widgets['contig'] = download_contig_button

    # Download metrics summary - Panel FileDownload widget
    download_metrics_button = pn.widgets.FileDownload(
        callback=make_metrics_download_callback,
        filename="metrics_summary.csv",
        label="DOWNLOAD METRICS SUMMARY",
        button_type="primary",
        height=30,
        visible=False
    )
    download_widgets['metrics'] = download_metrics_button

    # Download feature data - Panel FileDownload widget
    download_data_button = pn.widgets.FileDownload(
        callback=make_data_download_callback,
        filename="data.csv",
        label="DOWNLOAD DATA",
        button_type="primary",
        height=30,
        visible=False
    )
    download_widgets['data'] = download_data_button

    # Only Apply button in left panel now
    buttons_row = apply_button


    ## Initialize section titles with counts
    update_section_titles()

    ## Put together all DOM elements
    # Create visual separators (horizontal lines) - using 2px height for consistent rendering at all zoom levels
    separator_filtering = Div(text="", height=2, sizing_mode="stretch_width", styles={'background-color': '#333', 'margin-top': '10px', 'margin-bottom': '10px'})
    separator_samples = Div(text="", height=2, sizing_mode="stretch_width", styles={'background-color': '#333', 'margin-top': '10px', 'margin-bottom': '10px'})
    separator_contigs = Div(text="", height=2, sizing_mode="stretch_width", styles={'background-color': '#333', 'margin-top': '10px', 'margin-bottom': '10px'})
    separator_variables = Div(text="", height=2, sizing_mode="stretch_width", styles={'background-color': '#333', 'margin-top': '10px', 'margin-bottom': '10px'})
    
    # Gene map is now part of the Genome module's CheckboxButtonGroup
    # Build controls list conditionally based on whether samples exist
    # When no samples: hide views toggle, samples section, and variables section
    if widgets['has_samples']:
        controls_children = [logo, views, separator_filtering, filtering_header, filtering_content,
                             separator_contigs, contig_header, above_contig_content, widgets['contig_select'], below_contig_content,
                             separator_samples, sample_title, above_sample_content, widgets['sample_select'],
                             separator_variables,
                             variables_section_one,  # One Sample view (with module checkboxes)
                             variables_section_all,  # All Samples view (title headers only)
                             separator_plotting_params, plotting_params_header, plotting_params_content,
                             buttons_row]
        placeholder_text = "<i>No plot yet. Select one sample, one contig and at least one variable in \"One sample\" mode or one contig and one variable in \"All samples\" mode and click Apply.</i>"
    else:
        # No samples: simplified UI with only Filtering, Contigs, and Apply button
        controls_children = [logo, filtering_header, filtering_content,
                             separator_contigs, contig_header, above_contig_content, widgets['contig_select'], below_contig_content,
                             buttons_row]
        placeholder_text = "<i>No plot yet. Select one contig and click Apply to view the genome annotation.</i>"

    controls_column = pn.Column(*controls_children, sizing_mode="stretch_height", css_classes=["left-col"])

    peruse_button.visible = False  # Initially hidden
    main_placeholder = pn.Column(
        pn.pane.HTML(placeholder_text),
        sizing_mode="stretch_both"
    )

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