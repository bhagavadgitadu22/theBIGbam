import argparse
import os
import time
import duckdb
import traceback

import panel as pn

from bokeh.layouts import column, row
from bokeh.models import Div, InlineStyleSheet, Tooltip
from bokeh.models.widgets import CheckboxGroup, HelpButton, Button, RadioButtonGroup, CheckboxButtonGroup, Select, TextInput, Spinner, MultiChoice, ColorPicker

# Import the plotting function from the repo
from .plotting_data_per_sample import generate_bokeh_plot_per_sample, generate_bokeh_plot_mag_view, DEFAULT_GENEMAP_WINDOW, DEFAULT_SEQUENCE_WINDOW, _DEFAULT_MAX_BASE_RESOLUTION
from .plotting_data_all_samples import generate_bokeh_plot_all_samples
from ..database.database_getters import get_filtering_metadata, resolve_distinct_values, ANNOTATION_EXCLUDED_COLUMNS, is_mag_mode, get_mag_contig_map
from .searchable_select import SearchableSelect

def build_controls(conn):
    """Query DB and return widgets and helper mappings."""
    cur = conn.cursor()

    # Get annotation feature types from Annotated_types table
    cur.execute("SELECT Type_name FROM Annotated_types ORDER BY Frequency DESC")
    annotation_types = [r[0] for r in cur.fetchall()]

    # Detect MAG mode
    has_mags = is_mag_mode(conn)
    mag_to_contigs, contig_to_mag = get_mag_contig_map(conn)
    mags = sorted(mag_to_contigs.keys())

    # Widget Selector for MAGs (visible only in MAG-mode databases)
    mag_select = SearchableSelect(
        value=mags[0] if len(mags) == 1 else "",
        options=mags,
        placeholder="Type to search MAGs...",
        sizing_mode="stretch_width",
        margin=(0, 5, 0, 5),
        visible=has_mags,
    )

    # View toggle (MAG / Contig). Only relevant in MAG-mode databases.
    view_radio = RadioButtonGroup(
        labels=["MAG view", "Contig view"],
        active=1,
        visible=has_mags,
        sizing_mode="stretch_width",
        margin=(0, 5, 10, 5),
    )

    # Widget Selector for Contigs (autocomplete with max 20 suggestions)
    cur.execute("SELECT Contig_name, Contig_length FROM Contig ORDER BY Contig_name")
    rows = cur.fetchall()
    contigs = [r[0] for r in rows]
    contig_lengths = {r[0]: r[1] for r in rows}  # Dictionary mapping contig_name -> length

    # Pre-compute MAG-space offsets: {mag_name: {contig_name: cumulative_offset}}
    # Contigs are in longest-first order (matching get_mag_contigs / mag_to_contigs order).
    mag_to_contig_offsets = {}
    for _mag_name, _contigs in mag_to_contigs.items():
        off = 0
        d = {}
        for c in _contigs:
            d[c] = off
            off += contig_lengths.get(c, 0)
        mag_to_contig_offsets[_mag_name] = d

    # If only one contig in database, pre-fill the field
    contig_select = SearchableSelect(
        value=contigs[0] if len(contigs) == 1 else "",
        options=contigs,
        placeholder="Type to search contigs...",
        sizing_mode="stretch_width",
        margin=(0, 5, 0, 5)
    )

    # Widget Selector for Samples (autocomplete with max 20 suggestions)
    # Sample table only exists when BAM files were provided
    cur.execute("SELECT 1 FROM information_schema.tables WHERE table_name = 'Sample'")
    has_sample_table = cur.fetchone() is not None
    if has_sample_table:
        cur.execute("SELECT Sample_name FROM Sample ORDER BY Sample_name")
        samples = [r[0] for r in cur.fetchall()]
    else:
        samples = []
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
    sample_to_contigs = {}
    contig_to_samples = {}
    if has_sample_table:
        cur.execute("""
        SELECT Contig.Contig_name, Sample.Sample_name FROM Coverage
          JOIN Contig ON Coverage.Contig_id = Contig.Contig_id
          JOIN Sample ON Coverage.Sample_id = Sample.Sample_id
        """)
        for contig_name, sample_name in cur.fetchall():
            sample_to_contigs.setdefault(sample_name, set()).add(contig_name)
            contig_to_samples.setdefault(contig_name, set()).add(sample_name)

    # Pre-compute MAG→samples mapping: a sample belongs to a MAG if any member contig has coverage
    mag_to_samples = {}
    for _mag_name, _mag_contigs in mag_to_contigs.items():
        _s = set()
        for _c in _mag_contigs:
            _s |= contig_to_samples.get(_c, set())
        mag_to_samples[_mag_name] = _s

    # Get variables that have data (their feature table exists)
    cur.execute("SELECT DISTINCT Feature_table_name FROM Variable")
    tables_with_data = [r[0] for r in cur.fetchall()]

    # Get modules that have at least one variable with data
    # Define the display order for modules
    MODULE_ORDER = ["Genome", "Coverage", "Misalignment", "RNA", "Long-reads", "Paired-reads", "Phage termini"]

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
        'mag_select': mag_select,
        'view_radio': view_radio,
        'sample_to_contigs': sample_to_contigs,
        'contig_to_samples': contig_to_samples,
        'mag_to_contigs': mag_to_contigs,
        'contig_to_mag': contig_to_mag,
        'mag_to_samples': mag_to_samples,
        'mag_to_contig_offsets': mag_to_contig_offsets,
        'mags': mags,
        'has_mags': has_mags,
        'module_names': module_names,
        'module_widgets_one': module_widgets_one,
        'helps_widgets': helps_widgets,
        'variables_widgets_one': variables_widgets_one,
        'variables_widgets_all': variables_widgets_all,
        'contigs': contigs,
        'contig_lengths': contig_lengths,
        'samples': samples,
        'custom_contig_subplots': custom_contig_subplots,
        'annotation_types': annotation_types,
        'has_samples': has_samples  # True if database has any samples
    }
    return widgets

def create_layout(db_path, enable_timing=False):
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
            'MAG': 'MAG',
            'Contig': 'Contig',
            'Annotations': 'Contig_annotation',
            'Sample': 'Sample',
            'MAG Coverage': 'Explicit_coverage_per_MAG',
            'MAG misassembly': 'Explicit_misassembly_per_MAG',
            'MAG microdiversity': 'Explicit_microdiversity_per_MAG',
            'Coverage': 'Explicit_coverage',
            'Misassembly': 'Explicit_misassembly',
            'Microdiversity': 'Explicit_microdiversity',
            'Side misassembly': 'Explicit_side_misassembly',
            'Topology': 'Explicit_topology',
            'Termini': 'Explicit_phage_mechanisms',
        }
        # Categories that yield (MAG_name, Sample_name) pairs (or MAG-only) —
        # expanded to (Contig_name, Sample_name) via MAG_contigs_association.
        mag_categories = {'MAG', 'MAG Coverage', 'MAG misassembly', 'MAG microdiversity'}

        def get_pairs_for_condition(category, column_name, operator, value):
            """Query database for contig/sample pairs matching a single condition."""
            source = source_table_map.get(category)
            if not source:
                return set()

            # Translate "has"/"has not" to SQL LIKE/NOT LIKE
            if operator == "has":
                operator = "LIKE"
                value = f"%{value}%"
            elif operator == "has not":
                operator = "NOT LIKE"
                value = f"%{value}%"

            # Check if this column is from Contig_annotation / Annotation_qualifier / Contig_qualifier
            col_info = filtering_metadata.get(category, {}).get('columns', {}).get(column_name, {})
            col_source = col_info.get('source')
            qualifier_key = col_info.get('qualifier_key')
            params = [value]
            has_samples = widgets['has_samples']

            if col_source == 'Contig_annotation':
                # Pivoted column from the Contig_annotation view — join via Contig_annotation_core
                if has_samples:
                    query = f'''
                        SELECT DISTINCT c.Contig_name, s.Sample_name
                        FROM Contig_annotation ca
                        JOIN Contig c ON ca.Contig_id = c.Contig_id
                        LEFT JOIN Coverage p ON c.Contig_id = p.Contig_id
                        LEFT JOIN Sample s ON p.Sample_id = s.Sample_id
                        WHERE ca."{column_name}" {operator} ?
                    '''
                else:
                    query = f'''
                        SELECT DISTINCT c.Contig_name, NULL
                        FROM Contig_annotation ca
                        JOIN Contig c ON ca.Contig_id = c.Contig_id
                        WHERE ca."{column_name}" {operator} ?
                    '''
            elif col_source == 'Annotation_qualifier':
                # KV-lookup — filter by Key + Value against Annotation_qualifier
                params = [qualifier_key, value]
                if has_samples:
                    query = f'''
                        SELECT DISTINCT c.Contig_name, s.Sample_name
                        FROM Annotation_qualifier aq
                        JOIN Contig_annotation_core cac ON aq.Annotation_id = cac.Annotation_id
                        JOIN Contig c ON cac.Contig_id = c.Contig_id
                        LEFT JOIN Coverage p ON c.Contig_id = p.Contig_id
                        LEFT JOIN Sample s ON p.Sample_id = s.Sample_id
                        WHERE aq."Key" = ? AND aq."Value" {operator} ?
                    '''
                else:
                    query = f'''
                        SELECT DISTINCT c.Contig_name, NULL
                        FROM Annotation_qualifier aq
                        JOIN Contig_annotation_core cac ON aq.Annotation_id = cac.Annotation_id
                        JOIN Contig c ON cac.Contig_id = c.Contig_id
                        WHERE aq."Key" = ? AND aq."Value" {operator} ?
                    '''
            elif col_source == 'Contig_qualifier':
                # KV-lookup — filter by Key + Value against Contig_qualifier (contig-level)
                params = [qualifier_key, value]
                if has_samples:
                    query = f'''
                        SELECT DISTINCT c.Contig_name, s.Sample_name
                        FROM Contig_qualifier cq
                        JOIN Contig c ON cq.Contig_id = c.Contig_id
                        LEFT JOIN Coverage p ON c.Contig_id = p.Contig_id
                        LEFT JOIN Sample s ON p.Sample_id = s.Sample_id
                        WHERE cq."Key" = ? AND cq."Value" {operator} ?
                    '''
                else:
                    query = f'''
                        SELECT DISTINCT c.Contig_name, NULL
                        FROM Contig_qualifier cq
                        JOIN Contig c ON cq.Contig_id = c.Contig_id
                        WHERE cq."Key" = ? AND cq."Value" {operator} ?
                    '''
            elif category == 'Contig' and has_samples:
                # Contig table has no Sample_name, left-join to preserve contigs with 0 samples
                query = f'''
                    SELECT DISTINCT c.Contig_name, s.Sample_name
                    FROM Contig c
                    LEFT JOIN Coverage p ON c.Contig_id = p.Contig_id
                    LEFT JOIN Sample s ON p.Sample_id = s.Sample_id
                    WHERE c."{column_name}" {operator} ?
                '''
            elif category == 'Contig':
                # Contig table - no samples available
                query = f'''
                    SELECT DISTINCT c.Contig_name, NULL
                    FROM Contig c
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
            elif category == 'MAG':
                # MAG table — expand to all member contigs × samples present on them
                if has_samples:
                    query = f'''
                        SELECT DISTINCT c.Contig_name, s.Sample_name
                        FROM MAG mg
                        JOIN MAG_contigs_association mca ON mca.MAG_id = mg.MAG_id
                        JOIN Contig c ON c.Contig_id = mca.Contig_id
                        LEFT JOIN Coverage p ON c.Contig_id = p.Contig_id
                        LEFT JOIN Sample s ON p.Sample_id = s.Sample_id
                        WHERE mg."{column_name}" {operator} ?
                    '''
                else:
                    query = f'''
                        SELECT DISTINCT c.Contig_name, NULL
                        FROM MAG mg
                        JOIN MAG_contigs_association mca ON mca.MAG_id = mg.MAG_id
                        JOIN Contig c ON c.Contig_id = mca.Contig_id
                        WHERE mg."{column_name}" {operator} ?
                    '''
            elif category in mag_categories:
                # Explicit_*_per_MAG views — row grain is (MAG_name, Sample_name);
                # expand to member contigs via MAG_contigs_association.
                query = f'''
                    SELECT DISTINCT c.Contig_name, v.Sample_name
                    FROM {source} v
                    JOIN MAG mg ON mg.MAG_name = v.MAG_name
                    JOIN MAG_contigs_association mca ON mca.MAG_id = mg.MAG_id
                    JOIN Contig c ON c.Contig_id = mca.Contig_id
                    WHERE v."{column_name}" {operator} ?
                '''
            else:
                # Tables with both Contig_name and Sample_name (views)
                query = f'''
                    SELECT DISTINCT Contig_name, Sample_name
                    FROM {source}
                    WHERE "{column_name}" {operator} ?
                '''

            # Coerce value type to match column type (prevents DuckDB implicit cast errors)
            col_type = col_info.get('type', 'numeric')
            is_bool = col_info.get('is_bool', False)
            if is_bool:
                # Convert yes/no back to boolean for SQL
                value = value.lower() == 'yes' if isinstance(value, str) else bool(value)
            elif col_type == 'text' and not isinstance(value, str):
                value = str(value)
            elif col_type == 'numeric' and isinstance(value, str):
                try:
                    value = float(value)
                except (ValueError, TypeError):
                    return set()

            # value is always the last param; refresh it after coercion
            params[-1] = value

            try:
                cur.execute(query, params)
                return {(row[0], row[1]) for row in cur.fetchall()}
            except duckdb.Error as e:
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
        if widgets['has_mags']:
            # MAG mode: contig list is a child of the selected MAG, not of the sample.
            sel_mag = widgets['mag_select'].value
            if sel_mag and sel_mag in widgets['mag_to_contigs']:
                completions = list(widgets['mag_to_contigs'][sel_mag])
            else:
                completions = list(orig_contigs)

            # Apply Filtering query builder filters (MAG mode)
            filtered_pairs = get_filtering_filtered_pairs()
            if filtered_pairs is not None:
                if views.active == 0 and widgets['sample_select'].value:
                    sel_sample = widgets['sample_select'].value
                    allowed_contigs = {pair[0] for pair in filtered_pairs if pair[1] == sel_sample}
                else:
                    allowed_contigs = {pair[0] for pair in filtered_pairs}
                completions = [c for c in completions if c in allowed_contigs]
        else:
            # Contig mode: filter contigs by selected sample (ONE SAMPLE view).
            if views.active == 0 and widgets['sample_select'].value:
                sel_sample = widgets['sample_select'].value
                allowed = widgets['sample_to_contigs'].get(sel_sample, set())
                completions = [c for c in orig_contigs if c in allowed]
            else:
                completions = list(orig_contigs)

            # Apply Filtering query builder filters (contig mode)
            filtered_pairs = get_filtering_filtered_pairs()
            if filtered_pairs is not None:
                if views.active == 0 and widgets['sample_select'].value:
                    sel_sample = widgets['sample_select'].value
                    allowed_contigs = {pair[0] for pair in filtered_pairs if pair[1] == sel_sample}
                else:
                    allowed_contigs = {pair[0] for pair in filtered_pairs}
                completions = [c for c in completions if c in allowed_contigs]

        update_widget_completions(widgets['contig_select'], completions)

    def refresh_sample_options_unlocked():
        """Core logic — does NOT check the lock."""
        if widgets['has_mags']:
            # MAG mode: filter samples by selected MAG.
            sel_mag = widgets['mag_select'].value
            if views.active == 0 and sel_mag and sel_mag in widgets['mag_to_samples']:
                allowed = widgets['mag_to_samples'][sel_mag]
                completions = [s for s in orig_samples if s in allowed]
            else:
                completions = list(orig_samples)

            # Apply Filtering query builder filters (MAG mode)
            filtered_pairs = get_filtering_filtered_pairs()
            if filtered_pairs is not None:
                if views.active == 0 and widgets['contig_select'].value:
                    sel_contig = widgets['contig_select'].value
                    allowed_samples = {pair[1] for pair in filtered_pairs if pair[0] == sel_contig}
                else:
                    allowed_samples = {pair[1] for pair in filtered_pairs}
                completions = [s for s in completions if s in allowed_samples]
        else:
            # Contig mode: filter samples by selected contig (ONE SAMPLE view).
            if views.active == 0 and widgets['contig_select'].value:
                sel_contig = widgets['contig_select'].value
                allowed = widgets['contig_to_samples'].get(sel_contig, set())
                completions = [s for s in orig_samples if s in allowed]
            else:
                completions = list(orig_samples)

            # Apply Filtering query builder filters (contig mode)
            filtered_pairs = get_filtering_filtered_pairs()
            if filtered_pairs is not None:
                if views.active == 0 and widgets['contig_select'].value:
                    sel_contig = widgets['contig_select'].value
                    allowed_samples = {pair[1] for pair in filtered_pairs if pair[0] == sel_contig}
                else:
                    allowed_samples = {pair[1] for pair in filtered_pairs}
                completions = [s for s in completions if s in allowed_samples]
        update_widget_completions(widgets['sample_select'], completions)

    def refresh_mag_options_unlocked():
        """Recompute MAG dropdown based on selected sample (ONE SAMPLE view).
        A MAG is valid iff at least one of its member contigs has coverage for that sample.
        Only meaningful in MAG-mode DBs; no-op otherwise.
        """
        if not widgets['has_mags']:
            return
        mag_to_contigs = widgets['mag_to_contigs']
        sel_sample = widgets['sample_select'].value
        if views.active == 0 and sel_sample:
            valid_contigs = widgets['sample_to_contigs'].get(sel_sample, set())
            completions = [
                m for m in sorted(mag_to_contigs.keys())
                if any(c in valid_contigs for c in mag_to_contigs[m])
            ]
        else:
            completions = sorted(mag_to_contigs.keys())

        # Apply Filtering query builder filters: keep a MAG only if at least
        # one of its member contigs survives the filter.
        filtered_pairs = get_filtering_filtered_pairs()
        if filtered_pairs is not None:
            if views.active == 0 and sel_sample:
                allowed_contigs = {pair[0] for pair in filtered_pairs if pair[1] == sel_sample}
            else:
                allowed_contigs = {pair[0] for pair in filtered_pairs}
            completions = [
                m for m in completions
                if any(c in allowed_contigs for c in mag_to_contigs[m])
            ]
        update_widget_completions(widgets['mag_select'], completions)

    def update_section_titles():
        """Update Filtering, Contigs, Samples, and MAGs section titles with current counts."""
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
        if widgets['has_mags']:
            mags_count = len(set(widgets['mag_select'].options) - {""})
            mag_title.text = f"<span style='font-size: 1.2em;'><b>MAGs</b></span> ({mags_count} available)"

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
            from bokeh.io import curdoc
            doc = curdoc()
            doc.hold('combine')
            global_toggle_lock['locked'] = True
            try:
                for other in widgets['variables_widgets_all']:
                    if other is cbg:
                        target = [] if sel_index is None else [sel_index]
                        if other.active != target:
                            other.active = target
                    else:
                        if other.active:
                            other.active = []
            finally:
                global_toggle_lock['locked'] = False
                doc.unhold()
        return callback

    # Views (One sample / All samples) callback: show/hide sample-related controls
    def on_view_change(attr, old, new):
        is_all = (new == 1)  # True means All samples

        from bokeh.io import curdoc
        doc = curdoc()
        doc.hold('combine')
        # Lock callbacks during view change to prevent cascading updates
        global_toggle_lock['locked'] = True
        try:
            # Toggle Sample section - hide entirely in All Samples view
            separator_samples.visible = not is_all
            sample_title.visible = not is_all
            above_sample_content.visible = not is_all
            widgets['sample_select'].visible = not is_all

            # Toggle visibility between the two variables sections
            # Each section maintains its own state independently
            variables_section_one.visible = not is_all
            variables_section_all.visible = is_all
            sample_params_header.visible = is_all

            # Refresh options while still locked (suppresses cascading callbacks)
            # Don't invalidate filtering cache - filtering is shared between views and hasn't changed
            refresh_contig_options_unlocked()
            if not is_all:
                refresh_sample_options_unlocked()
            refresh_mag_options_unlocked()
        finally:
            global_toggle_lock['locked'] = False
        update_section_titles()
        doc.unhold()

    ## Apply button function
    def apply_clicked():
        try:
            if enable_timing:
                t_apply_start = time.perf_counter()
            contig = widgets['contig_select'].value
            has_samples = widgets['has_samples']
            
            # When no samples exist, treat as "One Sample" mode with no sample/variables
            is_all = (views.active == 1) if has_samples else False
            sample = widgets['sample_select'].value if has_samples else None

            # Genome module is shared between views (in Contigs section)
            # Gene map is shown if at least one feature type is selected in the multichoice
            selected_feature_types = feature_type_multichoice.value if feature_type_multichoice is not None else None
            genbank_path = db_path if (selected_feature_types and len(selected_feature_types) > 0) else None
            plot_isoforms = (0 in plot_isoforms_cbg.active) if (plot_isoforms_cbg is not None and genbank_path) else True
            feature_label_key = feature_label_select.value if feature_label_select is not None else None

            # Read custom color rules from color rows
            OP_TO_MODE = {
                '=': 'exact',
                '!=': 'not_equal',
                'has': 'has',
                'has not': 'has_not',
                '<': 'lt',
                '>': 'gt',
                'Use random colors': 'random',
            }
            custom_colors = []
            for row_data in custom_color_rows:
                key = row_data['qualifier_select'].value
                if not key:
                    continue
                operator = row_data['operator_select'].value
                mode = OP_TO_MODE.get(operator, 'exact')
                if mode == 'random':
                    custom_colors.append({'qualifier_key': key, 'match_mode': 'random'})
                    continue
                widget = row_data['input_ref']['widget']
                raw_val = widget.value
                color = row_data['color_picker'].color
                if raw_val is None or raw_val == "" or not color:
                    continue
                # Spinner values come through as int/float; keep as number so
                # downstream lt/gt comparisons don't have to re-parse.
                val = float(raw_val) if isinstance(raw_val, (int, float)) else raw_val
                custom_colors.append({
                    'qualifier_key': key, 'value': val, 'color': color,
                    'match_mode': mode,
                })

            # Select the correct widget set based on current view
            active_variables_widgets = widgets['variables_widgets_all'] if is_all else widgets['variables_widgets_one']

            # Read display/sizing parameters (needed by both MAG view and Contig view paths)
            max_genemap_window = int(max_genemap_window_input.value)
            same_y_scale = (0 in same_y_scale_cbg.active)
            subplot_size = int(subplot_height_input.value)
            genemap_size = int(genemap_height_input.value)
            sequence_size = int(sequence_height_input.value)
            translated_sequence_size = int(translated_sequence_height_input.value)
            max_binning = int(max_binning_window_input.value)
            min_coverage_freq = float(min_coverage_freq_input.value)

            # Parse and validate position inputs
            xstart = None
            xend = None

            # --- MAG view early path (runs instead of the contig-based path below) ---
            is_mag_view = widgets['has_mags'] and widgets['view_radio'].active == 0
            if is_mag_view:
                active_mag = widgets['mag_select'].value
                if not active_mag:
                    peruse_button.visible = False
                    main_placeholder.objects = [pn.pane.HTML("<pre>Error: Please select a MAG.</pre>")]
                    return

                # Parse position inputs (None → use full MAG extent inside plotting function)
                try:
                    xstart = int(from_position_input.value) if from_position_input.value.strip() else 1
                    xend = int(to_position_input.value) if to_position_input.value.strip() else None
                except ValueError:
                    peruse_button.visible = False
                    main_placeholder.objects = [pn.pane.HTML("<pre>Error: Invalid position range - positions must be integers.</pre>")]
                    return

                if xend is not None and xstart >= xend:
                    peruse_button.visible = False
                    main_placeholder.objects = [pn.pane.HTML("<pre>Error: Invalid position range - start must be less than end.</pre>")]
                    return

                mag_window = (xend - xstart) if xend is not None else float('inf')
                plot_genemap = mag_window <= max_genemap_window

                # Sequence / translated-sequence visibility — same threshold logic
                # as the contig-view path (see below).
                mag_max_sequence_window = int(max_sequence_window_input.value)
                mag_plot_sequence = False
                if sequence_cbg is not None and 0 in sequence_cbg.active:
                    if mag_window <= mag_max_sequence_window:
                        mag_plot_sequence = True
                    else:
                        print(f"Warning: Sequence will not be plotted for regions larger than {mag_max_sequence_window} bp.", flush=True)
                mag_plot_translated_sequence = False
                if translated_sequence_cbg is not None and 0 in translated_sequence_cbg.active:
                    if mag_window <= mag_max_sequence_window:
                        mag_plot_translated_sequence = True
                    else:
                        print(f"Warning: Translated sequence will not be plotted for regions larger than {mag_max_sequence_window} bp.", flush=True)

                # Collect requested features
                if is_all:
                    # ALL SAMPLES mode: pick one sample-level variable (same logic as contig-view all-samples)
                    selected_var = None
                    for cbg in active_variables_widgets:
                        if cbg.active and selected_var is None:
                            selected_var = cbg.labels[cbg.active[-1]]
                    mag_requested_features = []
                    if combined_features_cbg is not None:
                        for idx in combined_features_cbg.active:
                            mag_requested_features.append(combined_features_cbg.labels[idx])
                    if selected_var:
                        mag_requested_features.append(selected_var)
                    mag_requested_features = [f for f in mag_requested_features if f != "Gene map"]
                    if not mag_requested_features:
                        raise ValueError("In 'All samples' view you must select at least one variable.")
                    # Compute filtered samples for this MAG (union across all its contigs)
                    mag_contigs = widgets['mag_to_contigs'].get(active_mag, [])
                    filtered_samples = [s for s in orig_samples
                                        if any(s in widgets['contig_to_samples'].get(c, set()) for c in mag_contigs)]
                    filtering_pairs = get_filtering_filtered_pairs()
                    if filtering_pairs is not None:
                        allowed_s = {pair[1] for pair in filtering_pairs}
                        filtered_samples = [s for s in filtered_samples if s in allowed_s]
                    mag_allowed_samples = set(filtered_samples)
                else:
                    # ONE SAMPLE mode: collect all selected features (same logic as the One-sample Contig path below)
                    mag_requested_features = []
                    if combined_features_cbg is not None:
                        for idx in combined_features_cbg.active:
                            mag_requested_features.append(combined_features_cbg.labels[idx])
                    for cbg in active_variables_widgets:
                        for idx in cbg.active:
                            mag_requested_features.append(cbg.labels[idx])
                    # Gene map is handled via genbank_path, not as a feature name
                    mag_requested_features = [f for f in mag_requested_features if f != "Gene map"]
                    mag_allowed_samples = None

                print(f"[start_bokeh_server] MAG view: mag={active_mag}, is_all={is_all}, sample={sample}, features={mag_requested_features}", flush=True)
                if enable_timing:
                    t_params = time.perf_counter()
                    print(f"[timing] Parameter parsing: {t_params - t_apply_start:.3f}s", flush=True)

                # Preserve x-range when re-plotting same MAG/sample/range
                mag_preserve_xrange = (
                    current_plot_state['shared_xrange'] is not None
                    and current_plot_state['contig'] == active_mag
                    and current_plot_state['is_all'] == is_all
                    and (is_all or current_plot_state['sample'] == sample)
                    and current_plot_state['data_xstart'] == xstart
                    and current_plot_state['data_xend'] == xend
                )
                if mag_preserve_xrange:
                    mag_prev_xstart = current_plot_state['shared_xrange'].start
                    mag_prev_xend = current_plot_state['shared_xrange'].end

                if enable_timing:
                    t_plot = time.perf_counter()
                grid = generate_bokeh_plot_mag_view(
                    conn, mag_requested_features, active_mag, sample,
                    xstart=xstart, xend=xend,
                    genbank_path=genbank_path if plot_genemap else None,
                    feature_types=selected_feature_types, plot_isoforms=plot_isoforms,
                    plot_sequence=mag_plot_sequence,
                    plot_translated_sequence=mag_plot_translated_sequence,
                    same_y_scale=same_y_scale, subplot_size=subplot_size,
                    genemap_size=genemap_size, sequence_size=sequence_size,
                    translated_sequence_size=translated_sequence_size,
                    max_base_resolution=max_binning,
                    max_genemap_window=max_genemap_window,
                    max_sequence_window=mag_max_sequence_window,
                    min_relative_value=min_coverage_freq,
                    feature_label_key=feature_label_key,
                    custom_colors=custom_colors if custom_colors else None,
                    is_all=is_all,
                    allowed_samples=mag_allowed_samples,
                )
                if enable_timing:
                    print(f"[timing] generate_bokeh_plot_mag_view (DB queries + plotting): {time.perf_counter() - t_plot:.3f}s", flush=True)

                new_xrange = _get_shared_xrange(grid)
                if mag_preserve_xrange and new_xrange is not None:
                    new_xrange.start = mag_prev_xstart
                    new_xrange.end = mag_prev_xend

                current_plot_state['contig'] = active_mag
                current_plot_state['sample'] = sample
                current_plot_state['is_all'] = is_all
                current_plot_state['shared_xrange'] = new_xrange
                current_plot_state['data_xstart'] = xstart
                current_plot_state['data_xend'] = xend

                if new_xrange is not None:
                    def _sync_from_mag(attr, old, new_val):
                        from_position_input.value = str(int(new_val))
                    def _sync_to_mag(attr, old, new_val):
                        to_position_input.value = str(int(new_val))
                    new_xrange.on_change('start', _sync_from_mag)
                    new_xrange.on_change('end', _sync_to_mag)

                toolbar_row = pn.Row(
                    pn.Spacer(sizing_mode="stretch_width"),
                    peruse_button, download_mag_metrics_button, download_data_button,
                    margin=(0, 0, 5, 0)
                )
                peruse_button.visible = True
                download_mag_metrics_button.visible = bool(has_samples)
                download_data_button.visible = True
                command_hint_pane.visible = False
                main_placeholder.objects = [pn.Column(toolbar_row, command_hint_pane, grid, sizing_mode="stretch_both")]
                if enable_timing:
                    print(f"[timing] Total APPLY (MAG view): {time.perf_counter() - t_apply_start:.3f}s", flush=True)
                return
            # --- end of MAG view early path ---

            # Validate contig is selected
            if not contig:
                peruse_button.visible = False
                main_placeholder.objects = [pn.pane.HTML("<pre>Error: Please select a contig.</pre>")]
                return
            
            # Get contig length for validation
            contig_length = widgets['contig_lengths'].get(contig, 0)
            
            # Parse position inputs
            try:
                xstart = int(from_position_input.value) if from_position_input.value.strip() else 1
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

            # Only plot genome map if window is <= threshold from spinner
            plot_genemap = True
            if (xend - xstart) > max_genemap_window:
                plot_genemap = False
                print(f"Warning: Genome map will not be plotted for regions larger than {max_genemap_window} bp.", flush=True)

            # Only plot sequence if window is <= threshold from spinner
            plot_sequence = False
            if sequence_cbg is not None and 0 in sequence_cbg.active:
                max_seq_window = int(max_sequence_window_input.value)
                if (xend - xstart) <= max_seq_window:
                    plot_sequence = True
                else:
                    print(f"Warning: Sequence will not be plotted for regions larger than {max_seq_window} bp.", flush=True)

            # Only plot translated sequence if window is <= threshold from spinner
            plot_translated_sequence = False
            if translated_sequence_cbg is not None and 0 in translated_sequence_cbg.active:
                max_seq_window = int(max_sequence_window_input.value)
                if (xend - xstart) <= max_seq_window:
                    plot_translated_sequence = True
                else:
                    print(f"Warning: Translated sequence will not be plotted for regions larger than {max_seq_window} bp.", flush=True)

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

            if enable_timing:
                t_params = time.perf_counter()
                print(f"[timing] Parameter parsing: {t_params - t_apply_start:.3f}s", flush=True)

            grid = None
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

                # Get selected ordering column (map UI label "Sample name" to DB column "Sample_name")
                order_by = "Sample_name" if sample_order_select.value == "Sample name" else sample_order_select.value

                print(f"[start_bokeh_server] Generating plot for all samples with variable={selected_var}, contig={contig}, genome_features={genome_features}, filtered_samples={len(filtered_samples)}")
                # Pass plot_genemap to plotting function if supported, else filter genome_features
                if not plot_genemap and genome_features:
                    # Remove "Gene map" from genome_features if present
                    genome_features = [f for f in genome_features if f != "Gene map"]
                if enable_timing:
                    t_plot = time.perf_counter()
                grid = generate_bokeh_plot_all_samples(
                    conn, selected_var, contig, xstart=xstart, xend=xend, genbank_path=genbank_path,
                    genome_features=genome_features if genome_features else None, allowed_samples=set(filtered_samples),
                    feature_types=selected_feature_types, plot_isoforms=plot_isoforms, plot_sequence=plot_sequence,
                    plot_translated_sequence=plot_translated_sequence, same_y_scale=same_y_scale, subplot_size=subplot_size, genemap_size=genemap_size,
                    sequence_size=sequence_size, translated_sequence_size=translated_sequence_size, order_by_column=order_by, max_base_resolution=max_binning,
                    max_genemap_window=max_genemap_window, min_relative_value=min_coverage_freq,
                    feature_label_key=feature_label_key,
                    custom_colors=custom_colors if custom_colors else None
                )
                if enable_timing:
                    print(f"[timing] generate_bokeh_plot_all_samples (DB queries + plotting): {time.perf_counter() - t_plot:.3f}s", flush=True)
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
                # Remove "Gene map" from requested_features if window too large
                if not plot_genemap and requested_features:
                    requested_features = [f for f in requested_features if f != "Gene map"]
                # MAG view is handled by the early path above; MAG track is never shown in Contig view.
                active_mag = None
                if enable_timing:
                    t_plot = time.perf_counter()
                grid = generate_bokeh_plot_per_sample(
                    conn, requested_features, contig, sample, xstart=xstart, xend=xend, genbank_path=genbank_path,
                    feature_types=selected_feature_types, plot_isoforms=plot_isoforms,
                    plot_sequence=plot_sequence, plot_translated_sequence=plot_translated_sequence,
                    same_y_scale=False, subplot_size=subplot_size, genemap_size=genemap_size,
                    sequence_size=sequence_size, translated_sequence_size=translated_sequence_size, max_base_resolution=max_binning,
                    max_genemap_window=int(max_genemap_window_input.value),
                    max_sequence_window=int(max_sequence_window_input.value),
                    min_relative_value=min_coverage_freq,
                    feature_label_key=feature_label_key,
                    custom_colors=custom_colors if custom_colors else None,
                    mag_name=active_mag
                )
                if enable_timing:
                    print(f"[timing] generate_bokeh_plot_per_sample (DB queries + plotting): {time.perf_counter() - t_plot:.3f}s", flush=True)

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
                download_mag_metrics_button,
                download_metrics_button,
                download_data_button,
                margin=(0, 0, 5, 0)
            )
            # Show buttons when plot exists (but some are hidden in 0-sample mode)
            peruse_button.visible = True  # Always show — contig summary available even without samples
            download_data_button.visible = True  # Always show — contig features available even without samples
            download_metrics_button.visible = bool(has_samples)
            download_mag_metrics_button.visible = bool(has_samples and widgets['has_mags'])

            # Hide any previous command hint when plot is refreshed
            command_hint_pane.visible = False

            # Display the plot
            main_placeholder.objects = [pn.Column(toolbar_row, command_hint_pane, grid, sizing_mode="stretch_both")]
            if enable_timing:
                view_name = "all samples" if is_all else "one sample"
                print(f"[timing] Total APPLY ({view_name}): {time.perf_counter() - t_apply_start:.3f}s", flush=True)

        except Exception as e:
            peruse_button.visible = False
            download_metrics_button.visible = False
            download_mag_metrics_button.visible = False
            download_data_button.visible = False
            tb = traceback.format_exc()
            print(f"[start_bokeh_server] Exception: {tb}", flush=True)
            main_placeholder.objects = [pn.pane.HTML(f"<pre>Error building plot:\n{tb}</pre>")]

    ## Peruse button callback function
    def peruse_clicked():
        """Generate and open summary tables in a new browser window."""
        from .perusing_data import generate_and_open_peruse_html
        if enable_timing:
            t_peruse = time.perf_counter()

        is_mag_view = widgets['has_mags'] and widgets['view_radio'].active == 0

        if is_mag_view:
            # MAG view: show MAG characs + sample characs + MAG metrics
            mag = widgets['mag_select'].value
            if not mag:
                print("[start_bokeh_server] Peruse: No MAG selected", flush=True)
                return
            sample = widgets['sample_select'].value
            sample_names = [sample] if sample else []
            generate_and_open_peruse_html(conn, None, sample_names, mag_name=mag, is_mag_view=True)
            if enable_timing:
                print(f"[timing] SHOW SUMMARY (MAG view): {time.perf_counter() - t_peruse:.3f}s", flush=True)
            return

        # Contig view
        contig = widgets['contig_select'].value
        if not contig:
            print("[start_bokeh_server] Peruse: No contig selected", flush=True)
            return

        parent_mag = widgets['contig_to_mag'].get(contig)

        has_samples = widgets['has_samples']
        if not has_samples:
            sample_names = []
        else:
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
        generate_and_open_peruse_html(conn, contig, sample_names, mag_name=parent_mag, is_mag_view=False)
        if enable_timing:
            print(f"[timing] SHOW SUMMARY (contig view, {len(sample_names)} samples): {time.perf_counter() - t_peruse:.3f}s", flush=True)

    ## Download functionality using Panel FileDownload widgets
    import io

    # Store references to download widgets (created later, after widgets dict exists)
    download_widgets = {'contig_metrics': None, 'mag_metrics': None, 'data': None}

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

    def make_contig_metrics_download_callback():
        """Create callback for contig metrics download (contig-level metrics across samples)."""
        from .downloading_data import download_contig_metrics_csv
        if enable_timing:
            t_dl = time.perf_counter()

        contig = widgets['contig_select'].value
        if not contig:
            print("[start_bokeh_server] Download contig metrics: No contig selected", flush=True)
            return io.StringIO("")

        is_all = (views.active == 1)

        if is_all:
            filtered_samples = [s for s in orig_samples if s in widgets['contig_to_samples'].get(contig, set())]
            filtering_pairs = get_filtering_filtered_pairs()
            if filtering_pairs is not None:
                allowed_samples = {pair[1] for pair in filtering_pairs}
                filtered_samples = [s for s in filtered_samples if s in allowed_samples]

            if not filtered_samples:
                print("[start_bokeh_server] Download contig metrics: No samples match filters", flush=True)
                return io.StringIO("")

            sample_names = filtered_samples
            safe_contig = "".join(c if c.isalnum() or c in "-_" else "_" for c in contig)
            if download_widgets['contig_metrics']:
                download_widgets['contig_metrics'].filename = f"{safe_contig}_in_all_samples_contig_metrics.csv"
        else:
            sample = widgets['sample_select'].value
            if not sample:
                print("[start_bokeh_server] Download contig metrics: No sample selected", flush=True)
                return io.StringIO("")
            sample_names = [sample]
            safe_contig = "".join(c if c.isalnum() or c in "-_" else "_" for c in contig)
            safe_sample = "".join(c if c.isalnum() or c in "-_" else "_" for c in sample)
            if download_widgets['contig_metrics']:
                download_widgets['contig_metrics'].filename = f"{safe_contig}_in_{safe_sample}_contig_metrics.csv"

        csv_content = download_contig_metrics_csv(db_path, contig, sample_names)
        if enable_timing:
            print(f"[timing] DOWNLOAD CONTIG METRICS ({len(sample_names)} samples): {time.perf_counter() - t_dl:.3f}s", flush=True)
        if csv_content:
            return io.StringIO(csv_content)
        return io.StringIO("")

    def make_mag_metrics_download_callback():
        """Create callback for MAG metrics download (MAG-level metrics from Explicit_*_per_MAG views)."""
        from .downloading_data import download_mag_metrics_csv
        if enable_timing:
            t_dl = time.perf_counter()

        is_mag_view = widgets['has_mags'] and widgets['view_radio'].active == 0
        if is_mag_view:
            mag = widgets['mag_select'].value
        else:
            contig = widgets['contig_select'].value
            mag = widgets['contig_to_mag'].get(contig) if contig else None

        if not mag:
            print("[start_bokeh_server] Download MAG metrics: No MAG selected/found", flush=True)
            return io.StringIO("")

        sample = widgets['sample_select'].value
        sample_names = [sample] if sample else []

        csv_content = download_mag_metrics_csv(db_path, mag, sample_names)
        if enable_timing:
            print(f"[timing] DOWNLOAD MAG METRICS: {time.perf_counter() - t_dl:.3f}s", flush=True)
        if csv_content:
            safe_mag = "".join(c if c.isalnum() or c in "-_" else "_" for c in mag)
            if download_widgets['mag_metrics']:
                if sample_names:
                    safe_sample = "".join(c if c.isalnum() or c in "-_" else "_" for c in sample_names[0])
                    download_widgets['mag_metrics'].filename = f"{safe_mag}_in_{safe_sample}_mag_metrics.csv"
                else:
                    download_widgets['mag_metrics'].filename = f"{safe_mag}_mag_metrics.csv"
            return io.StringIO(csv_content)
        return io.StringIO("")

    ### Creating all DOM elements
    # Open DuckDB database connection to build widgets depending on data
    if enable_timing:
        t_init = time.perf_counter()
    conn = duckdb.connect(db_path, read_only=True)
    widgets = build_controls(conn)
    if enable_timing:
        print(f"[timing] build_controls (initial DB queries): {time.perf_counter() - t_init:.3f}s", flush=True)

    # Build subplot → variable_name(s) mapping for inspect command generation
    _subplot_to_varnames = {}
    _cur = conn.cursor()
    _cur.execute("SELECT Variable_name, Subplot FROM Variable WHERE Variable_name IS NOT NULL AND Subplot IS NOT NULL")
    for _vname, _subplot in _cur.fetchall():
        _subplot_to_varnames.setdefault(_subplot, []).append(_vname)
    _cur.close()

    if enable_timing:
        t_ui = time.perf_counter()
        t_section = time.perf_counter()

    # Load the CSS and logo
    static_path = os.path.join(os.path.dirname(__file__), "..", "static")
    css_path = os.path.join(static_path, "bokeh_styles.css")
    with open(css_path) as f:
        css_text = f.read()
    stylesheet = InlineStyleSheet(css=css_text)

    # Pink buttons css
    pink_buttons_css_path = os.path.join(static_path, "pink_buttons.css")
    with open(pink_buttons_css_path) as f:
        pink_buttons_css_text = f.read()
    pink_buttons_stylesheet = InlineStyleSheet(css=pink_buttons_css_text)
    widgets['view_radio'].stylesheets = [pink_buttons_stylesheet]

    # Separate stylesheet for toggle buttons (minimal styling)
    toggle_css_path = os.path.join(static_path, "toggle_styles.css")
    with open(toggle_css_path) as f:
        toggle_css_text = f.read()
    toggle_stylesheet = InlineStyleSheet(css=toggle_css_text)

    # Create main elements
    ## Views section
    logo_url = "https://raw.githubusercontent.com/bhagavadgitadu22/theBIGbam/master/thebigbam/static/LOGO.png"
    logo = Div(text=f"""<img src="{logo_url}" style="width:100%; max-width:800px; padding: 0 25%;">""")
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
    _filter_help_tooltip = Tooltip(
        content=(
            'Filtering rows are independent from each other. '
            'For example, "Number of samples" is based on the total number of presences '
            'of a genome in the database, and is not recomputed based on the other filtering rows'
        ),
        position="right",
    )
    _filter_help_btn = HelpButton(
        tooltip=_filter_help_tooltip,
        width=20, height=20, align="center",
        button_type="light",
        stylesheets=[toggle_stylesheet],
    )
    filtering_header = row(filtering_toggle_btn, filtering_title, _filter_help_btn, sizing_mode="stretch_width", align="center")

    # Cache filtering metadata once when document loads
    if enable_timing:
        t_fmeta = time.perf_counter()
    filtering_metadata = get_filtering_metadata(db_path)
    if enable_timing:
        print(f"[timing]   get_filtering_metadata: {time.perf_counter() - t_fmeta:.3f}s", flush=True)

    # Get Sample table columns for ordering dropdown (exclude ID and name columns)
    sample_order_columns = ["Sample name"]  # Default option
    if 'Sample' in filtering_metadata:
        sample_columns = list(filtering_metadata['Sample']['columns'].keys())
        sample_order_columns.extend(sample_columns)

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
        from bokeh.io import curdoc
        doc = curdoc()
        doc.hold('combine')
        _filtering_cache['valid'] = False
        global_toggle_lock['locked'] = True
        try:
            refresh_contig_options_unlocked()
            refresh_sample_options_unlocked()
            refresh_mag_options_unlocked()
        finally:
            global_toggle_lock['locked'] = False
        update_section_titles()
        doc.unhold()

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
            options=["=", "!=", "has", "has not"] if initial_is_text else ["=", ">", "<", "!="],
            value="=" if initial_is_text else ">",
            width=50,
            margin=(0, 2, 0, 0)
        )

        # Container for the dynamic input widget
        input_container = pn.Column(width=90, margin=(0, 2, 0, 0))

        # Create initial input widget based on column type
        if initial_is_text:
            distinct_values = resolve_distinct_values(db_path, filtering_metadata, initial_category, initial_column)
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
                comparison_select.options = ["=", "!=", "has", "has not"]
                if comparison_select.value not in comparison_select.options:
                    comparison_select.value = "="
            else:
                comparison_select.options = ["=", ">", "<", "!="]
                if comparison_select.value not in comparison_select.options:
                    comparison_select.value = "="

            # Create new input widget
            if is_text:
                current_op = comparison_select.value
                if current_op in ("has", "has not"):
                    new_input = TextInput(value="", placeholder="Search...", width=90, margin=(0, 2, 0, 0))
                    input_container.objects = [new_input]
                    current_input_ref['widget'] = new_input
                    current_input_ref['is_panel'] = False
                    new_input.on_change('value', lambda attr, old, new: refresh_on_filter_change())
                else:
                    distinct_values = resolve_distinct_values(db_path, filtering_metadata, category, col_name)
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

        def update_input_on_operator_change(attr, old, new):
            """Swap input widget between TextInput and SearchableSelect based on operator."""
            category = category_select.value
            col_name = subcategory_select.value
            col_info = filtering_metadata.get(category, {}).get('columns', {}).get(col_name, {})
            is_text = col_info.get('type') == 'text'

            if is_text:
                if new in ("has", "has not"):
                    new_input = TextInput(value="", placeholder="Search...", width=90, margin=(0, 2, 0, 0))
                    input_container.objects = [new_input]
                    current_input_ref['widget'] = new_input
                    current_input_ref['is_panel'] = False
                    new_input.on_change('value', lambda attr, old, new: refresh_on_filter_change())
                else:
                    distinct_values = resolve_distinct_values(db_path, filtering_metadata, category, col_name)
                    new_input = SearchableSelect(
                        value="", options=distinct_values,
                        placeholder="Search...", width=90
                    )
                    input_container.objects = [new_input]
                    current_input_ref['widget'] = new_input
                    current_input_ref['is_panel'] = True
                    new_input.param.watch(lambda event: refresh_on_filter_change(), 'value')

            refresh_on_filter_change()

        category_select.on_change('value', update_subcategories)
        subcategory_select.on_change('value', update_input_on_column_change)
        comparison_select.on_change('value', update_input_on_operator_change)

        query_row = pn.Row(category_select, subcategory_select, comparison_select, input_container, minus_btn,
                       sizing_mode="stretch_width", margin=(2, 0, 2, 0))

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
            # Add AND/OR div before each row except the first
            if i > 0:
                # Reuse existing widget to preserve user's AND/OR selection
                if row_data['and_div'] is not None:
                    select_widget = row_data['and_div']
                else:
                    select_widget = Select(
                        options=["AND", "OR"],
                        value="AND",
                        margin=(2, 0, 2, 0)
                    )
                    select_widget.on_change('value', lambda attr, old, new: refresh_on_filter_change())
                    row_data['and_div'] = select_widget
                section_children.append(select_widget)
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
            margin=(2, 0, 2, 0),
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
        # Preserve existing inter-section AND/OR values before clearing
        old_values = [s.value for s in inter_section_selects]
        inter_section_selects.clear()

        for i, section_data in enumerate(or_sections):
            # Add AND/OR div before each section except the first
            if i > 0:
                old_val = old_values[i - 1] if i - 1 < len(old_values) else "AND"
                select_widget = Select(
                    options=["AND", "OR"],
                    value=old_val,
                    margin=(2, 0, 2, 0)
                )
                select_widget.on_change('value', lambda attr, old, new: refresh_on_filter_change())
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

    if enable_timing:
        print(f"[timing]   Filtering section: {time.perf_counter() - t_section:.3f}s", flush=True)
        t_section = time.perf_counter()


    ## Build Sample section
    sample_title = Div(text="<b>Samples</b>")

    above_sample_children = []
    above_sample_content = column(
        *above_sample_children,
        visible=True, sizing_mode="stretch_width"
    )

    def _on_sample_change(event):
        if global_toggle_lock['locked']:
            return
        from bokeh.io import curdoc
        doc = curdoc()
        doc.hold('combine')
        global_toggle_lock['locked'] = True
        try:
            if widgets['has_mags']:
                # MAG mode: sample → MAGs → contigs chain.
                # Filter MAGs first (may clear an invalid MAG selection),
                # then update contig list based on the (now-current) MAG value.
                refresh_mag_options_unlocked()
                refresh_contig_options_unlocked()
            else:
                refresh_contig_options_unlocked()
        finally:
            global_toggle_lock['locked'] = False
        update_section_titles()
        doc.unhold()
    widgets['sample_select'].param.watch(_on_sample_change, 'value')


    ## Build Contig section
    # Keep original full lists so we can restore when filters are off
    orig_contigs = list(widgets['contigs'])
    orig_samples = list(widgets['samples'])

    # MAGs section header (visible only in MAG-mode databases)
    mag_title = Div(
        text=f"<span style='font-size: 1.2em;'><b>MAGs</b></span> ({len(widgets['mags'])} available)",
        align="center",
        visible=widgets['has_mags'],
    )
    mag_header = row(mag_title, sizing_mode="stretch_width", align="center", margin=(0, 0, 0, 0), visible=widgets['has_mags'])
    separator_mags = Div(text="", height=2, sizing_mode="stretch_width", styles={'background-color': '#333', 'margin-top': '10px', 'margin-bottom': '10px'}, visible=widgets['has_mags'])

    # Create Contigs section header — no outer collapse; the inner
    # "Genomic annotations" and "Other genomic features" subsections each
    # have their own toggle below.
    contig_title = Div(text="<b>Contigs</b>", align="center")
    contig_header = row(contig_title, sizing_mode="stretch_width", align="center", margin=(0, 0, 0, 0))

    def on_contig_change(event):
        if global_toggle_lock['locked']:
            return
        from bokeh.io import curdoc
        doc = curdoc()
        doc.hold('combine')
        new = event.new
        global_toggle_lock['locked'] = True
        try:
            # In MAG mode the sample filter is driven by the MAG, not the contig directly.
            # on_contig_sync_mag will fire next and update the MAG + refresh samples.
            if not widgets['has_mags']:
                refresh_sample_options_unlocked()
        finally:
            global_toggle_lock['locked'] = False
        update_section_titles()
        # Update position inputs when contig changes
        if widgets['has_mags'] and widgets['view_radio'].active == 0:
            # MAG view: zoom From/To to the selected contig's range in MAG space
            selected_mag = widgets['mag_select'].value
            offsets = widgets['mag_to_contig_offsets'].get(selected_mag, {})
            if new and selected_mag and new in offsets:
                off = offsets[new]
                c_len = widgets['contig_lengths'].get(new, 0)
                from_position_input.value = str(off + 1)
                to_position_input.value = str(off + c_len)
        elif new and new in widgets['contig_lengths']:
            from_position_input.value = "1"
            to_position_input.value = str(widgets['contig_lengths'][new])
        else:
            from_position_input.value = "1"
            to_position_input.value = ""
        doc.unhold()

    widgets['contig_select'].param.watch(on_contig_change, 'value')

    # MAG cross-filter callbacks (only attached when DB is MAG-mode)
    if widgets['has_mags']:
        def on_mag_change(event):
            if global_toggle_lock['locked']:
                return
            from bokeh.io import curdoc
            doc = curdoc()
            doc.hold('combine')
            new = event.new
            global_toggle_lock['locked'] = True
            try:
                # Contigs: children of selected MAG
                refresh_contig_options_unlocked()
                # Samples: filter to those containing this MAG
                refresh_sample_options_unlocked()
            finally:
                global_toggle_lock['locked'] = False
            update_section_titles()
            # In MAG view: set From/To to the full new MAG length
            if widgets['view_radio'].active == 0 and new:
                total = sum(
                    widgets['contig_lengths'].get(c, 0)
                    for c in widgets['mag_to_contigs'].get(new, [])
                )
                from_position_input.value = "1"
                to_position_input.value = str(total)
            doc.unhold()

        widgets['mag_select'].param.watch(on_mag_change, 'value')

        # When contig changes, auto-select its parent MAG (user may clear to override),
        # then refresh sample list for that MAG.
        def on_contig_sync_mag(event):
            if global_toggle_lock['locked']:
                return
            from bokeh.io import curdoc
            doc = curdoc()
            doc.hold('combine')
            contig_to_mag = widgets['contig_to_mag']
            new = event.new
            mag_select = widgets['mag_select']
            if new and new in contig_to_mag:
                parent = contig_to_mag[new]
                changed_mag = (mag_select.value != parent)
                global_toggle_lock['locked'] = True
                try:
                    if changed_mag:
                        mag_select.value = parent
                    # Refresh samples for the (possibly new) parent MAG
                    refresh_sample_options_unlocked()
                finally:
                    global_toggle_lock['locked'] = False
                if changed_mag:
                    update_section_titles()
            elif not new:
                # Contig cleared: restore full sample list in MAG mode
                global_toggle_lock['locked'] = True
                try:
                    refresh_sample_options_unlocked()
                finally:
                    global_toggle_lock['locked'] = False
                update_section_titles()
            doc.unhold()

        widgets['contig_select'].param.watch(on_contig_sync_mag, 'value')

        def on_view_change(attr, old, new):
            is_mag_view = (new == 0)
            if is_mag_view:
                selected_mag = widgets['mag_select'].value
                selected_contig = widgets['contig_select'].value
                offsets = widgets['mag_to_contig_offsets'].get(selected_mag, {})
                if selected_mag and selected_contig and selected_contig in offsets:
                    # Contig already selected: zoom to its range in MAG space
                    off = offsets[selected_contig]
                    c_len = widgets['contig_lengths'].get(selected_contig, 0)
                    from_position_input.value = str(off + 1)
                    to_position_input.value = str(off + c_len)
                elif selected_mag:
                    total = sum(
                        widgets['contig_lengths'].get(c, 0)
                        for c in widgets['mag_to_contigs'].get(selected_mag, [])
                    )
                    from_position_input.value = "1"
                    to_position_input.value = str(total)
            else:
                selected_contig = widgets['contig_select'].value
                if selected_contig and selected_contig in widgets['contig_lengths']:
                    from_position_input.value = "1"
                    to_position_input.value = str(widgets['contig_lengths'][selected_contig])

        widgets['view_radio'].on_change('active', on_view_change)


    if enable_timing:
        print(f"[timing]   Sample/Contig/MAG sections: {time.perf_counter() - t_section:.3f}s", flush=True)
        t_section = time.perf_counter()

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

    # Load color templates from database
    color_templates = {}  # {template_name: [{qualifier_name, operator, qualifier_value, color}, ...]}
    try:
        template_rows = conn.execute(
            "SELECT t.Template_name, r.Qualifier_name, r.Operator, r.Qualifier_value, r.Color "
            "FROM Color_templates t JOIN Color_rules r ON t.Template_id = r.Template_id "
            "ORDER BY t.Template_name, r.Rule_id"
        ).fetchall()
        for tname, qname, op, qvalue, color in template_rows:
            color_templates.setdefault(tname, []).append({
                'qualifier_name': qname, 'operator': op,
                'qualifier_value': qvalue, 'color': color,
            })
    except Exception:
        pass

    template_select = None
    if color_templates:
        template_options = ["(none)"] + list(color_templates.keys())
        template_select = Select(
            title="Use a color template:",
            value="(none)",
            options=template_options,
            sizing_mode="stretch_width"
        )

    # Build qualifier key options for custom color rows — reuse the exact same
    # list the Filtering section iterates over (text + numeric alike).
    annotation_meta = filtering_metadata.get('Annotations', {}).get('columns', {})
    color_qualifier_options = list(annotation_meta.keys())

    # Custom color row system (Panel widgets, same pattern as filtering rows)
    custom_color_rows = []
    add_color_btn = pn.widgets.Button(
        name="+ Add coloring rule",
        margin=(2, 0, 2, 0),
        button_type="success",
        stylesheets=[stylesheet]
    )
    custom_color_column = pn.Column(
        add_color_btn, sizing_mode="stretch_width",
        styles={'border-left': '3px solid #00b17c', 'padding-left': '10px', 'margin-left': '5px',
                'max-height': '300px', 'overflow-y': 'auto'}
    )

    TEXT_OPS = ["=", "!=", "has", "has not", "Use random colors"]
    NUMERIC_OPS = ["=", ">", "<", "!=", "Use random colors"]

    def _build_color_value_widget(is_text, operator, distinct_values):
        """Return the value widget that matches (type, operator).

        Widgets live directly inside the row (no wrapper container), so each
        one carries its own stretch_width sizing and the 2px right margin
        that separates it from the color picker.
        """
        if is_text:
            if operator in ("has", "has not"):
                return TextInput(value="", placeholder="Search...",
                                 sizing_mode="stretch_width",
                                 margin=(0, 2, 0, 0)), False
            # = / != / Use random colors → SearchableSelect (kept even in random
            # mode so the widget shape is consistent when user flips back).
            return SearchableSelect(
                value="", options=[str(v) for v in distinct_values],
                placeholder="Search...", sizing_mode="stretch_width",
                margin=(0, 2, 0, 0),
            ), True
        # Numeric columns always get a Spinner regardless of operator.
        return Spinner(value=0, placeholder="Value...",
                       sizing_mode="stretch_width",
                       margin=(0, 2, 0, 0)), False

    def create_color_row():
        """Create a single custom color row with qualifier / operator / value / color / remove widgets.

        Mirrors create_query_row(): text columns get SearchableSelect (or
        TextInput under has/has not), numeric columns get Spinner. Both types
        gain a 'Use random colors' operator that hides the value container and
        the color picker.
        """
        initial_key = color_qualifier_options[0] if color_qualifier_options else ""
        initial_info = annotation_meta.get(initial_key, {})
        initial_is_text = initial_info.get('type') == 'text'
        initial_distinct = resolve_distinct_values(db_path, filtering_metadata, 'Annotations', initial_key) if initial_is_text else []

        qualifier_select = Select(
            options=[(k, k.replace("_", " ").replace("percentage", "(%)")) for k in color_qualifier_options],
            value=initial_key,
            width=100,
            margin=(0, 2, 0, 0)
        )

        operator_select = Select(
            options=TEXT_OPS if initial_is_text else NUMERIC_OPS,
            value="=",
            width=50,
            margin=(0, 2, 0, 0),
        )

        # Dynamic value widget sits directly at index 2 of the Row — swapped
        # in place when the qualifier type or operator changes. No wrapper
        # pn.Column, because stretch_width on a wrapper introduces extra
        # vertical/horizontal padding that pushed the widget out of line.
        initial_input, initial_is_panel = _build_color_value_widget(
            initial_is_text, "=", initial_distinct
        )
        current_input_ref = {'widget': initial_input, 'is_panel': initial_is_panel}

        color_picker = ColorPicker(color="#cccccc", width=60, height=30, margin=(0, 2, 0, 0))

        minus_btn = pn.widgets.Button(name="\u2212", width=30, height=30, margin=(0, 10, 0, 0), stylesheets=[stylesheet])

        row_widget = pn.Row(qualifier_select, operator_select, initial_input,
                            color_picker, minus_btn,
                            sizing_mode="stretch_width", margin=(2, 0, 2, 0))

        VALUE_IDX = 2  # position of the value widget inside row_widget

        row_data = {
            'qualifier_select': qualifier_select,
            'operator_select': operator_select,
            'input_ref': current_input_ref,
            'color_picker': color_picker,
            'minus_btn': minus_btn,
            'row_widget': row_widget,
        }

        def _swap_value_widget(new_widget, is_panel):
            row_widget[VALUE_IDX] = new_widget
            current_input_ref['widget'] = new_widget
            current_input_ref['is_panel'] = is_panel

        def _apply_random_visibility():
            use_random = (operator_select.value == "Use random colors")
            current_input_ref['widget'].visible = not use_random
            color_picker.visible = not use_random
            # In random mode, the value widget and color picker are hidden;
            # let the operator Select grow to fill the freed space. Revert to
            # a fixed 50 px when any other operator is picked.
            if use_random:
                operator_select.sizing_mode = "stretch_width"
            else:
                operator_select.sizing_mode = "fixed"
                operator_select.width = 50

        def update_color_input_widget(col_name):
            """Rebuild operator options and the value widget for the new column type."""
            col_info = annotation_meta.get(col_name, {})
            is_text = col_info.get('type') == 'text'
            distinct = resolve_distinct_values(db_path, filtering_metadata, 'Annotations', col_name) if is_text else []

            # Swap operator options for the new type; preserve current operator
            # if still valid, otherwise default back to "=".
            new_ops = TEXT_OPS if is_text else NUMERIC_OPS
            operator_select.options = new_ops
            if operator_select.value not in new_ops:
                operator_select.value = "="

            new_widget, is_panel = _build_color_value_widget(
                is_text, operator_select.value, distinct
            )
            _swap_value_widget(new_widget, is_panel)

            # Respect random-mode visibility even after a qualifier swap.
            _apply_random_visibility()

        def on_qualifier_change(attr, old, new):
            update_color_input_widget(new)

        qualifier_select.on_change('value', on_qualifier_change)

        def on_operator_change(attr, old, new):
            col_info = annotation_meta.get(qualifier_select.value, {})
            is_text = col_info.get('type') == 'text'
            use_random = (new == "Use random colors")

            # For text columns, has/has not uses TextInput; everything else
            # uses SearchableSelect. Only swap when the current widget doesn't
            # already match, so user-typed values aren't clobbered needlessly.
            if is_text and not use_random:
                want_text_input = new in ("has", "has not")
                cur = current_input_ref['widget']
                if want_text_input and not isinstance(cur, TextInput):
                    new_input, is_panel = _build_color_value_widget(True, new, [])
                    _swap_value_widget(new_input, is_panel)
                elif not want_text_input and not isinstance(cur, SearchableSelect):
                    col_name_val = qualifier_select.value
                    distinct = resolve_distinct_values(db_path, filtering_metadata, 'Annotations', col_name_val)
                    new_input, is_panel = _build_color_value_widget(True, new, distinct)
                    _swap_value_widget(new_input, is_panel)

            _apply_random_visibility()

        operator_select.on_change('value', on_operator_change)

        def remove_row_callback(event):
            if row_data in custom_color_rows:
                custom_color_rows.remove(row_data)
                rebuild_color_rows()

        minus_btn.on_click(remove_row_callback)
        return row_data

    def rebuild_color_rows():
        children = [rd['row_widget'] for rd in custom_color_rows]
        children.append(add_color_btn)
        custom_color_column.objects = children

    def add_color_callback(event):
        new_row = create_color_row()
        custom_color_rows.append(new_row)
        rebuild_color_rows()

    add_color_btn.on_click(add_color_callback)

    # Template selection callback - populates color rows from template rules
    if template_select is not None:
        def on_template_change(attr, old, new):
            custom_color_rows.clear()
            if new != "(none)" and new in color_templates:
                for rule in color_templates[new]:
                    row_data = create_color_row()
                    # Order matters: qualifier first (triggers widget rebuild),
                    # then operator (may swap SearchableSelect ↔ TextInput),
                    # then the value on whichever widget is now current.
                    row_data['qualifier_select'].value = rule['qualifier_name']
                    row_data['operator_select'].value = rule['operator']
                    if rule['operator'] != "Use random colors":
                        widget = row_data['input_ref']['widget']
                        if isinstance(widget, SearchableSelect):
                            # Seed values (e.g. "DNA, RNA and nucleotide metabolism")
                            # may not be in the current distinct list; append so
                            # the Select can display them.
                            opts = list(widget.options)
                            if rule['qualifier_value'] not in opts:
                                widget.options = opts + [rule['qualifier_value']]
                        widget.value = rule['qualifier_value']
                        row_data['color_picker'].color = rule['color']
                    custom_color_rows.append(row_data)
            rebuild_color_rows()

        template_select.on_change('value', on_template_change)

    # Feature label dropdown — same qualifier list the Filtering Annotations
    # category and the coloring rules use, so every annotation attribute (KV
    # keys + direct Contig_annotation columns) can be picked as the tooltip.
    feature_label_select = None
    label_keys = color_qualifier_options
    if label_keys:
        feature_label_select = Select(
            title="Label annotations with:",
            value="product" if "product" in label_keys else label_keys[0],
            options=label_keys,
            sizing_mode="stretch_width",
            margin=(5, 5, 5, 5),
        )

    # Plot isoforms checkbox - only show if at least one locus_tag appears more than once
    plot_isoforms_cbg = None
    cur = conn.cursor()
    result = cur.execute(
        "SELECT Status FROM Constants_for_plotting WHERE Constant = 'isoforms'"
    ).fetchone()
    has_isoforms = result[0] if result else False

    if has_isoforms:
        plot_isoforms_cbg = CheckboxGroup(
            labels=["Plot isoforms"],
            active=[]  # Unchecked by default
        )

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
    from_position_input = TextInput(value="1", placeholder="Start position", sizing_mode="stretch_width", margin=(0, 0, 0, 0))
    to_position_input = TextInput(value="", placeholder="End position", sizing_mode="stretch_width", margin=(0, 0, 0, 0))
    
    position_label_from = Div(text="From", width=40, margin=(5, 0, 5, 5))
    position_label_to = Div(text="to", width=25, margin=(5, 0, 5, 5))
    
    # Create Reset button to reset position inputs
    position_reset_button = Button(label="Reset", stylesheets=[stylesheet], margin=(0, 5, 0, 5))
    
    def reset_position_inputs():
        from_position_input.value = "1"
        if widgets['has_mags'] and widgets['view_radio'].active == 0:
            # MAG view: clear contig selection and show full MAG extent
            widgets['contig_select'].value = ""
            selected_mag = widgets['mag_select'].value
            if selected_mag:
                total = sum(
                    widgets['contig_lengths'].get(c, 0)
                    for c in widgets['mag_to_contigs'].get(selected_mag, [])
                )
                to_position_input.value = str(total)
            else:
                to_position_input.value = ""
        elif widgets['contig_select'].value and widgets['contig_select'].value in widgets['contig_lengths']:
            to_position_input.value = str(widgets['contig_lengths'][widgets['contig_select'].value])
        else:
            to_position_input.value = ""
    
    position_reset_button.on_click(lambda event: reset_position_inputs())
    
    position_row = row(
        position_label_from, from_position_input, 
        position_label_to, to_position_input,
        position_reset_button,
        sizing_mode="stretch_width",
        margin=(10, 0, 5, 0)
    )
    
    below_contig_children.append(position_row)

    # Check if sequence data is available in the database
    cur = conn.cursor()
    cur.execute("SELECT 1 FROM information_schema.tables WHERE table_name = 'Contig_sequence'")
    has_sequence_data = cur.fetchone() is not None

    sequence_cbg = None
    sequence_row = None
    if has_sequence_data:
        sequence_cbg = CheckboxGroup(labels=["Plot sequence"], active=[])
        sequence_row = row(sequence_cbg, sizing_mode="stretch_width")

    # Check if translated annotation data is available
    cur = conn.cursor()
    cur.execute("SELECT 1 FROM information_schema.columns WHERE table_name = 'Contig_annotation' AND column_name = 'Protein_sequence'")
    has_translated_data = cur.fetchone() is not None

    translated_sequence_cbg = None
    translated_sequence_row = None
    if has_translated_data:
        translated_sequence_cbg = CheckboxGroup(labels=["Plot translated sequence"], active=[])
        translated_sequence_row = row(translated_sequence_cbg, sizing_mode="stretch_width")

    if feature_type_multichoice is not None or combined_features_cbg is not None:
        # --- Subsection 1: "Genomic annotations to plot" (collapsible) ---
        # Covers the feature type multichoice, isoform toggle, and both
        # sequence/translated-sequence rows. Coloring and labelling widgets
        # live in the separate "Customise genomic annotations" subsection.
        annotations_toggle_btn = Button(label="▶", width=20, height=20,
                                        button_type="primary", align="center",
                                        margin=0, stylesheets=[toggle_stylesheet])
        annotations_toggle_btn.styles = {'padding': '0px', 'line-height': '20px'}
        annotations_title = Div(text="Genomic annotations to plot", align="center")
        annotations_header = row(annotations_toggle_btn, annotations_title,
                                 sizing_mode="stretch_width", align="center")

        annotations_children = []
        if feature_type_multichoice is not None:
            annotations_children.append(feature_type_multichoice)
        if plot_isoforms_cbg is not None:
            annotations_children.append(plot_isoforms_cbg)
        if sequence_row is not None:
            annotations_children.append(sequence_row)
        if translated_sequence_row is not None:
            annotations_children.append(translated_sequence_row)
        annotations_content = pn.Column(*annotations_children, visible=False,
                                        sizing_mode="stretch_width",
                                        margin=(0, 0, 0, 0))
        annotations_toggle_btn.on_click(
            make_toggle_callback(annotations_toggle_btn, annotations_content))

        # --- Subsection 2: "Customise genomic annotations" (collapsible) ---
        # Independent sibling of the annotations-to-plot section. Bundles
        # the coloring rules (template + custom color rows) and the label
        # dropdown. Starts collapsed.
        customise_header = None
        customise_content = None
        has_customise = (
            template_select is not None
            or color_qualifier_options
            or feature_label_select is not None
        )
        if has_customise:
            customise_toggle_btn = Button(label="▶", width=20, height=20,
                                          button_type="primary", align="center",
                                          margin=0, stylesheets=[toggle_stylesheet])
            customise_toggle_btn.styles = {'padding': '0px', 'line-height': '20px'}
            customise_title = Div(text="Customise genomic annotations", align="center")
            customise_header = row(customise_toggle_btn, customise_title,
                                   sizing_mode="stretch_width", align="center")

            customise_children = []
            if feature_label_select is not None:
                customise_children.append(feature_label_select)
            if template_select is not None:
                customise_children.append(template_select)
            if template_select is not None or color_qualifier_options:
                customise_children.append(Div(text="Color annotations with:"))
            if color_qualifier_options:
                customise_children.append(custom_color_column)
            customise_content = pn.Column(*customise_children, visible=False,
                                          sizing_mode="stretch_width",
                                          margin=(0, 0, 5, 0))
            customise_toggle_btn.on_click(
                make_toggle_callback(customise_toggle_btn, customise_content))

        # --- Subsection 3: "Other genomic features to plot" (collapsible) ---
        other_features_toggle_btn = Button(label="▶", width=20, height=20,
                                           button_type="primary", align="center",
                                           margin=0, stylesheets=[toggle_stylesheet])
        other_features_toggle_btn.styles = {'padding': '0px', 'line-height': '20px'}
        # Master checkbox whose label doubles as the section title — same
        # pattern Variables modules use (CheckboxGroup with [module_name]).
        # Toggling it checks/unchecks every feature in combined_features_cbg
        # at once; bidirectional sync is wired below.
        genome_master_cbg = CheckboxGroup(labels=["Other genomic features to plot"], active=[])

        genome_help_tooltip = widgets['helps_widgets'][genome_index_one] if genome_index_one is not None else None
        if genome_help_tooltip is not None:
            help_btn = HelpButton(tooltip=genome_help_tooltip, width=20, height=20,
                                  align="center", button_type="light",
                                  stylesheets=[toggle_stylesheet])
            other_features_header = row(other_features_toggle_btn, genome_master_cbg, help_btn,
                                        sizing_mode="stretch_width", align="center")
        else:
            other_features_header = row(other_features_toggle_btn, genome_master_cbg,
                                        sizing_mode="stretch_width", align="center")

        other_features_content = None
        if combined_features_cbg is not None:
            combined_features_cbg.visible = True
            other_features_content = pn.Column(combined_features_cbg, visible=False,
                                               sizing_mode="stretch_width",
                                               margin=(0, 0, 0, 0))
            other_features_toggle_btn.on_click(
                make_toggle_callback(other_features_toggle_btn, other_features_content))

            # Bidirectional sync: master ⇄ individual feature toggles.
            # Mirrors make_module_callback / make_variable_callback used by
            # the Variables subsections at ~line 2025.
            genome_master_lock = {"locked": False}

            def _on_genome_master_change(attr, old, new):
                if genome_master_lock["locked"] or global_toggle_lock.get("locked", False):
                    return
                genome_master_lock["locked"] = True
                try:
                    master_on = 0 in genome_master_cbg.active
                    if master_on:
                        combined_features_cbg.active = list(range(len(combined_features_cbg.labels)))
                    else:
                        combined_features_cbg.active = []
                finally:
                    genome_master_lock["locked"] = False

            def _on_combined_features_change(attr, old, new):
                if genome_master_lock["locked"] or global_toggle_lock.get("locked", False):
                    return
                total = len(combined_features_cbg.labels)
                active_count = len(combined_features_cbg.active)
                genome_master_lock["locked"] = True
                try:
                    if total > 0 and active_count == total:
                        genome_master_cbg.active = [0]
                    else:
                        genome_master_cbg.active = []
                finally:
                    genome_master_lock["locked"] = False

            genome_master_cbg.on_change("active", _on_genome_master_change)
            combined_features_cbg.on_change("active", _on_combined_features_change)

        # Assemble the full genome section from the independent collapsibles.
        genome_section_children = [annotations_header, annotations_content]
        if customise_header is not None:
            genome_section_children.append(customise_header)
            if customise_content is not None:
                genome_section_children.append(customise_content)
        genome_section_children.append(other_features_header)
        if other_features_content is not None:
            genome_section_children.append(other_features_content)
        genome_section = pn.Column(*genome_section_children, visible=True,
                                   sizing_mode="stretch_width", margin=(0, 0, 0, 0))

    if genome_section is not None:
        below_contig_children = list(below_contig_children) + [genome_section]

    # Fallback: if no genome section exists, append sequence/translated rows directly
    if genome_section is None and sequence_row is not None:
        below_contig_children.append(sequence_row)
    if genome_section is None and translated_sequence_row is not None:
        below_contig_children.append(translated_sequence_row)

    below_contig_content = pn.Column(
        *below_contig_children,
        visible=True, sizing_mode="stretch_width", margin=(0, 0, 0, 0)
    )

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

    if enable_timing:
        print(f"[timing]   Variables + Genome/Annotations sections: {time.perf_counter() - t_section:.3f}s", flush=True)
        t_section = time.perf_counter()

    ## Plotting parameters section
    separator_plotting_params = Div(text="", height=2, sizing_mode="stretch_width",
        styles={'background-color': '#333', 'margin-top': '10px', 'margin-bottom': '10px'})
    plotting_params_title = Div(text="<span style='font-size: 1.2em;'><b>Plotting parameters</b></span>", align="center")
    plotting_params_header = row(plotting_params_title, sizing_mode="stretch_width", align="center")

    # Sample paramaters (only useful in All Samples view)
    sample_order_label = Div(text="Order samples by:", margin=(5, 5, 5, 0))
    sample_order_select = Select(value="Sample name", options=sample_order_columns, sizing_mode="stretch_width", margin=(0, 5, 0, 5))
    sample_order_row = row(sample_order_label, sample_order_select, sizing_mode="stretch_width", margin=(5, 0, 5, 0))

    same_y_scale_cbg = CheckboxGroup(labels=["Use same y scale for all samples"], active=[], margin=(5, 0, 5, 0))
    same_y_scale_row = row(same_y_scale_cbg, sizing_mode="stretch_width")
    
    ## Plotting parameters useful in both views
    min_coverage_freq_input = Spinner(value=0.0, low=0.0, high=1.0, step=0.01, width=100, margin=(0, 2, 0, 0))
    min_coverage_freq_label = Div(text="Minimum frequency for coverage-related features", margin=(5, 0, 5, 5))
    min_coverage_freq_row = row(min_coverage_freq_input, min_coverage_freq_label, sizing_mode="stretch_width", margin=(0, 0, 5, 0))

    # Subsection: Max window sizes for plotting
    max_genemap_window_input = Spinner(value=DEFAULT_GENEMAP_WINDOW, low=10, high=1000000, step=1000, width=100, margin=(0, 2, 0, 0))
    max_genemap_window_label = Div(text="Gene map (bp)", margin=(5, 0, 5, 5))
    max_genemap_window_row = row(max_genemap_window_input, max_genemap_window_label, sizing_mode="stretch_width", margin=(5, 0, 5, 0))

    max_sequence_window_input = Spinner(value=DEFAULT_SEQUENCE_WINDOW, low=10, high=1000000, step=100, width=100, margin=(0, 2, 0, 0))
    max_sequence_window_label = Div(text="Sequence plots (bp)", margin=(5, 0, 5, 5))
    max_sequence_window_row = row(max_sequence_window_input, max_sequence_window_label, sizing_mode="stretch_width", margin=(0, 0, 5, 0))

    max_binning_window_input = Spinner(value=_DEFAULT_MAX_BASE_RESOLUTION, low=10, high=1000000, step=1000, width=100, margin=(0, 2, 0, 0))
    max_binning_window_label = Div(text="Feature plots without binning (bp)", margin=(5, 0, 5, 5))
    max_binning_window_row = row(max_binning_window_input, max_binning_window_label, sizing_mode="stretch_width", margin=(0, 0, 5, 0))

    max_window_toggle_btn = Button(label="▶", width=20, height=20, button_type="primary", align="center", margin=0, stylesheets=[toggle_stylesheet])
    max_window_toggle_btn.styles = {'padding': '0px', 'line-height': '20px'}
    max_window_title = Div(text="Max window size for plotting", align="center")
    max_window_header = row(max_window_toggle_btn, max_window_title, sizing_mode="stretch_width", align="center", margin=(5, 0, 0, 0))
    max_window_content = pn.Column(
        max_genemap_window_row, max_sequence_window_row, max_binning_window_row,
        sizing_mode="stretch_width", visible=False
    )
    max_window_toggle_btn.on_click(make_toggle_callback(max_window_toggle_btn, max_window_content))

    # Subsection: Plot heights
    genemap_height_input = Spinner(value=100, low=10, high=1000, step=10, width=80, margin=(0, 2, 0, 0))
    genemap_height_label = Div(text="Of gene map (px)", margin=(5, 0, 5, 5))
    genemap_height_row = row(genemap_height_input, genemap_height_label, sizing_mode="stretch_width", margin=(0, 0, 5, 0))

    sequence_height_input = Spinner(value=50, low=10, high=1000, step=10, width=80, margin=(0, 2, 0, 0))
    sequence_height_label = Div(text="Of nucleotide sequence (px)", margin=(5, 0, 5, 5))
    sequence_height_row = row(sequence_height_input, sequence_height_label, sizing_mode="stretch_width", margin=(0, 0, 5, 0))

    translated_sequence_height_input = Spinner(value=50, low=10, high=1000, step=10, width=80, margin=(0, 2, 0, 0))
    translated_sequence_height_label = Div(text="Of translated sequence (px)", margin=(5, 0, 5, 5))
    translated_sequence_height_row = row(translated_sequence_height_input, translated_sequence_height_label, sizing_mode="stretch_width", margin=(0, 0, 5, 0))

    subplot_height_input = Spinner(value=100, low=10, high=1000, step=10, width=80, margin=(0, 2, 0, 0))
    subplot_height_label = Div(text="Per feature plot (px)", margin=(5, 0, 5, 5))
    subplot_height_row = row(subplot_height_input, subplot_height_label, sizing_mode="stretch_width")

    plot_heights_toggle_btn = Button(label="▶", width=20, height=20, button_type="primary", align="center", margin=0, stylesheets=[toggle_stylesheet])
    plot_heights_toggle_btn.styles = {'padding': '0px', 'line-height': '20px'}
    plot_heights_title = Div(text="Plot heights", align="center")
    plot_heights_header = row(plot_heights_toggle_btn, plot_heights_title, sizing_mode="stretch_width", align="center", margin=(5, 0, 0, 0))
    plot_heights_content = pn.Column(
        genemap_height_row, sequence_height_row, translated_sequence_height_row, subplot_height_row,
        sizing_mode="stretch_width", visible=False
    )
    plot_heights_toggle_btn.on_click(make_toggle_callback(plot_heights_toggle_btn, plot_heights_content))

    # Subsection: Sample parameters (All Samples view only)
    sample_params_toggle_btn = Button(label="▶", width=20, height=20, button_type="primary", align="center", margin=0, stylesheets=[toggle_stylesheet])
    sample_params_toggle_btn.styles = {'padding': '0px', 'line-height': '20px'}
    sample_params_title = Div(text="Sample parameters", align="center")
    sample_params_header = row(sample_params_toggle_btn, sample_params_title, sizing_mode="stretch_width", align="center", margin=(5, 0, 0, 0))
    sample_params_header.visible = False  # Only shown in All Samples mode
    sample_params_content = pn.Column(
        sample_order_row, same_y_scale_row,
        sizing_mode="stretch_width", visible=False
    )
    sample_params_toggle_btn.on_click(make_toggle_callback(sample_params_toggle_btn, sample_params_content))

    plotting_params_content = pn.Column(
        min_coverage_freq_row,
        max_window_header, max_window_content,
        plot_heights_header, plot_heights_content,
        sample_params_header, sample_params_content,
        sizing_mode="stretch_width"
    )

    # Hide samples-only parameters when no BAM files were provided
    if not widgets['has_samples']:
        min_coverage_freq_row.visible = False
        # sample_params_header is already visible=False by default (All Samples mode only)

    ## Create final Apply and Peruse data buttons
    apply_button = Button(label="APPLY", align="center", stylesheets=[stylesheet], css_classes=["apply-btn"], margin=(5, 0, 0, 0))
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

    # Download contig metrics - Panel FileDownload widget (contig view only)
    download_metrics_button = pn.widgets.FileDownload(
        callback=make_contig_metrics_download_callback,
        filename="contig_metrics.csv",
        label="DOWNLOAD CONTIG METRICS",
        button_type="primary",
        height=30,
        visible=False
    )
    download_widgets['contig_metrics'] = download_metrics_button

    # Download MAG metrics - Panel FileDownload widget (both MAG and contig views)
    download_mag_metrics_button = pn.widgets.FileDownload(
        callback=make_mag_metrics_download_callback,
        filename="mag_metrics.csv",
        label="DOWNLOAD MAG METRICS",
        button_type="primary",
        height=30,
        visible=False,
    )
    download_widgets['mag_metrics'] = download_mag_metrics_button

    # Download data — show CLI command hint instead of generating huge CSV in browser
    # Hidden pane to display the inspect command (shown on click, hidden on next APPLY)
    command_hint_pane = pn.pane.HTML("", visible=False, sizing_mode="stretch_width")

    def _show_inspect_command(event):
        """Build and display the thebigbam inspect CLI command from current widget state."""
        # Toggle off if already visible
        if command_hint_pane.visible:
            command_hint_pane.visible = False
            return

        is_mag_view = widgets['has_mags'] and widgets['view_radio'].active == 0
        is_all = (views.active == 1) if widgets['has_samples'] else False
        has_samples = widgets['has_samples']

        # Collect active feature subplots → variable names
        active_vars = widgets['variables_widgets_all'] if is_all else widgets['variables_widgets_one']
        feature_names = []
        if combined_features_cbg is not None:
            for idx in combined_features_cbg.active:
                subplot = combined_features_cbg.labels[idx]
                feature_names.extend(_subplot_to_varnames.get(subplot, [subplot]))
        if is_all:
            for cbg in active_vars:
                if cbg.active:
                    subplot = cbg.labels[cbg.active[-1]]
                    feature_names.extend(_subplot_to_varnames.get(subplot, [subplot]))
                    break
        else:
            for cbg in active_vars:
                for idx in cbg.active:
                    subplot = cbg.labels[idx]
                    feature_names.extend(_subplot_to_varnames.get(subplot, [subplot]))

        if not feature_names:
            command_hint_pane.object = '<div style="color:#c00;padding:6px">No features selected.</div>'
            command_hint_pane.visible = True
            return

        from thebigbam.database.blob_decoder import is_contig_blob_feature

        # --- MAG view: generate --mag command ---
        if is_mag_view:
            mag = widgets['mag_select'].value
            if not mag:
                command_hint_pane.object = '<div style="color:#c00;padding:6px">No MAG selected.</div>'
                command_hint_pane.visible = True
                return

            # Region (position inputs apply to the MAG layout, not a contig position)
            mag_length = sum(
                widgets['contig_lengths'].get(c, 0)
                for c in widgets['mag_to_contigs'].get(mag, [])
            )
            try:
                xstart = int(from_position_input.value) if from_position_input.value.strip() else 1
                xend = int(to_position_input.value) if to_position_input.value.strip() else mag_length
            except ValueError:
                xstart = 1
                xend = mag_length
            if mag_length > 0:
                xstart = max(1, min(xstart, mag_length))
                xend = max(1, min(xend, mag_length))
            region_arg = ""
            if mag_length > 0 and xstart < xend and (xstart > 1 or xend < mag_length):
                region_arg = f" --region {xstart}-{xend}"

            sample_arg = ""
            if has_samples:
                sample = widgets['sample_select'].value
                if sample:
                    sample_arg = f" --sample {sample}"

            contig_features = [f for f in feature_names if is_contig_blob_feature(f)]
            sample_features = [f for f in feature_names if not is_contig_blob_feature(f)]

            commands = []
            if contig_features:
                cmd = f"thebigbam inspect -d {db_path} --mag {mag} --feature {','.join(contig_features)}{region_arg} > output.tsv"
                commands.append(cmd)
            if sample_features and sample_arg:
                cmd = f"thebigbam inspect -d {db_path} --mag {mag}{sample_arg} --feature {','.join(sample_features)}{region_arg} > output.tsv"
                commands.append(cmd)
            elif sample_features and not sample_arg:
                cmd = f"thebigbam inspect -d {db_path} --mag {mag} --sample &lt;SAMPLE&gt; --feature {','.join(sample_features)}{region_arg} > output.tsv"
                commands.append(cmd)

        # --- Contig view: generate --contig command ---
        else:
            contig = widgets['contig_select'].value
            if not contig:
                command_hint_pane.object = '<div style="color:#c00;padding:6px">No contig selected.</div>'
                command_hint_pane.visible = True
                return

            contig_length = widgets['contig_lengths'].get(contig, 0)
            try:
                xstart = int(from_position_input.value) if from_position_input.value.strip() else 1
                xend = int(to_position_input.value) if to_position_input.value.strip() else contig_length
            except ValueError:
                xstart = 1
                xend = contig_length
            if contig_length > 0:
                xstart = max(1, min(xstart, contig_length))
                xend = max(1, min(xend, contig_length))
            region_arg = ""
            if contig_length > 0 and xstart < xend and (xstart > 1 or xend < contig_length):
                region_arg = f" --region {xstart}-{xend}"

            sample_arg = ""
            if has_samples:
                if is_all:
                    filtered_samples = [s for s in orig_samples if s in widgets['contig_to_samples'].get(contig, set())]
                    filtering_pairs = get_filtering_filtered_pairs()
                    if filtering_pairs is not None:
                        allowed = {pair[1] for pair in filtering_pairs}
                        filtered_samples = [s for s in filtered_samples if s in allowed]
                    if filtered_samples:
                        sample_arg = f" --sample {','.join(filtered_samples)}"
                else:
                    sample = widgets['sample_select'].value
                    if sample:
                        sample_arg = f" --sample {sample}"

            contig_features = [f for f in feature_names if is_contig_blob_feature(f)]
            sample_features = [f for f in feature_names if not is_contig_blob_feature(f)]

            commands = []
            if contig_features:
                cmd = f"thebigbam inspect -d {db_path} --contig {contig} --feature {','.join(contig_features)}{region_arg} > output.tsv"
                commands.append(cmd)
            if sample_features and sample_arg:
                cmd = f"thebigbam inspect -d {db_path} --contig {contig}{sample_arg} --feature {','.join(sample_features)}{region_arg} > output.tsv"
                commands.append(cmd)
            elif sample_features and not sample_arg:
                cmd = f"thebigbam inspect -d {db_path} --contig {contig} --sample &lt;SAMPLE&gt; --feature {','.join(sample_features)}{region_arg} > output.tsv"
                commands.append(cmd)

        cmd_html = "".join(
            f'<pre style="background:#1a1a2e;color:#e0e0e0;padding:8px 12px;border-radius:4px;'
            f'font-size:13px;white-space:pre-wrap;word-break:break-all;margin:4px 0;'
            f'user-select:all;cursor:pointer" title="Click to select, then Ctrl+C to copy">{c}</pre>'
            for c in commands
        )
        command_hint_pane.object = (
            f'<div style="padding:6px 0">'
            f'<div style="font-size:12px;color:#888;margin-bottom:4px">'
            f'Run in your terminal to download data as TSV (click command to select):</div>'
            f'{cmd_html}</div>'
        )
        command_hint_pane.visible = True

    download_data_button = pn.widgets.Button(
        name="DOWNLOAD DATA",
        height=30,
        button_type="primary",
        visible=False
    )
    download_data_button.on_click(_show_inspect_command)
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
                             separator_mags, mag_header, widgets['view_radio'], widgets['mag_select'],
                             separator_contigs, contig_header, widgets['contig_select'], below_contig_content,
                             separator_samples, sample_title, above_sample_content, widgets['sample_select'],
                             separator_variables,
                             variables_section_one,  # One Sample view (with module checkboxes)
                             variables_section_all,  # All Samples view (title headers only)
                             separator_plotting_params, plotting_params_header, plotting_params_content,
                             buttons_row]
        placeholder_text = "<i>No plot yet. Select one sample, one contig and at least one variable in \"One sample\" mode or one contig and one variable in \"All samples\" mode and click Apply.</i>"
    else:
        # No samples: simplified UI with Filtering, Contigs, Plotting parameters, and Apply button
        controls_children = [logo, filtering_header, filtering_content,
                             separator_mags, mag_header, widgets['view_radio'], widgets['mag_select'],
                             separator_contigs, contig_header, widgets['contig_select'], below_contig_content,
                             separator_plotting_params, plotting_params_header, plotting_params_content,
                             buttons_row]
        placeholder_text = "<i>No plot yet. Select one contig and click Apply to view the genome annotation.</i>"

    if enable_timing:
        print(f"[timing]   Plotting params + layout assembly: {time.perf_counter() - t_section:.3f}s", flush=True)

    controls_column = pn.Column(*controls_children, sizing_mode="fixed", width=400, css_classes=["left-col"])

    peruse_button.visible = False  # Initially hidden
    main_placeholder = pn.Column(
        pn.pane.HTML(placeholder_text),
        sizing_mode="stretch_both",
        css_classes=["main-right"],
    )

    # Wrap everything in a Flex container
    layout = pn.Row(controls_column, main_placeholder, sizing_mode="stretch_both", css_classes=["main-layout"])
    layout.stylesheets = [stylesheet]

    if enable_timing:
        print(f"[timing] UI construction (widgets + layout): {time.perf_counter() - t_ui:.3f}s", flush=True)
        print(f"[timing] Total initial load: {time.perf_counter() - t_init:.3f}s", flush=True)

    return layout

def add_serve_args(parser):
    parser.add_argument("--db", required=True, help="Path to DuckDB database")
    parser.add_argument("--port", type=int, default=5006, help="Port to serve Panel app")
    parser.add_argument('--time', action='store_true', default=False, help=argparse.SUPPRESS)

def run_serve(args):
    # Print database metadata if available
    import duckdb as _duckdb
    _conn = _duckdb.connect(args.db, read_only=True)
    try:
        rows = _conn.execute("SELECT Key, Value FROM Database_metadata").fetchall()
        meta = dict(rows)
        db_name = os.path.basename(args.db)
        params = []
        for key in ['Modules', 'Min_aligned_fraction', 'Min_coverage_depth',
                     'Coverage_percentage']:
            if key in meta:
                params.append(f"{key}={meta[key]}")
        print(f"Database '{db_name}': "
              f"\n Created on {meta.get('Date_of_creation', '?')} "
              f"(v{meta.get('Tool_version_used_for_creation', '?')}) "
              f"\n Last modified on {meta.get('Date_of_last_modification', '?')} "
              f"(v{meta.get('Tool_version_used_for_last_modification', '?')})")
        if params:
            params_str = '\n '.join(params)
            print(f"Calculate parameters used:\n {params_str}")
    finally:
        _conn.close()

    # Create a factory function that Panel will call for each session
    enable_timing = getattr(args, 'time', False)
    def create_app():
        return create_layout(args.db, enable_timing=enable_timing)

    static_path = os.path.join(os.path.dirname(__file__), "..", "static")
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