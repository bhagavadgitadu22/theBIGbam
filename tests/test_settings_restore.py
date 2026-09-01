from dataclasses import fields
from types import SimpleNamespace

from thebigbam.plotting.settings.restoration import SettingsRestoreBindings, restore_settings


def _widget(value=None, options=None, active=0):
    return SimpleNamespace(value=value, options=[] if options is None else options, active=active)


def test_restore_settings_applies_view_and_valid_subject_selection():
    sample = _widget("old", ["old", "new"])
    widgets = {
        "view_radio": _widget(active=1),
        "sample_select": sample,
        "contig_select": _widget("c1", ["c1"]),
        "mag_select": _widget("m1", ["m1"]),
        "module_names": [],
        "module_widgets_one": [],
        "variables_widgets_one": [],
        "variables_widgets_all": [],
    }
    values = {field.name: None for field in fields(SettingsRestoreBindings)}
    values.update(
        widgets=widgets,
        views=_widget(active=0),
        or_sections=[],
        inter_section_selects=[],
        custom_color_rows=[],
        mag_track_color_rows=[],
        filtering_metadata={},
        color_qualifier_options=[],
        mag_params_category_select=_widget(options=[]),
        mag_params_metric_select=_widget(options=[]),
        sample_order_category_select=_widget(options=[]),
        sample_order_metric_select=_widget(options=[]),
    )
    bindings = SettingsRestoreBindings(**values)

    restore_settings(
        {
            "view_mode": {"mag_or_contig": 0, "one_or_all_samples": 1},
            "selection": {"sample": "new", "contig": "missing"},
        },
        bindings,
    )

    assert widgets["view_radio"].active == 0
    assert bindings.views.active == 1
    assert sample.value == "new"
    assert widgets["contig_select"].value == "c1"


def test_explicit_empty_filtering_clears_existing_filter_builder():
    widgets = {
        "view_radio": _widget(active=0),
        "sample_select": _widget("", []),
        "contig_select": _widget("", []),
        "mag_select": _widget("", []),
        "module_names": [],
        "module_widgets_one": [],
        "variables_widgets_one": [],
        "variables_widgets_all": [],
    }
    existing_sections = [{"rows": ["old"]}]
    rebuilt = []
    values = {field.name: None for field in fields(SettingsRestoreBindings)}
    values.update(
        widgets=widgets,
        views=_widget(active=0),
        or_sections=existing_sections,
        inter_section_selects=[],
        custom_color_rows=[],
        mag_track_color_rows=[],
        filtering_metadata={},
        color_qualifier_options=[],
        create_or_section=lambda: {"rows": ["empty"]},
        rebuild_filtering_content=lambda: rebuilt.append(True),
        mag_params_category_select=_widget(options=[]),
        mag_params_metric_select=_widget(options=[]),
        sample_order_category_select=_widget(options=[]),
        sample_order_metric_select=_widget(options=[]),
    )

    restore_settings({"filtering": []}, SettingsRestoreBindings(**values))

    assert existing_sections == [{"rows": ["empty"]}]
    assert rebuilt == [True]
