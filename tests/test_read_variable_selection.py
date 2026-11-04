"""ALL SAMPLES switches between signal selection and read-class selection."""

from types import SimpleNamespace

from bokeh.models import CheckboxButtonGroup

from thebigbam.plotting.application.read_track_pipeline import read_track_config
from thebigbam.plotting.application.scope_transition import ScopeTransitionController
from thebigbam.plotting.models.read_plot import ReadPlotOptions


def test_all_sample_read_toggle_allows_multiple_variables_then_clears_them():
    groups = [CheckboxButtonGroup(labels=["Primary", "Other"]), CheckboxButtonGroup(labels=["Clipping"])]
    reads = CheckboxButtonGroup(labels=["Read plot"])
    select = SimpleNamespace(value="sample", options=["sample"])
    bindings = SimpleNamespace(
        interaction_lock={"locked": False}, sample_scope=SimpleNamespace(active=1),
        widgets={"variables_widgets_all": groups, "read_plot_controls": SimpleNamespace(all_reads=reads),
                 "has_mags": False, "view_radio": SimpleNamespace(active=1),
                 "sample_select": select, "contig_select": select},
    )
    controller = ScopeTransitionController(bindings)
    for group in groups:
        group.on_change("active", controller.variable_callback(group))
    reads.on_change("active", controller.read_plot_changed)

    groups[0].active = [0]
    groups[1].active = [0]
    assert groups[0].active == []
    reads.active = [0]
    groups[0].active = [0, 1]
    assert groups[1].active == [0]
    assert groups[0].active == [0, 1]
    reads.active = []
    assert all(group.active == [] for group in groups)
    groups[0].active = [0]
    groups[0].active = [0, 1]
    assert groups[0].active == [1]
    groups[1].active = [0]
    assert groups[0].active == []


def test_all_sample_read_classes_follow_selected_variables():
    widgets = {"bam_catalog": object(), "read_colors": {}, "subplot_to_varnames": {
        "Primary": ["primary_reads"], "Other": ["secondary_reads", "supplementary_reads"],
    }}
    inputs = SimpleNamespace(is_all=True, read_plot=ReadPlotOptions(enabled=True, selected_subplots=("Other",)))
    config = read_track_config(inputs, widgets)
    assert config.classes == frozenset({"secondary", "supplementary"})
    inputs.read_plot = ReadPlotOptions(enabled=True)
    assert read_track_config(inputs, widgets).classes == frozenset()
