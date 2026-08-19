from bokeh.models.widgets import CheckboxButtonGroup, CheckboxGroup

from thebigbam.plotting.controls.variables import attach_module_synchronization, build_variable_controls


def test_variable_controls_separate_genome_and_synchronize_module_selection():
    module = CheckboxGroup(labels=["Coverage"], active=[])
    variables = CheckboxButtonGroup(labels=["Depth", "Breadth"], active=[])
    genome = CheckboxButtonGroup(labels=["GC content"], active=[])
    widgets = {
        "module_names": ["Coverage", "Genome"],
        "module_widgets_one": [module, CheckboxGroup(labels=["Genome"], active=[])],
        "variables_widgets_one": [variables, genome],
        "variables_widgets_all": [
            CheckboxButtonGroup(labels=["Depth", "Breadth"], active=[]),
            CheckboxButtonGroup(labels=["GC content"], active=[]),
        ],
        "helps_widgets": [None, None],
    }

    controls = build_variable_controls(widgets, "", lambda button, content: lambda event: None)
    attach_module_synchronization(widgets, {"locked": False})
    module.active = [0]

    assert controls.genome_checkbox is genome
    assert controls.genome_index == 1
    assert variables.active == [0, 1]
    assert controls.one_sample.visible
    assert not controls.all_samples.visible
