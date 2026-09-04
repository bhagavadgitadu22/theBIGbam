import panel as pn
from bokeh.models import Div

from thebigbam.plotting.application.layout import LayoutParts, assemble_layout, separator


def _parts():
    values = {name: Div(text=name) for name in LayoutParts.__dataclass_fields__}
    return LayoutParts(**values)


def _history():
    return pn.Column(width=250)


def _bokeh_models(controls):
    models = []

    def visit(item):
        if hasattr(item, "object"):
            models.append(item.object)
        for child in getattr(item, "objects", ()):  # Panel containers
            visit(child)

    visit(controls)
    return models


def test_sample_layout_contains_scope_and_sample_specific_sections():
    parts = _parts()
    assembled = assemble_layout(
        parts, has_samples=True, summary_carrier=Div(), stylesheet="", history_drawer=_history()
    )

    models = _bokeh_models(assembled.controls)
    assert parts.sample_scope in models
    assert parts.sample_section in models
    assert parts.variables_all in models
    assert assembled.controls.objects[0].css_classes == ["sidebar-content"]
    bottom_space = assembled.controls.objects[0].objects[-1]
    assert bottom_space.height == 20
    assert bottom_space.sizing_mode == "fixed"
    assert bottom_space.css_classes == ["sidebar-bottom-space"]
    assert isinstance(assembled.layout, pn.Row)
    assert assembled.controls.visible is True
    assert assembled.history.visible is False


def test_genome_only_layout_omits_sample_specific_sections():
    parts = _parts()
    assembled = assemble_layout(
        parts, has_samples=False, summary_carrier=Div(), stylesheet="", history_drawer=_history()
    )

    models = _bokeh_models(assembled.controls)
    assert parts.sample_scope not in models
    assert parts.sample_section not in models
    assert parts.variables_one not in models


def test_drawer_toggles_stay_on_real_layout_rails():
    assembled = assemble_layout(
        _parts(), has_samples=True, summary_carrier=Div(), stylesheet="", history_drawer=_history()
    )
    left_toggle, right_toggle = assembled.left_toggle, assembled.right_toggle

    assert "drawer-toggle-left" in left_toggle.css_classes
    assert "drawer-toggle-right" in right_toggle.css_classes
    assert left_toggle.margin == 0
    assert right_toggle.margin == 0
    assert assembled.layout.objects[0] is assembled.left_shell
    assert assembled.layout.objects[3] is assembled.right_shell
    assert assembled.left_shell.objects[-1] is assembled.left_rail
    assert assembled.right_shell.objects[-1] is assembled.right_rail
    assert assembled.left_rail.objects[-1] is left_toggle
    assert assembled.right_rail.objects[-1] is right_toggle
    assert assembled.placeholder is assembled.layout.objects[1]
    assert assembled.placeholder.css_classes == ["main-right", "plot-area"]

    left_toggle.clicks += 1
    right_toggle.clicks += 1

    assert assembled.controls.visible is False
    assert assembled.history.visible is True
    assert assembled.left_resizer.enabled is False
    assert assembled.right_resizer.enabled is True
    assert assembled.left_shell.width == 0
    assert assembled.right_shell.width == 250


def test_drawer_resizers_are_client_side_bounded_and_persistent():
    assembled = assemble_layout(
        _parts(), has_samples=True, summary_carrier=Div(), stylesheet="", history_drawer=_history()
    )

    assert (assembled.left_resizer.minimum, assembled.left_resizer.maximum) == (100, 700)
    assert assembled.left_resizer.default_width == 400
    assert assembled.controls.width == 400
    assert (assembled.right_resizer.minimum, assembled.right_resizer.maximum) == (100, 700)
    assert assembled.right_resizer.default_width == 250
    assert assembled.history.width == 250
    assert assembled.left_rail.width == 8
    assert assembled.right_rail.width == 8
    assert assembled.left_rail.styles["right"] == "-8px"
    assert assembled.right_rail.styles["left"] == "-8px"
    assert assembled.left_shell.width == 400
    assert assembled.right_shell.width == 0
    assert "localStorage.setItem" in assembled.left_resizer._esm
    assert "dblclick" in assembled.left_resizer._esm

    assembled.left_resizer.value = 510
    assembled.right_resizer.value = 310
    assert assembled.controls.width == 510
    assert assembled.history.width == 310
    assert assembled.left_shell.width == 510
    assert assembled.right_shell.width == 0


def test_shared_separator_has_fixed_thickness_and_explicit_spacing():
    divider = separator()
    assert divider.margin == (10, 0)
    assert divider.height == 2
    assert divider.min_height == 2
    assert divider.max_height == 2
    assert divider.css_classes == ["section-separator"]
    assert divider.styles == {"background-color": "#333", "box-sizing": "border-box"}
