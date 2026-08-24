import panel as pn
from bokeh.models import Div

from thebigbam.plotting.application.layout import LayoutParts, assemble_layout


def _parts():
    values = {name: Div(text=name) for name in LayoutParts.__dataclass_fields__}
    return LayoutParts(**values)


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
        parts, has_samples=True, summary_carrier=Div(), stylesheet=""
    )

    models = _bokeh_models(assembled.controls)
    assert parts.sample_scope in models
    assert parts.sample_section in models
    assert parts.variables_all in models
    assert assembled.controls.objects[0].css_classes == ["sidebar-content"]
    assert isinstance(assembled.layout, pn.Row)


def test_genome_only_layout_omits_sample_specific_sections():
    parts = _parts()
    assembled = assemble_layout(
        parts, has_samples=False, summary_carrier=Div(), stylesheet=""
    )

    models = _bokeh_models(assembled.controls)
    assert parts.sample_scope not in models
    assert parts.sample_section not in models
    assert parts.variables_one not in models
