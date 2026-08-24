from dataclasses import fields
from pathlib import Path
from types import SimpleNamespace

from thebigbam.plotting.application.apply_controller import ApplyBindings, ApplyController


class Diagnostics:
    enabled = False

    def next_generation(self):
        pass


class Lifecycle:
    def __init__(self):
        self.prepared = False

    def prepare_replacement(self):
        self.prepared = True


def test_apply_controller_reports_missing_contig_and_completes_lifecycle():
    values = {field.name: None for field in fields(ApplyBindings)}
    placeholder = SimpleNamespace(objects=[object()], loading=True)
    lifecycle = Lifecycle()
    operations = []
    spinner = SimpleNamespace(value=100)
    values.update(
        widgets={
            "contig_select": SimpleNamespace(value=""),
            "has_samples": False,
            "has_mags": False,
            "view_radio": SimpleNamespace(active=1),
            "variables_widgets_one": [],
            "variables_widgets_all": [],
        },
        views=SimpleNamespace(active=0),
        custom_color_rows=[],
        mag_track_color_rows=[],
        mag_track_max_dots_input=spinner,
        max_genemap_window_input=spinner,
        same_y_scale_cbg=SimpleNamespace(active=[]),
        subplot_height_input=spinner,
        genemap_height_input=spinner,
        sequence_height_input=spinner,
        translated_sequence_height_input=spinner,
        max_binning_window_input=spinner,
        min_coverage_freq_input=SimpleNamespace(value=0),
        diagnostics=Diagnostics(),
        plot_lifecycle=lifecycle,
        main_placeholder=placeholder,
        peruse_button=SimpleNamespace(visible=True),
        download_metrics_button=SimpleNamespace(visible=True),
        download_mag_metrics_button=SimpleNamespace(visible=True),
        download_data_button=SimpleNamespace(visible=True),
        command_hint_pane=SimpleNamespace(visible=True),
        enable_timing=False,
        set_operation=operations.append,
    )

    ApplyController(ApplyBindings(**values)).apply()

    # Invalid requests leave the active plot and its callbacks untouched.
    assert not lifecycle.prepared
    assert "Please select a contig" in placeholder.objects[0].object
    assert placeholder.loading is False
    assert operations == ["apply/param_parse", "idle"]


def test_apply_controller_does_not_clear_current_plot_before_replacement():
    source = Path("thebigbam/plotting/application/apply_controller.py").read_text(encoding="utf-8")

    assert "main_placeholder.objects = []" not in source


def test_apply_controller_is_thin_and_registers_four_typed_handlers():
    source = Path("thebigbam/plotting/application/apply_controller.py").read_text(encoding="utf-8")

    assert len(source.splitlines()) <= 100
    assert ".execute(" not in source
    assert "compose_" not in source
    for handler_name in ("ContigOneHandler", "ContigAllHandler", "MagOneHandler", "MagAllHandler"):
        assert handler_name in source
