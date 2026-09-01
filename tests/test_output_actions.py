from types import SimpleNamespace

import panel as pn

from thebigbam.plotting.application import output_controls as module


def test_output_controls_record_summary_download_and_command_actions(monkeypatch):
    monkeypatch.setattr(module, "make_contig_metrics_callback", lambda *args: lambda: "contig")
    monkeypatch.setattr(module, "make_mag_metrics_callback", lambda *args: lambda: "mag")

    class FakeInspect:
        def __init__(self, _bindings):
            self.button = pn.widgets.Button(name="DOWNLOAD DATA")
            self.pane = pn.pane.HTML("")

    monkeypatch.setattr(module, "InspectCommandController", FakeInspect)
    actions = []
    summaries = []
    controls = module.build_output_controls(
        db_path="example.db",
        connection=object(),
        widgets={},
        sample_scope=object(),
        filtered_samples=lambda _contig: [],
        combined_features=None,
        subplot_variables={},
        from_position=SimpleNamespace(value="1"),
        to_position=SimpleNamespace(value="2"),
        stylesheet="",
        enable_timing=False,
        timing=None,
        report_timing=None,
        show_summary=lambda _event: summaries.append(True),
        record_action=lambda action, details: actions.append((action, details)),
    )

    controls.peruse_button.clicks += 1
    assert controls.contig_metrics_button.callback() == "contig"
    assert controls.mag_metrics_button.callback() == "mag"
    controls.data_button.clicks += 1

    assert summaries == [True]
    assert actions == [
        ("show_summary", {}),
        ("download_contig_metrics", {}),
        ("download_mag_metrics", {}),
        ("show_download_command", {}),
    ]
