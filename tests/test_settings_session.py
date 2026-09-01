from types import SimpleNamespace

import panel as pn

from thebigbam.plotting.application import settings_session as module


def test_save_settings_keeps_draft_collector_while_history_uses_applied(monkeypatch):
    collector = SimpleNamespace(
        collect=lambda: {"filtering": ["draft"]},
        collect_applied=lambda _expression: {"filtering": ["committed"]},
    )
    captured = {}

    class FakeSaveControls:
        def __init__(self, _db_path, collect, _stylesheet, _record_action=None):
            captured["collect"] = collect
            self.button = pn.widgets.Button(name="SAVE SETTINGS")
            self.confirmation = pn.pane.Markdown("")

    monkeypatch.setattr(module, "SettingsCollector", lambda _bindings: collector)
    monkeypatch.setattr(module, "SettingsSaveControls", FakeSaveControls)
    monkeypatch.setattr(module, "make_settings_collector_bindings", lambda *args: object())
    monkeypatch.setattr(module, "make_settings_restore_bindings", lambda **kwargs: object())
    projection = SimpleNamespace(applied_expression=lambda: "committed-expression")

    session = module.build_settings_session(
        db_path="example.db",
        widgets={},
        sample_scope=object(),
        genome=object(),
        parameters=object(),
        filtering=object(),
        filtering_metadata={},
        create_query_row=lambda: None,
        apply_button=pn.widgets.Button(name="APPLY"),
        filter_projection=projection,
        stylesheet="",
    )

    assert captured["collect"]() == {"filtering": ["draft"]}
    assert session.collect_applied() == {"filtering": ["committed"]}
    assert [item.name for item in session.buttons_row] == ["APPLY"]
