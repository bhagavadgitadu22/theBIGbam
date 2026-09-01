import datetime as dt
import io
import json
from types import SimpleNamespace

from thebigbam.plotting.downloads import callbacks as download_callbacks
from thebigbam.plotting.downloads.callbacks import make_contig_metrics_callback, make_mag_metrics_callback
from thebigbam.plotting.settings.controls import SettingsSaveControls
from thebigbam.plotting.settings.persistence import (
    load_settings_document,
    save_session_document,
    save_settings_document,
    serialize_color_rows,
    serialize_filter_sections,
    serialize_variable_selection,
)
from thebigbam.plotting.shared.timing import BrowserTimingRelay


def test_settings_serializers_use_stable_names():
    widgets = {
        "module_names": ["coverage"],
        "module_widgets_one": [SimpleNamespace(active=[0])],
        "variables_widgets_one": [SimpleNamespace(labels=["depth", "breadth"], active=[1])],
        "variables_widgets_all": [SimpleNamespace(labels=["depth", "breadth"], active=[0])],
    }
    row = {
        "category_select": SimpleNamespace(value="Contigs"),
        "subcategory_select": SimpleNamespace(value="Length"),
        "comparison_select": SimpleNamespace(value=">"),
        "input_ref": {"widget": SimpleNamespace(value=1000)},
        "and_div": None,
    }
    color_row = {
        "qualifier_select": SimpleNamespace(value="product"),
        "operator_select": SimpleNamespace(value="contains"),
        "input_ref": {"widget": SimpleNamespace(value="kinase")},
        "color_picker": SimpleNamespace(color="#abcdef"),
    }

    assert serialize_variable_selection(widgets) == {
        "coverage": {
            "module_enabled": True,
            "selected_one": ["breadth"],
            "selected_all": ["depth"],
        }
    }
    assert serialize_filter_sections([{"rows": [row]}], []) == [
        {
            "section_and_or": None,
            "rows": [
                {
                    "category": "Contigs",
                    "column": "Length",
                    "operator": ">",
                    "value": 1000,
                    "row_and_or": None,
                }
            ],
        }
    ]
    assert serialize_color_rows([color_row])[0]["color"] == "#abcdef"


def test_settings_document_round_trip(tmp_path):
    def fixed_time():
        return dt.datetime(2026, 8, 18, 12, 34, 56)

    path = save_settings_document(
        {"view_mode": {"one_or_all_samples": 1}},
        "example.db",
        tmp_path,
        fixed_time,
    )

    assert path.name == "example_settings_20260818_123456.json"
    assert load_settings_document(path) == {"view_mode": {"one_or_all_samples": 1}}
    assert json.loads(path.read_text(encoding="utf-8"))["view_mode"]


def test_session_document_uses_explicit_timestamped_name(tmp_path):
    path = save_session_document(
        {"format": "thebigbam-session-history"},
        "example.db",
        tmp_path,
        lambda: dt.datetime(2026, 8, 18, 12, 34, 56),
    )

    assert path.name == "example_session_20260818_123456.json"
    assert json.loads(path.read_text(encoding="utf-8"))["format"] == "thebigbam-session-history"


def test_settings_document_rejects_non_object(tmp_path):
    path = tmp_path / "bad.json"
    path.write_text("[]", encoding="utf-8")

    try:
        load_settings_document(path)
    except ValueError as error:
        assert "object" in str(error)
    else:
        raise AssertionError("non-object settings should be rejected")


def test_settings_save_controls_persist_collected_document(monkeypatch, tmp_path):
    saved = []

    def fake_save(document, db_path):
        saved.append((document, db_path))
        return tmp_path / "saved.json"

    monkeypatch.setattr("thebigbam.plotting.settings.controls.save_settings_document", fake_save)
    controls = SettingsSaveControls("example.db", lambda: {"selection": {"contig": "c1"}}, "")

    controls._save()

    assert saved == [({"selection": {"contig": "c1"}}, "example.db")]
    assert controls.confirmation.text == "1"


def test_settings_save_records_explicit_scenario_action(monkeypatch, tmp_path):
    monkeypatch.setattr(
        "thebigbam.plotting.settings.controls.save_settings_document",
        lambda *_args: tmp_path / "saved.json",
    )
    actions = []
    controls = SettingsSaveControls(
        "example.db",
        dict,
        "",
        lambda action, details: actions.append((action, details)),
    )

    controls._save()
    controls.button.clicks += 1

    assert actions == [("save_settings", {})]


class Value:
    def __init__(self, value="", active=0):
        self.value = value
        self.active = active


class Download:
    filename = ""


def test_contig_download_callback_uses_scope_and_updates_filename(monkeypatch):
    monkeypatch.setattr(
        download_callbacks, "download_contig_metrics_csv", lambda _db, contig, samples: f"{contig}:{','.join(samples)}"
    )
    widgets = {"contig_select": Value("c 1"), "sample_select": Value("s1")}
    downloads = {"contig_metrics": Download()}
    callback = make_contig_metrics_callback("db", widgets, Value(active=1), lambda _contig: ["s1", "s2"], downloads)

    result = callback()

    assert isinstance(result, io.StringIO)
    assert result.getvalue() == "c 1:s1,s2"
    assert downloads["contig_metrics"].filename == "c_1_in_all_samples_contig_metrics.csv"


def test_mag_download_callback_resolves_mag_from_contig(monkeypatch):
    monkeypatch.setattr(
        download_callbacks, "download_mag_metrics_csv", lambda _db, mag, samples: f"{mag}:{','.join(samples)}"
    )
    widgets = {
        "has_mags": True,
        "view_radio": Value(active=1),
        "mag_select": Value("m1"),
        "contig_select": Value("c1"),
        "contig_to_mag": {"c1": "m from contig"},
        "sample_select": Value("s1"),
    }
    downloads = {"mag_metrics": Download()}

    result = make_mag_metrics_callback("db", widgets, downloads)()

    assert result.getvalue() == "m from contig:s1"
    assert downloads["mag_metrics"].filename == "m_from_contig_in_s1_mag_metrics.csv"


def test_browser_timing_relay_is_inert_when_disabled():
    relay = BrowserTimingRelay(False)
    relay.send("view_change")
    assert relay.state == {}
    assert relay.ping.value == ""
