import argparse
import json

from thebigbam.plotting.application.server import add_serve_args
from thebigbam.plotting.settings.scenario import ScenarioPathAllocator, ScenarioRecorder


def _settings(*, sample="s1", saved_at="first"):
    return {
        "_meta": {"saved_at": saved_at},
        "view_mode": {"mag_or_contig": 0, "one_or_all_samples": 0},
        "selection": {"sample": sample, "contig": "c1", "mag": None},
        "filtering": [],
    }


def test_recorder_keeps_settings_ontology_and_ignores_volatile_metadata(tmp_path):
    path = tmp_path / "scenario.json"
    recorder = ScenarioRecorder(path, "example.db", _settings())
    recorder.flush()

    assert recorder.record_state(_settings(saved_at="second")) is False
    assert recorder.record_state(_settings(sample="s2")) is True
    recorder.close()

    document = json.loads(path.read_text(encoding="utf-8"))
    assert document["_meta"]["format"] == "thebigbam-scenario"
    assert "complete" not in document["_meta"]
    assert "started_at" not in document["_meta"]
    assert document["initial_state"]["selection"]["sample"] == "s1"
    assert document["steps"] == [
        {
            "sequence": 1,
            "action": "state_change",
            "changes": {"selection": {"sample": "s2"}},
        }
    ]
    assert document["final_state"]["selection"]["sample"] == "s2"


def test_atomic_failure_preserves_previous_valid_document(tmp_path, monkeypatch):
    path = tmp_path / "scenario.json"
    recorder = ScenarioRecorder(path, "example.db", _settings())
    recorder.flush()
    before = json.loads(path.read_text(encoding="utf-8"))

    def fail_replace(_source, _destination):
        raise OSError("interrupted")

    monkeypatch.setattr("thebigbam.plotting.settings.scenario.os.replace", fail_replace)
    recorder.record_state(_settings(sample="s2"))
    recorder.flush()

    assert json.loads(path.read_text(encoding="utf-8")) == before
    assert recorder.error is not None
    recorder.close()


def test_allocator_avoids_concurrent_session_overwrites(tmp_path):
    allocator = ScenarioPathAllocator(tmp_path / "workflow.json")

    assert allocator.next_path() == tmp_path / "workflow.json"
    assert allocator.next_path() == tmp_path / "workflow.session-2.json"


def test_action_is_recorded_even_when_settings_do_not_change(tmp_path):
    path = tmp_path / "scenario.json"
    recorder = ScenarioRecorder(path, "example.db", _settings())

    apply_step = recorder.record_action("apply_plot", _settings())
    recorder.close()

    document = json.loads(path.read_text(encoding="utf-8"))
    assert document["steps"][0]["action"] == "apply_plot"
    assert apply_step == document["steps"][0]["sequence"]
    assert "changes" not in document["steps"][0]


def test_semantic_action_records_intent_without_result_timing_or_status(tmp_path):
    path = tmp_path / "scenario.json"
    recorder = ScenarioRecorder(path, "example.db", _settings())

    recorder.record_action("apply_filters", _settings())
    recorder.close()

    step = json.loads(path.read_text(encoding="utf-8"))["steps"][0]
    assert step["action"] == "apply_filters"
    assert "status" not in step
    assert "duration_seconds" not in step


def test_semantic_action_preserves_json_safe_details(tmp_path):
    path = tmp_path / "scenario.json"
    recorder = ScenarioRecorder(path, "example.db", _settings())

    recorder.record_action(
        "filter_lookup",
        details={"category": "Coverage", "column": "Coverage_mean"},
    )
    recorder.close()

    step = json.loads(path.read_text(encoding="utf-8"))["steps"][0]
    assert step["details"] == {"category": "Coverage", "column": "Coverage_mean"}


def test_scenario_option_is_documented_after_json_in_help():
    parser = argparse.ArgumentParser()
    add_serve_args(parser)

    args = parser.parse_args(["--db", "example.db", "--scenario", "workflow.json"])

    assert args.scenario == "workflow.json"
    help_text = parser.format_help()
    assert "--scenario" in help_text
    assert "continuously recorded" in help_text
    assert help_text.index("--json") < help_text.index("--no-browser") < help_text.index("--scenario")
    assert "(for developers) Path where user actions" in help_text
