import json

import pytest

from thebigbam.plotting.application.composition_root import _build_scenario_action_carrier
from thebigbam.plotting.benchmark import replay
from thebigbam.plotting.settings.scenario import describe_scenario_document


def _document(step):
    return {
        "_meta": {"format": "thebigbam-scenario", "version": 1, "source_db": "example.db"},
        "initial_state": {},
        "steps": [step],
    }


def test_history_scenario_actions_validate_and_replay_semantically(tmp_path, monkeypatch):
    details = {
        "history_sequence": 2,
        "history_action": "apply_plot",
        "apply_step": 1,
        "settings": {"selection": {"contig": "c1"}},
    }
    path = tmp_path / "scenario.json"
    path.write_text(
        json.dumps(
            {
                **_document({"sequence": 2, "action": "restore_history", "details": details}),
                "steps": [
                    {"sequence": 1, "action": "apply_plot"},
                    {"sequence": 2, "action": "restore_history", "details": details},
                ],
            }
        ),
        encoding="utf-8",
    )
    actions = []
    monkeypatch.setattr(
        replay,
        "perform_semantic_action",
        lambda page, action, action_details, timeout: actions.append((action, action_details)),
    )

    step = replay.load_scenario(path)["steps"][1]
    replay.replay_action(object(), step, {}, 1000)

    assert actions == [("restore_history", details)]


@pytest.mark.parametrize(
    "details",
    [
        {"history_sequence": 0, "history_action": "apply_plot", "settings": {}},
        {"history_sequence": 1, "history_action": "unknown", "settings": {}},
        {"history_sequence": 1, "history_action": "apply_plot", "settings": []},
        {
            "history_sequence": 1,
            "history_action": "apply_plot",
            "apply_step": 99,
            "settings": {},
        },
    ],
)
def test_restore_history_scenario_rejects_invalid_details(tmp_path, details):
    path = tmp_path / "scenario.json"
    path.write_text(
        json.dumps(_document({"sequence": 1, "action": "restore_history", "details": details})),
        encoding="utf-8",
    )

    with pytest.raises(replay.ScenarioError):
        replay.load_scenario(path)


def test_async_semantic_action_acknowledges_only_after_completion():
    completions = []
    carrier, status = _build_scenario_action_carrier(
        lambda _action, _details, complete: completions.append(complete) or True
    )
    carrier.value = json.dumps(
        {"nonce": "restore-1", "action": "restore_history", "details": {}}
    )

    assert status.value == ""
    completions[0]()
    assert json.loads(status.value) == {"nonce": "restore-1", "status": "completed"}


def test_history_scenario_actions_have_explicit_descriptions():
    document = _document({"sequence": 1, "action": "apply_filters"})
    document["steps"].append(
        {
            "sequence": 2,
            "action": "remove_history",
            "details": {
                "history_sequence": 4,
                "history_action": "apply_filters",
                "apply_step": 1,
            },
        }
    )

    assert describe_scenario_document(document) == (
        "1. Apply filters",
        "2. Remove apply filters from step 1",
    )


def test_saved_session_history_action_needs_no_synthetic_apply_step(tmp_path):
    details = {
        "history_sequence": 4,
        "history_action": "apply_plot",
        "apply_step": None,
        "settings": {"selection": {"contig": "c1"}},
    }
    document = _document(
        {"sequence": 1, "action": "restore_history", "details": details}
    )
    path = tmp_path / "saved-session-restore.json"
    path.write_text(json.dumps(document), encoding="utf-8")

    assert replay.load_scenario(path)["steps"][0]["details"] == details
    assert describe_scenario_document(document) == (
        "1. Restore saved-session apply plot",
    )

@pytest.mark.parametrize(
    ("action", "label"),
    [
        ("show_summary", "Show summary"),
        ("download_contig_metrics", "Download contig metrics"),
        ("download_mag_metrics", "Download MAG metrics"),
        ("show_download_command", "Show download command"),
        ("save_settings", "Save settings"),
        ("save_session", "Save session"),
        ("reset_position", "Reset genomic position"),
    ],
)
def test_output_scenario_actions_have_explicit_descriptions(action, label):
    document = _document({"sequence": 1, "action": action, "details": {}})
    assert describe_scenario_document(document) == (f"1. {label}",)
