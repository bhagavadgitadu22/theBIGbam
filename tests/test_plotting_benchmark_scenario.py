import json

import pytest

from benchmarks.plotting.run_workflow import (
    ScenarioError,
    StepResult,
    augment_scenario,
    choose_filter_lookup_index,
    load_scenario,
    merge_changes,
)
from thebigbam.plotting.application.composition_root import _build_scenario_restore_carrier


def _scenario():
    return {
        "_meta": {"format": "thebigbam-scenario", "version": 1, "source_db": "example.db"},
        "initial_state": {"selection": {"sample": "S1"}, "filtering": []},
        "steps": [
            {
                "sequence": 4,
                "action": "state_change",
                "changes": {"selection": {"sample": "S2"}},
            },
            {"sequence": 7, "action": "apply_plot"},
        ],
    }


def test_load_scenario_validates_database_and_replay_contract(tmp_path):
    path = tmp_path / "scenario.json"
    path.write_text(json.dumps(_scenario()), encoding="utf-8")

    assert load_scenario(path, tmp_path / "example.db")["steps"][0]["sequence"] == 4
    with pytest.raises(ScenarioError, match="expects database"):
        load_scenario(path, tmp_path / "other.db")


def test_preflight_accepts_complete_settings_collector_ontology(tmp_path):
    document = _scenario()
    document["steps"] = [
        {
            "sequence": 1,
            "action": "state_change",
            "changes": {
                "view_mode": {"mag_or_contig": 0, "one_or_all_samples": 1},
                "selection": {"sample": "S1", "contig": "C1", "mag": "M1"},
                "filtering": [{"section_and_or": None, "rows": []}],
                "variables": {
                    "Coverage": {"module_enabled": True, "selected_one": [], "selected_all": ["Depth"]}
                },
                "contig": {
                    "position_range": {"from": "1", "to": "100"},
                    "feature_widgets": {"feature_types": ["CDS"], "sequence": True},
                    "coloring": {
                        "custom_color_rows": [],
                        "mag_track_color_rows": [],
                        "apply_annotation_rules_to_mag_track": False,
                    },
                },
                "plotting_params": {
                    "min_coverage_freq": 0.1,
                    "max_genemap_window": 1000,
                    "max_sequence_window": 500,
                    "max_binning_window": 100,
                    "genemap_height": 100,
                    "sequence_height": 50,
                    "translated_sequence_height": 50,
                    "subplot_height": 100,
                    "mag_params": {
                        "category": "Coverage",
                        "metric": "Coverage_mean",
                        "direction": 1,
                        "sort_sample": "S1",
                        "max_dots": 1000,
                    },
                    "sample_params": {
                        "max_samples": 20,
                        "order_category": "Coverage",
                        "order_metric": "Coverage_mean",
                        "order_direction": 0,
                        "same_y_scale": True,
                    },
                },
            },
        }
    ]
    path = tmp_path / "scenario.json"
    path.write_text(json.dumps(document), encoding="utf-8")

    assert load_scenario(path)["steps"][0]["sequence"] == 1


@pytest.mark.parametrize(
    ("change", "message"),
    [
        ({"steps": []}, "at least one step"),
        ({"steps": [{"sequence": 1, "action": "copy_contig"}]}, "unsupported action"),
        ({"steps": [{"sequence": 1, "action": "state_change"}]}, "requires.*changes"),
        (
            {
                "steps": [
                    {"sequence": 2, "action": "apply_plot"},
                    {"sequence": 1, "action": "apply_filters"},
                ]
            },
            "not in increasing order",
        ),
        (
            {
                "steps": [
                    {
                        "sequence": 1,
                        "action": "state_change",
                        "changes": {"unknown": {"value": 1}},
                    }
                ]
            },
            "unsupported settings section",
        ),
    ],
)
def test_load_scenario_rejects_non_replayable_documents(tmp_path, change, message):
    document = _scenario() | change
    path = tmp_path / "scenario.json"
    path.write_text(json.dumps(document), encoding="utf-8")

    with pytest.raises(ScenarioError, match=message):
        load_scenario(path)


def test_merge_changes_reconstructs_settings_state_without_mutating_source():
    state = {"selection": {"sample": "S1", "mag": "M1"}, "variables": ["Coverage"]}

    merged = merge_changes(state, {"selection": {"sample": "S2", "mag": None}})

    assert merged == {"selection": {"sample": "S2"}, "variables": ["Coverage"]}
    assert state["selection"] == {"sample": "S1", "mag": "M1"}


def test_result_keeps_scenario_shape_and_augments_matching_numbered_steps():
    scenario = _scenario()
    results = [
        StepResult(4, "state_change", "completed", 0.125),
        StepResult(7, "apply_plot", "failed", 4.218, error="Timeout"),
    ]

    result = augment_scenario(scenario, results, {"run": "cold"})

    assert scenario["steps"][0].get("status") is None
    assert result["initial_state"] == scenario["initial_state"]
    assert result["steps"][0]["duration_seconds"] == 0.125
    assert result["steps"][1]["status"] == "failed"
    assert result["steps"][1]["error"] == "Timeout"
    assert result["_result"] == {"run": "cold"}


def test_hidden_restore_carrier_uses_settings_restoration_boundary():
    restored = []

    class Session:
        def restore(self, settings):
            restored.append(settings)

    carrier, status = _build_scenario_restore_carrier(Session())

    assert carrier.name == "benchmark-scenario-restore"
    assert carrier.visible is False
    carrier.value = json.dumps(
        {"nonce": "request-1", "settings": {"contig": {"coloring": {"custom_color_rows": []}}}}
    )
    assert restored == [{"contig": {"coloring": {"custom_color_rows": []}}}]
    assert json.loads(status.value) == {"nonce": "request-1", "status": "completed"}


def test_filter_lookup_uses_live_monotonic_indices_and_explicit_occurrence():
    assert choose_filter_lookup_index([14, 3, 9, 3], 1) == 3
    assert choose_filter_lookup_index([14, 3, 9], 2) == 9

    with pytest.raises(ScenarioError, match="occurrence 4 does not exist"):
        choose_filter_lookup_index([3, 9], 4)
    with pytest.raises(ScenarioError, match="positive integer"):
        choose_filter_lookup_index([3], 0)
