import json

import pytest

from thebigbam import cli
from thebigbam.plotting.application.composition_root import _build_scenario_restore_carrier
from thebigbam.plotting.benchmark.replay import (
    ScenarioError,
    StepResult,
    augment_scenario,
    choose_filter_lookup_index,
    find_filter_lookup_index,
    load_scenario,
    merge_changes,
    set_model,
)


def test_replay_scenario_cli_delegates_to_packaged_runner(monkeypatch):
    captured = []
    monkeypatch.setattr(cli, "run_replay", lambda args: captured.append(args) or 0)

    result = cli.main(["replay-scenario", "--scenario", "workflow.json", "--db", "example.db"])

    assert result == 0
    assert captured[0].scenario.name == "workflow.json"
    assert captured[0].db.name == "example.db"


def test_set_model_uses_keyword_playwright_argument():
    calls = []

    class Page:
        def wait_for_function(self, expression, *, arg=None, timeout=None):
            calls.append(("wait", expression, arg, timeout))

        def evaluate(self, expression, arg):
            calls.append(("evaluate", expression, arg))

    set_model(Page(), "target", "value", "new")

    assert calls[0][0:1] == ("wait",)
    assert calls[0][2:] == ("target", None)
    assert calls[1][2] == ["target", "value", "new"]


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

    carrier, status = _build_scenario_restore_carrier(Session(), lambda: {"filters_pending": True})

    assert carrier.name == "benchmark-scenario-restore"
    assert carrier.visible is False
    carrier.value = json.dumps(
        {"nonce": "request-1", "settings": {"contig": {"coloring": {"custom_color_rows": []}}}}
    )
    assert restored == [{"contig": {"coloring": {"custom_color_rows": []}}}]
    curdoc = __import__("bokeh.io", fromlist=["curdoc"]).curdoc()
    callbacks = list(curdoc.session_callbacks)
    assert status.value == ""
    callbacks[-1].callback()
    assert json.loads(status.value) == {
        "nonce": "request-1",
        "status": "completed",
        "filters_pending": True,
    }


def test_filter_lookup_uses_live_monotonic_indices_and_explicit_occurrence():
    assert choose_filter_lookup_index([14, 3, 9, 3], 1) == 3
    assert choose_filter_lookup_index([14, 3, 9], 2) == 9

    with pytest.raises(ScenarioError, match="occurrence 4 does not exist"):
        choose_filter_lookup_index([3, 9], 4)
    with pytest.raises(ScenarioError, match="positive integer"):
        choose_filter_lookup_index([3], 0)


def test_filter_lookup_ignores_semantically_matching_stale_rows():
    calls = []

    class Locator:
        def evaluate_all(self, expression):
            calls.append(("locator", expression))
            return [19]

    class Page:
        def locator(self, selector):
            calls.append(("selector", selector))
            return Locator()

        def evaluate(self, expression, arg):
            calls.append(("evaluate", expression, arg))
            # Model 7 is stale; only rendered row 19 may be returned.
            assert arg["liveIndices"] == [19]
            return [19]

    index = find_filter_lookup_index(Page(), {"category": "MAG coverage", "column": "Coverage_mean"})

    assert index == 19
