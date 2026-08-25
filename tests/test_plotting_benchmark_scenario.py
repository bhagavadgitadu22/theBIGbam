import json

import pytest

from thebigbam import cli
from thebigbam.plotting.application.composition_root import (
    _build_scenario_action_carrier,
    _build_scenario_restore_carrier,
)
from thebigbam.plotting.benchmark.replay import (
    ScenarioError,
    StepResult,
    augment_scenario,
    install_safe_model_lookup,
    load_scenario,
    measure_memory,
    merge_changes,
    record_blocked_step,
    set_model,
    write_results,
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

        def evaluate(self, expression, arg=None):
            calls.append(("evaluate", expression, arg))

    set_model(Page(), "target", "value", "new")

    assert calls[0][0:1] == ("wait",)
    assert calls[2][0:1] == ("wait",)
    assert calls[2][2:] == ("target", None)
    assert calls[3][2] == ["target", "value", "new"]


def test_safe_model_lookup_skips_models_with_unset_names():
    scripts = []

    class Page:
        def wait_for_function(self, expression, *, arg=None, timeout=None):
            pass

        def evaluate(self, expression):
            scripts.append(expression)

    install_safe_model_lookup(Page())

    assert "try" in scripts[0]
    assert "model.name" in scripts[0]
    assert "continue" in scripts[0]


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
        StepResult(
            4,
            "state_change",
            "completed",
            0.125,
            memory={"server_rss_mb": 1024.5, "browser_heap_mb": 64.25},
        ),
        StepResult(7, "apply_plot", "failed", 4.218, error="Timeout"),
    ]

    result = augment_scenario(scenario, results, {"run": "cold"})

    assert scenario["steps"][0].get("status") is None
    assert result["initial_state"] == scenario["initial_state"]
    assert result["steps"][0]["duration_seconds"] == 0.125
    assert result["steps"][0]["memory"]["server_rss_mb"] == 1024.5
    assert result["steps"][1]["status"] == "failed"
    assert result["steps"][1]["error"] == "Timeout"
    assert result["_result"] == {"run": "cold"}


def test_memory_probe_is_best_effort(capsys):
    class Page:
        def evaluate(self, expression):
            raise RuntimeError("unsupported")

    assert measure_memory(Page(), None) == {"server_rss_mb": None, "browser_heap_mb": None}
    assert "Memory probe unavailable" in capsys.readouterr().out


def test_blocked_step_does_not_wait_for_or_click_control(tmp_path):
    class Page:
        def screenshot(self, **kwargs):
            pass

        def evaluate(self, expression):
            return None

    result = record_blocked_step(
        Page(),
        None,
        "cold",
        tmp_path,
        {"sequence": 4, "action": "apply_filters"},
        "filter restoration failed at step 3",
    )

    assert result.status == "blocked"
    assert result.duration_seconds == 0
    assert result.error == "filter restoration failed at step 3"


def test_results_tsv_flattens_memory_columns(tmp_path):
    result = StepResult(
        4,
        "state_change",
        "completed",
        0.125,
        memory={"server_rss_mb": 1024.5, "browser_heap_mb": None},
    )

    second_result = StepResult(
        7,
        "apply_plot",
        "completed",
        0.25,
        memory={"server_rss_mb": 1030, "browser_heap_mb": 65},
    )

    write_results(tmp_path, _scenario(), {"cold": [result, second_result]}, {})

    lines = (tmp_path / "results.tsv").read_text(encoding="utf-8").splitlines()
    assert "server_rss_mb\tbrowser_heap_mb" in lines[0]
    assert "1024.5\t" in lines[1]


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


def test_semantic_action_carrier_dispatches_before_acknowledgement():
    handled = []
    carrier, status = _build_scenario_action_carrier(
        lambda action, details: handled.append((action, details))
    )

    carrier.value = json.dumps(
        {
            "nonce": "action-1",
            "action": "filter_lookup",
            "details": {"category": "Termini", "column": "Packaging_mechanism", "opening": True},
        }
    )

    assert handled == [
        (
            "filter_lookup",
            {"category": "Termini", "column": "Packaging_mechanism", "opening": True},
        )
    ]
    assert status.value == ""
    curdoc = __import__("bokeh.io", fromlist=["curdoc"]).curdoc()
    list(curdoc.session_callbacks)[-1].callback()
    assert json.loads(status.value) == {"nonce": "action-1", "status": "completed"}
