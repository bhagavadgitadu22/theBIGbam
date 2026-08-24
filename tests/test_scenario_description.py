import json

from thebigbam import cli
from thebigbam.plotting.settings.scenario import ScenarioFormatError, describe_scenario_document


def _document(steps):
    return {"_meta": {"format": "thebigbam-scenario", "version": 1}, "steps": steps}


def test_description_uses_stable_sequence_and_sorted_settings_paths():
    document = _document(
        [
            {
                "sequence": 7,
                "action": "state_change",
                "changes": {
                    "selection": {"sample": "S2"},
                    "filtering": [{"rows": [{"category": "Coverage", "value": 100}]}],
                },
            },
            {
                "sequence": 9,
                "action": "filter_lookup",
                "details": {"column": "Coverage_mean", "category": "Coverage"},
                "status": "completed",
                "duration_seconds": 1.25,
            },
            {
                "sequence": 10,
                "action": "apply_filters",
                "status": "failed",
                "duration_seconds": 2,
            },
            {"sequence": 11, "action": "apply_plot"},
        ]
    )

    assert describe_scenario_document(document) == (
        '7. Change settings: filtering=[{"rows":[{"category":"Coverage","value":100}]}]; '
        'selection.sample="S2"',
        "9. Look up filter values for Coverage.Coverage_mean [completed, 1.25 s]",
        "10. Apply filters [failed, 2 s]",
        "11. Apply plot",
    )


def test_plain_scenario_description_needs_no_result_fields():
    assert describe_scenario_document(
        _document(
            [
                {"sequence": 1, "action": "filter_lookup", "details": {"category": "Coverage"}},
                {"sequence": 2, "action": "apply_filters"},
                {"sequence": 3, "action": "apply_plot"},
            ]
        )
    ) == (
        "1. Look up filter values for Coverage",
        "2. Apply filters",
        "3. Apply plot",
    )


def test_description_rejects_malformed_documents():
    try:
        describe_scenario_document(_document([{"sequence": 1}]))
    except ScenarioFormatError as error:
        assert str(error) == "step 1 has an invalid action"
    else:
        raise AssertionError("malformed scenario was accepted")


def test_cli_prints_lines_and_reports_invalid_json(tmp_path, capsys):
    valid = tmp_path / "valid.json"
    valid.write_text(json.dumps(_document([{"sequence": 4, "action": "apply_plot"}])), encoding="utf-8")

    assert cli.main(["describe-scenario", str(valid)]) == 0
    assert capsys.readouterr().out == "4. Apply plot\n"

    invalid = tmp_path / "invalid.json"
    invalid.write_text("{", encoding="utf-8")
    assert cli.main(["describe-scenario", str(invalid)]) == 2
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err.startswith("Error describing scenario: invalid JSON at line 1, column 2")


def test_describe_scenario_command_help_is_hidden():
    assert "describe-scenario" not in cli.build_argparser().format_help()
