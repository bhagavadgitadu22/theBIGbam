import pytest

from thebigbam import cli


def test_calculate_time_is_documented_as_final_developer_option(capsys):
    with pytest.raises(SystemExit) as exit_info:
        cli.main(["calculate", "--help"])

    assert exit_info.value.code == 0
    help_text = capsys.readouterr().out
    assert "--time" in help_text
    assert "(for developers) Print timing diagnostics" in help_text
    assert help_text.rindex("--time") > help_text.rindex("--blast")


def test_top_level_version_reports_installed_version(monkeypatch, capsys):
    monkeypatch.setattr(cli, "_get_version", lambda: "1.2.3")

    with pytest.raises(SystemExit) as exit_info:
        cli.main(["--version"])

    assert exit_info.value.code == 0
    assert capsys.readouterr().out == "thebigbam 1.2.3\n"


def test_unknown_calculate_option_is_rejected(capsys):
    with pytest.raises(SystemExit) as exit_info:
        cli.main(["calculate", "-o", "output.db", "--annotation_tool", "pharokka"])

    assert exit_info.value.code == 2
    captured = capsys.readouterr()
    assert captured.out == ""
    assert "unrecognized arguments: --annotation_tool pharokka" in captured.err
