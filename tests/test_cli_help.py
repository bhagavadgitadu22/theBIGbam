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
