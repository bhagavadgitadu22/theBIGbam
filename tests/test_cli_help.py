import pytest

from thebigbam import cli


def test_calculate_diagnose_is_documented_as_final_option(capsys):
    with pytest.raises(SystemExit) as exit_info:
        cli.main(["calculate", "--help"])

    assert exit_info.value.code == 0
    help_text = capsys.readouterr().out
    assert "--diagnose" in help_text
    assert "Write calculation timing diagnostics" in help_text
    assert help_text.rindex("--diagnose") > help_text.rindex("--blast")


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


@pytest.mark.parametrize("mode", ["none", "termini", "mags", "all"])
def test_calculate_blast_modes(mode):
    args = cli.build_argparser().parse_args(["calculate", "-o", "output.db", "--blast", mode])
    assert args.blast == mode


def test_calculate_blast_defaults_to_none():
    args = cli.build_argparser().parse_args(["calculate", "-o", "output.db", "--extend"])
    assert args.blast == "none"


def test_blast_help_explains_large_mag_dataset_recommendation():
    help_text = cli.build_argparser()._subparsers._group_actions[0].choices["calculate"].format_help()
    assert "mags" in help_text
    assert "prefer mags over all" in " ".join(help_text.split())
    assert "not recommended for large datasets" in " ".join(help_text.split())


def test_blast_mags_requires_mag_view(capsys):
    from thebigbam.database.calculating_data import run_calculate_args
    args = cli.build_argparser().parse_args(["calculate", "-o", "output.db", "--blast", "mags"])
    with pytest.raises(SystemExit) as error:
        run_calculate_args(args)
    assert "--blast mags requires --view mag" in str(error.value)


def test_blast_all_help_includes_cross_mag_comparisons():
    help_text = cli.build_argparser()._subparsers._group_actions[0].choices["calculate"].format_help()
    assert "including different MAGs" in " ".join(help_text.split())


def test_blast_all_help_includes_intercontig_search_in_both_views():
    help_text = cli.build_argparser()._subparsers._group_actions[0].choices["calculate"].format_help()
    assert "all reference contigs in either view" in " ".join(help_text.split())
