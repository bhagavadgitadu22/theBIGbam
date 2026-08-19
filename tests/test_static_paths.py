from pathlib import Path

from thebigbam.plotting.shared.paths import static_directory


def test_packaged_static_directory_contains_server_assets():
    assets = Path(static_directory())

    assert assets.name == "static"
    assert assets.parent.name == "thebigbam"
    assert (assets / "LOGO.png").is_file()
    assert (assets / "bokeh_styles.css").is_file()
