import re
from pathlib import Path

STATIC = Path("thebigbam/static")


def test_sidebar_design_tokens_and_semantic_classes_share_one_stylesheet():
    css = (STATIC / "bokeh_styles.css").read_text(encoding="utf-8")

    for contract in (
        "--sidebar-gutter",
        "--control-gap",
        ".sidebar-content",
        ".sidebar-busy",
        ".control-row",
        ".nested-section",
        ".action-row",
        ".sidebar-actions",
        ":host(.action-primary)",
        ":host(.action-add)",
        ":host(.action-muted)",
        ":host(.mode-switch)",
    ):
        assert contract in css

    assert ":root,\n:host {" in css
    assert re.search(r"var\([^)]*,", css) is None

    assert not (STATIC / "pink_buttons.css").exists()
    assert not (STATIC / "grey_buttons.css").exists()
