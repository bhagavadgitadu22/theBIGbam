import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MARKDOWN_LINK = re.compile(r"\[[^]]*\]\(([^)]+)\)")


def test_local_markdown_links_resolve():
    markdown_files = [ROOT / "README.md", *sorted((ROOT / "docs").glob("*.md"))]
    broken_links = []

    for markdown_file in markdown_files:
        text = markdown_file.read_text(encoding="utf-8")
        for match in MARKDOWN_LINK.finditer(text):
            target = match.group(1).split("#", 1)[0]
            if not target or "://" in target or target.startswith(("mailto:", "<")):
                continue

            destination = markdown_file.parent / target
            if not destination.exists():
                line = text.count("\n", 0, match.start()) + 1
                broken_links.append(f"{markdown_file.relative_to(ROOT)}:{line}: {target}")

    assert not broken_links, "Broken local Markdown links:\n" + "\n".join(broken_links)
