"""Native vertical scrolling for dynamically inserted read plots."""

from html import escape
from math import ceil


def scrollbar_html(range_id: str, height: int, row_count: int, rows_visible: float) -> str:
    """Bind on the scroll element itself, including plots added after document ready."""
    maximum = max(0.0, row_count + 1 - rows_visible)
    content_height = max(height, ceil(height * (row_count + 1) / rows_visible))
    handler = f"""
        const range = Bokeh.documents.map(doc => doc.get_model_by_id('{range_id}')).find(Boolean);
        if (range) {{
            const offset = this.scrollTop / Math.max(1, this.scrollHeight - this.clientHeight) * {maximum};
            range.end = {row_count + 1} - offset;
            range.start = range.end - {rows_visible};
        }}
    """
    return (
        '<style>.read-scrollbar::-webkit-scrollbar{width:14px}'
        '.read-scrollbar::-webkit-scrollbar-track{background:#eee}'
        '.read-scrollbar::-webkit-scrollbar-thumb{background:#888;border-radius:7px}</style>'
        f'<div class="read-scrollbar" style="height:{height}px;width:24px;box-sizing:border-box;'
        'overflow-y:scroll;overflow-x:hidden;scrollbar-gutter:stable;scrollbar-color:#888 #eee" '
        f'onscroll="{escape(handler, quote=True)}">'
        f'<div style="height:{content_height}px;width:1px"></div></div>'
    )
