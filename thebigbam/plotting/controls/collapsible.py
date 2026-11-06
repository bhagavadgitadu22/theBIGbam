"""Shared sidebar collapse buttons and content containers."""

import panel as pn
from bokeh.layouts import row
from bokeh.models import Button, Div


def build_collapsible(
    title,
    children,
    toggle_stylesheet="",
    make_toggle_callback=None,
    *,
    margin=(0, 0, 5, 0),
    slug=None,
    initially_open=False,
):
    title_text = title.text if isinstance(title, Div) else title
    benchmark_slug = slug or title_text.split(" (")[0].lower().replace(" ", "-")
    toggle = Button(label="▼" if initially_open else "▶", width=20, height=20, button_type="primary", align="center", margin=0,
                    stylesheets=[toggle_stylesheet], css_classes=[f"benchmark-toggle-{benchmark_slug}"])
    toggle.styles = {"padding": "0px", "line-height": "20px"}
    title_model = title if isinstance(title, Div) else Div(text=title, align="center")
    header = row(toggle, title_model, sizing_mode="stretch_width", align="center")
    content = pn.Column(*children, visible=initially_open, sizing_mode="stretch_width", margin=margin)

    def default_toggle():
        content.visible = not content.visible
        toggle.label = "▼" if content.visible else "▶"

    toggle.on_click(make_toggle_callback(toggle, content) if make_toggle_callback is not None else default_toggle)
    return header, content
