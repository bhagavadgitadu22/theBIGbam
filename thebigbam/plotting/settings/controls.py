"""Bokeh controls for saving plotting settings."""

from __future__ import annotations

from typing import Any, Callable, Mapping

import panel as pn
from bokeh.models import Div
from bokeh.models.callbacks import CustomJS

from ..shared.styles import right_panel_tooltip
from .persistence import save_settings_document


def save_confirmation_js(button_class: str) -> str:
    """Return the shared browser confirmation behavior for server-side saves."""
    script = """
if (!cb_obj.text) return;
var btn = document.querySelector('.__BUTTON_CLASS__');
if (btn) {
    var rect = btn.getBoundingClientRect();
    var tip = document.createElement('div');
    tip.textContent = 'Saved!';
    Object.assign(tip.style, {
        position: 'fixed', left: (rect.left + rect.width / 2) + 'px',
        top: (rect.top - 28) + 'px', transform: 'translateX(-50%)',
        background: '#333', color: '#fff', padding: '3px 8px',
        borderRadius: '4px', fontSize: '12px', pointerEvents: 'none',
        zIndex: '9999', opacity: '1', transition: 'opacity 0.4s ease'
    });
    document.body.appendChild(tip);
    setTimeout(function() { tip.style.opacity = '0'; }, 700);
    setTimeout(function() { if (tip.parentNode) tip.parentNode.removeChild(tip); }, 1100);
}
cb_obj.text = '';
"""
    return script.replace("__BUTTON_CLASS__", button_class)


class SettingsSaveControls:
    """Own the save button, confirmation carrier, and persistence callback."""

    def __init__(
        self,
        db_path: str,
        collect: Callable[[], Mapping[str, Any]],
        stylesheet: Any,
        record_action: Callable[[str, Mapping[str, Any]], Any] | None = None,
    ) -> None:
        self.db_path = db_path
        self.collect = collect
        self.confirmation = Div(text="", visible=False)
        self.confirmation.js_on_change("text", CustomJS(code=save_confirmation_js("save-settings-btn")))
        self.button = pn.widgets.Button(
            name="SAVE SETTINGS",
            description=right_panel_tooltip(
                "Save current settings only, including unapplied ones"
            ),
            align="center",
            button_type="primary",
            stylesheets=[stylesheet],
            css_classes=["action-primary", "apply-btn", "save-settings-btn"],
            margin=(5, 0, 0, 0),
        )
        self.button.on_click(self._save)
        if record_action is not None:
            self.button.on_click(lambda _event: record_action("save_settings", {}))

    def _save(self, event: Any = None) -> None:
        path = save_settings_document(self.collect(), self.db_path)
        print(f"[settings] Saved to {path}", flush=True)
        self.confirmation.text = "1"
