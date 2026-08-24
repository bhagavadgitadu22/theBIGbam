"""Bokeh controls for saving plotting settings."""

from __future__ import annotations

from typing import Any, Callable, Mapping

from bokeh.models import Div
from bokeh.models.callbacks import CustomJS
from bokeh.models.widgets import Button

from .persistence import save_settings_document

_SAVE_CONFIRM_JS = """
if (!cb_obj.text) return;
var btn = document.querySelector('.save-settings-btn');
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


class SettingsSaveControls:
    """Own the save button, confirmation carrier, and persistence callback."""

    def __init__(
        self,
        db_path: str,
        collect: Callable[[], Mapping[str, Any]],
        stylesheet: Any,
    ) -> None:
        self.db_path = db_path
        self.collect = collect
        self.confirmation = Div(text="", visible=False)
        self.confirmation.js_on_change("text", CustomJS(code=_SAVE_CONFIRM_JS))
        self.button = Button(
            label="SAVE SETTINGS",
            align="center",
            button_type="primary",
            stylesheets=[stylesheet],
            css_classes=["action-primary", "apply-btn", "save-settings-btn"],
            margin=(5, 0, 0, 5),
        )
        self.button.on_click(self._save)

    def _save(self, event: Any = None) -> None:
        path = save_settings_document(self.collect(), self.db_path)
        print(f"[settings] Saved to {path}", flush=True)
        self.confirmation.text = "1"
