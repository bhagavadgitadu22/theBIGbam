"""Client-side drag surface for resizing a serve sidebar."""

import param
from panel.custom import JSComponent


class PanelResizer(JSComponent):
    """Project pointer drag deltas into a bounded Panel width parameter."""

    side = param.Selector(default="left", objects=["left", "right"])
    storage_key = param.String()
    minimum = param.Integer(default=180, bounds=(1, None))
    maximum = param.Integer(default=700, bounds=(1, None))
    default_width = param.Integer(default=300, bounds=(1, None))
    value = param.Integer(default=300, bounds=(1, None))
    enabled = param.Boolean(default=True)

    _stylesheets = [
        """
        :host {
          cursor: col-resize;
          height: 100%;
          margin: 0;
          min-height: 40px;
          width: 8px;
        }
        .panel-resizer {
          background: transparent;
          height: 100%;
          touch-action: none;
          width: 8px;
        }
        """
    ]

    _esm = """
    export function render({ model }) {
      const handle = document.createElement("div");
      handle.className = `panel-resizer panel-resizer-${model.side}`;
      handle.title = "Drag to resize; double-click to reset";
      let startX = 0;
      let startWidth = model.value;

      function clamp(value) {
        return Math.max(model.minimum, Math.min(model.maximum, Math.round(value)));
      }

      function applyWidth(value, persist=true) {
        const width = clamp(value);
        if (model.value !== width) model.value = width;
        if (persist) localStorage.setItem(model.storage_key, String(width));
      }

      const saved = Number.parseInt(localStorage.getItem(model.storage_key), 10);
      applyWidth(Number.isFinite(saved) ? saved : model.default_width, false);

      handle.addEventListener("pointerdown", (event) => {
        if (!model.enabled) return;
        event.preventDefault();
        startX = event.clientX;
        startWidth = model.value;
        handle.setPointerCapture(event.pointerId);
        handle.classList.add("dragging");
      });
      handle.addEventListener("pointermove", (event) => {
        if (!handle.hasPointerCapture(event.pointerId)) return;
        const delta = event.clientX - startX;
        applyWidth(startWidth + (model.side === "left" ? delta : -delta));
      });
      handle.addEventListener("pointerup", (event) => {
        if (handle.hasPointerCapture(event.pointerId)) handle.releasePointerCapture(event.pointerId);
        handle.classList.remove("dragging");
      });
      handle.addEventListener("dblclick", () => {
        if (model.enabled) applyWidth(model.default_width);
      });
      return handle;
    }
    """
