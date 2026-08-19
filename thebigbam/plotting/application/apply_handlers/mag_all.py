from dataclasses import dataclass

from ..apply_pipeline import MagAllApplyRequest, PlotResult
from .rendering import ApplyRenderEngine


@dataclass(frozen=True)
class MagAllHandler:
    engine: ApplyRenderEngine

    def render(self, request: MagAllApplyRequest, started_at: float) -> PlotResult:
        return self.engine.render_mag(request.inputs, started_at)
