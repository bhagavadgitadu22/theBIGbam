from dataclasses import dataclass

from ..apply_pipeline import MagOneApplyRequest, PlotResult
from .rendering import ApplyRenderEngine


@dataclass(frozen=True)
class MagOneHandler:
    engine: ApplyRenderEngine

    def render(self, request: MagOneApplyRequest, started_at: float) -> PlotResult:
        return self.engine.render_mag(request.inputs, started_at)
