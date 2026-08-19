from dataclasses import dataclass

from ..apply_pipeline import ContigOneApplyRequest, PlotResult
from .rendering import ApplyRenderEngine


@dataclass(frozen=True)
class ContigOneHandler:
    engine: ApplyRenderEngine

    def render(self, request: ContigOneApplyRequest, started_at: float) -> PlotResult:
        return self.engine.render_contig(request.inputs, started_at)
