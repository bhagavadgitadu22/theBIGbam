from dataclasses import dataclass

from ..apply_pipeline import ContigAllApplyRequest, PlotResult
from .rendering import ApplyRenderEngine


@dataclass(frozen=True)
class ContigAllHandler:
    engine: ApplyRenderEngine

    def render(self, request: ContigAllApplyRequest, started_at: float) -> PlotResult:
        return self.engine.render_contig(request.inputs, started_at)
