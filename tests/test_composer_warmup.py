import sys
import types

if "thebigbam_rs" not in sys.modules:
    rust_stub = types.ModuleType("thebigbam_rs")
    rust_stub.decode_dense_chunk = lambda data: []
    rust_stub.decode_sparse_chunk = lambda data: ([], [])
    sys.modules["thebigbam_rs"] = rust_stub

from thebigbam.plotting.application.apply_render_handlers import warm_plot_pipeline_imports


def test_pipeline_warmup_loads_all_apply_pipeline_modules():
    warm_plot_pipeline_imports()

    assert "thebigbam.plotting.application.sample_mag_pipeline" in sys.modules
    assert "thebigbam.plotting.application.all_samples_pipeline" in sys.modules
