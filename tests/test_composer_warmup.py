import sys
import types

if "thebigbam_rs" not in sys.modules:
    rust_stub = types.ModuleType("thebigbam_rs")
    rust_stub.decode_dense_chunk = lambda data: []
    rust_stub.decode_sparse_chunk = lambda data: ([], [])
    sys.modules["thebigbam_rs"] = rust_stub

from thebigbam.plotting.application.apply_render_handlers import warm_composer_imports


def test_composer_warmup_loads_all_apply_composer_modules():
    warm_composer_imports()

    assert "thebigbam.plotting.composers.sample_mag" in sys.modules
    assert "thebigbam.plotting.composers.all_samples" in sys.modules
