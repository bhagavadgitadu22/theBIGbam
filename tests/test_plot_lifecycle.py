from types import SimpleNamespace

from thebigbam.plotting.shared.lifecycle import PlotLifecycle


class Model:
    def __init__(self):
        self.removed = []

    def remove_on_change(self, attribute, callback):
        self.removed.append((attribute, callback))


def test_prepare_replacement_detaches_callbacks_and_releases_range():
    lifecycle = PlotLifecycle()
    model = Model()
    callback = object()
    lifecycle.state["range_callbacks"] = ((model, "start", callback),)
    lifecycle.state["shared_xrange"] = object()

    lifecycle.prepare_replacement()

    assert model.removed == [("start", callback)]
    assert lifecycle.state["range_callbacks"] == ()
    assert lifecycle.state["shared_xrange"] is None


def test_shared_xrange_finds_first_figure():
    expected = object()
    grid = SimpleNamespace(children=[(SimpleNamespace(x_range=expected), 0, 0)])
    assert PlotLifecycle.shared_xrange(grid) is expected
