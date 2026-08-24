from types import SimpleNamespace

from thebigbam.plotting.application.interactions import InteractionCoordinator


def test_interaction_coordinator_restores_each_control_state():
    enabled = SimpleNamespace(disabled=False)
    already_disabled = SimpleNamespace(disabled=True)
    root = SimpleNamespace(
        objects=[enabled, SimpleNamespace(objects=[already_disabled])], css_classes=["left-col"]
    )
    coordinator = InteractionCoordinator()
    coordinator.attach(root)

    assert coordinator.begin("plot")
    assert enabled.disabled
    assert already_disabled.disabled
    assert not coordinator.begin("controls")
    assert root.css_classes == ["left-col", "sidebar-busy"]

    coordinator.end()

    assert not enabled.disabled
    assert already_disabled.disabled
    assert not coordinator.busy
    assert root.css_classes == ["left-col"]
