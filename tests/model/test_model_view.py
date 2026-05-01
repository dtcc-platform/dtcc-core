from types import ModuleType

from dtcc_core.model import Point
import dtcc_core.model.model as model_module


def test_model_view_warns_when_viewer_unavailable(monkeypatch):
    messages = []

    def fake_import_module(name):
        assert name == "dtcc_viewer"
        raise ImportError("missing viewer")

    monkeypatch.setattr(model_module.importlib, "import_module", fake_import_module)
    monkeypatch.setattr(model_module, "warning", messages.append)

    result = Point().view()

    assert result is None
    assert len(messages) == 1
    assert "dtcc-viewer" in messages[0]
    assert "Point" in messages[0]


def test_model_view_delegates_after_viewer_import(monkeypatch):
    calls = []

    def fake_view(self, *args, **kwargs):
        calls.append((self, args, kwargs))
        return "viewed"

    def fake_import_module(name):
        assert name == "dtcc_viewer"
        monkeypatch.setattr(Point, "view", fake_view, raising=False)
        return ModuleType("dtcc_viewer")

    monkeypatch.setattr(model_module.importlib, "import_module", fake_import_module)

    point = Point()
    result = point.view("arg", option=True)

    assert result == "viewed"
    assert calls == [(point, ("arg",), {"option": True})]
