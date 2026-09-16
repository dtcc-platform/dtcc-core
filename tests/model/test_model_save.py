from types import ModuleType

import pytest

from dtcc_core.model import Point
import dtcc_core.model.model as model_module


def test_model_save_raises_when_io_unavailable(monkeypatch):
    def fake_import_module(name):
        assert name == "dtcc_core.io"
        raise ImportError("missing io")

    monkeypatch.setattr(model_module.importlib, "import_module", fake_import_module)

    with pytest.raises(AttributeError) as excinfo:
        Point().save("point.pb")

    assert "Point" in str(excinfo.value)
    assert "dtcc_core.io" in str(excinfo.value)


def test_model_save_delegates_after_io_import(monkeypatch):
    calls = []

    def fake_save(self, *args, **kwargs):
        calls.append((self, args, kwargs))
        return "saved"

    def fake_import_module(name):
        assert name == "dtcc_core.io"
        monkeypatch.setattr(Point, "save", fake_save, raising=False)
        return ModuleType("dtcc_core.io")

    monkeypatch.setattr(model_module.importlib, "import_module", fake_import_module)

    point = Point()
    result = point.save("point.pb", option=True)

    assert result == "saved"
    assert calls == [(point, ("point.pb",), {"option": True})]
