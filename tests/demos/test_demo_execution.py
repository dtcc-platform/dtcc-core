"""Execute public examples with real plotting and meshing, without windows.

Live workflows are opt-in: ``pytest tests/demos --run-demos``. The small
synthetic height example also runs offline in the ordinary suite.
"""

import runpy
from pathlib import Path

import pytest

from dtcc_core.model.model import Model
from tests.demos.demo_inventory import WORKFLOW_DEMOS

REPO_ROOT = Path(__file__).resolve().parents[2]
DEMO_DIR = REPO_ROOT / "demos"


@pytest.fixture
def preview_without_windows(monkeypatch, tmp_path):
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    monkeypatch.setattr(Model, "view", lambda *args, **kwargs: None)
    monkeypatch.setattr(plt, "show", lambda *args, **kwargs: None)
    monkeypatch.chdir(REPO_ROOT)
    save = Model.save

    def save_in_temporary_directory(model, path, *args, **kwargs):
        return save(model, tmp_path / path, *args, **kwargs)

    monkeypatch.setattr(Model, "save", save_in_temporary_directory)
    yield
    plt.close("all")


@pytest.mark.demo
@pytest.mark.parametrize("demo_name", WORKFLOW_DEMOS + (
    "dataset_city_surface_mesh.py", "dataset_city_flat_mesh.py",
))
def test_demo_executes_cleanly(demo_name, preview_without_windows):
    runpy.run_path(str(DEMO_DIR / demo_name), run_name="__main__")


def test_height_demo_keeps_three_roofs_in_one_block(preview_without_windows):
    result = runpy.run_path(str(DEMO_DIR / "build_lod1_buildings.py"))
    assert len(result["city"].buildings) == 3
    assert len(result["cleaned"].buildings) == 1
    parts = result["cleaned"].buildings[0].building_parts
    assert sorted(part.estimated_height for part in parts) == [8, 12, 16]
    assert all(part.lod1 is not None for part in parts)
    assert sorted(part.lod1.zmax for part in parts) == [8, 12, 16]
    shell = result["shell"]
    z = shell.vertices[shell.faces][:, :, 2]
    horizontal = z.max(axis=1) == z.min(axis=1)
    assert sorted(set(z[horizontal, 0])) == [0, 8, 12, 16]
