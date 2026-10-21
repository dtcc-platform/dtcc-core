"""Execute public examples with real plotting and meshing, without windows.

Live workflows are opt-in: ``pytest tests/demos --run-demos``. The small
synthetic height example also runs offline in the ordinary suite.
"""

import runpy
from pathlib import Path
import numpy as np

import pytest

from dtcc_core.model.model import Model
from dtcc_core.model import City, Mesh, VolumeMesh, PointCloud, Raster
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
    # IO registers concrete save methods lazily; capture them before patching.
    import dtcc_core.io
    for cls in (City, Mesh, VolumeMesh, PointCloud, Raster):
        save = cls.save
        def save_in_temporary_directory(model, path, *args, _save=save, **kwargs):
            return _save(model, tmp_path / path, *args, **kwargs)
        monkeypatch.setattr(cls, 'save', save_in_temporary_directory)
    yield tmp_path
    plt.close("all")


@pytest.mark.demo
@pytest.mark.parametrize("demo_name", WORKFLOW_DEMOS + (
    "dataset_city_surface_mesh.py", "dataset_city_flat_mesh.py", "dataset_city.py",
))
def test_demo_executes_cleanly(demo_name, preview_without_windows):
    result = runpy.run_path(str(DEMO_DIR / demo_name), run_name="__main__")
    if demo_name in ('dataset_city.py', 'create_cityjson.py'):
        buildings = result['city'].buildings
        assert buildings and all(b.lod2 is not None or b.lod1 is not None for b in buildings)
        lod2 = [b.lod2 for b in buildings if b.lod2 is not None]
        assert lod2, 'The LOD2 demo must reconstruct roofs, not only fall back to LOD1'
        assert any(.01 < abs(s.calculate_normal()[2]) < .99 for g in lod2 for s in g.surfaces)
        if demo_name == 'create_cityjson.py':
            from dtcc_core.io import load_city
            loaded = load_city(preview_without_windows/'output/city.city.json')
            assert len(loaded.buildings) == len(buildings)
            assert sum(b.lod2 is not None for b in loaded.buildings) == len(lod2)
    elif demo_name == 'build_cityjson_surface_mesh.py':
        mesh = result['mesh']
        assert len(mesh.faces) and np.isfinite(mesh.vertices).all()
        xyz = mesh.vertices[mesh.faces]
        normals = np.cross(xyz[:, 1]-xyz[:, 0], xyz[:, 2]-xyz[:, 0])
        area2 = np.linalg.norm(normals, axis=1)
        assert np.all(area2 > 0)
        vertical_fraction = abs(normals[:, 2])/area2
        assert np.any((vertical_fraction > .01) & (vertical_fraction < .99))
        assert (preview_without_windows/'output/cityjson_lod2_surface_mesh.vtu').exists()


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
