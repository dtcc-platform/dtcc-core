from pathlib import Path

import numpy as np
import pytest

from dtcc_core.model.object.building import Building
from dtcc_core.model.object.object import GeometryType
from dtcc_core.model.geometry.surface import MultiSurface, Surface
from dtcc_core.builder.evaluation.vtk_export import write_building_vtk


def _square(z: float = 0.0) -> Surface:
    return Surface(vertices=np.array([[0,0,z],[1,0,z],[1,1,z],[0,1,z]], dtype=float))


def test_write_building_vtk_creates_file(tmp_path: Path):
    b = Building(id="b1")
    ms = MultiSurface()
    ms.surfaces.append(_square())
    b.add_geometry(ms, GeometryType.LOD2)
    out = tmp_path / "b1.vtk"
    write_building_vtk(b, out)
    assert out.exists()
    assert out.stat().st_size > 0


def test_write_building_vtk_no_lod2_raises(tmp_path: Path):
    b = Building(id="b2")
    out = tmp_path / "b2.vtk"
    with pytest.raises(ValueError):
        write_building_vtk(b, out)
