import numpy as np

from dtcc_core.model.object.building import Building
from dtcc_core.model.geometry.surface import Surface, MultiSurface
from dtcc_core.model.object.object import GeometryType
from dtcc_core.builder.evaluation.metrics import is_watertight


def test_is_watertight_none_lod2_returns_false():
    b = Building(id="b1")
    assert is_watertight(b) is False


def test_is_watertight_simple_box_returns_true():
    ms = MultiSurface()
    verts = np.array([
        [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
        [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1],
    ], dtype=float)
    faces = [
        [0, 1, 2, 3],
        [4, 5, 6, 7],
        [0, 1, 5, 4],
        [1, 2, 6, 5],
        [2, 3, 7, 6],
        [3, 0, 4, 7],
    ]
    for face in faces:
        ms.surfaces.append(Surface(vertices=verts[face]))
    b = Building(id="b2")
    b.add_geometry(ms, GeometryType.LOD2)
    assert is_watertight(b) is True
