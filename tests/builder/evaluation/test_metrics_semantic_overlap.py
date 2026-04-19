import numpy as np
import pytest

from dtcc_core.model.geometry.surface import MultiSurface, Surface
from dtcc_core.model.enums import SurfaceSemantic
from dtcc_core.builder.evaluation.metrics import semantic_overlap


def _unit_square_at_z(z: float) -> Surface:
    return Surface(vertices=np.array([[0,0,z],[1,0,z],[1,1,z],[0,1,z]], dtype=float))


def test_identical_shapes_overlap_is_one():
    pred = MultiSurface()
    pred.surfaces.append(_unit_square_at_z(0))
    pred.semantics = [SurfaceSemantic.ROOF]

    truth = MultiSurface()
    truth.surfaces.append(_unit_square_at_z(0))
    truth.semantics = [SurfaceSemantic.ROOF]

    overlaps = semantic_overlap(pred, truth)
    assert overlaps[SurfaceSemantic.ROOF] == pytest.approx(1.0, abs=1e-9)


def test_no_overlap_is_zero():
    pred = MultiSurface()
    pred.surfaces.append(Surface(vertices=np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0]], dtype=float)))
    pred.semantics = [SurfaceSemantic.ROOF]

    truth = MultiSurface()
    truth.surfaces.append(Surface(vertices=np.array([[10,10,0],[11,10,0],[11,11,0],[10,11,0]], dtype=float)))
    truth.semantics = [SurfaceSemantic.ROOF]

    overlaps = semantic_overlap(pred, truth)
    assert overlaps[SurfaceSemantic.ROOF] == 0.0


def test_half_overlap():
    pred = MultiSurface()
    pred.surfaces.append(Surface(vertices=np.array([[0,0,0],[2,0,0],[2,1,0],[0,1,0]], dtype=float)))
    pred.semantics = [SurfaceSemantic.ROOF]

    truth = MultiSurface()
    truth.surfaces.append(Surface(vertices=np.array([[1,0,0],[3,0,0],[3,1,0],[1,1,0]], dtype=float)))
    truth.semantics = [SurfaceSemantic.ROOF]

    overlaps = semantic_overlap(pred, truth)
    assert overlaps[SurfaceSemantic.ROOF] == pytest.approx(1.0 / 3.0, abs=1e-9)


def test_missing_semantic_in_pred_is_zero():
    pred = MultiSurface()
    pred.surfaces.append(_unit_square_at_z(0))
    pred.semantics = [SurfaceSemantic.WALL]

    truth = MultiSurface()
    truth.surfaces.append(_unit_square_at_z(0))
    truth.semantics = [SurfaceSemantic.ROOF]

    overlaps = semantic_overlap(pred, truth)
    assert overlaps[SurfaceSemantic.ROOF] == 0.0
