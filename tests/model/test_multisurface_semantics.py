import numpy as np

from dtcc_core.model.geometry.surface import Surface, MultiSurface
from dtcc_core.model.enums import SurfaceSemantic


def test_multisurface_semantics_default_none():
    s = Surface(
        vertices=np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], dtype=float)
    )
    ms = MultiSurface(surfaces=[s])
    assert ms.semantics is None


def test_multisurface_semantics_assigned():
    s = Surface(
        vertices=np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], dtype=float)
    )
    ms = MultiSurface(surfaces=[s], semantics=[SurfaceSemantic.ROOF])
    assert ms.semantics is not None
    assert len(ms.semantics) == len(ms.surfaces)
    assert ms.semantics[0] == SurfaceSemantic.ROOF


def test_multisurface_semantics_none_backward_compatible():
    s = Surface(
        vertices=np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], dtype=float)
    )
    ms = MultiSurface(surfaces=[s])
    assert ms.semantics is None


def test_multisurface_copy_preserves_semantics():
    s = Surface(
        vertices=np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], dtype=float)
    )
    ms = MultiSurface(surfaces=[s], semantics=[SurfaceSemantic.WALL])
    ms2 = ms.copy(geometry_only=True)
    assert ms2.semantics is not None
    assert ms2.semantics[0] == SurfaceSemantic.WALL


def test_multisurface_copy_no_semantics():
    s = Surface(
        vertices=np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], dtype=float)
    )
    ms = MultiSurface(surfaces=[s])
    ms2 = ms.copy(geometry_only=True)
    assert ms2.semantics is None
