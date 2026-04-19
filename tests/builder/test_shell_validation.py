import numpy as np

from dtcc_core.model.geometry.surface import Surface, MultiSurface
from dtcc_core.model.enums import SurfaceSemantic
from dtcc_core.builder.geometry.shell_validation import validate_shell


def _make_box_multisurface():
    bottom = Surface(vertices=np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0]], dtype=float))
    top = Surface(vertices=np.array([[0,0,1],[1,0,1],[1,1,1],[0,1,1]], dtype=float))
    front = Surface(vertices=np.array([[0,0,0],[1,0,0],[1,0,1],[0,0,1]], dtype=float))
    back = Surface(vertices=np.array([[0,1,0],[1,1,0],[1,1,1],[0,1,1]], dtype=float))
    left = Surface(vertices=np.array([[0,0,0],[0,1,0],[0,1,1],[0,0,1]], dtype=float))
    right = Surface(vertices=np.array([[1,0,0],[1,1,0],[1,1,1],[1,0,1]], dtype=float))
    semantics = [
        SurfaceSemantic.GROUND, SurfaceSemantic.ROOF,
        SurfaceSemantic.WALL, SurfaceSemantic.WALL,
        SurfaceSemantic.WALL, SurfaceSemantic.WALL,
    ]
    return MultiSurface(surfaces=[bottom, top, front, back, left, right], semantics=semantics)


def test_valid_closed_shell():
    ms = _make_box_multisurface()
    is_valid, issues = validate_shell(ms, tolerance=0.01)
    assert is_valid
    assert len(issues) == 0


def test_open_shell_detected():
    ms = _make_box_multisurface()
    ms.surfaces = ms.surfaces[:5]
    ms.semantics = ms.semantics[:5]
    is_valid, issues = validate_shell(ms, tolerance=0.01)
    assert not is_valid
    assert len(issues) > 0
