"""Verify existing LoD1 pipeline is unaffected by LoD2 changes."""
import numpy as np

from dtcc_core.model.geometry.surface import Surface, MultiSurface
from dtcc_core.builder.geometry_builders.surface import extrude_surface


def test_lod1_extrusion_still_works():
    fp = Surface(
        vertices=np.array([[0,0,5],[10,0,5],[10,10,5],[0,10,5]], dtype=float)
    )
    result = extrude_surface(fp, 0)
    assert isinstance(result, MultiSurface)
    assert len(result.surfaces) >= 6
    assert result.semantics is None  # LoD1 should NOT have semantics


def test_multisurface_without_semantics_backward_compatible():
    s = Surface(vertices=np.array([[0,0,0],[1,0,0],[1,1,0]], dtype=float))
    ms = MultiSurface(surfaces=[s])
    assert ms.semantics is None
    ms2 = ms.copy()
    assert ms2.semantics is None
