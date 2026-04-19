import numpy as np
import pytest

from dtcc_core.model.geometry.surface import Surface, MultiSurface
from dtcc_core.model.enums import SurfaceSemantic
from dtcc_core.builder.geometry_builders.roof_geometry import (
    build_flat_geometry,
    build_fallback_geometry,
    build_gabled_geometry,
    build_hipped_geometry,
)
from dtcc_core.builder.geometry.shell_validation import validate_shell


def _make_lod1_box(ground_z=0.0, roof_z=5.0, size=10.0):
    bottom = Surface(vertices=np.array(
        [[0,0,ground_z],[size,0,ground_z],[size,size,ground_z],[0,size,ground_z]], dtype=float))
    top = Surface(vertices=np.array(
        [[0,0,roof_z],[size,0,roof_z],[size,size,roof_z],[0,size,roof_z]], dtype=float))
    w1 = Surface(vertices=np.array(
        [[0,0,ground_z],[size,0,ground_z],[size,0,roof_z],[0,0,roof_z]], dtype=float))
    w2 = Surface(vertices=np.array(
        [[size,0,ground_z],[size,size,ground_z],[size,size,roof_z],[size,0,roof_z]], dtype=float))
    w3 = Surface(vertices=np.array(
        [[size,size,ground_z],[0,size,ground_z],[0,size,roof_z],[size,size,roof_z]], dtype=float))
    w4 = Surface(vertices=np.array(
        [[0,size,ground_z],[0,0,ground_z],[0,0,roof_z],[0,size,roof_z]], dtype=float))
    return MultiSurface(surfaces=[bottom, top, w1, w2, w3, w4])


# --- Flat geometry tests ---


def test_build_flat_geometry_has_correct_semantics():
    lod1 = _make_lod1_box()
    result = build_flat_geometry(lod1)
    assert result.semantics is not None
    assert len(result.semantics) == len(result.surfaces)
    assert result.semantics[0] == SurfaceSemantic.GROUND
    assert result.semantics[-1] == SurfaceSemantic.ROOF
    for sem in result.semantics[1:-1]:
        assert sem == SurfaceSemantic.WALL


def test_build_fallback_geometry_same_as_flat():
    lod1 = _make_lod1_box()
    result = build_fallback_geometry(lod1)
    assert result.semantics is not None
    assert result.semantics[0] == SurfaceSemantic.GROUND
    assert result.semantics[-1] == SurfaceSemantic.ROOF


def test_flat_geometry_surface_count_matches_lod1():
    lod1 = _make_lod1_box()
    result = build_flat_geometry(lod1)
    assert len(result.surfaces) == len(lod1.surfaces)


# --- Gabled geometry tests ---


def test_build_gabled_geometry_has_roof_surfaces():
    lod1 = _make_lod1_box()
    ridge_line = np.array([[5, 0, 7], [5, 10, 7]])
    footprint = Surface(vertices=np.array(
        [[0,0,0],[10,0,0],[10,10,0],[0,10,0]], dtype=float))

    result = build_gabled_geometry(
        lod1=lod1, footprint=footprint,
        ridge_line=ridge_line, eave_height=5.0, ridge_height=7.0,
    )

    assert result.semantics is not None
    assert len(result.semantics) == len(result.surfaces)
    roof_count = sum(1 for s in result.semantics if s == SurfaceSemantic.ROOF)
    assert roof_count == 2
    ground_count = sum(1 for s in result.semantics if s == SurfaceSemantic.GROUND)
    assert ground_count == 1


def test_gabled_fallback_on_bad_ridge():
    lod1 = _make_lod1_box()
    # Ridge completely outside footprint -- split should fail
    ridge_line = np.array([[100, 100, 7], [200, 200, 7]])
    footprint = Surface(vertices=np.array(
        [[0,0,0],[10,0,0],[10,10,0],[0,10,0]], dtype=float))

    result = build_gabled_geometry(
        lod1=lod1, footprint=footprint,
        ridge_line=ridge_line, eave_height=5.0, ridge_height=7.0,
    )

    # Should fall back to flat
    assert result.semantics is not None
    assert result.semantics[-1] == SurfaceSemantic.ROOF


# --- Hipped geometry tests ---


def test_build_hipped_geometry_has_4_roof_surfaces():
    lod1 = _make_lod1_box()
    ridge_line = np.array([[3, 5, 7], [7, 5, 7]])
    footprint = Surface(vertices=np.array(
        [[0,0,0],[10,0,0],[10,10,0],[0,10,0]], dtype=float))

    result = build_hipped_geometry(
        lod1=lod1, footprint=footprint,
        ridge_line=ridge_line, eave_height=5.0, ridge_height=7.0,
    )

    assert result.semantics is not None
    assert len(result.semantics) == len(result.surfaces)
    roof_count = sum(1 for s in result.semantics if s == SurfaceSemantic.ROOF)
    assert roof_count == 4
    ground_count = sum(1 for s in result.semantics if s == SurfaceSemantic.GROUND)
    assert ground_count == 1
