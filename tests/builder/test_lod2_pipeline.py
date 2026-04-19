import numpy as np
import pytest

from dtcc_core.model.geometry.pointcloud import PointCloud
from dtcc_core.model.geometry.surface import Surface, MultiSurface
from dtcc_core.model.object.building import Building
from dtcc_core.model.object.object import GeometryType
from dtcc_core.model.enums import SurfaceSemantic
from dtcc_core.builder.geometry_builders.buildings import build_lod2_buildings
from dtcc_core.builder.geometry_builders.surface import extrude_surface


def _make_building_with_flat_roof():
    b = Building()
    fp = Surface(
        vertices=np.array([[0,0,0],[10,0,0],[10,10,0],[0,10,0]], dtype=float)
    )
    b.add_geometry(fp, GeometryType.LOD0)
    b.attributes["ground_height"] = 0.0
    b.attributes["height"] = 5.0

    # Build LoD1
    lod1 = extrude_surface(fp, 5.0)
    b.add_geometry(lod1, GeometryType.LOD1)

    # Flat roof point cloud
    pts = np.column_stack([
        np.random.rand(100) * 10,
        np.random.rand(100) * 10,
        np.full(100, 5.0) + np.random.randn(100) * 0.1,
    ])
    b.add_geometry(PointCloud(points=pts), GeometryType.POINT_CLOUD)
    return b


def test_build_lod2_flat_building():
    b = _make_building_with_flat_roof()
    build_lod2_buildings([b])

    lod2 = b.lod2
    assert lod2 is not None
    assert lod2.semantics is not None
    assert SurfaceSemantic.ROOF in lod2.semantics
    assert SurfaceSemantic.GROUND in lod2.semantics
    assert b.attributes.get("roof_type") is not None


def test_build_lod2_no_pointcloud_falls_back():
    b = Building()
    fp = Surface(
        vertices=np.array([[0,0,5],[10,0,5],[10,10,5],[0,10,5]], dtype=float)
    )
    b.add_geometry(fp, GeometryType.LOD0)
    lod1 = extrude_surface(fp, 0)
    b.add_geometry(lod1, GeometryType.LOD1)

    build_lod2_buildings([b])

    assert b.attributes.get("fallback_reason") == "insufficient_points"


def test_build_lod2_no_pointcloud_no_fallback():
    b = Building()
    fp = Surface(
        vertices=np.array([[0,0,5],[10,0,5],[10,10,5],[0,10,5]], dtype=float)
    )
    b.add_geometry(fp, GeometryType.LOD0)
    lod1 = extrude_surface(fp, 0)
    b.add_geometry(lod1, GeometryType.LOD1)

    build_lod2_buildings([b], fallback_to_flat=False)

    assert b.lod2 is None
    assert b.attributes.get("fallback_reason") == "insufficient_points"


def test_build_lod2_attributes_are_json_serializable():
    import json
    b = _make_building_with_flat_roof()
    build_lod2_buildings([b])

    # All LoD2 attributes should be JSON-serializable
    attrs = {
        k: v for k, v in b.attributes.items()
        if k in ("roof_type", "roof_confidence", "fallback_reason")
    }
    serialized = json.dumps(attrs)
    assert isinstance(serialized, str)


def test_cityjson_semantic_mapping():
    from dtcc_core.io.cityjson.converters import semantic_type_for_surface
    assert semantic_type_for_surface(SurfaceSemantic.GROUND) == "GroundSurface"
    assert semantic_type_for_surface(SurfaceSemantic.WALL) == "WallSurface"
    assert semantic_type_for_surface(SurfaceSemantic.ROOF) == "RoofSurface"
    assert semantic_type_for_surface(None) == "WallSurface"


def test_build_lod2_rebuild_false_skips_existing():
    b = _make_building_with_flat_roof()
    build_lod2_buildings([b])
    assert b.lod2 is not None

    # Mark the existing LoD2 so we can tell if it was replaced
    original_lod2 = b.lod2
    build_lod2_buildings([b], rebuild=False)
    assert b.lod2 is original_lod2  # should not have been replaced
