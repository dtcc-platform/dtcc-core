import numpy as np

from dtcc_core.model.geometry.pointcloud import PointCloud
from dtcc_core.model.geometry.surface import Surface
from dtcc_core.model.object.building import Building
from dtcc_core.model.object.object import GeometryType
from dtcc_core.builder.geometry_builders.buildings import build_lod2_buildings
from dtcc_core.builder.geometry_builders.surface import extrude_surface


def _flat_roof_building(bid="b001"):
    b = Building(id=bid)
    fp = Surface(vertices=np.array([[0,0,0],[10,0,0],[10,10,0],[0,10,0]], dtype=float))
    b.add_geometry(fp, GeometryType.LOD0)
    b.attributes["ground_height"] = 0.0
    b.attributes["height"] = 5.0
    b.add_geometry(extrude_surface(fp, 5.0), GeometryType.LOD1)
    np.random.seed(0)
    pts = np.column_stack([
        np.random.rand(200) * 10,
        np.random.rand(200) * 10,
        np.full(200, 5.0) + np.random.randn(200) * 0.1,
    ])
    b.add_geometry(PointCloud(points=pts), GeometryType.POINT_CLOUD)
    return b


def test_diagnostics_records_stage_timings():
    b = _flat_roof_building()
    diagnostics = {}
    build_lod2_buildings([b], diagnostics_out=diagnostics)
    assert "b001" in diagnostics
    timings = diagnostics["b001"]["timings"]
    for stage in ("filter", "detect", "merge", "classify", "geometry", "validate"):
        assert stage in timings, f"missing stage: {stage}"
        assert timings[stage] >= 0.0


def test_diagnostics_records_planes_for_successful_build():
    b = _flat_roof_building()
    diagnostics = {}
    build_lod2_buildings([b], diagnostics_out=diagnostics)
    planes = diagnostics["b001"]["planes"]
    assert planes is not None
    assert len(planes) >= 1


def test_diagnostics_records_classification_result():
    b = _flat_roof_building()
    diagnostics = {}
    build_lod2_buildings([b], diagnostics_out=diagnostics)
    clsr = diagnostics["b001"]["classification"]
    assert clsr is not None
    assert clsr.roof_type.name == "FLAT"


def test_diagnostics_records_filtered_point_count():
    b = _flat_roof_building()
    diagnostics = {}
    build_lod2_buildings([b], diagnostics_out=diagnostics)
    assert diagnostics["b001"]["filtered_point_count"] > 0


def test_diagnostics_none_is_noop():
    b = _flat_roof_building()
    build_lod2_buildings([b])
    assert b.lod2 is not None
