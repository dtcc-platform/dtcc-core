import subprocess
import sys
import textwrap

import numpy as np
import pytest
from shapely.geometry import Point, Polygon, box

from dtcc_core.builder import (
    _dtcc_builder,
    build_city_flat_mesh,
    build_city_surface_mesh,
    build_city_volume_mesh,
    fix_building_footprint_clearance,
    merge_building_footprints,
    simplify_building_footprints,
)
from dtcc_core.builder.building.modify import clean_building_footprints
from dtcc_core.builder.geometry_builders import meshes as meshes_module
from dtcc_core.builder.model_conversion import builder_mesh_to_mesh, create_builder_polygon
from dtcc_core.builder.meshing import dtcc_mesher_backend as dtcc_mesher_backend_module
from dtcc_core.builder.meshing import flat_mesh_backends as flat_mesh_backends_module
from dtcc_core.builder.meshing.tetgen import is_tetgen_available
from dtcc_core.model import Building, Bounds, City, GeometryType, Mesh, Raster, Surface, Terrain, VolumeMesh


def make_surface(polygon: Polygon, z: float) -> Surface:
    surface = Surface()
    surface.from_polygon(polygon, z)
    return surface


def make_building(
    polygon: Polygon,
    *,
    roof_z: float = 10.0,
    lod1_polygon: Polygon | None = None,
    lod1_z: float | None = None,
    lod2_polygon: Polygon | None = None,
    lod2_z: float | None = None,
) -> Building:
    building = Building()
    building.add_geometry(make_surface(polygon, roof_z), GeometryType.LOD0)
    if lod1_polygon is not None:
        building.add_geometry(
            make_surface(lod1_polygon, lod1_z if lod1_z is not None else roof_z),
            GeometryType.LOD1,
        )
    if lod2_polygon is not None:
        building.add_geometry(
            make_surface(lod2_polygon, lod2_z if lod2_z is not None else roof_z),
            GeometryType.LOD2,
        )
    building.attributes["height"] = roof_z
    building.attributes["ground_height"] = 0.0
    return building


def make_flat_city(buildings: list[Building]) -> City:
    city = City()
    raster = Raster()
    raster.data = np.zeros((8, 8), dtype=float)
    raster.set_bounds(Bounds(0.0, 0.0, 80.0, 80.0))
    terrain = Terrain()
    terrain.add_geometry(raster, GeometryType.RASTER)
    city.add_terrain(terrain)
    city.add_buildings(buildings)
    return city


def _nearest_nonadjacent_boundary_vertex_distance(polygon: Polygon) -> float:
    exterior = np.asarray(polygon.exterior.coords, dtype=np.float64)
    best = np.inf

    for i in range(len(exterior) - 1):
        for j in range(i + 2, len(exterior) - 1):
            if i == 0 and j == len(exterior) - 2:
                continue
            best = min(best, float(np.linalg.norm(exterior[i] - exterior[j])))

    return best


def _minimum_boundary_edge_length(polygon: Polygon) -> float:
    best = np.inf
    for ring in [polygon.exterior, *polygon.interiors]:
        coords = np.asarray(ring.coords, dtype=np.float64)
        if len(coords) < 2:
            continue
        lengths = np.hypot(
            np.diff(coords[:, 0]),
            np.diff(coords[:, 1]),
        )
        if lengths.size:
            best = min(best, float(lengths.min()))
    return best


def test_build_city_flat_mesh_handles_pathological_footprints(monkeypatch):
    calls = []
    original = meshes_module._condition_meshing_footprints

    def wrapped(*args, **kwargs):
        calls.append(kwargs.get("lod"))
        return original(*args, **kwargs)

    monkeypatch.setattr(meshes_module, "_condition_meshing_footprints", wrapped)

    buildings = [
        make_building(Polygon([(8, 8), (18, 8), (18, 18), (8, 18), (8, 8)]), roof_z=12.0),
        make_building(Polygon([(18.1, 8), (28.1, 8), (28.1, 18), (18.1, 18), (18.1, 8)]), roof_z=8.0),
        make_building(Polygon([(35, 8), (42, 15), (35, 15), (42, 8), (35, 8)]), roof_z=6.0),
    ]
    city = make_flat_city(buildings)

    mesh = build_city_flat_mesh(
        city,
        lod=GeometryType.LOD0,
        merge_buildings=True,
        min_building_detail=0.3,
        min_building_area=1.0,
        merge_tolerance=0.25,
        max_mesh_size=5.0,
        min_mesh_angle=20.0,
        report_mesh_quality=False,
        mesher="spade",
    )

    assert calls
    assert mesh.vertices.shape[0] > 0
    assert mesh.faces.shape[0] > 0


def test_build_city_flat_mesh_forwards_cleaning_diagnostics_flag(monkeypatch):
    captured = {}
    original = meshes_module._condition_meshing_footprints

    def wrapped(*args, **kwargs):
        captured["cleaning_diagnostics"] = kwargs.get("cleaning_diagnostics")
        return original(*args, **kwargs)

    monkeypatch.setattr(meshes_module, "_condition_meshing_footprints", wrapped)

    city = make_flat_city(
        [make_building(box(8, 8, 18, 18), roof_z=10.0)]
    )

    mesh = build_city_flat_mesh(
        city,
        lod=GeometryType.LOD0,
        merge_buildings=True,
        min_building_detail=0.3,
        min_building_area=1.0,
        merge_tolerance=0.25,
        max_mesh_size=5.0,
        min_mesh_angle=20.0,
        report_mesh_quality=False,
        cleaning_diagnostics=False,
        mesher="spade",
    )

    assert captured["cleaning_diagnostics"] is False
    assert mesh.vertices.shape[0] > 0
    assert mesh.faces.shape[0] > 0


def test_condition_meshing_footprints_regularizes_touching_holes():
    touching_holes = Polygon(
        [(8, 8), (28, 8), (28, 28), (8, 28), (8, 8)],
        [
            [(12, 12), (16, 12), (16, 16), (12, 16), (12, 12)],
            [(16, 16), (22, 16), (22, 24), (16, 24), (16, 16)],
        ],
    )
    assert touching_holes.is_valid

    surfaces, source_map, resolutions, diagnostics = meshes_module._condition_meshing_footprints(
        [make_building(touching_holes, roof_z=10.0)],
        lod=GeometryType.LOD0,
        min_building_detail=0.5,
        min_building_area=1.0,
        merge_tolerance=0.5,
        merge_buildings=False,
        max_mesh_size=5.0,
        cleaning_diagnostics=False,
    )

    assert len(surfaces) == 1
    assert source_map == [[0]]
    assert resolutions == [5.0]
    normalized = surfaces[0].to_polygon(simplify=0.0)
    assert not meshes_module._polygon_has_ring_boundary_contacts(normalized)
    assert diagnostics.get("mesher_regularized_polygon_count", 0) == 0


def test_regularize_flat_mesh_ground_polygons_splits_case55_style_pinch():
    polygon = Polygon(
        [
            (675263.1249638698, 6581229.500071362),
            (675274.8125201109, 6581225.406301662),
            (675282.3125684819, 6581222.187530903),
            (675278.5937848733, 6581210.874928664),
            (675263.4062842835, 6581216.312428876),
            (675261.562570048, 6581210.906283745),
            (675276.0313192445, 6581205.718784033),
            (675274.656285481, 6581201.468679672),
            (675268.9687844442, 6581203.37493002),
            (675268.375039339, 6581201.65619419),
            (675265.3125405194, 6581201.687444177),
            (675263.0312874243, 6581194.437434341),
            (675245.0000353431, 6581198.656184828),
            (675241.6875341476, 6581189.499931524),
            (675232.7499284449, 6581192.249964047),
            (675235.6562160432, 6581200.437569969),
            (675239.9999676079, 6581198.96881944),
            (675241.9374480231, 6581203.937519215),
            (675248.093715854, 6581221.531320173),
            (675263.4062151447, 6581216.3125704145),
            (675265.8124304013, 6581223.499966635),
            (675261.7499317214, 6581224.999966148),
        ]
    )

    direct = meshes_module._regularize_flat_mesh_ground_polygons(
        [polygon],
        cleanup_scale=0.05,
        cleaning_diagnostics=False,
    )

    assert _nearest_nonadjacent_boundary_vertex_distance(polygon) < 1e-3
    assert len(direct) == 2
    assert min(_nearest_nonadjacent_boundary_vertex_distance(part) for part in direct) > 1.0


def test_build_city_flat_mesh_dtcc_mesher_handles_touching_holes():
    pytest.importorskip("dtcc_mesher")

    touching_holes = Polygon(
        [(8, 8), (28, 8), (28, 28), (8, 28), (8, 8)],
        [
            [(12, 12), (16, 12), (16, 16), (12, 16), (12, 12)],
            [(16, 16), (22, 16), (22, 24), (16, 24), (16, 16)],
        ],
    )
    city = make_flat_city([make_building(touching_holes, roof_z=10.0)])

    mesh = build_city_flat_mesh(
        city,
        lod=GeometryType.LOD0,
        merge_buildings=False,
        min_building_detail=0.5,
        min_building_area=1.0,
        merge_tolerance=0.5,
        max_mesh_size=5.0,
        min_mesh_angle=20.0,
        report_mesh_quality=False,
        mesher="dtcc_mesher",
    )

    assert mesh.vertices.shape[0] > 0
    assert mesh.faces.shape[0] > 0
    markers = set(np.asarray(mesh.markers, dtype=int))
    assert -2 in markers
    assert 0 in markers


def test_build_city_flat_mesh_dtcc_mesher_handles_nested_building_in_courtyard():
    pytest.importorskip("dtcc_mesher")

    shell = Polygon(
        [(8, 8), (32, 8), (32, 32), (8, 32), (8, 8)],
        [[(14, 14), (26, 14), (26, 26), (14, 26), (14, 14)]],
    )
    nested = box(18, 18, 22, 22)
    city = make_flat_city(
        [
            make_building(shell, roof_z=12.0),
            make_building(nested, roof_z=8.0),
        ]
    )

    mesh = build_city_flat_mesh(
        city,
        lod=GeometryType.LOD0,
        merge_buildings=False,
        min_building_detail=0.5,
        min_building_area=1.0,
        merge_tolerance=0.5,
        max_mesh_size=5.0,
        min_mesh_angle=20.0,
        report_mesh_quality=False,
        mesher="dtcc_mesher",
    )

    assert mesh.vertices.shape[0] > 0
    assert mesh.faces.shape[0] > 0
    assert set(np.unique(mesh.markers)).issuperset({-2, 0, 1})


def test_condition_flat_mesh_building_regions_removes_subscale_edges():
    short_edge_hole = [
        (675446.375, 6581300.46875),
        (675444.28125, 6581294.59375),
        (675440.6875, 6581295.875),
        (675434.5, 6581277.5),
        (675427.6875, 6581279.84375),
        (675426.3125, 6581275.875),
        (675433.09375, 6581273.5),
        (675434.4375, 6581277.375),
        (675434.53125, 6581277.40625),
        (675438.25, 6581276.09375),
        (675432.0, 6581258.4375),
        (675428.21875, 6581259.75),
        (675429.65625, 6581263.84375),
        (675423.0625, 6581266.15625),
        (675420.09375, 6581257.59375),
        (675426.8125, 6581255.25),
        (675427.59375, 6581253.5),
        (675421.875, 6581236.96875),
        (675420.375, 6581236.125),
        (675409.78125, 6581240.03125),
        (675424.125, 6581281.375),
        (675425.5, 6581280.875),
        (675430.5, 6581294.75),
        (675428.9375, 6581295.3125),
        (675433.75, 6581309.15625),
        (675437.53125, 6581307.8125),
        (675444.75, 6581328.625),
        (675454.8125, 6581325.21875),
        (675447.71875, 6581304.25),
        (675443.96875, 6581305.5625),
        (675442.6875, 6581301.78125),
        (675446.375, 6581300.46875),
    ]
    xs = [point[0] for point in short_edge_hole]
    ys = [point[1] for point in short_edge_hole]
    tiny_notch = Polygon(
        [
            (min(xs) - 20.0, min(ys) - 20.0),
            (max(xs) + 20.0, min(ys) - 20.0),
            (max(xs) + 20.0, max(ys) + 20.0),
            (min(xs) - 20.0, max(ys) + 20.0),
            (min(xs) - 20.0, min(ys) - 20.0),
        ],
        [short_edge_hole],
    )

    polygons, markers = meshes_module._condition_flat_mesh_building_regions(
        building_polygons=[tiny_notch],
        building_markers=[7],
        footprint_diagnostics={"output_grid": 0.03125},
        max_mesh_size=10.0,
        min_building_detail=0.5,
        cleaning_diagnostics=False,
    )

    assert markers == [7]
    assert len(polygons) == 1
    assert _minimum_boundary_edge_length(tiny_notch) < 0.2
    assert _minimum_boundary_edge_length(polygons[0]) > 1.0


def test_condition_flat_mesh_building_regions_preserves_mesher_ready_holes():
    touching_holes = Polygon(
        [(8, 8), (28, 8), (28, 28), (8, 28), (8, 8)],
        [
            [(12, 12), (16, 12), (16, 16), (12, 16), (12, 12)],
            [(16, 16), (22, 16), (22, 24), (16, 24), (16, 16)],
        ],
    )

    polygons, markers = meshes_module._condition_flat_mesh_building_regions(
        building_polygons=[touching_holes],
        building_markers=[3],
        footprint_diagnostics={"output_grid": 0.03125},
        max_mesh_size=10.0,
        min_building_detail=0.5,
        cleaning_diagnostics=False,
    )

    assert markers == [3]
    assert polygons
    assert all(not meshes_module._polygon_has_ring_boundary_contacts(polygon) for polygon in polygons)


def test_condition_flat_mesh_ground_polygons_returns_explicit_courtyards():
    building = Polygon(
        [(10, 10), (30, 10), (30, 30), (10, 30), (10, 10)],
        [[(16, 16), (24, 16), (24, 24), (16, 24), (16, 16)]],
    )

    ground_polygons = meshes_module._condition_flat_mesh_ground_polygons(
        bounds=(0.0, 0.0, 40.0, 40.0),
        building_polygons=[building],
        hole_polygons=[],
        max_mesh_size=10.0,
        footprint_diagnostics={"output_grid": 0.03125},
        cleaning_diagnostics=False,
    )

    courtyard = Polygon(building.interiors[0])
    assert any(polygon.covers(courtyard) and courtyard.covers(polygon) for polygon in ground_polygons)


def test_condition_flat_mesh_ground_polygons_courtyards_exclude_nested_buildings():
    shell = Polygon(
        [(10, 10), (34, 10), (34, 34), (10, 34), (10, 10)],
        [[(14, 14), (30, 14), (30, 30), (14, 30), (14, 14)]],
    )
    nested = box(18, 18, 22, 22)

    ground_polygons = meshes_module._condition_flat_mesh_ground_polygons(
        bounds=(0.0, 0.0, 40.0, 40.0),
        building_polygons=[shell, nested],
        hole_polygons=[],
        max_mesh_size=10.0,
        footprint_diagnostics={"output_grid": 0.03125},
        cleaning_diagnostics=False,
    )

    assert any(polygon.covers(Point(16.0, 16.0)) for polygon in ground_polygons)
    assert all(polygon.intersection(nested).area == 0.0 for polygon in ground_polygons)


def test_prepare_surface_ground_regions_preserves_building_holes_for_courtyards():
    building = Polygon(
        [(10, 10), (30, 10), (30, 30), (10, 30), (10, 10)],
        [[(16, 16), (24, 16), (24, 24), (16, 24), (16, 16)]],
    )
    surface = make_surface(building, 12.0)

    (
        active_surfaces,
        directives,
        region_polygons,
        region_markers,
        region_triangle_sizes,
        region_points,
    ) = meshes_module._prepare_surface_ground_regions(
        conditioned_surfaces=[surface],
        conditioned_resolution=[5.0],
        target_lods=[GeometryType.LOD1],
        bounds=(0.0, 0.0, 40.0, 40.0),
        max_mesh_size=10.0,
        min_building_detail=0.5,
        footprint_diagnostics={"output_grid": 0.03125},
        cleaning_diagnostics=False,
        treat_lod0_as_holes=False,
    )

    assert len(active_surfaces) == 1
    assert directives == [1]
    assert region_markers.count(-2) == 2
    assert region_markers.count(0) == 1
    assert region_triangle_sizes == {0: 5.0}

    building_index = region_markers.index(0)
    ground_indices = [index for index, marker in enumerate(region_markers) if marker == -2]
    courtyard = Polygon(building.interiors[0])
    assert len(region_polygons[building_index].interiors) == 1
    assert building.contains(Point(region_points[building_index]))
    assert not courtyard.covers(Point(region_points[building_index]))
    assert all(
        region_polygons[index].intersection(region_polygons[building_index]).area == 0.0
        for index in ground_indices
    )


def test_build_city_flat_mesh_dtcc_mesher_uses_single_coverage_call(monkeypatch):
    calls = {}

    class DummyRawMesh:
        def __init__(self):
            self.points = np.array(
                [
                    [0.0, 0.0],
                    [1.0, 0.0],
                    [0.0, 1.0],
                    [1.0, 1.0],
                ],
                dtype=np.float64,
            )
            self.triangles = np.array([[0, 1, 2], [1, 3, 2]], dtype=np.uint32)
            self.segments = np.empty((0, 2), dtype=np.uint32)
            self.markers = np.array([-2, 0], dtype=np.int32)

    class DummyMesher:
        class Coverage:
            def __init__(self, polygons, markers, *, tolerance=1e-9):
                self.polygons = tuple(polygons)
                self.markers = tuple(markers)
                self.tolerance = tolerance

        class MeshingOptions:
            def __init__(self, *, min_angle, max_edge_length, refine):
                self.min_angle = min_angle
                self.max_edge_length = max_edge_length
                self.refine = refine

        def mesh(self, geometry, *, options):
            calls["polygon_count"] = len(geometry.polygons)
            calls["markers"] = list(geometry.markers)
            calls["min_angle"] = options.min_angle
            calls["max_edge_length"] = options.max_edge_length
            calls["refine"] = options.refine
            return DummyRawMesh()

    monkeypatch.setattr(
        dtcc_mesher_backend_module,
        "_load_dtcc_mesher",
        lambda: DummyMesher(),
    )

    city = make_flat_city([make_building(box(8, 8, 18, 18), roof_z=10.0)])
    mesh = build_city_flat_mesh(
        city,
        lod=GeometryType.LOD0,
        merge_buildings=False,
        min_building_detail=0.5,
        min_building_area=1.0,
        merge_tolerance=0.5,
        max_mesh_size=5.0,
        min_mesh_angle=20.0,
        report_mesh_quality=False,
        mesher="dtcc_mesher",
    )

    assert calls["polygon_count"] == 2
    assert set(calls["markers"]) == {-2, 0}
    assert calls["min_angle"] == 20.0
    assert calls["max_edge_length"] == 5.0
    assert calls["refine"] is True
    assert set(np.asarray(mesh.markers, dtype=int)) == {-1, 0}


def test_build_city_flat_mesh_dtcc_mesher_unrestricted_enables_refinement(monkeypatch):
    calls = {}

    class DummyRawMesh:
        def __init__(self):
            self.points = np.array(
                [
                    [0.0, 0.0],
                    [1.0, 0.0],
                    [0.0, 1.0],
                ],
                dtype=np.float64,
            )
            self.triangles = np.array([[0, 1, 2]], dtype=np.uint32)
            self.segments = np.empty((0, 2), dtype=np.uint32)
            self.markers = np.array([-2], dtype=np.int32)

    class DummyMesher:
        class Coverage:
            def __init__(self, polygons, markers, *, tolerance=1e-9):
                self.polygons = tuple(polygons)
                self.markers = tuple(markers)
                self.tolerance = tolerance

        class MeshingOptions:
            def __init__(self, *, min_angle, max_edge_length, refine):
                self.min_angle = min_angle
                self.max_edge_length = max_edge_length
                self.refine = refine

        def mesh(self, geometry, *, options):
            calls["max_edge_length"] = options.max_edge_length
            calls["refine"] = options.refine
            return DummyRawMesh()

    monkeypatch.setattr(
        dtcc_mesher_backend_module,
        "_load_dtcc_mesher",
        lambda: DummyMesher(),
    )

    city = make_flat_city([make_building(box(8, 8, 18, 18), roof_z=10.0)])
    build_city_flat_mesh(
        city,
        lod=GeometryType.LOD0,
        merge_buildings=False,
        min_building_detail=0.5,
        min_building_area=1.0,
        merge_tolerance=0.5,
        max_mesh_size=None,
        min_mesh_angle=20.0,
        report_mesh_quality=False,
        mesher="dtcc_mesher",
    )

    assert calls["max_edge_length"] is None
    assert calls["refine"] is True


def test_build_city_flat_mesh_dtcc_mesher_uses_explicit_region_points(monkeypatch):
    calls = {}

    class DummyRawMesh:
        def __init__(self):
            self.points = np.array(
                [
                    [0.0, 0.0],
                    [1.0, 0.0],
                    [0.0, 1.0],
                ],
                dtype=np.float64,
            )
            self.triangles = np.array([[0, 1, 2]], dtype=np.uint32)
            self.segments = np.empty((0, 2), dtype=np.uint32)
            self.markers = np.array([-2], dtype=np.int32)

    class DummyMesher:
        class Coverage:
            def __init__(self, polygons, markers, *, tolerance=1e-9):
                self.polygons = tuple(polygons)
                self.markers = tuple(markers)
                self.tolerance = tolerance

            def graph(self, *, max_edge_length=None):
                calls["coverage_graph_max_edge_length"] = max_edge_length
                return DummyMesher.CoverageGraph(
                    np.array([[0.0, 0.0], [1.0, 0.0]], dtype=np.float64),
                    np.array([[0, 1]], dtype=np.uint32),
                    np.array([[0.25, 0.25], [0.75, 0.75]], dtype=np.float64),
                    np.array([-2, 0], dtype=np.int32),
                )

        class CoverageGraph:
            def __init__(self, points, segments, region_points, region_markers):
                self.points = points
                self.segments = segments
                self.region_points = region_points
                self.region_markers = region_markers

        class MeshingOptions:
            def __init__(self, *, min_angle, max_edge_length, refine):
                self.min_angle = min_angle
                self.max_edge_length = max_edge_length
                self.refine = refine

        def mesh(self, geometry, *, options):
            calls["geometry_type"] = type(geometry).__name__
            calls["region_points"] = np.asarray(geometry.region_points, dtype=np.float64)
            calls["region_markers"] = np.asarray(geometry.region_markers, dtype=np.int32)
            return DummyRawMesh()

    monkeypatch.setattr(
        dtcc_mesher_backend_module,
        "_load_dtcc_mesher",
        lambda: DummyMesher(),
    )

    mesh = dtcc_mesher_backend_module.build_city_flat_mesh_with_dtcc_mesher(
        region_polygons=[box(0, 0, 10, 10), box(2, 2, 4, 4)],
        region_markers=[-2, 0],
        region_points=[
            np.array([1.0, 1.0], dtype=np.float64),
            np.array([3.0, 3.0], dtype=np.float64),
        ],
        max_mesh_size=5.0,
        min_mesh_angle=20.0,
    )

    assert calls["coverage_graph_max_edge_length"] == 5.0
    assert calls["geometry_type"] == "CoverageGraph"
    assert np.allclose(calls["region_points"], np.array([[1.0, 1.0], [3.0, 3.0]]))
    assert calls["region_markers"].tolist() == [-2, 0]
    assert mesh.faces.shape[0] == 1


def test_build_city_flat_mesh_triangle_uses_conditioned_coverage(monkeypatch):
    captured = {}

    def fake_condition(*args, **kwargs):
        return [make_surface(box(8, 8, 18, 18), 10.0)], [[0]], [], {"output_grid": 0.03125}

    def fake_condition_coverage_regions(**kwargs):
        captured["coverage_building_count"] = len(kwargs["building_polygons"])
        captured["coverage_markers"] = list(kwargs["building_markers"])
        return [box(0, 0, 80, 80), box(8, 8, 18, 18)], [-2, 7]

    def fake_build_from_coverage(
        *,
        region_polygons,
        region_markers,
        region_points=None,
        bounds,
        max_mesh_size,
        min_mesh_angle,
        backend,
        sort_triangles,
        region_triangle_sizes,
    ):
        captured["region_polygon_count"] = len(region_polygons)
        captured["region_markers"] = list(region_markers)
        captured["bounds"] = bounds
        captured["backend"] = backend
        captured["sort_triangles"] = sort_triangles
        captured["region_triangle_sizes"] = region_triangle_sizes
        return Mesh(
            vertices=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
            faces=np.array([[0, 1, 2]], dtype=int),
            markers=np.array([-2], dtype=int),
        )

    monkeypatch.setattr(meshes_module, "_condition_meshing_footprints", fake_condition)
    monkeypatch.setattr(meshes_module, "_condition_flat_mesh_coverage_regions", fake_condition_coverage_regions)
    monkeypatch.setattr(meshes_module, "build_city_flat_mesh_from_coverage", fake_build_from_coverage)
    monkeypatch.setattr(
        meshes_module,
        "resolve_2d_mesher",
        lambda mesher=None: mesher or "triangle",
    )

    city = make_flat_city([make_building(box(8, 8, 18, 18), roof_z=10.0)])
    mesh = build_city_flat_mesh(
        city,
        lod=GeometryType.LOD0,
        merge_buildings=False,
        min_building_detail=0.5,
        min_building_area=1.0,
        merge_tolerance=0.5,
        max_mesh_size=5.0,
        min_mesh_angle=20.0,
        report_mesh_quality=False,
        mesher="triangle",
    )

    assert captured["coverage_building_count"] == 1
    assert captured["coverage_markers"] == [0]
    assert captured["region_polygon_count"] == 2
    assert captured["region_markers"] == [-2, 7]
    assert captured["bounds"] == (0.0, 0.0, 80.0, 80.0)
    assert captured["backend"] == "triangle"
    assert captured["sort_triangles"] is True
    assert captured["region_triangle_sizes"] is None
    assert mesh.faces.shape[0] >= 1


def test_build_city_flat_mesh_propagates_runtime_mesher_errors(monkeypatch):
    def fake_condition(*args, **kwargs):
        return [make_surface(box(8, 8, 18, 18), 10.0)], [[0]], [], {"output_grid": 0.03125}

    def fake_condition_coverage_regions(**kwargs):
        return [box(0, 0, 80, 80), box(8, 8, 18, 18)], [-2, 7]

    def fake_build_from_coverage(**kwargs):
        raise RuntimeError("internal triangulation error")

    monkeypatch.setattr(meshes_module, "_condition_meshing_footprints", fake_condition)
    monkeypatch.setattr(meshes_module, "_condition_flat_mesh_coverage_regions", fake_condition_coverage_regions)
    monkeypatch.setattr(meshes_module, "build_city_flat_mesh_from_coverage", fake_build_from_coverage)
    monkeypatch.setattr(
        meshes_module,
        "resolve_2d_mesher",
        lambda mesher=None: mesher or "triangle",
    )

    city = make_flat_city([make_building(box(8, 8, 18, 18), roof_z=10.0)])
    with pytest.raises(RuntimeError, match="internal triangulation error"):
        build_city_flat_mesh(
            city,
            lod=GeometryType.LOD0,
            merge_buildings=False,
            min_building_detail=0.5,
            min_building_area=1.0,
            merge_tolerance=0.5,
            max_mesh_size=None,
            min_mesh_angle=20.0,
            report_mesh_quality=False,
            mesher="dtcc_mesher",
        )


def test_build_city_flat_mesh_allows_empty_conditioned_footprints(monkeypatch):
    captured = {}

    def fake_condition(*args, **kwargs):
        return [], [], [], {"output_count": 0}

    class DummyCppMesh:
        def from_cpp(self):
            return Mesh(
                vertices=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
                faces=np.array([[0, 1, 2]], dtype=int),
            )

    def fake_build_flat_mesh(
        building_polygons,
        holes,
        subdomain_resolution,
        xmin,
        ymin,
        xmax,
        ymax,
        max_mesh_size,
        min_mesh_angle,
        sort_triangles,
        backend,
    ):
        captured["building_polygons"] = building_polygons
        captured["subdomain_resolution"] = subdomain_resolution
        captured["backend"] = backend
        return DummyCppMesh()

    monkeypatch.setattr(meshes_module, "_condition_meshing_footprints", fake_condition)
    monkeypatch.setattr(meshes_module._dtcc_builder, "build_city_flat_mesh", fake_build_flat_mesh)

    city = make_flat_city([])
    mesh = build_city_flat_mesh(
        city,
        lod=GeometryType.LOD0,
        merge_buildings=True,
        min_building_detail=0.5,
        min_building_area=1.0,
        merge_tolerance=0.5,
        max_mesh_size=5.0,
        min_mesh_angle=20.0,
        report_mesh_quality=False,
        mesher="spade",
    )

    assert captured["building_polygons"] == []
    assert captured["subdomain_resolution"] == []
    assert captured["backend"] == "spade"
    assert mesh.faces.shape[0] == 1


def test_build_city_flat_mesh_forwards_triangle_backend(monkeypatch):
    captured = {}

    def fake_condition(*args, **kwargs):
        return [], [], [], {"output_grid": 0.03125}

    class DummyCppMesh:
        def from_cpp(self):
            return Mesh(
                vertices=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
                faces=np.array([[0, 1, 2]], dtype=int),
            )

    def fake_build_flat_mesh(
        building_polygons,
        holes,
        subdomain_resolution,
        xmin,
        ymin,
        xmax,
        ymax,
        max_mesh_size,
        min_mesh_angle,
        sort_triangles,
        backend,
    ):
        captured["backend"] = backend
        return DummyCppMesh()

    monkeypatch.setattr(meshes_module, "_condition_meshing_footprints", fake_condition)
    monkeypatch.setattr(meshes_module, "resolve_2d_mesher", lambda mesher=None: mesher or "triangle")
    monkeypatch.setattr(meshes_module._dtcc_builder, "build_city_flat_mesh", fake_build_flat_mesh)

    city = make_flat_city([])
    mesh = build_city_flat_mesh(
        city,
        lod=GeometryType.LOD0,
        merge_buildings=True,
        min_building_detail=0.5,
        min_building_area=1.0,
        merge_tolerance=0.5,
        max_mesh_size=5.0,
        min_mesh_angle=20.0,
        report_mesh_quality=False,
        mesher="triangle",
    )

    assert captured["backend"] == "triangle"
    assert mesh.faces.shape[0] == 1


def test_build_city_flat_mesh_marks_halos_for_triangle_backend(monkeypatch):
    def fake_condition(*args, **kwargs):
        return [make_surface(box(8, 8, 18, 18), 10.0)], [[0]], [], {"output_grid": 0.03125}

    def fake_condition_coverage_regions(**kwargs):
        return [box(0, 0, 80, 80), box(8, 8, 18, 18)], [-2, 7]

    def fake_build_from_coverage(**kwargs):
        return Mesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [1.0, 1.0, 0.0],
                ]
            ),
            faces=np.array([[0, 1, 2], [1, 3, 2]], dtype=int),
            markers=np.array([-2, 7], dtype=int),
        )

    monkeypatch.setattr(meshes_module, "_condition_meshing_footprints", fake_condition)
    monkeypatch.setattr(meshes_module, "_condition_flat_mesh_coverage_regions", fake_condition_coverage_regions)
    monkeypatch.setattr(meshes_module, "build_city_flat_mesh_from_coverage", fake_build_from_coverage)
    monkeypatch.setattr(
        meshes_module,
        "resolve_2d_mesher",
        lambda mesher=None: mesher or "triangle",
    )

    city = make_flat_city([make_building(box(8, 8, 18, 18), roof_z=10.0)])
    mesh = build_city_flat_mesh(
        city,
        lod=GeometryType.LOD0,
        merge_buildings=False,
        min_building_detail=0.5,
        min_building_area=1.0,
        merge_tolerance=0.5,
        max_mesh_size=5.0,
        min_mesh_angle=20.0,
        report_mesh_quality=False,
        mesher="triangle",
    )

    assert np.array_equal(np.asarray(mesh.markers, dtype=int), np.array([-1, 7], dtype=int))


def test_builder_flat_mesh_backend_remaps_canonical_region_markers(monkeypatch):
    captured = {}

    class DummyCppMesh:
        def from_cpp(self):
            return Mesh(
                vertices=np.array(
                    [
                        [0.0, 0.0, 0.0],
                        [1.0, 0.0, 0.0],
                        [0.0, 1.0, 0.0],
                    ]
                ),
                faces=np.array([[0, 1, 2], [0, 2, 1], [1, 2, 0]], dtype=int),
                markers=np.array([0, -1, -2], dtype=int),
            )

    def fake_create_builder_polygon(polygon):
        return polygon

    def fake_build_flat_mesh(
        building_polygons,
        holes,
        subdomain_resolution,
        xmin,
        ymin,
        xmax,
        ymax,
        max_mesh_size,
        min_mesh_angle,
        sort_triangles,
        backend,
    ):
        captured["building_polygon_count"] = len(building_polygons)
        captured["subdomain_resolution"] = list(subdomain_resolution)
        captured["backend"] = backend
        return DummyCppMesh()

    monkeypatch.setattr(flat_mesh_backends_module, "create_builder_polygon", fake_create_builder_polygon)
    monkeypatch.setattr(flat_mesh_backends_module._dtcc_builder, "build_city_flat_mesh", fake_build_flat_mesh)

    mesh = flat_mesh_backends_module.build_city_flat_mesh_with_builder_backend(
        region_polygons=[box(0, 0, 10, 10), box(2, 2, 4, 4)],
        region_markers=[-2, 7],
        bounds=(0.0, 0.0, 10.0, 10.0),
        max_mesh_size=5.0,
        min_mesh_angle=20.0,
        backend="triangle",
    )

    assert captured["building_polygon_count"] == 1
    assert captured["subdomain_resolution"] == []
    assert captured["backend"] == "triangle"
    assert np.array_equal(np.asarray(mesh.markers, dtype=int), np.array([7, -1, -2], dtype=int))


def test_native_triangle_shared_vertex_regions_survive_python_conversion():
    if "triangle" not in _dtcc_builder.triangulation_backends():
        pytest.skip("Triangle backend not built")

    script = textwrap.dedent(
        """
        from shapely.geometry import Polygon
        from dtcc_core.builder import _dtcc_builder
        from dtcc_core.builder.model_conversion import create_builder_polygon

        polygon_a = create_builder_polygon(
            Polygon([(2, 2), (10, 2), (10, 10), (2, 10), (2, 2)])
        )
        polygon_b = create_builder_polygon(
            Polygon([(10, 10), (18, 10), (18, 18), (10, 18), (10, 10)])
        )
        raw_mesh = _dtcc_builder.build_city_flat_mesh(
            [polygon_a, polygon_b],
            [],
            [],
            0.0,
            0.0,
            20.0,
            20.0,
            5.0,
            20.0,
            False,
            "triangle",
        )
        mesh = raw_mesh.from_cpp()
        print(mesh.vertices.shape, mesh.faces.shape, sorted(set(mesh.markers.tolist())))
        """
    )

    result = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
        cwd="/Users/logg/scratch/dtcc/dtcc-core",
    )

    assert result.returncode == 0, result.stdout + "\n" + result.stderr


def test_native_triangle_build_city_flat_mesh_handles_internal_loops_and_holes():
    if "triangle" not in _dtcc_builder.triangulation_backends():
        pytest.skip("Triangle backend not built")

    building = create_builder_polygon(
        Polygon(
            [(10, 10), (25, 10), (25, 25), (10, 25), (10, 10)],
            [[(14, 14), (20, 14), (20, 20), (14, 20), (14, 14)]],
        )
    )
    explicit_hole = create_builder_polygon(
        Polygon([(30, 30), (40, 30), (40, 40), (30, 40), (30, 30)])
    )

    raw_mesh = _dtcc_builder.build_city_flat_mesh(
        [building],
        [explicit_hole],
        [5.0],
        0.0,
        0.0,
        100.0,
        100.0,
        10.0,
        20.0,
        False,
        "triangle",
    )
    mesh = builder_mesh_to_mesh(raw_mesh)

    assert mesh.vertices.shape[0] > 0
    assert mesh.faces.shape[0] > 0
    assert set(np.asarray(mesh.markers, dtype=int)) >= {-2, 0}


def test_build_city_surface_mesh_reduces_lod_from_source_map(monkeypatch):
    surfaces = [
        make_surface(box(10, 10, 20, 20), 10.0),
        make_surface(box(30, 10, 40, 20), 12.0),
    ]
    source_map = [[0, 1], [2]]
    resolutions = [4.0, 5.0]

    captured = {}

    def fake_condition(*args, **kwargs):
        return surfaces, source_map, resolutions, {"output_count": 2}

    class DummyCppMesh:
        def from_cpp(self):
            return Mesh(
                vertices=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
                faces=np.array([[0, 1, 2]], dtype=int),
            )

    def fake_condition_building_regions(**kwargs):
        return kwargs["building_polygons"], kwargs["building_markers"], [[0], [1]]

    def fake_condition_ground_polygons(**kwargs):
        return [box(0, 0, 80, 80)]

    def fake_build_from_coverage(**kwargs):
        captured["backend"] = kwargs["backend"]
        captured["region_markers"] = kwargs["region_markers"]
        captured["region_triangle_sizes"] = kwargs["region_triangle_sizes"]
        return Mesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                ]
            ),
            faces=np.array([[0, 1, 2]], dtype=int),
            markers=np.array([0], dtype=int),
        )

    def fake_build_terrain_surface_mesh_from_ground_mesh(*args, **kwargs):
        return object()

    def fake_build_surface_mesh_from_terrain_mesh(
        building_surfaces,
        building_lod_switches,
        terrain_mesh,
        smoothing,
        merge_meshes,
    ):
        captured["lod_switches"] = building_lod_switches
        return [DummyCppMesh()]

    monkeypatch.setattr(meshes_module, "_condition_meshing_footprints", fake_condition)
    monkeypatch.setattr(
        meshes_module,
        "_condition_flat_mesh_building_regions_with_sources",
        fake_condition_building_regions,
    )
    monkeypatch.setattr(
        meshes_module,
        "_condition_flat_mesh_ground_polygons",
        fake_condition_ground_polygons,
    )
    monkeypatch.setattr(meshes_module, "build_city_flat_mesh_from_coverage", fake_build_from_coverage)
    monkeypatch.setattr(
        meshes_module,
        "resolve_2d_mesher",
        lambda mesher=None: mesher or "triangle",
    )
    monkeypatch.setattr(
        meshes_module._dtcc_builder,
        "build_terrain_surface_mesh_from_ground_mesh",
        fake_build_terrain_surface_mesh_from_ground_mesh,
    )
    monkeypatch.setattr(
        meshes_module._dtcc_builder,
        "build_city_surface_mesh_from_terrain_mesh",
        fake_build_surface_mesh_from_terrain_mesh,
    )

    city = make_flat_city(
        [
            make_building(box(10, 10, 20, 20), roof_z=10.0, lod1_polygon=box(10, 10, 20, 20)),
            make_building(box(20, 10, 30, 20), roof_z=11.0),
            make_building(box(30, 10, 40, 20), roof_z=12.0, lod2_polygon=box(30, 10, 40, 20)),
        ]
    )

    mesh = build_city_surface_mesh(
        city,
        lod=[GeometryType.LOD1, GeometryType.LOD0, GeometryType.LOD2],
        merge_buildings=True,
        min_building_detail=0.5,
        min_building_area=1.0,
        merge_tolerance=0.5,
        building_mesh_triangle_size=6.0,
        max_mesh_size=8.0,
        min_mesh_angle=20.0,
        merge_meshes=True,
        report_mesh_quality=False,
        mesher="triangle",
    )

    assert mesh.faces.shape[0] >= 1
    assert captured["lod_switches"] == [0, 2]
    assert captured["backend"] == "triangle"
    assert captured["region_markers"] == [-2, 0, 1]
    assert captured["region_triangle_sizes"] == {0: 4.0, 1: 5.0}


def test_build_city_surface_mesh_uses_raw_ground_markers(monkeypatch):
    city = make_flat_city([make_building(box(10, 10, 20, 20), roof_z=10.0)])
    terrain = city.terrain
    terrain_raster = terrain.raster
    conditioned_surface = make_surface(box(10, 10, 20, 20), 10.0)
    diagnostics = {"output_grid": 0.25}
    captured = {}

    def fake_prepare(*args, **kwargs):
        return terrain, terrain_raster, [conditioned_surface], [[0]], [4.0], diagnostics

    def fake_prepare_regions(**kwargs):
        return (
            [conditioned_surface],
            [1],
            [box(0, 0, 80, 80), box(10, 10, 20, 20)],
            [-2, 0],
            {0: 4.0},
            [
                np.array([40.0, 40.0], dtype=np.float64),
                np.array([15.0, 15.0], dtype=np.float64),
            ],
        )

    def fake_build_ground(**kwargs):
        captured["add_halo_markers"] = kwargs["add_halo_markers"]
        return (
            Mesh(
                vertices=np.array(
                    [
                        [0.0, 0.0, 0.0],
                        [1.0, 0.0, 0.0],
                        [0.0, 1.0, 0.0],
                    ]
                ),
                faces=np.array([[0, 1, 2]], dtype=int),
                markers=np.array([0], dtype=int),
            ),
            "dtcc_mesher",
        )

    def fake_build_surface_from_ground(**kwargs):
        return Mesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                ]
            ),
            faces=np.array([[0, 1, 2]], dtype=int),
            markers=np.array([0], dtype=int),
        )

    monkeypatch.setattr(meshes_module, "_prepare_city_meshing_inputs", fake_prepare)
    monkeypatch.setattr(meshes_module, "_prepare_surface_ground_regions", fake_prepare_regions)
    monkeypatch.setattr(meshes_module, "_build_ground_mesh_from_coverage", fake_build_ground)
    monkeypatch.setattr(
        meshes_module,
        "_build_city_surface_mesh_from_ground_mesh",
        fake_build_surface_from_ground,
    )

    mesh = build_city_surface_mesh(
        city,
        lod=GeometryType.LOD0,
        merge_buildings=False,
        min_building_detail=0.0,
        min_building_area=1.0,
        merge_tolerance=0.0,
        max_mesh_size=6.0,
        min_mesh_angle=20.0,
        report_mesh_quality=False,
        mesher="dtcc_mesher",
    )

    assert mesh.faces.shape[0] == 1
    assert captured["add_halo_markers"] is False


def test_split_ground_mesh_building_components_splits_vertex_touching_patches():
    ground_mesh = Mesh(
        vertices=np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [1.0, 1.0, 0.0],
                [2.0, 1.0, 0.0],
            ]
        ),
        faces=np.array(
            [
                [0, 1, 2],
                [2, 3, 4],
            ],
            dtype=int,
        ),
        markers=np.array([0, 0], dtype=int),
    )
    building_surface = make_surface(box(0, 0, 2, 1), 10.0)

    split_mesh, split_surfaces, split_directives = (
        meshes_module._split_ground_mesh_building_components(
            ground_mesh=ground_mesh,
            building_surfaces=[building_surface],
            meshing_directives=[1],
        )
    )

    assert split_mesh.markers.tolist() == [0, 1]
    assert len(split_surfaces) == 2
    assert split_directives == [1, 1]


def test_stabilize_shell_building_polygon_removes_low_clearance_hole():
    polygon = Polygon(
        [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)],
        [
            [
                (2.0, 2.0),
                (8.0, 2.0),
                (8.0, 8.0),
                (5.0002, 8.0),
                (5.0, 7.9998),
                (2.0, 8.0),
            ]
        ],
    )

    stabilized, removed = meshes_module._stabilize_shell_building_polygon(
        polygon,
        min_hole_clearance=0.01,
    )

    assert removed == 1
    assert len(stabilized.interiors) == 0


def test_split_weakly_pinched_shell_polygon_splits_exterior_revisit():
    polygon = Polygon(
        [
            (674614.125, 6580137.40625),
            (674617.6249773939, 6580141.031226587),
            (674613.874953289, 6580144.875001294),
            (674617.2812514388, 6580148.000045467),
            (674621.0625469746, 6580144.06249805),
            (674617.6250215588, 6580141.031225638),
            (674627.0, 6580130.4375),
            (674651.09375, 6580152.65625),
            (674634.6875, 6580170.53125),
            (674620.5625, 6580186.0625),
            (674598.21875, 6580165.59375),
            (674586.84375, 6580154.3125),
            (674600.3125, 6580139.65625),
            (674616.78125, 6580121.0),
            (674623.28125, 6580127.03125),
        ]
    )

    parts = meshes_module._split_weakly_pinched_shell_polygon(
        polygon,
        pinch_tolerance=0.05,
    )

    assert len(parts) == 2
    assert all(part.minimum_clearance > 1.0 for part in parts)
    assert all(not meshes_module._polygon_has_ring_boundary_contacts(part) for part in parts)


def test_build_city_volume_mesh_uses_shared_surface_pipeline(monkeypatch):
    city = make_flat_city([make_building(box(10, 10, 20, 20), roof_z=10.0)])
    terrain = city.terrain
    terrain_raster = terrain.raster
    conditioned_surface = make_surface(box(10, 10, 20, 20), 10.0)
    diagnostics = {"output_grid": 0.25}
    captured = {}

    def fake_prepare(*args, **kwargs):
        return terrain, terrain_raster, [conditioned_surface], [[0]], [4.0], diagnostics

    def fake_prepare_regions(**kwargs):
        captured["target_lods"] = kwargs["target_lods"]
        captured["conditioned_resolution"] = kwargs["conditioned_resolution"]
        return (
            [conditioned_surface],
            [1],
            [box(0, 0, 80, 80), box(10, 10, 20, 20)],
            [-2, 0],
            {0: 4.0},
            [
                np.array([40.0, 40.0], dtype=np.float64),
                np.array([15.0, 15.0], dtype=np.float64),
            ],
        )

    def fake_build_ground(**kwargs):
        captured["mesher"] = kwargs["mesher"]
        captured["add_halo_markers"] = kwargs["add_halo_markers"]
        return (
            Mesh(
                vertices=np.array(
                    [
                        [0.0, 0.0, 0.0],
                        [1.0, 0.0, 0.0],
                        [0.0, 1.0, 0.0],
                    ]
                ),
                faces=np.array([[0, 1, 2]], dtype=int),
                markers=np.array([0], dtype=int),
            ),
            "spade",
        )

    def fake_build_surface_from_ground(**kwargs):
        captured["surface_called"] = True
        return Mesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                ]
            ),
            faces=np.array([[0, 1, 2]], dtype=int),
            markers=np.array([0], dtype=int),
        )

    def fake_tetgen_build(**kwargs):
        captured["tetgen_mesh_faces"] = len(kwargs["mesh"].faces)
        captured["closure_mesh_faces"] = len(kwargs["closure_mesh"].faces)
        return VolumeMesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                ]
            ),
            cells=np.array([[0, 1, 2, 3]], dtype=int),
        )

    monkeypatch.setattr(meshes_module, "_prepare_city_meshing_inputs", fake_prepare)
    monkeypatch.setattr(meshes_module, "_prepare_surface_ground_regions", fake_prepare_regions)
    monkeypatch.setattr(meshes_module, "_build_ground_mesh_from_coverage", fake_build_ground)
    monkeypatch.setattr(
        meshes_module,
        "_build_city_surface_mesh_from_ground_mesh",
        fake_build_surface_from_ground,
    )
    monkeypatch.setattr(meshes_module, "is_tetgen_available", lambda: True)
    monkeypatch.setattr(meshes_module, "tetgen_build_volume_mesh", fake_tetgen_build)

    volume_mesh = build_city_volume_mesh(
        city,
        lod=GeometryType.LOD0,
        merge_buildings=False,
        min_building_detail=0.0,
        min_building_area=1.0,
        merge_tolerance=0.0,
        max_mesh_size=6.0,
        min_mesh_angle=20.0,
        report_mesh_quality=False,
        mesher="spade",
    )

    assert captured["mesher"] == "spade"
    assert captured["add_halo_markers"] is False
    assert captured["target_lods"] == [GeometryType.LOD1]
    assert captured["conditioned_resolution"] == [4.0]
    assert captured["surface_called"] is True
    assert captured["tetgen_mesh_faces"] == 1
    assert captured["closure_mesh_faces"] == 1
    assert volume_mesh.cells.shape[0] == 1


def test_build_city_surface_mesh_from_ground_mesh_snaps_boundary_vertices(monkeypatch):
    raster = Raster()
    raster.data = np.zeros((80, 80), dtype=float)
    raster.set_bounds(Bounds(0.0, 0.0, 80.0, 80.0))
    ground_mesh = Mesh(
        vertices=np.array(
            [
                [0.05, 0.0, 0.0],
                [80.0, 0.04, 0.0],
                [79.96, 80.0, 0.0],
                [0.0, 79.97, 0.0],
                [40.0, 40.0, 0.0],
            ]
        ),
        faces=np.array([[0, 1, 4], [1, 2, 4], [2, 3, 4]], dtype=int),
        markers=np.array([0, 0, 0], dtype=int),
    )
    captured = {}

    class DummyBuilderMesh:
        def __init__(self, mesh):
            self._mesh = mesh

        def from_cpp(self):
            return self._mesh

    def fake_mesh_to_builder_mesh(mesh):
        captured["aligned_vertices"] = np.array(mesh.vertices, copy=True)
        return DummyBuilderMesh(mesh)

    def fake_build_terrain(builder_mesh, builder_dem, smoothing):
        captured["terrain_builder_input"] = np.array(builder_mesh._mesh.vertices, copy=True)
        return DummyBuilderMesh(
            Mesh(
                vertices=np.array(
                    [
                        [0.0, 0.0, 0.0],
                        [80.0, 0.0, 0.0],
                        [0.0, 80.0, 0.0],
                    ]
                ),
                faces=np.array([[0, 1, 2]], dtype=int),
                markers=np.array([0], dtype=int),
            )
        )

    def fail_city_surface(*args, **kwargs):
        raise AssertionError("Terrain-only surface mesh should not build building surfaces.")

    monkeypatch.setattr(meshes_module, "mesh_to_builder_mesh", fake_mesh_to_builder_mesh)
    monkeypatch.setattr(meshes_module, "raster_to_builder_gridfield", lambda raster: "grid")
    monkeypatch.setattr(
        meshes_module._dtcc_builder,
        "build_terrain_surface_mesh_from_ground_mesh",
        fake_build_terrain,
    )
    monkeypatch.setattr(
        meshes_module._dtcc_builder,
        "build_city_surface_mesh_from_terrain_mesh",
        fail_city_surface,
    )

    mesh = meshes_module._build_city_surface_mesh_from_ground_mesh(
        ground_mesh=ground_mesh,
        terrain_raster=raster,
        building_surfaces=[],
        meshing_directives=[],
        smoothing=0,
        merge_meshes=True,
    )

    assert np.allclose(captured["aligned_vertices"][0], [0.0, 0.0, 0.0])
    assert np.allclose(captured["aligned_vertices"][1], [80.0, 0.0, 0.0])
    assert np.allclose(captured["aligned_vertices"][2], [80.0, 80.0, 0.0])
    assert np.allclose(captured["aligned_vertices"][3], [0.0, 80.0, 0.0])
    assert np.allclose(captured["aligned_vertices"][4], [40.0, 40.0, 0.0])
    assert np.allclose(captured["terrain_builder_input"], captured["aligned_vertices"])
    assert mesh.faces.shape[0] == 1


def test_build_city_volume_mesh_keeps_requested_tetgen_switches(monkeypatch):
    city = make_flat_city([make_building(box(10, 10, 20, 20), roof_z=10.0)])
    terrain = city.terrain
    terrain_raster = terrain.raster
    conditioned_surface = make_surface(box(10, 10, 20, 20), 10.0)
    diagnostics = {"output_grid": 0.25}
    captured = {}

    def fake_prepare(*args, **kwargs):
        return terrain, terrain_raster, [conditioned_surface], [[0]], [4.0], diagnostics

    def fake_prepare_regions(**kwargs):
        return (
            [conditioned_surface],
            [1],
            [box(0, 0, 80, 80), box(10, 10, 20, 20)],
            [-2, 0],
            {0: 4.0},
            [
                np.array([40.0, 40.0], dtype=np.float64),
                np.array([15.0, 15.0], dtype=np.float64),
            ],
        )

    def fake_build_ground(**kwargs):
        return (
            Mesh(
                vertices=np.array(
                    [
                        [0.0, 0.0, 0.0],
                        [10.0, 0.0, 0.0],
                        [10.0, 10.0, 0.0],
                        [0.0, 10.0, 0.0],
                    ]
                ),
                faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
                markers=np.array([0, 0], dtype=int),
            ),
            "dtcc_mesher",
        )

    def fake_build_surface_from_ground(**kwargs):
        return Mesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [10.0, 0.0, 0.0],
                    [10.0, 10.0, 0.0],
                    [0.0, 10.0, 0.0],
                ]
            ),
            faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
            markers=np.array([0, 0], dtype=int),
        )

    def fake_tetgen_build(**kwargs):
        captured["tetgen_switches"] = dict(kwargs["switches_params"])
        captured["top_cap_backend"] = kwargs["top_cap_backend"]
        captured["top_cap_max_mesh_size"] = kwargs["top_cap_max_mesh_size"]
        captured["top_cap_min_mesh_angle"] = kwargs["top_cap_min_mesh_angle"]
        return VolumeMesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                ]
            ),
            cells=np.array([[0, 1, 2, 3]], dtype=int),
        )

    monkeypatch.setattr(meshes_module, "_prepare_city_meshing_inputs", fake_prepare)
    monkeypatch.setattr(meshes_module, "_prepare_surface_ground_regions", fake_prepare_regions)
    monkeypatch.setattr(meshes_module, "_build_ground_mesh_from_coverage", fake_build_ground)
    monkeypatch.setattr(
        meshes_module,
        "_build_city_surface_mesh_from_ground_mesh",
        fake_build_surface_from_ground,
    )
    monkeypatch.setattr(meshes_module, "is_tetgen_available", lambda: True)
    monkeypatch.setattr(meshes_module, "tetgen_build_volume_mesh", fake_tetgen_build)

    build_city_volume_mesh(
        city,
        lod=GeometryType.LOD0,
        merge_buildings=False,
        min_building_detail=0.0,
        min_building_area=1.0,
        merge_tolerance=0.0,
        max_mesh_size=6.0,
        min_mesh_angle=20.0,
        report_mesh_quality=False,
        mesher="dtcc_mesher",
        tetgen_switches={"quality": (1.6, 25.0)},
    )

    assert captured["tetgen_switches"]["quality"] == (1.6, 25.0)
    assert captured["tetgen_switches"]["preserve_surface"] is False
    assert captured["tetgen_switches"]["max_added_points"] is None


def test_build_city_volume_mesh_allows_empty_conditioned_footprints(monkeypatch):
    city = make_flat_city([])
    terrain = city.terrain
    terrain_raster = terrain.raster
    captured = {}

    def fake_prepare(*args, **kwargs):
        return terrain, terrain_raster, [], [], [], {"output_grid": 0.25}

    def fake_prepare_regions(**kwargs):
        captured["conditioned_surfaces"] = kwargs["conditioned_surfaces"]
        return (
            [],
            [],
            [box(0, 0, 80, 80)],
            [-2],
            {},
            [np.array([40.0, 40.0], dtype=np.float64)],
        )

    def fake_build_ground(**kwargs):
        captured["region_markers"] = list(kwargs["region_markers"])
        return (
            Mesh(
                vertices=np.array(
                    [
                        [0.0, 0.0, 0.0],
                        [80.0, 0.0, 0.0],
                        [80.0, 80.0, 0.0],
                        [0.0, 80.0, 0.0],
                    ]
                ),
                faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
                markers=np.array([-2, -2], dtype=int),
            ),
            "dtcc_mesher",
        )

    def fake_build_surface_from_ground(**kwargs):
        captured["building_surfaces"] = list(kwargs["building_surfaces"])
        return Mesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [80.0, 0.0, 0.0],
                    [80.0, 80.0, 0.0],
                    [0.0, 80.0, 0.0],
                ]
            ),
            faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
            markers=np.array([-1, -1], dtype=int),
        )

    def fake_tetgen_build(**kwargs):
        captured["tetgen_faces"] = len(kwargs["mesh"].faces)
        return VolumeMesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                ]
            ),
            cells=np.array([[0, 1, 2, 3]], dtype=int),
        )

    monkeypatch.setattr(meshes_module, "_prepare_city_meshing_inputs", fake_prepare)
    monkeypatch.setattr(meshes_module, "_prepare_surface_ground_regions", fake_prepare_regions)
    monkeypatch.setattr(meshes_module, "_build_ground_mesh_from_coverage", fake_build_ground)
    monkeypatch.setattr(
        meshes_module,
        "_build_city_surface_mesh_from_ground_mesh",
        fake_build_surface_from_ground,
    )
    monkeypatch.setattr(meshes_module, "is_tetgen_available", lambda: True)
    monkeypatch.setattr(meshes_module, "tetgen_build_volume_mesh", fake_tetgen_build)

    volume_mesh = build_city_volume_mesh(
        city,
        lod=GeometryType.LOD0,
        merge_buildings=False,
        min_building_detail=0.0,
        min_building_area=1.0,
        merge_tolerance=0.0,
        max_mesh_size=6.0,
        min_mesh_angle=20.0,
        report_mesh_quality=False,
        mesher="dtcc_mesher",
    )

    assert captured["conditioned_surfaces"] == []
    assert captured["region_markers"] == [-2]
    assert captured["building_surfaces"] == []
    assert captured["tetgen_faces"] == 2
    assert volume_mesh.cells.shape[0] == 1


def test_build_city_volume_mesh_respects_explicit_tetgen_switches_for_dtcc_mesher(monkeypatch):
    city = make_flat_city([make_building(box(10, 10, 20, 20), roof_z=10.0)])
    terrain = city.terrain
    terrain_raster = terrain.raster
    conditioned_surface = make_surface(box(10, 10, 20, 20), 10.0)
    diagnostics = {"output_grid": 0.25}
    captured = {}

    def fake_prepare(*args, **kwargs):
        return terrain, terrain_raster, [conditioned_surface], [[0]], [4.0], diagnostics

    def fake_prepare_regions(**kwargs):
        return (
            [conditioned_surface],
            [1],
            [box(0, 0, 80, 80), box(10, 10, 20, 20)],
            [-2, 0],
            {0: 4.0},
            [
                np.array([40.0, 40.0], dtype=np.float64),
                np.array([15.0, 15.0], dtype=np.float64),
            ],
        )

    def fake_build_ground(**kwargs):
        return (
            Mesh(
                vertices=np.array(
                    [
                        [0.0, 0.0, 0.0],
                        [10.0, 0.0, 0.0],
                        [10.0, 10.0, 0.0],
                        [0.0, 10.0, 0.0],
                    ]
                ),
                faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
                markers=np.array([0, 0], dtype=int),
            ),
            "dtcc_mesher",
        )

    def fake_build_surface_from_ground(**kwargs):
        return Mesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [10.0, 0.0, 0.0],
                    [10.0, 10.0, 0.0],
                    [0.0, 10.0, 0.0],
                ]
            ),
            faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
            markers=np.array([0, 0], dtype=int),
        )

    def fake_tetgen_build(**kwargs):
        captured["tetgen_switches"] = dict(kwargs["switches_params"])
        captured["top_cap_backend"] = kwargs["top_cap_backend"]
        captured["top_cap_max_mesh_size"] = kwargs["top_cap_max_mesh_size"]
        captured["top_cap_min_mesh_angle"] = kwargs["top_cap_min_mesh_angle"]
        return VolumeMesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                ]
            ),
            cells=np.array([[0, 1, 2, 3]], dtype=int),
        )

    monkeypatch.setattr(meshes_module, "_prepare_city_meshing_inputs", fake_prepare)
    monkeypatch.setattr(meshes_module, "_prepare_surface_ground_regions", fake_prepare_regions)
    monkeypatch.setattr(meshes_module, "_build_ground_mesh_from_coverage", fake_build_ground)
    monkeypatch.setattr(
        meshes_module,
        "_build_city_surface_mesh_from_ground_mesh",
        fake_build_surface_from_ground,
    )
    monkeypatch.setattr(meshes_module, "is_tetgen_available", lambda: True)
    monkeypatch.setattr(meshes_module, "tetgen_build_volume_mesh", fake_tetgen_build)

    build_city_volume_mesh(
        city,
        lod=GeometryType.LOD0,
        merge_buildings=False,
        min_building_detail=0.0,
        min_building_area=1.0,
        merge_tolerance=0.0,
        max_mesh_size=6.0,
        min_mesh_angle=20.0,
        report_mesh_quality=False,
        mesher="dtcc_mesher",
        tetgen_switches={
            "quality": (1.6, 25.0),
            "preserve_surface": False,
            "max_added_points": 1234,
        },
    )

    assert captured["tetgen_switches"]["quality"] == (1.6, 25.0)
    assert captured["tetgen_switches"]["preserve_surface"] is False
    assert captured["tetgen_switches"]["max_added_points"] == 1234
    assert captured["top_cap_backend"] == "dtcc_mesher"
    assert captured["top_cap_max_mesh_size"] == 6.0
    assert captured["top_cap_min_mesh_angle"] == 20.0


def test_build_city_volume_mesh_retries_internal_tetgen_error_with_preserve_surface(
    monkeypatch,
):
    city = make_flat_city([make_building(box(10, 10, 20, 20), roof_z=10.0)])
    terrain = city.terrain
    terrain_raster = terrain.raster
    conditioned_surface = make_surface(box(10, 10, 20, 20), 10.0)
    diagnostics = {"output_grid": 0.25}
    calls = []

    def fake_prepare(*args, **kwargs):
        return terrain, terrain_raster, [conditioned_surface], [[0]], [4.0], diagnostics

    def fake_prepare_regions(**kwargs):
        return (
            [conditioned_surface],
            [1],
            [box(0, 0, 80, 80), box(10, 10, 20, 20)],
            [-2, 0],
            {0: 4.0},
            [
                np.array([40.0, 40.0], dtype=np.float64),
                np.array([15.0, 15.0], dtype=np.float64),
            ],
        )

    def fake_build_ground(**kwargs):
        return (
            Mesh(
                vertices=np.array(
                    [
                        [0.0, 0.0, 0.0],
                        [10.0, 0.0, 0.0],
                        [10.0, 10.0, 0.0],
                        [0.0, 10.0, 0.0],
                    ]
                ),
                faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
                markers=np.array([0, 0], dtype=int),
            ),
            "dtcc_mesher",
        )

    def fake_build_surface_from_ground(**kwargs):
        return Mesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [10.0, 0.0, 0.0],
                    [10.0, 10.0, 0.0],
                    [0.0, 10.0, 0.0],
                ]
            ),
            faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
            markers=np.array([0, 0], dtype=int),
        )

    def fake_tetgen_build(**kwargs):
        calls.append(
            {
                "switches_params": dict(kwargs["switches_params"]),
                "switches_overrides": (
                    dict(kwargs["switches_overrides"])
                    if kwargs["switches_overrides"] is not None
                    else None
                ),
            }
        )
        if len(calls) == 1:
            raise RuntimeError("TetGen failed (code 2): internal error (report bug)")
        return VolumeMesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                ]
            ),
            cells=np.array([[0, 1, 2, 3]], dtype=int),
        )

    monkeypatch.setattr(meshes_module, "_prepare_city_meshing_inputs", fake_prepare)
    monkeypatch.setattr(meshes_module, "_prepare_surface_ground_regions", fake_prepare_regions)
    monkeypatch.setattr(meshes_module, "_build_ground_mesh_from_coverage", fake_build_ground)
    monkeypatch.setattr(
        meshes_module,
        "_build_city_surface_mesh_from_ground_mesh",
        fake_build_surface_from_ground,
    )
    monkeypatch.setattr(meshes_module, "is_tetgen_available", lambda: True)
    monkeypatch.setattr(meshes_module, "tetgen_build_volume_mesh", fake_tetgen_build)

    volume_mesh = build_city_volume_mesh(
        city,
        lod=GeometryType.LOD0,
        merge_buildings=False,
        min_building_detail=0.0,
        min_building_area=1.0,
        merge_tolerance=0.0,
        max_mesh_size=6.0,
        min_mesh_angle=20.0,
        report_mesh_quality=False,
        mesher="dtcc_mesher",
    )

    assert volume_mesh.cells.shape[0] == 1
    assert len(calls) == 2
    assert calls[0]["switches_params"]["preserve_surface"] is False
    assert calls[0]["switches_overrides"] is None
    assert calls[1]["switches_overrides"]["preserve_surface"] is True


def test_build_city_volume_mesh_saves_tetgen_debug_meshes(monkeypatch, tmp_path):
    city = make_flat_city([make_building(box(10, 10, 20, 20), roof_z=10.0)])
    terrain = city.terrain
    terrain_raster = terrain.raster
    conditioned_surface = make_surface(box(10, 10, 20, 20), 10.0)
    diagnostics = {"output_grid": 0.25}
    captured = {}

    def fake_prepare(*args, **kwargs):
        return terrain, terrain_raster, [conditioned_surface], [[0]], [4.0], diagnostics

    def fake_prepare_regions(**kwargs):
        return (
            [conditioned_surface],
            [1],
            [box(0, 0, 80, 80), box(10, 10, 20, 20)],
            [-2, 0],
            {0: 4.0},
            [
                np.array([40.0, 40.0], dtype=np.float64),
                np.array([15.0, 15.0], dtype=np.float64),
            ],
        )

    def fake_build_ground(**kwargs):
        return (
            Mesh(
                vertices=np.array(
                    [
                        [0.0, 0.0, 0.0],
                        [10.0, 0.0, 0.0],
                        [10.0, 10.0, 0.0],
                        [0.0, 10.0, 0.0],
                    ]
                ),
                faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
                markers=np.array([-2, -2], dtype=int),
            ),
            "dtcc_mesher",
        )

    def fake_build_surface_from_ground(**kwargs):
        return Mesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [10.0, 0.0, 0.0],
                    [10.0, 10.0, 0.0],
                    [0.0, 10.0, 0.0],
                ]
            ),
            faces=np.array([[0, 1, 2], [0, 2, 3]], dtype=int),
            markers=np.array([0, 0], dtype=int),
        )

    def fake_save_debug_meshes(**kwargs):
        captured["debug"] = kwargs
        return {
            "ground": str(tmp_path / "ground.xdmf"),
            "shell": str(tmp_path / "shell.xdmf"),
            "plc": str(tmp_path / "plc.xdmf"),
        }

    def fake_tetgen_build(**kwargs):
        return VolumeMesh(
            vertices=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                ]
            ),
            cells=np.array([[0, 1, 2, 3]], dtype=int),
        )

    monkeypatch.setattr(meshes_module, "_prepare_city_meshing_inputs", fake_prepare)
    monkeypatch.setattr(meshes_module, "_prepare_surface_ground_regions", fake_prepare_regions)
    monkeypatch.setattr(meshes_module, "_build_ground_mesh_from_coverage", fake_build_ground)
    monkeypatch.setattr(
        meshes_module,
        "_build_city_surface_mesh_from_ground_mesh",
        fake_build_surface_from_ground,
    )
    monkeypatch.setattr(meshes_module, "_save_tetgen_debug_meshes", fake_save_debug_meshes)
    monkeypatch.setattr(meshes_module, "is_tetgen_available", lambda: True)
    monkeypatch.setattr(meshes_module, "tetgen_build_volume_mesh", fake_tetgen_build)

    build_city_volume_mesh(
        city,
        lod=GeometryType.LOD0,
        merge_buildings=False,
        min_building_detail=0.0,
        min_building_area=1.0,
        merge_tolerance=0.0,
        max_mesh_size=6.0,
        min_mesh_angle=20.0,
        report_mesh_quality=False,
        mesher="dtcc_mesher",
        tetgen_debug_output_dir=tmp_path,
        tetgen_debug_output_stem="case_055",
    )

    assert captured["debug"]["output_dir"] == tmp_path
    assert captured["debug"]["stem"] == "case_055"
    assert captured["debug"]["top_cap_backend"] == "dtcc_mesher"
    assert captured["debug"]["top_cap_max_mesh_size"] == 6.0
    assert captured["debug"]["top_cap_min_mesh_angle"] == 20.0


def test_build_city_flat_mesh_runs_with_dtcc_mesher():
    pytest.importorskip("dtcc_mesher")

    buildings = [
        make_building(box(8, 8, 18, 18), roof_z=10.0),
        make_building(box(24, 8, 34, 18), roof_z=12.0),
    ]
    city = make_flat_city(buildings)

    mesh = build_city_flat_mesh(
        city,
        lod=GeometryType.LOD0,
        merge_buildings=False,
        min_building_detail=0.0,
        min_building_area=1.0,
        merge_tolerance=0.0,
        max_mesh_size=6.0,
        min_mesh_angle=20.0,
        report_mesh_quality=False,
        mesher="dtcc_mesher",
    )

    assert mesh.vertices.shape[0] > 0
    assert mesh.faces.shape[0] > 0
    assert np.any(mesh.markers >= 0)
    assert np.any(mesh.markers == -2)


def test_build_city_surface_mesh_runs_with_mixed_lod_directives():
    city = make_flat_city(
        [
            make_building(box(8, 8, 18, 18), roof_z=10.0, lod1_polygon=box(8, 8, 18, 18)),
            make_building(box(22, 8, 32, 18), roof_z=12.0, lod2_polygon=box(22, 8, 32, 18)),
        ]
    )

    mesh = build_city_surface_mesh(
        city,
        lod=[GeometryType.LOD1, GeometryType.LOD2],
        merge_buildings=False,
        min_building_detail=0.0,
        min_building_area=1.0,
        merge_tolerance=0.0,
        building_mesh_triangle_size=5.0,
        max_mesh_size=6.0,
        min_mesh_angle=20.0,
        merge_meshes=True,
        report_mesh_quality=False,
    )

    assert mesh.vertices.shape[0] > 0
    assert mesh.faces.shape[0] > 0


def test_legacy_building_wrappers_still_return_buildings():
    buildings = [
        make_building(box(0, 0, 4, 4), roof_z=8.0),
        make_building(box(4.1, 0, 8.1, 4), roof_z=10.0),
    ]

    merged, merged_map = merge_building_footprints(
        buildings,
        max_distance=0.2,
        min_area=0.0,
        return_index_map=True,
    )
    cleaned, cleaned_map = clean_building_footprints(
        buildings,
        clearance=0.2,
        smallest_hole_area=0.0,
        return_index_map=True,
    )
    fixed, fixed_map = fix_building_footprint_clearance(
        buildings,
        clearance=0.2,
        return_index_map=True,
    )
    simplified, simplified_map = simplify_building_footprints(
        buildings,
        tolerance=0.1,
        return_index_map=True,
    )

    assert merged and merged_map
    assert cleaned and cleaned_map
    assert fixed and fixed_map
    assert simplified and simplified_map


@pytest.mark.skipif(not is_tetgen_available(), reason="TetGen is not available")
def test_build_city_volume_mesh_smoke():
    city = make_flat_city(
        [
            make_building(box(10, 10, 18, 18), roof_z=10.0),
            make_building(box(24, 10, 32, 18), roof_z=12.0),
        ]
    )

    volume_mesh = build_city_volume_mesh(
        city,
        lod=GeometryType.LOD0,
        domain_height=40.0,
        max_mesh_size=8.0,
        min_mesh_angle=20.0,
        merge_buildings=True,
        min_building_detail=0.0,
        min_building_area=1.0,
        merge_tolerance=0.0,
        smoothing=0,
        boundary_face_markers=False,
        report_mesh_quality=False,
    )

    assert volume_mesh.vertices.shape[0] > 0
    assert volume_mesh.cells.shape[0] > 0
