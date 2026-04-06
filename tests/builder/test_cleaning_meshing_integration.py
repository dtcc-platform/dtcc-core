import subprocess
import sys
import textwrap

import numpy as np
import pytest
from shapely.geometry import Polygon, box

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
from dtcc_core.model import Building, Bounds, City, GeometryType, Mesh, Raster, Surface, Terrain


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
    assert diagnostics["mesher_regularized_polygon_count"] == 1
    normalized = surfaces[0].to_polygon(simplify=0.0)
    assert not meshes_module._polygon_has_ring_boundary_contacts(normalized)


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
        bounds,
        max_mesh_size,
        min_mesh_angle,
        backend,
    ):
        captured["region_polygon_count"] = len(region_polygons)
        captured["region_markers"] = list(region_markers)
        captured["bounds"] = bounds
        captured["backend"] = backend
        return Mesh(
            vertices=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
            faces=np.array([[0, 1, 2]], dtype=int),
            markers=np.array([-2], dtype=int),
        )

    monkeypatch.setattr(meshes_module, "_condition_meshing_footprints", fake_condition)
    monkeypatch.setattr(meshes_module, "_condition_flat_mesh_coverage_regions", fake_condition_coverage_regions)
    monkeypatch.setattr(meshes_module, "build_city_flat_mesh_from_coverage", fake_build_from_coverage)

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
    assert mesh.faces.shape[0] == 1


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

    def fake_build_surface_mesh(
        building_surfaces,
        hole_surfaces,
        building_lod_switches,
        building_resolution,
        builder_dem,
        max_mesh_size,
        min_mesh_angle,
        smoothing,
        merge_meshes,
        sort_triangles,
    ):
        captured["lod_switches"] = building_lod_switches
        captured["resolution"] = building_resolution
        return [DummyCppMesh()]

    monkeypatch.setattr(meshes_module, "_condition_meshing_footprints", fake_condition)
    monkeypatch.setattr(meshes_module._dtcc_builder, "build_city_surface_mesh", fake_build_surface_mesh)

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
    )

    assert mesh.faces.shape[0] == 1
    assert captured["lod_switches"] == [0, 2]
    assert captured["resolution"] == [4.0, 5.0]


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
