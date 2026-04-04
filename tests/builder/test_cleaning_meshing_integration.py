import numpy as np
import pytest
from shapely.geometry import Polygon, box

from dtcc_core.builder import (
    build_city_flat_mesh,
    build_city_surface_mesh,
    build_city_volume_mesh,
    fix_building_footprint_clearance,
    merge_building_footprints,
    simplify_building_footprints,
)
from dtcc_core.builder.building.modify import clean_building_footprints
from dtcc_core.builder.geometry_builders import meshes as meshes_module
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
    ):
        captured["building_polygons"] = building_polygons
        captured["subdomain_resolution"] = subdomain_resolution
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
    assert mesh.faces.shape[0] == 1


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
