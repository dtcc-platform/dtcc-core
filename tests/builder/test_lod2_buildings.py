from pathlib import Path

import numpy as np
import pytest
from shapely.geometry import Polygon

import dtcc_core.builder.geometry_builders.lod2 as lod2_module
from dtcc_core.model import Building, City, GeometryType, MultiSurface, PointCloud, Surface
from dtcc_core.builder.geometry_builders.lod2 import is_watertight
from dtcc_core.builder.geometry_builders.lod2 import _fit_plane, _ransac_planes
from dtcc_core.builder.geometry_builders.lod2 import _roof_surfaces_cover_footprint
from dtcc_core.builder.geometry_builders.lod2 import build_lod2_buildings


def _surface(coords):
    return Surface(vertices=np.array(coords, dtype=float))


def _closed_box():
    return MultiSurface(
        surfaces=[
            _surface([[0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]]),
            _surface([[0, 0, 0], [0, 1, 0], [1, 1, 0], [1, 0, 0]]),
            _surface([[0, 0, 0], [1, 0, 0], [1, 0, 1], [0, 0, 1]]),
            _surface([[1, 0, 0], [1, 1, 0], [1, 1, 1], [1, 0, 1]]),
            _surface([[1, 1, 0], [0, 1, 0], [0, 1, 1], [1, 1, 1]]),
            _surface([[0, 1, 0], [0, 0, 0], [0, 0, 1], [0, 1, 1]]),
        ]
    )


def test_watertight_validator_accepts_closed_box():
    assert is_watertight(_closed_box())


def test_watertight_validator_rejects_open_shell():
    shell = _closed_box()
    shell.surfaces.pop()

    assert not is_watertight(shell)


def test_watertight_validator_rejects_non_manifold_edge():
    shell = _closed_box()
    duplicate_top = shell.surfaces[0].copy(geometry_only=True)
    shell.surfaces.append(duplicate_top)

    assert not is_watertight(shell)


def test_fit_plane_predicts_roof_z_values():
    points = np.array(
        [
            [0, 0, 10],
            [1, 0, 11],
            [0, 1, 12],
            [1, 1, 13],
            [2, 1, 14],
        ],
        dtype=float,
    )

    plane = _fit_plane(points)

    assert np.isclose(plane.z_at(2, 2), 16.0)


def test_ransac_planes_is_deterministic_for_two_planes():
    left = np.array([[x, y, 10 + 0.3 * x] for x in range(6) for y in range(6)], dtype=float)
    right = np.array([[x + 20, y, 30 - 0.3 * x] for x in range(6) for y in range(6)], dtype=float)
    points = np.vstack([left, right])

    first = _ransac_planes(points, seed=7)
    second = _ransac_planes(points, seed=7)

    assert [len(plane.inliers) for plane in first] == [len(plane.inliers) for plane in second]
    assert np.allclose(
        [(plane.a, plane.b, plane.c) for plane in first],
        [(plane.a, plane.b, plane.c) for plane in second],
    )
    assert len(first) == 2


def _building_with_footprint(points, *, ground_height=0.0):
    building = Building()
    footprint = Surface()
    footprint.from_polygon(Polygon([(0, 0), (10, 0), (10, 10), (0, 10)]), 10.0)
    building.add_geometry(footprint, GeometryType.LOD0)
    building.add_geometry(PointCloud(points=np.array(points, dtype=float)), GeometryType.POINT_CLOUD)
    building.attributes["ground_height"] = ground_height
    building.attributes["height"] = 10.0 - ground_height
    return building


def _flat_roof_points(z=10.0):
    return [[x, y, z] for x in np.linspace(1, 9, 5) for y in np.linspace(1, 9, 5)]


def test_build_lod2_buildings_creates_watertight_flat_roof():
    building = _building_with_footprint(_flat_roof_points())

    result = build_lod2_buildings([building], build_lod1_fallback=False)

    assert result[0].lod2 is not None
    assert is_watertight(result[0].lod2)
    assert len(result[0].lod2.surfaces) == 6


def test_build_lod2_buildings_skips_sparse_roof_points():
    building = _building_with_footprint(_flat_roof_points()[:6])

    result = build_lod2_buildings([building], build_lod1_fallback=False)

    assert result[0].lod2 is None


def test_build_lod2_buildings_does_not_log_rejection_summary_by_default(monkeypatch):
    messages = []
    monkeypatch.setattr(lod2_module, "info", messages.append, raising=False)
    building = _building_with_footprint(_flat_roof_points()[:6])

    build_lod2_buildings([building])

    assert not any(message.startswith("LOD2 build summary:") for message in messages)
    assert not any(message.startswith("LOD2 rejection summary:") for message in messages)


def test_projected_roof_surfaces_must_cover_footprint():
    footprint = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
    half_roof = _surface([[0, 0, 10], [5, 0, 10], [5, 10, 10], [0, 10, 10]])
    other_half_roof = _surface([[5, 0, 10], [10, 0, 10], [10, 10, 10], [5, 10, 10]])

    assert not _roof_surfaces_cover_footprint(footprint, [half_roof])
    assert _roof_surfaces_cover_footprint(footprint, [half_roof, other_half_roof])


def _gable_roof_points():
    points = []
    for x in np.linspace(1, 4.5, 5):
        for y in np.linspace(1, 9, 5):
            points.append([x, y, 10.0 + 0.5 * x])
    for x in np.linspace(5.5, 9, 5):
        for y in np.linspace(1, 9, 5):
            points.append([x, y, 15.0 - 0.5 * x])
    return points


def _edge_keys_for_surface(surface, tolerance=1e-3):
    keys = []
    for index, vertex in enumerate(surface.vertices):
        start = tuple(np.round(vertex / tolerance).astype(int))
        end = tuple(np.round(surface.vertices[(index + 1) % len(surface.vertices)] / tolerance).astype(int))
        keys.append(tuple(sorted((start, end))))
    return keys


def test_build_lod2_buildings_creates_watertight_gable_with_shared_ridge():
    building = _building_with_footprint(_gable_roof_points())

    result = build_lod2_buildings([building], build_lod1_fallback=False)

    assert result[0].lod2 is not None
    roof_surfaces = result[0].lod2.surfaces[:2]
    shared_edges = set(_edge_keys_for_surface(roof_surfaces[0])).intersection(_edge_keys_for_surface(roof_surfaces[1]))
    assert is_watertight(result[0].lod2)
    assert len(roof_surfaces) == 2
    assert len(shared_edges) == 1


def test_rebuild_false_preserves_existing_lod2():
    building = _building_with_footprint(_flat_roof_points())
    existing = _closed_box()
    building.add_geometry(existing, GeometryType.LOD2)

    build_lod2_buildings([building], rebuild=False, build_lod1_fallback=False)

    assert building.lod2 is existing


def test_rebuild_true_replaces_existing_lod2_when_reconstruction_succeeds():
    building = _building_with_footprint(_flat_roof_points())
    existing = _closed_box()
    building.add_geometry(existing, GeometryType.LOD2)

    build_lod2_buildings([building], rebuild=True, build_lod1_fallback=False)

    assert building.lod2 is not existing
    assert is_watertight(building.lod2)


def test_rebuild_true_removes_existing_lod2_when_reconstruction_fails():
    building = _building_with_footprint(_flat_roof_points()[:6])
    building.add_geometry(_closed_box(), GeometryType.LOD2)

    build_lod2_buildings([building], rebuild=True, build_lod1_fallback=False)

    assert building.lod2 is None


def test_hole_footprint_skips_lod2_and_builds_lod1_fallback():
    building = Building()
    footprint = Surface()
    footprint.from_polygon(
        Polygon(
            [(0, 0), (10, 0), (10, 10), (0, 10)],
            holes=[[(4, 4), (6, 4), (6, 6), (4, 6)]],
        ),
        10.0,
    )
    building.add_geometry(footprint, GeometryType.LOD0)
    building.add_geometry(PointCloud(points=np.array(_flat_roof_points(), dtype=float)), GeometryType.POINT_CLOUD)
    building.attributes["ground_height"] = 0.0

    build_lod2_buildings([building], build_lod1_fallback=True)

    assert building.lod2 is None
    assert building.lod1 is not None


def test_build_lod2_buildings_preserves_existing_lod1_fallback():
    building = _building_with_footprint(_flat_roof_points()[:6])
    existing_lod1 = _closed_box()
    building.add_geometry(existing_lod1, GeometryType.LOD1)

    build_lod2_buildings([building], build_lod1_fallback=True)

    assert building.lod2 is None
    assert building.lod1 is existing_lod1


def test_public_builder_import_exposes_lod2_builder():
    import dtcc_core.builder as builder

    assert builder.build_lod2_buildings is build_lod2_buildings


def test_city_build_lod2_buildings_uses_existing_roof_points_without_recomputing_heights():
    city = City()
    city.add_building(_building_with_footprint(_flat_roof_points()))

    result = city.build_lod2_buildings(calculate_heights=False)

    assert result is city
    assert city.buildings[0].lod2 is not None
    assert is_watertight(city.buildings[0].lod2)


def test_city_build_lod2_buildings_logs_rejection_summary_when_requested(monkeypatch):
    messages = []
    monkeypatch.setattr(lod2_module, "info", messages.append, raising=False)
    city = City()
    city.add_building(_building_with_footprint(_flat_roof_points()))
    city.add_building(_building_with_footprint(_flat_roof_points()[:6]))

    city.build_lod2_buildings(calculate_heights=False, log_rejections=True)

    assert any(
        message == "LOD2 build summary: total=2 lod2=1 fallback=1 skipped_existing=0"
        for message in messages
    )
    assert any(
        message == "LOD2 rejection summary: insufficient_roof_points=1"
        for message in messages
    )


def test_build_lod2_buildings_logs_plane_count_summary_when_requested(monkeypatch):
    messages = []
    monkeypatch.setattr(lod2_module, "info", messages.append, raising=False)
    buildings = [
        _building_with_footprint(_flat_roof_points()),
        _building_with_footprint(_gable_roof_points()),
    ]

    build_lod2_buildings(buildings, build_lod1_fallback=False, log_rejections=True)

    assert any(message == "LOD2 plane count summary: planes_1=1 planes_2=1" for message in messages)


def test_build_lod2_buildings_logs_shell_rejection_sub_reason(monkeypatch):
    messages = []
    monkeypatch.setattr(lod2_module, "info", messages.append, raising=False)
    planes = [
        lod2_module.RoofPlane(0.0, 0.0, 10.0, np.arange(12)),
        lod2_module.RoofPlane(0.0, 0.0, 10.0, np.arange(12, 24)),
    ]
    monkeypatch.setattr(lod2_module, "_ransac_planes", lambda points: planes)
    building = _building_with_footprint(_flat_roof_points())

    build_lod2_buildings([building], build_lod1_fallback=False, log_rejections=True)

    assert any(
        message == "LOD2 rejection summary: split_failed=1"
        for message in messages
    )


def test_build_lod2_buildings_logs_max_plane_coplanar_pair_summary(monkeypatch):
    messages = []
    monkeypatch.setattr(lod2_module, "info", messages.append, raising=False)
    planes = [
        lod2_module.RoofPlane(0.0, 0.0, 10.0, np.arange(8)),
        lod2_module.RoofPlane(0.0, 0.0, 10.1, np.arange(8, 16)),
        lod2_module.RoofPlane(0.3, 0.0, 12.0, np.arange(16, 24)),
        lod2_module.RoofPlane(0.0, 0.3, 14.0, np.arange(24, 32)),
        lod2_module.RoofPlane(-0.3, 0.0, 16.0, np.arange(32, 40)),
        lod2_module.RoofPlane(0.0, -0.3, 18.0, np.arange(40, 48)),
    ]
    monkeypatch.setattr(lod2_module, "_ransac_planes", lambda points: planes)
    building = _building_with_footprint(_flat_roof_points())

    build_lod2_buildings([building], build_lod1_fallback=False, log_rejections=True)

    assert any(
        message == "LOD2 max-plane coplanar pair summary: coplanar_pairs_1=1"
        for message in messages
    )


def test_build_lod2_buildings_logs_inlier_share_recovery_simulation(monkeypatch):
    messages = []
    monkeypatch.setattr(lod2_module, "info", messages.append, raising=False)
    planes = [
        lod2_module.RoofPlane(0.0, 0.0, 10.0, np.arange(240)),
        lod2_module.RoofPlane(0.3, 0.0, 11.0, np.arange(240, 440)),
        lod2_module.RoofPlane(0.0, 0.3, 12.0, np.arange(440, 460)),
        lod2_module.RoofPlane(-0.3, 0.0, 13.0, np.arange(460, 483)),
    ]
    monkeypatch.setattr(lod2_module, "_ransac_planes", lambda points: planes)
    building = _building_with_footprint(_flat_roof_points())

    build_lod2_buildings([building], build_lod1_fallback=False, log_rejections=True)

    assert any(
        message
        == "LOD2 inlier-share recovery simulation: would_recover_at_10pct=1 would_recover_at_15pct=1 would_recover_at_20pct=1"
        for message in messages
    )
    assert any(
        message == "LOD2 trailing plane summary: trailing_planes_at_min_inliers=2"
        for message in messages
    )


@pytest.fixture
def minimal_case_dir():
    return Path(__file__).parent / ".." / "data" / "MinimalCase"


def test_city_build_lod2_buildings_calculates_heights_and_keeps_roof_points(minimal_case_dir):
    city = City()
    city.load_footprints(str(minimal_case_dir / "PropertyMap.shp"))
    city.load_pointcloud(str(minimal_case_dir / "pointcloud.las"))

    result = city.build_lod2_buildings(calculate_heights=True)

    assert result is city
    assert len(city.buildings) == 5
    for building in city.buildings:
        assert building.attributes.get("height") is not None
        assert building.attributes.get("ground_height") is not None
        assert building.point_cloud is not None
        assert building.lod1 is not None or building.lod2 is not None
