import numpy as np
from shapely.geometry import Polygon

from dtcc_core.model import Building, GeometryType, MultiSurface, PointCloud, Surface
from dtcc_core.builder.geometry_builders.lod2 import is_watertight
from dtcc_core.builder.geometry_builders.lod2 import _fit_plane, _ransac_planes
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
    left = np.array([[x, y, 10 + 0.2 * x] for x in range(5) for y in range(5)], dtype=float)
    right = np.array([[x + 6, y, 12 - 0.2 * x] for x in range(5) for y in range(5)], dtype=float)
    points = np.vstack([left, right])

    first = _ransac_planes(points, seed=7)
    second = _ransac_planes(points, seed=7)

    assert [len(plane.inliers) for plane in first] == [len(plane.inliers) for plane in second]
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
