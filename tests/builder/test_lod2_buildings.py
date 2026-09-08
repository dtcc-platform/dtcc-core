from collections import Counter
from pathlib import Path

import numpy as np
import pytest
from shapely.geometry import LineString, Point, Polygon

import dtcc_core.builder.geometry_builders.lod2 as lod2_module
from dtcc_core.model import Building, City, GeometryType, MultiSurface, PointCloud, Surface
from dtcc_core.builder.geometry_builders.lod2 import is_watertight
from dtcc_core.builder.geometry_builders.lod2 import _concave_vertex_indices
from dtcc_core.builder.geometry_builders.lod2 import _decompose_footprint
from dtcc_core.builder.geometry_builders.lod2 import _decomposition_shape_reason
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


def test_reconcile_shell_edge_vertices_splits_collinear_surface_edges():
    top = _surface([[0, 0, 1], [2, 0, 1], [2, 1, 1], [0, 1, 1]])
    ground = _surface([[0, 0, 0], [0, 1, 0], [2, 1, 0], [2, 0, 0]])
    split_front = _surface(
        [[0, 0, 1], [1, 0, 1], [2, 0, 1], [2, 0, 0], [1, 0, 0], [0, 0, 0]]
    )
    right = _surface([[2, 0, 0], [2, 1, 0], [2, 1, 1], [2, 0, 1]])
    back = _surface([[2, 1, 0], [0, 1, 0], [0, 1, 1], [2, 1, 1]])
    left = _surface([[0, 1, 0], [0, 0, 0], [0, 0, 1], [0, 1, 1]])
    surfaces = [top, ground, split_front, right, back, left]

    assert all(lod2_module._surface_is_simple(surface) for surface in surfaces)
    assert not is_watertight(MultiSurface(surfaces=surfaces))

    reconciled = lod2_module._reconcile_shell_edge_vertices(surfaces)

    assert is_watertight(MultiSurface(surfaces=reconciled))


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


def test_stepped_flat_height_levels_cluster_by_elevation():
    points = np.asarray(_stepped_flat_roof_points(), dtype=float)
    planes = lod2_module._ransac_planes(points)

    levels, reason = lod2_module._stepped_flat_height_levels(points, planes)

    assert reason is None
    assert levels is not None
    assert [round(level.plane.c, 6) for level in levels] == [10.0, 11.0, 13.0]
    assert all(level.plane.a == 0.0 and level.plane.b == 0.0 for level in levels)
    assert [len(level.point_indices) for level in levels] == [25, 25, 25]


def test_stepped_flat_height_levels_drop_minor_noise():
    points = np.asarray(_stepped_flat_noisy_roof_points(), dtype=float)
    planes = lod2_module._ransac_planes(points)

    levels, reason = lod2_module._stepped_flat_height_levels(points, planes)

    assert reason is None
    assert levels is not None
    assert [round(level.plane.c, 6) for level in levels] == [10.0, 11.0, 13.0]


def test_stepped_flat_height_levels_ignores_minor_sloped_artifacts():
    flat_levels = []
    for x0, z in [(1, 10.0), (11, 13.0), (21, 11.0)]:
        flat_levels.extend(
            [x, y, z]
            for x in np.linspace(x0, x0 + 8, 8)
            for y in np.linspace(1, 9, 10)
        )
    sloped_artifact = [
        [100.0 + index, 0.0, 20.0 + 0.3 * index]
        for index in range(20)
    ]
    points = np.asarray([*flat_levels, *sloped_artifact], dtype=float)
    planes = [
        lod2_module.RoofPlane(0.0, 0.0, 10.0, np.arange(0, 80)),
        lod2_module.RoofPlane(0.0, 0.0, 13.0, np.arange(80, 160)),
        lod2_module.RoofPlane(0.0, 0.0, 11.0, np.arange(160, 240)),
        lod2_module.RoofPlane(0.3, 0.0, 20.0, np.arange(240, 260)),
    ]

    levels, reason = lod2_module._stepped_flat_height_levels(points, planes)

    assert reason is None
    assert levels is not None
    assert [round(level.plane.c, 6) for level in levels] == [10.0, 11.0, 13.0]


def test_stepped_flat_height_levels_ignore_unassigned_roof_noise():
    flat_levels = []
    for x0, z in [(1, 10.0), (11, 13.0), (21, 11.0)]:
        flat_levels.extend(
            [x, y, z]
            for x in np.linspace(x0, x0 + 8, 8)
            for y in np.linspace(1, 9, 10)
        )
    unassigned_noise = [
        [60.0 + index % 14, 60.0 + index // 14, 20.0 + 0.25 * index]
        for index in range(140)
    ]
    points = np.asarray([*flat_levels, *unassigned_noise], dtype=float)
    planes = [
        lod2_module.RoofPlane(0.0, 0.0, 10.0, np.arange(0, 80)),
        lod2_module.RoofPlane(0.0, 0.0, 13.0, np.arange(80, 160)),
        lod2_module.RoofPlane(0.0, 0.0, 11.0, np.arange(160, 240)),
    ]

    levels, reason = lod2_module._stepped_flat_height_levels(points, planes)

    assert reason is None
    assert levels is not None
    assert [round(level.plane.c, 6) for level in levels] == [10.0, 11.0, 13.0]


def test_stepped_flat_height_levels_reject_sloped_major_planes():
    points = np.asarray(_gable_roof_points(), dtype=float)
    planes = lod2_module._ransac_planes(points)

    levels, reason = lod2_module._stepped_flat_height_levels(points, planes)

    assert levels is None
    assert reason == lod2_module.STEPPED_FLAT_SLOPED_MAJOR_EVIDENCE


def test_flat_collapse_plane_accepts_one_dominant_level_with_noise():
    points = np.asarray(_dominant_flat_with_minor_noise_points(), dtype=float)
    planes = [
        lod2_module.RoofPlane(0.0, 0.0, 10.0, np.arange(0, 60)),
        lod2_module.RoofPlane(0.0, 0.0, 10.0, np.arange(60, 100)),
        lod2_module.RoofPlane(0.0, 0.0, 8.5, np.arange(100, 108)),
    ]

    plane = lod2_module._flat_collapse_plane(points, planes)

    assert plane is not None
    assert plane.a == 0.0
    assert plane.b == 0.0
    assert round(plane.c, 6) == 10.0


def test_flat_collapse_plane_rejects_real_two_level_step():
    points = np.asarray(
        [
            *[
                [x, y, 10.0]
                for x in np.linspace(1, 9, 8)
                for y in np.linspace(1, 9, 10)
            ],
            *[
                [x, y, 13.0]
                for x in np.linspace(11, 19, 8)
                for y in np.linspace(1, 9, 10)
            ],
        ],
        dtype=float,
    )
    planes = [
        lod2_module.RoofPlane(0.0, 0.0, 10.0, np.arange(0, 80)),
        lod2_module.RoofPlane(0.0, 0.0, 13.0, np.arange(80, 160)),
        lod2_module.RoofPlane(0.0, 0.0, 10.0, np.arange(0, 20)),
    ]

    assert lod2_module._flat_collapse_plane(points, planes) is None


def test_flat_collapse_plane_rejects_continuous_near_flat_slope():
    points = np.asarray(
        [
            [x, y, 10.0 + 1.2 * (x - 1.0) / 28.0]
            for x in np.linspace(1, 29, 10)
            for y in np.linspace(1, 9, 10)
        ],
        dtype=float,
    )
    planes = [
        lod2_module.RoofPlane(0.04, 0.0, 10.0, np.arange(0, 40)),
        lod2_module.RoofPlane(0.04, 0.0, 10.0, np.arange(40, 80)),
        lod2_module.RoofPlane(0.04, 0.0, 10.0, np.arange(80, 100)),
    ]

    assert lod2_module._flat_collapse_plane(points, planes) is None


def test_stepped_flat_patches_cover_rectangular_strip_levels():
    footprint = Polygon([(0, 0), (30, 0), (30, 10), (0, 10)])
    points = np.asarray(_stepped_flat_roof_points(), dtype=float)
    levels, reason = lod2_module._stepped_flat_height_levels(
        points,
        lod2_module._ransac_planes(points),
    )

    patches, reason = lod2_module._stepped_flat_level_patches(points, footprint, levels)

    assert reason is None
    assert patches is not None
    assert len(patches) == 3
    assert lod2_module._patches_cover_footprint(footprint, patches)
    assert [round(patch.area, 6) for patch in patches] == [100.0, 100.0, 100.0]


def test_stepped_flat_region_patches_cover_corner_levels():
    points, footprint, low_indices, high_indices = _corner_stepped_flat_roof_fixture()
    levels = [
        lod2_module.SteppedFlatLevel(
            lod2_module.RoofPlane(0.0, 0.0, 10.0, low_indices),
            low_indices,
        ),
        lod2_module.SteppedFlatLevel(
            lod2_module.RoofPlane(0.0, 0.0, 13.0, high_indices),
            high_indices,
        ),
    ]

    patches, reason = lod2_module._stepped_flat_level_patches(points, footprint, levels)

    assert reason is None
    assert patches is not None
    assert len(patches) == 2
    assert all(patch.is_valid and patch.area > lod2_module.MIN_PATCH_AREA for patch in patches)
    assert lod2_module._patches_cover_footprint(footprint, patches)
    assert patches[0].covers(Point(3.5, 2.5))
    assert patches[1].covers(Point(4.0, 9.0))
    assert patches[1].covers(Point(16.5, 3.0))


def test_stepped_flat_region_patches_reject_scattered_blob_level():
    points, footprint, low_indices, high_indices = _scattered_blob_stepped_flat_roof_fixture()
    levels = [
        lod2_module.SteppedFlatLevel(
            lod2_module.RoofPlane(0.0, 0.0, 10.0, low_indices),
            low_indices,
        ),
        lod2_module.SteppedFlatLevel(
            lod2_module.RoofPlane(0.0, 0.0, 13.0, high_indices),
            high_indices,
        ),
    ]

    patches, reason = lod2_module._stepped_flat_level_patches(points, footprint, levels)

    assert patches is None
    assert reason == lod2_module.STEPPED_FLAT_PATCH_FAILED


def test_stepped_flat_region_assignment_rejects_overlapping_patches():
    footprint = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
    low_points = [
        [x, y, 10.0]
        for x in np.linspace(1, 2, 5)
        for y in np.linspace(1, 9, 5)
    ]
    high_points = [
        [x, y, 13.0]
        for x in np.linspace(8, 9, 5)
        for y in np.linspace(1, 9, 5)
    ]
    points = np.asarray([*low_points, *high_points], dtype=float)
    levels = [
        lod2_module.SteppedFlatLevel(
            lod2_module.RoofPlane(0.0, 0.0, 10.0, np.arange(0, 25)),
            np.arange(0, 25),
        ),
        lod2_module.SteppedFlatLevel(
            lod2_module.RoofPlane(0.0, 0.0, 13.0, np.arange(25, 50)),
            np.arange(25, 50),
        ),
    ]
    patches = lod2_module._assign_candidate_regions_to_levels(
        points,
        footprint,
        [
            Polygon([(0, 0), (7, 0), (7, 10), (0, 10)]),
            Polygon([(3, 0), (10, 0), (10, 10), (3, 10)]),
        ],
        levels,
    )

    assert patches is None


def test_stepped_flat_patches_reject_unsupported_strip_fill_extent():
    footprint = Polygon([(0, 0), (20, 0), (20, 10), (0, 10)])
    low_points = [
        [x, y, 10.0]
        for x in np.linspace(1, 2, 3)
        for y in np.linspace(1, 2, 3)
    ]
    high_points = [
        [x, y, 13.0]
        for x in np.linspace(11, 19, 5)
        for y in np.linspace(1, 9, 5)
    ]
    points = np.asarray([*low_points, *high_points], dtype=float)
    low_indices = np.arange(0, len(low_points))
    high_indices = np.arange(len(low_points), len(points))
    levels = [
        lod2_module.SteppedFlatLevel(
            lod2_module.RoofPlane(0.0, 0.0, 10.0, low_indices),
            low_indices,
        ),
        lod2_module.SteppedFlatLevel(
            lod2_module.RoofPlane(0.0, 0.0, 13.0, high_indices),
            high_indices,
        ),
    ]

    patches, reason = lod2_module._stepped_flat_level_patches(points, footprint, levels)

    assert patches is None
    assert reason == lod2_module.STEPPED_FLAT_PATCH_FAILED


def test_stepped_flat_region_patches_cover_l_limb_levels():
    footprint = Polygon([(0, 0), (10, 0), (10, 4), (4, 4), (4, 10), (0, 10)])
    points = np.asarray(_l_limb_stepped_flat_points(), dtype=float)
    levels = [
        lod2_module.SteppedFlatLevel(
            lod2_module.RoofPlane(0.0, 0.0, 10.0, np.arange(0, 36)),
            np.arange(0, 36),
        ),
        lod2_module.SteppedFlatLevel(
            lod2_module.RoofPlane(0.0, 0.0, 13.0, np.arange(36, len(points))),
            np.arange(36, len(points)),
        ),
    ]

    patches, reason = lod2_module._stepped_flat_level_patches(points, footprint, levels)

    assert reason is None
    assert patches is not None
    assert len(patches) == 2
    assert lod2_module._patches_cover_footprint(footprint, patches)


def test_stepped_flat_patches_reject_unstriped_levels():
    footprint = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
    points = np.asarray(_mixed_height_unstriped_points(), dtype=float)
    levels, reason = lod2_module._stepped_flat_height_levels(
        points,
        [
            lod2_module.RoofPlane(0.0, 0.0, 10.0, np.arange(50)),
            lod2_module.RoofPlane(0.0, 0.0, 13.0, np.arange(50, 100)),
        ],
    )

    patches, reason = lod2_module._stepped_flat_level_patches(points, footprint, levels)

    assert patches is None
    assert reason == lod2_module.STEPPED_FLAT_PATCH_FAILED


def test_stepped_flat_patches_reject_partly_mixed_strip_support():
    footprint = Polygon([(0, 0), (20, 0), (20, 10), (0, 10)])
    low_points = [
        [x, y, 10.0]
        for x in np.linspace(1, 9, 5)
        for y in np.linspace(1, 9, 5)
    ]
    high_points = [
        [x, y, 13.0]
        for x in np.linspace(11, 19, 5)
        for y in np.linspace(1, 9, 5)
    ]
    mixed_high_points = [
        [x, y, 13.0]
        for x in np.linspace(2, 6, 5)
        for y in np.linspace(2, 8, 2)
    ]
    points = np.asarray([*low_points, *high_points, *mixed_high_points], dtype=float)
    levels = [
        lod2_module.SteppedFlatLevel(
            lod2_module.RoofPlane(0.0, 0.0, 10.0, np.arange(0, 25)),
            np.arange(0, 25),
        ),
        lod2_module.SteppedFlatLevel(
            lod2_module.RoofPlane(0.0, 0.0, 13.0, np.arange(25, 60)),
            np.arange(25, 60),
        ),
    ]

    patches, reason = lod2_module._stepped_flat_level_patches(points, footprint, levels)

    assert patches is None
    assert reason == lod2_module.STEPPED_FLAT_PATCH_FAILED


def test_build_stepped_flat_shell_from_patches_is_watertight():
    footprint = Polygon([(0, 0), (30, 0), (30, 10), (0, 10)])
    points = np.asarray(_stepped_flat_roof_points(), dtype=float)
    levels, reason = lod2_module._stepped_flat_height_levels(
        points,
        lod2_module._ransac_planes(points),
    )
    patches, reason = lod2_module._stepped_flat_level_patches(points, footprint, levels)

    shell, reason = lod2_module._build_stepped_flat_shell(footprint, patches, levels, 0.0)

    assert reason == lod2_module.STEPPED_FLAT_SUCCESS
    assert shell is not None
    assert is_watertight(shell)


def test_build_stepped_flat_shell_nodes_partial_shared_boundaries():
    footprint = Polygon([(0, 0), (20, 0), (20, 20), (0, 20)])
    patches = [
        Polygon([(0, 0), (10, 0), (10, 20), (0, 20)]),
        Polygon([(10, 0), (20, 0), (20, 10), (10, 10)]),
        Polygon([(10, 10), (20, 10), (20, 20), (10, 20)]),
    ]
    levels = [
        lod2_module.SteppedFlatLevel(
            lod2_module.RoofPlane(0.0, 0.0, 10.0, np.array([], dtype=int)),
            np.array([], dtype=int),
        ),
        lod2_module.SteppedFlatLevel(
            lod2_module.RoofPlane(0.0, 0.0, 12.0, np.array([], dtype=int)),
            np.array([], dtype=int),
        ),
        lod2_module.SteppedFlatLevel(
            lod2_module.RoofPlane(0.0, 0.0, 12.0, np.array([], dtype=int)),
            np.array([], dtype=int),
        ),
    ]

    shell, reason = lod2_module._build_stepped_flat_shell(footprint, patches, levels, 0.0)

    assert reason == lod2_module.STEPPED_FLAT_SUCCESS
    assert shell is not None
    assert is_watertight(shell)


def test_build_stepped_flat_shell_rejects_three_height_t_junction():
    footprint = Polygon([(0, 0), (20, 0), (20, 20), (0, 20)])
    patches = [
        Polygon([(0, 0), (10, 0), (10, 20), (0, 20)]),
        Polygon([(10, 0), (20, 0), (20, 10), (10, 10)]),
        Polygon([(10, 10), (20, 10), (20, 20), (10, 20)]),
    ]
    levels = [
        lod2_module.SteppedFlatLevel(
            lod2_module.RoofPlane(0.0, 0.0, 10.0, np.array([], dtype=int)),
            np.array([], dtype=int),
        ),
        lod2_module.SteppedFlatLevel(
            lod2_module.RoofPlane(0.0, 0.0, 12.0, np.array([], dtype=int)),
            np.array([], dtype=int),
        ),
        lod2_module.SteppedFlatLevel(
            lod2_module.RoofPlane(0.0, 0.0, 14.0, np.array([], dtype=int)),
            np.array([], dtype=int),
        ),
    ]

    shell, reason = lod2_module._build_stepped_flat_shell(footprint, patches, levels, 0.0)

    assert shell is None
    assert reason == lod2_module.STEPPED_FLAT_STEP_WALL_FAILED


def test_build_stepped_flat_shell_rejects_non_simple_surfaces(monkeypatch):
    footprint = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
    patches = [
        Polygon([(0, 0), (5, 0), (5, 10), (0, 10)]),
        Polygon([(5, 0), (10, 0), (10, 10), (5, 10)]),
    ]
    levels = [
        lod2_module.SteppedFlatLevel(
            lod2_module.RoofPlane(0.0, 0.0, 10.0, np.array([], dtype=int)),
            np.array([], dtype=int),
        ),
        lod2_module.SteppedFlatLevel(
            lod2_module.RoofPlane(0.0, 0.0, 12.0, np.array([], dtype=int)),
            np.array([], dtype=int),
        ),
    ]
    bowtie_wall = _surface([[0, 0, 0], [1, 1, 0], [0, 1, 0], [1, 0, 0]])
    monkeypatch.setattr(
        lod2_module,
        "_wall_surfaces",
        lambda candidate_footprint, roof_surfaces, ground_height: [bowtie_wall],
    )
    monkeypatch.setattr(lod2_module, "is_watertight", lambda shell: True)

    shell, reason = lod2_module._build_stepped_flat_shell(footprint, patches, levels, 0.0)

    assert shell is None
    assert reason == lod2_module.STEPPED_FLAT_WATERTIGHT_FAILED


def _building_with_footprint(points, *, ground_height=0.0):
    building = Building()
    footprint = Surface()
    footprint.from_polygon(Polygon([(0, 0), (10, 0), (10, 10), (0, 10)]), 10.0)
    building.add_geometry(footprint, GeometryType.LOD0)
    building.add_geometry(PointCloud(points=np.array(points, dtype=float)), GeometryType.POINT_CLOUD)
    building.attributes["ground_height"] = ground_height
    building.attributes["height"] = 10.0 - ground_height
    return building


def _building_with_polygon(points, polygon, *, ground_height=0.0):
    building = Building()
    footprint = Surface()
    footprint.from_polygon(polygon, 10.0)
    building.add_geometry(footprint, GeometryType.LOD0)
    building.add_geometry(PointCloud(points=np.array(points, dtype=float)), GeometryType.POINT_CLOUD)
    building.attributes["ground_height"] = ground_height
    building.attributes["height"] = 10.0 - ground_height
    return building


def _flat_roof_points(z=10.0):
    return [[x, y, z] for x in np.linspace(1, 9, 5) for y in np.linspace(1, 9, 5)]


def _region_points(count):
    return np.asarray(
        [[1.0 + 0.01 * (index % 10), 1.0 + 0.01 * (index // 10), 10.0] for index in range(count)],
        dtype=float,
    )


def _planes_with_inlier_counts(*counts):
    start = 0
    planes = []
    for index, count in enumerate(counts):
        planes.append(lod2_module.RoofPlane(0.1 * index, 0.0, 10.0 + index, np.arange(start, start + count)))
        start += count
    return planes


def _decomposition_summary_values(messages):
    summary = next(message for message in messages if message.startswith("LOD2 decomposition summary:"))
    return {
        part.split("=")[0]: int(part.split("=")[1])
        for part in summary.removeprefix("LOD2 decomposition summary: ").split()
    }


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
    assert not any(message.startswith("LOD2 template gate summary:") for message in messages)
    assert not any(message.startswith("LOD2 decomposition summary:") for message in messages)


def test_projected_roof_surfaces_must_cover_footprint():
    footprint = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
    half_roof = _surface([[0, 0, 10], [5, 0, 10], [5, 10, 10], [0, 10, 10]])
    other_half_roof = _surface([[5, 0, 10], [10, 0, 10], [10, 10, 10], [5, 10, 10]])

    assert not _roof_surfaces_cover_footprint(footprint, [half_roof])
    assert _roof_surfaces_cover_footprint(footprint, [half_roof, other_half_roof])


def test_concave_vertex_indices_classifies_basic_footprints():
    rectangle = Polygon([(0, 0), (10, 0), (10, 4), (0, 4)])
    l_shape = Polygon([(0, 0), (10, 0), (10, 4), (4, 4), (4, 10), (0, 10)])
    u_shape = Polygon([(0, 0), (10, 0), (10, 10), (7, 10), (7, 3), (3, 3), (3, 10), (0, 10)])
    plus_like = Polygon([(3, 0), (7, 0), (7, 3), (10, 3), (10, 7), (7, 7), (7, 10), (3, 10), (3, 7), (0, 7), (0, 3), (3, 3)])

    assert len(_concave_vertex_indices(rectangle)) == 0
    assert len(_concave_vertex_indices(l_shape)) == 1
    assert len(_concave_vertex_indices(u_shape)) == 2
    assert len(_concave_vertex_indices(plus_like)) == 4


def test_decomposition_shape_reason_uses_simplified_footprint():
    l_shape = Polygon([(0, 0), (10, 0), (10, 4), (4, 4), (4.1, 5), (4, 6), (4, 10), (0, 10)])
    u_shape = Polygon([(0, 0), (10, 0), (10, 10), (7, 10), (7, 3), (3, 3), (3, 10), (0, 10)])

    assert _decomposition_shape_reason(l_shape) == lod2_module.DECOMPOSITION_L_LIKE
    assert _decomposition_shape_reason(u_shape) == lod2_module.DECOMPOSITION_T_OR_U_LIKE


def test_decompose_footprint_splits_l_shape_into_two_regions():
    l_shape = Polygon([(0, 0), (10, 0), (10, 4), (4, 4), (4, 10), (0, 10)])

    decomposition = _decompose_footprint(l_shape)

    assert decomposition is not None
    assert decomposition.family_reason == lod2_module.DECOMPOSITION_L_LIKE
    assert len(decomposition.pieces) == 2
    assert len(decomposition.slice_lines) == 1
    assert lod2_module._patches_cover_footprint(l_shape, decomposition.pieces)


def test_decompose_footprint_splits_u_shape_into_three_regions():
    u_shape = Polygon([(0, 0), (10, 0), (10, 10), (7, 10), (7, 3), (3, 3), (3, 10), (0, 10)])

    decomposition = _decompose_footprint(u_shape)

    assert decomposition is not None
    assert decomposition.family_reason == lod2_module.DECOMPOSITION_T_OR_U_LIKE
    assert len(decomposition.pieces) == 3
    assert len(decomposition.slice_lines) == 2
    assert lod2_module._patches_cover_footprint(u_shape, decomposition.pieces)


def test_decompose_footprint_rejects_complex_shape():
    plus_like = Polygon([(3, 0), (7, 0), (7, 3), (10, 3), (10, 7), (7, 7), (7, 10), (3, 10), (3, 7), (0, 7), (0, 3), (3, 3)])

    assert _decompose_footprint(plus_like) is None


def test_decompose_footprint_rejects_axis_misaligned_concavity():
    skewed = Polygon([(0, 0), (10, 0), (10, 4), (4, 7), (4, 10), (0, 10)])

    assert _decompose_footprint(skewed) is None


def test_decompose_footprint_tiebreaking_is_deterministic():
    u_shape = Polygon([(0, 0), (12, 0), (12, 10), (8, 10), (8, 4), (4, 4), (4, 10), (0, 10)])

    first = _decompose_footprint(u_shape)
    second = _decompose_footprint(u_shape)

    assert first is not None
    assert second is not None
    assert [round(piece.area, 6) for piece in first.pieces] == [round(piece.area, 6) for piece in second.pieces]
    assert [tuple(np.round(line.bounds, 6)) for line in first.slice_lines] == [tuple(np.round(line.bounds, 6)) for line in second.slice_lines]


def test_decomposition_region_roof_surfaces_use_local_ransac():
    left = [[x, y, 10.0] for x in np.linspace(1, 3, 5) for y in np.linspace(1, 9, 5)]
    right = [[x, y, 14.0] for x in np.linspace(5, 9, 5) for y in np.linspace(1, 3, 5)]
    points = np.asarray([*left, *right], dtype=float)
    region = Polygon([(0, 0), (4, 0), (4, 10), (0, 10)])

    roof_surfaces, reason = lod2_module._region_roof_surfaces(points, region)

    assert reason is None
    assert roof_surfaces is not None
    assert len(roof_surfaces) == 1
    assert np.allclose(roof_surfaces[0].vertices[:, 2], 10.0)


@pytest.mark.parametrize("counts", [(200, 50, 20), (100, 90, 80)])
def test_region_roof_surfaces_caps_oversegmented_planes(monkeypatch, counts):
    region = Polygon([(0, 0), (4, 0), (4, 4), (0, 4)])
    points = _region_points(sum(counts))
    planes = _planes_with_inlier_counts(*counts)
    captured = []
    monkeypatch.setattr(lod2_module, "_ransac_planes", lambda region_points: planes)

    def roof_surfaces_for_planes(region_points, footprint, kept_planes):
        captured.append(kept_planes)
        return [_surface([[0, 0, 10], [4, 0, 10], [4, 4, 10], [0, 4, 10]])], None

    monkeypatch.setattr(lod2_module, "_roof_surfaces_for_planes", roof_surfaces_for_planes)

    roof_surfaces, reason = lod2_module._region_roof_surfaces(points, region)

    assert reason is None
    assert roof_surfaces is not None
    assert captured == [planes[:2]]


def test_region_roof_surfaces_rejects_balanced_oversegmentation(monkeypatch):
    region = Polygon([(0, 0), (4, 0), (4, 4), (0, 4)])
    points = _region_points(200)
    planes = _planes_with_inlier_counts(50, 50, 50, 50)
    decomposition_counts = Counter()
    monkeypatch.setattr(lod2_module, "_ransac_planes", lambda region_points: planes)

    roof_surfaces, reason = lod2_module._region_roof_surfaces(
        points,
        region,
        decomposition_counts,
    )

    assert roof_surfaces is None
    assert reason == lod2_module.REGION_ROOF_DROPPED_TOO_MANY
    assert decomposition_counts[lod2_module.REGION_ROOF_DROPPED_TOO_MANY] == 1


def _gable_roof_points():
    points = []
    for x in np.linspace(1, 4.5, 5):
        for y in np.linspace(1, 9, 5):
            points.append([x, y, 10.0 + 0.5 * x])
    for x in np.linspace(5.5, 9, 5):
        for y in np.linspace(1, 9, 5):
            points.append([x, y, 15.0 - 0.5 * x])
    return points


def _l_shape_flat_points():
    left = [[x, y, 10.0] for x in np.linspace(1, 3, 5) for y in np.linspace(1, 9, 5)]
    bottom = [[x, y, 12.0] for x in np.linspace(5, 9, 5) for y in np.linspace(1, 3, 5)]
    return [*left, *bottom]


def _u_shape_flat_points():
    left = [[x, y, 10.0] for x in np.linspace(1, 2.5, 5) for y in np.linspace(3.5, 9, 5)]
    bottom = [[x, y, 12.0] for x in np.linspace(4, 6, 5) for y in np.linspace(1, 2.5, 5)]
    right = [[x, y, 14.0] for x in np.linspace(7.5, 9, 5) for y in np.linspace(3.5, 9, 5)]
    return [*left, *bottom, *right]


def _l_limb_stepped_flat_points():
    low = [
        [x, y, 10.0]
        for x in np.linspace(1, 3, 4)
        for y in np.linspace(1, 9, 9)
    ]
    high = [
        [x, y, 13.0]
        for x in np.linspace(4.5, 9, 7)
        for y in np.linspace(1, 3, 4)
    ]
    return [*low, *high]


def _stepped_flat_roof_points():
    points = []
    for x0, x1, z in [(1, 9, 10.0), (11, 19, 13.0), (21, 29, 11.0)]:
        for x in np.linspace(x0, x1, 5):
            for y in np.linspace(1, 9, 5):
                points.append([x, y, z])
    return points


def _stepped_flat_noisy_roof_points():
    points = _stepped_flat_roof_points()
    points.extend(
        [[15.0 + 0.1 * index, 5.0, 18.0 + 0.02 * index] for index in range(6)]
    )
    return points


def _corner_stepped_flat_roof_fixture():
    footprint = Polygon([(0, 0), (20, 0), (20, 12), (0, 12)])
    low_points = [
        [x, y, 10.0]
        for x in np.linspace(1, 6, 5)
        for y in np.linspace(1, 4, 5)
    ]
    high_upper = [
        [x, y, 13.0]
        for x in np.linspace(1, 19, 7)
        for y in np.linspace(6, 11, 4)
    ]
    high_right = [
        [x, y, 13.0]
        for x in np.linspace(9, 19, 5)
        for y in np.linspace(1, 4, 4)
    ]
    points = np.asarray([*low_points, *high_upper, *high_right], dtype=float)
    low_indices = np.arange(0, len(low_points))
    high_indices = np.arange(len(low_points), len(points))
    return points, footprint, low_indices, high_indices


def _scattered_blob_stepped_flat_roof_fixture():
    footprint = Polygon([(0, 0), (20, 0), (20, 12), (0, 12)])
    high_points = [
        [x, y, 13.0]
        for x in np.linspace(1, 19, 7)
        for y in np.linspace(1, 11, 6)
    ]
    low_points = [
        [cx + dx, cy + dy, 10.0]
        for cx, cy in [(3.0, 3.0), (17.0, 3.0), (4.0, 9.0), (16.0, 10.0)]
        for dx in (-0.2, 0.2, 0.0)
        for dy in (-0.2, 0.2)
    ]
    points = np.asarray([*high_points, *low_points], dtype=float)
    high_indices = np.arange(0, len(high_points))
    low_indices = np.arange(len(high_points), len(points))
    return points, footprint, low_indices, high_indices


def _dominant_flat_with_minor_noise_points():
    flat = [
        [x, y, 10.0]
        for x in np.linspace(1, 29, 10)
        for y in np.linspace(1, 9, 10)
    ]
    lower_edge = [
        [x, y, 8.5]
        for x in np.linspace(1, 29, 4)
        for y in np.linspace(10.5, 11.5, 2)
    ]
    return [*flat, *lower_edge]


def _mixed_height_unstriped_points():
    points = []
    for x in np.linspace(1, 9, 10):
        for y in np.linspace(1, 9, 10):
            z = 10.0 if int(x + y) % 2 == 0 else 13.0
            points.append([x, y, z])
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


def test_candidate_from_parts_preserves_flat_and_gable_paths():
    flat_building = _building_with_footprint(_flat_roof_points())
    gable_building = _building_with_footprint(_gable_roof_points())
    flat_candidate = lod2_module._candidate_lod2_from_parts(
        flat_building.lod0.to_polygon(),
        np.asarray(flat_building.point_cloud.points, dtype=float),
        0.0,
    )
    gable_candidate = lod2_module._candidate_lod2_from_parts(
        gable_building.lod0.to_polygon(),
        np.asarray(gable_building.point_cloud.points, dtype=float),
        0.0,
    )

    assert flat_candidate is not None
    assert gable_candidate is not None
    assert is_watertight(flat_candidate)
    assert is_watertight(gable_candidate)


def _horizontal_roof_heights(multisurface):
    heights = []
    for surface in multisurface.surfaces:
        z_values = np.asarray(surface.vertices, dtype=float)[:, 2]
        if len(z_values) >= 3 and np.allclose(z_values, z_values[0]):
            if z_values[0] > 0.0:
                heights.append(round(float(z_values[0]), 6))
    return sorted(set(heights))


def _horizontal_surface_at_height_covers(multisurface, height, xy):
    point = Point(xy)
    for surface in multisurface.surfaces:
        z_values = np.asarray(surface.vertices, dtype=float)[:, 2]
        if len(z_values) < 3 or not np.allclose(z_values, height):
            continue
        polygon = Polygon(np.asarray(surface.vertices, dtype=float)[:, :2])
        if polygon.is_valid and polygon.covers(point):
            return True
    return False


def test_build_lod2_buildings_creates_watertight_stepped_flat_roof():
    building = _building_with_polygon(
        _stepped_flat_roof_points(),
        Polygon([(0, 0), (30, 0), (30, 10), (0, 10)]),
    )

    result = build_lod2_buildings([building], build_lod1_fallback=False)

    assert result[0].lod2 is not None
    assert is_watertight(result[0].lod2)
    assert _horizontal_roof_heights(result[0].lod2) == [10.0, 11.0, 13.0]


def test_build_lod2_buildings_creates_corner_region_stepped_flat_roof(monkeypatch):
    points, footprint, low_indices, high_indices = _corner_stepped_flat_roof_fixture()
    planes = [
        lod2_module.RoofPlane(0.0, 0.0, 10.0, low_indices),
        lod2_module.RoofPlane(0.0, 0.0, 13.0, high_indices),
    ]
    monkeypatch.setattr(lod2_module, "_ransac_planes", lambda roof_points: planes)
    building = _building_with_polygon(points, footprint)

    result = build_lod2_buildings([building], build_lod1_fallback=False)

    assert result[0].lod2 is not None
    assert is_watertight(result[0].lod2)
    assert _horizontal_roof_heights(result[0].lod2) == [10.0, 13.0]
    assert _horizontal_surface_at_height_covers(result[0].lod2, 10.0, (3.5, 2.5))
    assert _horizontal_surface_at_height_covers(result[0].lod2, 13.0, (4.0, 9.0))
    assert _horizontal_surface_at_height_covers(result[0].lod2, 13.0, (16.5, 3.0))


def test_build_lod2_buildings_rejects_scattered_region_stepped_flat_roof(monkeypatch):
    points, footprint, low_indices, high_indices = _scattered_blob_stepped_flat_roof_fixture()
    planes = [
        lod2_module.RoofPlane(0.0, 0.0, 10.0, low_indices),
        lod2_module.RoofPlane(0.0, 0.0, 13.0, high_indices),
    ]
    monkeypatch.setattr(lod2_module, "_ransac_planes", lambda roof_points: planes)
    building = _building_with_polygon(points, footprint)

    result = build_lod2_buildings([building], build_lod1_fallback=False)

    assert result[0].lod2 is None


def test_build_lod2_buildings_drops_minor_noise_for_stepped_flat_roof():
    building = _building_with_polygon(
        _stepped_flat_noisy_roof_points(),
        Polygon([(0, 0), (30, 0), (30, 10), (0, 10)]),
    )

    result = build_lod2_buildings([building], build_lod1_fallback=False)

    assert result[0].lod2 is not None
    assert is_watertight(result[0].lod2)
    assert _horizontal_roof_heights(result[0].lod2) == [10.0, 11.0, 13.0]


def test_build_lod2_buildings_flat_collapses_dominant_level_with_noise(monkeypatch):
    planes = [
        lod2_module.RoofPlane(0.0, 0.0, 10.0, np.arange(0, 60)),
        lod2_module.RoofPlane(0.0, 0.0, 10.0, np.arange(60, 100)),
        lod2_module.RoofPlane(0.0, 0.0, 8.5, np.arange(100, 108)),
    ]
    monkeypatch.setattr(lod2_module, "_ransac_planes", lambda points: planes)
    building = _building_with_polygon(
        _dominant_flat_with_minor_noise_points(),
        Polygon([(0, 0), (30, 0), (30, 12), (0, 12)]),
    )

    result = build_lod2_buildings([building], build_lod1_fallback=False)

    assert result[0].lod2 is not None
    assert is_watertight(result[0].lod2)
    assert _horizontal_roof_heights(result[0].lod2) == [10.0]


def test_build_lod2_buildings_does_not_flat_collapse_real_step(monkeypatch):
    planes = [
        lod2_module.RoofPlane(0.0, 0.0, 10.0, np.arange(0, 80)),
        lod2_module.RoofPlane(0.0, 0.0, 13.0, np.arange(80, 160)),
        lod2_module.RoofPlane(0.0, 0.0, 10.0, np.arange(0, 20)),
    ]
    monkeypatch.setattr(lod2_module, "_ransac_planes", lambda points: planes)
    building = _building_with_polygon(
        [
            *[
                [x, y, 10.0]
                for x in np.linspace(1, 9, 8)
                for y in np.linspace(1, 9, 10)
            ],
            *[
                [x, y, 13.0]
                for x in np.linspace(11, 19, 8)
                for y in np.linspace(1, 9, 10)
            ],
        ],
        Polygon([(0, 0), (20, 0), (20, 10), (0, 10)]),
    )

    result = build_lod2_buildings([building], build_lod1_fallback=False)

    assert result[0].lod2 is not None
    assert is_watertight(result[0].lod2)
    assert _horizontal_roof_heights(result[0].lod2) == [10.0, 13.0]


def test_build_lod2_buildings_rejects_unstriped_stepped_flat_points():
    building = _building_with_polygon(
        _mixed_height_unstriped_points(),
        Polygon([(0, 0), (10, 0), (10, 10), (0, 10)]),
    )

    result = build_lod2_buildings([building], build_lod1_fallback=False)

    assert result[0].lod2 is None


def test_build_lod2_buildings_decomposes_l_shape_flat_regions(monkeypatch):
    l_shape = Polygon([(0, 0), (10, 0), (10, 4), (4, 4), (4, 10), (0, 10)])
    building = _building_with_polygon(_l_shape_flat_points(), l_shape)
    monkeypatch.setattr(
        lod2_module,
        "_ransac_planes",
        lambda points: [
            lod2_module.RoofPlane(0.0, 0.0, 10.0, np.arange(24)),
            lod2_module.RoofPlane(0.0, 0.0, 12.0, np.arange(24, len(points))),
            lod2_module.RoofPlane(0.2, 0.0, 11.0, np.arange(20)),
        ] if len(points) > 30 else [lod2_module._fit_plane(points)],
    )

    result = build_lod2_buildings([building], build_lod1_fallback=False)

    assert result[0].lod2 is not None
    assert is_watertight(result[0].lod2)


def test_build_lod2_buildings_decomposes_u_shape_flat_regions(monkeypatch):
    u_shape = Polygon([(0, 0), (10, 0), (10, 10), (7, 10), (7, 3), (3, 3), (3, 10), (0, 10)])
    building = _building_with_polygon(_u_shape_flat_points(), u_shape)
    monkeypatch.setattr(
        lod2_module,
        "_ransac_planes",
        lambda points: [
            lod2_module.RoofPlane(0.0, 0.0, 10.0, np.arange(25)),
            lod2_module.RoofPlane(0.0, 0.0, 12.0, np.arange(25, 50)),
            lod2_module.RoofPlane(0.0, 0.0, 14.0, np.arange(50, len(points))),
        ] if len(points) > 50 else [lod2_module._fit_plane(points)],
    )

    result = build_lod2_buildings([building], build_lod1_fallback=False)

    assert result[0].lod2 is not None
    assert is_watertight(result[0].lod2)


def test_split_failed_irregular_footprint_attempts_decomposition(monkeypatch):
    messages = []
    monkeypatch.setattr(lod2_module, "info", messages.append, raising=False)
    l_shape = Polygon([(0, 0), (10, 0), (10, 4), (4, 4), (4, 10), (0, 10)])
    building = _building_with_polygon(_l_shape_flat_points(), l_shape)
    two_planes = [
        lod2_module.RoofPlane(0.2, 0.0, 10.0, np.arange(25)),
        lod2_module.RoofPlane(0.2, 0.0, 12.0, np.arange(25, 50)),
    ]
    monkeypatch.setattr(
        lod2_module,
        "_ransac_planes",
        lambda points: two_planes if len(points) > 30 else [lod2_module._fit_plane(points)],
    )

    build_lod2_buildings([building], build_lod1_fallback=False, log_rejections=True)

    assert building.lod2 is not None
    assert any("decomposition_candidate=1" in message for message in messages)
    assert any("decomposition_success=1" in message for message in messages)


def test_irregular_flat_levels_use_decomposition_not_stepped_flat(monkeypatch):
    messages = []
    monkeypatch.setattr(lod2_module, "info", messages.append, raising=False)
    l_shape = Polygon([(0, 0), (10, 0), (10, 4), (4, 4), (4, 10), (0, 10)])
    building = _building_with_polygon(_l_shape_flat_points(), l_shape)

    build_lod2_buildings([building], build_lod1_fallback=False, log_rejections=True)

    assert building.lod2 is not None
    assert is_watertight(building.lod2)
    assert any("decomposition_success=1" in message for message in messages)
    assert not any("stepped_flat_success" in message for message in messages)


def test_decomposition_counter_invariant(monkeypatch):
    messages = []
    monkeypatch.setattr(lod2_module, "info", messages.append, raising=False)
    unsupported_shape = Polygon([
        (3, 0), (7, 0), (7, 3), (10, 3), (10, 7), (7, 7),
        (7, 10), (3, 10), (3, 7), (0, 7), (0, 3), (3, 3),
    ])
    planes = [
        lod2_module.RoofPlane(0.0, 0.0, 10.0, np.arange(50)),
        lod2_module.RoofPlane(0.2, 0.0, 11.0, np.arange(50, 80)),
        lod2_module.RoofPlane(-0.2, 0.0, 13.0, np.arange(80, 110)),
    ]
    monkeypatch.setattr(lod2_module, "_ransac_planes", lambda points: planes)

    build_lod2_buildings(
        [_building_with_polygon(_flat_roof_points(), unsupported_shape)],
        build_lod1_fallback=False,
        log_rejections=True,
    )

    summary = next(message for message in messages if message.startswith("LOD2 decomposition summary:"))
    values = {
        part.split("=")[0]: int(part.split("=")[1])
        for part in summary.removeprefix("LOD2 decomposition summary: ").split()
    }
    failures = sum(values.get(reason, 0) for reason in lod2_module.DECOMPOSITION_FAILURE_REASONS)
    assert values["decomposition_candidate"] == values.get("decomposition_success", 0) + failures


def test_decomposition_region_rejections_are_counted(monkeypatch):
    messages = []
    monkeypatch.setattr(lod2_module, "info", messages.append, raising=False)
    l_shape = Polygon([(0, 0), (10, 0), (10, 4), (4, 4), (4, 10), (0, 10)])
    building = _building_with_polygon(_l_shape_flat_points(), l_shape)
    planes = [
        lod2_module.RoofPlane(0.0, 0.0, 10.0, np.arange(50)),
        lod2_module.RoofPlane(0.2, 0.0, 11.0, np.arange(50, 80)),
        lod2_module.RoofPlane(-0.2, 0.0, 13.0, np.arange(80, 110)),
    ]
    monkeypatch.setattr(lod2_module, "_ransac_planes", lambda points: planes)

    build_lod2_buildings([building], build_lod1_fallback=False, log_rejections=True)

    summary = next(message for message in messages if message.startswith("LOD2 decomposition summary:"))
    assert "decomposition_candidate=1" in summary
    assert "decomposition_region_roof_failed=1" in summary


def test_decomposition_no_valid_slices_logs_axis_misaligned_subcounter(monkeypatch):
    messages = []
    monkeypatch.setattr(lod2_module, "info", messages.append, raising=False)
    skewed = Polygon([(0, 0), (10, 0), (10, 4), (4, 7), (4, 10), (0, 10)])
    planes = [
        lod2_module.RoofPlane(0.0, 0.0, 10.0, np.arange(50)),
        lod2_module.RoofPlane(0.2, 0.0, 11.0, np.arange(50, 80)),
        lod2_module.RoofPlane(-0.2, 0.0, 13.0, np.arange(80, 110)),
    ]
    monkeypatch.setattr(lod2_module, "_ransac_planes", lambda points: planes)

    build_lod2_buildings(
        [_building_with_polygon(_flat_roof_points(), skewed)],
        build_lod1_fallback=False,
        log_rejections=True,
    )

    values = _decomposition_summary_values(messages)
    assert values["decomposition_no_valid_slices"] == 1
    assert values["no_valid_slices_axis_misaligned"] == 1
    assert values["decomposition_no_valid_slices"] == sum(
        values.get(reason, 0)
        for reason in lod2_module.NO_VALID_SLICE_REASONS
    )


def test_decomposition_axis_misalignment_logs_angle_histogram(monkeypatch):
    messages = []
    monkeypatch.setattr(lod2_module, "info", messages.append, raising=False)
    mostly_aligned = Polygon([(0, 0), (10, 0), (10, 4), (4, 6), (4, 10), (0, 10)])
    strongly_misaligned = Polygon([(0, 0), (10, 0), (10, 4), (4, 7), (4, 10), (0, 10)])
    planes = [
        lod2_module.RoofPlane(0.0, 0.0, 10.0, np.arange(50)),
        lod2_module.RoofPlane(0.2, 0.0, 11.0, np.arange(50, 80)),
        lod2_module.RoofPlane(-0.2, 0.0, 13.0, np.arange(80, 110)),
    ]
    monkeypatch.setattr(lod2_module, "_ransac_planes", lambda points: planes)

    build_lod2_buildings(
        [
            _building_with_polygon(_flat_roof_points(), mostly_aligned),
            _building_with_polygon(_flat_roof_points(), strongly_misaligned),
        ],
        build_lod1_fallback=False,
        log_rejections=True,
    )

    values = _decomposition_summary_values(messages)
    assert values["no_valid_slices_axis_misaligned"] == 2
    assert values["axis_misaligned_angle_15_20"] == 1
    assert values["axis_misaligned_angle_25_30"] == 1
    assert values["no_valid_slices_axis_misaligned"] == sum(
        values.get(reason, 0)
        for reason in lod2_module.AXIS_MISALIGNED_ANGLE_REASONS
    )


def test_decomposition_region_roof_failure_logs_subcounter(monkeypatch):
    messages = []
    monkeypatch.setattr(lod2_module, "info", messages.append, raising=False)
    l_shape = Polygon([(0, 0), (10, 0), (10, 4), (4, 4), (4, 10), (0, 10)])
    building = _building_with_polygon(_l_shape_flat_points(), l_shape)
    planes = [
        lod2_module.RoofPlane(0.0, 0.0, 10.0, np.arange(50)),
        lod2_module.RoofPlane(0.2, 0.0, 11.0, np.arange(50, 80)),
        lod2_module.RoofPlane(-0.2, 0.0, 13.0, np.arange(80, 110)),
    ]
    monkeypatch.setattr(lod2_module, "_ransac_planes", lambda points: planes)

    build_lod2_buildings([building], build_lod1_fallback=False, log_rejections=True)

    values = _decomposition_summary_values(messages)
    assert values["decomposition_region_roof_failed"] == 1
    assert values["region_roof_split_failed"] == 1
    assert values["decomposition_region_roof_failed"] == sum(
        values.get(reason, 0)
        for reason in lod2_module.REGION_ROOF_REASONS
    )


def test_decomposition_region_dropped_too_many_preserves_counter_invariant(monkeypatch):
    messages = []
    monkeypatch.setattr(lod2_module, "info", messages.append, raising=False)
    l_shape = Polygon([(0, 0), (10, 0), (10, 4), (4, 4), (4, 10), (0, 10)])
    planes = _planes_with_inlier_counts(50, 50, 50, 50)
    monkeypatch.setattr(lod2_module, "_ransac_planes", lambda points: planes)

    build_lod2_buildings(
        [_building_with_polygon(_l_shape_flat_points(), l_shape)],
        build_lod1_fallback=False,
        log_rejections=True,
    )

    values = _decomposition_summary_values(messages)
    assert values["decomposition_region_roof_failed"] == 1
    assert values["region_roof_dropped_too_many"] == 1
    assert values["decomposition_region_roof_failed"] == sum(
        values.get(reason, 0)
        for reason in lod2_module.REGION_ROOF_REASONS
    )


def test_decomposed_watertight_failure_detects_unpaired_ridge_on_slice():
    footprint = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
    slice_lines = [LineString([(5, 0), (5, 10)])]
    roof_surfaces = [
        _surface([[0, 0, 10], [5, 0, 10], [5, 10, 10], [0, 10, 10]]),
        _surface([[5, 0, 12], [10, 0, 12], [10, 10, 12], [5, 10, 12]]),
        _surface([[5, 4, 14], [7, 4, 14], [7, 6, 14], [5, 6, 14]]),
    ]

    reason = lod2_module._decomposed_watertight_failure_reason(
        footprint,
        roof_surfaces,
        slice_lines,
        MultiSurface(surfaces=[]),
    )

    assert reason == lod2_module.WATERTIGHT_UNPAIRED_RIDGE_ON_SLICE


def test_decomposed_watertight_failure_detects_too_few_paired_vertices():
    footprint = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
    slice_lines = [LineString([(5, 0), (5, 10)])]
    roof_surfaces = [
        _surface([[0, 0, 10], [5, 0, 10], [5, 10, 10], [0, 10, 10]]),
        _surface([[5, 0, 12], [10, 0, 12], [10, 10, 12], [8, 10, 12]]),
    ]

    reason = lod2_module._decomposed_watertight_failure_reason(
        footprint,
        roof_surfaces,
        slice_lines,
        MultiSurface(surfaces=[]),
    )

    assert reason == lod2_module.WATERTIGHT_TOO_FEW_PAIRED_VERTICES


def test_edge_count_mismatch_reason_classifies_failed_edges():
    footprint = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
    slice_lines = [LineString([(5, 0), (5, 10)])]

    slice_shell = MultiSurface(surfaces=[_surface([[5, 2, 0], [5, 8, 0], [6, 8, 0], [6, 2, 0]])])
    exterior_shell = MultiSurface(surfaces=[_surface([[2, 0, 0], [4, 0, 0], [4, 1, 0], [2, 1, 0]])])
    interior_shell = MultiSurface(surfaces=[_surface([[2, 2, 0], [4, 2, 0], [4, 4, 0], [2, 4, 0]])])
    outside_shell = MultiSurface(surfaces=[_surface([[12, 2, 0], [14, 2, 0], [14, 4, 0], [12, 4, 0]])])
    excess_shell = _closed_box()
    excess_shell.surfaces.append(excess_shell.surfaces[0].copy(geometry_only=True))

    assert (
        lod2_module._edge_count_mismatch_reason(footprint, slice_lines, slice_shell)
        == lod2_module.EDGE_COUNT_UNMATCHED_ON_SLICE
    )
    assert (
        lod2_module._edge_count_mismatch_reason(footprint, slice_lines, exterior_shell)
        == lod2_module.EDGE_COUNT_UNMATCHED_ON_FOOTPRINT_EXTERIOR
    )
    assert (
        lod2_module._edge_count_mismatch_reason(footprint, slice_lines, interior_shell)
        == lod2_module.EDGE_COUNT_UNMATCHED_INTERIOR
    )
    assert (
        lod2_module._edge_count_mismatch_reason(footprint, slice_lines, outside_shell)
        == lod2_module.EDGE_COUNT_OTHER
    )
    assert (
        lod2_module._edge_count_mismatch_reason(footprint, slice_lines, excess_shell)
        == lod2_module.EDGE_COUNT_EXCESS_COUNT
    )


def test_decomposition_edge_count_mismatch_preserves_counter_invariant(monkeypatch):
    footprint = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
    left_region = Polygon([(0, 0), (5, 0), (5, 10), (0, 10)])
    right_region = Polygon([(5, 0), (10, 0), (10, 10), (5, 10)])
    decomposition = lod2_module.FootprintDecomposition(
        pieces=[left_region, right_region],
        slice_lines=[LineString([(5, 0), (5, 10)])],
        family_reason=lod2_module.DECOMPOSITION_L_LIKE,
        concave_vertex_count=1,
    )
    left_roof = _surface([[0, 0, 10], [5, 0, 10], [5, 10, 10], [0, 10, 10]])
    right_roof = _surface([[5, 0, 12], [10, 0, 12], [10, 10, 12], [5, 10, 12]])
    region_surfaces = iter([([left_roof], None), ([right_roof], None)])
    decomposition_counts = Counter()

    monkeypatch.setattr(
        lod2_module,
        "_decompose_footprint",
        lambda candidate_footprint, decomposition_counts=None: decomposition,
    )
    monkeypatch.setattr(
        lod2_module,
        "_region_roof_surfaces",
        lambda points, region, decomposition_counts=None: next(region_surfaces),
    )
    monkeypatch.setattr(lod2_module, "_internal_junction_surfaces", lambda *args: [])

    shell, reason = lod2_module._build_decomposed_shell(
        footprint,
        np.empty((0, 3)),
        0.0,
        decomposition_counts,
    )

    assert shell is None
    assert reason == lod2_module.DECOMPOSITION_WATERTIGHT_FAILED
    assert decomposition_counts[lod2_module.WATERTIGHT_EDGE_COUNT_MISMATCH] == 1
    assert decomposition_counts[lod2_module.EDGE_COUNT_UNMATCHED_ON_SLICE] == 1
    assert decomposition_counts[lod2_module.WATERTIGHT_EDGE_COUNT_MISMATCH] == sum(
        decomposition_counts[subreason]
        for subreason in lod2_module.EDGE_COUNT_MISMATCH_REASONS
    )


def test_decomposition_watertight_failure_preserves_counter_invariant(monkeypatch):
    messages = []
    monkeypatch.setattr(lod2_module, "info", messages.append, raising=False)
    l_shape = Polygon([(0, 0), (10, 0), (10, 4), (4, 4), (4, 10), (0, 10)])
    planes = _planes_with_inlier_counts(50, 40, 30)

    def build_failed_shell(footprint, roof_points, ground_height, decomposition_counts=None):
        decomposition_counts[lod2_module.WATERTIGHT_EDGE_COUNT_MISMATCH] += 1
        return None, lod2_module.DECOMPOSITION_WATERTIGHT_FAILED

    monkeypatch.setattr(lod2_module, "_ransac_planes", lambda points: planes)
    monkeypatch.setattr(lod2_module, "_build_decomposed_shell", build_failed_shell)

    build_lod2_buildings(
        [_building_with_polygon(_l_shape_flat_points(), l_shape)],
        build_lod1_fallback=False,
        log_rejections=True,
    )

    values = _decomposition_summary_values(messages)
    assert values["decomposition_watertight_failed"] == 1
    assert values["watertight_edge_count_mismatch"] == 1
    assert values["decomposition_watertight_failed"] == sum(
        values.get(reason, 0)
        for reason in lod2_module.WATERTIGHT_FAILURE_REASONS
    )


def test_decomposed_shell_reconciles_unpaired_ridge_on_slice(monkeypatch):
    footprint = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
    left_region = Polygon([(0, 0), (5, 0), (5, 10), (0, 10)])
    right_region = Polygon([(5, 0), (10, 0), (10, 10), (5, 10)])
    decomposition = lod2_module.FootprintDecomposition(
        pieces=[left_region, right_region],
        slice_lines=[LineString([(5, 0), (5, 10)])],
        family_reason=lod2_module.DECOMPOSITION_L_LIKE,
        concave_vertex_count=1,
    )
    left_roof = _surface(
        [[0, 0, 10], [5, 0, 10], [5, 5, 14], [5, 10, 10], [0, 10, 10]]
    )
    right_roof = _surface([[5, 0, 8], [10, 0, 8], [10, 10, 8], [5, 10, 8]])
    region_surfaces = iter([([left_roof], None), ([right_roof], None)])

    monkeypatch.setattr(
        lod2_module,
        "_decompose_footprint",
        lambda candidate_footprint, decomposition_counts=None: decomposition,
    )
    monkeypatch.setattr(
        lod2_module,
        "_region_roof_surfaces",
        lambda points, region, decomposition_counts=None: next(region_surfaces),
    )

    shell, reason = lod2_module._build_decomposed_shell(footprint, np.empty((0, 3)), 0.0)

    assert reason == lod2_module.DECOMPOSITION_SUCCESS
    assert shell is not None
    assert is_watertight(shell)


def test_decomposed_shell_reconciles_slice_vertex_when_height_order_changes(monkeypatch):
    footprint = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
    left_region = Polygon([(0, 0), (5, 0), (5, 10), (0, 10)])
    right_region = Polygon([(5, 0), (10, 0), (10, 10), (5, 10)])
    decomposition = lod2_module.FootprintDecomposition(
        pieces=[left_region, right_region],
        slice_lines=[LineString([(5, 0), (5, 10)])],
        family_reason=lod2_module.DECOMPOSITION_L_LIKE,
        concave_vertex_count=1,
    )
    left_roof = _surface(
        [[0, 0, 10], [5, 0, 10], [5, 5, 14], [5, 10, 10], [0, 10, 10]]
    )
    right_roof = _surface([[5, 0, 12], [10, 0, 12], [10, 10, 12], [5, 10, 12]])
    region_surfaces = iter([([left_roof], None), ([right_roof], None)])

    monkeypatch.setattr(
        lod2_module,
        "_decompose_footprint",
        lambda candidate_footprint, decomposition_counts=None: decomposition,
    )
    monkeypatch.setattr(
        lod2_module,
        "_region_roof_surfaces",
        lambda points, region, decomposition_counts=None: next(region_surfaces),
    )

    shell, reason = lod2_module._build_decomposed_shell(footprint, np.empty((0, 3)), 0.0)

    assert reason == lod2_module.DECOMPOSITION_SUCCESS
    assert shell is not None
    assert is_watertight(shell)


def _surface_is_simple(surface):
    vertices = np.asarray(surface.vertices, dtype=float)
    centered = vertices - vertices.mean(axis=0)
    _, _, rotation = np.linalg.svd(centered)
    projected = centered @ rotation[:2].T
    return Polygon(projected).is_valid


def test_decomposed_shell_splits_crossing_junction_walls_into_simple_triangles(monkeypatch):
    footprint = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
    left_region = Polygon([(0, 0), (5, 0), (5, 10), (0, 10)])
    right_region = Polygon([(5, 0), (10, 0), (10, 10), (5, 10)])
    decomposition = lod2_module.FootprintDecomposition(
        pieces=[left_region, right_region],
        slice_lines=[LineString([(5, 0), (5, 10)])],
        family_reason=lod2_module.DECOMPOSITION_L_LIKE,
        concave_vertex_count=1,
    )
    left_roof = _surface(
        [[0, 0, 10], [5, 0, 10], [5, 5, 14], [5, 10, 10], [0, 10, 10]]
    )
    right_roof = _surface([[5, 0, 12], [10, 0, 12], [10, 10, 12], [5, 10, 12]])
    region_surfaces = iter([([left_roof], None), ([right_roof], None)])

    monkeypatch.setattr(
        lod2_module,
        "_decompose_footprint",
        lambda candidate_footprint, decomposition_counts=None: decomposition,
    )
    monkeypatch.setattr(
        lod2_module,
        "_region_roof_surfaces",
        lambda points, region, decomposition_counts=None: next(region_surfaces),
    )

    shell, reason = lod2_module._build_decomposed_shell(footprint, np.empty((0, 3)), 0.0)

    assert reason == lod2_module.DECOMPOSITION_SUCCESS
    assert shell is not None
    assert is_watertight(shell)
    # The left ridge crosses the flat right roof twice along the slice line;
    # each crossing must yield two simple triangles, never a bowtie quad.
    assert all(_surface_is_simple(surface) for surface in shell.surfaces)
    assert sum(len(surface.vertices) == 3 for surface in shell.surfaces) == 4


def test_decomposed_shell_handles_near_coplanar_crossing_regions(monkeypatch):
    footprint = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
    left_region = Polygon([(0, 0), (5, 0), (5, 10), (0, 10)])
    right_region = Polygon([(5, 0), (10, 0), (10, 10), (5, 10)])
    decomposition = lod2_module.FootprintDecomposition(
        pieces=[left_region, right_region],
        slice_lines=[LineString([(5, 0), (5, 10)])],
        family_reason=lod2_module.DECOMPOSITION_L_LIKE,
        concave_vertex_count=1,
    )
    # Nearly coincident planes whose centimetre-scale height gap changes sign
    # along the slice line (the Kungsbacka regression).
    left_roof = _surface(
        [[0, 0, 10.0], [5, 0, 10.0], [5, 10, 10.08], [0, 10, 10.08]]
    )
    right_roof = _surface(
        [[5, 0, 10.03], [10, 0, 10.03], [10, 10, 10.0], [5, 10, 10.0]]
    )
    region_surfaces = iter([([left_roof], None), ([right_roof], None)])

    monkeypatch.setattr(
        lod2_module,
        "_decompose_footprint",
        lambda candidate_footprint, decomposition_counts=None: decomposition,
    )
    monkeypatch.setattr(
        lod2_module,
        "_region_roof_surfaces",
        lambda points, region, decomposition_counts=None: next(region_surfaces),
    )

    shell, reason = lod2_module._build_decomposed_shell(footprint, np.empty((0, 3)), 0.0)

    assert reason == lod2_module.DECOMPOSITION_SUCCESS
    assert shell is not None
    assert is_watertight(shell)
    assert all(_surface_is_simple(surface) for surface in shell.surfaces)
    assert sum(len(surface.vertices) == 3 for surface in shell.surfaces) == 2


def test_decomposed_shell_welds_near_equal_boundary_endpoint_heights(monkeypatch):
    footprint = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
    left_region = Polygon([(0, 0), (5, 0), (5, 10), (0, 10)])
    right_region = Polygon([(5, 0), (10, 0), (10, 10), (5, 10)])
    decomposition = lod2_module.FootprintDecomposition(
        pieces=[left_region, right_region],
        slice_lines=[LineString([(5, -20), (5, 30)])],
        family_reason=lod2_module.DECOMPOSITION_L_LIKE,
        concave_vertex_count=1,
    )
    left_roof = _surface([[0, 0, 10.0], [5, 0, 10.0], [5, 10, 10.0], [0, 10, 10.0]])
    right_roof = _surface(
        [[5, 0, 10.009], [10, 0, 10.009], [10, 10, 10.0], [5, 10, 10.0]]
    )
    region_surfaces = iter([([left_roof], None), ([right_roof], None)])

    monkeypatch.setattr(
        lod2_module,
        "_decompose_footprint",
        lambda candidate_footprint, decomposition_counts=None: decomposition,
    )
    monkeypatch.setattr(
        lod2_module,
        "_region_roof_surfaces",
        lambda points, region, decomposition_counts=None: next(region_surfaces),
    )

    shell, reason = lod2_module._build_decomposed_shell(footprint, np.empty((0, 3)), 0.0)

    assert reason == lod2_module.DECOMPOSITION_SUCCESS
    assert shell is not None
    assert is_watertight(shell)
    boundary_heights = [
        vertex[2]
        for surface in shell.surfaces
        for vertex in surface.vertices
        if vertex[2] > 0.0 and np.linalg.norm(vertex[:2] - np.array([5.0, 0.0])) <= 1e-3
    ]
    assert max(boundary_heights) - min(boundary_heights) <= 1e-3


def test_decomposed_shell_rejects_crossing_too_close_to_region_corner(monkeypatch):
    footprint = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
    quadrants = [
        Polygon([(0, 0), (5, 0), (5, 5), (0, 5)]),
        Polygon([(5, 0), (10, 0), (10, 5), (5, 5)]),
        Polygon([(0, 5), (5, 5), (5, 10), (0, 10)]),
        Polygon([(5, 5), (10, 5), (10, 10), (5, 10)]),
    ]
    decomposition = lod2_module.FootprintDecomposition(
        pieces=quadrants,
        slice_lines=[LineString([(5, 0), (5, 10)]), LineString([(0, 5), (10, 5)])],
        family_reason=lod2_module.DECOMPOSITION_T_OR_U_LIKE,
        concave_vertex_count=2,
    )
    # The SW and NW planes cross on the y=5 slice line 0.75 mm west of the
    # shared corner (5, 5), too close to insert a split vertex. The shell
    # must be rejected; the unrelated SE/NE corner heights (3 m apart from
    # the crossing pair) must never be averaged into the crossing repair.
    sw_roof = _surface([[0, 0, 20.0], [5, 0, 10.0], [5, 5, 10.0], [0, 5, 20.0]])
    se_roof = _surface([[5, 0, 3.0], [10, 0, 3.0], [10, 5, 3.0], [5, 5, 3.0]])
    nw_roof = _surface(
        [[0, 5, 10.0015], [5, 5, 10.0015], [5, 10, 10.0015], [0, 10, 10.0015]]
    )
    ne_roof = _surface([[5, 5, 15.0], [10, 5, 15.0], [10, 10, 15.0], [5, 10, 15.0]])
    region_surfaces = iter(
        [([sw_roof], None), ([se_roof], None), ([nw_roof], None), ([ne_roof], None)]
    )

    monkeypatch.setattr(
        lod2_module,
        "_decompose_footprint",
        lambda candidate_footprint, decomposition_counts=None: decomposition,
    )
    monkeypatch.setattr(
        lod2_module,
        "_region_roof_surfaces",
        lambda points, region, decomposition_counts=None: next(region_surfaces),
    )

    shell, reason = lod2_module._build_decomposed_shell(footprint, np.empty((0, 3)), 0.0)

    assert shell is None
    assert reason == lod2_module.DECOMPOSITION_JUNCTION_FAILED


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
        lod2_module.RoofPlane(0.2, 0.0, 10.0, np.arange(12)),
        lod2_module.RoofPlane(0.2, 0.0, 10.0, np.arange(12, 24)),
    ]
    monkeypatch.setattr(lod2_module, "_ransac_planes", lambda points: planes)
    building = _building_with_footprint(_flat_roof_points())

    build_lod2_buildings([building], build_lod1_fallback=False, log_rejections=True)

    assert any(
        message == "LOD2 rejection summary: split_failed=1"
        for message in messages
    )


def test_stepped_flat_success_is_logged_when_rejections_enabled(monkeypatch):
    messages = []
    monkeypatch.setattr(lod2_module, "info", messages.append, raising=False)
    building = _building_with_polygon(
        _stepped_flat_roof_points(),
        Polygon([(0, 0), (30, 0), (30, 10), (0, 10)]),
    )

    build_lod2_buildings([building], build_lod1_fallback=False, log_rejections=True)

    assert building.lod2 is not None
    assert any(
        message == "LOD2 stepped-flat summary: stepped_flat_candidate=1 stepped_flat_success=1"
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


def test_build_lod2_buildings_logs_template_gate_summary(monkeypatch):
    messages = []
    monkeypatch.setattr(lod2_module, "info", messages.append, raising=False)
    four_planes = [
        lod2_module.RoofPlane(0.30, 0.00, 10.0, np.arange(100)),
        lod2_module.RoofPlane(-0.30, 0.00, 16.0, np.arange(100, 180)),
        lod2_module.RoofPlane(0.00, 0.30, 10.0, np.arange(180, 220)),
        lod2_module.RoofPlane(0.00, -0.30, 16.0, np.arange(220, 250)),
    ]
    three_planes = four_planes[:3]
    monkeypatch.setattr(
        lod2_module,
        "_ransac_planes",
        lambda points: three_planes if len(points) == 30 else four_planes,
    )
    l_shape = Polygon([(0, 0), (10, 0), (10, 4), (4, 4), (4, 10), (0, 10)])
    wide_rect = Polygon([(0, 0), (30, 0), (30, 10), (0, 10)])
    large_square = Polygon([(0, 0), (20, 0), (20, 20), (0, 20)])
    rectangular_few_dominant_points = [*_flat_roof_points(), *[[1, 1, 10]] * 5]
    buildings = [
        _building_with_footprint(_flat_roof_points()),
        _building_with_polygon(_flat_roof_points(), wide_rect),
        _building_with_polygon(_flat_roof_points(), large_square),
        _building_with_polygon(_flat_roof_points(), l_shape),
        _building_with_footprint(rectangular_few_dominant_points),
    ]

    build_lod2_buildings(buildings, build_lod1_fallback=False, log_rejections=True)

    summary_lines = [
        message for message in messages
        if message.startswith("LOD2 template gate summary:")
    ]
    assert len(summary_lines) == 1
    assert "unsupported_rect_4plane=3" in summary_lines[0]
    assert "unsupported_near_square_4plane=1" in summary_lines[0]
    assert "unsupported_irregular_or_other=2" in summary_lines[0]
    assert "unsupported_irregular_footprint=1" in summary_lines[0]
    assert "unsupported_rect_few_dominant=1" in summary_lines[0]


def test_build_lod2_buildings_logs_decomposition_summary_when_requested(monkeypatch):
    messages = []
    monkeypatch.setattr(lod2_module, "info", messages.append, raising=False)
    planes = [
        lod2_module.RoofPlane(0.0, 0.0, 10.0, np.arange(100)),
        lod2_module.RoofPlane(0.2, 0.0, 11.0, np.arange(100, 140)),
        lod2_module.RoofPlane(-0.2, 0.0, 13.0, np.arange(140, 180)),
    ]
    monkeypatch.setattr(lod2_module, "_ransac_planes", lambda points: planes)
    l_shape = Polygon([(0, 0), (10, 0), (10, 4), (4, 4), (4, 10), (0, 10)])
    building = _building_with_polygon(_flat_roof_points(), l_shape)

    build_lod2_buildings([building], build_lod1_fallback=False, log_rejections=True)

    summary_lines = [
        message for message in messages
        if message.startswith("LOD2 decomposition summary:")
    ]
    assert len(summary_lines) == 1
    assert "decomposition_candidate=1" in summary_lines[0]
    assert "decomposition_l_like=1" in summary_lines[0]
    assert "decomposition_sparse_region_points=1" in summary_lines[0]


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
