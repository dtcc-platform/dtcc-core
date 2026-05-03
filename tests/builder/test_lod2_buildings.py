from collections import Counter
from pathlib import Path

import numpy as np
import pytest
from shapely.geometry import LineString, Polygon

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
        lod2_module.RoofPlane(0.0, 0.0, 10.0, np.arange(25)),
        lod2_module.RoofPlane(0.0, 0.0, 12.0, np.arange(25, 50)),
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
