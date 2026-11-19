"""Private prism decomposition: occupancy, shell closure and absolute budgets."""

from copy import deepcopy
from math import pi

import numpy as np
import pytest
import shapely
from shapely.affinity import translate
from shapely.geometry import MultiPolygon, Polygon, box
from shapely.geometry.polygon import orient
from shapely.ops import triangulate

from dtcc_core.builder.geometry_builders.prisms import (
    _clusters,
    _decompose_prisms,
    _node_ring,
    _PrismError,
    _PrismFace,
    _volume_budget,
)


def _parts(polygon):
    return [polygon] if isinstance(polygon, Polygon) else list(polygon.geoms)


def _cap(polygon, z, down=False):
    polygon = orient(polygon, sign=-1 if down else 1)
    return _PrismFace(
        tuple(
            np.column_stack([np.asarray(r.coords)[:-1], np.full(len(r.coords) - 1, z)])
            for r in (polygon.exterior, *polygon.interiors)
        )
    )


def _prism(polygon, ground=0, roof=10):
    top, bottom = _cap(polygon, roof), _cap(polygon, ground)
    faces = [top, _cap(polygon, ground, down=True)]
    for upper, lower in zip(top.rings, bottom.rings, strict=True):
        for i in range(len(upper)):
            j = (i + 1) % len(upper)
            faces.append(
                _PrismFace((np.array([upper[j], upper[i], lower[i], lower[j]]),))
            )
    return faces


def _stack(lower, upper, ground=0, step=6, roof=30):
    """Two horizontal slabs, retaining only their exposed caps."""
    faces = []
    for polygon in _parts(lower):
        faces += _prism(polygon, ground, step)[1:]
    for polygon in _parts(upper):
        solid = _prism(polygon, step, roof)
        faces += [solid[0], *solid[2:]]
    faces += [_cap(p, step) for p in _parts(lower.difference(upper)) if not p.is_empty]
    faces += [
        _cap(p, step, down=True)
        for p in _parts(upper.difference(lower))
        if not p.is_empty
    ]
    return faces


def _triangles(faces, alternate=False):
    result = []
    for face in faces:
        ring = face.rings[0]
        normal = np.cross(ring - ring[0], np.roll(ring, -1, axis=0) - ring[0]).sum(
            axis=0
        )
        omitted = int(np.argmax(abs(normal)))
        axes = [i for i in range(3) if i != omitted]
        polygon = Polygon(ring[:, axes], [h[:, axes] for h in face.rings[1:]])
        if alternate and len(ring) == 4 and len(face.rings) == 1:
            pieces = [ring[[0, 1, 3]], ring[[1, 2, 3]]]
        else:
            pieces = []
            for triangle in triangulate(polygon):
                if not polygon.covers(triangle):
                    continue
                values = np.zeros((3, 3))
                values[:, axes] = np.asarray(triangle.exterior.coords)[:-1]
                values[:, omitted] = (
                    ring[0, omitted]
                    - ((values[:, axes] - ring[0, axes]) @ normal[axes])
                    / normal[omitted]
                )
                pieces.append(values)
        result += [_PrismFace((values,), face.source_index) for values in pieces]
    return result


def _input(faces, kind):
    if kind == "rings":
        return faces
    if kind == "triangles":
        return _triangles(faces)
    if kind == "alternate":
        return _triangles(faces, alternate=True)
    from dtcc_core.model import Mesh, MultiSurface, Surface

    if kind == "surfaces":
        return MultiSurface(
            surfaces=[
                Surface(
                    vertices=f.rings[0].copy(), holes=[r.copy() for r in f.rings[1:]]
                )
                for f in faces
            ]
        )
    triangles = _triangles(faces)
    vertices = np.vstack([f.rings[0] for f in triangles])
    indices = np.arange(len(vertices)).reshape(-1, 3)
    if kind == "mesh":
        vertices, inverse = np.unique(vertices, axis=0, return_inverse=True)
        indices = inverse[indices]
    return Mesh(vertices=vertices, faces=indices)


BASE = box(0, 0, 20, 20)
TOWER = box(7, 7, 13, 13)
COURTYARD = BASE.difference(TOWER)
L_SHAPE = Polygon([(0, 0), (10, 0), (10, 3), (3, 3), (3, 10), (0, 10)])
KINDS = ("rings", "surfaces", "triangles", "alternate", "mesh", "soup")


def _signature(results):
    return [
        (c.ground_z, c.base.wkb_hex, [(z, p.wkb_hex) for z, p in c.levels])
        for c in results
    ]


def _assert_cells(component, base, levels, volume):
    assert component.base.equals(base)
    assert [z for z, _ in component.levels] == sorted(levels)
    polygons = []
    for z, polygon in component.levels:
        assert polygon.is_valid and polygon.area > 0
        assert polygon.equals(levels[z])
        polygons += _parts(polygon)
    assert shapely.union_all(polygons).equals(base)
    assert sum(p.area for p in polygons) == pytest.approx(base.area, abs=1e-10, rel=0)
    for i, polygon in enumerate(polygons):
        for other in polygons[i + 1 :]:
            assert polygon.intersection(other).area == 0
    diagnostics = component.diagnostics
    assert diagnostics["shell_volume"] == pytest.approx(volume, abs=1e-8, rel=0)
    assert diagnostics["source_shell_volume"] == pytest.approx(volume, abs=1e-8, rel=0)
    assert diagnostics["prism_volume"] == pytest.approx(volume, abs=1e-8, rel=0)
    assert diagnostics["source_volume_residual"] <= diagnostics["volume_budget"]
    assert diagnostics["atomic_edges"] > 0
    assert diagnostics["max_displacement"] == 0
    assert diagnostics["boundary_residual"] <= 1e-10
    if hasattr(shapely, "coverage_is_valid") and shapely.geos_version >= (3, 12, 0):
        assert shapely.coverage_is_valid(polygons)


@pytest.mark.parametrize("kind", KINDS)
@pytest.mark.parametrize("polygon", [BASE, L_SHAPE, COURTYARD, box(0, 0, 0.05, 0.05)])
def test_flat_prisms_preserve_holes_concavity_and_small_areas(kind, polygon):
    results = _decompose_prisms(_input(_prism(polygon, 2, 12), kind))
    assert len(results) == 1 and results[0].ground_z == 2
    _assert_cells(results[0], polygon, {12: polygon}, polygon.area * 10)


@pytest.mark.parametrize("kind", KINDS)
@pytest.mark.parametrize(
    "tower",
    [TOWER, box(0, 7, 6, 13), MultiPolygon([box(2, 2, 6, 6), box(14, 14, 18, 18)])],
)
def test_stepped_shells_share_boundaries_and_keep_absolute_roof_heights(kind, tower):
    faces = _stack(BASE, tower)
    results = _decompose_prisms(_input(faces, kind))
    assert len(results) == 1 and results[0].ground_z == 0
    _assert_cells(
        results[0],
        BASE,
        {6: BASE.difference(tower), 30: tower},
        (BASE.area - tower.area) * 6 + tower.area * 30,
    )
    _assert_shared_segmentation(results[0])
    assert _signature(results) == _signature(_decompose_prisms(_input(faces, "rings")))


def _assert_shared_segmentation(component):
    by_level = []
    for _, polygon in component.levels:
        segments = []
        for piece in _parts(polygon):
            for ring in (piece.exterior, *piece.interiors):
                points = np.asarray(ring.coords)
                segments += [
                    shapely.LineString([a, b])
                    for a, b in zip(points[:-1], points[1:], strict=True)
                ]
        by_level.append(segments)
    for i, segments in enumerate(by_level):
        for others in by_level[i + 1 :]:
            for a in segments:
                for b in others:
                    overlap = a.intersection(b)
                    if overlap.length > 0:
                        assert overlap.equals(a) and overlap.equals(b)


@pytest.mark.parametrize("kind", KINDS)
def test_mixed_winding_face_order_and_triangulation_do_not_change_result(kind):
    faces = _stack(BASE, TOWER)
    expected = _signature(_decompose_prisms(_input(faces, kind)))
    rng = np.random.default_rng(42)
    shuffled = []
    for face in faces:
        shuffled.append(
            _PrismFace(
                tuple(
                    r[::-1].copy() if rng.random() < 0.5 else np.roll(r, 1, axis=0)
                    for r in face.rings
                )
            )
        )
    rng.shuffle(shuffled)
    assert _signature(_decompose_prisms(_input(shuffled, kind))) == expected


def test_duplicate_faces_with_different_closures_and_subdivisions_are_counted():
    faces = _prism(COURTYARD)
    duplicate = deepcopy(faces[0])
    duplicate.rings = tuple(
        np.vstack([r[0], (r[0] + r[1]) / 2, r[1:], r[0]]) for r in duplicate.rings
    )
    result = _decompose_prisms([*faces, duplicate, deepcopy(faces[-1])])[0]
    _assert_cells(result, COURTYARD, {10: COURTYARD}, 3640)
    assert result.diagnostics["duplicate_faces"] == 2


@pytest.mark.parametrize("kind", KINDS)
def test_disconnected_components_with_different_grounds_are_ordered(kind):
    low, high = box(0, 0, 3, 3), box(30, 0, 35, 5)
    faces = _prism(high, 5, 15) + _prism(low, 2, 6)
    result = _decompose_prisms(_input(faces, kind))
    assert [c.ground_z for c in result] == [2, 5]
    _assert_cells(result[0], low, {6: low}, 36)
    _assert_cells(result[1], high, {15: high}, 250)


@pytest.mark.parametrize("kind", KINDS)
def test_isolated_box_above_zero_has_inferred_ground(kind):
    result = _decompose_prisms(_input(_prism(BASE, 3, 4), kind))[0]
    assert result.ground_z == 3
    _assert_cells(result, BASE, {4: BASE}, 400)


@pytest.mark.parametrize("kind", KINDS)
def test_closed_shells_touching_on_cap_interiors_are_not_disconnected_components(kind):
    faces = _prism(BASE, 0, 6) + _prism(TOWER, 6, 30)
    with pytest.raises(_PrismError, match="not_2_5d"):
        _decompose_prisms(_input(faces, kind))


@pytest.mark.parametrize("missing", [0, 1, 2, 6])
@pytest.mark.parametrize("kind", KINDS)
def test_missing_roof_base_outside_or_courtyard_wall_is_rejected(kind, missing):
    faces = _prism(COURTYARD)
    del faces[missing]
    with pytest.raises(_PrismError):
        _decompose_prisms(_input(faces, kind))


def test_missing_wall_reports_source_face_and_location():
    faces = _prism(BASE)[:-1]
    for i, face in enumerate(faces):
        face.source_index = 100 + i
    with pytest.raises(_PrismError, match="open_shell") as caught:
        _decompose_prisms(faces)
    assert caught.value.component == 0
    assert caught.value.face_index >= 100
    assert len(caught.value.location) == 3


@pytest.mark.parametrize("kind", KINDS)
@pytest.mark.parametrize(
    "case",
    [
        "internal_cap",
        "shared_wall",
        "point_contact",
        "cavity",
        "cantilever",
        "bridge",
        "stepped_base",
    ],
)
def test_unsupported_occupancy_is_rejected_explicitly(kind, case):
    cases = {
        "internal_cap": _stack(BASE, TOWER) + [_cap(BASE, 6)],
        "shared_wall": _prism(box(0, 0, 10, 10)) + _prism(box(10, 0, 20, 10)),
        "point_contact": _prism(box(0, 0, 10, 10)) + _prism(box(10, 10, 20, 20)),
        "cavity": _prism(BASE) + _prism(TOWER, 2, 8),
        "cantilever": _stack(box(0, 0, 10, 10), box(-2, 0, 12, 10)),
        "bridge": _stack(
            MultiPolygon([box(0, 0, 3, 10), box(7, 0, 10, 10)]), box(0, 0, 10, 10)
        ),
        "stepped_base": _stack(TOWER, BASE),
    }
    with pytest.raises(_PrismError):
        _decompose_prisms(_input(cases[case], kind))


@pytest.mark.parametrize("kind", KINDS)
def test_translated_equal_area_roof_cannot_pass_coverage(kind):
    faces = _prism(BASE)
    faces[0] = _cap(translate(BASE, 1, 0), 10)
    with pytest.raises(_PrismError):
        _decompose_prisms(_input(faces, kind))


@pytest.mark.parametrize("kind", ["rings", "soup"])
def test_tiny_distant_island_cannot_pass_an_area_only_check(kind):
    faces = _prism(BASE) + [_cap(box(100, 100, 100.01, 100.01), 10)]
    with pytest.raises(_PrismError):
        _decompose_prisms(_input(faces, kind))


@pytest.mark.parametrize("kind", KINDS)
@pytest.mark.parametrize("case", ["roof", "taper", "tall_lean"])
def test_slopes_and_displacement_beyond_budget_fail(kind, case):
    faces = _prism(BASE, 0, 1000 if case == "tall_lean" else 10)
    for face in faces:
        for ring in face.rings:
            if case == "roof":
                ring[:, 2] += np.where(ring[:, 2] > 0, ring[:, 0] * 0.01, 0)
            else:
                ring[:, 0] += np.where(ring[:, 2] > 0, 0.1, 0)
    with pytest.raises(_PrismError):
        _decompose_prisms(_input(faces, kind))


@pytest.mark.parametrize("triangulated", [False, True])
@pytest.mark.parametrize("dtype", [np.float64, np.float32])
@pytest.mark.parametrize("amplitude", [0.0002, 0.0005])
def test_noise_is_bounded_and_reported_with_independent_volume(
    dtype, triangulated, amplitude
):
    faces = _stack(BASE, TOWER)
    if triangulated:
        faces = _triangles(faces)
    rng = np.random.default_rng(53)
    for face in faces:
        face.rings = tuple(
            (r + rng.uniform(-amplitude, amplitude, r.shape)).astype(dtype)
            for r in face.rings
        )
    before = deepcopy(faces)
    result = _decompose_prisms(faces)[0]
    diagnostics = result.diagnostics
    assert 0 < diagnostics["max_displacement"] <= 0.002
    assert abs(diagnostics["prism_volume"] - 3264) <= diagnostics["volume_budget"]
    assert diagnostics["source_volume_residual"] <= diagnostics["volume_budget"]
    assert diagnostics["volume_residual"] < 1e-8
    assert max(diagnostics["height_cluster_spans"]) <= 0.002
    for actual, expected in zip(faces, before, strict=True):
        for a, e in zip(actual.rings, expected.rings, strict=True):
            np.testing.assert_array_equal(a, e)
    assert _signature(_decompose_prisms(faces[::-1])) == _signature((result,))


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("dimension", [2, 3])
def test_edge_noding_emits_each_endpoint_once(dtype, dimension):
    ring = np.array(
        [
            [-7.768348693847656, 16.1204891204834, -12.347917556762695],
            [17.943483352661133, -4.046733379364014, 4.249589920043945],
            [3.9026734828948975, -3.9385478496551514, 17.218843460083008],
            [5.956786632537842, 3.5636281967163086, 8.662971496582031],
        ],
        dtype=dtype,
    )[:, :dimension]
    vertices = np.unique(ring, axis=0)
    np.testing.assert_array_equal(_node_ring(ring, vertices, 1e-5), ring)
    midpoint = (ring[0] + ring[1]) / 2
    vertices = np.vstack([vertices, midpoint])
    expected = np.insert(ring, 1, midpoint, axis=0)
    np.testing.assert_array_equal(_node_ring(ring, vertices[::-1], 1e-5), expected)


def test_close_heights_across_rounding_boundary_form_one_bounded_cluster():
    faces = _triangles(_prism(BASE, 0, 10.00005))
    roofs = [f for f in faces if np.all(f.rings[0][:, 2] > 10)]
    for face, z in zip(roofs, [10.00004999, 10.00005001], strict=True):
        face.rings[0][:, 2] = z
    result = _decompose_prisms(faces, tolerance=0.0001)[0]
    assert len(result.levels) == 1
    assert result.levels[0][0] == pytest.approx(10.00005, abs=1e-12)
    assert result.diagnostics["height_cluster_spans"][-1] == pytest.approx(
        0.00000002, abs=1e-12
    )


def test_height_chains_are_not_merged_transitively():
    faces = [_cap(box(i * 2, 0, i * 2 + 2, 2), 10 + i * 0.0015) for i in range(4)]
    # A common base connects the provisional component before level clustering.
    faces += _prism(box(0, 0, 8, 2), 0, 10)[1:]
    with pytest.raises(_PrismError, match="ambiguous_levels"):
        _decompose_prisms(faces)


def test_vertex_chains_are_not_collapsed_transitively():
    with pytest.raises(_PrismError, match="precision_loss"):
        _clusters(np.array([[0, 0, 0], [0.0015, 0, 0], [0.003, 0, 0]]), 0.002)


def test_normalization_failure_location_is_restored_to_source_coordinates():
    faces = _prism(BASE)
    ring = faces[0].rings[0]
    faces[0].rings = (
        np.vstack(
            [ring[0], ring[0] + [0, 0.0015, 0], ring[0] + [0, 0.003, 0], ring[1:]]
        ),
    )
    for face in faces:
        for values in face.rings:
            values += [100000, 200000, 100]
    with pytest.raises(_PrismError, match="precision_loss") as caught:
        _decompose_prisms(faces)
    assert caught.value.location[0] >= 100000
    assert caught.value.location[1] >= 200000
    assert caught.value.location[2] >= 100


def test_nonzero_feature_collapsed_by_normalization_is_not_dropped():
    with pytest.raises(_PrismError, match="precision_loss"):
        _decompose_prisms(_prism(box(0, 0, 0.0005, 1)))


@pytest.mark.parametrize("kind", KINDS)
def test_unequal_subdivisions_and_t_junctions_are_noded(kind):
    faces = _stack(BASE, box(0, 7, 6, 13))
    for face in faces[::2]:
        ring = face.rings[0]
        face.rings = (
            np.vstack([ring[0], (ring[0] + ring[1]) / 2, ring[1:]]),
            *face.rings[1:],
        )
    result = _decompose_prisms(_input(faces, kind))[0]
    _assert_cells(
        result, BASE, {6: BASE.difference(box(0, 7, 6, 13)), 30: box(0, 7, 6, 13)}, 3264
    )
    _assert_shared_segmentation(result)


def test_large_translation_keeps_absolute_physical_tolerances():
    faces = _stack(BASE, TOWER)
    offset = np.array([1_000_000, 6_000_000, 100])
    for face in faces:
        for ring in face.rings:
            ring += offset
    result = _decompose_prisms(faces)[0]
    assert result.ground_z == 100
    _assert_cells(
        result,
        translate(BASE, offset[0], offset[1]),
        {
            106: translate(COURTYARD, offset[0], offset[1]),
            130: translate(TOWER, offset[0], offset[1]),
        },
        3264,
    )
    faces[0].rings[0][0, 0] += 0.1
    with pytest.raises(_PrismError):
        _decompose_prisms(faces)


def test_volume_budget_formula_charges_boundary_and_cap_displacement():
    levels = [(6, COURTYARD), (30, TOWER)]
    d, epsilon = 0.001, 1e-10
    expected = sum(
        (p.length * d + pi * d**2 * ring_count) * (z - 2) + 2 * d * p.area
        for z, p, ring_count in [(6, COURTYARD, 2), (30, TOWER, 1)]
    )
    expected += epsilon * sum(p.area + p.length * (z - 2) for z, p in levels)
    assert _volume_budget(levels, 2, d, epsilon) == expected


def test_source_precision_loss_at_large_float32_coordinates_is_rejected():
    faces = _prism(translate(BASE, 100000, 200000))
    for face in faces:
        face.rings = tuple(r.astype(np.float32) for r in face.rings)
    with pytest.raises(_PrismError, match="precision_loss"):
        _decompose_prisms(faces)


def test_tall_leaning_wall_within_displacement_and_angle_budget_is_normalized():
    faces = _prism(BASE, 0, 1000)
    for face in faces:
        for ring in face.rings:
            ring[:, 0] += np.where(ring[:, 2] > 0, 0.001, 0)
    result = _decompose_prisms(faces)[0]
    assert 0 < result.diagnostics["max_displacement"] < 0.002
    assert result.levels[0][0] == 1000
    assert result.diagnostics["prism_volume"] == pytest.approx(400000, abs=1e-8, rel=0)


def test_noisy_disconnected_components_use_independent_local_volume_references():
    first = _prism(BASE)
    second = _prism(translate(BASE, 100000, 200000), 2, 12)
    rng = np.random.default_rng(17)
    for face in second:
        face.rings = tuple(
            r + rng.uniform(-0.0002, 0.0002, r.shape) for r in face.rings
        )
    result = _decompose_prisms(first + second)
    assert len(result) == 2
    for component in result:
        diagnostics = component.diagnostics
        assert diagnostics["source_volume_residual"] <= diagnostics["volume_budget"]
        assert diagnostics["prism_volume"] == pytest.approx(4000, abs=0.1, rel=0)


def test_derived_polygons_rebuild_closed_prisms_with_independent_mesh_volume():
    from collections import Counter

    from dtcc_core.builder import building_from_footprint

    component = _decompose_prisms(_stack(BASE, TOWER))[0]
    total_volume = 0.0
    for roof, polygon in component.levels:
        for piece in _parts(polygon):
            building = building_from_footprint(
                piece, roof - component.ground_z, ground_height=component.ground_z
            )
            mesh = building.lod1.mesh(weld=True)
            edges = Counter()
            for face in mesh.faces:
                for a, b in zip(face, np.roll(face, -1), strict=True):
                    edges[int(a), int(b)] += 1
            assert all(n == 1 and edges[b, a] == 1 for (a, b), n in edges.items())
            xyz = mesh.vertices[mesh.faces] - mesh.vertices[0]
            total_volume += (
                np.einsum("ij,ij->i", xyz[:, 0], np.cross(xyz[:, 1], xyz[:, 2])).sum()
                / 6
            )
    assert total_volume == pytest.approx(3264, abs=1e-8, rel=0)


def test_programming_errors_are_not_reported_as_invalid_geometry(monkeypatch):
    from dtcc_core.builder.geometry_builders import prisms

    def broken(*args):
        raise RuntimeError("programming bug")

    monkeypatch.setattr(prisms, "_validate_walls", broken)
    with pytest.raises(RuntimeError, match="programming bug"):
        _decompose_prisms(_prism(BASE))


def test_one_shell_solid_is_supported_and_interior_shells_are_rejected():
    from dtcc_core.model import Solid

    surfaces = _input(_prism(COURTYARD), "surfaces").surfaces
    solid = Solid(surfaces=surfaces, shells=[np.arange(len(surfaces))])
    _assert_cells(_decompose_prisms(solid)[0], COURTYARD, {10: COURTYARD}, 3640)
    solid.shells.append(np.array([0]))
    with pytest.raises(_PrismError, match="not_2_5d"):
        _decompose_prisms(solid)


@pytest.mark.parametrize("case", ["indices", "shape", "markers", "unused_nonfinite"])
def test_invalid_or_untransferable_mesh_data_is_rejected(case):
    mesh = _input(_prism(BASE), "mesh")
    if case == "indices":
        mesh.faces[0, 0] = len(mesh.vertices)
    elif case == "shape":
        mesh.faces = np.array([[0, 1, 2, 3]])
    elif case == "markers":
        mesh.markers = np.zeros(len(mesh.faces), dtype=int)
    else:
        mesh.vertices = np.vstack([mesh.vertices, [np.nan, 0, 0]])
    with pytest.raises(_PrismError):
        _decompose_prisms(mesh)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"tolerance": 0},
        {"tolerance": -1},
        {"tolerance": np.inf},
        {"tolerance": np.nan},
        {"wall_angle_tolerance": -1},
        {"wall_angle_tolerance": 90},
        {"wall_angle_tolerance": np.nan},
    ],
)
def test_parameter_errors_precede_geometry_processing(kwargs):
    with pytest.raises(ValueError) as caught:
        _decompose_prisms(None, **kwargs)
    assert not isinstance(caught.value, _PrismError)


@pytest.mark.parametrize(
    "values",
    [
        np.zeros((3, 2)),
        np.zeros((3, 3)),
        np.array([[0, 0, 0], [1, 0, 0], [0, 1, np.nan]]),
        np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype=complex),
    ],
)
def test_invalid_raw_faces_are_rejected(values):
    with pytest.raises(_PrismError, match="degenerate"):
        _decompose_prisms([_PrismFace((values,), 42)])


def test_common_translation_and_srs_are_returned_without_double_application():
    geometry = _input(_prism(COURTYARD, 2, 12), "surfaces")
    geometry.transform.set_translation(100, 200, 30)
    geometry.transform.srs = "EPSG:3006"
    for surface in geometry.surfaces:
        surface.transform.set_translation(5, 6, 7)
        surface.transform.srs = "EPSG:3006"
    before = deepcopy(geometry)
    result = _decompose_prisms(geometry)[0]
    assert result.translation == (105, 206, 37)
    assert result.srs == "EPSG:3006"
    assert result.ground_z == 2 and result.levels[0][0] == 12
    assert result.base.equals(COURTYARD)
    for actual, expected in zip(geometry.surfaces, before.surfaces, strict=True):
        np.testing.assert_array_equal(actual.vertices, expected.vertices)
        np.testing.assert_array_equal(
            actual.transform.affine, expected.transform.affine
        )


@pytest.mark.parametrize(
    "case",
    [
        "scale",
        "tilt",
        "mixed_translation",
        "mixed_srs",
        "fields",
        "regions",
        "schema",
        "context",
    ],
)
def test_untransferable_frames_and_metadata_are_rejected(case):
    from dtcc_core.model import Field, SemanticRegion

    geometry = _input(_prism(BASE), "surfaces")
    if case == "scale":
        geometry.transform.affine[0, 0] = 2
    elif case == "tilt":
        geometry.transform.affine[:3, :3] = [[1, 0, 0], [0, 0, -1], [0, 1, 0]]
    elif case == "mixed_translation":
        geometry.surfaces[0].transform.set_translation(1, 0, 0)
    elif case == "mixed_srs":
        geometry.transform.srs = "EPSG:3006"
        geometry.surfaces[0].transform.srs = "EPSG:28992"
    elif case == "fields":
        geometry.surfaces[0].fields = [Field(name="temperature")]
    elif case == "regions":
        geometry.regions = [
            SemanticRegion(semantic_type="RoofSurface", indices=np.array([0]))
        ]
    elif case == "schema":
        geometry.schema_id = "example"
    else:
        geometry.dataset_context = object()
    with pytest.raises(_PrismError, match="unsupported_"):
        _decompose_prisms(geometry)
