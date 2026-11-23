"""Footprint construction preserves rings, placement and extrusion ownership."""

from collections import Counter
from copy import deepcopy
from uuid import UUID

import numpy as np
import pytest
from shapely.geometry import MultiPolygon, Polygon, box

from dtcc_core import io
from dtcc_core.builder import (
    build_lod1_buildings,
    building_from_footprint,
    set_building_heights_from_attribute,
)
from dtcc_core.builder.cleaning.blocks import invalidate_lod1
from dtcc_core.model import City, Raster, Surface

SQUARE = box(0, 0, 10, 10)
CONCAVE = Polygon([(0, 0), (10, 0), (10, 3), (3, 3), (3, 10), (0, 10)])
COURTYARD = Polygon(SQUARE.exterior.coords, [box(3, 3, 7, 7).exterior.coords])


def _polygon(surface):
    return Polygon(surface.vertices[:, :2], [h[:, :2] for h in surface.holes])


def _assert_closed_volume(building, expected):
    # The mesher triangulates the cap with its concavities and holes intact.
    mesh = building.lod1.mesh(weld=True)
    vertices, faces = mesh.vertices, mesh.faces
    directed = Counter()
    for face in faces:
        for a, b in zip(face, np.roll(face, -1), strict=True):
            directed[int(a), int(b)] += 1
    assert all(n == 1 and directed[b, a] == 1 for (a, b), n in directed.items())
    triangles = vertices[faces] - vertices[0]
    volume = (
        np.einsum(
            "ij,ij->i", triangles[:, 0], np.cross(triangles[:, 1], triangles[:, 2])
        ).sum()
        / 6
    )
    assert volume == pytest.approx(expected)


@pytest.mark.parametrize("polygon", [SQUARE, CONCAVE, COURTYARD])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("kind", ["polygon", "surface"])
def test_rings_height_and_closed_outward_solid(polygon, reverse, kind):
    if reverse:
        polygon = Polygon(
            list(polygon.exterior.coords)[::-1],
            [list(h.coords)[::-1] for h in polygon.interiors],
        )
    footprint = polygon if kind == "polygon" else Surface().from_polygon(polygon, 97)
    building = building_from_footprint(footprint, 10, ground_height=2)
    assert UUID(building.id)
    assert building.attributes["ground_height"] == 2
    assert building.attributes["measured_height"] == 10
    assert "estimated_height" not in building.attributes
    assert building.lod0.vertices[:, 2] == pytest.approx(12)
    assert all(h[:, 2] == pytest.approx(12) for h in building.lod0.holes)
    assert _polygon(building.lod0).equals(polygon)
    assert len(building.lod0.holes) == len(polygon.interiors)
    assert (building.lod1.bounds.zmin, building.lod1.bounds.zmax) == (2, 12)
    _assert_closed_volume(building, polygon.area * 10)


@pytest.mark.parametrize("dimensions", [2, 3])
@pytest.mark.parametrize("closed", [False, True])
def test_array_closure_duplicates_and_ignored_z(dimensions, closed):
    xy = np.array([[0, 0], [10, 0], [10, 0], [10, 10], [0, 10]], dtype=float)
    if closed:
        xy = np.vstack([xy, xy[0], xy[0]])
    values = xy if dimensions == 2 else np.column_stack([xy, np.arange(len(xy))])
    before = values.copy()
    building = building_from_footprint(values, 10, ground_height=2)
    np.testing.assert_array_equal(values, before)
    assert len(building.lod0.vertices) == 4
    assert _polygon(building.lod0).equals(SQUARE)
    _assert_closed_volume(building, 1000)


def test_surface_hole_normalization_and_frame_preserved_once():
    source = Surface().from_polygon(COURTYARD, 123)
    source.vertices = np.vstack(
        [source.vertices, source.vertices[-1], source.vertices[0]]
    )
    source.holes[0] = np.vstack(
        [source.holes[0], source.holes[0][-1], source.holes[0][0]]
    )
    source.transform.set_translation(100000, 200000, 30)
    source.transform.srs = "EPSG:3006"
    before = deepcopy(source)
    building = building_from_footprint(source, 10, ground_height=2)
    assert len(building.lod0.vertices) == len(building.lod0.holes[0]) == 4
    assert _polygon(building.lod0).equals(COURTYARD)
    for face in [building.lod0, *building.lod1.surfaces]:
        assert face.transform is not source.transform
        np.testing.assert_array_equal(face.transform.affine, source.transform.affine)
        assert face.transform.srs == "EPSG:3006"
        world = face.transform(face.vertices)
        np.testing.assert_allclose(world, face.vertices + [100000, 200000, 30], rtol=0)
    np.testing.assert_array_equal(building.transform.affine, np.eye(4))
    np.testing.assert_array_equal(building.lod1.transform.affine, np.eye(4))
    np.testing.assert_array_equal(source.vertices, before.vertices)
    np.testing.assert_array_equal(source.holes[0], before.holes[0])
    building.lod0.transform.set_translation(0, 0, 0)
    np.testing.assert_array_equal(source.transform.affine, before.transform.affine)
    assert building.lod1.surfaces[0].transform.offset[0] == 100000


@pytest.mark.parametrize("height_attribute", ["measured_height", "estimated_height"])
def test_attributes_are_independent_and_reserved_values_are_removed(height_attribute):
    reserved = {
        "measured_height",
        "estimated_height",
        "ground_height",
        "lod1_extrusion",
        "merged_from",
        "source_footprint",
        "footprint_walls",
    }
    attributes = dict.fromkeys(reserved, "stale")
    attributes.update(source_id="upstream", custom={"values": [1]})
    before = deepcopy(attributes)
    building = building_from_footprint(
        SQUARE,
        10,
        ground_height=2,
        id="a",
        attributes=attributes,
        height_attribute=height_attribute,
        lod1=False,
    )
    assert building.id == "a"
    assert building.attributes == {
        "ground_height": 2,
        height_attribute: 10,
        "source_id": "upstream",
        "custom": {"values": [1]},
    }
    assert building.lod1 is None
    building.attributes["custom"]["values"].append(2)
    assert attributes == before
    build_lod1_buildings([building], rebuild=False)
    assert "lod1_extrusion" in building.attributes
    _assert_closed_volume(building, 1000)


@pytest.mark.parametrize(
    "footprint",
    [
        np.empty((0, 2)),
        np.array([0, 1, 2]),
        np.zeros((4, 4)),
        np.zeros((4, 1)),
        np.zeros((2, 3)),
        np.ones((4, 2)),
        np.array([[0, 0], [1, 0], [2, 0]]),
        np.array([[0, 0], [2, 2], [0, 2], [2, 0]]),
        np.array([[0, 0], [1, np.inf], [0, 1]]),
        np.array([[0, 0, 0], [1, 0, np.nan], [0, 1, 0]]),
        np.array([[0, 0], [1, 0], [0, 1]], dtype=complex),
        Polygon(),
        Surface(),
        Polygon(SQUARE.exterior.coords, [box(20, 20, 21, 21).exterior.coords]),
        Surface(
            vertices=np.array([[0, 0, 1], [1, 0, 1], [0, 1, 1]]),
            holes=[np.zeros((3, 3))],
        ),
    ],
)
def test_rejects_invalid_rings(footprint):
    with pytest.raises(ValueError):
        building_from_footprint(footprint, 10)


def test_invalid_polygon_includes_validity_reason():
    with pytest.raises(ValueError, match="Self-intersection"):
        building_from_footprint(Polygon([(0, 0), (2, 2), (0, 2), (2, 0)]), 10)


def test_preserves_small_positive_area_and_boundary_detail():
    polygon = Polygon([(0, 0), (0.01, 0), (0.01, 0.01), (0.009, 0.009), (0, 0.01)])
    building = building_from_footprint(polygon, 10)
    np.testing.assert_array_equal(
        building.lod0.vertices[:, :2], np.asarray(polygon.exterior.coords)[:-1]
    )
    assert len(building.lod1.surfaces) == 7
    assert _polygon(building.lod0).area == polygon.area


@pytest.mark.parametrize(
    "footprint", [MultiPolygon([SQUARE]), "square", [[0, 0], [1, 0], [0, 1]]]
)
def test_rejects_unsupported_types(footprint):
    with pytest.raises(TypeError, match="Polygon, Surface or numpy array"):
        building_from_footprint(footprint, 10)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"height": 0},
        {"height": -1},
        {"height": np.nan},
        {"height": np.inf},
        {"ground_height": np.nan},
        {"ground_height": -np.inf},
        {"height": 1e308, "ground_height": 1e308},
        {"height": 1, "ground_height": 1e30},
        {"height_attribute": "height"},
    ],
)
def test_rejects_invalid_heights_and_attribute(kwargs):
    with pytest.raises(ValueError):
        building_from_footprint(SQUARE, **({"height": 10} | kwargs))


@pytest.mark.parametrize("terrain", [None, Raster(data=np.full((20, 20), 2.0))])
def test_height_builder_invalidates_owned_extrusion(terrain):
    building = building_from_footprint(SQUARE, 10, ground_height=2)
    original = building.lod1
    build_lod1_buildings([building], rebuild=False)
    assert building.lod1 is original
    building.attributes["measured_height"] = 20
    set_building_heights_from_attribute(
        [building],
        terrain,
        always_use_default_ground=terrain is None,
        default_ground_height=2,
    )
    assert building.lod0.zmax == 22 and building.attributes["ground_height"] == 2
    assert building.lod1 is None and "lod1_extrusion" not in building.attributes
    build_lod1_buildings([building], rebuild=False)
    assert (building.lod1.bounds.zmin, building.lod1.bounds.zmax) == (2, 22)
    _assert_closed_volume(building, 2000)


@pytest.mark.parametrize("translated", [False, True])
def test_native_roundtrip_preserves_geometry_and_lifecycle(tmp_path, translated):
    footprint = Surface().from_polygon(COURTYARD, 100)
    if translated:
        footprint.transform.set_translation(100000, 200000, 30)
        footprint.transform.srs = "EPSG:3006"
    building = building_from_footprint(footprint, 10, ground_height=2, id="courtyard")
    path = tmp_path / "building.dtcc"
    io.save_model(building, path)
    restored = io.load_model(path)
    assert restored.id == building.id and restored.attributes == building.attributes
    np.testing.assert_array_equal(restored.lod0.vertices, building.lod0.vertices)
    np.testing.assert_array_equal(restored.lod0.holes[0], building.lod0.holes[0])
    np.testing.assert_array_equal(
        restored.lod0.transform.affine, footprint.transform.affine
    )
    assert restored.lod0.transform.srs == footprint.transform.srs
    for actual, expected in zip(
        restored.lod1.surfaces, building.lod1.surfaces, strict=True
    ):
        np.testing.assert_array_equal(actual.vertices, expected.vertices)
        np.testing.assert_array_equal(
            actual.transform.affine, footprint.transform.affine
        )
        assert actual.transform.srs == footprint.transform.srs
        for a, e in zip(actual.holes, expected.holes, strict=True):
            np.testing.assert_array_equal(a, e)
    original = restored.lod1
    build_lod1_buildings([restored], rebuild=False)
    assert restored.lod1 is original
    restored.attributes["ground_height"] = 1
    invalidate_lod1(restored)
    assert restored.lod1 is None


@pytest.mark.parametrize("strict", [False, True])
def test_cityjson_quantization_preserves_both_lods_and_requires_explicit_rebuild(
    tmp_path,
    strict,
):
    polygon = box(100000.1234, 200000.2345, 100010.3456, 200010.4567)
    building = building_from_footprint(
        polygon,
        10.1234,
        ground_height=2.1234,
        id="quantized",
        attributes={"source_id": "upstream", "custom": {"name": "a"}},
    )
    city = City()
    city.add_building(building)
    path = tmp_path / "city.city.json"
    city.save(path, strict=strict)
    restored = io.load_city(path, strict=strict).buildings[0]
    assert restored.id == building.id and restored.attributes == building.attributes
    assert restored.lod0 is not None
    np.testing.assert_allclose(
        restored.lod0.surfaces[0].vertices, building.lod0.vertices, atol=0.001, rtol=0
    )
    for a, e in zip(restored.lod1.surfaces, building.lod1.surfaces, strict=True):
        np.testing.assert_allclose(a.vertices, e.vertices, atol=0.001, rtol=0)
    original = restored.lod1
    representation_ids = set(restored.geometry)
    invalidate_lod1(restored)
    assert restored.lod1 is original and "lod1_extrusion" not in restored.attributes
    restored.attributes["measured_height"] = 20
    set_building_heights_from_attribute(
        [restored],
        None,
        always_use_default_ground=True,
        default_ground_height=2.1234,
    )
    build_lod1_buildings([restored], rebuild=False)
    assert restored.lod1 is original
    build_lod1_buildings([restored], rebuild=True)
    assert restored.lod1 is not original
    assert set(restored.geometry) == representation_ids
    assert restored.lod1.bounds.zmax == pytest.approx(22.1234)
    assert "lod1_extrusion" in restored.attributes


def test_cityjson_lod0_only_preserves_courtyard_within_quantization(tmp_path):
    building = building_from_footprint(
        COURTYARD, 10.1234, ground_height=2.1234, lod1=False
    )
    city = City()
    city.add_building(building)
    path = tmp_path / "footprint.city.json"
    city.save(path)
    restored = io.load_city(path).buildings[0]
    assert restored.lod1 is None
    actual = restored.lod0.surfaces[0]
    np.testing.assert_allclose(
        actual.vertices, building.lod0.vertices, rtol=0, atol=0.001
    )
    np.testing.assert_allclose(
        actual.holes[0], building.lod0.holes[0], rtol=0, atol=0.001
    )
