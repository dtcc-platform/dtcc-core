"""Tests for the GeoJSON web writer of Mesh."""

import json
import zipfile
from copy import deepcopy

import numpy as np
import pytest
from shapely.geometry import shape

import dtcc_core.datasets as datasets
from dtcc_core.datasets import attach_dataset_context, load_model_package
from dtcc_core.io.meshes import centered_copy
from dtcc_core.model import Mesh
from dtcc_core.model.geometry.transform import Transform

# Oracles: pyproj EPSG:3006 -> EPSG:4326 (always_xy) in Gothenburg.
SQUARE = [(320000, 6399000), (320100, 6399000), (320100, 6399100), (320000, 6399100)]
SQUARE_WGS84_FIRST = [11.979646369152894, 57.697102381277034]
NATIVE_CRS = {"type": "name", "properties": {"name": "EPSG:3006"}}


def _mesh(markers=(0, 0, 1, 1), z=0.0):
    """A flat square split at x=320050 into two triangle pairs, one per marker."""
    vertices = np.array(
        [
            [320000.0, 6399000.0, z],
            [320050.0, 6399000.0, z],
            [320100.0, 6399000.0, z],
            [320100.0, 6399100.0, z],
            [320050.0, 6399100.0, z],
            [320000.0, 6399100.0, z],
        ]
    )
    faces = np.array([[0, 1, 4], [0, 4, 5], [1, 2, 3], [1, 3, 4]])
    mesh = Mesh(vertices=vertices, faces=faces, markers=np.array(markers))
    mesh.transform.srs = "EPSG:3006"
    return mesh


def _context(name):
    dataset = getattr(datasets, name)
    return dataset.create_context(dataset.validate({"bounds": (0.0, 0.0, 1.0, 1.0)}))


def _triangle_area(mesh, face):
    (x0, y0, _), (x1, y1, _), (x2, y2, _) = mesh.vertices[face]
    return abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)) / 2


def _signed_area(ring):
    return (
        sum(a[0] * b[1] - b[0] * a[1] for a, b in zip(ring, ring[1:], strict=False)) / 2
    )


def _positions(payload):
    return np.array(
        [
            position
            for feature in payload["features"]
            for polygon in feature["geometry"]["coordinates"]
            for ring in polygon
            for position in ring
        ]
    )


def _geojson_payload(package):
    [artifact] = [a for a in package.artifacts if a.format == "geojson"]
    if package.package_format == "dtccpkg":
        with zipfile.ZipFile(package.path) as archive:
            return json.loads(archive.read(artifact.path))
    return json.loads((package.path / artifact.path).read_text())


def _tilted(mesh):
    angle = np.radians(30.0)
    rotation = np.eye(4)
    rotation[1, 1] = rotation[2, 2] = np.cos(angle)
    rotation[1, 2] = -np.sin(angle)
    rotation[2, 1] = np.sin(angle)
    mesh.transform = Transform(srs="EPSG:3006", affine=rotation)
    return mesh


def test_one_dissolved_multipolygon_per_marker_native():
    mesh = _mesh()
    payload = mesh.to_geojson()

    assert payload["crs"] == NATIVE_CRS
    assert [f["properties"] for f in payload["features"]] == [
        {"marker": 0, "face_count": 2},
        {"marker": 1, "face_count": 2},
    ]
    for feature in payload["features"]:
        assert feature["geometry"]["type"] == "MultiPolygon"
        assert "id" not in feature
        [polygon] = feature["geometry"]["coordinates"]
        assert len(polygon) == 1  # the shared edge dissolves, no interior ring
        exterior = polygon[0]
        assert exterior[0] == exterior[-1]
        assert _signed_area(exterior) > 0  # counter-clockwise
        assert all(len(position) == 2 for position in exterior)

    areas = [shape(f["geometry"]).area for f in payload["features"]]
    triangle_sum = sum(_triangle_area(mesh, face) for face in mesh.faces)
    assert sum(areas) == pytest.approx(triangle_sum, rel=0, abs=1e-6)
    assert triangle_sum == pytest.approx(100.0 * 100.0)
    assert {tuple(p) for p in _positions(payload)} == {
        (float(x), float(y)) for x, y in SQUARE
    } | {(320050.0, 6399000.0), (320050.0, 6399100.0)}


def test_wgs84_output_has_no_crs_member():
    payload = _mesh().to_geojson(crs="EPSG:4326")

    assert "crs" not in payload
    found = _positions(payload)
    assert (found[:, 0] > -180).all() and (found[:, 0] < 180).all()
    assert (found[:, 1] > -90).all() and (found[:, 1] < 90).all()
    assert any(
        np.allclose(position, SQUARE_WGS84_FIRST, rtol=0, atol=1e-9)
        for position in found
    )


def test_mesh_without_markers_yields_one_feature_with_null_marker():
    payload = _mesh(markers=()).to_geojson()

    [feature] = payload["features"]
    assert feature["properties"] == {"marker": None, "face_count": 4}
    assert shape(feature["geometry"]).area == pytest.approx(
        100.0 * 100.0, rel=0, abs=1e-6
    )


def test_empty_mesh_yields_empty_collection():
    mesh = Mesh()
    mesh.transform.srs = "EPSG:3006"

    assert mesh.to_geojson(crs="EPSG:4326")["features"] == []


def test_markers_of_wrong_length_raise():
    with pytest.raises(ValueError, match="markers"):
        _mesh(markers=(0, 1)).to_geojson()


def _grid_mesh(centre_marker=7, other_marker=-3):
    """A 5x5 grid of 10 m cells, two triangles each; the centre cell is separate."""
    n = 5
    xs, ys = np.meshgrid(np.arange(n + 1) * 10.0 + 320000, np.arange(n + 1) * 10.0)
    vertices = np.column_stack([xs.ravel(), ys.ravel() + 6399000, np.zeros(xs.size)])
    faces, markers = [], []
    for row in range(n):
        for col in range(n):
            a = row * (n + 1) + col
            b, c, d = a + 1, a + n + 2, a + n + 1
            faces += [[a, b, c], [a, c, d]]
            centre = (row, col) == (2, 2)
            markers += [centre_marker if centre else other_marker] * 2
    mesh = Mesh(vertices=vertices, faces=np.array(faces), markers=np.array(markers))
    mesh.transform.srs = "EPSG:3006"
    return mesh


@pytest.mark.parametrize("crs", [None, "EPSG:4326"])
def test_surrounding_marker_keeps_the_surrounded_one_as_an_interior_ring(crs):
    mesh = _grid_mesh()
    payload = mesh.to_geojson(crs=crs)

    outer, inner = payload["features"]  # ordered by marker: -3, then 7
    assert [outer["properties"]["marker"], inner["properties"]["marker"]] == [-3, 7]
    [outer_polygon] = outer["geometry"]["coordinates"]
    [inner_polygon] = inner["geometry"]["coordinates"]
    assert len(outer_polygon) == 2
    assert len(inner_polygon) == 1
    exterior, hole = outer_polygon
    assert _signed_area(exterior) > 0
    assert hole[0] == hole[-1]
    assert _signed_area(hole) < 0  # clockwise
    assert _signed_area(inner_polygon[0]) > 0

    if crs is None:
        areas = {
            marker: sum(
                _triangle_area(mesh, face)
                for face, m in zip(mesh.faces, mesh.markers, strict=True)
                if m == marker
            )
            for marker in (-3, 7)
        }
        assert areas == {-3: 2400.0, 7: 100.0}
        assert shape(outer["geometry"]).area == pytest.approx(2400.0, rel=0, abs=1e-6)
        assert shape(inner["geometry"]).area == pytest.approx(100.0, rel=0, abs=1e-6)


def test_group_of_one_degenerate_triangle_raises():
    mesh = Mesh(
        vertices=np.array(
            [
                [320000.0, 6399000.0, 0],
                [320010.0, 6399000.0, 0],
                [320005.0, 6399000.0, 0],
            ]
        ),
        faces=np.array([[0, 1, 2]]),
        markers=np.array([4]),
    )
    mesh.transform.srs = "EPSG:3006"
    with pytest.raises(ValueError, match="no area"):
        mesh.to_geojson()


def test_varying_z_raises():
    mesh = _mesh()
    mesh.vertices[0, 2] = 1.0
    with pytest.raises(ValueError, match="flat"):
        mesh.to_geojson(crs="EPSG:4326")


def test_flat_mesh_under_tilted_transform_raises():
    with pytest.raises(ValueError, match="flat"):
        _tilted(_mesh()).to_geojson(crs="EPSG:4326")


@pytest.mark.parametrize("make", [lambda m: m, _tilted], ids=["varying_z", "tilted"])
def test_refusals_hold_after_canonical_reload(tmp_path, make):
    mesh = _mesh()
    if make is _tilted:
        _tilted(mesh)
    else:
        mesh.vertices[0, 2] = 1.0
    attach_dataset_context(mesh, _context("city_flat_mesh"))
    package = mesh.export(tmp_path / "pkg", canonical=True)

    with pytest.raises(ValueError, match="flat"):
        load_model_package(package.path).to_geojson(crs="EPSG:4326")


def test_centered_mesh_exports_original_positions_without_mutation():
    mesh = _mesh()
    centered = centered_copy(mesh)
    assert not np.array_equal(centered.transform.affine, np.eye(4))
    before = centered.vertices.copy()

    expected = _positions(mesh.to_geojson(crs="EPSG:4326"))
    found = _positions(centered.to_geojson(crs="EPSG:4326"))

    np.testing.assert_allclose(found, expected, rtol=0, atol=1e-9)
    np.testing.assert_array_equal(centered.vertices, before)
    assert not np.array_equal(centered.transform.affine, np.eye(4))


@pytest.mark.parametrize("target", ["pkg", "pkg.dtccpkg"])
def test_canonical_geojson_package_round_trip(tmp_path, target):
    mesh = _mesh()
    attach_dataset_context(mesh, _context("city_flat_mesh"))
    package = mesh.export(tmp_path / target, canonical=True, format="geojson")
    again = load_model_package(package.path).export(
        tmp_path / f"again_{target}", canonical=True, format="geojson"
    )

    payload = _geojson_payload(package)
    assert payload == _geojson_payload(again)
    assert "crs" not in payload
    assert any(
        np.allclose(position, SQUARE_WGS84_FIRST, rtol=0, atol=1e-9)
        for position in _positions(payload)
    )
    [artifact] = [a for a in package.artifacts if a.format == "geojson"]
    assert artifact.crs == "EPSG:4326"
    assert artifact.bounds is not None


@pytest.mark.parametrize("target", ["pkg", "pkg.dtccpkg"])
def test_non_flat_mesh_package_export_publishes_nothing(tmp_path, target):
    mesh = _mesh()
    mesh.vertices[:, 2] = np.arange(len(mesh.vertices), dtype=float)
    attach_dataset_context(mesh, _context("city_surface_mesh"))

    with pytest.raises(ValueError, match="flat"):
        mesh.export(tmp_path / target, canonical=True, format="geojson")

    assert list(tmp_path.iterdir()) == []


def test_export_leaves_the_mesh_unmodified():
    mesh = _mesh()
    before = deepcopy(mesh)
    mesh.to_geojson(crs="EPSG:4326")
    np.testing.assert_array_equal(mesh.vertices, before.vertices)
    np.testing.assert_array_equal(mesh.transform.affine, before.transform.affine)
