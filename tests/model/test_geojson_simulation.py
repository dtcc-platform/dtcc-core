"""GeoJSON writers of FieldSlice and StreamlineCollection."""

from __future__ import annotations

import json
import zipfile

import numpy as np
import pytest
from pyproj import Transformer

import dtcc_core.datasets as datasets
from dtcc_core.datasets import attach_dataset_context, load_model_package
from dtcc_core.datasets.schema import (
    DatasetContext,
    DatasetIdentity,
    DatasetMetadata,
    DatasetPresentation,
    DatasetProvenance,
    DatasetRequest,
)
from dtcc_core.model import (
    Bounds,
    Field,
    FieldSlice,
    LineString,
    StreamlineCollection,
    exchange,
)

NATIVE_CRS = {"type": "name", "properties": {"name": "EPSG:3006"}}
BOUNDS = (320000.0, 6399000.0, 320100.0, 6399100.0)
DOMAIN = Bounds(
    xmin=320000.0, ymin=6399000.0, zmin=2.0, xmax=320100.0, ymax=6399100.0, zmax=40.0
)
# Oracle: pyproj EPSG:3006 -> EPSG:4326 of (320000, 6399000).
LONLAT = [11.979646369152894, 57.697102381277034]
UNITS = {"velocity": "m/s", "speed": "m/s", "pressure": "Pa", "temperature": "K"}


def _lonlat(x, y):
    return list(Transformer.from_crs(3006, 4326, always_xy=True).transform(x, y))


def _fields(count, association):
    values = np.arange(1.0, 1.0 + 3 * count).reshape((count, 3))
    columns = {
        "velocity": values,
        "speed": np.linalg.norm(values, axis=1).reshape((count, 1)),
        "pressure": values[:, :1] * 10.0,
        "temperature": values[:, :1] + 270.0,
    }
    return [
        Field(
            name=name,
            unit=UNITS[name],
            values=rows,
            dim=rows.shape[1],
            association=association,
        )
        for name, rows in columns.items()
    ]


def _slice(**kwargs):
    points = np.array(
        [
            [320000.0, 6399000.0, 10.0],
            [320100.0, 6399000.0, 10.0],
            [320000.0, 6399100.0, 10.0],
            [320100.0, 6399100.0, 10.0],
        ]
    )
    result = FieldSlice(
        points=points,
        resolution=2,
        time=1.5,
        period=8.0,
        crs="EPSG:3006",
        domain_bounds=DOMAIN,
        **kwargs,
    )
    result.fields = _fields(4, "sample")
    return result


def _streamlines(**kwargs):
    lines = []
    for index in range(2):
        base = np.array([320000.0 + 10 * index, 6399000.0, 10.0])
        line = LineString(
            vertices=base + np.array([[0, 0, 0], [10, 20, 1], [20, 40, 2]])
        )
        line.fields = _fields(3, "vertex")
        lines.append(line)
    return StreamlineCollection(
        lines=lines,
        seed_axis="z",
        seed_position=0.5,
        requested_line_count=2,
        streamline_steps=3,
        time=1.5,
        period=8.0,
        crs="EPSG:3006",
        domain_bounds=DOMAIN,
        **kwargs,
    )


PRODUCTS = pytest.mark.parametrize("make", [_slice, _streamlines])


def _positions(payload):
    """Every position of a payload as an (N, d) array."""
    rows = []
    for feature in payload["features"]:
        coordinates = feature["geometry"]["coordinates"]
        rows.extend(
            coordinates
            if feature["geometry"]["type"] == "LineString"
            else [coordinates]
        )
    return np.array(rows)


def _context(crs):
    return DatasetContext(
        identity=DatasetIdentity(name="example", title="Example"),
        metadata=DatasetMetadata(crs=crs),
        provenance=DatasetProvenance(sources=["example"]),
        presentation=DatasetPresentation(),
        request=DatasetRequest(dataset_name="example", bounds=[0.0, 0.0, 1.0, 1.0]),
        health={"status": "ok"},
        warnings=[],
    )


@PRODUCTS
def test_native_output_names_the_source_crs(make):
    model = make()
    payload = model.to_geojson()

    assert payload["crs"] == NATIVE_CRS
    assert payload["metadata"]["crs"] == "EPSG:3006"
    assert payload["metadata"]["bounds"] == [
        320000.0, 6399000.0, 2.0, 320100.0, 6399100.0, 40.0
    ]  # fmt: skip
    assert payload["metadata"]["time"] == 1.5
    assert payload["metadata"]["time_period"] == 8.0
    assert payload["name"] == model.name
    stored = model.points if make is _slice else model.lines[0].vertices
    np.testing.assert_array_equal(_positions(payload)[: len(stored)], stored)
    json.dumps(payload, allow_nan=False)


@PRODUCTS
def test_wgs84_output_is_reprojected(make):
    model = make()
    payload = model.to_geojson(crs="EPSG:4326")

    assert "crs" not in payload
    assert payload["metadata"]["crs"] == "EPSG:4326"
    positions = _positions(payload)
    assert positions.shape[1] == 3
    assert np.all(np.abs(positions[:, 0]) <= 180) and np.all(
        np.abs(positions[:, 1]) <= 90
    )
    assert _lonlat(320000.0, 6399000.0) == pytest.approx(LONLAT, abs=1e-12)
    assert positions[0, :2] == pytest.approx(LONLAT, abs=1e-9)

    xmin, ymin, zmin, xmax, ymax, zmax = payload["metadata"]["bounds"]
    assert (zmin, zmax) == (2.0, 40.0)
    assert xmin <= positions[:, 0].min() and positions[:, 0].max() <= xmax
    assert ymin <= positions[:, 1].min() and positions[:, 1].max() <= ymax
    assert [xmin, ymin] == pytest.approx(_lonlat(320000.0, 6399000.0), abs=1e-3)
    assert [xmax, ymax] == pytest.approx(_lonlat(320100.0, 6399100.0), abs=1e-3)
    json.dumps(payload, allow_nan=False)


def test_other_target_is_reprojected_and_named():
    payload = _slice().to_geojson(crs="EPSG:3857")

    assert payload["crs"] == {"type": "name", "properties": {"name": "EPSG:3857"}}
    assert payload["metadata"]["crs"] == "EPSG:3857"
    expected = Transformer.from_crs(3006, 3857, always_xy=True).transform(
        320000.0, 6399000.0
    )
    assert _positions(payload)[0, :2] == pytest.approx(expected, abs=1e-6)


def test_slice_properties_carry_every_field_with_units():
    payload = _slice().to_geojson(crs="EPSG:4326")

    first, second = payload["features"][:2]
    assert first["properties"] == {
        "sample_index": 0,
        "velocity": [1.0, 2.0, 3.0],
        "velocity_unit": "m/s",
        "velocity_components": ["u", "v", "w"],
        "speed": pytest.approx(np.sqrt(14.0)),
        "speed_unit": "m/s",
        "pressure": 10.0,
        "pressure_unit": "Pa",
        "temperature": 271.0,
        "temperature_unit": "K",
    }
    assert second["properties"]["velocity"] == [4.0, 5.0, 6.0]
    assert not {"u", "v", "w"} & set(first["properties"])
    assert payload["metadata"]["fields"] == list(UNITS)


def test_streamline_properties_carry_every_field_with_units():
    payload = _streamlines().to_geojson(crs="EPSG:4326")

    [first, second] = payload["features"]
    assert (first["properties"]["seed_index"], second["properties"]["seed_index"]) == (
        0,
        1,
    )
    properties = first["properties"]
    assert properties["num_points"] == 3
    assert properties["velocity"] == [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
    assert properties["velocity_components"] == ["u", "v", "w"]
    assert properties["temperature"] == [271.0, 274.0, 277.0]
    assert properties["speed"] == pytest.approx(
        np.linalg.norm([[1, 2, 3], [4, 5, 6], [7, 8, 9]], axis=1)
    )
    assert {name: properties[f"{name}_unit"] for name in UNITS} == UNITS
    assert "speed_components" not in properties


def test_short_lines_are_skipped_and_keep_their_seed_index():
    model = _streamlines()
    short = LineString(vertices=np.array([[320000.0, 6399000.0, 1.0]]))
    short.fields = _fields(1, "vertex")
    model.lines.insert(0, short)

    payload = model.to_geojson()

    assert [f["properties"]["seed_index"] for f in payload["features"]] == [1, 2]


def test_nan_values_become_null():
    model = _slice()
    model.fields[0].values = np.array([[np.nan, 1.0, 2.0]] * 4)
    payload = model.to_geojson(crs="EPSG:4326")
    assert payload["features"][0]["properties"]["velocity"] == [None, 1.0, 2.0]

    model = _slice()
    model.fields[3].values = np.array([[np.nan], [1.0], [2.0], [3.0]])
    payload = model.to_geojson(crs="EPSG:4326")
    assert payload["features"][0]["properties"]["temperature"] is None
    json.dumps(payload, allow_nan=False)

    model = _streamlines()
    model.lines[0].fields[3].values = np.array([[np.nan], [1.0], [2.0]])
    payload = model.to_geojson(crs="EPSG:4326")
    assert payload["features"][0]["properties"]["temperature"] == [None, 1.0, 2.0]
    json.dumps(payload, allow_nan=False)


@PRODUCTS
def test_include_z(make):
    model = make()

    flat = _positions(model.to_geojson(include_z=False, crs="EPSG:4326"))
    full = _positions(model.to_geojson(include_z=True, crs="EPSG:4326"))

    assert flat.shape[1] == 2 and full.shape[1] == 3
    np.testing.assert_array_equal(flat, full[:, :2])
    stored = model.points if make is _slice else model.lines[0].vertices
    np.testing.assert_array_equal(full[: len(stored), 2], stored[:, 2])


def test_slice_default_follows_its_include_z():
    assert _positions(_slice(include_z=False).to_geojson()).shape[1] == 2
    assert _positions(_slice().to_geojson()).shape[1] == 3


def test_the_slice_transform_is_applied_once_without_modifying_it():
    model = _slice()
    model.transform.srs = "EPSG:3006"
    model.transform.set_rotation(
        np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    )
    model.transform.set_translation(1000.0, 2000.0, 3.0)
    stored = model.points.copy()
    affine = model.transform.affine.copy()

    payload = model.to_geojson()

    np.testing.assert_allclose(_positions(payload), model.transform(stored))
    np.testing.assert_array_equal(model.points, stored)
    np.testing.assert_array_equal(model.transform.affine, affine)
    wgs84 = _positions(model.to_geojson(crs="EPSG:4326"))
    expected = _lonlat(*model.transform(stored)[0, :2])
    assert wgs84[0, :2] == pytest.approx(expected, abs=1e-9)


def test_the_slice_bounds_follow_its_transform():
    model = _slice()
    model.transform.srs = "EPSG:3006"
    model.transform.set_rotation(
        np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    )
    model.transform.set_translation(1000.0, 2000.0, 3.0)
    corners = model.transform(
        np.array(
            [
                [x, y, z]
                for x in (DOMAIN.xmin, DOMAIN.xmax)
                for y in (DOMAIN.ymin, DOMAIN.ymax)
                for z in (DOMAIN.zmin, DOMAIN.zmax)
            ]
        )
    )
    envelope = [*corners.min(axis=0), *corners.max(axis=0)]
    assert envelope[2] == 5.0 and envelope[5] == 43.0

    native = model.to_geojson()
    assert native["metadata"]["bounds"] == pytest.approx(envelope)

    wgs84 = model.to_geojson(crs="EPSG:4326")
    xmin, ymin, zmin, xmax, ymax, zmax = wgs84["metadata"]["bounds"]
    assert (zmin, zmax) == (5.0, 43.0)
    assert [xmin, ymin] == pytest.approx(_lonlat(*corners[:, :2].min(axis=0)), abs=1e-3)
    for payload in (native, wgs84):
        x, y, z = _positions(payload).T
        low, high = payload["metadata"]["bounds"][:3], payload["metadata"]["bounds"][3:]
        assert low[0] <= x.min() and x.max() <= high[0]
        assert low[1] <= y.min() and y.max() <= high[1]
        assert low[2] <= z.min() and z.max() <= high[2]


@pytest.mark.parametrize("crs", [None, "EPSG:4326"])
def test_vector_fields_follow_the_slice_transform(crs):
    model = _slice()
    model.transform.srs = "EPSG:3006"
    model.transform.set_rotation(
        np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    )

    properties = model.to_geojson(crs=crs)["features"][1]["properties"]

    # Rotating (4, 5, 6) by 90 degrees about z gives (-5, 4, 6).
    assert properties["velocity"] == pytest.approx([-5.0, 4.0, 6.0])
    assert properties["speed"] == pytest.approx(np.linalg.norm([4.0, 5.0, 6.0]))
    assert properties["pressure"] == 40.0
    assert model.fields[0].values[1].tolist() == [4.0, 5.0, 6.0]


def _rotated(rotation, **kwargs):
    model = _slice(**kwargs)
    model.transform.srs = "EPSG:3006"
    model.transform.set_rotation(rotation)
    return model


Z_ROTATION = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
X_ROTATION = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]])
SCALE = 1000.0 * np.eye(3)


def test_vector_fields_need_an_orthogonal_transform():
    model = _rotated(SCALE)
    with pytest.raises(ValueError, match="orthogonal"):
        model.to_geojson()

    model.fields = [f for f in model.fields if f.dim == 1]
    payload = model.to_geojson()
    assert payload["features"][1]["properties"]["pressure"] == 40.0
    assert _positions(payload)[1, 0] == pytest.approx(1000.0 * 320100.0)

    reflection = np.diag([1.0, 1.0, -1.0])
    properties = _rotated(reflection).to_geojson()["features"][1]["properties"]
    assert properties["velocity"] == pytest.approx([4.0, 5.0, -6.0])


def test_two_component_fields_need_the_xy_plane_kept_separate_from_z():
    field = Field(
        name="flow",
        unit="m/s",
        values=np.array([[0.0, 1.0]] * 4),
        dim=2,
        association="sample",
    )
    model = _rotated(Z_ROTATION)
    model.fields.append(field)
    properties = model.to_geojson()["features"][0]["properties"]
    assert properties["flow"] == pytest.approx([-1.0, 0.0])
    assert properties["flow_components"] == ["u", "v"]

    model = _rotated(X_ROTATION)
    model.fields.append(field)
    with pytest.raises(ValueError, match="flow.*x/y plane"):
        model.to_geojson()


def test_wide_fields_are_labelled_by_index_and_never_transformed():
    model = _rotated(Z_ROTATION)
    model.fields.append(
        Field(
            name="modes",
            unit="1",
            values=np.arange(16.0).reshape((4, 4)),
            dim=4,
            association="sample",
        )
    )

    properties = model.to_geojson()["features"][1]["properties"]

    assert properties["modes"] == [4.0, 5.0, 6.0, 7.0]
    assert properties["modes_components"] == ["c0", "c1", "c2", "c3"]


def test_a_geometry_field_cannot_overwrite_reserved_keys():
    for name in ("bounds", "time"):
        model = _slice()
        model.fields.append(
            Field(name=name, unit="m", values=np.array([1.0]), association="geometry")
        )
        with pytest.raises(ValueError, match=name):
            model.to_geojson()

    model = _streamlines()
    model.lines[0].fields.append(
        Field(
            name="seed_index", unit="1", values=np.array([7.0]), association="geometry"
        )
    )
    with pytest.raises(ValueError, match="seed_index"):
        model.to_geojson()


def test_a_row_aligned_field_cannot_overwrite_reserved_or_derived_keys():
    model = _slice()
    model.fields.append(
        Field(
            name="sample_index", unit="1", values=np.ones((4, 1)), association="sample"
        )
    )
    with pytest.raises(ValueError, match="sample_index"):
        model.to_geojson()

    model = _slice()
    model.fields.append(
        Field(name="speed_unit", unit="1", values=np.ones((4, 1)), association="sample")
    )
    with pytest.raises(ValueError, match="speed_unit"):
        model.to_geojson()

    for name in ("seed_index", "num_points"):
        model = _streamlines()
        model.lines[0].fields.append(
            Field(name=name, unit="1", values=np.ones((3, 1)), association="vertex")
        )
        with pytest.raises(ValueError, match=name):
            model.to_geojson()


def _empty_slice(shape):
    return FieldSlice(points=np.empty(shape), crs="EPSG:3006")


@pytest.mark.parametrize("shape", [(0,), (0, 3)])
@pytest.mark.parametrize("restored", [False, True], ids=["original", "restored"])
def test_an_empty_slice_exports_no_features(shape, restored):
    model = _empty_slice(shape)
    if restored:
        model = exchange.loads(exchange.dumps(model))

    for crs in (None, "EPSG:4326"):
        payload = model.to_geojson(crs=crs)
        assert payload["features"] == []
        assert payload["metadata"]["sample_count"] == 0
        json.dumps(payload, allow_nan=False)


@pytest.mark.usefixtures("no_network")
@pytest.mark.parametrize("shape", [(0,), (0, 3)])
def test_canonical_geojson_exports_an_empty_slice(tmp_path, shape):
    model = _empty_slice(shape)
    attach_dataset_context(model, _smoke("slice").dataset_context)

    package = model.export(tmp_path / "pkg", canonical=True, format="geojson")
    again = load_model_package(package.path).export(
        tmp_path / "again", canonical=True, format="geojson"
    )

    payload = _payload(package)
    assert payload == _payload(again)
    assert payload["features"] == []
    assert "crs" not in payload
    [artifact] = [a for a in package.artifacts if a.format == "geojson"]
    assert artifact.bounds is None


def test_a_line_with_a_transform_is_rejected():
    model = _streamlines()
    model.lines[1].transform.set_translation(1.0, 2.0, 3.0)

    with pytest.raises(ValueError, match="non-identity transform"):
        model.to_geojson()


def test_a_field_with_the_wrong_row_count_is_rejected():
    model = _slice()
    model.fields[1].values = np.ones((3, 1))
    with pytest.raises(ValueError, match="speed"):
        model.to_geojson()

    model = _streamlines()
    model.lines[1].fields[0].values = np.ones((2, 3))
    with pytest.raises(ValueError, match="velocity"):
        model.to_geojson()


def test_a_geometry_field_is_one_value_per_line_and_slice_metadata():
    model = _streamlines()
    model.lines[0].fields.append(
        Field(
            name="length",
            unit="m",
            values=np.array([12.5]),
            association="geometry",
        )
    )
    model.lines[1].fields.append(
        Field(
            name="drift",
            unit="m/s",
            values=np.array([[1.0, 2.0, 3.0]]),
            dim=3,
            association="geometry",
        )
    )

    [first, second] = model.to_geojson(crs="EPSG:4326")["features"]

    assert first["properties"]["length"] == 12.5
    assert first["properties"]["length_unit"] == "m"
    assert second["properties"]["drift"] == [1.0, 2.0, 3.0]
    assert second["properties"]["drift_components"] == ["u", "v", "w"]

    model = _slice()
    model.fields.append(
        Field(
            name="mean",
            unit="m/s",
            values=np.array([[2.5]]),
            association="geometry",
        )
    )
    model.fields.append(
        Field(
            name="drift",
            unit="m/s",
            values=np.array([[1.0, 2.0, 3.0]]),
            dim=3,
            association="geometry",
        )
    )

    payload = model.to_geojson()

    assert payload["metadata"]["mean"] == 2.5
    assert payload["metadata"]["mean_unit"] == "m/s"
    assert payload["metadata"]["drift"] == [1.0, 2.0, 3.0]
    assert payload["metadata"]["drift_components"] == ["u", "v", "w"]
    assert not {"mean", "mean_unit", "drift"} & set(
        payload["features"][0]["properties"]
    )


@PRODUCTS
def test_the_payload_frame_wins_over_the_metadata_payload(make):
    model = make(metadata_payload={"bounds": [0, 0, 0, 1, 1, 1], "crs": "EPSG:1"})

    for crs in (None, "EPSG:4326"):
        metadata = model.to_geojson(crs=crs)["metadata"]
        assert metadata["bounds"][2] == 2.0 and metadata["bounds"][5] == 40.0
        assert metadata["crs"] == (crs or "EPSG:3006")


def test_streamlines_follow_their_include_z_by_default():
    flat = _streamlines(include_z=False)

    assert _positions(flat.to_geojson()).shape[1] == 2
    assert _positions(flat.to_geojson(include_z=True)).shape[1] == 3
    assert flat.to_geojson()["metadata"]["include_z"] is False
    assert _positions(_streamlines().to_geojson()).shape[1] == 3


@PRODUCTS
def test_the_source_crs_comes_from_a_declaration_or_a_single_context_crs(make):
    model = make()
    model.crs = None

    with pytest.raises(ValueError, match="Cannot resolve the source CRS"):
        model.to_geojson()
    with pytest.raises(ValueError, match="Cannot resolve the source CRS"):
        model.to_geojson(crs="EPSG:4326")

    attach_dataset_context(model, _context(["EPSG:3006"]))
    assert model.to_geojson()["crs"] == NATIVE_CRS
    assert model.to_geojson(crs="EPSG:4326")["features"][0]["geometry"]

    attach_dataset_context(model, _context(["EPSG:3006", "EPSG:4326"]))
    with pytest.raises(ValueError, match="Cannot resolve the source CRS"):
        model.to_geojson()


def _smoke(product, **kwargs):
    return datasets.smoke(
        bounds=BOUNDS, product=product, resolution=4, streamline_count=4,
        streamline_steps=20, **kwargs,
    )  # fmt: skip


def _payload(package):
    [artifact] = [a for a in package.artifacts if a.format == "geojson"]
    if package.package_format == "dtccpkg":
        with zipfile.ZipFile(package.path) as archive:
            return json.loads(archive.read(artifact.path))
    return json.loads((package.path / artifact.path).read_text())


@pytest.mark.usefixtures("no_network")
@pytest.mark.parametrize("target", ["pkg", "pkg.dtccpkg"])
@pytest.mark.parametrize("product", ["slice", "streamlines"])
def test_canonical_geojson_round_trips(tmp_path, product, target):
    result = _smoke(product)

    package = result.export(tmp_path / target, canonical=True, format="geojson")
    restored = load_model_package(package.path)
    again = restored.export(
        tmp_path / f"again_{target}", canonical=True, format="geojson"
    )

    payload = _payload(package)
    assert payload == _payload(again)
    assert "crs" not in payload
    assert payload["metadata"]["crs"] == "EPSG:4326"
    assert "time" in payload["metadata"]
    for copy in (payload, _payload(again)):
        properties = copy["features"][0]["properties"]
        for name in ("velocity", "speed", "pressure"):
            assert properties[f"{name}_unit"] == UNITS[name]
    stored = result.points[0] if product == "slice" else result.lines[0].vertices[0]
    expected = _lonlat(stored[0], stored[1])
    positions = _positions(payload)
    assert positions[0, :2] == pytest.approx(expected, abs=1e-9)
    assert positions[0, 2] == pytest.approx(stored[2])
    if product == "slice":
        assert stored[:2].tolist() == [320000.0, 6399000.0]
        assert positions[0, :2] == pytest.approx(LONLAT, abs=1e-9)

    [artifact] = [a for a in package.artifacts if a.format == "geojson"]
    x, y, z = positions.T
    assert artifact.crs == "EPSG:4326"
    assert artifact.bounds == pytest.approx(
        [x.min(), y.min(), z.min(), x.max(), y.max(), z.max()]
    )
    xmin, ymin, zmin, xmax, ymax, zmax = payload["metadata"]["bounds"]
    assert xmin <= x.min() and x.max() <= xmax
    assert ymin <= y.min() and y.max() <= ymax
    assert zmin <= z.min() and z.max() <= zmax


@pytest.mark.usefixtures("no_network")
@pytest.mark.parametrize("target", ["pkg", "pkg.dtccpkg"])
def test_canonical_streamlines_keep_geometry_fields_and_include_z(tmp_path, target):
    result = _smoke("streamlines")
    result.include_z = False
    result.lines[0].fields.append(
        Field(name="length", unit="m", values=np.array([12.5]), association="geometry")
    )

    package = result.export(tmp_path / target, canonical=True, format="geojson")
    restored = load_model_package(package.path)
    again = restored.export(
        tmp_path / f"again_{target}", canonical=True, format="geojson"
    )

    payload = _payload(package)
    assert payload == _payload(again)
    assert payload["metadata"]["include_z"] is False
    assert _positions(payload).shape[1] == 2
    assert payload["features"][0]["properties"]["length"] == 12.5
    assert payload["features"][0]["properties"]["length_unit"] == "m"


@pytest.mark.usefixtures("no_network")
def test_canonical_slice_keeps_geometry_fields_in_metadata(tmp_path):
    result = _smoke("slice")
    result.fields.append(
        Field(name="mean", unit="m/s", values=np.array([[2.5]]), association="geometry")
    )

    package = result.export(tmp_path / "pkg", canonical=True, format="geojson")
    again = load_model_package(package.path).export(
        tmp_path / "again", canonical=True, format="geojson"
    )

    payload = _payload(package)
    assert payload == _payload(again)
    assert payload["metadata"]["mean"] == 2.5
    assert payload["metadata"]["mean_unit"] == "m/s"


@pytest.mark.parametrize("crs", [None, "EPSG:4326"])
def test_a_vertical_slice_exports_real_easting_and_northing(crs):
    result = _smoke("slice", slice_axis="x", slice_position=0.5)

    positions = _positions(result.to_geojson(crs=crs))

    assert len(positions) == len(result.points)
    stored = result.points
    assert np.ptp(stored[:, 0]) == 0.0 and np.ptp(stored[:, 1]) > 0.0
    if crs is None:
        np.testing.assert_array_equal(positions, stored)
    else:
        expected = np.array([_lonlat(x, y) for x, y in stored[:, :2]])
        np.testing.assert_allclose(positions[:, :2], expected, atol=1e-9)
        assert np.ptp(positions[:, 1]) > 0.0


def _integer_field(name, values, association, dim=1):
    return Field(
        name=name,
        unit="m",
        values=np.asarray(values, dtype=np.int64).reshape(-1, dim),
        dim=dim,
        association=association,
    )


def _with_integer(product, value, kind):
    if kind == "geometry":
        field = _integer_field("count", [value], "geometry")
        target = product if isinstance(product, FieldSlice) else product.lines[0]
    elif kind == "vector":
        rows = np.zeros((4, 3), dtype=np.int64)
        rows[2, 1] = value
        field = _integer_field("flow", rows, "sample", dim=3)
        target = product
    else:
        association = (
            "vertex" if isinstance(product, StreamlineCollection) else "sample"
        )
        count = 3 if association == "vertex" else 4
        values = np.zeros(count, dtype=np.int64)
        values[-1] = value
        field = _integer_field("count", values, association)
        target = product if isinstance(product, FieldSlice) else product.lines[0]
    target.fields.append(field)
    return product


@pytest.mark.parametrize("kind", ["scalar", "vector", "geometry"])
def test_a_slice_field_float64_would_round_is_refused(kind):
    model = _with_integer(_slice(), 2**53 + 1, kind)

    with pytest.raises(ValueError, match=rf"Field '\w+'.*{2**53 + 1}"):
        model.to_geojson()


def test_a_slice_field_float64_holds_exactly_is_exported():
    properties = _with_integer(_slice(), 2**60, "scalar").to_geojson()["features"][-1][
        "properties"
    ]

    assert properties["count"] == float(2**60)


@pytest.mark.parametrize("kind", ["scalar", "geometry"])
def test_a_streamline_field_float64_would_round_is_refused(kind):
    model = _with_integer(_streamlines(), 2**53 + 1, kind)

    with pytest.raises(ValueError, match=rf"Field 'count'.*{2**53 + 1}"):
        model.to_geojson()


@pytest.mark.usefixtures("no_network")
def test_a_canonical_package_keeps_the_integer_float64_would_round(tmp_path):
    model = _with_integer(_slice(), 2**53 + 1, "scalar")
    attach_dataset_context(model, _context(["EPSG:3006"]))

    package = model.export(tmp_path / "pkg", canonical=True)
    restored = load_model_package(package.path)

    assert restored.field("count").values.dtype == np.int64
    assert int(restored.field("count").values[-1, 0]) == 2**53 + 1
    target = tmp_path / "geojson"
    with pytest.raises(ValueError, match=rf"Field 'count'.*{2**53 + 1}"):
        model.export(target, canonical=True, format="geojson")
    assert not target.exists()
