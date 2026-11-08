"""Tests for the GeoJSON web profile helpers in dtcc_core.io.geojson."""

import json
import math
import re

import numpy as np
import pytest

import dtcc_core.datasets as datasets
from dtcc_core.datasets import attach_dataset_context
from dtcc_core.io.geojson import (
    exported_coordinates,
    feature_collection,
    polygon_rings,
    positions,
    resolve_source_crs,
)
from dtcc_core.model import (
    CalibrationGrid,
    Field,
    Mesh,
    MultiLineString,
    Object,
    Point,
    RoadNetwork,
    SensorCollection,
    Transform,
    VehicleCollection,
)

# Oracle: pyproj EPSG:3006 -> EPSG:4326 for a point in Gothenburg.
SWEREF = (319950.0, 6398950.0)
WGS84 = (11.9788460825461, 57.696633876500854)


def _signed_area(ring):
    return sum(a[0] * b[1] - b[0] * a[1] for a, b in zip(ring, ring[1:], strict=False))


def _station(station_id, x, y, fields=()):
    station = Object(attributes={"station_id": station_id})
    point = Point(x=x, y=y, z=10.0)
    point.fields = list(fields)
    station.add_geometry(point, "location")
    return station


def _field(name, value, unit):
    return Field(
        name=name,
        unit=unit,
        association="sample",
        values=np.array([value], dtype=np.float32),
    )


def _context(obj, crs):
    context = datasets.weather.create_context(
        datasets.weather.validate({"bounds": (0.0, 0.0, 1.0, 1.0)})
    )
    context.metadata.crs = list(crs)
    attach_dataset_context(obj, context)
    return obj


def test_positions_reproject_sweref_to_wgs84_within_1e9():
    [position] = positions(np.array([SWEREF]), "EPSG:3006", "EPSG:4326")

    assert len(position) == 2  # 2D input is accepted and stays 2D
    assert position[0] == pytest.approx(WGS84[0], abs=1e-9)
    assert position[1] == pytest.approx(WGS84[1], abs=1e-9)


def test_positions_keep_z_and_native_coordinates():
    [position] = positions([[*SWEREF, 12.5]], "EPSG:3006", None)

    assert position == [*SWEREF, 12.5]
    assert all(type(value) is float for value in position)


def test_positions_reject_non_finite_coordinate():
    with pytest.raises(ValueError, match="finite"):
        positions([[SWEREF[0], math.nan]], "EPSG:3006", "EPSG:4326")


def test_polygon_rings_close_and_orient_exterior_ccw_hole_cw():
    exterior_cw = [[0, 0], [0, 10], [10, 10], [10, 0]]
    hole_ccw = [[2, 2], [4, 2], [4, 4], [2, 4]]

    rings = polygon_rings(exterior_cw, [hole_ccw])

    assert len(rings) == 2
    for ring in rings:
        assert ring[0] == ring[-1]
        assert len(ring) == 5
    assert _signed_area(rings[0]) > 0
    assert _signed_area(rings[1]) < 0


def test_feature_collection_properties_are_json_safe():
    feature = {
        "type": "Feature",
        "geometry": {"type": "Point", "coordinates": list(SWEREF)},
        "properties": {
            "missing": math.nan,
            "infinite": np.float64(np.inf),
            "count": np.int64(3),
            "values": np.array([1.5, np.nan]),
        },
    }

    payload = feature_collection([feature], source_crs="EPSG:3006", crs=None)

    properties = payload["features"][0]["properties"]
    assert properties == {
        "missing": None,
        "infinite": None,
        "count": 3,
        "values": [1.5, None],
    }
    assert type(properties["count"]) is int
    assert payload["crs"] == {"type": "name", "properties": {"name": "EPSG:3006"}}


def test_feature_collection_wgs84_has_no_crs_member():
    payload = feature_collection([], source_crs="EPSG:3006", crs="EPSG:4326")

    assert payload == {"type": "FeatureCollection", "features": []}


def test_sensor_collection_projected_mode_keeps_native_coordinates():
    collection = SensorCollection(attributes={"crs": "EPSG:3006"})
    collection.add_station(_station("s1", *SWEREF))

    payload = collection.to_geojson()

    assert payload["crs"]["properties"]["name"] == "EPSG:3006"
    feature = payload["features"][0]
    assert feature["geometry"]["coordinates"] == [*SWEREF, 10.0]

    web = collection.to_geojson(crs="EPSG:4326")
    assert "crs" not in web
    lon, lat, z = web["features"][0]["geometry"]["coordinates"]
    assert (lon, lat) == pytest.approx(WGS84, abs=1e-9)
    assert z == 10.0


def test_sensor_collection_properties_carry_units_and_null_for_missing():
    collection = SensorCollection(attributes={"crs": "EPSG:3006"})
    collection.add_station(
        _station(
            "s1",
            *SWEREF,
            [_field("air_temperature", -3.5, "celsius"), _field("wind", 4.0, "m/s")],
        )
    )
    collection.add_station(
        _station("s2", *SWEREF, [_field("air_temperature", math.nan, "celsius")])
    )

    features = collection.to_geojson(crs="EPSG:4326")["features"]

    assert [f["id"] for f in features] == ["s1", "s2"]
    first, second = (f["properties"] for f in features)
    assert first["air_temperature"] == -3.5
    assert first["air_temperature_unit"] == "celsius"
    assert first["wind"] == 4.0
    assert first["wind_unit"] == "m/s"
    assert second["station_id"] == "s2"
    assert second["air_temperature"] is None
    assert second["wind"] is None
    assert second["wind_unit"] == "m/s"


def test_empty_sensor_collection_exports_empty_feature_collection():
    collection = SensorCollection(attributes={"crs": "EPSG:3006"})

    assert collection.to_geojson(crs="EPSG:4326") == {
        "type": "FeatureCollection",
        "features": [],
    }


def test_resolve_crs_from_calibration_grid_crs():
    assert resolve_source_crs(CalibrationGrid(crs="EPSG:3006"), []) == "EPSG:3006"


def test_resolve_crs_from_sensor_attributes_with_blank_point_srs():
    collection = SensorCollection(attributes={"crs": "EPSG:3006"})
    collection.add_station(_station("s1", *SWEREF))
    [point] = collection.stations()[0].get_geometries()
    assert point.transform.srs == ""
    assert np.array_equal(point.transform.affine, np.eye(4))

    assert resolve_source_crs(collection, [point]) == "EPSG:3006"


def test_resolve_crs_from_road_network_representation_srs():
    road = RoadNetwork()
    lines = MultiLineString(transform=Transform(srs="EPSG:3006"))
    road.add_geometry(lines, id="lines")

    assert resolve_source_crs(road, road.get_geometries()) == "EPSG:3006"


def test_resolve_crs_aliases_agree():
    collection = SensorCollection(attributes={"crs": "epsg:3006"})
    point = Point(transform=Transform(srs="EPSG:3006"))

    assert resolve_source_crs(collection, [point]) == "epsg:3006"


def test_resolve_crs_conflicting_geometries_raise():
    first = Point(transform=Transform(srs="EPSG:3006"))
    second = Point(transform=Transform(srs="EPSG:4326"))

    with pytest.raises(ValueError, match="'EPSG:3006' and 'EPSG:4326'"):
        resolve_source_crs(SensorCollection(), [first, second])


def test_resolve_crs_context_is_ignored_when_a_declaration_exists():
    collection = _context(
        SensorCollection(attributes={"crs": "EPSG:3006"}), ["EPSG:4326"]
    )

    assert resolve_source_crs(collection, []) == "EPSG:3006"


@pytest.mark.parametrize("crs", [["EPSG:3006", "EPSG:4326"], []])
def test_resolve_crs_context_must_hold_exactly_one_value(crs):
    collection = _context(SensorCollection(), crs)

    with pytest.raises(ValueError, match=f"SensorCollection.*{re.escape(repr(crs))}"):
        resolve_source_crs(collection, [])


def test_resolve_crs_without_declaration_or_context_raises():
    with pytest.raises(ValueError, match="SensorCollection.*\\[\\]"):
        resolve_source_crs(SensorCollection(), [])


def test_root_translation_frames_export_the_same_positions():
    vertices = np.array(
        [
            [319950.0, 6398950.0, 0.0],
            [320050.0, 6398950.0, 0.0],
            [319950.0, 6399050.0, 5.0],
        ]
    )
    t = np.array([100.0, 100.0, 0.0])
    faces = np.array([[0, 1, 2]])
    reference = Mesh(vertices=vertices.copy(), faces=faces)
    reference.transform.srs = "EPSG:3006"
    shifted = Mesh(vertices=vertices - t, faces=faces)
    shifted.transform.srs = "EPSG:3006"
    shifted.transform.set_translation(*t)
    stored_reference, stored_shifted = (
        reference.vertices.copy(),
        shifted.vertices.copy(),
    )

    expected = positions(
        exported_coordinates(reference).vertices, "EPSG:3006", "EPSG:4326"
    )
    exported = exported_coordinates(shifted)
    actual = positions(exported.vertices, "EPSG:3006", "EPSG:4326")

    np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-9)
    np.testing.assert_array_equal(reference.vertices, stored_reference)
    np.testing.assert_array_equal(shifted.vertices, stored_shifted)
    assert not np.array_equal(shifted.transform.affine, np.eye(4))
    assert np.array_equal(exported.transform.affine, np.eye(4))
    assert exported.transform.srs == "EPSG:3006"


def test_non_identity_transform_below_the_root_is_rejected():
    obj = Object()
    point = Point(x=1.0, y=2.0)
    point.transform.set_translation(10.0, 0.0, 0.0)
    obj.add_geometry(point, "location")

    with pytest.raises(ValueError, match="Point below the root has a non-identity"):
        exported_coordinates(obj)
    assert point.x == 1.0


def test_shared_point_is_transformed_once():
    obj = Object()
    obj.transform.set_translation(10.0, 0.0, 0.0)
    point = Point(x=1.0, y=2.0)
    obj.add_geometry(point, id="first")
    obj.add_geometry(point, id="second")

    exported = exported_coordinates(obj)

    [first, second] = exported.get_geometries()
    assert first is second
    assert first.x == 11.0
    assert point.x == 1.0


def test_sensor_float32_field_and_numpy_id_export_exactly():
    collection = SensorCollection(attributes={"crs": "EPSG:3006"})
    collection.add_station(_station(np.int64(7), *SWEREF, [_field("wind", 7.3, "m/s")]))
    collection.add_station(_station("s2", *SWEREF, [_field("wind", math.nan, "m/s")]))

    payload = collection.to_geojson()

    first, second = payload["features"]
    assert first["properties"]["wind"] == 7.3
    assert second["properties"]["wind"] is None
    assert first["id"] == 7
    assert type(first["id"]) is int
    json.dumps(payload)


def test_sensor_field_units_are_per_station_with_first_seen_fallback():
    collection = SensorCollection(attributes={"crs": "EPSG:3006"})
    collection.add_station(_station("s1", *SWEREF, [_field("temperature", 20.0, "C")]))
    collection.add_station(_station("s2", *SWEREF, [_field("temperature", 68.0, "F")]))
    collection.add_station(_station("s3", *SWEREF))

    features = collection.to_geojson()["features"]

    first, second, third = (f["properties"] for f in features)
    assert first["temperature"] == 20.0
    assert first["temperature_unit"] == "C"
    assert second["temperature"] == 68.0
    assert second["temperature_unit"] == "F"
    assert third["temperature"] is None
    assert third["temperature_unit"] == "C"


def test_sensor_numpy_array_attributes_export_like_scalars():
    station = _station("s1", *SWEREF)
    station.attributes["temps"] = np.array([7.3, 1.1], dtype=np.float32)
    station.attributes["exact"] = np.array([0.1, 1 / 3])
    station.attributes["scalar"] = np.array(7.3, dtype=np.float32)
    collection = SensorCollection(attributes={"crs": "EPSG:3006"})
    collection.add_station(station)

    payload = collection.to_geojson()

    properties = payload["features"][0]["properties"]
    assert properties["temps"] == [7.3, 1.1]
    assert properties["exact"] == [0.1, 1 / 3]
    assert properties["scalar"] == 7.3
    json.dumps(payload)


def _vehicle(vehicle_id, x, y):
    vehicle = Object(attributes={"vehicle_id": vehicle_id})
    vehicle.add_geometry(Point(x=x, y=y, z=10.0), "position")
    return vehicle


def _point_collection(kind, member_srs=None, crs=None, extra_member_srs=None):
    attributes = {} if crs is None else {"crs": crs}
    if kind == "sensor":
        collection = SensorCollection(attributes=attributes)
        member = _station("s1", *SWEREF)
        collection.add_station(member)
    else:
        collection = VehicleCollection(attributes=attributes)
        member = _vehicle("v1", *SWEREF)
        collection.add_vehicle(member)
    if member_srs is not None:
        member.transform.srs = member_srs
    if extra_member_srs is not None:
        if kind == "sensor":
            extra = _station("s2", *SWEREF)
            collection.add_station(extra)
        else:
            extra = _vehicle("v2", *SWEREF)
            collection.add_vehicle(extra)
        extra.transform.srs = extra_member_srs
    return collection


@pytest.mark.parametrize("kind", ["sensor", "vehicle"])
def test_point_collection_resolves_crs_declared_only_on_members(kind):
    collection = _point_collection(kind, member_srs="EPSG:3006")

    payload = collection.to_geojson()
    web = collection.to_geojson(crs="EPSG:4326")

    assert payload["crs"]["properties"]["name"] == "EPSG:3006"
    assert payload["features"][0]["geometry"]["coordinates"] == [*SWEREF, 10.0]
    assert "crs" not in web
    lon, lat, z = web["features"][0]["geometry"]["coordinates"]
    assert (lon, lat) == pytest.approx(WGS84, abs=1e-9)
    assert z == 10.0


@pytest.mark.parametrize("kind", ["sensor", "vehicle"])
def test_point_collection_member_crs_conflicting_with_collection_raises(kind):
    collection = _point_collection(kind, member_srs="EPSG:4326", crs="EPSG:3006")

    with pytest.raises(ValueError, match="Conflicting CRS"):
        collection.to_geojson(crs="EPSG:4326")


@pytest.mark.parametrize("kind", ["sensor", "vehicle"])
def test_point_collection_member_crs_equal_to_collection_is_unchanged(kind):
    plain = _point_collection(kind, crs="EPSG:3006")
    declared = _point_collection(kind, member_srs="epsg:3006", crs="EPSG:3006")

    assert declared.to_geojson() == plain.to_geojson()
    assert declared.to_geojson()["crs"]["properties"]["name"] == "EPSG:3006"
    assert declared.to_geojson(crs="EPSG:4326") == plain.to_geojson(crs="EPSG:4326")


@pytest.mark.parametrize("kind", ["sensor", "vehicle"])
def test_point_collection_members_with_conflicting_crs_raise(kind):
    collection = _point_collection(
        kind, member_srs="EPSG:3006", extra_member_srs="EPSG:4326"
    )

    with pytest.raises(ValueError, match="Conflicting CRS declarations"):
        collection.to_geojson()
