"""Tests for the GeoJSON web writers of RoadNetwork and DeSO."""

import json
import math

import numpy as np
import pytest

import dtcc_core.datasets as datasets
from dtcc_core.builder.roadnetwork.space_syntax import analyze_space_syntax
from dtcc_core.datasets import attach_dataset_context, load_model_package
from dtcc_core.io import load_roadnetwork
from dtcc_core.model import (
    DeSO,
    Field,
    GeometryType,
    LineString,
    MultiLineString,
    MultiSurface,
    Object,
    RoadNetwork,
    Surface,
)

# Oracles: pyproj EPSG:3006 -> EPSG:4326 (always_xy) in Gothenburg.
A = (319950.0, 6398950.0)
B = (320050.0, 6399000.0)
C = (320150.0, 6398950.0)
A_WGS84 = [11.9788460825461, 57.696633876500854]
B_WGS84 = [11.980484024390593, 57.697122380732424]
C_WGS84 = [11.982196664423405, 57.69671386186802]
SQUARE = [(320000, 6399000), (320100, 6399000), (320100, 6399100), (320000, 6399100)]
HOLE = [(320040, 6399040), (320060, 6399040), (320060, 6399060), (320040, 6399060)]
SQUARE_WGS84_FIRST = [11.979646369152894, 57.697102381277034]
HOLE_WGS84_FIRST = [11.980286613302345, 57.697477181573014]
SQUARE_WGS84_THIRD = [11.981247004300723, 57.698039376141686]


def _road(curved=True, fields=(), srs="EPSG:3006"):
    road = RoadNetwork()
    road.vertices = np.array([A, C])
    road.edges = np.array([[0, 1]])
    road.length = np.array([223.6])
    road.attributes = {"highway": ["primary"], "id": [77]}
    road.transform.srs = srs
    if curved:
        line = MultiLineString(linestrings=[LineString(vertices=np.array([A, B, C]))])
        line.fields = list(fields)
        road.add_geometry(line, GeometryType.MULTILINESTRING)
    return road


def _flat(coordinates):
    return [value for position in coordinates for value in position]


def _context(obj, name, crs="EPSG:3006"):
    dataset = getattr(datasets, name)
    context = dataset.create_context(dataset.validate({"bounds": (0.0, 0.0, 1.0, 1.0)}))
    context.metadata.crs = [crs]
    return attach_dataset_context(obj, context)


def _deso(hole=True, nan=False):
    deso = DeSO()
    deso.transform.srs = "EPSG:3006"
    for code, offset in (("1480A0001", 0.0), ("1480A0002", 200.0)):
        area = Object(attributes={"desokod": code})
        area.id = code
        surface = Surface(vertices=np.array([[x + offset, y, 0.0] for x, y in SQUARE]))
        if hole and offset == 0.0:
            surface.holes = [np.array([[x, y, 0.0] for x, y in HOLE])]
        area.add_geometry(MultiSurface(surfaces=[surface]), GeometryType.LOD0)
        deso.add_child(area)
    deso.attach_field(
        Field(
            name="population",
            unit="persons",
            values=np.array([1200.0, math.nan if nan else 800.0]),
            association="sample",
        )
    )
    return deso


def _round_trip(obj, tmp_path):
    package = obj.export(tmp_path / "p", canonical=True, format="geojson")
    restored = load_model_package(package.path)
    again = restored.export(tmp_path / "q", canonical=True, format="geojson")
    return [
        json.loads((p.path / a.path).read_text())
        for p in (package, again)
        for a in p.artifacts
        if a.format == "geojson"
    ]


def test_road_exports_full_line_when_representation_exists():
    payload = _road().to_geojson(crs="EPSG:4326")

    assert payload["metadata"] == {"geometry": "full"}
    [feature] = payload["features"]
    assert feature["geometry"]["type"] == "LineString"
    assert len(feature["geometry"]["coordinates"]) == 3
    assert _flat(feature["geometry"]["coordinates"]) == pytest.approx(
        _flat([A_WGS84, B_WGS84, C_WGS84]), abs=1e-9
    )
    assert feature["id"] == 77
    assert feature["properties"]["highway"] == "primary"


def test_road_exports_chord_without_representation():
    payload = _road(curved=False).to_geojson(crs="EPSG:4326")

    assert payload["metadata"] == {"geometry": "chord"}
    [feature] = payload["features"]
    assert len(feature["geometry"]["coordinates"]) == 2
    assert _flat(feature["geometry"]["coordinates"]) == pytest.approx(
        _flat([A_WGS84, C_WGS84]), abs=1e-9
    )


def test_road_chord_when_representation_is_not_aligned_with_edges():
    road = _road()
    road.vertices = np.array([A, B, C])
    road.edges = np.array([[0, 1], [1, 2]])
    road.length = np.array([111.8, 111.8])
    road.attributes = {}

    payload = road.to_geojson()

    assert payload["metadata"] == {"geometry": "chord"}
    assert [f["id"] for f in payload["features"]] == [0, 1]
    assert [len(f["geometry"]["coordinates"]) for f in payload["features"]] == [2, 2]


def test_road_native_export_names_crs_and_keeps_coordinates():
    payload = _road().to_geojson()

    assert payload["crs"] == {"type": "name", "properties": {"name": "EPSG:3006"}}
    assert payload["features"][0]["geometry"]["coordinates"] == [
        list(A),
        list(B),
        list(C),
    ]


def test_road_root_transform_applied_once_without_mutation():
    road = _road()
    road.transform.affine = np.array(
        [[1, 0, 0, 100], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]], dtype=float
    )

    payload = road.to_geojson()

    assert payload["features"][0]["geometry"]["coordinates"][0] == [A[0] + 100, A[1]]
    assert road.linestrings[0].vertices[0].tolist() == list(A)


def test_road_length_is_stored_value_with_source_crs_unit():
    metric = _road().to_geojson(crs="EPSG:4326")["features"][0]["properties"]
    degree = _road(srs="EPSG:4326")
    degree.length = np.array([0.0021])
    degree = degree.to_geojson(crs="EPSG:4326")["features"][0]["properties"]

    assert metric["length"] == 223.6
    assert metric["length_unit"] == "metre"
    assert degree["length"] == 0.0021
    assert degree["length_unit"] == "degree"


def test_loaded_networks_export_length_unit_of_their_crs(tmp_path):
    source = tmp_path / "roads.geojson"
    source.write_text(
        json.dumps(
            {
                "type": "FeatureCollection",
                "features": [
                    {
                        "type": "Feature",
                        "properties": {"id": "way-9", "highway": "service"},
                        "geometry": {
                            "type": "LineString",
                            "coordinates": [
                                [11.0, 57.0],
                                [11.001, 57.001],
                                [11.002, 57],
                            ],
                        },
                    }
                ],
            }
        )
    )

    native = load_roadnetwork(source, target_crs="EPSG:4326")
    projected = load_roadnetwork(source, target_crs="EPSG:3006")

    degree = native.to_geojson()["features"][0]
    metre = projected.to_geojson(crs="EPSG:4326")["features"][0]
    assert degree["id"] == "way-9"
    assert degree["properties"]["length"] == float(native.length[0])
    assert degree["properties"]["length_unit"] == "degree"
    assert len(degree["geometry"]["coordinates"]) == 3
    assert metre["properties"]["length"] == float(projected.length[0])
    assert metre["properties"]["length_unit"] == "metre"
    assert _flat(metre["geometry"]["coordinates"]) == pytest.approx(
        _flat(degree["geometry"]["coordinates"]), abs=1e-6
    )


def test_road_carrier_field_is_exported_per_edge_with_unit():
    speed = Field(
        name="speed",
        unit="km/h",
        values=np.array([[50.0]]),
        association="sample",
    )

    properties = _road(fields=[speed]).to_geojson()["features"][0]["properties"]

    assert properties["speed"] == 50.0
    assert properties["speed_unit"] == "km/h"


def test_road_carrier_field_must_have_one_row_per_edge():
    speed = Field(name="speed", unit="km/h", values=np.array([1.0, 2.0]))
    speed.association = "sample"

    with pytest.raises(ValueError, match="speed"):
        _road(fields=[speed]).to_geojson()


def test_road_crs_declared_only_on_a_line_is_used():
    road = _road(srs="")
    road.linestrings[0].transform.srs = "EPSG:3006"

    payload = road.to_geojson(crs="EPSG:4326")

    assert _flat(payload["features"][0]["geometry"]["coordinates"]) == pytest.approx(
        _flat([A_WGS84, B_WGS84, C_WGS84]), abs=1e-9
    )
    assert road.to_geojson()["crs"]["properties"]["name"] == "EPSG:3006"


def test_road_conflicting_line_crs_raises():
    road = _road()
    road.linestrings[0].transform.srs = "EPSG:4326"

    with pytest.raises(ValueError, match="Conflicting CRS"):
        road.to_geojson(crs="EPSG:4326")


@pytest.mark.parametrize("crs", [None, "EPSG:4326"])
def test_empty_road_network_with_field_exports_no_features(crs):
    road = RoadNetwork()
    road.transform.srs = "EPSG:3006"
    line = MultiLineString()
    line.fields = [
        Field(name="speed", unit="km/h", values=np.empty((0, 1)), association="sample")
    ]
    road.add_geometry(line, GeometryType.MULTILINESTRING)

    payload = road.to_geojson(crs=crs)

    assert payload["features"] == []
    assert payload["metadata"] == {"geometry": "chord"}


def test_road_length_must_have_one_value_per_edge():
    road = _road(curved=False)
    road.vertices = np.array([A, B, C])
    road.edges = np.array([[0, 1], [1, 2]])
    road.length = np.array([1.0])
    road.attributes = {}

    with pytest.raises(ValueError, match="length"):
        road.to_geojson()


def test_road_without_stored_lengths_exports_no_length_properties():
    road = _road(curved=False)
    road.vertices = np.array([A, B, C])
    road.edges = np.array([[0, 1], [1, 2]])
    road.length = np.empty(0)
    road.attributes = {}

    features = road.to_geojson()["features"]

    assert len(features) == 2
    for feature in features:
        assert "length" not in feature["properties"]
        assert "length_unit" not in feature["properties"]


def test_road_space_syntax_units():
    road = _road(curved=False)
    road.vertices = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]])
    road.edges = np.array([[0, 1], [1, 2]])
    road.length = np.array([1.0, 1.0])
    road.attributes = {}

    for cost, unit in (("topological", "steps"), ("metric", "meters")):
        result = analyze_space_syntax(road, cost=cost)
        result.transform.srs = "EPSG:3006"

        properties = result.to_geojson()["features"][1]["properties"]

        assert properties["space_syntax_connectivity"] == 1.0
        assert properties["space_syntax_mean_depth"] == 1.0
        assert properties["space_syntax_mean_depth_unit"] == unit
        assert properties["space_syntax_integration_unit"] == f"1/{unit}"
        for dimensionless in ("connectivity", "reach", "choice"):
            assert f"space_syntax_{dimensionless}" in properties
            assert f"space_syntax_{dimensionless}_unit" not in properties


def test_road_space_syntax_units_only_for_present_measures():
    result = analyze_space_syntax(_road(curved=False), measures=("connectivity",))
    result.transform.srs = "EPSG:3006"

    properties = result.to_geojson()["features"][0]["properties"]

    assert "space_syntax_mean_depth_unit" not in properties
    assert "space_syntax_integration_unit" not in properties


def test_road_canonical_round_trip(tmp_path):
    speed = Field(
        name="speed", unit="km/h", values=np.array([[50.0]]), association="sample"
    )
    road = _context(_road(fields=[speed]), "roads")

    payload, again = _round_trip(road, tmp_path)

    assert payload == again
    assert "crs" not in payload
    assert payload["metadata"] == {"geometry": "full"}
    [feature] = payload["features"]
    assert feature["geometry"]["coordinates"][1] == pytest.approx(B_WGS84, abs=1e-9)
    assert feature["id"] == 77
    assert feature["properties"]["speed"] == 50.0
    assert feature["properties"]["speed_unit"] == "km/h"
    assert feature["properties"]["length_unit"] == "metre"


def test_deso_keeps_hole_after_reprojection():
    payload = _deso().to_geojson(crs="EPSG:4326")

    first, second = payload["features"]
    assert first["geometry"]["type"] == "Polygon"
    rings = first["geometry"]["coordinates"]
    assert len(rings) == 2
    assert rings[0][0] == pytest.approx(SQUARE_WGS84_FIRST, abs=1e-9)
    assert rings[0][2] == pytest.approx(SQUARE_WGS84_THIRD, abs=1e-9)
    assert rings[1][0] == pytest.approx(HOLE_WGS84_FIRST, abs=1e-9)
    assert len(second["geometry"]["coordinates"]) == 1


def test_deso_crs_declared_only_on_a_surface_is_used():
    deso = _deso()
    deso.transform.srs = ""
    for area in deso.areas:
        for surface in area.get_geometry(GeometryType.LOD0).surfaces:
            surface.transform.srs = "EPSG:3006"

    payload = deso.to_geojson(crs="EPSG:4326")

    ring = payload["features"][0]["geometry"]["coordinates"][0]
    assert ring[0] == pytest.approx(SQUARE_WGS84_FIRST, abs=1e-9)
    assert deso.to_geojson()["crs"]["properties"]["name"] == "EPSG:3006"


def test_deso_crs_declared_only_on_an_area_is_used():
    deso = _deso()
    deso.transform.srs = ""
    for area in deso.areas:
        area.transform.srs = "EPSG:3006"

    payload = deso.to_geojson(crs="EPSG:4326")

    ring = payload["features"][0]["geometry"]["coordinates"][0]
    assert ring[0] == pytest.approx(SQUARE_WGS84_FIRST, abs=1e-9)
    assert deso.to_geojson()["crs"]["properties"]["name"] == "EPSG:3006"


def test_deso_conflicting_area_crs_raises():
    deso = _deso()
    deso.areas[0].transform.srs = "EPSG:4326"

    with pytest.raises(ValueError, match="Conflicting CRS"):
        deso.to_geojson(crs="EPSG:4326")


def test_deso_conflicting_surface_crs_raises():
    deso = _deso()
    deso.areas[0].get_geometry(GeometryType.LOD0).surfaces[
        0
    ].transform.srs = "EPSG:4326"

    with pytest.raises(ValueError, match="Conflicting CRS"):
        deso.to_geojson(crs="EPSG:4326")


def test_deso_multipart_area_is_multipolygon():
    deso = _deso()
    area = deso.areas[0]
    geometry = area.get_geometry(GeometryType.LOD0)
    geometry.surfaces.append(
        Surface(vertices=np.array([[x + 500, y, 0.0] for x, y in SQUARE]))
    )

    feature = deso.to_geojson()["features"][0]

    assert feature["geometry"]["type"] == "MultiPolygon"
    assert [len(polygon) for polygon in feature["geometry"]["coordinates"]] == [2, 1]


def test_deso_statistics_ids_units_and_nan():
    payload = _deso(nan=True).to_geojson()

    first, second = payload["features"]
    assert [first["id"], second["id"]] == ["1480A0001", "1480A0002"]
    assert first["properties"]["population"] == 1200.0
    assert first["properties"]["population_unit"] == "persons"
    assert second["properties"]["population"] is None
    assert second["properties"]["population_unit"] == "persons"
    assert first["properties"]["desokod"] == "1480A0001"
    json.dumps(payload, allow_nan=False)


def test_deso_native_export_names_crs():
    payload = _deso().to_geojson()

    assert payload["crs"] == {"type": "name", "properties": {"name": "EPSG:3006"}}
    ring = payload["features"][0]["geometry"]["coordinates"][0]
    assert ring[0] == [320000.0, 6399000.0]
    assert ring[0] == ring[-1]


def test_deso_canonical_round_trip(tmp_path):
    deso = _context(_deso(nan=True), "deso")

    payload, again = _round_trip(deso, tmp_path)

    assert payload == again
    assert "crs" not in payload
    first, second = payload["features"]
    assert len(first["geometry"]["coordinates"]) == 2
    assert first["geometry"]["coordinates"][0][0] == pytest.approx(
        SQUARE_WGS84_FIRST, abs=1e-9
    )
    assert second["properties"]["population"] is None
    assert second["properties"]["population_unit"] == "persons"


def test_deso_area_field_float64_would_round_is_refused():
    deso = _deso()
    geometry = deso.areas[1].get_geometry(GeometryType.LOD0)
    geometry.fields.append(
        Field(
            name="households",
            unit="count",
            values=np.array([2**53 + 1], dtype=np.int64),
            association="geometry",
        )
    )

    with pytest.raises(ValueError, match=rf"Field 'households'.*{2**53 + 1}"):
        deso.to_geojson()


def test_deso_area_without_the_statistic_still_exports_null():
    deso = _deso()
    deso.areas[1].get_geometry(GeometryType.LOD0).fields.clear()

    first, second = deso.to_geojson()["features"]

    assert first["properties"]["population"] == 1200.0
    assert second["properties"]["population"] is None
