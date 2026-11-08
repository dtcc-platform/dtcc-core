"""Tests for the GeoJSON web writers of the collection types and City."""

import copy
from unittest.mock import patch

import numpy as np
import pytest
from shapely.geometry import LinearRing

from dtcc_core.io.geojson import exported_coordinates
from dtcc_core.model import (
    Bounds,
    Building,
    BuildingCollection,
    BuildingPart,
    CalibrationGrid,
    City,
    FootprintCollection,
    GeometryType,
    MultiSurface,
    Surface,
    Tree,
    TreeCollection,
)

# Oracles: pyproj EPSG:3006 -> EPSG:4326 (always_xy) in Gothenburg.
SQUARE = [(320000, 6399000), (320100, 6399000), (320100, 6399100), (320000, 6399100)]
HOLE = [(320040, 6399040), (320060, 6399040), (320060, 6399060), (320040, 6399060)]
SQUARE_WGS84_FIRST = [11.979646369152894, 57.697102381277034]
HOLE_WGS84_FIRST = [11.980286613302345, 57.697477181573014]
SQUARE_WGS84_SECOND = [11.981321680736151, 57.69714237464992]
NATIVE_CRS = {"type": "name", "properties": {"name": "EPSG:3006"}}


def _surface(points=SQUARE, hole=None, srs="EPSG:3006"):
    surface = Surface(vertices=np.array([[x, y, 0.0] for x, y in points]))
    if hole is not None:
        surface.holes = [np.array([[x, y, 0.0] for x, y in hole])]
    surface.transform.srs = srs
    return surface


def _building(building_id="b1", surface=None, **attributes):
    building = Building(id=building_id, attributes=attributes)
    building.add_geometry(surface or _surface(hole=HOLE), GeometryType.LOD0)
    return building


def _assert_wgs84(payload):
    assert "crs" not in payload
    for feature in payload["features"]:
        geometry = feature["geometry"]
        coordinates = geometry["coordinates"]
        if geometry["type"] == "Polygon":
            coordinates = [position for ring in coordinates for position in ring]
        elif geometry["type"] == "Point":
            coordinates = [coordinates]
        for lon, lat, *_ in coordinates:
            assert -180 <= lon <= 180 and -90 <= lat <= 90


def test_footprint_collection_reprojects_and_keeps_hole():
    collection = FootprintCollection.from_buildings([_building("b1", name="x")])

    payload = collection.to_geojson(crs="EPSG:4326")

    _assert_wgs84(payload)
    [feature] = payload["features"]
    rings = feature["geometry"]["coordinates"]
    assert len(rings) == 2
    assert rings[0][0] == pytest.approx(SQUARE_WGS84_FIRST, abs=1e-9)
    assert rings[1][0] == pytest.approx(HOLE_WGS84_FIRST, abs=1e-9)
    assert feature["id"] == "b1"
    assert feature["properties"] == {
        "index": 0,
        "source_index": 0,
        "source_id": "b1",
    }


def test_footprint_collection_without_source_ids_has_no_feature_id():
    payload = FootprintCollection(footprints=[_surface()]).to_geojson(crs="EPSG:4326")

    [feature] = payload["features"]
    assert "id" not in feature


def test_footprint_collection_native_names_crs_and_keeps_coordinates():
    collection = FootprintCollection.from_buildings([_building()])

    payload = collection.to_geojson()

    assert payload["crs"] == NATIVE_CRS
    rings = payload["features"][0]["geometry"]["coordinates"]
    assert len(rings) == 2
    assert sorted(map(tuple, rings[0][:-1])) == sorted(SQUARE)


def test_building_collection_features_carry_ids_attributes_and_skips():
    keep = _building("b1", storeys=3, name="Hall")
    bare = Building(id="b2")
    collection = BuildingCollection([keep, bare])

    payload = collection.to_geojson(crs="EPSG:4326")

    _assert_wgs84(payload)
    assert payload["metadata"] == {"skipped": 1}
    [feature] = payload["features"]
    assert feature["id"] == "b1"
    assert feature["properties"] == {"storeys": 3, "name": "Hall"}
    rings = feature["geometry"]["coordinates"]
    assert len(rings) == 2
    assert rings[0][0] == pytest.approx(SQUARE_WGS84_FIRST, abs=1e-9)
    assert rings[1][0] == pytest.approx(HOLE_WGS84_FIRST, abs=1e-9)


def test_building_collection_native_names_crs():
    payload = BuildingCollection([_building()]).to_geojson()

    assert payload["crs"] == NATIVE_CRS
    assert payload["metadata"] == {"skipped": 0}
    assert payload["features"][0]["geometry"]["coordinates"][0][0] in (
        [float(x), float(y)] for x, y in SQUARE
    )


def test_city_writes_building_footprints_and_applies_root_transform():
    city = City()
    city.transform.srs = "EPSG:3006"
    city.add_building(_building("b1", surface=_surface(srs="")))
    city.transform.affine = np.array(
        [[1, 0, 0, 100], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]], dtype=float
    )

    payload = city.to_geojson(crs="EPSG:4326")

    _assert_wgs84(payload)
    [feature] = payload["features"]
    assert feature["id"] == "b1"
    # The root translation moves the first corner to (320100, 6399000).
    assert feature["geometry"]["coordinates"][0][0] == pytest.approx(
        SQUARE_WGS84_SECOND, abs=1e-9
    )
    assert city.transform.affine[0, 3] == 100
    assert city.buildings[0].get_geometry(GeometryType.LOD0).vertices[0, 0] == 320000


def test_city_native_names_crs():
    city = City()
    city.add_building(_building("b1"))

    payload = city.to_geojson()

    assert payload["crs"] == NATIVE_CRS
    assert payload["metadata"] == {"skipped": 0}


def test_tree_collection_points_carry_height_and_crown_radius_with_units():
    tree = Tree(
        id="t1",
        position=np.array([320000.0, 6399000.0, 12.0]),
        height=9.5,
        crown_radius=2.5,
    )
    tree.transform.srs = "EPSG:3006"

    payload = TreeCollection([tree]).to_geojson(crs="EPSG:4326")

    _assert_wgs84(payload)
    [feature] = payload["features"]
    assert feature["id"] == "t1"
    assert feature["geometry"]["type"] == "Point"
    assert feature["geometry"]["coordinates"][:2] == pytest.approx(
        SQUARE_WGS84_FIRST, abs=1e-9
    )
    assert feature["geometry"]["coordinates"][2] == 12.0
    assert feature["properties"] == {
        "height": 9.5,
        "height_unit": "m",
        "crown_radius": 2.5,
        "crown_radius_unit": "m",
    }


def test_tree_collection_native_names_crs():
    tree = Tree(id="t1", position=np.array([320000.0, 6399000.0, 12.0]))
    tree.transform.srs = "EPSG:3006"

    payload = TreeCollection([tree]).to_geojson()

    assert payload["crs"] == NATIVE_CRS
    assert payload["features"][0]["geometry"]["coordinates"] == [
        320000.0,
        6399000.0,
        12.0,
    ]


def _grid(crs="EPSG:3006"):
    line = {
        "type": "Feature",
        "geometry": {
            "type": "LineString",
            "coordinates": [[320000.0, 6399000.0], [320100.0, 6399000.0]],
        },
        "properties": {"axis": "x"},
    }
    return CalibrationGrid(features=[line], divisions=1, crs=crs)


def test_calibration_grid_reprojects_lines_and_keeps_metadata():
    grid = _grid()
    native = grid.to_geojson()

    payload = grid.to_geojson(crs="EPSG:4326")

    _assert_wgs84(payload)
    [feature] = payload["features"]
    assert feature["geometry"]["coordinates"][0] == pytest.approx(
        SQUARE_WGS84_FIRST, abs=1e-9
    )
    assert feature["geometry"]["coordinates"][1] == pytest.approx(
        SQUARE_WGS84_SECOND, abs=1e-9
    )
    assert feature["properties"] == {"axis": "x"}
    assert payload["name"] == "calibration_grid"
    frame = {"crs", "bounds", "spacing"}
    assert {k: v for k, v in payload["metadata"].items() if k not in frame} == {
        k: v for k, v in native["metadata"].items() if k not in frame
    }
    assert grid.features[0]["geometry"]["coordinates"][0] == [320000.0, 6399000.0]


def _framed_grid(**kwargs):
    grid = _grid()
    grid.bounds = Bounds(320000.0, 6399000.0, 320100.0, 6399100.0)
    grid.divisions = 2
    for name, value in kwargs.items():
        setattr(grid, name, value)
    return grid


def _positions(payload):
    return [p for f in payload["features"] for p in f["geometry"]["coordinates"]]


def test_calibration_grid_web_payload_names_its_own_frame():
    grid = _framed_grid()
    payload = grid.to_geojson(crs="EPSG:4326")
    metadata = payload["metadata"]

    assert metadata["crs"] == "EPSG:4326"
    xmin, ymin, xmax, ymax = metadata["bounds"]
    for lon, lat in _positions(payload):
        assert xmin <= lon <= xmax and ymin <= lat <= ymax
    assert 11 < xmin < xmax < 12 and 57 < ymin < ymax < 58
    assert metadata["spacing"] == pytest.approx([(xmax - xmin) / 2, (ymax - ymin) / 2])
    assert metadata["divisions"] == 2
    assert metadata["line_count"] == 1


def test_calibration_grid_web_payload_round_trips():
    grid = _framed_grid()
    payload = grid.to_geojson(crs="EPSG:4326")

    reloaded = CalibrationGrid.from_geojson(payload)

    assert reloaded.crs == "EPSG:4326"
    assert reloaded.to_geojson(crs="EPSG:4326") == payload
    assert reloaded.bounds.xmin <= min(p[0] for p in _positions(payload))
    assert reloaded.bounds.xmax >= max(p[0] for p in _positions(payload))
    assert reloaded.bounds.ymin <= min(p[1] for p in _positions(payload))
    assert reloaded.bounds.ymax >= max(p[1] for p in _positions(payload))


def test_calibration_grid_native_round_trip_is_unchanged():
    grid = _framed_grid()

    reloaded = CalibrationGrid.from_geojson(grid.to_geojson())

    assert reloaded.crs == grid.crs
    assert reloaded.bounds == grid.bounds
    assert reloaded.divisions == grid.divisions
    assert reloaded.features == grid.features
    assert reloaded.to_geojson() == grid.to_geojson()


@pytest.mark.parametrize("divisions", [0, 4])
def test_calibration_grid_empty_web_payload_bounds_from_corners(divisions):
    from pyproj import Transformer

    grid = CalibrationGrid(
        bounds=Bounds(320000.0, 6399000.0, 320100.0, 6399100.0),
        divisions=divisions,
        crs="EPSG:3006",
    )
    payload = grid.to_geojson(crs="EPSG:4326")

    transformer = Transformer.from_crs("EPSG:3006", "EPSG:4326", always_xy=True)
    lons, lats = transformer.transform(
        [320000.0, 320100.0, 320100.0, 320000.0],
        [6399000.0, 6399000.0, 6399100.0, 6399100.0],
    )
    metadata = payload["metadata"]
    assert metadata["bounds"] == pytest.approx(
        [min(lons), min(lats), max(lons), max(lats)], abs=1e-9
    )
    xmin, ymin, xmax, ymax = metadata["bounds"]
    expected = [0.0, 0.0]
    if divisions:
        expected = [(xmax - xmin) / divisions, (ymax - ymin) / divisions]
    assert metadata["spacing"] == pytest.approx(expected)
    assert payload == grid.to_geojson(crs="EPSG:4326")
    reloaded = CalibrationGrid.from_geojson(payload)
    assert reloaded.crs == "EPSG:4326"
    assert reloaded.to_geojson(crs="EPSG:4326") == payload


def test_calibration_grid_web_export_leaves_the_grid_unchanged():
    grid = _framed_grid(metadata_payload={"note": "kept", "bounds": [1, 2, 3, 4]})
    before = copy.deepcopy(
        (grid.features, grid.bounds, grid.crs, grid.metadata_payload)
    )

    grid.to_geojson(crs="EPSG:4326")

    assert (grid.features, grid.bounds, grid.crs, grid.metadata_payload) == before


def test_calibration_grid_web_payload_keeps_custom_metadata():
    grid = _framed_grid(metadata_payload={"note": "kept", "dataset": "x"})

    metadata = grid.to_geojson(crs="EPSG:4326")["metadata"]

    assert metadata["note"] == "kept"
    assert metadata["dataset"] == "x"
    assert metadata["crs"] == "EPSG:4326"


def test_calibration_grid_native_names_crs_and_keeps_coordinates():
    payload = _grid().to_geojson()

    assert payload["crs"] == NATIVE_CRS
    assert payload["features"][0]["geometry"]["coordinates"][0] == [
        320000.0,
        6399000.0,
    ]


def test_calibration_grid_without_crs_has_no_crs_member():
    assert "crs" not in _grid(crs=None).to_geojson()


@pytest.mark.parametrize(
    "collection",
    [
        FootprintCollection(),
        BuildingCollection(),
        TreeCollection(),
        CalibrationGrid(),
        City(),
    ],
    ids=lambda c: type(c).__name__,
)
def test_collection_types_declare_web_geojson(collection):
    assert type(collection).web_geojson is True


def _part_building(own_geometry):
    building = Building(id="b1")
    if own_geometry:
        building.add_geometry(_surface(), GeometryType.LOD0)
    part = BuildingPart()
    if not own_geometry:
        part.add_geometry(_surface(), GeometryType.LOD0)
    building.add_child(part)
    return building


def test_footprint_keeps_srs_when_own_geometry_coexists_with_a_part():
    building = _part_building(own_geometry=True)

    assert building.footprint().transform.srs == "EPSG:3006"


def test_footprint_keeps_srs_when_geometry_is_stored_on_a_part():
    building = _part_building(own_geometry=False)

    assert building.footprint().transform.srs == "EPSG:3006"


def test_footprint_srs_comes_from_the_declaring_part():
    building = Building(id="b1")
    blank = BuildingPart()
    blank.add_geometry(_surface(srs=""), GeometryType.LOD0)
    declared = BuildingPart()
    declared.add_geometry(_surface(srs="EPSG:3006"), GeometryType.LOD0)
    building.add_child(blank)
    building.add_child(declared)

    assert building.footprint().transform.srs == "EPSG:3006"


def _clockwise_exterior_building():
    # Clockwise exterior and counter-clockwise hole: the opposite of RFC 7946.
    return _building("b1", surface=_surface(SQUARE[::-1], hole=HOLE))


@pytest.mark.parametrize("crs", [None, "EPSG:4326"], ids=["native", "wgs84"])
@pytest.mark.parametrize(
    "make",
    [
        lambda b: FootprintCollection.from_buildings([b]),
        lambda b: BuildingCollection([b]),
    ],
    ids=["footprints", "buildings"],
)
def test_polygon_rings_are_closed_exterior_ccw_and_hole_cw(make, crs):
    source = _clockwise_exterior_building().footprint()
    assert not LinearRing(source.vertices[:, :2]).is_ccw
    assert LinearRing(source.holes[0][:, :2]).is_ccw

    payload = make(_clockwise_exterior_building()).to_geojson(crs=crs)

    exterior, hole = payload["features"][0]["geometry"]["coordinates"]
    assert exterior[0] == exterior[-1] and hole[0] == hole[-1]
    assert LinearRing(exterior).is_ccw
    assert not LinearRing(hole).is_ccw


def _two_part_building(first_srs, second_srs, part_srs=("", "")):
    building = Building(id="b1")
    for srs, object_srs in zip((first_srs, second_srs), part_srs, strict=True):
        part = BuildingPart()
        part.transform.srs = object_srs
        part.add_geometry(_surface(srs=srs), GeometryType.LOD0)
        building.add_child(part)
    return building


def test_parts_declaring_different_crs_are_refused():
    building = _two_part_building("EPSG:3006", "EPSG:3008")

    with pytest.raises(ValueError, match="EPSG:3006.*EPSG:3008"):
        BuildingCollection([building]).to_geojson(crs="EPSG:4326")
    with pytest.raises(ValueError, match="EPSG:3006.*EPSG:3008"):
        FootprintCollection.from_buildings([building])


def test_parts_declaring_the_same_crs_under_different_spelling_export():
    building = _two_part_building("EPSG:3006", "epsg:3006")

    [feature] = BuildingCollection([building]).to_geojson(crs="EPSG:4326")["features"]

    assert any(
        position == pytest.approx(SQUARE_WGS84_FIRST, abs=1e-9)
        for position in feature["geometry"]["coordinates"][0]
    )


def test_part_object_crs_is_checked_against_the_building():
    building = _two_part_building("EPSG:3006", "EPSG:3006", ("EPSG:4326", ""))
    building.transform.srs = "EPSG:3006"

    with pytest.raises(ValueError, match="Conflicting CRS"):
        BuildingCollection([building]).to_geojson(crs="EPSG:4326")


def test_own_lod0_geometry_decides_the_footprint_srs_over_parts():
    building = _building("b1", surface=_surface(srs="EPSG:3006"))
    part = BuildingPart()
    part.add_geometry(_surface(srs="EPSG:3008"), GeometryType.LOD0)
    building.add_child(part)

    for geom_type in (None, GeometryType.LOD0, "LOD0"):
        assert building.footprint(geom_type).transform.srs == "EPSG:3006"


def _translation():
    affine = np.eye(4)
    affine[0, 3], affine[1, 3] = 320000, 6399000
    return affine


def _shifted_square_surface(srs=""):
    return _surface([(x - 320000, y - 6399000) for x, y in SQUARE], srs=srs)


def _translated_building():
    building = Building(id="b1")
    building.add_geometry(_shifted_square_surface(), GeometryType.LOD0)
    building.transform.affine = _translation()
    building.transform.srs = "EPSG:3006"
    return building


def test_building_owned_translation_is_applied_to_the_footprint():
    building = _translated_building()

    footprint = building.footprint()
    payload = FootprintCollection.from_buildings([building]).to_geojson(crs="EPSG:4326")

    assert sorted(map(tuple, footprint.vertices[:, :2].tolist())) == sorted(
        (float(x), float(y)) for x, y in SQUARE
    )
    assert np.array_equal(footprint.transform.affine, np.eye(4))
    assert footprint.transform.srs == "EPSG:3006"
    assert any(
        position == pytest.approx(SQUARE_WGS84_FIRST, abs=1e-9)
        for position in payload["features"][0]["geometry"]["coordinates"][0]
    )
    assert building.transform.affine[0, 3] == 320000
    assert building.get_geometry(GeometryType.LOD0).vertices[0, 0] == 0


def test_footprint_copies_the_building_only_when_a_transform_needs_applying():
    with patch(
        "dtcc_core.io.geojson.exported_coordinates", wraps=exported_coordinates
    ) as helper:
        _building().footprint()
        assert helper.call_count == 0
        _translated_building().footprint()
        assert helper.call_count == 1


def _geometry_translated_building():
    surface = _shifted_square_surface(srs="EPSG:3006")
    surface.transform.affine = _translation()
    building = Building(id="b1")
    building.add_geometry(surface, GeometryType.LOD0)
    return building


def test_geometry_owned_translation_is_refused():
    building = _geometry_translated_building()

    with pytest.raises(ValueError, match="non-identity transform"):
        building.footprint()
    with pytest.raises(ValueError, match="non-identity transform"):
        FootprintCollection.from_buildings([building])
    with pytest.raises(ValueError, match="non-identity transform"):
        BuildingCollection([building]).to_geojson(crs="EPSG:4326")


def test_part_owned_translation_is_refused():
    building = Building(id="b1")
    part = BuildingPart()
    part.add_geometry(_shifted_square_surface(srs="EPSG:3006"), GeometryType.LOD0)
    part.transform.affine = _translation()
    building.add_child(part)

    with pytest.raises(ValueError, match="non-identity transform"):
        building.footprint()


def test_building_srs_is_a_declaration():
    building = _building("b1", surface=_surface(srs=""))
    building.transform.srs = "EPSG:3006"

    assert building.footprint().transform.srs == "EPSG:3006"

    building = _building("b1", surface=_surface(srs="EPSG:3006"))
    building.transform.srs = "EPSG:3008"
    with pytest.raises(ValueError, match="Conflicting CRS declarations"):
        building.footprint()


def _multi_surface_building(*component_srs, carrier_srs=""):
    carrier = MultiSurface(surfaces=[_surface(srs=srs) for srs in component_srs])
    carrier.transform.srs = carrier_srs
    building = Building(id="b1")
    building.add_geometry(carrier, GeometryType.LOD0)
    return building


def test_component_surface_srs_is_a_declaration():
    building = _multi_surface_building("EPSG:3006")

    assert building.footprint().transform.srs == "EPSG:3006"
    for payload in (
        BuildingCollection([building]).to_geojson(crs="EPSG:4326"),
        FootprintCollection.from_buildings([building]).to_geojson(crs="EPSG:4326"),
    ):
        assert any(
            position == pytest.approx(SQUARE_WGS84_FIRST, abs=1e-9)
            for position in payload["features"][0]["geometry"]["coordinates"][0]
        )


def test_component_surface_srs_conflicting_with_the_building_is_refused():
    building = _multi_surface_building("EPSG:3006")
    building.transform.srs = "EPSG:4326"

    with pytest.raises(ValueError, match="Conflicting CRS declarations"):
        building.footprint()


def test_components_declaring_different_crs_are_refused():
    building = _multi_surface_building("EPSG:3006", "EPSG:3008")

    with pytest.raises(ValueError, match="Conflicting CRS declarations"):
        building.footprint()
