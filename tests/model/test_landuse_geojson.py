import geopandas as gpd
import pytest
from shapely.geometry import box

from dtcc_core.io.data import urban_atlas
from dtcc_core.model import Landuse

X0, Y0 = 319000.0, 6399000.0


def _urban_atlas(count):
    """Urban Atlas Landuse with ``count`` polygons and requested bounds."""
    codes = ["14100", "11100", "12220", "50000"]
    gdf = gpd.GeoDataFrame(
        [
            {
                "ua_code": codes[i],
                "ua_class": f"class {codes[i]}",
                "fua_name": "Göteborg",
                "geometry": box(X0 + 100 * i, Y0, X0 + 100 * (i + 1), Y0 + 100),
            }
            for i in range(count)
        ],
        crs="EPSG:3006",
    )
    landuse = urban_atlas.landuse_from_urban_atlas(gdf)
    landuse.attributes["requested_bounds"] = [X0, Y0, X0 + 400.0, Y0 + 100.0]
    return landuse


def test_each_polygon_is_a_feature_with_its_class_and_codes():
    payload = _urban_atlas(2).to_geojson()

    assert payload["crs"]["properties"]["name"] == "EPSG:3006"
    assert [f["properties"] for f in payload["features"]] == [
        {
            "landuse": "GRASS",
            "ua_code": "14100",
            "ua_class": "class 14100",
            "fua_name": "Göteborg",
        },
        {
            "landuse": "HEAVY_URBAN",
            "ua_code": "11100",
            "ua_class": "class 11100",
            "fua_name": "Göteborg",
        },
    ]
    ring = payload["features"][0]["geometry"]["coordinates"][0]
    assert ring[0] == ring[-1]
    assert {tuple(p) for p in ring} == {
        (X0, Y0),
        (X0 + 100, Y0),
        (X0 + 100, Y0 + 100),
        (X0, Y0 + 100),
    }


def test_collection_attributes_stay_metadata_when_they_match_the_polygon_count():
    # Four polygons and four requested-bounds values: the bounds stay metadata.
    payload = _urban_atlas(4).to_geojson()

    assert len(payload["features"]) == 4
    assert "requested_bounds" not in payload["features"][0]["properties"]
    assert payload["metadata"]["requested_bounds"] == [
        X0,
        Y0,
        X0 + 400.0,
        Y0 + 100.0,
    ]
    assert payload["metadata"]["year"] == 2012
    assert "ua_code" not in payload["metadata"]


def test_wgs84_output_is_rfc7946():
    payload = _urban_atlas(1).to_geojson(crs="EPSG:4326")

    assert "crs" not in payload
    lon, lat = payload["features"][0]["geometry"]["coordinates"][0][0]
    assert lon == pytest.approx(11.97, abs=0.1)
    assert lat == pytest.approx(57.7, abs=0.1)


def test_classes_are_null_when_they_do_not_match_the_surfaces():
    landuse = _urban_atlas(2)
    landuse.landuses = landuse.landuses[:1]

    payload = landuse.to_geojson()

    assert [f["properties"]["landuse"] for f in payload["features"]] == [None, None]


def test_empty_landuse_with_a_crs_has_no_features():
    payload = Landuse(attributes={"crs": "EPSG:3006"}).to_geojson()

    assert payload["features"] == []
    assert payload["metadata"] == {"crs": "EPSG:3006"}
