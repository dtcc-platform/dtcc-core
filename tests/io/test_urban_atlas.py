import pytest
import requests
from shapely.geometry import box, mapping

from dtcc_core.datasets.dataset import DatasetUpstreamError
from dtcc_core.io.data import urban_atlas
from dtcc_core.model import Bounds, Landuse, LanduseClasses

BOUNDS = Bounds(0.0, 0.0, 100.0, 100.0)


def _feature(code, label, geometry):
    return {
        "type": "Feature",
        "properties": {"code_2018": code, "class_2018": label, "fua_name": "Göteborg"},
        "geometry": mapping(geometry),
    }


class FakeResponse:
    def __init__(self, payload):
        self._payload = payload

    def raise_for_status(self):
        return None

    def json(self):
        return self._payload


@pytest.fixture
def fake_service(monkeypatch, tmp_path):
    """Serve two pages of non-road features and one road feature."""
    calls = []
    pages = [
        {
            "features": [
                _feature("11210", "Discontinuous dense urban fabric", box(0, 0, 50, 50))
            ],
            "properties": {"exceededTransferLimit": True},
        },
        {"features": [_feature("14100", "Green urban areas", box(50, 0, 150, 50))]},
    ]
    road = {"features": [_feature("12220", "Other roads", box(-500, 50, 500, 60))]}

    def fake_get(url, params, timeout):
        calls.append(dict(params))
        if params["where"].endswith("= '12220'"):
            return FakeResponse(road)
        return FakeResponse(pages[params["resultOffset"]])

    monkeypatch.setattr(urban_atlas.requests, "get", fake_get)
    monkeypatch.setattr(urban_atlas, "cache_dir", tmp_path)
    return calls


def test_download_paginates_clips_and_fetches_road_separately(fake_service):
    gdf = urban_atlas.download_urban_atlas_geodataframe(BOUNDS)

    assert gdf["ua_code"].tolist() == ["11210", "14100", "12220"]
    assert gdf.crs.to_epsg() == 3006
    # Everything is clipped to the requested bounds.
    assert gdf.total_bounds.tolist() == [0.0, 0.0, 100.0, 60.0]
    assert gdf.area.tolist() == pytest.approx([2500.0, 2500.0, 1000.0])

    non_road = [c for c in fake_service if c["where"].endswith("<> '12220'")]
    road = [c for c in fake_service if c["where"].endswith("= '12220'")]
    assert [c["resultOffset"] for c in non_road] == [0, 1]
    assert "maxAllowableOffset" not in non_road[0]
    assert road[0]["maxAllowableOffset"] == urban_atlas.URBAN_ATLAS_ROAD_SIMPLIFY_M
    assert non_road[0]["inSR"] == 3006 and non_road[0]["outSR"] == 3006


def test_download_uses_cache_on_second_call(fake_service):
    urban_atlas.download_urban_atlas_geodataframe(BOUNDS)
    call_count = len(fake_service)

    cached = urban_atlas.download_urban_atlas_geodataframe(BOUNDS)

    assert len(fake_service) == call_count
    assert cached["ua_code"].tolist() == ["11210", "14100", "12220"]


def test_download_urban_atlas_returns_landuse_with_codes(fake_service):
    landuse = urban_atlas.download_urban_atlas(BOUNDS)

    assert isinstance(landuse, Landuse)
    assert len(landuse.surfaces) == 3
    assert landuse.landuses == [
        LanduseClasses.URBAN,
        LanduseClasses.GRASS,
        LanduseClasses.ROAD,
    ]
    assert landuse.attributes["ua_code"] == ["11210", "14100", "12220"]
    assert landuse.attributes["ua_class"][1] == "Green urban areas"
    assert landuse.attributes["fua_name"] == ["Göteborg"] * 3
    assert landuse.attributes["year"] == 2018
    assert "Copernicus" in landuse.attributes["attribution"]


def test_landuse_keeps_holes():
    import geopandas as gpd

    block = box(0, 0, 10, 10).difference(box(4, 4, 6, 6))
    gdf = gpd.GeoDataFrame(
        [
            {
                "ua_code": "11100",
                "ua_class": "Continuous urban fabric",
                "fua_name": "X",
                "geometry": block,
            }
        ],
        crs="EPSG:3006",
    )

    landuse = urban_atlas.landuse_from_urban_atlas(gdf)

    assert landuse.landuses == [LanduseClasses.HEAVY_URBAN]
    assert len(landuse.surfaces[0].holes) == 1


def test_every_urban_atlas_class_is_mapped():
    assert len(urban_atlas.URBAN_ATLAS_LANDUSE_MAP) == 27
    assert LanduseClasses.UNKNOWN not in urban_atlas.URBAN_ATLAS_LANDUSE_MAP.values()


def test_empty_area_warns(monkeypatch, tmp_path):
    warnings = []
    monkeypatch.setattr(
        urban_atlas.requests, "get", lambda *a, **k: FakeResponse({"features": []})
    )
    monkeypatch.setattr(urban_atlas, "cache_dir", tmp_path)
    monkeypatch.setattr(urban_atlas, "warning", warnings.append)

    gdf = urban_atlas.download_urban_atlas_geodataframe(BOUNDS)

    assert len(gdf) == 0
    assert "Functional Urban Areas" in warnings[0]


class ErrorResponse(FakeResponse):
    def __init__(self, status_code):
        super().__init__({})
        self.status_code = status_code

    def raise_for_status(self):
        response = requests.Response()
        response.status_code = self.status_code
        raise requests.HTTPError(f"HTTP {self.status_code}", response=response)


class BadJsonResponse(FakeResponse):
    def json(self):
        raise ValueError("Expecting value")


def _raise(exc):
    def fail(*args, **kwargs):
        raise exc

    return fail


@pytest.mark.parametrize(
    "get, failure_class, status_code",
    [
        (
            lambda *a, **k: FakeResponse({"error": {"code": 500, "message": "boom"}}),
            "http_5xx",
            500,
        ),
        (
            lambda *a, **k: FakeResponse(
                {"error": {"code": 400, "message": "bad where"}}
            ),
            "invalid_payload",
            400,
        ),
        (_raise(requests.exceptions.Timeout("read timed out")), "timeout", None),
        (_raise(requests.exceptions.ConnectionError("no route")), "connection", None),
        (lambda *a, **k: ErrorResponse(503), "http_5xx", 503),
        (lambda *a, **k: ErrorResponse(404), "http_4xx", 404),
        (lambda *a, **k: BadJsonResponse({}), "invalid_json", None),
        (
            lambda *a, **k: FakeResponse({"type": "FeatureCollection"}),
            "invalid_payload",
            None,
        ),
        (lambda *a, **k: FakeResponse([]), "invalid_payload", None),
        (
            lambda *a, **k: FakeResponse(
                {"features": [{"properties": {}, "geometry": {}}]}
            ),
            "invalid_payload",
            None,
        ),
        (
            lambda *a, **k: FakeResponse(
                {
                    "features": [
                        {**_feature("11100", "X", box(0, 0, 1, 1)), "geometry": None}
                    ]
                }
            ),
            "invalid_payload",
            None,
        ),
    ],
)
def test_upstream_failures_raise_typed_errors_and_are_not_cached(
    monkeypatch, tmp_path, get, failure_class, status_code
):
    monkeypatch.setattr(urban_atlas.requests, "get", get)
    monkeypatch.setattr(urban_atlas, "cache_dir", tmp_path)

    with pytest.raises(DatasetUpstreamError) as error:
        urban_atlas.download_urban_atlas_geodataframe(BOUNDS)

    assert error.value.dataset == "urban_atlas"
    assert error.value.failure_class == failure_class
    assert error.value.status_code == status_code
    assert list(tmp_path.rglob("*.gpkg")) == []


def test_paging_stops_after_the_page_limit(monkeypatch, tmp_path):
    monkeypatch.setattr(urban_atlas, "URBAN_ATLAS_MAX_PAGES", 3)
    monkeypatch.setattr(urban_atlas, "cache_dir", tmp_path)
    page = {
        "features": [_feature("11100", "X", box(0, 0, 1, 1))],
        "exceededTransferLimit": True,
    }
    calls = []
    monkeypatch.setattr(
        urban_atlas.requests,
        "get",
        lambda *a, **k: calls.append(1) or FakeResponse(page),
    )

    with pytest.raises(DatasetUpstreamError, match="after 3 pages") as error:
        urban_atlas.download_urban_atlas_geodataframe(BOUNDS)

    assert error.value.failure_class == "invalid_payload"
    assert len(calls) == 3


def test_nearby_bounds_do_not_share_a_cache_entry(fake_service):
    # These round to the same whole metres; the cached clip must not be reused.
    first = urban_atlas.download_urban_atlas_geodataframe(Bounds(0.0, 0.0, 99.6, 100.0))
    calls = len(fake_service)

    second = urban_atlas.download_urban_atlas_geodataframe(
        Bounds(0.0, 0.0, 100.4, 100.0)
    )

    assert len(fake_service) > calls
    assert first.total_bounds[2] == pytest.approx(99.6)
    assert second.total_bounds[2] == pytest.approx(100.4)


def test_empty_result_round_trips_through_the_cache(monkeypatch, tmp_path):
    monkeypatch.setattr(
        urban_atlas.requests, "get", lambda *a, **k: FakeResponse({"features": []})
    )
    monkeypatch.setattr(urban_atlas, "cache_dir", tmp_path)
    urban_atlas.download_urban_atlas_geodataframe(BOUNDS)
    monkeypatch.setattr(
        urban_atlas.requests, "get", _raise(AssertionError("cache not used"))
    )

    cached = urban_atlas.download_urban_atlas_geodataframe(BOUNDS)
    landuse = urban_atlas.landuse_from_urban_atlas(cached)

    assert len(cached) == 0
    assert list(cached.columns) == urban_atlas.URBAN_ATLAS_COLUMNS
    assert landuse.surfaces == [] and landuse.attributes["ua_code"] == []


def test_landuse_records_requested_bounds(fake_service):
    landuse = urban_atlas.download_urban_atlas(BOUNDS)

    assert landuse.attributes["requested_bounds"] == [0.0, 0.0, 100.0, 100.0]


def test_bounds_accept_four_or_six_values():
    assert urban_atlas._ensure_bounds([0, 1, 2, 3]).tuple == (0, 1, 2, 3)
    assert urban_atlas._ensure_bounds([0, 1, -5, 2, 3, 50]).tuple == (0, 1, 2, 3)
    with pytest.raises(ValueError, match="got 5 values"):
        urban_atlas._ensure_bounds([0, 1, 2, 3, 4])


def test_unsupported_year_raises():
    with pytest.raises(ValueError, match="Unsupported Urban Atlas year"):
        urban_atlas.download_urban_atlas_geodataframe(BOUNDS, year=2021)
