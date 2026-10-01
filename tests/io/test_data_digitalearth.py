import numpy as np
import pytest
from rasterio.io import MemoryFile

from dtcc_core.io.data import digitalearth


# Central Gothenburg, 8 km x 8 km in EPSG:3006, large enough to need tiling.
BBOX = (317000, 6394000, 325000, 6402000)


class DummyResponse:
    """Minimal stand-in for a requests response."""

    def __init__(self, content=b"", status_code=200, content_type="image/png"):
        self.content = content
        self.status_code = status_code
        self.headers = {"Content-Type": content_type}
        self.text = ""


class JsonResponse:
    """Minimal stand-in for a requests response carrying JSON."""

    def __init__(self, payload):
        self.status_code = 200
        self.text = ""
        self._payload = payload

    def json(self):
        return self._payload


class RecordingSession:
    """Session that records the parameters it was called with."""

    def __init__(self, response=None):
        self.response = response or DummyResponse()
        self.calls = []

    def get(self, url, params=None, timeout=None):
        self.calls.append({"url": url, "params": params})
        return self.response


def make_png(bands, width=4, height=4, fill=200):
    """Build an in-memory PNG with the given number of bands."""
    data = np.full((bands, height, width), fill, dtype=np.uint8)
    with MemoryFile() as memfile:
        with memfile.open(
            driver="PNG", width=width, height=height, count=bands, dtype="uint8"
        ) as dst:
            dst.write(data)
        return memfile.read()


def test_to_wgs84_bbox():
    lon_min, lat_min, lon_max, lat_max = digitalearth._to_wgs84_bbox(BBOX)
    # Gothenburg sits near 11.9-12.1 E, 57.6-57.8 N.
    assert 11.5 < lon_min < lon_max < 12.5
    assert 57.5 < lat_min < lat_max < 58.0


def test_tile_grid_covers_bbox_without_gaps():
    width, height, tiles = digitalearth._tile_grid(BBOX, 10.0)
    assert (width, height) == (800, 800)
    assert len(tiles) == 4

    # No tile may exceed what the server is willing to render.
    for tile in tiles:
        assert tile["width"] <= digitalearth.MAX_TILE_PIXELS
        assert tile["height"] <= digitalearth.MAX_TILE_PIXELS

    # Together the tiles must cover every pixel of the mosaic exactly once.
    covered = np.zeros((height, width), dtype=int)
    for tile in tiles:
        rows = slice(tile["row_off"], tile["row_off"] + tile["height"])
        cols = slice(tile["col_off"], tile["col_off"] + tile["width"])
        covered[rows, cols] += 1
    assert (covered == 1).all()

    # And their bounding boxes must span the requested area.
    xmin, ymin, xmax, ymax = BBOX
    assert min(t["bbox"][0] for t in tiles) == pytest.approx(xmin)
    assert min(t["bbox"][1] for t in tiles) == pytest.approx(ymin)
    assert max(t["bbox"][2] for t in tiles) == pytest.approx(xmax)
    assert max(t["bbox"][3] for t in tiles) == pytest.approx(ymax)


def test_tile_grid_small_bbox_is_a_single_tile():
    width, height, tiles = digitalearth._tile_grid((0, 0, 1000, 1000), 10.0)
    assert (width, height) == (100, 100)
    assert len(tiles) == 1


def test_layer_for_date_switches_at_2022():
    assert digitalearth._layer_for_date("2021-12-31") == digitalearth.LAYER_BEFORE_2022
    assert digitalearth._layer_for_date("2022-01-01") == digitalearth.LAYER_FROM_2022
    assert digitalearth._layer_for_date("2026-05-05") == digitalearth.LAYER_FROM_2022


def test_acquisition_datetime_falls_back_to_start_datetime():
    # Every item in the live catalogue leaves 'datetime' null, so the
    # fallback is the normal path rather than an edge case.
    feature = {"properties": {"datetime": None, "start_datetime": "2026-05-05T10:27:01Z"}}
    assert digitalearth._acquisition_datetime(feature) == "2026-05-05T10:27:01Z"

    feature = {"properties": {"datetime": "2026-06-01T10:20:41Z"}}
    assert digitalearth._acquisition_datetime(feature) == "2026-06-01T10:20:41Z"

    assert digitalearth._acquisition_datetime({"properties": {}}) is None


def test_get_map_uses_northing_first_bbox():
    # EPSG:3006 declares northing before easting, and WMS 1.3.0 follows the
    # CRS. Sending easting first returns a blank image instead of an error,
    # so the ordering is pinned down here.
    session = RecordingSession()
    digitalearth.get_map(
        session, "sentinel2_reflectance", "rgb", (317000, 6394000, 321000, 6398000),
        512, 512, "2026-05-05",
    )
    params = session.calls[0]["params"]
    assert params["bbox"] == "6394000,317000,6398000,321000"
    assert params["crs"] == "EPSG:3006"
    assert params["time"] == "2026-05-05"
    assert params["width"] == 512


def test_get_map_raises_on_service_exception():
    # The server reports errors as XML with status 200.
    session = RecordingSession(
        DummyResponse(content=b"<ServiceExceptionReport/>", content_type="application/xml")
    )
    with pytest.raises(RuntimeError, match="WMS returned an error"):
        digitalearth.get_map(
            session, "sentinel2_reflectance", "rgb", BBOX, 32, 32, "2026-05-05"
        )


def test_get_map_raises_on_bad_status():
    session = RecordingSession(DummyResponse(status_code=500))
    with pytest.raises(RuntimeError, match="status 500"):
        digitalearth.get_map(
            session, "sentinel2_reflectance", "rgb", BBOX, 32, 32, "2026-05-05"
        )


def test_tile_to_array_treats_greyscale_as_transparent():
    # A date with no data comes back as a tiny greyscale image.
    data = digitalearth._tile_to_array(make_png(1), 4, 4)
    assert data.shape == (4, 4, 4)
    assert not data[3].any()


def test_tile_to_array_adds_alpha_to_rgb():
    data = digitalearth._tile_to_array(make_png(3), 4, 4)
    assert data.shape == (4, 4, 4)
    assert (data[3] == 255).all()


def test_tile_to_array_keeps_rgba():
    data = digitalearth._tile_to_array(make_png(4), 4, 4)
    assert data.shape == (4, 4, 4)
    assert (data[0] == 200).all()


def test_list_available_dates_filters_and_sorts(monkeypatch):
    features = [
        # One date spread over two satellite tiles, one clear and one not.
        {"properties": {"start_datetime": "2026-06-21T10:00:00Z", "eo:cloud_cover": 0.1}},
        {"properties": {"start_datetime": "2026-06-21T10:00:00Z", "eo:cloud_cover": 80.0}},
        {"properties": {"start_datetime": "2026-07-15T10:00:00Z", "eo:cloud_cover": 3.0}},
        {"properties": {"start_datetime": "2026-08-24T10:00:00Z", "eo:cloud_cover": 7.0}},
        {"properties": {"start_datetime": "2026-09-01T10:00:00Z", "eo:cloud_cover": 99.0}},
        # Items without usable metadata are skipped rather than crashing.
        {"properties": {"start_datetime": None, "eo:cloud_cover": 1.0}},
        {"properties": {"start_datetime": "2026-09-05T10:00:00Z"}},
    ]
    monkeypatch.setattr(digitalearth, "_stac_search_all", lambda session, payload: features)

    dates = digitalearth.list_available_dates(BBOX, max_cloud=10.0, session=object())

    # Newest first, and the part-cloudy date is judged by its worst tile.
    assert dates == [("2026-08-24", 7.0), ("2026-07-15", 3.0)]


def test_list_available_dates_respects_max_cloud(monkeypatch):
    features = [
        {"properties": {"start_datetime": "2026-07-15T10:00:00Z", "eo:cloud_cover": 25.0}}
    ]
    monkeypatch.setattr(digitalearth, "_stac_search_all", lambda session, payload: features)

    assert digitalearth.list_available_dates(BBOX, max_cloud=10.0, session=object()) == []
    assert digitalearth.list_available_dates(BBOX, max_cloud=30.0, session=object()) == [
        ("2026-07-15", 25.0)
    ]


def test_latest_clear_date_skips_dates_without_pixels(monkeypatch):
    # The catalogue lists acquisitions the imagery server has no pixels for,
    # so the newest clear date is not always the one that can be downloaded.
    monkeypatch.setattr(
        digitalearth,
        "list_available_dates",
        lambda *a, **k: [("2026-08-24", 7.0), ("2026-07-15", 3.0)],
    )
    probed = []

    def fake_has_imagery(session, layer, style, bbox, date):
        probed.append(date)
        return date == "2026-07-15"

    monkeypatch.setattr(digitalearth, "_has_imagery", fake_has_imagery)

    assert digitalearth.latest_clear_date(BBOX, session=object()) == "2026-07-15"
    assert probed == ["2026-08-24", "2026-07-15"]


def test_latest_clear_date_raises_when_catalogue_is_empty(monkeypatch):
    monkeypatch.setattr(digitalearth, "list_available_dates", lambda *a, **k: [])
    with pytest.raises(RuntimeError, match="cloud cover"):
        digitalearth.latest_clear_date(BBOX, session=object())


def test_latest_clear_date_raises_when_no_candidate_has_pixels(monkeypatch):
    monkeypatch.setattr(
        digitalearth, "list_available_dates", lambda *a, **k: [("2026-08-24", 7.0)]
    )
    monkeypatch.setattr(digitalearth, "_has_imagery", lambda *a, **k: False)
    with pytest.raises(RuntimeError, match="none of the"):
        digitalearth.latest_clear_date(BBOX, session=object())


def test_stac_search_all_follows_next_links():
    # The catalogue returns results oldest first and pages them, so the newest
    # acquisitions are only seen once every page has been read.
    pages = [
        {
            "features": [{"id": "a"}],
            "links": [{"rel": "next", "href": "https://example/stac/search", "body": {"_o": 1}}],
        },
        {"features": [{"id": "b"}], "links": []},
    ]

    class PagingSession:
        def __init__(self):
            self.bodies = []

        def post(self, url, json=None, timeout=None):
            self.bodies.append(json)
            return JsonResponse(pages[len(self.bodies) - 1])

    session = PagingSession()
    features = digitalearth._stac_search_all(session, {"collections": ["s2_msi_l2a"]})

    assert [f["id"] for f in features] == ["a", "b"]
    assert session.bodies[1]["_o"] == 1


def test_download_imagery_rejects_unknown_style():
    with pytest.raises(ValueError, match="Invalid style"):
        digitalearth.download_imagery(BBOX, date="2026-05-05", style="nonsense")


def test_download_imagery_rejects_inverted_bounds():
    with pytest.raises(ValueError, match="Invalid bounds"):
        digitalearth.download_imagery((325000, 6402000, 317000, 6394000), date="2026-05-05")


def test_stac_search_all_warns_when_page_limit_is_hit(monkeypatch):
    # Results arrive oldest first, so stopping early loses the newest
    # acquisitions silently unless the truncation is reported.
    page = {
        "features": [{"id": "x"}],
        "links": [{"rel": "next", "href": "https://example/stac/search", "body": {"_o": 1}}],
    }

    class EndlessSession:
        def post(self, url, json=None, timeout=None):
            return JsonResponse(page)

    warnings = []
    monkeypatch.setattr(digitalearth, "warning", warnings.append)

    features = digitalearth._stac_search_all(EndlessSession(), {}, max_pages=3)

    assert len(features) == 3
    assert len(warnings) == 1
    assert "3 pages" in warnings[0]


def test_stac_search_all_does_not_warn_when_results_run_out(monkeypatch):
    pages = [
        {
            "features": [{"id": "a"}],
            "links": [{"rel": "next", "href": "https://example/stac/search", "body": {"_o": 1}}],
        },
        {"features": [{"id": "b"}], "links": []},
    ]

    class PagingSession:
        def __init__(self):
            self.count = 0

        def post(self, url, json=None, timeout=None):
            self.count += 1
            return JsonResponse(pages[self.count - 1])

    warnings = []
    monkeypatch.setattr(digitalearth, "warning", warnings.append)

    digitalearth._stac_search_all(PagingSession(), {}, max_pages=10)
    assert warnings == []


def _stub_download(monkeypatch, tmp_path, alpha_fraction):
    """Wire download_imagery up to tiles with a given fraction of coverage."""
    monkeypatch.setattr(digitalearth, "CACHE_DIR", str(tmp_path))
    monkeypatch.setattr(digitalearth, "create_retry_session", lambda: object())
    monkeypatch.setattr(digitalearth, "get_map", lambda *a, **k: b"")

    def fake_tile(png_bytes, width, height):
        data = np.zeros((4, height, width), dtype=np.uint8)
        opaque_rows = int(round(height * alpha_fraction))
        data[:, :opaque_rows, :] = 255
        return data

    monkeypatch.setattr(digitalearth, "_tile_to_array", fake_tile)

    warnings = []
    monkeypatch.setattr(digitalearth, "warning", warnings.append)
    return warnings


def test_download_imagery_warns_on_partial_coverage(monkeypatch, tmp_path):
    warnings = _stub_download(monkeypatch, tmp_path, alpha_fraction=0.5)
    digitalearth.download_imagery((0, 0, 1000, 1000), date="2026-05-05")
    assert len(warnings) == 1
    assert "50%" in warnings[0]


def test_download_imagery_warns_when_empty(monkeypatch, tmp_path):
    warnings = _stub_download(monkeypatch, tmp_path, alpha_fraction=0.0)
    digitalearth.download_imagery((0, 0, 1000, 1000), date="2026-05-05")
    assert len(warnings) == 1
    assert "is empty" in warnings[0]


def test_download_imagery_quiet_on_full_coverage(monkeypatch, tmp_path):
    warnings = _stub_download(monkeypatch, tmp_path, alpha_fraction=1.0)
    digitalearth.download_imagery((0, 0, 1000, 1000), date="2026-05-05")
    assert warnings == []
