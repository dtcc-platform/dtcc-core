import datetime

import numpy as np
import pytest
import requests
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
        self.request = None


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


def test_acquisition_datetime_falls_back_to_start_datetime():
    # Every item in the live catalogue leaves 'datetime' null, so the
    # fallback is the normal path rather than an edge case.
    feature = {
        "properties": {"datetime": None, "start_datetime": "2026-05-05T10:27:01Z"}
    }
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
        session,
        "sentinel2_reflectance",
        "rgb",
        (317000, 6394000, 321000, 6398000),
        512,
        512,
        "2026-05-05",
    )
    params = session.calls[0]["params"]
    assert params["bbox"] == "6394000,317000,6398000,321000"
    assert params["crs"] == "EPSG:3006"
    assert params["time"] == "2026-05-05"
    assert params["width"] == 512


def test_get_map_raises_on_service_exception():
    # The server reports errors as XML with status 200.
    session = RecordingSession(
        DummyResponse(
            content=b"<ServiceExceptionReport/>", content_type="application/xml"
        )
    )
    with pytest.raises(digitalearth.ServiceError, match="WMS returned an error"):
        digitalearth.get_map(
            session, "sentinel2_reflectance", "rgb", BBOX, 32, 32, "2026-05-05"
        )


def test_get_map_rejects_date_outside_archive():
    response = DummyResponse(status_code=400, content_type="application/xml")
    response.text = '<ServiceException code="InvalidDimensionValue"/>'
    with pytest.raises(ValueError, match="outside the date range"):
        digitalearth.get_map(
            RecordingSession(response),
            "sentinel2_reflectance",
            "rgb",
            BBOX,
            32,
            32,
            "2099-01-01",
        )


def test_get_map_raises_on_bad_status():
    session = RecordingSession(DummyResponse(status_code=500))
    with pytest.raises(requests.HTTPError, match="status 500"):
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
        {
            "properties": {
                "start_datetime": "2026-06-21T10:00:00Z",
                "eo:cloud_cover": 0.1,
            }
        },
        {
            "properties": {
                "start_datetime": "2026-06-21T10:00:00Z",
                "eo:cloud_cover": 80.0,
            }
        },
        {
            "properties": {
                "start_datetime": "2026-07-15T10:00:00Z",
                "eo:cloud_cover": 3.0,
            }
        },
        {
            "properties": {
                "start_datetime": "2026-08-24T10:00:00Z",
                "eo:cloud_cover": 7.0,
            }
        },
        {
            "properties": {
                "start_datetime": "2026-09-01T10:00:00Z",
                "eo:cloud_cover": 99.0,
            }
        },
        # Items without usable metadata are skipped rather than crashing.
        {"properties": {"start_datetime": None, "eo:cloud_cover": 1.0}},
        {"properties": {"start_datetime": "2026-09-05T10:00:00Z"}},
    ]
    monkeypatch.setattr(
        digitalearth, "_stac_search_all", lambda session, payload: features
    )

    dates = digitalearth.list_available_dates(BBOX, max_cloud=10.0, session=object())

    # Newest first, and the part-cloudy date is judged by its worst tile.
    assert dates == [("2026-08-24", 7.0), ("2026-07-15", 3.0)]


def test_list_available_dates_respects_max_cloud(monkeypatch):
    features = [
        {
            "properties": {
                "start_datetime": "2026-07-15T10:00:00Z",
                "eo:cloud_cover": 25.0,
            }
        }
    ]
    monkeypatch.setattr(
        digitalearth, "_stac_search_all", lambda session, payload: features
    )

    assert (
        digitalearth.list_available_dates(BBOX, max_cloud=10.0, session=object()) == []
    )
    assert digitalearth.list_available_dates(
        BBOX, max_cloud=30.0, session=object()
    ) == [("2026-07-15", 25.0)]


def test_latest_clear_date_skips_dates_without_pixels(monkeypatch):
    # The catalogue lists acquisitions the imagery server has no pixels for,
    # so the newest clear date is not always the one that can be downloaded.
    monkeypatch.setattr(
        digitalearth,
        "list_available_dates",
        lambda *a, **k: [("2026-08-24", 7.0), ("2026-07-15", 3.0)],
    )
    probed = []

    def fake_coverage(session, bbox, date):
        probed.append(date)
        return 1.0 if date == "2026-07-15" else 0.0

    monkeypatch.setattr(digitalearth, "_imagery_coverage", fake_coverage)

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
    monkeypatch.setattr(digitalearth, "_imagery_coverage", lambda *a, **k: 0.0)
    with pytest.raises(RuntimeError, match="none of the"):
        digitalearth.latest_clear_date(BBOX, session=object())


def test_imagery_coverage_always_probes_rgb():
    session = RecordingSession(DummyResponse(content=make_png(4, 32, 32)))
    assert digitalearth._imagery_coverage(session, BBOX, "2026-05-05") == 1.0
    assert session.calls[0]["params"]["styles"] == "rgb"
    assert session.calls[0]["params"]["width"] == digitalearth.PROBE_PIXELS


def test_stac_search_all_follows_next_links():
    # The catalogue returns results oldest first and pages them, so the newest
    # acquisitions are only seen once every page has been read.
    pages = [
        {
            "features": [{"id": "a"}],
            "links": [
                {
                    "rel": "next",
                    "href": "https://example/stac/search",
                    "body": {"_o": 1},
                }
            ],
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
        digitalearth.download_imagery(
            (325000, 6402000, 317000, 6394000), date="2026-05-05"
        )


def test_stac_search_all_warns_when_page_limit_is_hit(monkeypatch):
    # Results arrive oldest first, so stopping early loses the newest
    # acquisitions silently unless the truncation is reported.
    page = {
        "features": [{"id": "x"}],
        "links": [
            {"rel": "next", "href": "https://example/stac/search", "body": {"_o": 1}}
        ],
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
            "links": [
                {
                    "rel": "next",
                    "href": "https://example/stac/search",
                    "body": {"_o": 1},
                }
            ],
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


def _stub_download(monkeypatch, cache_dir, alpha_fraction, probe_coverage=1.0):
    """Wire download_imagery up to tiles with a given fraction of coverage."""
    monkeypatch.setattr(digitalearth, "CACHE_DIR", str(cache_dir))
    monkeypatch.setattr(digitalearth, "create_session", lambda: object())
    monkeypatch.setattr(
        digitalearth, "_imagery_coverage", lambda *a, **k: probe_coverage
    )

    calls = []

    def fake_get_map(session, layer, style, bbox, width, height, date):
        calls.append({"layer": layer, "style": style, "bbox": bbox})
        return b""

    monkeypatch.setattr(digitalearth, "get_map", fake_get_map)

    def fake_tile(png_bytes, width, height):
        data = np.zeros((4, height, width), dtype=np.uint8)
        opaque_rows = int(round(height * alpha_fraction))
        data[:, :opaque_rows, :] = 255
        return data

    monkeypatch.setattr(digitalearth, "_tile_to_array", fake_tile)

    warnings = []
    monkeypatch.setattr(digitalearth, "warning", warnings.append)
    return warnings, calls


def _cached_files(cache_dir):
    return sorted(p.name for p in cache_dir.iterdir()) if cache_dir.exists() else []


def test_download_imagery_uses_one_layer_for_all_dates(monkeypatch, tmp_path):
    _, calls = _stub_download(monkeypatch, tmp_path, alpha_fraction=1.0)
    digitalearth.download_imagery((0, 0, 1000, 1000), date="2019-06-13")
    assert {call["layer"] for call in calls} == {"sentinel2_reflectance"}


def test_download_imagery_returns_rgba_raster(monkeypatch, tmp_path):
    _stub_download(monkeypatch, tmp_path, alpha_fraction=1.0)
    raster = digitalearth.download_imagery((0, 0, 1000, 2000), date="2026-05-05")
    assert raster.data.shape == (200, 100, 4)
    assert raster.crs == "EPSG:3006"
    assert raster.georef == digitalearth.from_origin(0, 2000, 10.0, 10.0)


def test_download_imagery_caches_complete_result(monkeypatch, tmp_path):
    warnings, _ = _stub_download(monkeypatch, tmp_path, alpha_fraction=1.0)
    first = digitalearth.download_imagery((0, 0, 1000, 1000), date="2026-05-05")
    assert warnings == []
    assert len(_cached_files(tmp_path)) == 1

    def no_network(*args, **kwargs):
        raise AssertionError("cache hit must not contact the server")

    monkeypatch.setattr(digitalearth, "get_map", no_network)
    monkeypatch.setattr(digitalearth, "_imagery_coverage", no_network)
    second = digitalearth.download_imagery((0, 0, 1000, 1000), date="2026-05-05")
    np.testing.assert_array_equal(second.data, first.data)
    assert second.georef == first.georef


def test_download_imagery_raises_when_explicit_date_has_no_imagery(
    monkeypatch, tmp_path
):
    _, calls = _stub_download(
        monkeypatch, tmp_path, alpha_fraction=1.0, probe_coverage=0.0
    )
    with pytest.raises(RuntimeError, match="no imagery for 2026-10-16"):
        digitalearth.download_imagery((0, 0, 1000, 1000), date="2026-10-16")
    assert calls == []
    assert _cached_files(tmp_path) == []


def test_download_imagery_raises_on_empty_mosaic_and_caches_nothing(
    monkeypatch, tmp_path
):
    _stub_download(monkeypatch, tmp_path, alpha_fraction=0.0)
    with pytest.raises(RuntimeError, match="is empty"):
        digitalearth.download_imagery((0, 0, 1000, 1000), date="2026-05-05")
    assert _cached_files(tmp_path) == []


def test_download_imagery_returns_partial_result_without_caching(monkeypatch, tmp_path):
    warnings, _ = _stub_download(monkeypatch, tmp_path, alpha_fraction=0.5)
    raster = digitalearth.download_imagery((0, 0, 1000, 1000), date="2026-05-05")
    assert raster.data.shape == (100, 100, 4)
    assert len(warnings) == 1
    assert "50%" in warnings[0] and "not cached" in warnings[0]
    assert _cached_files(tmp_path) == []


def test_download_imagery_leaves_no_file_after_interrupted_write(monkeypatch, tmp_path):
    _stub_download(monkeypatch, tmp_path, alpha_fraction=1.0)
    real_open = digitalearth.rasterio.open

    class InterruptingDataset:
        def __init__(self, dataset):
            self.dataset = dataset

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            self.dataset.close()
            return False

        def write(self, data):
            raise KeyboardInterrupt

    def interrupting_open(path, mode="r", **kwargs):
        dataset = real_open(path, mode, **kwargs)
        return InterruptingDataset(dataset) if mode == "w" else dataset

    monkeypatch.setattr(digitalearth.rasterio, "open", interrupting_open)
    with pytest.raises(KeyboardInterrupt):
        digitalearth.download_imagery((0, 0, 1000, 1000), date="2026-05-05")
    assert _cached_files(tmp_path) == []


def test_download_imagery_recreates_missing_cache_dir(monkeypatch, tmp_path):
    cache_dir = tmp_path / "removed" / "downloaded_des"
    _stub_download(monkeypatch, cache_dir, alpha_fraction=1.0)
    digitalearth.download_imagery((0, 0, 1000, 1000), date="2026-05-05")
    assert len(_cached_files(cache_dir)) == 1


def test_download_imagery_caches_masked_index_style(monkeypatch, tmp_path):
    # ndwi hides land, so a dry area is fully transparent but not missing.
    warnings, _ = _stub_download(
        monkeypatch, tmp_path, alpha_fraction=0.0, probe_coverage=1.0
    )
    raster = digitalearth.download_imagery(
        (0, 0, 1000, 1000), date="2026-05-05", style="ndwi"
    )
    assert not raster.data[..., 3].any()
    assert warnings == []
    assert len(_cached_files(tmp_path)) == 1


def test_download_imagery_index_style_uses_rgb_probe_for_coverage(
    monkeypatch, tmp_path
):
    warnings, _ = _stub_download(
        monkeypatch, tmp_path, alpha_fraction=1.0, probe_coverage=0.5
    )
    digitalearth.download_imagery((0, 0, 1000, 1000), date="2026-05-05", style="ndvi")
    assert len(warnings) == 1 and "50%" in warnings[0]
    assert _cached_files(tmp_path) == []


def test_download_imagery_reuses_date_probe_for_index_style(monkeypatch, tmp_path):
    _stub_download(monkeypatch, tmp_path, alpha_fraction=0.0)
    monkeypatch.setattr(
        digitalearth, "_pick_clear_date", lambda *a, **k: ("2026-05-05", 1.0)
    )

    def no_second_probe(*args, **kwargs):
        raise AssertionError("the automatic date was already probed")

    monkeypatch.setattr(digitalearth, "_imagery_coverage", no_second_probe)
    raster = digitalearth.download_imagery((0, 0, 1000, 1000), style="ndwi")
    assert raster.data.shape == (100, 100, 4)
    assert len(_cached_files(tmp_path)) == 1


def test_download_imagery_ignores_files_from_older_cache_version(monkeypatch, tmp_path):
    _, calls = _stub_download(monkeypatch, tmp_path, alpha_fraction=1.0)
    old_name = "des_rgb_2026-05-05_0_0_1000_1000_10m.tif"
    (tmp_path / old_name).write_bytes(b"written by an older version")

    raster = digitalearth.download_imagery((0, 0, 1000, 1000), date="2026-05-05")

    assert raster.data.shape == (100, 100, 4)
    assert calls, "a file from an older cache version must not be reused"
    assert (
        f"des_v{digitalearth.CACHE_VERSION}_rgb_2026-05-05_0.00_0.00_1000.00_1000.00_10m.tif"
        in (_cached_files(tmp_path))
    )


def test_download_imagery_works_after_empty_cache(monkeypatch, tmp_path):
    from dtcc_core.io.data import cache

    monkeypatch.setattr(cache, "cache_dir", tmp_path)
    cache_dir = tmp_path / "downloaded_des"
    _stub_download(monkeypatch, cache_dir, alpha_fraction=1.0)
    digitalearth.download_imagery((0, 0, 1000, 1000), date="2026-05-05")

    cache.empty_cache()
    assert not cache_dir.exists()

    digitalearth.download_imagery((0, 0, 1000, 1000), date="2026-05-05")
    assert len(_cached_files(cache_dir)) == 1


def _fail_on_network(monkeypatch):
    def no_network(*args, **kwargs):
        raise AssertionError("invalid bounds must be rejected before any request")

    monkeypatch.setattr(digitalearth, "create_session", no_network)
    monkeypatch.setattr(digitalearth, "_pick_clear_date", no_network)


@pytest.mark.parametrize(
    "bbox",
    [
        (317000, 6394000, 317005, 6396000),  # 5 m wide
        (317000, 6394000, 319000, 6394009),  # 9 m tall
    ],
)
def test_download_imagery_rejects_bounds_below_one_pixel(monkeypatch, bbox):
    _fail_on_network(monkeypatch)
    with pytest.raises(ValueError, match="at least one pixel"):
        digitalearth.download_imagery(bbox)


@pytest.mark.parametrize(
    ("bbox", "resolution"),
    [
        ((0, 0, 200_000, 200_000), 10.0),  # 20000 x 20000 px
        ((317000, 6394000, 319000, 6396000), 0.01),  # 200000 x 200000 px
    ],
)
def test_download_imagery_rejects_oversized_requests(monkeypatch, bbox, resolution):
    _fail_on_network(monkeypatch)
    with pytest.raises(ValueError, match="coarser resolution or a smaller area"):
        digitalearth.download_imagery(bbox, resolution=resolution)


def test_download_imagery_rejects_non_positive_resolution(monkeypatch):
    _fail_on_network(monkeypatch)
    with pytest.raises(ValueError, match="Resolution must be positive"):
        digitalearth.download_imagery((0, 0, 1000, 1000), resolution=0)


def test_mosaic_size_accepts_the_largest_allowed_request():
    side = int(digitalearth.MAX_PIXELS**0.5)
    assert digitalearth._mosaic_size((0, 0, side * 10, side * 10), 10.0) == (side, side)


def test_download_imagery_keeps_one_pixel_wide_raster_rgba(monkeypatch, tmp_path):
    _stub_download(monkeypatch, tmp_path, alpha_fraction=1.0)
    raster = digitalearth.download_imagery((0, 0, 10, 2000), date="2026-05-05")
    assert raster.data.shape == (200, 1, 4)
    assert raster.channels == 4


def test_create_session_retries_dropped_connections_and_timeouts():
    retry = digitalearth.create_session().get_adapter(digitalearth.WMS_URL).max_retries
    assert retry.connect == digitalearth.HTTP_RETRIES
    assert retry.read == digitalearth.HTTP_RETRIES
    assert 503 in retry.status_forcelist
    assert "POST" in retry.allowed_methods


def _stub_failing_tiles(monkeypatch, tmp_path, failing, style="rgb"):
    """Stub a 2x2 tile download in which the tiles listed in 'failing' fail."""
    _stub_download(monkeypatch, tmp_path, alpha_fraction=1.0)
    calls = []

    def flaky_get_map(session, layer, style, bbox, width, height, date):
        calls.append(bbox)
        if len(calls) in failing:
            raise requests.ConnectionError(f"tile {len(calls)} dropped")
        return b""

    monkeypatch.setattr(digitalearth, "get_map", flaky_get_map)


def test_download_imagery_raises_failed_tile_by_default(monkeypatch, tmp_path):
    _stub_failing_tiles(monkeypatch, tmp_path, failing={2})
    with pytest.raises(requests.ConnectionError, match="tile 2"):
        digitalearth.download_imagery((0, 0, 6000, 6000), date="2026-05-05")
    assert _cached_files(tmp_path) == []


@pytest.mark.parametrize("style", ["rgb", "ndvi"])
def test_download_imagery_collects_failed_tiles_when_asked(
    monkeypatch, tmp_path, style
):
    _stub_failing_tiles(monkeypatch, tmp_path, failing={2})
    tile_errors = []
    raster = digitalearth.download_imagery(
        (0, 0, 6000, 6000), date="2026-05-05", style=style, tile_errors=tile_errors
    )
    assert [str(e) for e in tile_errors] == ["tile 2 dropped"]
    assert raster.data.shape == (600, 600, 4)
    # The failed tile is the top right one: columns 512 onwards of rows 0-511.
    assert not raster.data[:512, 512:, 3].any()
    assert raster.data[:512, :512, 3].all()
    assert _cached_files(tmp_path) == []


def test_download_imagery_raises_when_every_tile_fails(monkeypatch, tmp_path):
    _stub_failing_tiles(monkeypatch, tmp_path, failing={1, 2, 3, 4})
    tile_errors = []
    with pytest.raises(requests.ConnectionError, match="tile 4"):
        digitalearth.download_imagery(
            (0, 0, 6000, 6000), date="2026-05-05", tile_errors=tile_errors
        )
    assert len(tile_errors) == 4


@pytest.mark.parametrize(
    "value",
    [
        "2026-09-22",
        datetime.date(2026, 9, 22),
        datetime.datetime(2026, 9, 22, 10, 30),
    ],
)
def test_normalize_date_accepts_strings_dates_and_datetimes(value):
    assert digitalearth._normalize_date(value) == "2026-09-22"


@pytest.mark.parametrize("value", ["2026-02-30", "22/09/2026", 20260922])
def test_normalize_date_rejects_invalid_dates(value):
    with pytest.raises(ValueError, match="Invalid date"):
        digitalearth._normalize_date(value)


def test_download_imagery_accepts_date_object(monkeypatch, tmp_path):
    _stub_download(monkeypatch, tmp_path, alpha_fraction=1.0)
    digitalearth.download_imagery((0, 0, 1000, 1000), date=datetime.date(2026, 5, 5))
    assert _cached_files(tmp_path) == [
        f"des_v{digitalearth.CACHE_VERSION}_rgb_2026-05-05_0.00_0.00_1000.00_1000.00_10m.tif"
    ]


def test_download_imagery_rejects_impossible_date_before_any_request(monkeypatch):
    _fail_on_network(monkeypatch)
    with pytest.raises(ValueError, match="Invalid date"):
        digitalearth.download_imagery((0, 0, 1000, 1000), date="2026-02-30")


def test_download_imagery_writes_georeferenced_geotiff(monkeypatch, tmp_path):
    _stub_download(monkeypatch, tmp_path, alpha_fraction=1.0)
    tiles = []

    def numbered_tile(png_bytes, width, height):
        # Fill each tile with its request order, so its position can be checked.
        tiles.append(len(tiles) + 1)
        return np.full((4, height, width), tiles[-1], dtype=np.uint8)

    monkeypatch.setattr(digitalearth, "_tile_to_array", numbered_tile)
    bbox = (317000, 6394000, 323000, 6400000)  # 600 x 600 px, 2 x 2 tiles
    digitalearth.download_imagery(bbox, date="2026-05-05")

    (name,) = _cached_files(tmp_path)
    with digitalearth.rasterio.open(tmp_path / name) as src:
        assert src.crs.to_epsg() == 3006
        assert src.transform == digitalearth.from_origin(317000, 6400000, 10.0, 10.0)
        assert (src.width, src.height, src.count) == (600, 600, 4)
        data = src.read()
    # Tiles are requested row by row from the north west corner.
    assert data[0, 0, 0] == 1  # top left
    assert data[0, 0, -1] == 2  # top right
    assert data[0, -1, 0] == 3  # bottom left
    assert data[0, -1, -1] == 4  # bottom right


def test_cache_key_separates_boxes_shifted_by_less_than_a_metre(monkeypatch, tmp_path):
    _, calls = _stub_download(monkeypatch, tmp_path, alpha_fraction=1.0)
    first = digitalearth.download_imagery((0, 0, 1000, 1000), date="2026-05-05")
    shifted = digitalearth.download_imagery((0.4, 0, 1000.4, 1000), date="2026-05-05")
    assert len(calls) == 2
    assert len(_cached_files(tmp_path)) == 2
    assert shifted.georef.c == 0.4 and first.georef.c == 0.0


def _record_stac_payload(monkeypatch):
    payloads = []

    def fake_search(session, payload):
        payloads.append(payload)
        return []

    monkeypatch.setattr(digitalearth, "_stac_search_all", fake_search)
    return payloads


def test_list_available_dates_defaults_to_recent_dates(monkeypatch):
    payloads = _record_stac_payload(monkeypatch)
    digitalearth.list_available_dates(BBOX, session=object())
    start = datetime.datetime.now(datetime.UTC).date() - datetime.timedelta(
        days=digitalearth.DEFAULT_SEARCH_DAYS
    )
    assert payloads[0]["datetime"] == f"{start.isoformat()}/.."


def test_list_available_dates_accepts_date_objects(monkeypatch):
    payloads = _record_stac_payload(monkeypatch)
    digitalearth.list_available_dates(
        BBOX,
        start=datetime.date(2018, 1, 1),
        end=datetime.date(2026, 9, 30),
        session=object(),
    )
    assert payloads[0]["datetime"] == "2018-01-01/2026-09-30"


def test_tile_to_array_does_not_warn_about_georeferencing():
    import warnings

    png = make_png(4)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        digitalearth._tile_to_array(png, 4, 4)


def test_create_session_sends_dtcc_user_agent():
    session = digitalearth.create_session()
    assert session.headers["User-Agent"] == digitalearth.USER_AGENT


def test_cache_lives_in_the_shared_dtcc_data_cache():
    from dtcc_core.io.data import cache

    assert digitalearth.CACHE_DIR == str(cache.cache_dir / "downloaded_des")


def _stub_candidates(monkeypatch, coverage_by_date):
    monkeypatch.setattr(
        digitalearth,
        "list_available_dates",
        lambda *a, **k: [(date, 5.0) for date in coverage_by_date],
    )
    probed = []

    def fake_coverage(session, bbox, date):
        probed.append(date)
        return coverage_by_date[date]

    monkeypatch.setattr(digitalearth, "_imagery_coverage", fake_coverage)
    return probed


def test_latest_clear_date_prefers_older_date_covering_whole_area(monkeypatch):
    probed = _stub_candidates(
        monkeypatch, {"2026-08-24": 0.02, "2026-07-15": 1.0, "2026-06-01": 1.0}
    )
    assert digitalearth.latest_clear_date(BBOX, session=object()) == "2026-07-15"
    assert probed == ["2026-08-24", "2026-07-15"]


def test_latest_clear_date_falls_back_to_best_partial_date(monkeypatch):
    _stub_candidates(
        monkeypatch, {"2026-08-24": 0.3, "2026-07-15": 0.8, "2026-06-01": 0.0}
    )
    assert digitalearth.latest_clear_date(BBOX, session=object()) == "2026-07-15"


def test_latest_clear_date_error_counts_dates_actually_probed(monkeypatch):
    _stub_candidates(monkeypatch, {"2026-08-24": 0.0, "2026-07-15": 0.0})
    with pytest.raises(RuntimeError, match="none of the 2 most recent"):
        digitalearth.latest_clear_date(BBOX, session=object())
