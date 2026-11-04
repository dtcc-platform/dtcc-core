import io
import json

import numpy as np
import pytest
import rasterio
import requests
from affine import Affine
from pyproj import Transformer

from dtcc_core.datasets.dataset import DatasetUpstreamError
from dtcc_core.io.data import imperviousness
from dtcc_core.model import Bounds, Raster

# 100 m x 100 m in central Gothenburg (EPSG:3006): a 10 x 10 grid of 10 m cells.
BOUNDS = Bounds(319720.0, 6397660.0, 319820.0, 6397760.0)
TO_3035 = Transformer.from_crs("EPSG:3006", "EPSG:3035", always_xy=True)


def _palette():
    colormap = {i: (i, i, i, 255) for i in range(256)}
    for index, rgb in imperviousness.IMPERVIOUSNESS_PALETTE_CHECKS.items():
        colormap[index] = (*rgb, 255)
    return colormap


def _geotiff(values, crs, transform, colormap=None):
    """Return GeoTIFF bytes for a (bands, rows, cols) or (rows, cols) uint8 array."""
    values = np.asarray(values, dtype=np.uint8)
    if values.ndim == 2:
        values = values[np.newaxis]
    buffer = io.BytesIO()
    with rasterio.open(
        buffer,
        "w",
        driver="GTiff",
        width=values.shape[2],
        height=values.shape[1],
        count=values.shape[0],
        dtype="uint8",
        crs=crs,
        transform=transform,
    ) as dst:
        dst.write(values)
        if colormap:
            dst.write_colormap(1, colormap)
    return buffer.getvalue()


def _gradient(transform, rows, cols):
    """Values 0-99 that differ between neighbouring 10 m cells in both axes."""
    col_ids, row_ids = np.meshgrid(np.arange(cols), np.arange(rows))
    xs, ys = transform * (col_ids + 0.5, row_ids + 0.5)
    return ((np.floor(xs / 10) * 7 + np.floor(ys / 10) * 3) % 100).astype(np.uint8)


def _request_grid(params):
    """Return the transform, rows, cols and EPSG code a request asks for."""
    if "size" in params:
        cols, rows = (int(v) for v in params["size"].split(","))
        xmin, ymin, xmax, ymax = (float(v) for v in params["bbox"].split(","))
        epsg = 3006
    else:
        cols, rows = int(params["WIDTH"]), int(params["HEIGHT"])
        ymin, xmin, ymax, xmax = (float(v) for v in params["BBOX"].split(","))
        epsg = 3035
    transform = Affine((xmax - xmin) / cols, 0, xmin, 0, -(ymax - ymin) / rows, ymax)
    return transform, rows, cols, epsg


class FakeResponse:
    def __init__(self, content, content_type="image/tiff", status=200):
        self.content = content
        self.headers = {"content-type": content_type}
        self.status_code = status

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError(f"{self.status_code} error", response=self)

    def iter_content(self, chunk_size):
        for start in range(0, len(self.content), chunk_size):
            yield self.content[start : start + chunk_size]


def image_response(params, values=None, colormap=None):
    """A correctly georeferenced response; ``values`` defaults to a gradient."""
    transform, rows, cols, epsg = _request_grid(params)
    if values is None:
        values = _gradient(transform, rows, cols)
    if colormap is None and epsg == 3035:
        colormap = _palette()
    return FakeResponse(_geotiff(values, f"EPSG:{epsg}", transform, colormap or None))


@pytest.fixture
def service(monkeypatch, tmp_path):
    """Record requests; each response is built by ``service.respond(params)``."""
    calls = []

    class Service:
        respond = staticmethod(image_response)

    svc = Service()

    class FakeSession:
        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

        def get(self, url, params, timeout, stream):
            calls.append((url, dict(params)))
            result = svc.respond(params)
            if isinstance(result, Exception):
                raise result
            return result

    monkeypatch.setattr(imperviousness, "_session", FakeSession)
    monkeypatch.setattr(imperviousness, "cache_dir", tmp_path)
    svc.calls = calls
    return svc


def test_2018_returns_gradient_on_the_requested_grid(service):
    raster = imperviousness.download_imperviousness(BOUNDS, year=2018)

    assert isinstance(raster, Raster)
    assert raster.crs == "EPSG:3006"
    assert raster.georef == Affine(10.0, 0, 319720.0, 0, -10.0, 6397760.0)
    assert raster.data.tolist() == _gradient(raster.georef, 10, 10).tolist()
    url, params = service.calls[0]
    assert url == imperviousness.IMPERVIOUSNESS_IMAGESERVER_URL
    assert params["bboxSR"] == 3006 and params["imageSR"] == 3006
    assert params["bbox"] == "319720.0,6397660.0,319820.0,6397760.0"
    assert params["size"] == "10,10"
    assert params["noData"] == 255


def test_2018_nodata_is_nan(service):
    def respond(params):
        values = np.full((10, 10), 30)
        values[0, :] = 255  # sea or outside the product
        return image_response(params, values)

    service.respond = respond

    raster = imperviousness.download_imperviousness(BOUNDS, year=2018)

    assert np.isnan(raster.data[0]).all()
    assert (raster.data[1:] == 30).all()


def test_2018_large_area_is_tiled_and_stitched(service, monkeypatch):
    monkeypatch.setattr(imperviousness, "_IMAGESERVER_MAX_PIXELS", (4, 6))

    raster = imperviousness.download_imperviousness(BOUNDS, year=2018)

    assert len(service.calls) == 6  # 3 tile columns x 2 tile rows
    assert raster.data.tolist() == _gradient(raster.georef, 10, 10).tolist()


@pytest.mark.parametrize("year", [2021, 2024])
def test_wms_years_reproject_gradient(service, year):
    raster = imperviousness.download_imperviousness(BOUNDS, year=year)

    # Each EPSG:3006 cell should hold the EPSG:3035 source cell under its centre.
    col_ids, row_ids = np.meshgrid(np.arange(10), np.arange(10))
    xs, ys = raster.georef * (col_ids + 0.5, row_ids + 0.5)
    sx, sy = TO_3035.transform(xs, ys)
    expected = (np.floor(sx / 10) * 7 + np.floor(sy / 10) * 3) % 100
    # Nearest-neighbour warping may pick a neighbour for centres on a cell edge.
    assert np.mean(raster.data == expected) > 0.9

    url, params = service.calls[0]
    workspace, layer = imperviousness.IMPERVIOUSNESS_GEOSERVER_LAYERS[year]
    assert workspace in url
    assert params["LAYERS"] == layer
    assert params["CRS"] == "EPSG:3035"
    # WMS 1.3.0 EPSG:3035 axis order is northing first; the request is snapped
    # to the product's 10 m grid.
    ymin, xmin, ymax, xmax = (float(v) for v in params["BBOX"].split(","))
    assert 3.8e6 < ymin < ymax < 3.9e6 and 4.4e6 < xmin < xmax < 4.5e6
    assert all(v % 10 == 0 for v in (ymin, xmin, ymax, xmax))


def test_wms_value_above_100_is_nan(service):
    service.respond = lambda params: image_response(
        params, np.full(_request_grid(params)[1:3], 101)
    )

    raster = imperviousness.download_imperviousness(BOUNDS, year=2024)

    assert np.isnan(raster.data).all()


def test_wms_palette_change_raises(service):
    restyled = _palette()
    restyled[50] = (0, 0, 255, 255)
    service.respond = lambda params: image_response(params, colormap=restyled)

    with pytest.raises(DatasetUpstreamError, match="map style has changed") as info:
        imperviousness.download_imperviousness(BOUNDS, year=2024)
    assert info.value.failure_class == "invalid_payload"


def test_wms_missing_palette_raises(service):
    def respond(params):
        transform, rows, cols, _ = _request_grid(params)
        return FakeResponse(_geotiff(np.zeros((rows, cols)), "EPSG:3035", transform))

    service.respond = respond

    with pytest.raises(DatasetUpstreamError, match="map style has changed"):
        imperviousness.download_imperviousness(BOUNDS, year=2021)


@pytest.mark.parametrize(
    ("change", "message"),
    [
        ("size", "size 1x1"),
        ("bands", "4 bands"),
        ("crs", "CRS"),
        ("offset", "grid"),
    ],
)
def test_mismatched_image_raises_before_reading(service, change, message):
    def respond(params):
        transform, rows, cols, epsg = _request_grid(params)
        values = np.full((rows, cols), 40)
        crs = f"EPSG:{epsg}"
        if change == "size":
            values = np.full((1, 1), 40)
        elif change == "bands":
            values = np.zeros((4, rows, cols))
        elif change == "crs":
            crs = "EPSG:4326"
        elif change == "offset":
            transform = transform * Affine.translation(0.5, 0)
        return FakeResponse(_geotiff(values, crs, transform))

    service.respond = respond

    with pytest.raises(DatasetUpstreamError, match=message) as info:
        imperviousness.download_imperviousness(BOUNDS, year=2018)
    assert info.value.failure_class == "invalid_payload"
    assert not list(imperviousness._cache_path(BOUNDS, 2018).parent.glob("*.tif"))


def test_oversized_response_is_refused(service, monkeypatch):
    monkeypatch.setattr(imperviousness, "_RESPONSE_HEADROOM_BYTES", 0)
    service.respond = lambda params: FakeResponse(b"\0" * 1000)

    with pytest.raises(DatasetUpstreamError, match="exceeds 400 bytes"):
        imperviousness.download_imperviousness(BOUNDS, year=2018)


def test_arcgis_json_error_in_image_response(service):
    body = {
        "error": {
            "code": 400,
            "message": "Invalid or missing input parameters.",
            "details": ["The requested image exceeds the size limit."],
        }
    }
    service.respond = lambda params: FakeResponse(json.dumps(body).encode())

    with pytest.raises(DatasetUpstreamError, match="exceeds the size limit") as info:
        imperviousness.download_imperviousness(BOUNDS, year=2018)
    assert info.value.failure_class == "http_4xx"
    assert info.value.status_code == 400


@pytest.mark.parametrize(
    ("error", "failure_class", "message"),
    [
        (
            {"code": 400, "message": "Invalid bbox", "details": "bad extent"},
            "http_4xx",
            "Invalid bbox bad extent",
        ),
        ({"code": 500, "message": "Server error"}, "http_5xx", "Server error"),
        ({"message": "No code"}, "invalid_payload", "No code"),
    ],
)
def test_arcgis_json_error_classification(service, error, failure_class, message):
    service.respond = lambda params: FakeResponse(json.dumps({"error": error}).encode())

    with pytest.raises(DatasetUpstreamError, match=message) as info:
        imperviousness.download_imperviousness(BOUNDS, year=2018)
    assert info.value.failure_class == failure_class


@pytest.mark.parametrize("year", [2018, 2024])
def test_unreadable_image_is_a_typed_error(service, year):
    service.respond = lambda params: FakeResponse(b"II*\0truncated", "image/tiff")

    with pytest.raises(DatasetUpstreamError, match="unreadable image") as info:
        imperviousness.download_imperviousness(BOUNDS, year=year)
    assert info.value.failure_class == "invalid_payload"
    assert not imperviousness._cache_path(BOUNDS, year).exists()


def test_session_ignores_retry_after():
    session = imperviousness._session()
    retry = session.get_adapter("https://example.test").max_retries

    assert retry.respect_retry_after_header is False
    assert 503 in retry.status_forcelist


def test_cache_write_uses_a_unique_temporary_file(tmp_path, monkeypatch):
    seen = []
    real_replace = imperviousness.os.replace

    def record(src, dst):
        seen.append(src)
        real_replace(src, dst)

    monkeypatch.setattr(imperviousness.os, "replace", record)
    path = tmp_path / "cache.tif"
    data = np.zeros((2, 2), dtype=np.float32)

    imperviousness._write_cache(path, data, Affine(10, 0, 0, 0, -10, 20))
    imperviousness._write_cache(path, data, Affine(10, 0, 0, 0, -10, 20))

    assert len(set(seen)) == 2
    assert path.exists() and list(tmp_path.iterdir()) == [path]


def test_non_image_response_raises(service):
    service.respond = lambda params: FakeResponse(
        b"<ServiceExceptionReport>boom</ServiceExceptionReport>", "text/xml"
    )

    with pytest.raises(DatasetUpstreamError, match="instead of an image"):
        imperviousness.download_imperviousness(BOUNDS, year=2024)


@pytest.mark.parametrize(
    ("result", "failure_class"),
    [
        (requests.ConnectTimeout("slow"), "timeout"),
        (requests.ConnectionError("down"), "connection"),
        (FakeResponse(b"", status=503), "http_5xx"),
        (FakeResponse(b"", status=404), "http_4xx"),
    ],
)
def test_request_failures_are_typed_and_not_cached(service, result, failure_class):
    service.respond = lambda params: result

    with pytest.raises(DatasetUpstreamError) as info:
        imperviousness.download_imperviousness(BOUNDS, year=2018)

    assert info.value.failure_class == failure_class
    assert not imperviousness._cache_path(BOUNDS, 2018).exists()


def test_download_uses_cache_on_second_call(service):
    first = imperviousness.download_imperviousness(BOUNDS, year=2018)
    call_count = len(service.calls)

    cached = imperviousness.download_imperviousness(BOUNDS, year=2018)

    assert len(service.calls) == call_count
    assert np.array_equal(cached.data, first.data)
    assert cached.georef == first.georef


def test_cache_is_keyed_by_exact_bounds(service):
    shifted = Bounds(BOUNDS.xmin + 0.4, BOUNDS.ymin, BOUNDS.xmax + 0.4, BOUNDS.ymax)

    imperviousness.download_imperviousness(BOUNDS, year=2018)
    imperviousness.download_imperviousness(shifted, year=2018)

    assert len(service.calls) == 2
    assert imperviousness._cache_path(BOUNDS, 2018) != imperviousness._cache_path(
        shifted, 2018
    )


def test_wider_bounds_get_their_own_grid(service):
    wider = Bounds(BOUNDS.xmin, BOUNDS.ymin, BOUNDS.xmax + 0.4, BOUNDS.ymax)

    imperviousness.download_imperviousness(BOUNDS, year=2018)
    raster = imperviousness.download_imperviousness(wider, year=2018)

    assert raster.data.shape == (10, 11)


def test_empty_area_warns(service, monkeypatch):
    warnings = []
    monkeypatch.setattr(imperviousness, "warning", warnings.append)
    service.respond = lambda params: image_response(
        params, np.full(_request_grid(params)[1:3], 255)
    )

    imperviousness.download_imperviousness(BOUNDS, year=2018)

    assert "No imperviousness 2018 data" in warnings[0]


def test_area_limit(service):
    huge = Bounds(0.0, 0.0, 60_000.0, 60_000.0)

    with pytest.raises(ValueError, match="limit"):
        imperviousness.download_imperviousness(huge, year=2018)
    assert service.calls == []


def test_progress_is_reported_per_tile(service, monkeypatch):
    reports = []
    monkeypatch.setattr(imperviousness, "_IMAGESERVER_MAX_PIXELS", (5, 5))
    monkeypatch.setattr(
        imperviousness, "report_progress", lambda **kwargs: reports.append(kwargs)
    )

    imperviousness.download_imperviousness(BOUNDS, year=2018)

    assert [(r["current"], r["total"]) for r in reports] == [
        (1, 4),
        (2, 4),
        (3, 4),
        (4, 4),
    ]


def test_empty_raster_matches_grid():
    raster = imperviousness.empty_imperviousness_raster(BOUNDS)

    assert raster.data.shape == (10, 10)
    assert np.isnan(raster.data).all()
    assert raster.georef == Affine(10.0, 0, 319720.0, 0, -10.0, 6397760.0)


def test_unsupported_year_raises():
    with pytest.raises(ValueError, match="Unsupported imperviousness year"):
        imperviousness.download_imperviousness(BOUNDS, year=2015)


# The CRS the 2018 ImageServer writes: EPSG:3006 with a SWEREF99 datum code
# that GDAL only matches at reduced confidence.
ARCGIS_SWEREF99_TM = (
    'PROJCS["SWEREF99 TM",GEOGCS["SWEREF99",DATUM["SWEREF99",SPHEROID["GRS 1980",'
    '6378137,298.257222101004,AUTHORITY["EPSG","7019"]],AUTHORITY["EPSG","6619"]],'
    'PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],'
    'AUTHORITY["EPSG","4619"]],PROJECTION["Transverse_Mercator"],'
    'PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",15],'
    'PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],'
    'PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],'
    'AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","3006"]]'
)


def test_image_check_accepts_the_imageserver_crs():
    from types import SimpleNamespace

    from rasterio.crs import CRS

    transform = Affine(10, 0, 319720, 0, -10, 6397760)
    src = SimpleNamespace(
        width=10,
        height=10,
        count=1,
        dtypes=("uint8",),
        crs=CRS.from_wkt(ARCGIS_SWEREF99_TM),
        transform=transform,
    )

    imperviousness._check_image(src, "url", 2018, 10, 10, transform, 3006)


def test_all_zero_wms_result_warns(service, monkeypatch):
    warnings = []
    monkeypatch.setattr(imperviousness, "warning", warnings.append)
    service.respond = lambda params: image_response(
        params, np.zeros(_request_grid(params)[1:3])
    )

    imperviousness.download_imperviousness(BOUNDS, year=2024)

    assert "reads 0 %" in warnings[0]
