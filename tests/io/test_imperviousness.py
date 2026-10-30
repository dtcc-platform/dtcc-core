import io

import numpy as np
import pytest
import rasterio
from affine import Affine

from dtcc_core.io.data import imperviousness
from dtcc_core.model import Bounds, Raster

# 100 m x 100 m in central Gothenburg (EPSG:3006): a 10 x 10 grid of 10 m cells.
BOUNDS = Bounds(319720.0, 6397660.0, 319820.0, 6397760.0)


def _geotiff(values, crs="EPSG:3006", transform=None, colormap=False):
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
        transform=transform or Affine.identity(),
    ) as dst:
        dst.write(values)
        if colormap:
            dst.write_colormap(1, {i: (i, i, i, 255) for i in range(256)})
    return buffer.getvalue()


class FakeResponse:
    def __init__(self, content, content_type="image/tiff"):
        self.content = content
        self.headers = {"content-type": content_type}
        self.text = content.decode("latin-1")[:200]

    def raise_for_status(self):
        return None


def _size(params):
    if "size" in params:
        cols, rows = (int(v) for v in params["size"].split(","))
    else:
        cols, rows = int(params["WIDTH"]), int(params["HEIGHT"])
    return rows, cols


@pytest.fixture
def service(monkeypatch, tmp_path):
    """Record requests; each response is built by ``service.respond(params)``."""
    calls = []

    class Service:
        def respond(self, params):
            rows, cols = _size(params)
            return FakeResponse(_geotiff(np.full((rows, cols), 40)))

    svc = Service()

    def fake_get(url, params, timeout):
        calls.append((url, dict(params)))
        return svc.respond(params)

    monkeypatch.setattr(imperviousness.requests, "get", fake_get)
    monkeypatch.setattr(imperviousness, "cache_dir", tmp_path)
    svc.calls = calls
    return svc


def test_2018_returns_10m_raster_with_nodata_as_nan(service):
    values = np.full((10, 10), 30)
    values[0, :] = 255  # sea
    service.respond = lambda params: FakeResponse(_geotiff(values))

    raster = imperviousness.download_imperviousness(BOUNDS, year=2018)

    assert isinstance(raster, Raster)
    assert raster.data.shape == (10, 10)
    assert raster.crs == "EPSG:3006"
    assert raster.georef == Affine(10.0, 0, 319720.0, 0, -10.0, 6397760.0)
    assert np.isnan(raster.data[0]).all()
    assert (raster.data[1:] == 30).all()
    url, params = service.calls[0]
    assert url == imperviousness.IMPERVIOUSNESS_IMAGESERVER_URL
    assert params["bboxSR"] == 3006 and params["imageSR"] == 3006
    assert params["bbox"] == "319720.0,6397660.0,319820.0,6397760.0"
    assert params["size"] == "10,10"


def test_2018_large_area_is_tiled(service, monkeypatch):
    monkeypatch.setattr(imperviousness, "_IMAGESERVER_MAX_PIXELS", (4, 6))

    def respond(params):
        rows, cols = _size(params)
        xmin = float(params["bbox"].split(",")[0])
        # Tag each tile with its column offset so stitching can be checked.
        return FakeResponse(_geotiff(np.full((rows, cols), (xmin - BOUNDS.xmin) / 10)))

    service.respond = respond

    raster = imperviousness.download_imperviousness(BOUNDS, year=2018)

    assert len(service.calls) == 6  # 3 tile columns x 2 tile rows
    assert raster.data.shape == (10, 10)
    assert raster.data[:, 0:4].tolist() == [[0.0] * 4] * 10
    assert raster.data[:, 4:8].tolist() == [[4.0] * 4] * 10
    assert raster.data[:, 8:10].tolist() == [[8.0] * 2] * 10


@pytest.mark.parametrize("year", [2021, 2024])
def test_wms_years_reproject_palette_values(service, year):
    def respond(params):
        rows, cols = _size(params)
        ymin, xmin, ymax, xmax = (float(v) for v in params["BBOX"].split(","))
        transform = Affine((xmax - xmin) / cols, 0, xmin, 0, -(ymax - ymin) / rows, ymax)
        return FakeResponse(
            _geotiff(np.full((rows, cols), 55), crs="EPSG:3035", transform=transform, colormap=True)
        )

    service.respond = respond

    raster = imperviousness.download_imperviousness(BOUNDS, year=year)

    assert raster.data.shape == (10, 10)
    assert (raster.data == 55).all()
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
    def respond(params):
        rows, cols = _size(params)
        ymin, xmin, ymax, xmax = (float(v) for v in params["BBOX"].split(","))
        transform = Affine((xmax - xmin) / cols, 0, xmin, 0, -(ymax - ymin) / rows, ymax)
        return FakeResponse(_geotiff(np.full((rows, cols), 101), crs="EPSG:3035", transform=transform))

    service.respond = respond

    raster = imperviousness.download_imperviousness(BOUNDS, year=2024)

    assert np.isnan(raster.data).all()


def test_wms_style_change_raises(service):
    service.respond = lambda params: FakeResponse(_geotiff(np.zeros((4, 3, 3)), crs="EPSG:3035"))

    with pytest.raises(RuntimeError, match="style has changed"):
        imperviousness.download_imperviousness(BOUNDS, year=2024)


def test_non_image_response_raises(service):
    service.respond = lambda params: FakeResponse(
        b"<ServiceExceptionReport>boom</ServiceExceptionReport>", "text/xml"
    )

    with pytest.raises(RuntimeError, match="instead of an image"):
        imperviousness.download_imperviousness(BOUNDS, year=2018)


def test_download_uses_cache_on_second_call(service):
    first = imperviousness.download_imperviousness(BOUNDS, year=2018)
    call_count = len(service.calls)

    cached = imperviousness.download_imperviousness(BOUNDS, year=2018)

    assert len(service.calls) == call_count
    assert np.array_equal(cached.data, first.data)
    assert cached.georef == first.georef


def test_empty_area_warns(service, monkeypatch):
    warnings = []
    monkeypatch.setattr(imperviousness, "warning", warnings.append)
    service.respond = lambda params: FakeResponse(_geotiff(np.full(_size(params), 255)))

    imperviousness.download_imperviousness(BOUNDS, year=2018)

    assert "No imperviousness 2018 data" in warnings[0]


def test_unsupported_year_raises():
    with pytest.raises(ValueError, match="Unsupported imperviousness year"):
        imperviousness.download_imperviousness(BOUNDS, year=2015)
