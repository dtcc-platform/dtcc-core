import numpy as np
import pytest
from affine import Affine
from pyproj import Transformer
from rasterio.errors import RasterioIOError

from dtcc_core.io.data import ghsl
from dtcc_core.model import Bounds, Raster

# 1 km x 1 km in central Gothenburg (EPSG:3006): a 10 x 10 grid of 100 m cells.
BOUNDS = Bounds(319000.0, 6397000.0, 320000.0, 6398000.0)
MOLLWEIDE = "ESRI:54009"


def _mollweide_window(value=2500.0, margin=500.0):
    """A uniform 100 m Mollweide grid covering BOUNDS, ``value`` m² per cell."""
    to_moll = Transformer.from_crs("EPSG:3006", MOLLWEIDE, always_xy=True)
    xs, ys = to_moll.transform(
        [BOUNDS.xmin, BOUNDS.xmin, BOUNDS.xmax, BOUNDS.xmax],
        [BOUNDS.ymin, BOUNDS.ymax, BOUNDS.ymin, BOUNDS.ymax],
    )
    left = np.floor((min(xs) - margin) / 100) * 100
    top = np.ceil((max(ys) + margin) / 100) * 100
    cols = int(np.ceil((max(xs) + margin - left) / 100))
    rows = int(np.ceil((top - (min(ys) - margin)) / 100))
    data = np.full((rows, cols), value, dtype=np.float32)
    return data, Affine(100.0, 0, left, 0, -100.0, top), MOLLWEIDE


@pytest.fixture
def fake_read(monkeypatch, tmp_path):
    calls = []
    state = {"window": _mollweide_window()}

    def read(url, bounds):
        calls.append((url, bounds))
        data, transform, crs = state["window"]
        return data.copy(), transform, crs

    monkeypatch.setattr(ghsl, "_read_window", read)
    monkeypatch.setattr(ghsl, "cache_dir", tmp_path)
    read.calls = calls
    read.state = state
    return read


def test_returns_100m_raster_preserving_totals(fake_read):
    raster = ghsl.download_built_up(BOUNDS)

    assert isinstance(raster, Raster)
    assert raster.data.shape == (10, 10)
    assert raster.crs == "EPSG:3006"
    assert raster.georef == Affine(100.0, 0, 319000.0, 0, -100.0, 6398000.0)
    # A uniform 25% built-up share stays 2,500 m² per 100 m cell after reprojection.
    assert raster.data == pytest.approx(np.full((10, 10), 2500.0), rel=1e-4)
    assert np.nansum(raster.data) == pytest.approx(250_000.0, rel=1e-4)


@pytest.mark.parametrize(
    "epoch, component, name",
    [
        (2020, "total", "GHS_BUILT_S_E2020_GLOBE_R2023A_54009_100"),
        (1975, "total", "GHS_BUILT_S_E1975_GLOBE_R2023A_54009_100"),
        (2030, "non_residential", "GHS_BUILT_S_NRES_E2030_GLOBE_R2023A_54009_100"),
    ],
)
def test_selects_file_for_epoch_and_component(fake_read, epoch, component, name):
    ghsl.download_built_up(BOUNDS, epoch=epoch, component=component)

    url, bounds = fake_read.calls[0]
    assert url == (
        "zip+https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/"
        f"GHS_BUILT_S_GLOBE_R2023A/{name}/V1-0/{name}_V1_0.zip!{name}_V1_0.tif"
    )
    assert bounds == BOUNDS


def test_nodata_source_gives_nan_and_warns(fake_read, monkeypatch):
    warnings = []
    monkeypatch.setattr(ghsl, "warning", warnings.append)
    data, transform, crs = _mollweide_window()
    fake_read.state["window"] = (np.full_like(data, np.nan), transform, crs)

    raster = ghsl.download_built_up(BOUNDS)

    assert np.isnan(raster.data).all()
    assert "No GHSL built-up 2020 data" in warnings[0]


def test_uses_cache_on_second_call(fake_read):
    first = ghsl.download_built_up(BOUNDS, epoch=2000)

    cached = ghsl.download_built_up(BOUNDS, epoch=2000)

    assert len(fake_read.calls) == 1
    assert np.array_equal(cached.data, first.data)
    assert cached.georef == first.georef


def test_cache_is_separate_per_epoch_and_component(fake_read):
    ghsl.download_built_up(BOUNDS, epoch=2000)
    ghsl.download_built_up(BOUNDS, epoch=2005)
    ghsl.download_built_up(BOUNDS, epoch=2005, component="non_residential")

    assert len(fake_read.calls) == 3


def test_read_failure_raises_runtime_error(monkeypatch, tmp_path):
    def fail(*args, **kwargs):
        raise RasterioIOError("HTTP response code: 404")

    monkeypatch.setattr(ghsl.rasterio, "open", fail)
    monkeypatch.setattr(ghsl, "cache_dir", tmp_path)

    with pytest.raises(RuntimeError, match="Failed to read GHSL data"):
        ghsl.download_built_up(BOUNDS)


def test_unsupported_epoch_raises():
    with pytest.raises(ValueError, match="Unsupported GHSL epoch"):
        ghsl.download_built_up(BOUNDS, epoch=2022)


def test_unsupported_component_raises():
    with pytest.raises(ValueError, match="Unsupported GHSL component"):
        ghsl.download_built_up(BOUNDS, component="residential")
