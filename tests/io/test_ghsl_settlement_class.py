import numpy as np
import pytest
from affine import Affine
from pyproj import Transformer

from dtcc_core.io.data import ghsl
from dtcc_core.model import Bounds, Raster

# 3 km x 2 km in central Gothenburg (EPSG:3006): a 3 x 2 grid of 1 km cells.
BOUNDS = Bounds(318000.0, 6397000.0, 321000.0, 6399000.0)
MOLLWEIDE = "ESRI:54009"
SOURCE_NODATA = -200


def _mollweide_window(classes, margin=3000.0):
    """A 1 km Mollweide grid covering BOUNDS; ``classes(x, y)`` gives each cell's code."""
    to_moll = Transformer.from_crs("EPSG:3006", MOLLWEIDE, always_xy=True)
    xs, ys = to_moll.transform(
        [BOUNDS.xmin, BOUNDS.xmin, BOUNDS.xmax, BOUNDS.xmax],
        [BOUNDS.ymin, BOUNDS.ymax, BOUNDS.ymin, BOUNDS.ymax],
    )
    left = np.floor((min(xs) - margin) / 1000) * 1000
    top = np.ceil((max(ys) + margin) / 1000) * 1000
    cols = int(np.ceil((max(xs) + margin - left) / 1000))
    rows = int(np.ceil((top - (min(ys) - margin)) / 1000))
    transform = Affine(1000.0, 0, left, 0, -1000.0, top)
    data = np.empty((rows, cols), dtype=np.int16)
    for r in range(rows):
        for c in range(cols):
            x, y = transform * (c + 0.5, r + 0.5)
            data[r, c] = classes(x, y)
    return data, transform, MOLLWEIDE, SOURCE_NODATA


@pytest.fixture
def fake_read(monkeypatch, tmp_path):
    calls = []
    state = {"window": _mollweide_window(lambda x, y: 30)}

    def read(url, bounds):
        calls.append((url, bounds))
        data, transform, crs, nodata = state["window"]
        return data.copy(), transform, crs, nodata

    monkeypatch.setattr(ghsl, "_read_window", read)
    monkeypatch.setattr(ghsl, "cache_dir", tmp_path)
    read.calls = calls
    read.state = state
    return read


def test_returns_1km_class_raster(fake_read):
    raster = ghsl.download_settlement_class(BOUNDS)

    assert isinstance(raster, Raster)
    assert raster.data.shape == (2, 3)
    assert raster.data.dtype == np.int16
    assert raster.crs == "EPSG:3006"
    assert raster.nodata == ghsl.SETTLEMENT_CLASS_NODATA
    assert raster.georef == Affine(1000.0, 0, 318000.0, 0, -1000.0, 6399000.0)
    assert (raster.data == 30).all()


def test_classes_are_taken_from_nearest_cell_not_averaged(fake_read):
    # Split the source grid into urban centre (30) in the west and rural (11)
    # in the east; averaging would invent codes between them.
    to_moll = Transformer.from_crs("EPSG:3006", MOLLWEIDE, always_xy=True)
    split_x, _ = to_moll.transform(319500.0, 6398000.0)
    fake_read.state["window"] = _mollweide_window(lambda x, y: 30 if x < split_x else 11)

    raster = ghsl.download_settlement_class(BOUNDS)

    assert set(np.unique(raster.data)) == {11, 30}
    assert (raster.data[:, 0] == 30).all()
    assert (raster.data[:, -1] == 11).all()


def test_selects_v2_file_for_epoch(fake_read):
    ghsl.download_settlement_class(BOUNDS, epoch=1975)

    name = "GHS_SMOD_E1975_GLOBE_R2023A_54009_1000"
    url, bounds = fake_read.calls[0]
    assert url == (
        "zip+https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_SMOD_GLOBE_R2023A/"
        f"{name}/V2-0/{name}_V2_0.zip!{name}_V2_0.tif"
    )
    assert bounds == BOUNDS


def test_nodata_source_warns(fake_read, monkeypatch):
    warnings = []
    monkeypatch.setattr(ghsl, "warning", warnings.append)
    fake_read.state["window"] = _mollweide_window(lambda x, y: SOURCE_NODATA)

    raster = ghsl.download_settlement_class(BOUNDS)

    assert (raster.data == ghsl.SETTLEMENT_CLASS_NODATA).all()
    assert "No GHSL settlement class 2020 data" in warnings[0]


def test_uses_cache_on_second_call(fake_read):
    first = ghsl.download_settlement_class(BOUNDS, epoch=1990)

    cached = ghsl.download_settlement_class(BOUNDS, epoch=1990)

    assert len(fake_read.calls) == 1
    assert np.array_equal(cached.data, first.data)
    assert cached.data.dtype == np.int16


def test_every_class_has_a_name():
    assert sorted(ghsl.SETTLEMENT_CLASSES) == [10, 11, 12, 13, 21, 22, 23, 30]


def test_unsupported_epoch_raises():
    with pytest.raises(ValueError, match="Unsupported GHSL epoch"):
        ghsl.download_settlement_class(BOUNDS, epoch=2021)
