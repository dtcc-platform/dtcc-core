from pathlib import Path

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
    # Uniform input stays uniform, and the total matches the exact overlap
    # total. A 100 m Mollweide cell measures about 10,029 m² in EPSG:3006
    # units here, so 2,500 m² per source cell is about 2,493 m² per target cell.
    data, transform, _ = _mollweide_window()
    assert raster.data == pytest.approx(np.full((10, 10), raster.data[5, 5]), rel=1e-4)
    assert np.nansum(raster.data) == pytest.approx(
        _exact_total(data, transform), rel=1e-4
    )


def _cells_inside_bounds(transform, shape, margin=150.0):
    """Indices of source cells whose EPSG:3006 footprint lies well inside BOUNDS."""
    to_3006 = Transformer.from_crs(MOLLWEIDE, "EPSG:3006", always_xy=True)
    cells = []
    for row in range(shape[0]):
        for col in range(shape[1]):
            corners = [
                transform @ (col + dc, row + dr) for dc in (0, 1) for dr in (0, 1)
            ]
            xs, ys = to_3006.transform(*zip(*corners, strict=True))
            if (
                min(xs) > BOUNDS.xmin + margin
                and max(xs) < BOUNDS.xmax - margin
                and min(ys) > BOUNDS.ymin + margin
                and max(ys) < BOUNDS.ymax - margin
            ):
                cells.append((row, col))
    return cells


def _exact_total(data, transform):
    """Built-up m² inside BOUNDS, splitting each source cell by overlap area."""
    from shapely.geometry import Polygon, box

    to_3006 = Transformer.from_crs(MOLLWEIDE, "EPSG:3006", always_xy=True)
    target = box(*BOUNDS.tuple)
    total = 0.0
    for row, col in zip(*np.nonzero(data), strict=True):
        ring = [
            transform @ (col + dc, row + dr)
            for dc, dr in ((0, 0), (1, 0), (1, 1), (0, 1))
        ]
        cell = Polygon([to_3006.transform(x, y) for x, y in ring])
        total += float(data[row, col]) * cell.intersection(target).area / cell.area
    return total


def test_isolated_built_up_cells_keep_their_full_area(fake_read):
    _, transform, crs = _mollweide_window()
    shape = _mollweide_window()[0].shape
    cells = _cells_inside_bounds(transform, shape)
    assert len(cells) > 20

    for row, col in cells[:: max(1, len(cells) // 8)]:
        data = np.zeros(shape, dtype=np.float32)
        data[row, col] = 10_000.0
        fake_read.state["window"] = (data, transform, crs)

        raster = ghsl.download_built_up(BOUNDS, use_cache=False)

        assert np.nansum(raster.data) == pytest.approx(10_000.0, abs=1.0)


def test_non_uniform_built_up_keeps_the_exact_total(fake_read):
    _, transform, crs = _mollweide_window()
    rng = np.random.default_rng(118)
    data = rng.uniform(0.0, 10_000.0, _mollweide_window()[0].shape).astype(np.float32)
    data[rng.random(data.shape) < 0.5] = 0.0
    fake_read.state["window"] = (data, transform, crs)

    raster = ghsl.download_built_up(BOUNDS, use_cache=False)

    assert np.nansum(raster.data) == pytest.approx(
        _exact_total(data, transform), rel=1e-3
    )


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


def test_cache_separates_grids_that_differ_by_less_than_a_metre(fake_read):
    ten_columns = ghsl.download_built_up(BOUNDS)
    eleven_columns = ghsl.download_built_up(
        Bounds(BOUNDS.xmin, BOUNDS.ymin, BOUNDS.xmax + 0.4, BOUNDS.ymax)
    )
    shifted = ghsl.download_built_up(
        Bounds(BOUNDS.xmin + 0.3, BOUNDS.ymin, BOUNDS.xmax + 0.3, BOUNDS.ymax)
    )

    assert len(fake_read.calls) == 3
    assert ten_columns.data.shape == (10, 10)
    assert eleven_columns.data.shape == (10, 11)
    assert shifted.data.shape == (10, 10)
    assert shifted.georef.c == pytest.approx(BOUNDS.xmin + 0.3)


def test_cache_file_with_another_grid_is_ignored(fake_read, monkeypatch):
    warnings = []
    monkeypatch.setattr(ghsl, "warning", warnings.append)
    transform, width, height = ghsl._target_grid(BOUNDS)
    path = ghsl._cache_path(transform, width, height, 2020, "total")
    other = Affine(100.0, 0, BOUNDS.xmin + 50.0, 0, -100.0, BOUNDS.ymax)
    ghsl._write_cache(path, np.zeros((height, width), dtype=np.float32), other)

    raster = ghsl.download_built_up(BOUNDS)

    assert len(fake_read.calls) == 1
    assert raster.georef == transform
    assert any("does not match" in message for message in warnings)


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


def test_concurrent_cache_writers_each_publish_a_complete_file(monkeypatch, tmp_path):
    import threading
    from concurrent.futures import ThreadPoolExecutor

    import rasterio

    path = tmp_path / "built_up.tif"
    transform, width, height = ghsl._target_grid(BOUNDS)
    data = np.full((height, width), 2500.0, dtype=np.float32)
    replace = ghsl.os.replace
    barrier = threading.Barrier(2)
    sources = []

    def replace_together(src, dst):
        sources.append(Path(src))
        barrier.wait(timeout=5)
        replace(src, dst)

    monkeypatch.setattr(ghsl.os, "replace", replace_together)

    with ThreadPoolExecutor(max_workers=2) as pool:
        for future in [
            pool.submit(ghsl._write_cache, path, data, transform) for _ in range(2)
        ]:
            future.result()

    assert len(set(sources)) == 2
    with rasterio.open(path) as src:
        assert src.read(1).tolist() == data.tolist()
    assert sorted(p.name for p in tmp_path.iterdir()) == ["built_up.tif"]


def test_cache_writer_accepts_a_file_published_concurrently(monkeypatch, tmp_path):
    path = tmp_path / "built_up.tif"
    transform, width, height = ghsl._target_grid(BOUNDS)
    data = np.zeros((height, width), dtype=np.float32)
    ghsl._write_cache(path, data, transform)

    def refuse(src, dst):
        raise PermissionError("file in use")

    monkeypatch.setattr(ghsl.os, "replace", refuse)

    ghsl._write_cache(path, data, transform)

    assert sorted(p.name for p in tmp_path.iterdir()) == ["built_up.tif"]


def test_cache_writer_reraises_when_nothing_was_published(monkeypatch, tmp_path):
    path = tmp_path / "built_up.tif"
    transform, width, height = ghsl._target_grid(BOUNDS)

    def refuse(src, dst):
        raise PermissionError("access denied")

    monkeypatch.setattr(ghsl.os, "replace", refuse)

    with pytest.raises(PermissionError):
        ghsl._write_cache(path, np.zeros((height, width), np.float32), transform)
    assert list(tmp_path.iterdir()) == []
