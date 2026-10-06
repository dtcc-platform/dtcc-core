import math
import os
import struct
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pytest
import rasterio
import requests
from affine import Affine
from rasterio.enums import ColorInterp
from rasterio.env import get_gdal_config
from rasterio.errors import NotGeoreferencedWarning
from rasterio.transform import from_bounds

from dtcc_core.io import raster_tiles
from dtcc_core.io.raster_tiles import (
    DEFAULT_MAX_MEMORY_BYTES,
    LOAD_LAYOUTS,
    READ_OVERHEAD_BYTES,
    RasterTileLayoutError,
    RasterTileReadError,
    WorkingMemoryError,
    load_tile,
)
from dtcc_core.model import Bounds
from dtcc_core.model import exchange
from dtcc_core.model.values.raster_tiles import RasterTile

CELL = (672500.0, 6577500.0, 672502.0, 6577502.0)
FULL_CELL = (672500.0, 6577500.0, 675000.0, 6580000.0)
LM_COLORS = [ColorInterp.red, ColorInterp.green, ColorInterp.blue]
BANDS = {"rgb": 3, "rgbi": 4, "cir": 3}
# Fixtures write and check files with the real open, unseen by recording spies.
RASTERIO_OPEN = rasterio.open


def band_values(bands, height=4, width=4, dtype="uint8"):
    data = np.arange(1, bands * height * width + 1).reshape(bands, height, width)
    data = (data * 7 % 251).astype(dtype)
    if bands == 4:
        data[3, 0, :] = 0
    return data


def write_tiff(
    path,
    data,
    *,
    bounds=CELL,
    crs="EPSG:3006",
    nodata=None,
    transform=None,
    colorinterp=None,
    mask=False,
    band_masks=False,
    **profile,
):
    """Write a GeoTIFF labelled like the LM originals: R, G, B, then undefined.

    ``mask`` adds an internal per-dataset mask; ``band_masks`` adds a ``.msk``
    sidecar holding one mask per band, which GDAL reads with the file.
    """
    bands, height, width = data.shape
    colors = (colorinterp or LM_COLORS + [ColorInterp.undefined] * bands)[:bands]
    if bands > 3 and ColorInterp.alpha not in colors:
        profile.setdefault("alpha", "unspecified")
    path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.Env(GDAL_TIFF_INTERNAL_MASK=True):
        with RASTERIO_OPEN(
            path,
            "w",
            driver="GTiff",
            width=width,
            height=height,
            count=bands,
            dtype=data.dtype,
            crs=crs,
            transform=transform or from_bounds(*bounds, width, height),
            nodata=nodata,
            **profile,
        ) as dst:
            dst.write(data)
            dst.colorinterp = colors
            if mask:
                dst.write_mask(np.full((height, width), 255, dtype=np.uint8))
    if band_masks:
        masks = np.full((bands, height, width), 255, dtype=np.uint8)
        masks[:, 0, 0] = 0
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", NotGeoreferencedWarning)
            with RASTERIO_OPEN(
                f"{path}.msk",
                "w",
                driver="GTiff",
                width=width,
                height=height,
                count=bands,
                dtype="uint8",
            ) as dst:
                dst.write(masks)
                dst.update_tags(
                    **{f"INTERNAL_MASK_FLAGS_{b}": "0" for b in range(1, bands + 1)}
                )
    with RASTERIO_OPEN(path) as src:
        assert list(src.colorinterp) == colors
    return path


def make_tile(path, spektraltyp="rgb", extent=CELL, crs="EPSG:3006", size_bytes=None):
    return RasterTile(
        path=path,
        id="o65775_6725_25_mr25",
        collection="orto-o2-2025",
        datetime=datetime(2025, 5, 31, 8, 10, 56, tzinfo=timezone.utc),
        extent=Bounds(*extent),
        crs=crs,
        spektraltyp=spektraltyp,
        resolution=0.5,
        size_bytes=size_bytes,
    )


def written_tile(tmp_path, spektraltyp="rgb", bands=None, tile_kwargs=None, **kwargs):
    data = band_values(bands or BANDS[spektraltyp])
    path = write_tiff(tmp_path / "tile.tif", data, **kwargs)
    return make_tile(path, spektraltyp, **(tile_kwargs or {})), data


def load(tile, budget=DEFAULT_MAX_MEMORY_BYTES):
    return load_tile(tile, max_memory_bytes=budget)


def estimate(height, width, bands):
    return height * width * bands + height * width + READ_OVERHEAD_BYTES


@pytest.fixture
def no_band_reads(monkeypatch):
    def forbidden(self, *args, **kwargs):
        raise AssertionError("pixels must not be read")

    monkeypatch.setattr(rasterio.io.DatasetReader, "read", forbidden)


def forbid_allocation(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("arrays must not be allocated")

    monkeypatch.setattr(np, "empty", forbidden)
    monkeypatch.setattr(np, "zeros", forbidden)


@pytest.fixture
def opened(monkeypatch):
    """Record arguments, GDAL cache limit and driver of every Rasterio open."""
    calls = []

    def recording_open(*args, **kwargs):
        call = {
            "args": args,
            "kwargs": kwargs,
            "cache": get_gdal_config("GDAL_CACHEMAX"),
        }
        calls.append(call)
        dataset = RASTERIO_OPEN(*args, **kwargs)
        call["driver"] = dataset.driver
        return dataset

    monkeypatch.setattr(rasterio, "open", recording_open)
    return calls


# Constants


def test_constants():
    assert DEFAULT_MAX_MEMORY_BYTES == 2 * 1024**3
    assert READ_OVERHEAD_BYTES == 64 * 1024**2
    assert LOAD_LAYOUTS == {"rgb": 3, "rgbi": 4, "cir": 3}


# Loading


@pytest.mark.parametrize("spektraltyp", ["rgb", "rgbi", "cir"])
@pytest.mark.parametrize("nodata", [0, None])
def test_loads_every_band_in_hwc_order(tmp_path, spektraltyp, nodata):
    tile, data = written_tile(tmp_path, spektraltyp, nodata=nodata)
    raster = load(tile)
    assert raster.data.shape == (4, 4, BANDS[spektraltyp])
    assert raster.data.dtype == np.uint8
    assert raster.data.flags.c_contiguous
    np.testing.assert_array_equal(raster.data, data.transpose(1, 2, 0))
    with rasterio.open(tile.path) as src:
        assert raster.georef == src.transform
    assert raster.crs == "EPSG:3006"
    if nodata is None:
        assert math.isnan(raster.nodata)
    else:
        assert raster.nodata == 0
    exchange.validate(raster)


def test_zero_samples_in_fourth_band_are_kept(tmp_path):
    tile, data = written_tile(tmp_path, "rgbi", nodata=0)
    raster = load(tile)
    assert (raster.data[0, :, 3] == 0).all()
    assert (raster.data[1:, :, 3] == data[3, 1:, :]).all()


@pytest.mark.parametrize("height, width", [(1, 5), (5, 1), (1, 1)])
def test_spatial_singleton_axes_are_kept(tmp_path, height, width):
    data = band_values(3, height, width)
    path = write_tiff(tmp_path / "tile.tif", data)
    raster = load(make_tile(path))
    assert raster.data.shape == (height, width, 3)
    np.testing.assert_array_equal(raster.data, data.transpose(1, 2, 0))
    exchange.validate(raster)


def test_loading_is_local_and_leaves_the_file_unchanged(tmp_path, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("loading must not use the network")

    monkeypatch.setattr(requests, "get", forbidden)
    monkeypatch.setattr(requests.Session, "request", forbidden)
    tile, _ = written_tile(tmp_path, "rgbi", nodata=0)
    before = (tile.path.read_bytes(), tile.path.stat().st_mtime_ns)
    load(tile)
    assert (tile.path.read_bytes(), tile.path.stat().st_mtime_ns) == before


# Header admission


def test_pan_loading_is_not_supported(tmp_path, opened):
    tile, _ = written_tile(tmp_path, "pan", bands=1, colorinterp=[ColorInterp.gray])
    with pytest.raises(NotImplementedError) as info:
        load(tile)
    assert str(tile.path) in str(info.value)
    assert tile.id in str(info.value)
    assert opened == []


REFUSALS = {
    "unknown spektraltyp": (dict(spektraltyp="nir", bands=3), "spektraltyp"),
    "uint16 rgb": (dict(spektraltyp="rgb", dtype="uint16"), "uint8"),
    "rgb with 4 bands": (dict(spektraltyp="rgb", bands=4), "bands"),
    "rgbi with 3 bands": (dict(spektraltyp="rgbi", bands=3), "bands"),
    "cir with 4 bands": (dict(spektraltyp="cir", bands=4), "bands"),
    "fourth band alpha": (
        dict(spektraltyp="rgbi", colorinterp=LM_COLORS + [ColorInterp.alpha]),
        "alpha",
    ),
    "internal mask": (dict(spektraltyp="rgb", mask=True), "mask"),
    "per-band masks": (dict(spektraltyp="rgb", band_masks=True), "mask"),
    "other CRS": (dict(spektraltyp="rgb", crs="EPSG:3021"), "CRS"),
    "no CRS": (dict(spektraltyp="rgb", crs=None), "CRS"),
    "tile in other CRS": (
        dict(spektraltyp="rgb", tile_kwargs={"crs": "EPSG:3021"}),
        "CRS",
    ),
    "wrong cell": (
        dict(spektraltyp="rgb", bounds=tuple(v + 2500.0 for v in CELL)),
        "bounds",
    ),
    "NaN georeferencing": (
        dict(
            spektraltyp="rgb",
            transform=Affine(0.5, 0.0, math.nan, 0.0, -0.5, CELL[3]),
        ),
        "bounds",
    ),
}


@pytest.mark.parametrize("case", REFUSALS, ids=list(REFUSALS))
def test_unadmitted_layouts_are_refused_before_pixel_reads(
    tmp_path, no_band_reads, case
):
    kwargs, reason = REFUSALS[case]
    kwargs = dict(kwargs)
    spektraltyp = kwargs.pop("spektraltyp")
    bands = kwargs.pop("bands", None)
    dtype = kwargs.pop("dtype", "uint8")
    tile_kwargs = kwargs.pop("tile_kwargs", {})
    data = band_values(bands or BANDS.get(spektraltyp, 3), dtype=dtype)
    path = write_tiff(tmp_path / "tile.tif", data, **kwargs)
    tile = make_tile(path, spektraltyp, **tile_kwargs)
    with pytest.raises(RasterTileLayoutError) as info:
        load(tile)
    assert info.value.path == path
    assert info.value.spektraltyp == spektraltyp
    assert reason in info.value.reason
    assert str(path) in str(info.value)


@pytest.mark.parametrize(
    "nodatavals, reason",
    [
        ((0.0, 0.0, 1.0), "nodata differs between bands"),
        ((0.0, None, 0.0), "nodata declared on some bands only"),
        ((None, None, 0.0), "nodata declared on some bands only"),
    ],
)
def test_nodata_must_be_common_to_every_band(
    tmp_path, monkeypatch, no_band_reads, nodatavals, reason
):
    # GeoTIFF stores one nodata value per file, so the per-band declarations are
    # set on the header of a GDAL-written file.
    tile, _ = written_tile(tmp_path, "rgb", nodata=0)
    real_open = rasterio.open

    class Header:
        def __init__(self, dataset):
            self._dataset = dataset

        def __getattr__(self, name):
            return getattr(self._dataset, name)

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            self._dataset.close()

    Header.nodatavals = nodatavals
    monkeypatch.setattr(
        rasterio, "open", lambda *args, **kwargs: Header(real_open(*args, **kwargs))
    )
    with pytest.raises(RasterTileLayoutError) as info:
        load(tile)
    assert reason in info.value.reason


@pytest.mark.parametrize("edge", range(4))
@pytest.mark.parametrize("offset, accepted", [(0.0005, True), (0.002, False)])
def test_bounds_tolerance_is_absolute_per_edge(tmp_path, edge, offset, accepted):
    bounds = list(CELL)
    bounds[edge] += offset if edge >= 2 else -offset
    tile, _ = written_tile(tmp_path, "rgb", bounds=tuple(bounds))
    if accepted:
        assert load(tile).data.shape == (4, 4, 3)
    else:
        with pytest.raises(RasterTileLayoutError, match="bounds"):
            load(tile)


# Local and source failures


def test_missing_file_raises_file_not_found(tmp_path, opened):
    tile = make_tile(tmp_path / "missing.tif")
    with pytest.raises(FileNotFoundError) as info:
        load(tile)
    assert not isinstance(info.value, RasterTileReadError)
    assert opened == []


def test_directory_raises_is_a_directory(tmp_path, opened):
    path = tmp_path / "tile.tif"
    path.mkdir()
    with pytest.raises(IsADirectoryError) as info:
        load(make_tile(path))
    assert not isinstance(info.value, RasterTileReadError)
    assert opened == []


@pytest.mark.skipif(os.geteuid() == 0, reason="root ignores file permissions")
def test_unreadable_file_raises_permission_error(tmp_path, opened):
    tile, _ = written_tile(tmp_path)
    tile.path.chmod(0)
    try:
        with pytest.raises(PermissionError) as info:
            load(tile)
    finally:
        tile.path.chmod(0o644)
    assert not isinstance(info.value, RasterTileReadError)
    assert opened == []


def test_rasterio_gets_the_path_and_only_the_geotiff_driver(tmp_path, opened):
    tile, _ = written_tile(tmp_path)
    load(tile)
    assert len(opened) == 1
    (argument,) = opened[0]["args"]
    assert isinstance(argument, (str, Path)) and Path(argument) == tile.path
    assert opened[0]["kwargs"] == {"driver": "GTiff"}


def test_non_raster_file_is_a_read_error(tmp_path):
    path = tmp_path / "tile.tif"
    path.write_text("<html><body>Bad gateway</body></html>")
    with pytest.raises(RasterTileReadError) as info:
        load(make_tile(path))
    assert info.value.path == path
    assert str(path) in str(info.value)
    assert isinstance(info.value.__cause__, rasterio.errors.RasterioError)


def test_vrt_file_is_a_read_error_and_never_opened_as_vrt(tmp_path, opened):
    source = write_tiff(tmp_path / "src" / "valid.tif", band_values(3))
    path = tmp_path / "tile.tif"
    path.write_text(
        '<VRTDataset rasterXSize="4" rasterYSize="4"><SRS>EPSG:3006</SRS>'
        f"<GeoTransform>{CELL[0]},0.5,0,{CELL[3]},0,-0.5</GeoTransform>"
        + "".join(
            f'<VRTRasterBand dataType="Byte" band="{b}"><SimpleSource>'
            f'<SourceFilename relativeToVRT="0">{source}</SourceFilename>'
            f"<SourceBand>{b}</SourceBand></SimpleSource></VRTRasterBand>"
            for b in (1, 2, 3)
        )
        + "</VRTDataset>"
    )
    with pytest.raises(RasterTileReadError) as info:
        load(make_tile(path))
    assert info.value.path == path
    assert all(call.get("driver") != "VRT" for call in opened)


def damage_first_block(path):
    data = bytearray(path.read_bytes())
    assert data[:4] == b"II*\0"
    ifd = struct.unpack_from("<I", data, 4)[0]
    for entry in range(struct.unpack_from("<H", data, ifd)[0]):
        entry_offset = ifd + 2 + 12 * entry
        tag, kind, count, value = struct.unpack_from("<HHII", data, entry_offset)
        if tag == 324:  # TileOffsets
            assert kind == 4 and count > 1
            offset = struct.unpack_from("<I", data, value)[0]
    data[offset : offset + 32] = b"\xff" * 32
    path.write_bytes(data)


def test_damaged_block_is_a_read_error(tmp_path):
    data = np.random.default_rng(0).integers(0, 255, (3, 64, 64), dtype=np.uint8)
    path = write_tiff(
        tmp_path / "tile.tif",
        data,
        tiled=True,
        blockxsize=16,
        blockysize=16,
        compress="deflate",
    )
    damage_first_block(path)
    with pytest.raises(RasterTileReadError) as info:
        load(make_tile(path))
    assert info.value.path == path
    assert str(path) in str(info.value)
    assert isinstance(info.value.__cause__, rasterio.errors.RasterioError)


def test_memory_error_from_allocation_propagates(tmp_path, monkeypatch):
    tile, _ = written_tile(tmp_path)
    real_empty = np.empty

    def failing_empty(shape, *args, **kwargs):
        if tuple(np.atleast_1d(shape)) == (4, 4, 3):
            raise MemoryError("no memory")
        return real_empty(shape, *args, **kwargs)

    monkeypatch.setattr(np, "empty", failing_empty)
    with pytest.raises(MemoryError, match="no memory"):
        load(tile)


def test_unexpected_read_exception_propagates(tmp_path, monkeypatch):
    tile, _ = written_tile(tmp_path)

    def failing_read(self, *args, **kwargs):
        raise RuntimeError("unexpected")

    monkeypatch.setattr(rasterio.io.DatasetReader, "read", failing_read)
    with pytest.raises(RuntimeError, match="unexpected") as info:
        load(tile)
    assert type(info.value) is RuntimeError


# Budget


@pytest.mark.parametrize(
    "budget, error",
    [(True, TypeError), (1.5, TypeError), ("2", TypeError), (0, ValueError),
     (-1, ValueError)],
)
def test_invalid_budget_is_rejected_before_file_access(tmp_path, opened, budget, error):
    with pytest.raises(error):
        load(make_tile(tmp_path / "missing.tif"), budget)
    assert opened == []


def test_full_tile_estimate_is_refused_before_allocation(
    tmp_path, monkeypatch, no_band_reads
):
    path = tmp_path / "full.tif"
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        width=15625,
        height=15625,
        count=4,
        dtype="uint8",
        crs="EPSG:3006",
        transform=from_bounds(*FULL_CELL, 15625, 15625),
        tiled=True,
        blockxsize=512,
        blockysize=512,
        sparse_ok=True,
        alpha="unspecified",
    ) as dst:
        dst.colorinterp = LM_COLORS + [ColorInterp.undefined]
    assert path.stat().st_size < 1024**2
    tile = make_tile(path, "rgbi", extent=FULL_CELL, size_bytes=1)
    forbid_allocation(monkeypatch)
    with pytest.raises(WorkingMemoryError) as info:
        load(tile, 1_287_811_989 - 1)
    assert info.value.estimate_bytes == 1_287_811_989
    assert info.value.budget_bytes == 1_287_811_989 - 1
    assert "1287811989" in str(info.value) and "1287811988" in str(info.value)
    assert "max_memory_bytes" in str(info.value)


def test_budget_equal_to_the_estimate_is_accepted(tmp_path):
    tile, _ = written_tile(tmp_path, "rgbi")
    assert load(tile, estimate(4, 4, 4)).data.shape == (4, 4, 4)


def test_budget_below_the_estimate_is_refused(tmp_path, monkeypatch, no_band_reads):
    tile, _ = written_tile(tmp_path, "rgbi")
    forbid_allocation(monkeypatch)
    with pytest.raises(WorkingMemoryError) as info:
        load(tile, estimate(4, 4, 4) - 1)
    assert info.value.estimate_bytes == estimate(4, 4, 4)


# Memory discipline


def test_bands_are_read_into_one_buffer_with_gdal_cache_capped(
    tmp_path, monkeypatch, opened
):
    tile, data = written_tile(tmp_path, "rgbi")
    real_read = rasterio.io.DatasetReader.read
    reads = []

    def recording_read(self, *args, **kwargs):
        reads.append((args, kwargs, get_gdal_config("GDAL_CACHEMAX")))
        return real_read(self, *args, **kwargs)

    monkeypatch.setattr(rasterio.io.DatasetReader, "read", recording_read)
    cache_before = get_gdal_config("GDAL_CACHEMAX")
    raster = load(tile)
    assert [args for args, _, _ in reads] == [(1,), (2,), (3,), (4,)]
    buffers = [kwargs["out"] for _, kwargs, _ in reads]
    assert all(buffer is buffers[0] for buffer in buffers)
    assert buffers[0].shape == (4, 4) and buffers[0].dtype == np.uint8
    assert [kwargs.keys() for _, kwargs, _ in reads] == [{"out"}] * 4
    assert [cache for _, _, cache in reads] == [READ_OVERHEAD_BYTES] * 4
    assert opened[0]["cache"] == READ_OVERHEAD_BYTES
    assert get_gdal_config("GDAL_CACHEMAX") == cache_before
    np.testing.assert_array_equal(raster.data, data.transpose(1, 2, 0))
