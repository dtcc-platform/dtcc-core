import gc
import math
import os
import struct
import sys
import warnings
import tracemalloc
import weakref
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest
import rasterio
import requests
from affine import Affine
from rasterio.enums import ColorInterp, MaskFlags, Resampling
from rasterio.env import get_gdal_config
from rasterio.errors import NotGeoreferencedWarning
from rasterio.transform import from_bounds

from dtcc_core.io import raster_tiles
from dtcc_core.io.raster_tiles import (
    DEFAULT_MAX_MEMORY_BYTES,
    LOAD_LAYOUTS,
    MOSAIC_LAYOUTS,
    READ_OVERHEAD_BYTES,
    STRIP_BYTES,
    RasterTileLayoutError,
    RasterTileReadError,
    WorkingMemoryError,
    build_mosaic,
    load_tile,
    mosaic_grid,
    mosaic_minimum_bytes,
    without_tracebacks,
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
# chmod(0) cannot deny reads on Windows, and root ignores file permissions.
CANNOT_DENY_READ = sys.platform == "win32" or os.geteuid() == 0


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
    overviews=None,
    full_resolution=None,
    **profile,
):
    """Write a GeoTIFF labelled like the LM originals: R, G, B, then undefined.

    ``mask`` adds an internal per-dataset mask; ``band_masks`` adds a ``.msk``
    sidecar holding one mask per band, which GDAL reads with the file.
    ``overviews`` builds internal NEAREST overviews from ``data``; then
    ``full_resolution``, if given, replaces the full-resolution samples only, so
    overview samples differ from them.
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
            if overviews:
                dst.build_overviews(overviews, Resampling.nearest)
    if full_resolution is not None:
        with RASTERIO_OPEN(path, "r+") as dst:
            dst.write(full_resolution)
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


def make_tile(
    path, spektraltyp="rgb", extent=CELL, crs="EPSG:3006", size_bytes=None, **fields
):
    record = dict(
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
    record.update(fields)
    return RasterTile(**record)


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
            "readdir": get_gdal_config("GDAL_DISABLE_READDIR_ON_OPEN"),
        }
        calls.append(call)
        dataset = RASTERIO_OPEN(*args, **kwargs)
        call["driver"] = dataset.driver
        call["files"] = dataset.files
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
        "labelled alpha",
    ),
    "internal mask": (dict(spektraltyp="rgb", mask=True), "mask"),
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
    override_header(monkeypatch, nodatavals=nodatavals)
    with pytest.raises(RasterTileLayoutError) as info:
        load(tile)
    assert reason in info.value.reason


def override_header(monkeypatch, **attributes):
    """Make every Rasterio open return the real dataset with some attributes set."""

    class Header:
        def __init__(self, dataset):
            self._dataset = dataset

        def __getattr__(self, name):
            return getattr(self._dataset, name)

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            self._dataset.close()

    for name, value in attributes.items():
        setattr(Header, name, value)
    monkeypatch.setattr(
        rasterio, "open", lambda *args, **kwargs: Header(RASTERIO_OPEN(*args, **kwargs))
    )


@pytest.mark.parametrize(
    "flags",
    [
        ([], [], []),
        ([MaskFlags.all_valid], [MaskFlags.all_valid], []),
        ([MaskFlags.all_valid, MaskFlags.nodata],) * 3,
    ],
)
def test_only_nodata_or_no_masks_are_admitted(
    tmp_path, monkeypatch, no_band_reads, flags
):
    # GeoTIFF has no internal per-band masks (GDAL reports those with no
    # flags) and sidecars are ignored, so the flags are set on a real header.
    tile, _ = written_tile(tmp_path, "rgb")
    override_header(monkeypatch, mask_flag_enums=flags)
    with pytest.raises(RasterTileLayoutError, match="mask band"):
        load(tile)


def write_aux_nodata(path, value=7):
    Path(f"{path}.aux.xml").write_text(
        '<PAMDataset><PAMRasterBand band="1">'
        f"<NoDataValue>{value}</NoDataValue></PAMRasterBand></PAMDataset>"
    )


@pytest.mark.parametrize("sidecar", ["msk", "aux.xml"])
def test_sidecar_files_are_ignored(tmp_path, opened, sidecar):
    tile, data = written_tile(tmp_path, "rgb", band_masks=sidecar == "msk")
    if sidecar == "aux.xml":
        write_aux_nodata(tile.path)
    with RASTERIO_OPEN(tile.path) as src:
        # GDAL reads the sidecar unless told to ignore the directory.
        assert f"{tile.path}.{sidecar}" in src.files
    readdir_before = get_gdal_config("GDAL_DISABLE_READDIR_ON_OPEN")
    raster = load(tile)
    np.testing.assert_array_equal(raster.data, data.transpose(1, 2, 0))
    assert math.isnan(raster.nodata)
    assert [call["files"] for call in opened] == [[str(tile.path)]]
    assert opened[0]["readdir"] == "EMPTY_DIR"
    assert get_gdal_config("GDAL_DISABLE_READDIR_ON_OPEN") == readdir_before


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


@pytest.mark.skipif(CANNOT_DENY_READ, reason="file permissions cannot deny reads")
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


def damage_block(path, which=0):
    data = bytearray(path.read_bytes())
    assert data[:4] == b"II*\0"
    ifd = struct.unpack_from("<I", data, 4)[0]
    for entry in range(struct.unpack_from("<H", data, ifd)[0]):
        entry_offset = ifd + 2 + 12 * entry
        tag, kind, count, value = struct.unpack_from("<HHII", data, entry_offset)
        if tag == 324:  # TileOffsets
            assert kind == 4 and count > 1
            offsets = struct.unpack_from(f"<{count}I", data, value)
    offset = offsets[which]
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
    damage_block(path)
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


# Mosaic

ORIGIN = (672500.0, 6577508.0)  # top-left corner of the mosaic fixtures


def grid_tiff(
    path, width, height, pixel, origin=ORIGIN, *, bands=3, tag=1, data=None, **kwargs
):
    """GeoTIFF whose bands 1 and 2 hold each pixel's column and row, band 3 a tag."""
    if data is None:
        cols, rows = np.meshgrid(np.arange(width), np.arange(height))
        layers = [cols % 256, rows % 256, np.full_like(cols, tag)]
        layers += [np.full_like(cols, 77)] * (bands - 3)
        data = np.stack(layers).astype(np.uint8)
    x0, y0 = origin
    transform = kwargs.pop("transform", Affine(pixel, 0.0, x0, 0.0, -pixel, y0))
    return write_tiff(path, data, transform=transform, **kwargs)


def tile_at(path, spektraltyp="rgb", year=2020, **fields):
    """A tile record for a written file, with the file's own extent."""
    with RASTERIO_OPEN(path) as src:
        extent = tuple(src.bounds)
    fields.setdefault("id", path.stem)
    return make_tile(
        path,
        spektraltyp,
        extent=extent,
        datetime=datetime(year, 6, 1, tzinfo=timezone.utc),
        **fields,
    )


def level_arrays(path):
    """Samples and transform of the full resolution, then each overview level."""
    with RASTERIO_OPEN(path) as src:
        levels = [(src.read(), src.transform)]
        count = len(src.overviews(1))
    for level in range(count):
        with RASTERIO_OPEN(path, overview_level=level) as src:
            levels.append((src.read(), src.transform))
    return levels


def reference(tiles, bounds, resolution, level=None):
    """The mosaic computed per output pixel with exact rational arithmetic.

    ``level`` is the overview level each tile is read from (None: full
    resolution); the level rule itself is tested separately.
    """
    xmin, ymin, xmax, ymax = (Fraction(v) for v in bounds)
    r = Fraction(resolution)
    width = max(1, math.ceil((xmax - xmin) / r))
    height = max(1, math.ceil((ymax - ymin) / r))
    layers = []
    order = sorted(tiles, key=lambda t: (t.datetime, t.collection, t.id))
    for tile in reversed(order):
        levels = level_arrays(tile.path)
        base_samples, base = levels[0]
        samples, transform = levels[0 if level is None else level + 1]
        with RASTERIO_OPEN(tile.path) as src:
            nodata = src.nodatavals[0]
        layers.append((base, base_samples.shape[1:], samples, transform, nodata))
    out = np.zeros((height, width, 4), dtype=np.uint8)
    for i in range(height):
        y = ymax - (i + Fraction(1, 2)) * r
        for j in range(width):
            x = xmin + (j + Fraction(1, 2)) * r
            for base, (bh, bw), samples, t, nodata in layers:
                column = math.floor((x - Fraction(base.c)) / Fraction(base.a))
                row = math.floor((Fraction(base.f) - y) / Fraction(-base.e))
                if not (0 <= column < bw and 0 <= row < bh):
                    continue
                column = math.floor((x - Fraction(t.c)) / Fraction(t.a))
                row = math.floor((Fraction(t.f) - y) / Fraction(-t.e))
                column = min(max(column, 0), samples.shape[2] - 1)
                row = min(max(row, 0), samples.shape[1] - 1)
                sample = samples[:, row, column]
                if nodata is not None and (sample == nodata).all():
                    continue
                out[i, j, :3] = sample[:3]
                out[i, j, 3] = 255
                break
    return out


def mosaic(tiles, bounds, resolution=None, budget=DEFAULT_MAX_MEMORY_BYTES):
    return build_mosaic(
        tiles, bounds, resolution=resolution, max_memory_bytes=budget
    )


def extent_of(tile):
    e = tile.extent
    return (e.xmin, e.ymin, e.xmax, e.ymax)


def shifted(bounds, offsets, resolution):
    return tuple(v + d * resolution for v, d in zip(bounds, offsets))


def estimate_of(tiles, bounds, resolution=None):
    with pytest.raises(WorkingMemoryError) as info:
        mosaic(tiles, bounds, resolution, budget=1)
    return info.value.estimate_bytes


def block_bytes(rows, cols, bands, level_shape, ratio_x=1, ratio_y=1):
    """Documented bound on one block's arena share."""
    level_height, level_width = level_shape
    nh = min(math.ceil(rows * Fraction(ratio_y)) + 2, level_height)
    nw = min(math.ceil(cols * Fraction(ratio_x)) + 2, level_width)
    return (
        bands * nh * nw
        + bands * rows * nw
        + bands * rows * cols
        + 2 * rows * cols
        + 8 * (rows + cols)
    )


def arena_root(array):
    while array.base is not None:
        array = array.base
    return array


@pytest.fixture
def reads(monkeypatch):
    """Record the window, output buffer and GDAL settings of every pixel read."""
    calls = []
    real_read = rasterio.io.DatasetReader.read

    def recording_read(self, *args, **kwargs):
        calls.append(
            {
                "args": args,
                "kwargs": dict(kwargs),
                "cache": get_gdal_config("GDAL_CACHEMAX"),
                "path": self.name,
            }
        )
        return real_read(self, *args, **kwargs)

    monkeypatch.setattr(rasterio.io.DatasetReader, "read", recording_read)
    return calls


# Grid


def test_mosaic_layouts():
    assert MOSAIC_LAYOUTS == {"rgb": 3, "rgbi": 4}
    assert STRIP_BYTES == 32 * 1024**2


def test_grid_is_anchored_at_the_top_left_and_rounded_up():
    width, height, transform = mosaic_grid((100.0, 200.0, 1100.0, 1200.0), 0.16)
    assert (width, height) == (6250, 6250)
    assert transform == Affine(0.16, 0.0, 100.0, 0.0, -0.16, 1200.0)
    assert mosaic_grid((0.0, 0.0, 1000.001, 1.0), 0.16)[:2] == (6251, 7)


def test_grid_covers_a_small_excess_over_whole_pixels():
    assert mosaic_grid((0.0, 0.0, 1.0000005, 1.0), 1.0)[:2] == (2, 1)
    assert mosaic_grid((0.0, 0.0, 1.0, 1.0000005), 1)[:2] == (1, 2)


def test_grid_is_at_least_one_pixel():
    assert mosaic_grid((10.0, 10.0, 10.01, 10.02), 1.0)[:2] == (1, 1)


def test_minimum_bytes_counts_the_rgba_output_and_gdal_allowance():
    assert mosaic_minimum_bytes((0.0, 0.0, 10.0, 7.0), 0.5) == (
        20 * 14 * 4 + READ_OVERHEAD_BYTES
    )


# Arguments


BAD_BOUNDS = {
    "NaN": ((math.nan, 0.0, 1.0, 1.0), ValueError),
    "inf": ((0.0, 0.0, math.inf, 1.0), ValueError),
    "inverted": ((1.0, 0.0, 0.0, 1.0), ValueError),
    "zero span": ((0.0, 1.0, 1.0, 1.0), ValueError),
    "three values": ((0.0, 0.0, 1.0), ValueError),
    "six values": ((0.0, 0.0, 0.0, 1.0, 1.0, 1.0), ValueError),
    "text": (("0", 0.0, 1.0, 1.0), TypeError),
}


@pytest.mark.parametrize("case", BAD_BOUNDS, ids=list(BAD_BOUNDS))
def test_invalid_bounds_are_rejected_before_file_access(tmp_path, opened, case):
    bounds, error = BAD_BOUNDS[case]
    with pytest.raises(error):
        mosaic([make_tile(tmp_path / "missing.tif")], bounds)
    assert opened == []


@pytest.mark.parametrize(
    "resolution, error",
    [(True, TypeError), ("1", TypeError), (0, ValueError), (-1.0, ValueError),
     (math.nan, ValueError), (math.inf, ValueError)],
)
def test_invalid_resolution_is_rejected_before_file_access(
    tmp_path, opened, resolution, error
):
    with pytest.raises(error):
        mosaic([make_tile(tmp_path / "missing.tif")], CELL, resolution)
    assert opened == []


@pytest.mark.parametrize(
    "budget, error",
    [(True, TypeError), (1.5, TypeError), ("2", TypeError), (0, ValueError),
     (-1, ValueError)],
)
def test_invalid_mosaic_budget_is_rejected_before_file_access(
    tmp_path, opened, budget, error
):
    with pytest.raises(error):
        mosaic([make_tile(tmp_path / "missing.tif")], CELL, 1.0, budget)
    assert opened == []


# Sampling

EDGE_OFFSETS = {
    "aligned": (0.0, 0.0, 0.0, 0.0),
    "outside every edge": (-0.3, -0.7, 0.45, 0.2),
    "inside every edge": (0.3, 0.7, -0.45, -0.2),
    "half pixel inside": (0.5, 0.5, -0.5, -0.5),
    "top-left corner outside": (-0.4, 0.4, -0.3, 0.3),
    "bottom-right corner outside": (0.3, -0.4, 0.4, -0.3),
}


@pytest.mark.parametrize("resolution", [0.16, 0.4, 0.5, 1.0])
@pytest.mark.parametrize("case", EDGE_OFFSETS, ids=list(EDGE_OFFSETS))
def test_samples_follow_exact_pixel_centres_at_tile_edges(tmp_path, resolution, case):
    tile = tile_at(grid_tiff(tmp_path / "t.tif", 20, 20, 0.4))
    bounds = shifted(extent_of(tile), EDGE_OFFSETS[case], resolution)
    result = mosaic([tile], bounds, resolution)
    np.testing.assert_array_equal(
        result.raster.data, reference([tile], bounds, resolution)
    )
    exchange.validate(result.raster)


def test_seams_and_centres_on_edges_go_right_and_down(tmp_path):
    x0, y0 = ORIGIN
    tiles = []
    for tag, (dx, dy) in enumerate([(0, 0), (5, 0), (0, -5), (5, -5)], start=1):
        path = grid_tiff(
            tmp_path / f"t{tag}.tif", 10, 10, 0.5, (x0 + dx, y0 + dy), tag=tag
        )
        tiles.append(tile_at(path))
    # Output centres fall on whole metres: on tile edges and seams. The last
    # column and row have centres on the right and bottom edges, outside.
    bounds = (x0 - 0.5, y0 - 10.5, x0 + 10.5, y0 + 0.5)
    data = mosaic(tiles, bounds, 1.0).raster.data
    np.testing.assert_array_equal(data, reference(tiles, bounds, 1.0))
    assert data.shape == (11, 11, 4)
    assert (data[10, :, 3] == 0).all() and (data[:, 10, 3] == 0).all()
    assert (data[:10, :10, 3] == 255).all()
    assert data[0, 0, 2] == 1 and tuple(data[0, 0, :2]) == (0, 0)
    assert data[4, 4, 2] == 1 and tuple(data[4, 4, :2]) == (8, 8)
    assert data[0, 5, 2] == 2 and data[0, 5, 0] == 0
    assert data[5, 0, 2] == 3 and data[5, 0, 1] == 0
    assert data[5, 5, 2] == 4


@pytest.mark.parametrize("manifest_resolution", [None, 0.4, 1.0])
def test_default_resolution_is_the_finest_header(tmp_path, manifest_resolution):
    fine = tile_at(
        grid_tiff(tmp_path / "fine.tif", 25, 25, 0.16, tag=1),
        resolution=manifest_resolution,
    )
    x0, y0 = ORIGIN
    coarse = tile_at(grid_tiff(tmp_path / "coarse.tif", 10, 10, 0.4, (x0 + 4, y0)))
    bounds = (x0, y0 - 4.0, x0 + 8.0, y0)
    result = mosaic([coarse, fine], bounds)
    assert result.resolution == 0.16
    assert result.raster.georef == Affine(0.16, 0.0, x0, 0.0, -0.16, y0)
    np.testing.assert_array_equal(
        result.raster.data, reference([coarse, fine], bounds, 0.16)
    )


def test_result_records_requested_and_actual_grid(tmp_path):
    tile = tile_at(grid_tiff(tmp_path / "t.tif", 20, 20, 0.4))
    x0, y0 = ORIGIN
    bounds = (x0 + 0.1, y0 - 7.9, x0 + 7.75, y0 - 0.05)
    result = mosaic([tile], bounds, 0.5)
    width, height, transform = mosaic_grid(bounds, 0.5)
    assert result.requested_bounds == bounds
    assert result.bounds == (
        bounds[0], bounds[3] - height * 0.5, bounds[0] + width * 0.5, bounds[3]
    )
    assert result.resolution == 0.5
    assert result.raster.data.shape == (height, width, 4)
    assert result.raster.georef == transform
    assert result.raster.crs == "EPSG:3006"
    assert math.isnan(result.raster.nodata)


# Levels


def overview_fixture(path, width, height, pixel, factors):
    """Overviews built from tag 10 samples; full resolution rewritten with tag 200."""
    cols, rows = np.meshgrid(np.arange(width), np.arange(height))
    first = np.stack([cols % 256, rows % 256, np.full_like(cols, 10)]).astype(np.uint8)
    final = first.copy()
    final[2] = 200
    return grid_tiff(
        path,
        width,
        height,
        pixel,
        data=first,
        tiled=True,
        blockxsize=16,
        blockysize=16,
        overviews=factors,
        full_resolution=final,
    )


def level_spacing(path, level):
    with RASTERIO_OPEN(path, overview_level=level) as src:
        return src.transform.a, -src.transform.e


@pytest.mark.parametrize(
    "resolution, level",
    [
        (0.2, None),
        (0.4, None),
        (0.79, None),
        ("level 0", 0),
        (1.0, 0),
        ("level 1", 1),
        (2.0, 1),
        ("level 2", 2),
        (10.0, 2),
    ],
)
def test_level_is_the_coarsest_not_coarser_than_the_output(
    tmp_path, resolution, level
):
    path = overview_fixture(tmp_path / "t.tif", 50, 50, 0.4, [2, 4, 8])
    if isinstance(resolution, str):
        resolution = level_spacing(path, int(resolution[-1]))[0]
    tile = tile_at(path)
    bounds = extent_of(tile)
    data = mosaic([tile], bounds, resolution).raster.data
    opaque = data[..., 3] == 255
    assert set(np.unique(data[..., 2][opaque])) == {200 if level is None else 10}
    np.testing.assert_array_equal(
        data, reference([tile], bounds, resolution, level=level)
    )


@pytest.mark.parametrize("resolution, level", [(3.6, None), (3.9, 0)])
def test_level_needs_both_spacings_within_the_output(tmp_path, resolution, level):
    path = overview_fixture(tmp_path / "t.tif", 7, 50, 1.0, [4])
    assert level_spacing(path, 0) == pytest.approx((3.5, 50 / 13))
    tile = tile_at(path)
    bounds = extent_of(tile)
    data = mosaic([tile], bounds, resolution).raster.data
    opaque = data[..., 3] == 255
    assert set(np.unique(data[..., 2][opaque])) == {200 if level is None else 10}
    np.testing.assert_array_equal(
        data, reference([tile], bounds, resolution, level=level)
    )


def test_level_indices_are_clipped_to_the_level(tmp_path):
    # A 99 px tile's 13 px overview ends a few 1e-15 m before the tile's edge:
    # an output centre in that gap is inside the tile but past the level.
    path = grid_tiff(tmp_path / "t.tif", 99, 99, 0.4, overviews=[2, 4, 8])
    tile = tile_at(path)
    with RASTERIO_OPEN(path) as src:
        tile_edge = Fraction(src.transform.c) + 99 * Fraction(src.transform.a)
    with RASTERIO_OPEN(path, overview_level=2) as src:
        assert src.width == 13
        level_edge = Fraction(src.transform.c) + 13 * Fraction(src.transform.a)
    assert level_edge < tile_edge
    xmin = float(math.floor(level_edge) - 2)
    resolution = float(2 * (level_edge - Fraction(xmin)))
    for _ in range(64):
        if level_edge <= Fraction(xmin) + Fraction(resolution) / 2 < tile_edge:
            break
        resolution = np.nextafter(resolution, math.inf)
    else:
        raise AssertionError("no float resolution puts a centre in the gap")
    resolution = float(resolution)
    y0 = ORIGIN[1]
    bounds = (xmin, y0 - resolution, xmin + resolution, y0)
    result = mosaic([tile], bounds, resolution)
    assert result.failures == ()
    rows, cols, _ = result.raster.data.shape
    assert cols == 1 and result.valid_pixels == rows
    np.testing.assert_array_equal(
        result.raster.data, reference([tile], bounds, resolution, level=2)
    )


def test_indices_are_exact_where_float_arithmetic_rounds_down(tmp_path):
    origin = (6577500.0, 6577500.0)
    tile = tile_at(grid_tiff(tmp_path / "t.tif", 10, 10, 1e-5, origin))
    bounds = extent_of(tile)
    # The first output centre lies exactly on the edge of source pixels 0 and 1;
    # float arithmetic gives 0.99996... there.
    assert math.floor((origin[0] + 1e-5 - origin[0]) / 1e-5) == 0
    data = mosaic([tile], bounds, 2e-5).raster.data
    assert tuple(data[0, 0, :2]) == (1, 1)
    np.testing.assert_array_equal(data, reference([tile], bounds, 2e-5))


# Blocks and memory


def test_blocks_split_in_both_directions_with_the_same_result(
    tmp_path, monkeypatch, reads
):
    tile = tile_at(grid_tiff(tmp_path / "t.tif", 50, 40, 0.4))
    bounds = extent_of(tile)
    whole = mosaic([tile], bounds, 0.4).raster.data
    reads.clear()
    monkeypatch.setattr(raster_tiles, "STRIP_BYTES", 600)
    split = mosaic([tile], bounds, 0.4).raster.data
    np.testing.assert_array_equal(split, whole)
    windows = [call["kwargs"]["window"] for call in reads]
    assert len({w.col_off for w in windows}) > 1
    assert len({w.row_off for w in windows}) > 1
    for call in reads:
        out = call["kwargs"]["out"]
        assert out.nbytes <= 600


def test_arena_is_the_largest_single_block(tmp_path):
    x0, y0 = ORIGIN
    square = tile_at(grid_tiff(tmp_path / "square.tif", 100, 100, 0.4, tag=1))
    strip = tile_at(
        grid_tiff(tmp_path / "strip.tif", 2000, 1, 0.4, (x0, y0 - 40.0), tag=2)
    )
    bounds = (x0, y0 - 40.4, x0 + 800.0, y0)
    width, height, _ = mosaic_grid(bounds, 0.4)
    square_block = block_bytes(100, 100, 3, (100, 100))
    strip_block = block_bytes(1, 2000, 3, (1, 2000))
    # The square tile needs the most pixel memory, the strip the most indices.
    assert 8 * 2001 > 8 * 200 and square_block > strip_block
    arena = estimate_of([square, strip], bounds, 0.4) - width * height * 4
    assert arena - READ_OVERHEAD_BYTES == max(square_block, strip_block)


def test_numpy_and_python_memory_stays_within_the_estimate(tmp_path):
    rng = np.random.default_rng(0)
    data = rng.integers(1, 255, (4, 800, 800), dtype=np.uint8)
    path = grid_tiff(tmp_path / "t.tif", 800, 800, 0.4, data=data, nodata=0)
    tile = tile_at(path, "rgbi")
    bounds = extent_of(tile)
    estimate = estimate_of([tile], bounds, 0.64)
    width, height, _ = mosaic_grid(bounds, 0.64)
    arena = estimate - width * height * 4 - READ_OVERHEAD_BYTES
    assert arena > 4 * 1024**2
    tracemalloc.start()
    try:
        tracemalloc.reset_peak()
        result = mosaic([tile], bounds, 0.64)
        peak = tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()
    assert result.valid_pixels == width * height
    assert peak <= width * height * 4 + arena + 256 * 1024


@pytest.mark.parametrize("case", ["mid-read", "header", "no valid pixel"])
def test_recorded_failures_keep_no_mosaic_buffer_alive(tmp_path, monkeypatch, case):
    older = tile_at(grid_tiff(tmp_path / "old.tif", 40, 40, 0.4, tag=1), year=2010)
    path = grid_tiff(
        tmp_path / "new.tif", 40, 40, 0.4, tag=2, tiled=True, blockxsize=16,
        blockysize=16, compress="deflate",
    )
    if case == "mid-read":
        damage_block(path, -1)
        tiles = [older, tile_at(path, year=2020)]
        monkeypatch.setattr(raster_tiles, "STRIP_BYTES", 5000)
    elif case == "header":
        # Refused when opened, before the buffers exist; the error has a cause.
        tiles = [older, bad_tile(tmp_path, "non-raster", older)]
    else:
        damage_block(path, 0)
        tiles = [tile_at(path, year=2020)]
    buffers = []
    real_zeros = np.zeros
    real_read = rasterio.io.DatasetReader.read

    def recording_zeros(*args, **kwargs):
        array = real_zeros(*args, **kwargs)
        buffers.append(weakref.ref(array))
        return array

    def recording_read(self, *args, **kwargs):
        buffers.append(weakref.ref(arena_root(kwargs["out"])))
        return real_read(self, *args, **kwargs)

    monkeypatch.setattr(np, "zeros", recording_zeros)
    monkeypatch.setattr(rasterio.io.DatasetReader, "read", recording_read)
    # With cyclic collection off, a buffer outlives the result only if the
    # recorded failures still reach it.
    gc.disable()
    try:
        result = mosaic(tiles, extent_of(older))
        failures = result.failures
        del result
        assert len(failures) == 1 and len(buffers) >= 2
        assert [buffer for buffer in buffers if buffer() is not None] == []
    finally:
        gc.enable()


def test_without_tracebacks_clears_every_chained_traceback():
    def caught(error, cause=None):
        try:
            raise error from cause
        except Exception as raised:
            return raised

    cause = caught(OSError("cause"))
    try:
        raise KeyError("context")
    except KeyError:
        error = caught(ValueError("error"), cause)
    context = error.__context__
    assert isinstance(context, KeyError) and error.__cause__ is cause
    cause.__context__ = error  # a cycle
    chain = [error, cause, context]
    assert all(link.__traceback__ is not None for link in chain)
    assert without_tracebacks(error) is error
    assert [link.__traceback__ for link in chain] == [None, None, None]


# Validity and alpha


def test_validity_uses_every_band_and_nir_never_drives_opacity(tmp_path):
    data = np.array(
        [[(0, 0, 0, 50), (10, 20, 30, 0)], [(0, 0, 0, 0), (1, 2, 3, 4)]],
        dtype=np.uint8,
    ).transpose(2, 0, 1)
    # 0.5 m pixels keep the file's float extent an exact number of pixels.
    tile = tile_at(grid_tiff(tmp_path / "t.tif", 2, 2, 0.5, data=data, nodata=0),
                   "rgbi")
    result = mosaic([tile], extent_of(tile))
    expected = [
        [(0, 0, 0, 255), (10, 20, 30, 255)],
        [(0, 0, 0, 0), (1, 2, 3, 255)],
    ]
    np.testing.assert_array_equal(result.raster.data, np.array(expected))
    assert result.valid_pixels == 3


@pytest.mark.parametrize("nodata, alpha", [(0, 0), (None, 255)])
def test_black_rgb_is_missing_only_with_declared_nodata(tmp_path, nodata, alpha):
    data = np.array([[(0, 0, 0), (1, 2, 3)]], dtype=np.uint8).transpose(2, 0, 1)
    tile = tile_at(grid_tiff(tmp_path / "t.tif", 2, 1, 0.5, data=data, nodata=nodata))
    result = mosaic([tile], extent_of(tile))
    np.testing.assert_array_equal(
        result.raster.data, np.array([[(0, 0, 0, alpha), (1, 2, 3, 255)]])
    )


# Precedence


def same_cell(tmp_path, specs):
    """Tiles over one cell: (name, year, collection, tag[, data])."""
    tiles = []
    for name, year, collection, tag, *data in specs:
        path = grid_tiff(
            tmp_path / f"{name}.tif", 10, 10, 0.4, tag=tag, nodata=0,
            data=data[0] if data else None,
        )
        tiles.append(tile_at(path, year=year, collection=collection, id=name))
    return tiles


def test_newest_tile_wins_regardless_of_manifest_order(tmp_path):
    older, newer = same_cell(
        tmp_path, [("old", 2010, "orto-a", 1), ("new", 2020, "orto-a", 2)]
    )
    result = mosaic([older, newer], extent_of(older))
    assert set(np.unique(result.raster.data[..., 2])) == {2}
    assert result.sources == ((newer, 100),)
    assert result.valid_pixels == 100


@pytest.mark.parametrize(
    "first, second",
    [
        (("y", 2020, "orto-a", 1), ("x", 2020, "orto-b", 2)),
        (("a1", 2020, "orto-a", 1), ("a2", 2020, "orto-a", 2)),
    ],
)
def test_ties_fall_to_collection_then_id(tmp_path, first, second):
    tiles = same_cell(tmp_path, [second, first])
    result = mosaic(tiles, extent_of(tiles[0]))
    assert set(np.unique(result.raster.data[..., 2])) == {2}


def test_older_tile_fills_where_the_newer_is_empty(tmp_path):
    cols, rows = np.meshgrid(np.arange(10), np.arange(10))
    partial = np.stack([cols, rows, np.full_like(cols, 2)]).astype(np.uint8)
    partial[:, :, :5] = 0
    older, newer = same_cell(
        tmp_path,
        [("old", 2010, "orto-a", 1), ("new", 2020, "orto-a", 2, partial)],
    )
    result = mosaic([newer, older], extent_of(older))
    tags = result.raster.data[..., 2]
    assert (tags[:, :5] == 1).all() and (tags[:, 5:] == 2).all()
    assert result.sources == ((newer, 50), (older, 50))
    np.testing.assert_array_equal(
        result.raster.data, reference([newer, older], extent_of(older), 0.4)
    )


# Coverage


def test_area_outside_every_tile_is_transparent(tmp_path):
    tile = tile_at(grid_tiff(tmp_path / "t.tif", 10, 10, 0.4))
    x0, y0 = ORIGIN
    bounds = (x0, y0 - 4.0, x0 + 6.0, y0)
    result = mosaic([tile], bounds)
    data = result.raster.data
    assert data.shape == (10, 15, 4)
    assert (data[:, 10:] == 0).all() and (data[:, :10, 3] == 255).all()
    assert result.valid_pixels == 100


def test_tiles_outside_the_bounds_are_not_opened(tmp_path, opened):
    x0, y0 = ORIGIN
    inside = tile_at(grid_tiff(tmp_path / "in.tif", 10, 10, 0.4))
    touching = tile_at(grid_tiff(tmp_path / "touch.tif", 10, 10, 0.4, (x0 + 4, y0)))
    away = tile_at(grid_tiff(tmp_path / "away.tif", 10, 10, 0.4, (x0 + 100, y0)))
    mosaic([inside, touching, away], extent_of(inside))
    paths = {Path(call["args"][0]) for call in opened}
    assert paths == {inside.path}


def assert_empty(result):
    assert result.raster.data.shape == ()
    assert math.isnan(result.raster.data)
    assert result.raster.crs == "EPSG:3006"
    assert (result.bounds, result.resolution) == (None, None)
    assert result.sources == ()
    assert result.valid_pixels == 0


def test_no_tiles_give_the_empty_result():
    result = mosaic([], CELL, 1.0)
    assert_empty(result)
    assert result.failures == ()
    assert result.requested_bounds == CELL


def test_only_failing_headers_give_the_empty_result(tmp_path):
    tile = tile_at(grid_tiff(tmp_path / "t.tif", 10, 10, 0.4), "cir")
    result = mosaic([tile], extent_of(tile), 0.4)
    assert_empty(result)
    assert [(f.tile, f.pixels) for f in result.failures] == [(tile, 0)]


def test_failed_reads_give_the_empty_result(tmp_path):
    path = grid_tiff(
        tmp_path / "t.tif", 32, 32, 0.4, tiled=True, blockxsize=16, blockysize=16,
        compress="deflate",
    )
    damage_block(path)
    tile = tile_at(path)
    result = mosaic([tile], extent_of(tile))
    assert_empty(result)
    (failure,) = result.failures
    assert failure.tile == tile and failure.pixels == 0
    assert isinstance(failure.error, RasterTileReadError)
    assert path.is_file()


def test_all_nodata_samples_give_the_empty_result(tmp_path):
    data = np.zeros((3, 10, 10), dtype=np.uint8)
    tile = tile_at(grid_tiff(tmp_path / "t.tif", 10, 10, 0.4, data=data, nodata=0))
    result = mosaic([tile], extent_of(tile))
    assert_empty(result)
    assert result.failures == ()


# Failures


def bad_tile(tmp_path, kind, good):
    path = tmp_path / f"bad-{kind.replace(' ', '-')}.tif"
    x0, y0 = ORIGIN
    if kind == "cir":
        return tile_at(grid_tiff(path, 10, 10, 0.4), "cir")
    if kind == "pan":
        data = np.ones((1, 10, 10), dtype=np.uint8)
        path = grid_tiff(path, 10, 10, 0.4, data=data, colorinterp=[ColorInterp.gray])
        return tile_at(path, "pan")
    if kind == "non-raster":
        path.write_text("<html><body>Bad gateway</body></html>")
        return make_tile(path, extent=extent_of(good), id="bad")
    if kind == "VRT":
        path.write_text(
            '<VRTDataset rasterXSize="10" rasterYSize="10"><SRS>EPSG:3006</SRS>'
            f"<GeoTransform>{x0},0.4,0,{y0},0,-0.4</GeoTransform>"
            + "".join(
                f'<VRTRasterBand dataType="Byte" band="{b}"><SimpleSource>'
                f'<SourceFilename relativeToVRT="0">{good.path}</SourceFilename>'
                f"<SourceBand>{b}</SourceBand></SimpleSource></VRTRasterBand>"
                for b in (1, 2, 3)
            )
            + "</VRTDataset>"
        )
        return make_tile(path, extent=extent_of(good), id="bad")
    if kind == "rotated":
        transform = Affine(0.4, 1e-6, x0, 1e-6, -0.4, y0)
        return tile_at(grid_tiff(path, 10, 10, 0.4, transform=transform))
    if kind == "non-square pixels":
        transform = Affine(0.4, 0.0, x0, 0.0, -0.5, y0)
        return tile_at(grid_tiff(path, 10, 8, 0.4, transform=transform))
    if kind == "uint16":
        data = np.ones((3, 10, 10), dtype=np.uint16)
        return tile_at(grid_tiff(path, 10, 10, 0.4, data=data))
    if kind == "alpha band":
        colors = LM_COLORS + [ColorInterp.alpha]
        return tile_at(grid_tiff(path, 10, 10, 0.4, bands=4, colorinterp=colors),
                       "rgbi")
    if kind == "internal mask":
        return tile_at(grid_tiff(path, 10, 10, 0.4, mask=True))
    if kind == "other CRS":
        return tile_at(grid_tiff(path, 10, 10, 0.4, crs="EPSG:3021"))
    if kind == "wrong cell":
        tile = tile_at(grid_tiff(path, 10, 10, 0.4))
        return make_tile(path, extent=extent_of(good), id="bad",
                         datetime=tile.datetime)
    raise AssertionError(kind)


BAD_KINDS = [
    "cir", "pan", "non-raster", "VRT", "rotated", "non-square pixels", "uint16",
    "alpha band", "internal mask", "other CRS", "wrong cell",
]


@pytest.mark.parametrize("kind", BAD_KINDS)
def test_source_failures_are_recorded_and_the_rest_still_mosaics(tmp_path, kind):
    good = tile_at(grid_tiff(tmp_path / "good.tif", 10, 10, 0.4, tag=9), year=2010)
    bad = bad_tile(tmp_path, kind, good)
    if kind == "wrong cell":
        # The file's own extent is elsewhere; the record claims the good cell.
        x0, y0 = ORIGIN
        grid_tiff(bad.path, 10, 10, 0.4, (x0 + 1.0, y0))
    bounds = extent_of(good)
    result = mosaic([bad, good], bounds, 0.4)
    (failure,) = result.failures
    assert failure.tile == bad and failure.pixels == 0
    assert isinstance(failure.error, (RasterTileLayoutError, RasterTileReadError))
    assert str(bad.path) in str(failure.error)
    np.testing.assert_array_equal(
        result.raster.data, reference([good], bounds, 0.4)
    )
    assert result.sources == ((good, 100),)


@pytest.mark.parametrize(
    "kind, reason", [("rotated", "north-up"), ("non-square pixels", "square")]
)
def test_transform_refusals_name_the_reason(tmp_path, kind, reason):
    good = tile_at(grid_tiff(tmp_path / "good.tif", 10, 10, 0.4, tag=9))
    bad = bad_tile(tmp_path, kind, good)
    result = mosaic([bad], extent_of(bad), 0.4)
    assert reason in result.failures[0].error.reason


def test_failure_mid_read_keeps_composited_blocks(tmp_path, monkeypatch):
    older = tile_at(grid_tiff(tmp_path / "old.tif", 40, 40, 0.4, tag=1), year=2010)
    path = grid_tiff(
        tmp_path / "new.tif", 40, 40, 0.4, tag=2, tiled=True, blockxsize=16,
        blockysize=16, compress="deflate",
    )
    damage_block(path, -1)  # the bottom-right 16 x 16 block
    newer = tile_at(path, year=2020)
    monkeypatch.setattr(raster_tiles, "STRIP_BYTES", 5000)
    result = mosaic([older, newer], extent_of(older), 0.4)
    tags = result.raster.data[..., 2]
    (failure,) = result.failures
    assert failure.tile == newer
    assert isinstance(failure.error, RasterTileReadError)
    assert 0 < failure.pixels == np.count_nonzero(tags == 2)
    assert (tags[32:, 32:] == 1).all()
    assert (result.raster.data[..., 3] == 255).all()
    assert result.sources == ((newer, failure.pixels), (older, 1600 - failure.pixels))
    assert path.is_file()


@pytest.mark.parametrize(
    "replacement, reason",
    [
        ({"bands": 4}, "4 bands, expected 3"),
        # Same geometry and bands: only readmission sees the change.
        ({"crs": "EPSG:3021"}, "CRS is EPSG:3021"),
    ],
)
def test_layout_change_before_reading_is_a_failure(
    tmp_path, monkeypatch, replacement, reason
):
    older = tile_at(grid_tiff(tmp_path / "old.tif", 10, 10, 0.4, tag=1), year=2010)
    newer = tile_at(grid_tiff(tmp_path / "new.tif", 10, 10, 0.4, tag=2), year=2020)
    opens = []

    def replacing_open(path, *args, **kwargs):
        if Path(path) == newer.path and "overview_level" not in kwargs:
            opens.append(path)
            if len(opens) == 2:
                grid_tiff(newer.path, 10, 10, 0.4, tag=2, **replacement)
        return RASTERIO_OPEN(path, *args, **kwargs)

    monkeypatch.setattr(rasterio, "open", replacing_open)
    result = mosaic([older, newer], extent_of(older), 0.4)
    (failure,) = result.failures
    assert failure.tile == newer and failure.pixels == 0
    assert isinstance(failure.error, RasterTileLayoutError)
    assert reason in failure.error.reason
    assert set(np.unique(result.raster.data[..., 2])) == {1}


def test_same_bounds_replacement_before_reading_is_a_failure(tmp_path, monkeypatch):
    older = tile_at(grid_tiff(tmp_path / "old.tif", 10, 10, 0.4, tag=1), year=2010)
    newer = tile_at(grid_tiff(tmp_path / "new.tif", 10, 10, 0.4, tag=2), year=2020)
    opens = []

    def replacing_open(path, *args, **kwargs):
        if Path(path) == newer.path and "overview_level" not in kwargs:
            opens.append(path)
            if len(opens) == 2:
                # Admissible and with the same bounds, at a different resolution.
                grid_tiff(newer.path, 20, 20, 0.2, tag=2)
        return RASTERIO_OPEN(path, *args, **kwargs)

    monkeypatch.setattr(rasterio, "open", replacing_open)
    result = mosaic([older, newer], extent_of(older), 0.4)
    (failure,) = result.failures
    assert failure.tile == newer and failure.pixels == 0
    assert isinstance(failure.error, RasterTileLayoutError)
    assert "changed" in failure.error.reason
    assert set(np.unique(result.raster.data[..., 2])) == {1}


def test_replaced_overview_before_reading_is_a_failure(tmp_path, monkeypatch):
    older = tile_at(grid_tiff(tmp_path / "old.tif", 50, 50, 0.4, tag=1), year=2010)
    path = grid_tiff(tmp_path / "new.tif", 50, 50, 0.4, tag=2, overviews=[2, 4, 8])
    newer = tile_at(path, year=2020)
    # The same full-resolution geometry, but level 0 is 17 px, not 25.
    other = grid_tiff(tmp_path / "other.tif", 50, 50, 0.4, tag=2, overviews=[3])
    level_opens = []

    def replacing_open(path, *args, **kwargs):
        dataset = RASTERIO_OPEN(path, *args, **kwargs)
        if Path(path) == newer.path and "overview_level" in kwargs:
            level_opens.append(path)
            if len(level_opens) == 3:
                # The header pass has read every level: replace the file
                # before compositing takes its identity.
                os.replace(other, newer.path)
        return dataset

    monkeypatch.setattr(rasterio, "open", replacing_open)
    result = mosaic([older, newer], extent_of(older), 1.0)
    assert len(level_opens) == 4
    (failure,) = result.failures
    assert failure.tile == newer and failure.pixels == 0
    assert "changed" in failure.error.reason
    opaque = result.raster.data[..., 3] == 255
    assert set(np.unique(result.raster.data[..., 2][opaque])) == {1}


def test_overview_from_a_replaced_file_is_admitted_before_reading(
    tmp_path, monkeypatch
):
    older = tile_at(grid_tiff(tmp_path / "old.tif", 50, 50, 0.4, tag=1), year=2010)
    path = grid_tiff(tmp_path / "new.tif", 50, 50, 0.4, tag=2, overviews=[2, 4, 8])
    newer = tile_at(path, year=2020)
    level_opens = []

    def replacing_open(path, *args, **kwargs):
        if Path(path) == newer.path and "overview_level" in kwargs:
            level_opens.append(path)
            if len(level_opens) == 4:
                # Same geometry, bands and overview sizes; another CRS.
                grid_tiff(
                    newer.path, 50, 50, 0.4, tag=2, overviews=[2, 4, 8],
                    crs="EPSG:3021",
                )
        return RASTERIO_OPEN(path, *args, **kwargs)

    monkeypatch.setattr(rasterio, "open", replacing_open)
    result = mosaic([older, newer], extent_of(older), 1.0)
    assert len(level_opens) == 4
    (failure,) = result.failures
    assert failure.tile == newer and failure.pixels == 0
    assert "CRS is EPSG:3021" in failure.error.reason
    opaque = result.raster.data[..., 3] == 255
    assert set(np.unique(result.raster.data[..., 2][opaque])) == {1}


def non_square_with_the_same_overview(path):
    """49 x 50 pixels over the same 20 m cell: non-square full resolution whose
    factor 2 overview is 25 x 25 at 0.8 m, like a 50 x 50 0.4 m file's."""
    x0, y0 = ORIGIN
    transform = Affine(20 / 49, 0.0, x0, 0.0, -0.4, y0)
    grid_tiff(path, 49, 50, 0.4, tag=3, overviews=[2], transform=transform)


@pytest.mark.parametrize("moment", ["before the full resolution", "before the level"])
def test_overview_plans_read_one_admitted_file(tmp_path, monkeypatch, moment):
    older = tile_at(grid_tiff(tmp_path / "old.tif", 50, 50, 0.4, tag=1), year=2010)
    path = grid_tiff(tmp_path / "new.tif", 50, 50, 0.4, tag=2, overviews=[2])
    newer = tile_at(path, year=2020)
    with RASTERIO_OPEN(path, overview_level=0) as src:
        level_geometry = (src.transform, src.width, src.height)
    base_opens, level_opens = [], []

    def replacing_open(path, *args, **kwargs):
        if Path(path) == newer.path:
            opens = level_opens if "overview_level" in kwargs else base_opens
            opens.append(path)
            if (moment == "before the full resolution" and opens is base_opens
                    and len(base_opens) == 2) or (
                moment == "before the level" and opens is level_opens
                and len(level_opens) == 2
            ):
                non_square_with_the_same_overview(newer.path)
        return RASTERIO_OPEN(path, *args, **kwargs)

    monkeypatch.setattr(rasterio, "open", replacing_open)
    result = mosaic([older, newer], extent_of(older), 0.8)
    # Refused at the full resolution, the level is not opened again.
    assert len(base_opens) == 2
    assert len(level_opens) == (1 if moment == "before the full resolution" else 2)
    with RASTERIO_OPEN(newer.path, overview_level=0) as src:
        assert (src.transform, src.width, src.height) == level_geometry
    (failure,) = result.failures
    assert failure.tile == newer and failure.pixels == 0
    expected = "square" if moment == "before the full resolution" else "changed"
    assert expected in failure.error.reason
    opaque = result.raster.data[..., 3] == 255
    assert set(np.unique(result.raster.data[..., 2][opaque])) == {1}


@pytest.mark.parametrize("change", ["replaced, same size and mtime", "rewritten"])
def test_file_identity_covers_replacement_and_rewriting(tmp_path, monkeypatch, change):
    older = tile_at(grid_tiff(tmp_path / "old.tif", 50, 50, 0.4, tag=1), year=2010)
    path = grid_tiff(tmp_path / "new.tif", 50, 50, 0.4, tag=2, overviews=[2])
    newer = tile_at(path, year=2020)
    # An admissible file of the same size and geometry with other samples.
    other = grid_tiff(tmp_path / "other.tif", 50, 50, 0.4, tag=3, overviews=[2])
    assert other.stat().st_size == path.stat().st_size
    level_opens = []

    def changing_open(path, *args, **kwargs):
        if Path(path) == newer.path and "overview_level" in kwargs:
            level_opens.append(path)
            if len(level_opens) == 2:
                original = newer.path.stat()
                if change == "rewritten":
                    # Same inode and size; only the modification time changes.
                    with open(newer.path, "r+b") as file:
                        file.write(other.read_bytes())
                    assert newer.path.stat().st_ino == original.st_ino
                    assert newer.path.stat().st_mtime_ns != original.st_mtime_ns
                else:
                    os.utime(other, ns=(original.st_atime_ns, original.st_mtime_ns))
                    os.replace(other, newer.path)
                    assert newer.path.stat().st_ino != original.st_ino
                    assert newer.path.stat().st_mtime_ns == original.st_mtime_ns
                assert newer.path.stat().st_size == original.st_size
        return RASTERIO_OPEN(path, *args, **kwargs)

    monkeypatch.setattr(rasterio, "open", changing_open)
    result = mosaic([older, newer], extent_of(older), 0.8)
    (failure,) = result.failures
    assert failure.tile == newer and failure.pixels == 0
    assert "changed" in failure.error.reason
    opaque = result.raster.data[..., 3] == 255
    assert set(np.unique(result.raster.data[..., 2][opaque])) == {1}


def test_overview_of_a_file_with_declared_nodata_keeps_its_gaps(tmp_path):
    cols, rows = np.meshgrid(np.arange(40), np.arange(40))
    data = np.stack([cols, rows, np.full_like(cols, 5)]).astype(np.uint8)
    data[:, :, :20] = 0
    path = grid_tiff(tmp_path / "t.tif", 40, 40, 0.4, data=data, nodata=0,
                     overviews=[2])
    tile = tile_at(path)
    result = mosaic([tile], extent_of(tile), 0.8)
    assert result.failures == ()
    alpha = result.raster.data[..., 3]
    assert (alpha[:, :10] == 0).all() and (alpha[:, 10:] == 255).all()
    np.testing.assert_array_equal(
        result.raster.data, reference([tile], extent_of(tile), 0.8, level=0)
    )


@pytest.mark.parametrize("change", ["removed", "unreadable"])
def test_local_changes_after_the_header_pass_raise_their_own_error(
    tmp_path, monkeypatch, change
):
    if change == "unreadable" and CANNOT_DENY_READ:
        pytest.skip("file permissions cannot deny reads")
    tile = tile_at(grid_tiff(tmp_path / "t.tif", 10, 10, 0.4))

    def changing_open(path, *args, **kwargs):
        dataset = RASTERIO_OPEN(path, *args, **kwargs)
        if change == "removed":
            Path(path).unlink(missing_ok=True)
        else:
            Path(path).chmod(0)
        return dataset

    monkeypatch.setattr(rasterio, "open", changing_open)
    error = FileNotFoundError if change == "removed" else PermissionError
    try:
        with pytest.raises(error) as info:
            mosaic([tile], extent_of(tile))
    finally:
        if tile.path.exists():
            tile.path.chmod(0o644)
    assert not isinstance(info.value, RasterTileReadError)


@pytest.mark.parametrize("kind", ["missing", "directory", "unreadable"])
def test_local_file_errors_propagate(tmp_path, kind):
    good = tile_at(grid_tiff(tmp_path / "good.tif", 10, 10, 0.4))
    path = tmp_path / "local.tif"
    if kind == "directory":
        path.mkdir()
    elif kind == "unreadable":
        if CANNOT_DENY_READ:
            pytest.skip("file permissions cannot deny reads")
        grid_tiff(path, 10, 10, 0.4)
        path.chmod(0)
    tile = make_tile(path, extent=extent_of(good), id="local")
    error = {
        "missing": FileNotFoundError,
        "directory": IsADirectoryError,
        "unreadable": PermissionError,
    }[kind]
    try:
        with pytest.raises(error) as info:
            mosaic([good, tile], extent_of(good), 0.4)
    finally:
        if kind == "unreadable":
            path.chmod(0o644)
    assert not isinstance(info.value, RasterTileReadError)


def test_memory_error_from_output_allocation_propagates(tmp_path, monkeypatch):
    tile = tile_at(grid_tiff(tmp_path / "t.tif", 10, 10, 0.4))
    real_zeros = np.zeros

    def failing_zeros(shape, *args, **kwargs):
        if tuple(np.atleast_1d(shape)) == (10, 10, 4):
            raise MemoryError("no memory")
        return real_zeros(shape, *args, **kwargs)

    monkeypatch.setattr(np, "zeros", failing_zeros)
    with pytest.raises(MemoryError, match="no memory"):
        mosaic([tile], extent_of(tile))


def test_unexpected_mosaic_read_exception_propagates(tmp_path, monkeypatch):
    tile = tile_at(grid_tiff(tmp_path / "t.tif", 10, 10, 0.4))

    def failing_read(self, *args, **kwargs):
        raise RuntimeError("unexpected")

    monkeypatch.setattr(rasterio.io.DatasetReader, "read", failing_read)
    with pytest.raises(RuntimeError, match="unexpected") as info:
        mosaic([tile], extent_of(tile))
    assert type(info.value) is RuntimeError


# Budget


def test_over_budget_is_refused_before_allocation_and_reads(
    tmp_path, monkeypatch, no_band_reads
):
    tile = tile_at(grid_tiff(tmp_path / "t.tif", 20, 20, 0.4))
    bounds = extent_of(tile)
    estimate = estimate_of([tile], bounds)
    forbid_allocation(monkeypatch)
    with pytest.raises(WorkingMemoryError) as info:
        mosaic([tile], bounds, budget=estimate - 1)
    assert info.value.estimate_bytes == estimate
    assert info.value.budget_bytes == estimate - 1
    message = str(info.value)
    assert "resolution" in message and 'product="tiles"' in message
    assert "max_memory_bytes" in message


def test_estimate_counts_output_arena_and_gdal_allowance(tmp_path):
    tile = tile_at(grid_tiff(tmp_path / "t.tif", 20, 20, 0.4))
    bounds = extent_of(tile)
    assert estimate_of([tile], bounds) == (
        20 * 20 * 4 + block_bytes(20, 20, 3, (20, 20)) + READ_OVERHEAD_BYTES
    )


def test_budget_equal_to_the_mosaic_estimate_is_accepted(tmp_path):
    tile = tile_at(grid_tiff(tmp_path / "t.tif", 20, 20, 0.4))
    bounds = extent_of(tile)
    estimate = estimate_of([tile], bounds)
    assert mosaic([tile], bounds, budget=estimate).valid_pixels == 400


def test_small_request_in_a_large_tile_reads_small_windows(tmp_path, reads, opened):
    path = tmp_path / "full.tif"
    with RASTERIO_OPEN(
        path, "w", driver="GTiff", width=15625, height=15625, count=4, dtype="uint8",
        crs="EPSG:3006", transform=from_bounds(*FULL_CELL, 15625, 15625),
        tiled=True, blockxsize=512, blockysize=512, sparse_ok=True,
        alpha="unspecified",
    ) as dst:
        dst.colorinterp = LM_COLORS + [ColorInterp.undefined]
    tile = make_tile(path, "rgbi", extent=FULL_CELL, size_bytes=1)
    x, y = FULL_CELL[0] + 1000.0, FULL_CELL[1] + 1000.0
    bounds = (x, y, x + 10.0, y + 10.0)
    result = mosaic([tile], bounds, budget=READ_OVERHEAD_BYTES + 1024**2)
    assert result.resolution == 0.16
    assert result.raster.data.shape == (63, 63, 4)
    assert (result.raster.data[..., 3] == 255).all()
    assert reads
    for call in reads:
        window = call["kwargs"]["window"]
        assert window.width <= 65 and window.height <= 65
        assert call["cache"] == READ_OVERHEAD_BYTES
    assert all(call["readdir"] == "EMPTY_DIR" for call in opened)


def test_reads_reuse_one_arena_and_restore_gdal_settings(
    tmp_path, monkeypatch, reads, opened
):
    tile = tile_at(grid_tiff(tmp_path / "t.tif", 50, 40, 0.4))
    monkeypatch.setattr(raster_tiles, "STRIP_BYTES", 600)
    before = (
        get_gdal_config("GDAL_CACHEMAX"),
        get_gdal_config("GDAL_DISABLE_READDIR_ON_OPEN"),
    )
    mosaic([tile], extent_of(tile), 0.4)
    assert len(reads) > 1
    roots = {id(arena_root(call["kwargs"]["out"])) for call in reads}
    assert len(roots) == 1
    assert all(call["cache"] == READ_OVERHEAD_BYTES for call in reads)
    assert all(call["cache"] == READ_OVERHEAD_BYTES for call in opened)
    assert all(call["readdir"] == "EMPTY_DIR" for call in opened)
    assert all(call["kwargs"].get("driver") == "GTiff" for call in opened)
    assert (
        get_gdal_config("GDAL_CACHEMAX"),
        get_gdal_config("GDAL_DISABLE_READDIR_ON_OPEN"),
    ) == before


def test_mosaic_ignores_sidecar_overviews(tmp_path):
    tile = tile_at(grid_tiff(tmp_path / "t.tif", 40, 40, 0.4))
    with RASTERIO_OPEN(f"{tile.path}.ovr", "w", driver="GTiff", width=20, height=20,
                       count=3, dtype="uint8") as dst:
        dst.write(np.full((3, 20, 20), 99, dtype=np.uint8))
    with RASTERIO_OPEN(tile.path) as src:
        assert src.overviews(1) == [2]
    bounds = extent_of(tile)
    data = mosaic([tile], bounds, 0.8).raster.data
    assert 99 not in data[..., 2]


# Side effects


def test_mosaic_is_local_and_leaves_files_unchanged(tmp_path, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("mosaics must not use the network")

    monkeypatch.setattr(requests, "get", forbidden)
    monkeypatch.setattr(requests.Session, "request", forbidden)
    tile = tile_at(grid_tiff(tmp_path / "t.tif", 10, 10, 0.4, nodata=0))
    before = (tile.path.read_bytes(), tile.path.stat().st_mtime_ns)
    mosaic([tile], extent_of(tile))
    assert (tile.path.read_bytes(), tile.path.stat().st_mtime_ns) == before


# Progress


def mosaic_reporting(tiles, bounds, progress, resolution=None):
    return build_mosaic(
        tiles, bounds, resolution=resolution,
        max_memory_bytes=DEFAULT_MAX_MEMORY_BYTES, progress=progress,
    )


def test_header_progress_counts_only_candidates(tmp_path):
    older = tile_at(grid_tiff(tmp_path / "old.tif", 40, 40, 0.4, tag=1), year=2010)
    newer = tile_at(grid_tiff(tmp_path / "new.tif", 40, 40, 0.4, tag=2), year=2020)
    away = tile_at(
        grid_tiff(tmp_path / "away.tif", 10, 10, 0.4,
                  origin=(ORIGIN[0] + 1000.0, ORIGIN[1])),
        year=2015,
    )
    tiles = [older, newer, bad_tile(tmp_path, "cir", older), away]
    events = []
    result = mosaic_reporting(
        tiles, extent_of(older), lambda *event: events.append(event)
    )
    assert len(result.failures) == 1
    assert events[:3] == [("headers", 1, 3), ("headers", 2, 3), ("headers", 3, 3)]
    assert events[3][:2] == ("mosaic", 0)


def test_mosaic_progress_reports_each_planned_block_after_its_read(
    tmp_path, monkeypatch
):
    older = tile_at(grid_tiff(tmp_path / "old.tif", 40, 40, 0.4, tag=1), year=2010)
    newer = tile_at(grid_tiff(tmp_path / "new.tif", 40, 40, 0.4, tag=2), year=2020)
    # Below one 40-pixel row of buffers: blocks split in both directions.
    monkeypatch.setattr(raster_tiles, "STRIP_BYTES", 600)
    log = []
    real_read = rasterio.io.DatasetReader.read

    def recording_read(self, *args, **kwargs):
        log.append("read")
        return real_read(self, *args, **kwargs)

    monkeypatch.setattr(rasterio.io.DatasetReader, "read", recording_read)
    mosaic_reporting([older, newer], extent_of(older), lambda *event: log.append(event))
    blocks = [entry for entry in log if entry != "read" and entry[0] == "mosaic"]
    total = blocks[0][2]
    assert total > 2
    assert blocks == [("mosaic", done, total) for done in range(total + 1)]
    after_start = log[log.index(blocks[0]) + 1 :]
    assert after_start == [entry for pair in zip(["read"] * total, blocks[1:])
                           for entry in pair]


def test_mosaic_progress_skips_a_failed_tiles_remaining_blocks(tmp_path, monkeypatch):
    older = tile_at(grid_tiff(tmp_path / "old.tif", 40, 40, 0.4, tag=1), year=2010)
    path = grid_tiff(
        tmp_path / "new.tif", 40, 40, 0.4, tag=2, tiled=True, blockxsize=16,
        blockysize=16, compress="deflate",
    )
    damage_block(path, 0)  # read first, so the tile's later blocks are skipped
    newer = tile_at(path, year=2020)
    monkeypatch.setattr(raster_tiles, "STRIP_BYTES", 5000)
    reads = []
    real_read = rasterio.io.DatasetReader.read

    def counting_read(self, *args, **kwargs):
        reads.append(1)
        return real_read(self, *args, **kwargs)

    monkeypatch.setattr(rasterio.io.DatasetReader, "read", counting_read)
    events = []
    result = mosaic_reporting([older, newer], extent_of(older),
                              lambda *event: events.append(event))
    assert len(result.failures) == 1
    done = [event[1] for event in events if event[0] == "mosaic"]
    total = events[-1][2]
    assert done == sorted(done) and done[-1] == total
    assert len(reads) < total


def test_no_admitted_tile_reports_headers_only(tmp_path):
    refused = bad_tile(tmp_path, "cir", None)
    events = []
    result = mosaic_reporting([refused], extent_of(refused),
                              lambda *event: events.append(event))
    assert result.valid_pixels == 0
    assert events == [("headers", 1, 1)]


CALLBACK_ERRORS = {
    "runtime": lambda: RuntimeError("subscriber"),
    "layout": lambda: RasterTileLayoutError(Path("x.tif"), "rgb", "subscriber"),
    "read": lambda: RasterTileReadError(Path("x.tif"), "subscriber"),
    "rasterio": lambda: rasterio.errors.RasterioError("subscriber"),
}


@pytest.mark.parametrize("kind", CALLBACK_ERRORS, ids=list(CALLBACK_ERRORS))
@pytest.mark.parametrize("stage", ["headers", "mosaic"])
def test_progress_callback_errors_propagate_unchanged(tmp_path, stage, kind):
    tile = tile_at(grid_tiff(tmp_path / "t.tif", 40, 40, 0.4))
    error = CALLBACK_ERRORS[kind]()

    def progress(event_stage, done, total):
        if event_stage == stage and done >= 1:
            raise error

    with pytest.raises(type(error)) as info:
        mosaic_reporting([tile], extent_of(tile), progress)
    assert info.value is error
