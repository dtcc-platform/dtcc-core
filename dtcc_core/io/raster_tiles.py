"""Budgeted loading of file-backed raster tiles into Rasters."""

from pathlib import Path

import numpy as np
import rasterio
from rasterio.enums import ColorInterp, MaskFlags

from ..model import Raster
from ..model.values.raster_tiles import RasterTile

DEFAULT_MAX_MEMORY_BYTES = 2 * 1024**3
# Allowance for GDAL's block cache and read buffers, added to every estimate.
# GDAL_CACHEMAX is capped to it from opening a file to its last read.
READ_OVERHEAD_BYTES = 64 * 1024**2
LOAD_LAYOUTS = {"rgb": 3, "rgbi": 4, "cir": 3}
CRS = "EPSG:3006"
BOUNDS_TOLERANCE = 0.001


class RasterTileLayoutError(ValueError):
    """A tile file's layout is not admitted, or its header contradicts the tile."""

    def __init__(self, path, spektraltyp, reason):
        self.path = Path(path)
        self.spektraltyp = spektraltyp
        self.reason = reason
        super().__init__(f"Cannot load {spektraltyp!r} tile {path}: {reason}")


class WorkingMemoryError(RuntimeError):
    """The estimated working memory exceeds the budget."""

    def __init__(self, estimate_bytes, budget_bytes, path):
        self.estimate_bytes = estimate_bytes
        self.budget_bytes = budget_bytes
        super().__init__(
            f"Loading {path} needs an estimated {estimate_bytes} bytes of working "
            f"memory, above the budget of {budget_bytes} bytes. Pass a larger "
            "max_memory_bytes, or read windows of the tile's path with Rasterio."
        )


class RasterTileReadError(OSError):
    """Rasterio cannot open or decode a tile file that the filesystem can open."""

    def __init__(self, path, message):
        self.path = Path(path)
        super().__init__(f"Cannot read raster tile {path}: {message}")


def load_tile(tile: RasterTile, *, max_memory_bytes: int) -> Raster:
    """Load every band of a local tile file into an (H, W, C) Raster.

    The header is admitted and the working memory estimated before any pixel is
    read. Local filesystem errors propagate unchanged; Rasterio failures to open
    or decode the file raise RasterTileReadError. The file is never modified and
    no network access is made.
    """
    if not isinstance(max_memory_bytes, int) or isinstance(max_memory_bytes, bool):
        kind = type(max_memory_bytes).__name__
        raise TypeError(f"max_memory_bytes must be an integer, not {kind}")
    if max_memory_bytes <= 0:
        raise ValueError(f"max_memory_bytes must be positive, not {max_memory_bytes}")
    path = Path(tile.path)
    if tile.spektraltyp == "pan":
        raise NotImplementedError(
            f"Loading pan tile {tile.collection}/{tile.id} is not supported; "
            f"read {path} directly."
        )
    # Raise missing files, directories and permission problems as their own
    # OSError before Rasterio sees the path.
    with open(path, "rb"):
        pass
    with rasterio.Env(GDAL_CACHEMAX=READ_OVERHEAD_BYTES):
        try:
            source = rasterio.open(path, driver="GTiff")
        except rasterio.errors.RasterioError as error:
            raise RasterTileReadError(path, str(error)) from error
        with source:
            nodata = _admit(tile, path, source)
            height, width, bands = source.height, source.width, source.count
            itemsize = np.dtype(source.dtypes[0]).itemsize
            estimate = (
                height * width * bands * itemsize
                + height * width * itemsize
                + READ_OVERHEAD_BYTES
            )
            if estimate > max_memory_bytes:
                raise WorkingMemoryError(estimate, max_memory_bytes, path)
            data = np.empty((height, width, bands), dtype=source.dtypes[0])
            buffer = np.empty((height, width), dtype=source.dtypes[0])
            for band in range(1, bands + 1):
                try:
                    source.read(band, out=buffer)
                except rasterio.errors.RasterioError as error:
                    raise RasterTileReadError(path, str(error)) from error
                data[:, :, band - 1] = buffer
            transform = source.transform
    return Raster(data=data, georef=transform, crs=CRS, nodata=nodata)


def _admit(tile: RasterTile, path: Path, source) -> float:
    """Check the header against the admitted layouts and the tile; return nodata."""

    def refuse(reason):
        return RasterTileLayoutError(path, tile.spektraltyp, reason)

    expected = LOAD_LAYOUTS.get(tile.spektraltyp)
    if expected is None:
        raise refuse(
            f"unsupported spektraltyp, expected one of {sorted(LOAD_LAYOUTS)}"
        )
    if source.count != expected:
        raise refuse(f"{source.count} bands, expected {expected}")
    if any(dtype != "uint8" for dtype in source.dtypes):
        raise refuse(f"band types {source.dtypes}, expected uint8")
    if ColorInterp.alpha in source.colorinterp:
        raise refuse("a band is labelled alpha")
    # Only masks derived from nodata, or none, are admitted. GDAL reports an
    # explicit per-band mask (such as a .msk sidecar) with no flags at all.
    if any(
        list(flags) not in ([MaskFlags.all_valid], [MaskFlags.nodata])
        for flags in source.mask_flag_enums
    ):
        raise refuse(f"the file has a mask band: {source.mask_flag_enums}")
    nodatavals = tuple(source.nodatavals)
    if None in nodatavals and set(nodatavals) != {None}:
        raise refuse(f"nodata declared on some bands only: {nodatavals}")
    if len(set(nodatavals)) != 1:
        raise refuse(f"nodata differs between bands: {nodatavals}")
    epsg = source.crs.to_epsg() if source.crs is not None else None
    if epsg != 3006 or tile.crs != CRS:
        raise refuse(f"CRS is EPSG:{epsg} for a {tile.crs} tile, expected {CRS}")
    extent = tile.extent
    expected_bounds = (extent.xmin, extent.ymin, extent.xmax, extent.ymax)
    # Written as "all within" so that NaN bounds fail the check.
    if not all(
        abs(actual - wanted) <= BOUNDS_TOLERANCE
        for actual, wanted in zip(source.bounds, expected_bounds)
    ):
        raise refuse(
            f"bounds {tuple(source.bounds)} do not match the tile extent "
            f"{expected_bounds}"
        )
    return np.nan if nodatavals[0] is None else float(nodatavals[0])
