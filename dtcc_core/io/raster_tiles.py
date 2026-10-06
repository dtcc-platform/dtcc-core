"""Budgeted loading and mosaicking of file-backed raster tiles into Rasters."""

import math
from collections.abc import Sequence
from dataclasses import dataclass
from fractions import Fraction
from numbers import Real
from pathlib import Path

import numpy as np
import rasterio
from affine import Affine
from rasterio.enums import ColorInterp, MaskFlags
from rasterio.windows import Window

from ..model import Raster
from ..model.values.raster_tiles import RasterTile

DEFAULT_MAX_MEMORY_BYTES = 2 * 1024**3
# Allowance for GDAL's block cache and read buffers, added to every estimate.
# GDAL_CACHEMAX is capped to it from opening a file to its last read.
READ_OVERHEAD_BYTES = 64 * 1024**2
LOAD_LAYOUTS = {"rgb": 3, "rgbi": 4, "cir": 3}
MOSAIC_LAYOUTS = {"rgb": 3, "rgbi": 4}
# Target size of one mosaic block's working memory (see build_mosaic).
STRIP_BYTES = 32 * 1024**2
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

    def __init__(self, estimate_bytes, budget_bytes, action, advice):
        self.estimate_bytes = estimate_bytes
        self.budget_bytes = budget_bytes
        super().__init__(
            f"{action} needs an estimated {estimate_bytes} bytes of working "
            f"memory, above the budget of {budget_bytes} bytes. {advice}"
        )


class RasterTileReadError(OSError):
    """Rasterio cannot open or decode a tile file that the filesystem can open."""

    def __init__(self, path, message):
        self.path = Path(path)
        super().__init__(f"Cannot read raster tile {path}: {message}")


def load_tile(tile: RasterTile, *, max_memory_bytes: int) -> Raster:
    """Load every band of a local tile file into an (H, W, C) Raster.

    Only the file itself is read; GDAL sidecar files beside it are ignored. The
    header is admitted and the working memory estimated before any pixel is
    read. Local filesystem errors propagate unchanged; Rasterio failures to open
    or decode the file raise RasterTileReadError. The file is never modified and
    no network access is made.
    """
    _check_budget(max_memory_bytes)
    path = Path(tile.path)
    if tile.spektraltyp == "pan":
        raise NotImplementedError(
            f"Loading pan tile {tile.collection}/{tile.id} is not supported; "
            f"read {path} directly."
        )
    _preflight(path)
    with _reading():
        with _open_source(path) as source:
            nodata = _admit(tile, path, source, LOAD_LAYOUTS)
            height, width, bands = source.height, source.width, source.count
            itemsize = np.dtype(source.dtypes[0]).itemsize
            estimate = (
                height * width * bands * itemsize
                + height * width * itemsize
                + READ_OVERHEAD_BYTES
            )
            if estimate > max_memory_bytes:
                raise WorkingMemoryError(
                    estimate,
                    max_memory_bytes,
                    f"Loading {path}",
                    "Pass a larger max_memory_bytes, or read windows of the tile's "
                    "path with Rasterio.",
                )
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


def _check_budget(max_memory_bytes) -> None:
    if not isinstance(max_memory_bytes, int) or isinstance(max_memory_bytes, bool):
        kind = type(max_memory_bytes).__name__
        raise TypeError(f"max_memory_bytes must be an integer, not {kind}")
    if max_memory_bytes <= 0:
        raise ValueError(f"max_memory_bytes must be positive, not {max_memory_bytes}")


def _preflight(path: Path) -> None:
    """Raise missing files, directories and permission problems as their own
    OSError before Rasterio sees the path."""
    with open(path, "rb"):
        pass


def _reading():
    """GDAL settings for every tile open and read.

    GDAL's block cache is capped to READ_OVERHEAD_BYTES, and EMPTY_DIR keeps GDAL
    from reading sidecar files (.msk, .ovr, .aux.xml) next to a tile, so only the
    original file's bytes are read.
    """
    return rasterio.Env(
        GDAL_CACHEMAX=READ_OVERHEAD_BYTES, GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR"
    )


def _open_source(path: Path, **options):
    """Open with the GeoTIFF driver only; Rasterio failures are read errors."""
    try:
        return rasterio.open(path, driver="GTiff", **options)
    except rasterio.errors.RasterioError as error:
        raise RasterTileReadError(path, str(error)) from error


def _admit(tile: RasterTile, path: Path, source, layouts) -> float:
    """Check the header against the admitted layouts and the tile; return nodata."""

    def refuse(reason):
        return RasterTileLayoutError(path, tile.spektraltyp, reason)

    expected = layouts.get(tile.spektraltyp)
    if expected is None:
        raise refuse(f"unsupported spektraltyp, expected one of {sorted(layouts)}")
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


@dataclass(frozen=True)
class TileFailure:
    """A tile left out of a mosaic, or only partly used, because of its source.

    ``pixels`` counts the output pixels it contributed before failing.
    """

    tile: RasterTile
    error: RasterTileLayoutError | RasterTileReadError
    pixels: int


@dataclass(frozen=True)
class Mosaic:
    """An RGBA mosaic of local tiles and how it was made.

    ``bounds`` (the output grid) and ``resolution`` are None when no pixel is
    valid; ``raster`` is then the empty Raster. ``sources`` lists the tiles that
    contributed pixels, with their pixel counts, in precedence order.
    """

    raster: Raster
    requested_bounds: tuple[float, float, float, float]
    bounds: tuple[float, float, float, float] | None
    resolution: float | None
    sources: tuple[tuple[RasterTile, int], ...]
    failures: tuple[TileFailure, ...]
    valid_pixels: int


@dataclass(frozen=True)
class _Level:
    """The full resolution (index None) or an internal overview level."""

    index: int | None
    transform: Affine
    width: int
    height: int


@dataclass(frozen=True)
class _Axis:
    """Output pixel ``k`` samples level pixel ``(a + k * b) // d`` on one axis."""

    a: int
    b: int
    d: int

    @classmethod
    def of(cls, offset: Fraction, resolution: Fraction, spacing: Fraction):
        start = (offset + resolution / 2) / spacing
        step = resolution / spacing
        d = math.lcm(start.denominator, step.denominator)
        return cls(
            start.numerator * (d // start.denominator),
            step.numerator * (d // step.denominator),
            d,
        )

    def span(self, size: int, count: int) -> tuple[int, int]:
        """Output pixels in [0, count) whose index falls in [0, size)."""
        first = math.ceil(Fraction(-self.a, self.b))
        stop = math.ceil(Fraction(size * self.d - self.a, self.b))
        return max(first, 0), min(stop, count)

    def fill(self, out, first: int, size: int) -> None:
        """Write the indices of output pixels first, first + 1, ... into out."""
        last = size - 1
        for k in range(len(out)):
            index = (self.a + (first + k) * self.b) // self.d
            out[k] = 0 if index < 0 else last if index > last else index


@dataclass(frozen=True)
class _Plan:
    tile: RasterTile
    path: Path
    base: _Level
    level: _Level
    bands: int
    rows: tuple[int, int]
    cols: tuple[int, int]
    row_axis: _Axis
    col_axis: _Axis
    block: tuple[int, int]
    block_bytes: int


def _check_request(bounds, resolution, required: bool):
    if len(bounds) != 4:
        raise ValueError(f"bounds must have four values, not {len(bounds)}")
    for value in bounds:
        if isinstance(value, bool) or not isinstance(value, Real):
            raise TypeError(f"bounds must be real numbers, not {value!r}")
    xmin, ymin, xmax, ymax = (float(value) for value in bounds)
    if not all(math.isfinite(value) for value in (xmin, ymin, xmax, ymax)):
        raise ValueError(f"bounds must be finite, not {bounds}")
    if not (xmin < xmax and ymin < ymax):
        raise ValueError(f"bounds must have xmin < xmax and ymin < ymax: {bounds}")
    if resolution is not None or required:
        if isinstance(resolution, bool) or not isinstance(resolution, Real):
            raise TypeError(f"resolution must be a number, not {resolution!r}")
        if not (math.isfinite(resolution) and resolution > 0):
            raise ValueError(f"resolution must be finite and positive: {resolution}")
    return xmin, ymin, xmax, ymax


def mosaic_grid(bounds, resolution) -> tuple[int, int, Affine]:
    """Width, height and transform of the output grid for bounds and resolution.

    The grid is anchored at the top-left corner of the bounds and rounded up to
    whole pixels (exactly, over the given float values), so it always covers
    the bounds and has at least one pixel.
    """
    xmin, ymin, xmax, ymax = _check_request(bounds, resolution, required=True)
    r = Fraction(resolution)
    width = math.ceil((Fraction(xmax) - Fraction(xmin)) / r)
    height = math.ceil((Fraction(ymax) - Fraction(ymin)) / r)
    r = float(resolution)
    return width, height, Affine(r, 0.0, xmin, 0.0, -r, ymax)


def mosaic_minimum_bytes(bounds, resolution) -> int:
    """A lower bound on the working memory of any mosaic of this grid."""
    width, height, _ = mosaic_grid(bounds, resolution)
    return width * height * 4 + READ_OVERHEAD_BYTES


def build_mosaic(
    tiles: Sequence[RasterTile],
    bounds: tuple[float, float, float, float],
    *,
    resolution: float | None,
    max_memory_bytes: int,
) -> Mosaic:
    """Mosaic local rgb and rgbi tiles into an RGBA Raster over bounds.

    The newest tile, by (datetime, collection, id), gives each pixel; older tiles
    fill pixels it leaves empty. A pixel of a tile with declared nodata is empty
    when every band, the fourth (NIR) included, equals it; without declared
    nodata every pixel is valid. RGB comes from bands 1-3 and alpha is 255 or 0.

    Each output pixel takes the sample of the source pixel containing its
    centre (a centre on a pixel edge goes to the right or lower pixel), from the
    coarsest of the full resolution and the internal overview levels whose
    spacing is at most the output resolution. The resolution is the finest
    admitted file's unless given; the mosaic never falls back to a coarser one.

    Headers are admitted and the working memory (the RGBA output, the largest
    block of read buffers and READ_OVERHEAD_BYTES) checked before any pixel is
    read. Tiles whose files are not admitted or cannot be decoded are recorded
    as failures; local filesystem errors propagate. Without any valid pixel the
    result is the empty Raster.
    """
    _check_request(bounds, resolution, required=False)
    _check_budget(max_memory_bytes)
    requested = tuple(float(value) for value in bounds)
    xmin, ymin, xmax, ymax = requested
    candidates = sorted(
        (
            tile
            for tile in tiles
            if tile.extent.xmin < xmax
            and tile.extent.xmax > xmin
            and tile.extent.ymin < ymax
            and tile.extent.ymax > ymin
        ),
        key=lambda tile: (tile.datetime, tile.collection, tile.id),
        reverse=True,
    )
    failures = []
    with _reading():
        admitted = []
        for tile in candidates:
            path = Path(tile.path)
            _preflight(path)
            try:
                with _open_source(path) as source:
                    _admit_mosaic(tile, path, source)
                    levels = _levels(path, source)
                    bands = source.count
            except (RasterTileLayoutError, RasterTileReadError) as error:
                failures.append(TileFailure(tile, error, 0))
                continue
            admitted.append((tile, path, levels, bands))
        if not admitted:
            return _empty(requested, failures)
        if resolution is None:
            spacings = (levels[0].transform.a for _, _, levels, _ in admitted)
            resolution = float(min(spacings, key=Fraction))
        width, height, transform = mosaic_grid(requested, resolution)
        plans = [
            plan
            for plan in (
                _plan(tile, path, levels, bands, requested, resolution, width, height)
                for tile, path, levels, bands in admitted
            )
            if plan is not None
        ]
        arena_bytes = max((plan.block_bytes for plan in plans), default=0)
        estimate = width * height * 4 + arena_bytes + READ_OVERHEAD_BYTES
        if estimate > max_memory_bytes:
            raise WorkingMemoryError(
                estimate,
                max_memory_bytes,
                f"A {width} x {height} mosaic at {resolution} m",
                'Use a coarser resolution, product="tiles", or a larger '
                "max_memory_bytes.",
            )
        output = np.zeros((height, width, 4), dtype=np.uint8)
        arena = np.empty(arena_bytes, dtype=np.uint8)
        sources = []
        for plan in plans:
            counter = [0]
            try:
                _composite(plan, output, arena, counter)
            except (RasterTileLayoutError, RasterTileReadError) as error:
                failures.append(TileFailure(plan.tile, error, counter[0]))
            if counter[0]:
                sources.append((plan.tile, counter[0]))
    valid_pixels = sum(pixels for _, pixels in sources)
    if not valid_pixels:
        return _empty(requested, failures)
    r = float(resolution)
    return Mosaic(
        raster=Raster(data=output, georef=transform, crs=CRS, nodata=np.nan),
        requested_bounds=requested,
        bounds=(xmin, ymax - height * r, xmin + width * r, ymax),
        resolution=r,
        sources=tuple(sources),
        failures=tuple(failures),
        valid_pixels=valid_pixels,
    )


def _empty(requested, failures) -> Mosaic:
    return Mosaic(
        raster=Raster(crs=CRS),
        requested_bounds=requested,
        bounds=None,
        resolution=None,
        sources=(),
        failures=tuple(failures),
        valid_pixels=0,
    )


def _admit_mosaic(tile: RasterTile, path: Path, source) -> float:
    nodata = _admit(tile, path, source, MOSAIC_LAYOUTS)
    t = source.transform
    if not (t.b == 0 and t.d == 0 and t.a > 0 and t.e < 0):
        raise RasterTileLayoutError(
            path, tile.spektraltyp, f"transform {tuple(t)[:6]} is not north-up"
        )
    if t.e != -t.a:
        raise RasterTileLayoutError(
            path, tile.spektraltyp, f"pixels are not square: {t.a} x {-t.e}"
        )
    return nodata


def _levels(path: Path, source) -> list[_Level]:
    levels = [_Level(None, source.transform, source.width, source.height)]
    for index in range(len(source.overviews(1))):
        with _open_source(path, overview_level=index) as overview:
            levels.append(
                _Level(index, overview.transform, overview.width, overview.height)
            )
    return levels


def _plan(tile, path, levels, bands, bounds, resolution, width, height):
    """The tile's output window, read level and block size; None if empty."""
    xmin, _, _, ymax = (Fraction(value) for value in bounds)
    r = Fraction(resolution)
    base = levels[0]
    cols = _Axis.of(xmin - Fraction(base.transform.c), r, Fraction(base.transform.a))
    rows = _Axis.of(
        Fraction(base.transform.f) - ymax, r, Fraction(-base.transform.e)
    )
    col_span = cols.span(base.width, width)
    row_span = rows.span(base.height, height)
    if col_span[0] >= col_span[1] or row_span[0] >= row_span[1]:
        return None
    fits = [
        level
        for level in levels[1:]
        if Fraction(level.transform.a) <= r and Fraction(-level.transform.e) <= r
    ]
    level = max(fits, key=lambda level: level.transform.a) if fits else base
    pixel_x, pixel_y = Fraction(level.transform.a), Fraction(-level.transform.e)
    ratio_x, ratio_y = r / pixel_x, r / pixel_y

    def block_bytes(block_rows, block_cols):
        nh = min(math.ceil(block_rows * ratio_y) + 2, level.height)
        nw = min(math.ceil(block_cols * ratio_x) + 2, level.width)
        return (
            bands * nh * nw
            + bands * block_rows * nw
            + bands * block_rows * block_cols
            + 2 * block_rows * block_cols
            + 8 * (block_rows + block_cols)
        )

    window_rows = row_span[1] - row_span[0]
    window_cols = col_span[1] - col_span[0]
    block_cols = _largest(window_cols, lambda n: block_bytes(1, n))
    block_rows = _largest(window_rows, lambda n: block_bytes(n, block_cols))
    return _Plan(
        tile=tile,
        path=path,
        base=base,
        level=level,
        bands=bands,
        rows=row_span,
        cols=col_span,
        row_axis=_Axis.of(Fraction(level.transform.f) - ymax, r, pixel_y),
        col_axis=_Axis.of(xmin - Fraction(level.transform.c), r, pixel_x),
        block=(block_rows, block_cols),
        block_bytes=block_bytes(block_rows, block_cols),
    )


def _largest(limit: int, cost) -> int:
    """The largest n in [1, limit] with cost(n) <= STRIP_BYTES, at least 1."""
    low, high = 1, limit
    while low < high:
        middle = (low + high + 1) // 2
        if cost(middle) <= STRIP_BYTES:
            low = middle
        else:
            high = middle - 1
    return low


def _carve(arena, offset: int, shape, dtype):
    """A contiguous array of shape and dtype in the arena at offset."""
    size = math.prod(shape) * np.dtype(dtype).itemsize
    return arena[offset : offset + size].view(dtype).reshape(shape), offset + size


def _composite(plan: _Plan, output, arena, counter) -> None:
    """Fill output pixels that are still empty from one tile, block by block.

    The handle that is read is admitted and checked against the plan first: the
    full resolution with the mosaic rules, or an overview level with the same
    rules except square pixels (overview spacing may differ between X and Y).
    counter[0] counts the pixels filled, including when a later block fails.
    """
    _preflight(plan.path)
    if plan.level.index is None:
        with _open_source(plan.path) as source:
            nodata = _admit_mosaic(plan.tile, plan.path, source)
            _check_unchanged(plan, source, plan.base)
            _blend(plan, source, _sample_value(nodata), output, arena, counter)
        return
    with _open_source(plan.path, overview_level=plan.level.index) as level:
        nodata = _admit(plan.tile, plan.path, level, MOSAIC_LAYOUTS)
        _check_unchanged(plan, level, plan.level)
        _blend(plan, level, _sample_value(nodata), output, arena, counter)


def _sample_value(nodata: float):
    """The declared nodata as a uint8 sample, or None if no sample can equal it."""
    if math.isnan(nodata) or nodata != int(nodata) or not 0 <= nodata <= 255:
        return None
    return np.uint8(nodata)


def _check_unchanged(plan: _Plan, dataset, expected: _Level) -> None:
    """Refuse a file whose geometry differs from the one the plan was made for."""
    actual = (dataset.transform, dataset.width, dataset.height, dataset.count)
    if actual != (expected.transform, expected.width, expected.height, plan.bands):
        raise RasterTileLayoutError(
            plan.path,
            plan.tile.spektraltyp,
            "the file changed while the mosaic was built: transform, size and "
            f"bands {actual}, expected {expected}",
        )


def _blend(plan: _Plan, source, nodata, output, arena, counter) -> None:
    block_rows, block_cols = plan.block
    for col in range(*plan.cols, block_cols):
        ncols = min(block_cols, plan.cols[1] - col)
        col_index, offset = _carve(arena, 0, (ncols,), np.intp)
        plan.col_axis.fill(col_index, col, plan.level.width)
        col_first = int(col_index[0])
        window_cols = int(col_index[-1]) - col_first + 1
        np.subtract(col_index, col_first, out=col_index)
        for row in range(*plan.rows, block_rows):
            nrows = min(block_rows, plan.rows[1] - row)
            row_index, at = _carve(arena, offset, (nrows,), np.intp)
            plan.row_axis.fill(row_index, row, plan.level.height)
            row_first = int(row_index[0])
            window_rows = int(row_index[-1]) - row_first + 1
            np.subtract(row_index, row_first, out=row_index)
            bands = plan.bands
            window, at = _carve(
                arena, at, (bands, window_rows, window_cols), np.uint8
            )
            gathered, at = _carve(arena, at, (bands, nrows, window_cols), np.uint8)
            samples, at = _carve(arena, at, (bands, nrows, ncols), np.uint8)
            valid, at = _carve(arena, at, (nrows, ncols), np.bool_)
            scratch, at = _carve(arena, at, (nrows, ncols), np.bool_)
            try:
                source.read(
                    out=window,
                    window=Window(col_first, row_first, window_cols, window_rows),
                )
            except rasterio.errors.RasterioError as error:
                raise RasterTileReadError(plan.path, str(error)) from error
            # mode="clip" with out= gathers without a temporary copy.
            np.take(window, row_index, axis=1, out=gathered, mode="clip")
            np.take(gathered, col_index, axis=2, out=samples, mode="clip")
            if nodata is None:
                valid.fill(True)
            else:
                np.not_equal(samples[0], nodata, out=valid)
                for band in range(1, bands):
                    np.not_equal(samples[band], nodata, out=scratch)
                    np.logical_or(valid, scratch, out=valid)
            target = output[row : row + nrows, col : col + ncols]
            np.equal(target[..., 3], 0, out=scratch)
            np.logical_and(valid, scratch, out=valid)
            for band in range(3):
                np.copyto(target[..., band], samples[band], where=valid)
            np.copyto(target[..., 3], 255, where=valid)
            counter[0] += int(np.count_nonzero(valid))
