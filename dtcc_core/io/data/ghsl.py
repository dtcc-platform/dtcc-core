"""Read GHSL built-up surface and settlement class grids from the JRC open data server."""

from __future__ import annotations

import math
import os
from pathlib import Path

import numpy as np
from affine import Affine

from ...model import Bounds, Raster
from .cache import cache_dir
from .logging import info, warning

try:
    import rasterio
    from rasterio.errors import RasterioIOError
    from rasterio.warp import Resampling, reproject, transform_bounds
    from rasterio.windows import from_bounds

    HAS_RASTERIO = True
except ImportError:
    HAS_RASTERIO = False


GHSL_EPOCHS = tuple(range(1975, 2031, 5))
GHSL_DEFAULT_EPOCH = 2020

BUILT_UP_EPOCHS = GHSL_EPOCHS
BUILT_UP_DEFAULT_EPOCH = GHSL_DEFAULT_EPOCH
BUILT_UP_COMPONENTS = ("total", "non_residential")
BUILT_UP_CELL_SIZE = 100.0
BUILT_UP_DOI = "https://doi.org/10.2905/9F06F36F-4B11-47EC-ABB0-4F8B7B1D72EA"
BUILT_UP_CITATION = (
    "Pesaresi, M., Politis, P. (2023): GHS-BUILT-S R2023A - GHS built-up surface "
    "grid, derived from Sentinel2 composite and Landsat, multitemporal (1975-2030). "
    "European Commission, Joint Research Centre (JRC). " + BUILT_UP_DOI
)
# Global 100 m Mollweide grids, one zipped GeoTIFF per epoch and component.
# GDAL reads the requested window directly from the remote zip.
BUILT_UP_URL = (
    "zip+https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_BUILT_S_GLOBE_R2023A/"
    "{name}/V1-0/{name}_V1_0.zip!{name}_V1_0.tif"
)
_BUILT_UP_NAMES = {
    "total": "GHS_BUILT_S_E{epoch}_GLOBE_R2023A_54009_100",
    "non_residential": "GHS_BUILT_S_NRES_E{epoch}_GLOBE_R2023A_54009_100",
}

SETTLEMENT_CLASS_EPOCHS = GHSL_EPOCHS
SETTLEMENT_CLASS_DEFAULT_EPOCH = GHSL_DEFAULT_EPOCH
SETTLEMENT_CLASS_CELL_SIZE = 1000.0
SETTLEMENT_CLASS_NODATA = -200
# Degree of Urbanisation level 2 classes used by GHS-SMOD.
SETTLEMENT_CLASSES = {
    30: "urban centre",
    23: "dense urban cluster",
    22: "semi-dense urban cluster",
    21: "suburban or peri-urban",
    13: "rural cluster",
    12: "low density rural",
    11: "very low density rural",
    10: "water",
}
SETTLEMENT_CLASS_DOI = "https://doi.org/10.2905/A0DF7A6F-49DE-46EA-9BDE-563437A6E2BA"
SETTLEMENT_CLASS_CITATION = (
    "Schiavina, M., Melchiorri, M., Pesaresi, M. (2023): GHS-SMOD R2023A - GHS "
    "settlement layers, application of the Degree of Urbanisation methodology "
    "(stage I) to GHS-POP R2023A and GHS-BUILT-S R2023A, multitemporal (1975-2030). "
    "European Commission, Joint Research Centre (JRC). " + SETTLEMENT_CLASS_DOI
)
# Global 1 km Mollweide grids, one zipped GeoTIFF per epoch. Version 2.0
# replaces the obsolete 1.0 files.
SETTLEMENT_CLASS_URL = (
    "zip+https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_SMOD_GLOBE_R2023A/"
    "{name}/V2-0/{name}_V2_0.zip!{name}_V2_0.tif"
)
_SETTLEMENT_CLASS_NAME = "GHS_SMOD_E{epoch}_GLOBE_R2023A_54009_1000"

_GDAL_ENV = {
    "GDAL_DISABLE_READDIR_ON_OPEN": "EMPTY_DIR",
    "CPL_VSIL_CURL_ALLOWED_EXTENSIONS": ".zip,.tif",
    "GDAL_HTTP_TIMEOUT": os.environ.get("DTCC_GHSL_HTTP_TIMEOUT", "120"),
    "GDAL_HTTP_MAX_RETRY": "3",
    "GDAL_HTTP_RETRY_DELAY": "2",
}
# Source cells read around the requested area so reprojection covers its edges.
_WINDOW_MARGIN_CELLS = 2


def download_built_up(
    bounds: Bounds,
    epoch: int = BUILT_UP_DEFAULT_EPOCH,
    component: str = "total",
    use_cache: bool = True,
) -> Raster:
    """Read GHSL built-up surface for ``bounds`` as a :class:`Raster`.

    Parameters
    ----------
    bounds : Bounds
        Area of interest in EPSG:3006.
    epoch : int, default 2020
        Epoch from 1975 to 2030 in steps of 5. 2025 and 2030 are projections.
    component : {"total", "non_residential"}, default "total"
        All built-up surface, or only the non-residential part.
    use_cache : bool, default True
        Reuse a cached result for the same bounds, epoch and component.

    Returns
    -------
    Raster
        Built-up surface in square metres per 100 m cell, on a grid anchored at
        the top-left corner of ``bounds`` in EPSG:3006. Reprojection from the
        Mollweide source grid preserves area totals. Cells without data are NaN.
    """
    if not HAS_RASTERIO:
        raise ImportError("rasterio is required to read GHSL data.")
    _check_epoch(epoch)
    if component not in BUILT_UP_COMPONENTS:
        raise ValueError(
            f"Unsupported GHSL component {component!r}. Supported: {BUILT_UP_COMPONENTS}."
        )
    bounds = _ensure_bounds(bounds)
    transform, width, height = _target_grid(bounds, BUILT_UP_CELL_SIZE)

    cache_path = _cache_path(f"built_up_{component}", bounds, epoch)
    data = _read_cache(cache_path) if use_cache else None
    if data is None:
        name = _BUILT_UP_NAMES[component].format(epoch=epoch)
        info(f"Reading GHSL {name} from the JRC server")
        values, source_transform, source_crs, nodata = _read_window(
            BUILT_UP_URL.format(name=name), bounds
        )
        source = values.astype(np.float32)
        if nodata is not None:
            source[values == nodata] = np.nan
        data = _reproject_preserving_totals(
            source, source_transform, source_crs, transform, width, height
        )
        _write_cache(cache_path, data, transform, np.nan)
    else:
        info(f"Using cached GHSL built-up {epoch} data from {cache_path}")

    if np.isnan(data).all():
        warning(f"No GHSL built-up {epoch} data in bounds {bounds.tuple}.")
    return Raster(data=data, georef=transform, nodata=np.nan, crs="EPSG:3006")


def download_settlement_class(
    bounds: Bounds, epoch: int = SETTLEMENT_CLASS_DEFAULT_EPOCH, use_cache: bool = True
) -> Raster:
    """Read the GHSL settlement class (GHS-SMOD) for ``bounds`` as a :class:`Raster`.

    Parameters
    ----------
    bounds : Bounds
        Area of interest in EPSG:3006.
    epoch : int, default 2020
        Epoch from 1975 to 2030 in steps of 5. 2025 and 2030 are projections.
    use_cache : bool, default True
        Reuse a cached result for the same bounds and epoch.

    Returns
    -------
    Raster
        Degree of Urbanisation class codes (see :data:`SETTLEMENT_CLASSES`) per
        1 km cell, on a grid anchored at the top-left corner of ``bounds`` in
        EPSG:3006, taken from the nearest Mollweide source cell. Cells without
        data hold :data:`SETTLEMENT_CLASS_NODATA`.
    """
    if not HAS_RASTERIO:
        raise ImportError("rasterio is required to read GHSL data.")
    _check_epoch(epoch)
    bounds = _ensure_bounds(bounds)
    transform, width, height = _target_grid(bounds, SETTLEMENT_CLASS_CELL_SIZE)

    cache_path = _cache_path("settlement_class", bounds, epoch)
    data = _read_cache(cache_path) if use_cache else None
    if data is None:
        name = _SETTLEMENT_CLASS_NAME.format(epoch=epoch)
        info(f"Reading GHSL {name} from the JRC server")
        values, source_transform, source_crs, nodata = _read_window(
            SETTLEMENT_CLASS_URL.format(name=name), bounds
        )
        source = values.astype(np.int16)
        if nodata is not None:
            source[values == nodata] = SETTLEMENT_CLASS_NODATA
        data = np.full((height, width), SETTLEMENT_CLASS_NODATA, dtype=np.int16)
        reproject(
            source,
            data,
            src_transform=source_transform,
            src_crs=source_crs,
            src_nodata=SETTLEMENT_CLASS_NODATA,
            dst_transform=transform,
            dst_crs="EPSG:3006",
            dst_nodata=SETTLEMENT_CLASS_NODATA,
            resampling=Resampling.nearest,
        )
        _write_cache(cache_path, data, transform, SETTLEMENT_CLASS_NODATA)
    else:
        info(f"Using cached GHSL settlement class {epoch} data from {cache_path}")

    if (data == SETTLEMENT_CLASS_NODATA).all():
        warning(f"No GHSL settlement class {epoch} data in bounds {bounds.tuple}.")
    return Raster(
        data=data, georef=transform, nodata=SETTLEMENT_CLASS_NODATA, crs="EPSG:3006"
    )


def _read_window(url: str, bounds: Bounds):
    """Return (values, transform, crs, nodata) of the source window covering ``bounds``."""
    try:
        with rasterio.Env(**_GDAL_ENV), rasterio.open(url) as src:
            left, bottom, right, top = transform_bounds(
                "EPSG:3006", src.crs, *bounds.tuple, densify_pts=21
            )
            margin = _WINDOW_MARGIN_CELLS * abs(src.res[0])
            window = from_bounds(
                left - margin, bottom - margin, right + margin, top + margin, src.transform
            ).round_offsets().round_lengths()
            values = src.read(1, window=window, boundless=True, fill_value=src.nodata)
            return values, src.window_transform(window), src.crs, src.nodata
    except RasterioIOError as exc:
        raise RuntimeError(f"Failed to read GHSL data from {url}: {exc}") from exc


def _reproject_preserving_totals(
    source, source_transform, source_crs, transform: Affine, width: int, height: int
):
    """Reproject m² per cell so that area totals are preserved.

    Values are converted to a built-up share, averaged onto the target grid
    and converted back to m² per target cell.
    """
    source_cell_area = abs(source_transform.a * source_transform.e)
    share = np.full((height, width), np.nan, dtype=np.float32)
    reproject(
        source / source_cell_area,
        share,
        src_transform=source_transform,
        src_crs=source_crs,
        src_nodata=np.nan,
        dst_transform=transform,
        dst_crs="EPSG:3006",
        dst_nodata=np.nan,
        resampling=Resampling.average,
    )
    return share * abs(transform.a * transform.e)


def _check_epoch(epoch: int):
    if epoch not in GHSL_EPOCHS:
        raise ValueError(f"Unsupported GHSL epoch {epoch}. Supported epochs: {GHSL_EPOCHS}.")


def _target_grid(bounds: Bounds, cell: float):
    width = max(1, math.ceil(round((bounds.xmax - bounds.xmin) / cell, 6)))
    height = max(1, math.ceil(round((bounds.ymax - bounds.ymin) / cell, 6)))
    transform = Affine(cell, 0, bounds.xmin, 0, -cell, bounds.ymax)
    return transform, width, height


def _read_cache(path: Path):
    if not path.exists():
        return None
    with rasterio.open(path) as src:
        return src.read(1)


def _write_cache(path: Path, data, transform: Affine, nodata):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_name(path.name + ".part.tif")
    with rasterio.open(
        tmp_path,
        "w",
        driver="GTiff",
        width=data.shape[1],
        height=data.shape[0],
        count=1,
        dtype=data.dtype,
        crs="EPSG:3006",
        transform=transform,
        nodata=nodata,
        compress="deflate",
    ) as dst:
        dst.write(data, 1)
    os.replace(tmp_path, path)


def _ensure_bounds(bounds) -> Bounds:
    if isinstance(bounds, Bounds):
        return bounds
    return Bounds(*bounds)


def _cache_path(product: str, bounds: Bounds, epoch: int) -> Path:
    xmin, ymin, xmax, ymax = (round(value) for value in bounds.tuple)
    return (
        Path(cache_dir)
        / "downloaded_ghsl"
        / f"{product}_{epoch}_{xmin}_{ymin}_{xmax}_{ymax}.tif"
    )


__all__ = [
    "BUILT_UP_CITATION",
    "BUILT_UP_COMPONENTS",
    "BUILT_UP_DEFAULT_EPOCH",
    "BUILT_UP_DOI",
    "BUILT_UP_EPOCHS",
    "SETTLEMENT_CLASSES",
    "SETTLEMENT_CLASS_CITATION",
    "SETTLEMENT_CLASS_DEFAULT_EPOCH",
    "SETTLEMENT_CLASS_DOI",
    "SETTLEMENT_CLASS_EPOCHS",
    "SETTLEMENT_CLASS_NODATA",
    "download_built_up",
    "download_settlement_class",
]
