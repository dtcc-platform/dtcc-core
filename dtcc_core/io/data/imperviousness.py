"""Download Copernicus HRL Imperviousness Density rasters from anonymous services."""

from __future__ import annotations

import io
import math
import os
from pathlib import Path

import numpy as np
import requests
from affine import Affine
from pyproj import Transformer

from ...model import Bounds, Raster
from .cache import cache_dir
from .logging import info, warning

try:
    import rasterio
    from rasterio.warp import Resampling, reproject

    HAS_RASTERIO = True
except ImportError:
    HAS_RASTERIO = False


IMPERVIOUSNESS_SUPPORTED_YEARS = (2018, 2021, 2024)
IMPERVIOUSNESS_DEFAULT_YEAR = 2024
IMPERVIOUSNESS_CELL_SIZE = 10.0
IMPERVIOUSNESS_DOI = {
    2018: "https://doi.org/10.2909/3bf542bd-eebd-4d73-b53c-a0243f2ed862",
    2021: "https://doi.org/10.2909/34ef6334-d432-4041-a3da-67e156d6501d",
    2024: "https://doi.org/10.2909/f0bcb6f0-d775-4218-bec1-13adea01ff8f",
}
IMPERVIOUSNESS_ATTRIBUTION = (
    "Generated using European Union's Copernicus Land Monitoring Service information; {doi}"
)
# 2018 is served as raw values by the EEA ArcGIS ImageServer, in any CRS.
IMPERVIOUSNESS_IMAGESERVER_URL = (
    "https://image.discomap.eea.europa.eu/arcgis/rest/services/GioLandPublic/"
    "HRL_ImperviousnessDensity_2018/ImageServer/exportImage"
)
# 2021 and 2024 are only served as WMS map images in EPSG:3035. The default
# style is a single-band palette whose pixel index is the density in percent.
IMPERVIOUSNESS_GEOSERVER_URL = "https://geoserver.geoville.com/geoserver/{workspace}/wms"
IMPERVIOUSNESS_GEOSERVER_LAYERS = {
    2021: ("nvlcc_2021", "HRL_NVLCC_IMD_10m"),
    2024: ("nvlcc_2024", "HRL_NVLCC_IMD_2024_10m"),
}
IMPERVIOUSNESS_CONNECT_TIMEOUT_SECONDS = float(
    os.environ.get("DTCC_IMPERVIOUSNESS_CONNECT_TIMEOUT", "10")
)
IMPERVIOUSNESS_READ_TIMEOUT_SECONDS = float(
    os.environ.get("DTCC_IMPERVIOUSNESS_READ_TIMEOUT", "120")
)
# Largest request per service, in pixels (ImageServer limit 15000 x 4100).
_IMAGESERVER_MAX_PIXELS = (4000, 4000)
_GEOSERVER_MAX_PIXELS = (2048, 2048)
# Valid densities are 0-100. Sea and areas outside the product are 255 (2018)
# or 101 (2021, 2024).
_MAX_DENSITY = 100
# Margin in metres around the EPSG:3035 request so reprojection covers the
# whole EPSG:3006 grid.
_REPROJECT_MARGIN_M = 50.0


def download_imperviousness(
    bounds: Bounds, year: int = IMPERVIOUSNESS_DEFAULT_YEAR, use_cache: bool = True
) -> Raster:
    """Download HRL Imperviousness Density for ``bounds`` as a :class:`Raster`.

    Parameters
    ----------
    bounds : Bounds
        Area of interest in EPSG:3006.
    year : int, default 2024
        Reference year: 2018, 2021 or 2024.
    use_cache : bool, default True
        Reuse a cached download for the same bounds and year.

    Returns
    -------
    Raster
        Sealed share of each 10 m cell in percent (0-100) on a grid anchored at
        the top-left corner of ``bounds`` in EPSG:3006. Cells without data are
        NaN.
    """
    if not HAS_RASTERIO:
        raise ImportError("rasterio is required to download imperviousness data.")
    if year not in IMPERVIOUSNESS_SUPPORTED_YEARS:
        raise ValueError(
            f"Unsupported imperviousness year {year}. "
            f"Supported years: {IMPERVIOUSNESS_SUPPORTED_YEARS}."
        )
    bounds = _ensure_bounds(bounds)
    transform, width, height = _target_grid(bounds)

    cache_path = _cache_path(bounds, year)
    if use_cache and cache_path.exists():
        info(f"Using cached imperviousness {year} data from {cache_path}")
        with rasterio.open(cache_path) as src:
            data = src.read(1)
    else:
        info(f"Downloading HRL Imperviousness Density {year}")
        if year == 2018:
            data = _download_imageserver(bounds, transform, width, height)
        else:
            data = _download_geoserver(year, transform, width, height)
        _write_cache(cache_path, data, transform)

    if np.isnan(data).all():
        warning(
            f"No imperviousness {year} data in bounds {bounds.tuple}. The product "
            "covers the EEA38 countries and the United Kingdom."
        )
    return Raster(data=data, georef=transform, nodata=np.nan, crs="EPSG:3006")


def _download_imageserver(bounds: Bounds, transform: Affine, width: int, height: int):
    data = np.full((height, width), np.nan, dtype=np.float32)
    for row, col, rows, cols in _tiles(width, height, _IMAGESERVER_MAX_PIXELS):
        xmin, ymax = transform * (col, row)
        xmax, ymin = transform * (col + cols, row + rows)
        params = {
            "bbox": f"{xmin},{ymin},{xmax},{ymax}",
            "bboxSR": 3006,
            "imageSR": 3006,
            "size": f"{cols},{rows}",
            "format": "tiff",
            "pixelType": "U8",
            "interpolation": "RSP_NearestNeighbor",
            "f": "image",
        }
        content = _get(IMPERVIOUSNESS_IMAGESERVER_URL, params, 2018)
        with rasterio.open(io.BytesIO(content)) as src:
            values = src.read(1)
        data[row : row + rows, col : col + cols] = _density(values)
    return data


def _download_geoserver(year: int, transform: Affine, width: int, height: int):
    workspace, layer = IMPERVIOUSNESS_GEOSERVER_LAYERS[year]
    url = IMPERVIOUSNESS_GEOSERVER_URL.format(workspace=workspace)
    cell = IMPERVIOUSNESS_CELL_SIZE

    # Request the EPSG:3035 extent of the target grid, snapped to the product's
    # native 10 m cells, then reproject onto the EPSG:3006 grid.
    xmin, ymax = transform * (0, 0)
    xmax, ymin = transform * (width, height)
    to_3035 = Transformer.from_crs("EPSG:3006", "EPSG:3035", always_xy=True)
    xs, ys = to_3035.transform(
        [xmin, xmin, xmax, xmax], [ymin, ymax, ymin, ymax]
    )
    sxmin = math.floor((min(xs) - _REPROJECT_MARGIN_M) / cell) * cell
    symin = math.floor((min(ys) - _REPROJECT_MARGIN_M) / cell) * cell
    sxmax = math.ceil((max(xs) + _REPROJECT_MARGIN_M) / cell) * cell
    symax = math.ceil((max(ys) + _REPROJECT_MARGIN_M) / cell) * cell
    src_width = round((sxmax - sxmin) / cell)
    src_height = round((symax - symin) / cell)
    src_transform = Affine(cell, 0, sxmin, 0, -cell, symax)

    source = np.full((src_height, src_width), np.nan, dtype=np.float32)
    for row, col, rows, cols in _tiles(src_width, src_height, _GEOSERVER_MAX_PIXELS):
        txmin, tymax = src_transform * (col, row)
        txmax, tymin = src_transform * (col + cols, row + rows)
        params = {
            "SERVICE": "WMS",
            "VERSION": "1.3.0",
            "REQUEST": "GetMap",
            "LAYERS": layer,
            "STYLES": "",
            "CRS": "EPSG:3035",
            # WMS 1.3.0 uses northing, easting axis order for EPSG:3035.
            "BBOX": f"{tymin},{txmin},{tymax},{txmax}",
            "WIDTH": cols,
            "HEIGHT": rows,
            "FORMAT": "image/geotiff",
        }
        content = _get(url, params, year)
        source[row : row + rows, col : col + cols] = _read_palette_density(content, year)

    data = np.full((height, width), np.nan, dtype=np.float32)
    reproject(
        source,
        data,
        src_transform=src_transform,
        src_crs="EPSG:3035",
        src_nodata=np.nan,
        dst_transform=transform,
        dst_crs="EPSG:3006",
        dst_nodata=np.nan,
        resampling=Resampling.nearest,
    )
    return data


def _read_palette_density(content: bytes, year: int):
    """Return densities from a WMS palette GeoTIFF, checking its layout."""
    with rasterio.open(io.BytesIO(content)) as src:
        if src.count != 1 or src.dtypes[0] != "uint8":
            raise RuntimeError(
                f"Unexpected imperviousness {year} image layout ({src.count} bands, "
                f"{src.dtypes[0]}). The map service style has changed, so pixel "
                "values no longer give the density."
            )
        values = src.read(1)
    return _density(values)


def _density(values):
    data = values.astype(np.float32)
    data[values > _MAX_DENSITY] = np.nan
    return data


def _get(url: str, params: dict, year: int) -> bytes:
    try:
        response = requests.get(
            url,
            params=params,
            timeout=(IMPERVIOUSNESS_CONNECT_TIMEOUT_SECONDS, IMPERVIOUSNESS_READ_TIMEOUT_SECONDS),
        )
        response.raise_for_status()
    except requests.exceptions.RequestException as exc:
        raise RuntimeError(f"Failed to download imperviousness {year} data: {exc}") from exc
    content_type = response.headers.get("content-type", "")
    if not content_type.startswith("image/"):
        raise RuntimeError(
            f"Imperviousness {year} service returned {content_type or 'no content type'} "
            f"instead of an image: {response.text[:200]}"
        )
    return response.content


def _tiles(width: int, height: int, max_pixels: tuple[int, int]):
    max_cols, max_rows = max_pixels
    for row in range(0, height, max_rows):
        for col in range(0, width, max_cols):
            yield row, col, min(max_rows, height - row), min(max_cols, width - col)


def _target_grid(bounds: Bounds):
    cell = IMPERVIOUSNESS_CELL_SIZE
    width = max(1, math.ceil(round((bounds.xmax - bounds.xmin) / cell, 6)))
    height = max(1, math.ceil(round((bounds.ymax - bounds.ymin) / cell, 6)))
    transform = Affine(cell, 0, bounds.xmin, 0, -cell, bounds.ymax)
    return transform, width, height


def _write_cache(path: Path, data, transform: Affine):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_name(path.name + ".part.tif")
    with rasterio.open(
        tmp_path,
        "w",
        driver="GTiff",
        width=data.shape[1],
        height=data.shape[0],
        count=1,
        dtype="float32",
        crs="EPSG:3006",
        transform=transform,
        nodata=np.nan,
        compress="deflate",
    ) as dst:
        dst.write(data, 1)
    os.replace(tmp_path, path)


def _ensure_bounds(bounds) -> Bounds:
    if isinstance(bounds, Bounds):
        return bounds
    return Bounds(*bounds)


def _cache_path(bounds: Bounds, year: int) -> Path:
    xmin, ymin, xmax, ymax = (round(value) for value in bounds.tuple)
    return (
        Path(cache_dir)
        / "downloaded_imperviousness"
        / f"imperviousness_{year}_{xmin}_{ymin}_{xmax}_{ymax}.tif"
    )


__all__ = [
    "IMPERVIOUSNESS_ATTRIBUTION",
    "IMPERVIOUSNESS_DEFAULT_YEAR",
    "IMPERVIOUSNESS_DOI",
    "IMPERVIOUSNESS_SUPPORTED_YEARS",
    "download_imperviousness",
]
