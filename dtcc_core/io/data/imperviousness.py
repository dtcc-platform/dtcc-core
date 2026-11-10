"""Download Copernicus HRL Imperviousness Density rasters from anonymous services."""

from __future__ import annotations

import hashlib
import io
import json
import math
import os
import uuid
from pathlib import Path

import numpy as np
import requests
from affine import Affine
from pyproj import Transformer
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from ...common.progress import report_progress
from ...model import Bounds, Raster
from .cache import cache_dir
from .logging import info, warning

try:
    import rasterio
    import rasterio.errors
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
# Largest request in 10 m cells: 50 km x 50 km, about 100 MB as float32.
IMPERVIOUSNESS_MAX_CELLS = 25_000_000
IMPERVIOUSNESS_CACHE_VERSION = 2
# 2018 is served as raw values by the EEA ArcGIS ImageServer, in any CRS.
IMPERVIOUSNESS_IMAGESERVER_URL = (
    "https://image.discomap.eea.europa.eu/arcgis/rest/services/GioLandPublic/"
    "HRL_ImperviousnessDensity_2018/ImageServer/exportImage"
)
# 2021 and 2024 are only served as WMS map images in EPSG:3035. The default
# style is a single-band palette whose pixel index is the density in percent.
# The service lists this style as CLMS_NVLCC_IMD but rejects that name in
# GetMap, so the request uses the default style and the colormap is checked.
IMPERVIOUSNESS_GEOSERVER_URL = (
    "https://geoserver.geoville.com/geoserver/{workspace}/wms"
)
IMPERVIOUSNESS_GEOSERVER_LAYERS = {
    2021: ("nvlcc_2021", "HRL_NVLCC_IMD_10m"),
    2024: ("nvlcc_2024", "HRL_NVLCC_IMD_2024_10m"),
}
# Palette colours (index: RGB) of the default style. A different colour at
# any of these indices means the style no longer maps index to density.
IMPERVIOUSNESS_PALETTE_CHECKS = {
    0: (255, 255, 255),
    50: (175, 74, 51),
    100: (113, 12, 2),
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
# Valid densities are 0-100. 2018 marks sea and areas outside the product as
# 255 (requested with noData=255); 2021 and 2024 mark sea as 101.
_MAX_DENSITY = 100
_NODATA_2018 = 255
# Margin in metres around the EPSG:3035 request so reprojection covers the
# whole EPSG:3006 grid.
_REPROJECT_MARGIN_M = 50.0
# A response may hold at most this many bytes per requested pixel, plus
# headroom for the GeoTIFF header.
_MAX_BYTES_PER_PIXEL = 4
_RESPONSE_HEADROOM_BYTES = 1 << 20
# Allowed difference between the returned and requested grid origin, in metres.
_ORIGIN_TOLERANCE_M = 0.01


def download_imperviousness(
    bounds: Bounds, year: int = IMPERVIOUSNESS_DEFAULT_YEAR, use_cache: bool = True
) -> Raster:
    """Download HRL Imperviousness Density for ``bounds`` as a :class:`Raster`.

    Parameters
    ----------
    bounds : Bounds
        Area of interest in EPSG:3006, at most ``IMPERVIOUSNESS_MAX_CELLS``
        10 m cells.
    year : int, default 2024
        Reference year: 2018, 2021 or 2024.
    use_cache : bool, default True
        Reuse a cached download for the same exact bounds and year.

    Returns
    -------
    Raster
        Sealed share of each 10 m cell in percent (0-100) on a grid anchored at
        the top-left corner of ``bounds`` in EPSG:3006. Sea is NaN. Outside
        the product, 2018 is NaN, but 2021 and 2024 read as 0 because their
        map images draw it the same as unsealed ground.

    Raises
    ------
    DatasetUpstreamError
        If a service cannot be reached, returns an error, or returns an image
        that does not match the request or the expected palette. Its
        ``failure_class`` tells transient failures (connection, timeout,
        HTTP 5xx) from others. Nothing is cached on failure.
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
    if width * height > IMPERVIOUSNESS_MAX_CELLS:
        raise ValueError(
            f"Imperviousness bounds {bounds.tuple} cover {width * height} cells of "
            f"10 m; the limit is {IMPERVIOUSNESS_MAX_CELLS} (50 km x 50 km)."
        )

    cache_path = _cache_path(bounds, year)
    if use_cache and cache_path.exists():
        info(f"Using cached imperviousness {year} data from {cache_path}")
        with rasterio.open(cache_path) as src:
            data = src.read(1)
    else:
        info(f"Downloading HRL Imperviousness Density {year}")
        with _session() as session:
            if year == 2018:
                data = _download_imageserver(session, transform, width, height)
            else:
                data = _download_geoserver(session, year, transform, width, height)
        _write_cache(cache_path, data, transform)

    if np.isnan(data).all():
        warning(
            f"No imperviousness {year} data in bounds {bounds.tuple}. The product "
            "covers the EEA38 countries and the United Kingdom."
        )
    elif year != 2018 and np.nanmax(data) == 0:
        warning(
            f"Every imperviousness {year} cell in bounds {bounds.tuple} reads 0 %. "
            f"Outside the product, {year} also reads 0 %; the product covers the "
            "EEA38 countries and the United Kingdom."
        )
    return Raster(data=data, georef=transform, nodata=np.nan, crs="EPSG:3006")


def empty_imperviousness_raster(bounds: Bounds) -> Raster:
    """Return an all-NaN imperviousness Raster on the 10 m grid of ``bounds``."""
    transform, width, height = _target_grid(_ensure_bounds(bounds))
    data = np.full((height, width), np.nan, dtype=np.float32)
    return Raster(data=data, georef=transform, nodata=np.nan, crs="EPSG:3006")


def _download_imageserver(session, transform: Affine, width: int, height: int):
    data = np.full((height, width), np.nan, dtype=np.float32)
    tiles = list(_tiles(width, height, _IMAGESERVER_MAX_PIXELS))
    for index, (row, col, rows, cols) in enumerate(tiles):
        xmin, ymax = transform * (col, row)
        xmax, ymin = transform * (col + cols, row + rows)
        params = {
            "bbox": f"{xmin},{ymin},{xmax},{ymax}",
            "bboxSR": 3006,
            "imageSR": 3006,
            "size": f"{cols},{rows}",
            "format": "tiff",
            "pixelType": "U8",
            "noData": _NODATA_2018,
            "interpolation": "RSP_NearestNeighbor",
            "f": "image",
        }
        url = IMPERVIOUSNESS_IMAGESERVER_URL
        content = _get(session, url, params, 2018, cols * rows)
        expected = Affine(transform.a, 0, xmin, 0, transform.e, ymax)
        values = _read_image(content, url, 2018, cols, rows, expected, 3006)
        data[row : row + rows, col : col + cols] = _density(values)
        _report_tile(index, len(tiles), 2018)
    return data


def _download_geoserver(session, year: int, transform: Affine, width: int, height: int):
    workspace, layer = IMPERVIOUSNESS_GEOSERVER_LAYERS[year]
    url = IMPERVIOUSNESS_GEOSERVER_URL.format(workspace=workspace)
    cell = IMPERVIOUSNESS_CELL_SIZE

    # Request the EPSG:3035 extent of the target grid, snapped to the product's
    # native 10 m cells, then reproject onto the EPSG:3006 grid.
    xmin, ymax = transform * (0, 0)
    xmax, ymin = transform * (width, height)
    to_3035 = Transformer.from_crs("EPSG:3006", "EPSG:3035", always_xy=True)
    xs, ys = to_3035.transform([xmin, xmin, xmax, xmax], [ymin, ymax, ymin, ymax])
    sxmin = math.floor((min(xs) - _REPROJECT_MARGIN_M) / cell) * cell
    symin = math.floor((min(ys) - _REPROJECT_MARGIN_M) / cell) * cell
    sxmax = math.ceil((max(xs) + _REPROJECT_MARGIN_M) / cell) * cell
    symax = math.ceil((max(ys) + _REPROJECT_MARGIN_M) / cell) * cell
    src_width = round((sxmax - sxmin) / cell)
    src_height = round((symax - symin) / cell)
    src_transform = Affine(cell, 0, sxmin, 0, -cell, symax)

    source = np.full((src_height, src_width), np.nan, dtype=np.float32)
    tiles = list(_tiles(src_width, src_height, _GEOSERVER_MAX_PIXELS))
    for index, (row, col, rows, cols) in enumerate(tiles):
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
        content = _get(session, url, params, year, cols * rows)
        expected = Affine(cell, 0, txmin, 0, -cell, tymax)
        values = _read_image(
            content, url, year, cols, rows, expected, 3035, palette=True
        )
        source[row : row + rows, col : col + cols] = _density(values)
        _report_tile(index, len(tiles), year)

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


def _read_image(
    content: bytes,
    url: str,
    year: int,
    cols: int,
    rows: int,
    expected,
    epsg: int,
    palette: bool = False,
):
    """Check a response image against the request, then return its band."""
    try:
        with rasterio.open(io.BytesIO(content)) as src:
            _check_image(src, url, year, cols, rows, expected, epsg)
            if palette:
                _check_palette(src, url, year)
            return src.read(1)
    except rasterio.errors.RasterioError as exc:
        raise _upstream_error(
            year, url, "invalid_payload", f"unreadable image: {exc}"
        ) from exc


def _check_image(src, url: str, year: int, cols: int, rows: int, expected, epsg: int):
    """Raise unless a response image has the requested size, layout and grid."""
    problems = []
    if (src.width, src.height) != (cols, rows):
        problems.append(f"size {src.width}x{src.height}, requested {cols}x{rows}")
    if src.count != 1 or src.dtypes[0] != "uint8":
        problems.append(
            f"{src.count} bands of {src.dtypes[0]}, expected 1 band of uint8"
        )
    # The 2018 service labels EPSG:3006 with a SWEREF99 datum code that GDAL
    # only matches at reduced confidence.
    if src.crs is None or src.crs.to_epsg(confidence_threshold=20) != epsg:
        problems.append(f"CRS {src.crs}, expected EPSG:{epsg}")
    if not src.transform.almost_equals(expected, precision=_ORIGIN_TOLERANCE_M):
        problems.append(
            f"grid {tuple(src.transform)[:6]}, expected {tuple(expected)[:6]}"
        )
    if problems:
        raise _upstream_error(
            year, url, "invalid_payload", "unexpected image: " + "; ".join(problems)
        )


def _check_palette(src, url: str, year: int):
    """Raise unless the WMS image palette still maps pixel index to density."""
    try:
        colormap = src.colormap(1)
    except ValueError:
        colormap = {}
    changed = {
        index: colormap.get(index)
        for index, rgb in IMPERVIOUSNESS_PALETTE_CHECKS.items()
        if tuple(colormap.get(index, ()))[:3] != rgb
    }
    if changed:
        raise _upstream_error(
            year,
            url,
            "invalid_payload",
            f"the map style has changed (palette entries {changed}), so pixel "
            "values no longer give the density",
        )


def _density(values):
    data = values.astype(np.float32)
    data[values > _MAX_DENSITY] = np.nan
    return data


def _session() -> requests.Session:
    """Return a session that retries rate limits and transient server errors."""
    session = requests.Session()
    retry = Retry(
        total=3,
        connect=2,
        read=0,
        status=3,
        backoff_factor=1.0,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=["GET"],
        raise_on_status=False,
        # A long Retry-After would block the download for that long.
        respect_retry_after_header=False,
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("http://", adapter)
    session.mount("https://", adapter)
    return session


def _get(session, url: str, params: dict, year: int, pixels: int) -> bytes:
    """Download one image, raising DatasetUpstreamError for any failure."""
    from ...datasets.dataset import DatasetDescriptor

    max_bytes = pixels * _MAX_BYTES_PER_PIXEL + _RESPONSE_HEADROOM_BYTES
    try:
        response = session.get(
            url,
            params=params,
            timeout=(
                IMPERVIOUSNESS_CONNECT_TIMEOUT_SECONDS,
                IMPERVIOUSNESS_READ_TIMEOUT_SECONDS,
            ),
            stream=True,
        )
        with response:
            response.raise_for_status()
            content = bytearray()
            for chunk in response.iter_content(chunk_size=1 << 16):
                content.extend(chunk)
                if len(content) > max_bytes:
                    raise _upstream_error(
                        year,
                        url,
                        "invalid_payload",
                        f"the response exceeds {max_bytes} bytes for {pixels} pixels",
                    )
    except requests.exceptions.RequestException as exc:
        raise DatasetDescriptor.build_upstream_error(
            "imperviousness", "download_image", url, exc
        ) from exc
    content = bytes(content)

    # ArcGIS answers errors with HTTP 200, labelled image/tiff, and a JSON body.
    if content[:64].lstrip()[:1] == b"{":
        raise _json_error(content, year, url)
    content_type = response.headers.get("content-type", "")
    if not content_type.startswith("image/"):
        text = content[:300].decode("utf-8", errors="replace")
        raise _upstream_error(
            year,
            url,
            "invalid_payload",
            f"the service returned {content_type or 'no content type'} instead of "
            f"an image: {text}",
        )
    return content


def _json_error(content: bytes, year: int, url: str):
    try:
        error = json.loads(content).get("error") or {}
    except (ValueError, AttributeError):
        error = {}
    if not isinstance(error, dict):
        error = {}
    code = error.get("code")
    code = code if isinstance(code, int) else None
    parts = [str(error.get("message") or "unreadable JSON instead of an image")]
    details = error.get("details")
    if isinstance(details, list):
        parts += [str(detail) for detail in details]
    elif details:
        parts.append(str(details))
    if code is not None and 500 <= code <= 599:
        failure_class = "http_5xx"
    elif code is not None and 400 <= code <= 499:
        failure_class = "http_4xx"
    else:
        failure_class = "invalid_payload"
    return _upstream_error(year, url, failure_class, " ".join(parts), status_code=code)


def _upstream_error(
    year: int, url: str, failure_class: str, message: str, status_code=None
):
    from ...datasets.dataset import DatasetUpstreamError

    return DatasetUpstreamError(
        dataset="imperviousness",
        operation="download_image",
        target=url,
        failure_class=failure_class,
        status_code=status_code,
        message=f"Imperviousness {year} download failed: {message}",
    )


def _report_tile(index: int, total: int, year: int):
    report_progress(
        current=index + 1,
        total=total,
        message=f"Downloading imperviousness {year} tiles ({index + 1}/{total})",
    )


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
    # A unique temporary name per writer, so parallel downloads of the same
    # area cannot interleave their writes before the rename.
    tmp_path = path.with_name(f"{path.name}.{uuid.uuid4().hex}.part.tif")
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
    exact = ",".join(repr(float(value)) for value in bounds.tuple)
    digest = hashlib.sha256(exact.encode("ascii")).hexdigest()[:16]
    return (
        Path(cache_dir)
        / "downloaded_imperviousness"
        / f"imperviousness_v{IMPERVIOUSNESS_CACHE_VERSION}_{year}_{digest}.tif"
    )


__all__ = [
    "IMPERVIOUSNESS_DEFAULT_YEAR",
    "IMPERVIOUSNESS_DOI",
    "IMPERVIOUSNESS_MAX_CELLS",
    "IMPERVIOUSNESS_SUPPORTED_YEARS",
    "download_imperviousness",
    "empty_imperviousness_raster",
]
