#!/usr/bin/env python3

"""
Access to Digital Earth Sweden, the Sentinel-2 satellite imagery platform
operated by RISE together with the Swedish National Space Agency.

Two public services are used, both of which are anonymous:

  * an OGC WMS server (datacube-ows) that renders the imagery, and
  * a STAC catalogue that describes which acquisitions exist and how
    cloudy each one is.

The catalogue is used to pick a date and the WMS server to fetch pixels, so
by default the user gets the most recent *clear* image rather than simply the
most recent one, which over Sweden is very often solid cloud.
"""

import math
import os
import warnings
from datetime import UTC, datetime, timedelta
from datetime import date as Date

import numpy as np
import pyproj
import rasterio
import requests
from rasterio.errors import NotGeoreferencedWarning
from rasterio.io import MemoryFile
from rasterio.transform import from_origin
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from dtcc_core.common.progress import report_progress
from dtcc_core.model import Raster

from . import cache
from .logging import debug, info, warning

# ------------------------------------------------------------------------
# 1) Global constants/paths
# ------------------------------------------------------------------------
CACHE_DIR = os.path.join(cache.cache_dir, "downloaded_des")

# Part of every cache file name. Raise it when the stored imagery changes, so
# files written under an earlier version are downloaded again.
CACHE_VERSION = 2

WMS_URL = "https://ows.digitalearth.se"
STAC_URL = "https://explorer.digitalearth.se/stac"
STAC_COLLECTION = "s2_msi_l2a"

# The server advertises MaxWidth/MaxHeight of 512 in its WMS capabilities,
# so anything larger has to be requested as a grid of tiles and stitched.
MAX_TILE_PIXELS = 512

# Native ground resolution of the Sentinel-2 visual bands, in metres.
DEFAULT_RESOLUTION = 10.0

# Largest mosaic, in pixels, that one call may request. At the default
# resolution this is a 100 km square: about 400 tiles and 400 MB of memory.
MAX_PIXELS = 10_000 * 10_000

# Attempts per request after a dropped connection, a timeout or a 429/5xx.
HTTP_RETRIES = 3

USER_AGENT = "dtcc-core (https://github.com/dtcc-platform/dtcc-core)"

# Acquisitions cloudier than this are skipped when picking a date.
DEFAULT_MAX_CLOUD = 10.0

# How far back to look for a clear acquisition, in days.
DEFAULT_SEARCH_DAYS = 120

# Size of the throwaway request used to check that a date really has pixels.
PROBE_PIXELS = 32

# How many catalogue dates to probe before giving up.
MAX_DATE_CANDIDATES = 12

# Warn when less than this fraction of the requested area has imagery.
COVERAGE_WARNING_THRESHOLD = 0.99

# WMS layer serving the whole Sentinel-2 archive.
LAYER = "sentinel2_reflectance"

VALID_STYLES = [
    "rgb",
    "false_color",
    "false_color_urban",
    "swir",
    "ndvi",
    "ndmi",
    "ndwi",
    "nbi",
    "msi",
]

# Styles that render every pixel the acquisition covers, so a transparent
# pixel means no data. The index styles also hide pixels by content, for
# example ndwi hides land and ndvi hides water.
OPAQUE_STYLES = ("rgb", "false_color", "false_color_urban", "swir")


def create_session(retries=HTTP_RETRIES, backoff_factor=1.0):
    """
    Create a requests session for the Digital Earth Sweden services.

    Dropped connections and timeouts are retried as well as 429 and 5xx
    responses, so one network fault does not discard a long tiled download.
    Both services are read-only, so the POSTed catalogue searches are safe
    to repeat.
    """
    retry = Retry(
        total=retries,
        connect=retries,
        read=retries,
        status=retries,
        backoff_factor=backoff_factor,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=["GET", "POST"],
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry)
    session = requests.Session()
    session.headers["User-Agent"] = USER_AGENT
    session.mount("http://", adapter)
    session.mount("https://", adapter)
    return session


class ServiceError(requests.RequestException):
    """The imagery server answered with an error report instead of an image."""


# ------------------------------------------------------------------------
# 2) Catalogue queries (STAC)
# ------------------------------------------------------------------------
def _to_wgs84_bbox(bbox):
    """
    Convert a bounding box from EPSG:3006 to EPSG:4326.

    The STAC catalogue expects longitude/latitude, while the rest of
    this package works in SWEREF 99 TM. All four corners are transformed rather
    than just two, since the reprojected box is not axis aligned.

    'bbox' is a tuple (xmin, ymin, xmax, ymax) in EPSG:3006.
    Returns (lon_min, lat_min, lon_max, lat_max).
    """
    transformer = pyproj.Transformer.from_crs("EPSG:3006", "EPSG:4326", always_xy=True)
    xmin, ymin, xmax, ymax = bbox
    corners = ((xmin, ymin), (xmin, ymax), (xmax, ymin), (xmax, ymax))
    lons, lats = zip(*(transformer.transform(x, y) for x, y in corners), strict=True)
    return (min(lons), min(lats), max(lons), max(lats))


def _acquisition_datetime(feature):
    """
    Return the acquisition timestamp of a STAC feature as a string.

    Some items carry the timestamp in 'datetime' and others, notably the
    mosaicked ones, leave that null and use 'start_datetime' instead.
    """
    props = feature.get("properties", {})
    return props.get("datetime") or props.get("start_datetime")


def _stac_search_all(session, payload, max_pages=20):
    """
    Run a STAC search and follow its 'next' links until the results run out.

    The catalogue returns results oldest first and its sort extension cannot
    order by time here, because every item leaves 'datetime' null and carries
    the timestamp in 'start_datetime' instead. So the only way to be sure the
    newest acquisitions are seen is to read every page.
    """
    features = []
    url = f"{STAC_URL}/search"
    body = dict(payload)
    truncated = False
    for page_number in range(max_pages):
        debug(f"[POST] to {url} with payload={body}")
        resp = session.post(url, json=body, timeout=60)
        if resp.status_code != 200:
            raise requests.HTTPError(
                f"STAC search failed with status {resp.status_code}:\n{resp.text}",
                request=resp.request,
                response=resp,
            )
        data = resp.json()
        page = data.get("features", [])
        features.extend(page)
        next_link = next(
            (link for link in data.get("links", []) if link.get("rel") == "next"), None
        )
        if not page or not next_link:
            break
        if page_number == max_pages - 1:
            truncated = True
            break
        url = next_link.get("href", url)
        body = {**payload, **(next_link.get("body") or {})}

    if truncated:
        # Results arrive oldest first, so stopping early drops the *newest*
        # acquisitions, which are exactly the ones a caller usually wants.
        warning(
            f"Stopped reading the catalogue after {max_pages} pages, so the most "
            f"recent acquisitions were not seen. Narrow the search with a shorter "
            f"search_days or a smaller area, or raise max_pages."
        )
    return features


def list_available_dates(
    bbox,
    max_cloud=DEFAULT_MAX_CLOUD,
    start=None,
    end=None,
    limit=500,
    session=None,
):
    """
    List the Sentinel-2 acquisition dates that cover 'bbox'.

    'bbox' is a tuple (xmin, ymin, xmax, ymax) in EPSG:3006. 'start' and
    'end' bound the search, as datetime.date objects or 'YYYY-MM-DD'
    strings. 'start' defaults to DEFAULT_SEARCH_DAYS days ago, because the
    whole archive takes minutes to read; pass an early date such as
    '2018-01-01' to search all of it. Items with a cloud cover above
    'max_cloud' percent are discarded.

    Returns a list of (date, cloud_cover) tuples sorted newest first, where
    'date' is a 'YYYY-MM-DD' string. A date appearing on several satellite
    tiles is reported once, with the cloudiest of those tiles, so that a
    date is only called clear when the whole area is clear.
    """
    if start is None:
        start = datetime.now(UTC).date() - timedelta(days=DEFAULT_SEARCH_DAYS)
    start = _normalize_date(start)
    end = _normalize_date(end)

    session = session or create_session()
    payload = {
        "collections": [STAC_COLLECTION],
        "bbox": list(_to_wgs84_bbox(bbox)),
        "limit": limit,
        "datetime": f"{start}/{end or '..'}",
    }

    worst_cloud_per_date = {}
    for feature in _stac_search_all(session, payload):
        timestamp = _acquisition_datetime(feature)
        if not timestamp:
            continue
        date = timestamp[:10]
        cloud = feature.get("properties", {}).get("eo:cloud_cover")
        if cloud is None:
            continue
        # Keep the worst tile for the date; a date is only as clear as its
        # cloudiest part, otherwise half the requested area comes back white.
        worst_cloud_per_date[date] = max(worst_cloud_per_date.get(date, 0.0), cloud)

    dates = [
        (date, cloud)
        for date, cloud in worst_cloud_per_date.items()
        if cloud <= max_cloud
    ]
    dates.sort(reverse=True)
    debug(f"Found {len(dates)} dates below {max_cloud}% cloud cover")
    return dates


def latest_clear_date(
    bbox,
    max_cloud=DEFAULT_MAX_CLOUD,
    search_days=DEFAULT_SEARCH_DAYS,
    session=None,
    max_candidates=MAX_DATE_CANDIDATES,
):
    """
    Return the most recent date on which 'bbox' was imaged without much cloud.

    Only the last 'search_days' days are considered, to keep the catalogue
    query small. Each candidate is then probed against the imagery server,
    because the catalogue and the server are not perfectly in step: the
    catalogue lists acquisitions that the server has no pixels for, and
    trusting it blindly yields a blank image. The newest candidate that covers
    the whole area is returned; if none does, the one covering most of it.

    Raises RuntimeError if nothing clear enough was found, since silently
    returning a cloudy date would hand the caller a white image.
    """
    return _pick_clear_date(bbox, max_cloud, search_days, session, max_candidates)[0]


def _pick_clear_date(bbox, max_cloud, search_days, session, max_candidates):
    """Do the work of latest_clear_date and return (date, rgb coverage)."""
    session = session or create_session()
    start = (datetime.now(UTC) - timedelta(days=search_days)).strftime("%Y-%m-%d")
    candidates = list_available_dates(
        bbox, max_cloud=max_cloud, start=start, session=session
    )
    if not candidates:
        raise RuntimeError(
            f"No Sentinel-2 acquisition below {max_cloud}% cloud cover in the last "
            f"{search_days} days for this area. Pass an explicit date, raise "
            f"max_cloud, or widen search_days."
        )

    # A date at the edge of a satellite pass covers only part of the area, so
    # an older date covering all of it is preferred.
    probed = candidates[:max_candidates]
    best = None
    for date, cloud in probed:
        coverage = _imagery_coverage(session, bbox, date)
        if coverage >= COVERAGE_WARNING_THRESHOLD:
            info(f"Using acquisition from {date} ({cloud:.1f}% cloud cover)")
            return date, coverage
        debug(f"{date} has imagery for {coverage:.0%} of the area")
        if coverage > 0.0 and (best is None or coverage > best[2]):
            best = (date, cloud, coverage)

    if best is not None:
        date, cloud, coverage = best
        info(
            f"Using acquisition from {date} ({cloud:.1f}% cloud cover), which "
            f"covers {coverage:.0%} of the area; no candidate covered all of it"
        )
        return date, coverage

    raise RuntimeError(
        f"The catalogue lists {len(candidates)} acquisition(s) below {max_cloud}% "
        f"cloud cover for this area, but none of the {len(probed)} most recent "
        f"had imagery on the server. Pass an explicit date or widen search_days."
    )


# ------------------------------------------------------------------------
# 3) Imagery requests (WMS)
# ------------------------------------------------------------------------
def _normalize_date(value):
    """
    Return 'value' as a 'YYYY-MM-DD' string, or None if it is None.

    Accepts a datetime.date, a datetime.datetime or an ISO date string, and
    raises ValueError for anything else, including impossible dates such as
    '2026-02-30'.
    """
    if value is None:
        return None
    if isinstance(value, datetime):
        value = value.date()
    if isinstance(value, Date):
        return value.isoformat()
    try:
        return Date.fromisoformat(value).isoformat()
    except (TypeError, ValueError):
        raise ValueError(
            f"Invalid date {value!r}; expected a datetime.date or 'YYYY-MM-DD'."
        ) from None


def _mosaic_size(bbox, resolution):
    """
    Return the (width, height) in pixels of a mosaic covering 'bbox'.

    Raises ValueError if 'bbox' is narrower than one pixel in either
    direction, or the mosaic would exceed MAX_PIXELS.
    """
    if resolution <= 0:
        raise ValueError(f"Resolution must be positive, got {resolution}.")
    xmin, ymin, xmax, ymax = bbox
    if xmax - xmin < resolution or ymax - ymin < resolution:
        raise ValueError(
            f"Bounds must be at least one pixel ({resolution:g} m) wide and tall; "
            f"got {xmax - xmin:g} m x {ymax - ymin:g} m."
        )
    width = math.ceil((xmax - xmin) / resolution)
    height = math.ceil((ymax - ymin) / resolution)
    if width * height > MAX_PIXELS:
        raise ValueError(
            f"Request is {width}x{height} px, more than the limit of {MAX_PIXELS:,}. "
            "Use a coarser resolution or a smaller area."
        )
    return width, height


def _tile_grid(bbox, resolution):
    """
    Split 'bbox' into tiles no larger than MAX_TILE_PIXELS on a side.

    Returns (width, height, tiles) where width and height are the pixel
    dimensions of the whole mosaic and each tile is a dict with its pixel
    offset, pixel size and its own bounding box in EPSG:3006.
    """
    xmin, ymin, xmax, ymax = bbox
    width, height = _mosaic_size(bbox, resolution)

    tiles = []
    for row_off in range(0, height, MAX_TILE_PIXELS):
        tile_height = min(MAX_TILE_PIXELS, height - row_off)
        for col_off in range(0, width, MAX_TILE_PIXELS):
            tile_width = min(MAX_TILE_PIXELS, width - col_off)
            # Rows count downwards from the top of the mosaic, which is ymax.
            tiles.append(
                {
                    "col_off": col_off,
                    "row_off": row_off,
                    "width": tile_width,
                    "height": tile_height,
                    "bbox": (
                        xmin + col_off * resolution,
                        ymax - (row_off + tile_height) * resolution,
                        xmin + (col_off + tile_width) * resolution,
                        ymax - row_off * resolution,
                    ),
                }
            )
    debug(f"Requesting {width}x{height} px as {len(tiles)} tile(s)")
    return width, height, tiles


def get_map(session, layer, style, bbox, width, height, date):
    """
    Fetch a single rendered tile from the WMS server and return the raw PNG.

    'bbox' is a tuple (xmin, ymin, xmax, ymax) in EPSG:3006.
    """
    xmin, ymin, xmax, ymax = bbox
    params = {
        "service": "WMS",
        "version": "1.3.0",
        "request": "GetMap",
        "layers": layer,
        "styles": style,
        "crs": "EPSG:3006",
        # WMS 1.3.0 orders the bounding box the way the CRS declares its axes.
        # EPSG:3006 is northing first, so this is ymin,xmin,ymax,xmax and not
        # the usual xmin,ymin,xmax,ymax. Getting this backwards does not fail,
        # it quietly returns an empty image.
        "bbox": f"{ymin},{xmin},{ymax},{xmax}",
        "width": width,
        "height": height,
        "format": "image/png",
        "time": date,
    }
    debug(f"[GET] {WMS_URL} bbox={params['bbox']} time={date}")
    resp = session.get(WMS_URL, params=params, timeout=120)
    if resp.status_code == 400 and "InvalidDimensionValue" in resp.text:
        raise ValueError(
            f"{date} is outside the date range served by Digital Earth Sweden."
        )
    if resp.status_code != 200:
        raise requests.HTTPError(
            f"WMS request failed with status {resp.status_code}:\n{resp.text}",
            request=resp.request,
            response=resp,
        )
    # Errors come back as an XML ServiceExceptionReport with status 200.
    if "xml" in resp.headers.get("Content-Type", ""):
        raise ServiceError(
            f"WMS returned an error:\n{resp.text}",
            request=resp.request,
            response=resp,
        )
    return resp.content


def _tile_to_array(png_bytes, width, height):
    """
    Decode a PNG tile into an (4, height, width) uint8 RGBA array.

    Tiles with no data for the requested date come back as a small greyscale
    image rather than RGBA, and are treated as fully transparent.
    """
    # Tiles are plain PNGs; their position comes from the request, so the
    # warning that they carry no georeferencing does not apply.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", NotGeoreferencedWarning)
        with MemoryFile(png_bytes) as memfile:
            with memfile.open() as src:
                data = src.read()

    if data.shape[0] < 3:
        return np.zeros((4, height, width), dtype=np.uint8)
    if data.shape[0] == 3:
        alpha = np.full((1, data.shape[1], data.shape[2]), 255, dtype=np.uint8)
        data = np.concatenate([data, alpha], axis=0)
    return data[:4].astype(np.uint8)


def _imagery_coverage(session, bbox, date):
    """
    Return the fraction of 'bbox' that the WMS server has pixels for on 'date'.

    This is a deliberately tiny request, PROBE_PIXELS on a side, because an
    empty answer comes back as a few hundred bytes and a full sized one would
    waste a download only to be thrown away. It always uses the rgb style,
    since the index styles leave pixels transparent where they have data.
    """
    png_bytes = get_map(session, LAYER, "rgb", bbox, PROBE_PIXELS, PROBE_PIXELS, date)
    alpha = _tile_to_array(png_bytes, PROBE_PIXELS, PROBE_PIXELS)[3]
    return float((alpha > 0).mean())


# ------------------------------------------------------------------------
# 4) Cache
# ------------------------------------------------------------------------
def _read_geotiff(path):
    """Load a cached GeoTIFF as an (H, W, 4) RGBA Raster."""
    with rasterio.open(path) as src:
        data = src.read()
        transform = src.transform
    return Raster(data=np.moveaxis(data, 0, -1), georef=transform, crs="EPSG:3006")


def _write_geotiff(path, mosaic, transform):
    """
    Write an RGBA mosaic to 'path' as a GeoTIFF.

    The file is written under a temporary name and renamed when complete, so
    an interrupted write never leaves a truncated file at 'path'.
    """
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp_path = path + ".part"
    try:
        with rasterio.open(
            tmp_path,
            "w",
            driver="GTiff",
            height=mosaic.shape[1],
            width=mosaic.shape[2],
            count=4,
            dtype="uint8",
            crs="EPSG:3006",
            transform=transform,
            photometric="RGB",
            compress="deflate",
        ) as dst:
            dst.write(mosaic)
            dst.colorinterp = [
                rasterio.enums.ColorInterp.red,
                rasterio.enums.ColorInterp.green,
                rasterio.enums.ColorInterp.blue,
                rasterio.enums.ColorInterp.alpha,
            ]
        os.replace(tmp_path, path)
    except BaseException:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        raise


def _download_mosaic(session, bbox, resolution, style, date, tile_errors):
    """
    Download 'bbox' tile by tile and stitch the tiles into one mosaic.

    Returns (mosaic, failed_fraction): a (4, height, width) uint8 RGBA array
    and the fraction of its pixels whose tile failed. Failed tiles are handled
    as described for download_imagery's 'tile_errors'.
    """
    width, height, tiles = _tile_grid(bbox, resolution)
    mosaic = np.zeros((4, height, width), dtype=np.uint8)

    failed_pixels = 0
    for index, tile in enumerate(tiles, start=1):
        message = f"Downloading tile {index}/{len(tiles)} for {date}"
        info(message)
        report_progress(current=index - 1, total=len(tiles), message=message)
        try:
            png_bytes = get_map(
                session, LAYER, style, tile["bbox"], tile["width"], tile["height"], date
            )
        except requests.RequestException as exc:
            if tile_errors is None:
                raise
            tile_errors.append(exc)
            failed_pixels += tile["width"] * tile["height"]
            warning(f"Tile {index}/{len(tiles)} failed and was left transparent: {exc}")
            continue
        data = _tile_to_array(png_bytes, tile["width"], tile["height"])
        rows = slice(tile["row_off"], tile["row_off"] + tile["height"])
        cols = slice(tile["col_off"], tile["col_off"] + tile["width"])
        mosaic[:, rows, cols] = data
    report_progress(current=len(tiles), total=len(tiles), message="Tiles downloaded")

    if failed_pixels == width * height:
        raise tile_errors[-1]
    return mosaic, failed_pixels / (width * height)


# ------------------------------------------------------------------------
# 5) Entry point
# ------------------------------------------------------------------------
def download_imagery(
    bbox,
    date=None,
    style="rgb",
    max_cloud=DEFAULT_MAX_CLOUD,
    resolution=DEFAULT_RESOLUTION,
    search_days=DEFAULT_SEARCH_DAYS,
    session=None,
    tile_errors=None,
):
    """
    Download Sentinel-2 imagery covering 'bbox' as an RGBA Raster.

    'bbox' is a tuple (xmin, ymin, xmax, ymax) in EPSG:3006. 'date' is a
    datetime.date or a 'YYYY-MM-DD' string; if omitted the most recent
    acquisition with less than 'max_cloud' percent cloud is used. 'style'
    selects the rendering, see VALID_STYLES, and 'resolution' is the ground
    sample distance in metres.

    A tile that still fails after retries raises its requests exception. If
    'tile_errors' is a list, the exception is appended to it instead and the
    tile is left transparent, unless every tile failed.

    Returns a Raster with (height, width, 4) uint8 data. Complete results are
    cached as GeoTIFFs, so asking for the same image twice only downloads it
    once. A result that covers only part of the area is returned with a
    warning and not cached.

    Raises RuntimeError if the server has no imagery for the area on 'date',
    and ValueError if 'date' is invalid or outside the served archive.
    """
    if style not in VALID_STYLES:
        raise ValueError(f"Invalid style '{style}'. Must be one of {VALID_STYLES}.")

    xmin, ymin, xmax, ymax = bbox
    if xmax <= xmin or ymax <= ymin:
        raise ValueError(f"Invalid bounds {bbox}; expected (xmin, ymin, xmax, ymax).")
    _mosaic_size(bbox, resolution)
    date = _normalize_date(date)

    session = session or create_session()
    explicit_date = date is not None
    probe_coverage = None
    if not explicit_date:
        report_progress(percent=0, message="Finding a clear acquisition date...")
        date, probe_coverage = _pick_clear_date(
            bbox, max_cloud, search_days, session, MAX_DATE_CANDIDATES
        )

    out_path = os.path.join(
        CACHE_DIR,
        f"des_v{CACHE_VERSION}_{style}_{date}_{xmin:.2f}_{ymin:.2f}_{xmax:.2f}_{ymax:.2f}"
        f"_{resolution:g}m.tif",
    )
    if os.path.exists(out_path):
        info(f"Imagery for {date} already in cache, skipping download.")
        return _read_geotiff(out_path)

    # An explicit date is probed after the cache lookup, so cached results
    # still load offline, and before downloading every tile of a date with no
    # pixels. Automatic dates were probed while being picked.
    if explicit_date:
        probe_coverage = _imagery_coverage(session, bbox, date)
        if probe_coverage == 0.0:
            raise RuntimeError(
                f"Digital Earth Sweden has no imagery for {date} in this area. "
                "Use list_available_dates() to pick a date, or omit date."
            )

    mosaic, failed_fraction = _download_mosaic(
        session, bbox, resolution, style, date, tile_errors
    )

    # For opaque styles the full mosaic also reveals a single tile that came
    # back blank; index styles rely on the rgb probe and the failed tiles.
    if style in OPAQUE_STYLES:
        covered = float((mosaic[3] > 0).mean())
    else:
        covered = probe_coverage * (1.0 - failed_fraction)
    if covered == 0.0:
        raise RuntimeError(
            f"The imagery returned for {date} is empty. The area may lie outside "
            f"the Digital Earth Sweden coverage, or have no acquisition that day."
        )

    # The mosaic is anchored at its top left corner, which is (xmin, ymax).
    transform = from_origin(xmin, ymax, resolution, resolution)
    raster = Raster(data=np.moveaxis(mosaic, 0, -1), georef=transform, crs="EPSG:3006")

    if covered < COVERAGE_WARNING_THRESHOLD:
        # A partly covered mosaic may come from a passing server fault, so it
        # is returned for this call only and the next call downloads again.
        warning(
            f"Only {covered:.0%} of the requested area has imagery for {date}; the "
            f"rest is transparent. The result was not cached."
        )
        return raster

    _write_geotiff(out_path, mosaic, transform)
    info(f"Saved imagery to {out_path}")
    return raster
