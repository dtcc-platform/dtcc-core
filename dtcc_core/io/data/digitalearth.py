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
from datetime import datetime, timedelta, timezone

import numpy as np
import pyproj
import rasterio
from platformdirs import user_cache_dir
from rasterio.io import MemoryFile
from rasterio.transform import from_origin

from .logging import info, warning, debug, error
from .overpass import create_retry_session

# ------------------------------------------------------------------------
# 1) Global constants/paths
# ------------------------------------------------------------------------
BASE_CACHE_DIR = user_cache_dir(appname="dtcc-data")
CACHE_DIR = os.path.join(BASE_CACHE_DIR, "downloaded_des")
os.makedirs(CACHE_DIR, exist_ok=True)

WMS_URL = "https://ows.digitalearth.se"
STAC_URL = "https://explorer.digitalearth.se/stac"
STAC_COLLECTION = "s2_msi_l2a"

# The server advertises MaxWidth/MaxHeight of 512 in its WMS capabilities,
# so anything larger has to be requested as a grid of tiles and stitched.
MAX_TILE_PIXELS = 512

# Native ground resolution of the Sentinel-2 visual bands, in metres.
DEFAULT_RESOLUTION = 10.0

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

# Two layers cover the same archive but are tuned differently. Despite its
# name, 'sentinel2_reflectance2022' is the one for imagery *before* 2022.
LAYER_FROM_2022 = "sentinel2_reflectance"
LAYER_BEFORE_2022 = "sentinel2_reflectance2022"

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
    lons, lats = zip(*(transformer.transform(x, y) for x, y in corners))
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
            raise RuntimeError(
                f"STAC search failed with status {resp.status_code}:\n{resp.text}"
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
    'end' are ISO date strings bounding the search, both optional. Items
    with a cloud cover above 'max_cloud' percent are discarded.

    Returns a list of (date, cloud_cover) tuples sorted newest first, where
    'date' is a 'YYYY-MM-DD' string. A date appearing on several satellite
    tiles is reported once, with the cloudiest of those tiles, so that a
    date is only called clear when the whole area is clear.
    """
    session = session or create_retry_session()
    payload = {
        "collections": [STAC_COLLECTION],
        "bbox": list(_to_wgs84_bbox(bbox)),
        "limit": limit,
    }
    if start or end:
        payload["datetime"] = f"{start or '..'}/{end or '..'}"

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
    style="rgb",
    session=None,
    max_candidates=MAX_DATE_CANDIDATES,
):
    """
    Return the most recent date on which 'bbox' was imaged without much cloud.

    Only the last 'search_days' days are considered, to keep the catalogue
    query small. Each candidate is then probed against the imagery server,
    because the catalogue and the server are not perfectly in step: the
    catalogue lists acquisitions that the server has no pixels for, and
    trusting it blindly yields a blank image. The first candidate that really
    has pixels is returned.

    Raises RuntimeError if nothing clear enough was found, since silently
    returning a cloudy date would hand the caller a white image.
    """
    session = session or create_retry_session()
    start = (datetime.now(timezone.utc) - timedelta(days=search_days)).strftime(
        "%Y-%m-%d"
    )
    candidates = list_available_dates(
        bbox, max_cloud=max_cloud, start=start, session=session
    )
    if not candidates:
        raise RuntimeError(
            f"No Sentinel-2 acquisition below {max_cloud}% cloud cover in the last "
            f"{search_days} days for this area. Pass an explicit date, raise "
            f"max_cloud, or widen search_days."
        )

    for date, cloud in candidates[:max_candidates]:
        if _has_imagery(session, _layer_for_date(date), style, bbox, date):
            info(f"Using acquisition from {date} ({cloud:.1f}% cloud cover)")
            return date
        debug(f"{date} is in the catalogue but the imagery server has no pixels for it")

    raise RuntimeError(
        f"The catalogue lists {len(candidates)} acquisition(s) below {max_cloud}% "
        f"cloud cover for this area, but none of the {max_candidates} most recent "
        f"had imagery on the server. Pass an explicit date or widen search_days."
    )


# ------------------------------------------------------------------------
# 3) Imagery requests (WMS)
# ------------------------------------------------------------------------
def _layer_for_date(date):
    """
    Pick the WMS layer tuned for the year of 'date'.

    The two layers span the same archive but are styled for different
    periods; see LAYER_FROM_2022 above for the naming trap.
    """
    return LAYER_FROM_2022 if int(date[:4]) >= 2022 else LAYER_BEFORE_2022


def _tile_grid(bbox, resolution):
    """
    Split 'bbox' into tiles no larger than MAX_TILE_PIXELS on a side.

    Returns (width, height, tiles) where width and height are the pixel
    dimensions of the whole mosaic and each tile is a dict with its pixel
    offset, pixel size and its own bounding box in EPSG:3006.
    """
    xmin, ymin, xmax, ymax = bbox
    width = max(1, int(math.ceil((xmax - xmin) / resolution)))
    height = max(1, int(math.ceil((ymax - ymin) / resolution)))

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
    if resp.status_code != 200:
        raise RuntimeError(
            f"WMS request failed with status {resp.status_code}:\n{resp.text}"
        )
    # Errors come back as an XML ServiceExceptionReport with status 200.
    if "xml" in resp.headers.get("Content-Type", ""):
        raise RuntimeError(f"WMS returned an error:\n{resp.text}")
    return resp.content


def _tile_to_array(png_bytes, width, height):
    """
    Decode a PNG tile into an (4, height, width) uint8 RGBA array.

    Tiles with no data for the requested date come back as a small greyscale
    image rather than RGBA, and are treated as fully transparent.
    """
    with MemoryFile(png_bytes) as memfile:
        with memfile.open() as src:
            data = src.read()

    if data.shape[0] < 3:
        return np.zeros((4, height, width), dtype=np.uint8)
    if data.shape[0] == 3:
        alpha = np.full((1, data.shape[1], data.shape[2]), 255, dtype=np.uint8)
        data = np.concatenate([data, alpha], axis=0)
    return data[:4].astype(np.uint8)


def _has_imagery(session, layer, style, bbox, date):
    """
    Report whether the WMS server holds any pixels for 'bbox' on 'date'.

    This is a deliberately tiny request, PROBE_PIXELS on a side, because an
    empty answer comes back as a few hundred bytes and a full sized one would
    waste a download only to be thrown away.
    """
    png_bytes = get_map(session, layer, style, bbox, PROBE_PIXELS, PROBE_PIXELS, date)
    return bool(_tile_to_array(png_bytes, PROBE_PIXELS, PROBE_PIXELS)[3].any())


# ------------------------------------------------------------------------
# 4) Entry point
# ------------------------------------------------------------------------
def download_imagery(
    bbox,
    date=None,
    style="rgb",
    max_cloud=DEFAULT_MAX_CLOUD,
    resolution=DEFAULT_RESOLUTION,
    search_days=DEFAULT_SEARCH_DAYS,
    session=None,
):
    """
    Download Sentinel-2 imagery covering 'bbox' and save it as a GeoTIFF.

    'bbox' is a tuple (xmin, ymin, xmax, ymax) in EPSG:3006. If 'date' is
    omitted the most recent acquisition with less than 'max_cloud' percent
    cloud is used. 'style' selects the rendering, see VALID_STYLES, and
    'resolution' is the ground sample distance in metres.

    Returns the path to the written GeoTIFF, which is cached so that asking
    for the same image twice only downloads it once.
    """
    if style not in VALID_STYLES:
        raise ValueError(f"Invalid style '{style}'. Must be one of {VALID_STYLES}.")

    xmin, ymin, xmax, ymax = bbox
    if xmax <= xmin or ymax <= ymin:
        raise ValueError(f"Invalid bounds {bbox}; expected (xmin, ymin, xmax, ymax).")

    session = session or create_retry_session()
    if date is None:
        date = latest_clear_date(
            bbox,
            max_cloud=max_cloud,
            search_days=search_days,
            style=style,
            session=session,
        )

    out_path = os.path.join(
        CACHE_DIR,
        f"des_{style}_{date}_{xmin:.0f}_{ymin:.0f}_{xmax:.0f}_{ymax:.0f}"
        f"_{resolution:g}m.tif",
    )
    if os.path.exists(out_path):
        info(f"Imagery for {date} already in cache, skipping download.")
        return out_path

    layer = _layer_for_date(date)
    width, height, tiles = _tile_grid(bbox, resolution)
    mosaic = np.zeros((4, height, width), dtype=np.uint8)

    for index, tile in enumerate(tiles, start=1):
        info(f"Downloading tile {index}/{len(tiles)} for {date}")
        png_bytes = get_map(
            session, layer, style, tile["bbox"], tile["width"], tile["height"], date
        )
        data = _tile_to_array(png_bytes, tile["width"], tile["height"])
        rows = slice(tile["row_off"], tile["row_off"] + tile["height"])
        cols = slice(tile["col_off"], tile["col_off"] + tile["width"])
        mosaic[:, rows, cols] = data

    covered = float((mosaic[3] > 0).mean())
    if covered == 0.0:
        warning(
            f"The imagery returned for {date} is empty. The area may lie outside "
            f"the Digital Earth Sweden coverage, or have no acquisition that day."
        )
    elif covered < COVERAGE_WARNING_THRESHOLD:
        # A partly covered mosaic looks like a normal image with transparent
        # patches, so say so rather than letting it pass for a complete one.
        warning(
            f"Only {covered:.0%} of the requested area has imagery for {date}; the "
            f"rest is transparent. The acquisition may not span the whole area."
        )

    # The mosaic is anchored at its top left corner, which is (xmin, ymax).
    transform = from_origin(xmin, ymax, resolution, resolution)
    with rasterio.open(
        out_path,
        "w",
        driver="GTiff",
        height=height,
        width=width,
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

    info(f"Saved imagery to {out_path}")
    return out_path
