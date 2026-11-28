#!/usr/bin/env python3

import json
import os
import re
import time
from datetime import UTC, datetime
from email.utils import parsedate_to_datetime

import geopandas as gpd
import requests
from platformdirs import user_cache_dir
from requests.adapters import HTTPAdapter
from shapely.geometry import LineString, Polygon, box
from shapely.ops import polygonize, unary_union
from urllib3.util.retry import Retry

from ...reproject.bounds import bounds_to_wgs84

# ``warning`` is re-exported for tests that assert nothing warns here.
from .logging import debug, error, info, warning  # noqa: F401

# ------------------------------------------------------------------------
# 1) Global constants/paths
# ------------------------------------------------------------------------
BASE_CACHE_DIR = user_cache_dir(appname="dtcc-data")
os.makedirs(BASE_CACHE_DIR, exist_ok=True)
CACHE_DIR = os.path.join(BASE_CACHE_DIR, "downloaded_osm")
os.makedirs(CACHE_DIR, exist_ok=True)
# Where we store local cache metadata
CACHE_METADATA_FILE = os.path.join(BASE_CACHE_DIR, "cache_metadata.json")

# Each entry is a separate server. Do not add DNS aliases of a listed server
# (overpass.kumi.systems points at overpass.private.coffee): failing over to
# the same server only adds waiting.
OVERPASS_ENDPOINTS = [
    "https://overpass-api.de/api/interpreter",
    "https://overpass.private.coffee/api/interpreter",
]
ROAD_CACHE_VERSION = 4
BUILDING_CACHE_VERSION = 3
OVERPASS_ENDPOINT_METADATA_KEY = "_dtcc_overpass_endpoint"
OVERPASS_SOURCE_ENDPOINT_COLUMN = "source_endpoint"
OVERPASS_CONNECT_TIMEOUT_SECONDS = float(
    os.environ.get("DTCC_OVERPASS_CONNECT_TIMEOUT", "8")
)
# A fixed read timeout set through the environment. When unset, the read
# timeout follows each query's own [timeout:N] (see _overpass_read_timeout).
_READ_TIMEOUT_ENV = os.environ.get("DTCC_OVERPASS_READ_TIMEOUT")
OVERPASS_READ_TIMEOUT_SECONDS = float(_READ_TIMEOUT_ENV) if _READ_TIMEOUT_ENV else None
# Extra seconds on top of the query's server timeout for waiting in the
# server's queue and transferring the response.
OVERPASS_QUEUE_ALLOWANCE_SECONDS = float(
    os.environ.get("DTCC_OVERPASS_QUEUE_ALLOWANCE", "60")
)
# Overpass applies this server timeout when a query sets no [timeout:N].
OVERPASS_DEFAULT_QUERY_TIMEOUT_SECONDS = 180
# How many more times to ask the same server after a quick 429 or 5xx answer.
OVERPASS_RETRIES = max(0, int(os.environ.get("DTCC_OVERPASS_RETRIES", "2")))
# An error answer that arrives within this many seconds counts as quick. A
# slower error, or a timeout, moves on to the next endpoint without retrying,
# so a struggling server is not waited on several times.
OVERPASS_QUICK_FAILURE_SECONDS = float(
    os.environ.get("DTCC_OVERPASS_QUICK_FAILURE", "10")
)
OVERPASS_RETRY_STATUSES = (429, 500, 502, 503, 504)
# Wait before the first retry; it doubles before each later one.
OVERPASS_RETRY_BACKOFF_SECONDS = 1.0
# The longest Retry-After wait honoured before a retry. A server asking for a
# longer pause is not retried; the next endpoint is tried instead.
OVERPASS_RETRY_AFTER_CAP_SECONDS = float(
    os.environ.get("DTCC_OVERPASS_RETRY_AFTER_CAP", "10")
)
# A query may take one full read timeout plus this many seconds across all
# endpoints, retries and waits, before the failover gives up.
OVERPASS_FAILOVER_ALLOWANCE_SECONDS = float(
    os.environ.get("DTCC_OVERPASS_FAILOVER_ALLOWANCE", "60")
)
OVERPASS_SERVER_TIMEOUT_SECONDS = int(
    os.environ.get("DTCC_OVERPASS_SERVER_TIMEOUT", "25")
)
# The settings statement that starts an Overpass QL query, such as
# "[out:json][timeout:25];", after any leading whitespace and comments.
# Possessive repeats keep a failed match linear in the query length.
_QL_SPACE = r"(?:\s|/\*.*?\*/|//[^\n]*+(?:\n|\Z))*+"
_QUERY_SETTINGS_PATTERN = re.compile(
    rf"\A{_QL_SPACE}((?:\[[^\]]*+\]{_QL_SPACE})++);", re.DOTALL
)
# Time left that is worth sending a request in; with less, the failover stops.
_MIN_REQUEST_SECONDS = 1.0
_QUERY_TIMEOUT_PATTERN = re.compile(r"\[\s*timeout\s*:\s*(\d+)\s*\]")
OVERPASS_USER_AGENT = os.environ.get(
    "DTCC_OVERPASS_USER_AGENT",
    "dtcc-core/0.9.8 (https://github.com/dtcc-platform/dtcc-core)",
)


def create_retry_session(
    retries=0,
    backoff_factor=1.0,
    status_forcelist=OVERPASS_RETRY_STATUSES,
):
    """
    Create a requests session for Overpass, optionally retrying HTTP errors.

    Parameters
    ----------
    retries
        Number of retry attempts for transient HTTP status responses. The
        default is 0 because query_overpass_with_failover retries quick
        failures itself. Connect and read failures are never retried here;
        endpoint failover handles those so DNS/connect issues do not block
        for minutes.
    backoff_factor
        Multiplier for exponential backoff (1.0 = 1s, 2s, 4s, 8s, 16s)
    status_forcelist
        HTTP status codes that trigger a retry
    """
    session = requests.Session()
    session.headers.update(
        {
            "User-Agent": OVERPASS_USER_AGENT,
            "Accept": "application/json",
        }
    )
    retry = Retry(
        total=retries,
        read=0,
        connect=0,
        status=retries,
        other=0,
        backoff_factor=backoff_factor,
        status_forcelist=status_forcelist,
        allowed_methods=["GET", "POST"],
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("http://", adapter)
    session.mount("https://", adapter)
    return session


def _overpass_read_timeout(query):
    """
    Return how long to wait for an Overpass response to ``query``, in seconds.

    The server may spend up to the query's ``[timeout:N]`` working, after
    waiting for a free slot, so the client waits that long plus
    ``OVERPASS_QUEUE_ALLOWANCE_SECONDS``. ``[timeout:N]`` is read only from the
    settings statement that starts the query, after any whitespace and
    comments, as Overpass reads it. A query without one gets the server
    default of 180 s. ``DTCC_OVERPASS_READ_TIMEOUT`` replaces the result with
    a fixed value.
    """
    if OVERPASS_READ_TIMEOUT_SECONDS is not None:
        return OVERPASS_READ_TIMEOUT_SECONDS
    settings = _QUERY_SETTINGS_PATTERN.match(query)
    match = settings and _QUERY_TIMEOUT_PATTERN.search(settings.group(1))
    server_timeout = (
        int(match.group(1)) if match else OVERPASS_DEFAULT_QUERY_TIMEOUT_SECONDS
    )
    return server_timeout + OVERPASS_QUEUE_ALLOWANCE_SECONDS


def _retry_after_seconds(response):
    """
    Return the wait a ``Retry-After`` header asks for, in seconds, or None.

    The header holds either a whole number of seconds or an HTTP date; a date
    without a zone, as in the obsolete asctime form, is UTC. A missing or
    unreadable header gives None; a date in the past gives 0.
    """
    value = response.headers.get("Retry-After")
    if value is None:
        return None
    value = value.strip()
    if value.isascii() and value.isdigit():
        return float(value)
    try:
        when = parsedate_to_datetime(value)
    except (TypeError, ValueError, OverflowError):
        return None
    if when.tzinfo is None:
        when = when.replace(tzinfo=UTC)
    return max(0.0, (when - datetime.now(UTC)).total_seconds())


def _quick_retry_delay(response, elapsed, attempt):
    """
    Return the wait before asking the same endpoint again, and an error note.

    The wait is None when the answer should move on to the next endpoint: it
    is not a 429 or 5xx, it took ``OVERPASS_QUICK_FAILURE_SECONDS`` or longer,
    or its ``Retry-After`` asks for more than
    ``OVERPASS_RETRY_AFTER_CAP_SECONDS``. Otherwise it is the exponential
    backoff for ``attempt``, or the ``Retry-After`` wait when that is longer.
    """
    if (
        response.status_code not in OVERPASS_RETRY_STATUSES
        or elapsed >= OVERPASS_QUICK_FAILURE_SECONDS
    ):
        return None, ""
    delay = OVERPASS_RETRY_BACKOFF_SECONDS * 2**attempt
    retry_after = _retry_after_seconds(response)
    if retry_after is None:
        return delay, ""
    if retry_after > OVERPASS_RETRY_AFTER_CAP_SECONDS:
        return None, f" (asked to retry after {retry_after:g}s)"
    return max(delay, retry_after), ""


def _ask_endpoint(session, endpoint, query, request_timeout, attempt):
    """
    Send ``query`` to ``endpoint`` once.

    Returns ``(data, None, None)`` for a usable answer, otherwise
    ``(None, error, delay)``, where ``delay`` is the wait before asking the
    same endpoint again, or None to move on (see :func:`_quick_retry_delay`).
    """
    started = time.monotonic()
    try:
        resp = session.post(endpoint, data={"data": query}, timeout=request_timeout)
    except requests.exceptions.RequestException as e:
        debug(f"Endpoint {endpoint} failed: {e}")
        return None, str(e), None
    elapsed = time.monotonic() - started
    if resp.status_code == 200:
        try:
            data = resp.json()
        except ValueError as e:
            debug(f"Endpoint {endpoint} returned invalid JSON: {e}")
            return None, f"invalid JSON from {endpoint}: {e}", None
        debug(f"Overpass query successful via {endpoint}")
        data[OVERPASS_ENDPOINT_METADATA_KEY] = endpoint
        return data, None, None
    debug(
        f"Endpoint {endpoint} returned status {resp.status_code} after {elapsed:.1f}s"
    )
    delay, note = _quick_retry_delay(resp, elapsed, attempt)
    return None, f"HTTP {resp.status_code}{note}", delay


def query_overpass_with_failover(
    query,
    timeout=None,
    connect_timeout=OVERPASS_CONNECT_TIMEOUT_SECONDS,
    read_timeout=None,
):
    """
    Query Overpass API with automatic failover between endpoints.

    Each endpoint is asked up to ``1 + OVERPASS_RETRIES`` times. A 429 or 5xx
    answer that arrives within ``OVERPASS_QUICK_FAILURE_SECONDS`` is retried on
    the same endpoint after a short backoff, or after the answer's
    ``Retry-After`` wait when that is longer and at most
    ``OVERPASS_RETRY_AFTER_CAP_SECONDS``; a longer requested wait moves on to
    the next endpoint. A slower error, a timeout, a connection failure or
    invalid JSON also moves on.

    The whole call is limited to the read timeout plus
    ``OVERPASS_FAILOVER_ALLOWANCE_SECONDS``: each request's connect and read
    timeouts are cut to the time left, no request starts with less than a
    second left, a retry whose wait would leave less than that is skipped, and
    endpoints not reached in time are listed in the error. The read timeout
    limits each wait for data, not a whole slowly streamed answer, so a call
    can run a little past the limit.

    Parameters
    ----------
    query
        The Overpass QL query string
    timeout
        Backward-compatible read timeout override in seconds. Prefer
        connect_timeout/read_timeout for new code.
    connect_timeout
        Timeout for DNS/TCP/TLS connection setup in seconds.
    read_timeout
        Timeout while waiting for the Overpass response in seconds. When
        None, it follows the query's ``[timeout:N]``; see
        :func:`_overpass_read_timeout`.

    Returns
    -------
    dict
        Parsed JSON response from the API

    Raises
    ------
    RuntimeError
        If all endpoints fail
    """
    if timeout is not None:
        read_timeout = timeout
    if read_timeout is None:
        read_timeout = _overpass_read_timeout(query)

    time_limit = read_timeout + OVERPASS_FAILOVER_ALLOWANCE_SECONDS
    deadline = time.monotonic() + time_limit
    last_error = None
    # Each endpoint's error, so callers can tell a broken query from busy mirrors.
    errors = []
    not_tried = []
    for endpoint in OVERPASS_ENDPOINTS:
        if deadline - time.monotonic() < _MIN_REQUEST_SECONDS:
            not_tried.append(endpoint)
            continue
        info(
            "Trying Overpass endpoint "
            f"{endpoint} (connect={connect_timeout:g}s, read={read_timeout:g}s)"
        )
        session = create_retry_session()
        delay = None
        for attempt in range(OVERPASS_RETRIES + 1):
            if attempt > 0:
                if time.monotonic() + delay + _MIN_REQUEST_SECONDS > deadline:
                    info(
                        f"Not retrying {endpoint}: waiting {delay:g}s would pass "
                        f"the {time_limit:g}s time limit"
                    )
                    break
                info(f"Retrying {endpoint} in {delay:g}s after {last_error}")
                time.sleep(delay)
            # A sleep can overrun slightly; never send a timeout below the minimum.
            remaining = max(deadline - time.monotonic(), _MIN_REQUEST_SECONDS)
            request_timeout = (
                min(connect_timeout, remaining),
                min(read_timeout, remaining),
            )
            data, last_error, delay = _ask_endpoint(
                session, endpoint, query, request_timeout, attempt
            )
            if data is not None:
                return data
            if delay is None:
                break
        errors.append((endpoint, last_error))

    by_endpoint = "\n".join(f"- {endpoint}: {problem}" for endpoint, problem in errors)
    skipped = (
        f"\nNot tried within the {time_limit:g}s time limit: {', '.join(not_tried)}"
        if not_tried
        else ""
    )
    error(
        f"All Overpass endpoints failed. Last error: {last_error}\n"
        f"Errors by endpoint:\n{by_endpoint}{skipped}"
    )
    raise RuntimeError(f"All Overpass endpoints failed: {last_error}")


def _overpass_source_endpoint(payload):
    return payload.get(OVERPASS_ENDPOINT_METADATA_KEY)


def _set_overpass_source_endpoint(gdf, source_endpoint):
    if source_endpoint:
        gdf[OVERPASS_SOURCE_ENDPOINT_COLUMN] = source_endpoint
        gdf.attrs[OVERPASS_SOURCE_ENDPOINT_COLUMN] = source_endpoint
    return gdf


def _gdf_source_endpoint(gdf):
    source_endpoint = gdf.attrs.get(OVERPASS_SOURCE_ENDPOINT_COLUMN)
    if source_endpoint:
        return source_endpoint
    if OVERPASS_SOURCE_ENDPOINT_COLUMN in gdf.columns and len(gdf) > 0:
        source_endpoint = gdf[OVERPASS_SOURCE_ENDPOINT_COLUMN].iloc[0]
        return source_endpoint if source_endpoint else None
    return None


# ------------------------------------------------------------------------
# 2) Utilities for bounding boxes
# ------------------------------------------------------------------------
def is_superset_bbox(bbox_sup, bbox_sub):
    """
    Return True if bbox_sub is fully contained within bbox_sup.
    Each bbox is (xmin, ymin, xmax, ymax).
    """
    xminS, yminS, xmaxS, ymaxS = bbox_sup
    xminT, yminT, xmaxT, ymaxT = bbox_sub
    return xminS <= xminT and yminS <= yminT and xmaxS >= xmaxT and ymaxS >= ymaxT


def filter_gdf_to_bbox(gdf, bbox_3006):
    """
    Filter a GeoDataFrame (already in EPSG:3006) to the specified bounding box by intersection.
    Returns a copy of the subset.
    """
    minx, miny, maxx, maxy = bbox_3006
    bbox_poly = box(minx, miny, maxx, maxy)  # shapely
    return gdf[gdf.geometry.intersects(bbox_poly)].copy()


# Extra OSM way tags kept as road attributes, mapped to column names.
ROAD_EXTRA_TAG_COLUMNS = {
    "surface": "surface",
    "lit": "lit",
    "sidewalk": "sidewalk",
    "sidewalk:left": "sidewalk_left",
    "sidewalk:right": "sidewalk_right",
    "sidewalk:both": "sidewalk_both",
    "cycleway": "cycleway",
    "cycleway:left": "cycleway_left",
    "cycleway:right": "cycleway_right",
    "cycleway:both": "cycleway_both",
    "footway": "footway",
    "width": "width",
    "access": "access",
    "service": "service",
    "foot": "foot",
    "bicycle": "bicycle",
    "vehicle": "vehicle",
    "motor_vehicle": "motor_vehicle",
    "motorcar": "motorcar",
    "motorroad": "motorroad",
}

# Point features on road nodes kept as segment end attributes.
ROAD_NODE_HIGHWAYS = ("crossing", "traffic_signals")

# OSM building tags kept as building attributes, mapped to column names.
# ``height`` is stored as ``osm_height`` so the raw tag text is not mistaken
# for the numeric building height used elsewhere in dtcc-core.
BUILDING_TAG_COLUMNS = {
    "building": "building",
    "building:levels": "building_levels",
    "roof:shape": "roof_shape",
    "height": "osm_height",
    "addr:street": "addr_street",
    "addr:housenumber": "addr_housenumber",
    "addr:postcode": "addr_postcode",
    "addr:city": "addr_city",
}


def _parse_bool_tag(value):
    if value is None:
        return False
    return str(value).strip().lower() in {"yes", "true", "1"}


def _is_reverse_oneway(value):
    if value is None:
        return False
    return str(value).strip().lower() in {"-1", "reverse"}


def _is_oneway(tags):
    oneway = tags.get("oneway")
    if _parse_bool_tag(oneway) or _is_reverse_oneway(oneway):
        return True
    return tags.get("junction") == "roundabout"


def _parse_maxspeed_kmh(value):
    if value is None:
        return None
    text = str(value).strip().lower()
    if not text:
        return None
    match = re.search(r"(\d+(?:\.\d+)?)", text)
    if match is None:
        return None
    speed = float(match.group(1))
    if "mph" in text:
        speed *= 1.609344
    return speed


# ------------------------------------------------------------------------
# 3) Metadata I/O
# ------------------------------------------------------------------------
def load_cache_metadata(meta_path=CACHE_METADATA_FILE):
    """Load the Overpass download cache index.

    Parameters
    ----------
    meta_path : str, optional
        Path of the JSON index file. Defaults to the cache's metadata file.

    Returns
    -------
    list[dict]
        Cache records, with keys such as ``type``, ``bbox``, ``filepath`` and
        ``layer``. Empty when the file does not exist.
    """
    if not os.path.exists(meta_path):
        return []
    with open(meta_path, encoding="utf-8") as f:
        return json.load(f)


def save_cache_metadata(records, meta_path=CACHE_METADATA_FILE):
    """Write the Overpass download cache index.

    Parameters
    ----------
    records : list[dict]
        Cache records to store, as returned by :func:`load_cache_metadata`.
    meta_path : str, optional
        Path of the JSON index file. Defaults to the cache's metadata file.
    """
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(records, f, indent=2)


def to_relative_path(filepath):
    """Convert absolute path to relative path from BASE_CACHE_DIR."""
    if os.path.isabs(filepath) and filepath.startswith(BASE_CACHE_DIR):
        return os.path.relpath(filepath, BASE_CACHE_DIR)
    return filepath


def to_absolute_path(filepath):
    """Convert relative path to absolute path using BASE_CACHE_DIR.
    Handles both relative and already-absolute paths for backward compatibility.
    """
    if os.path.isabs(filepath):
        return filepath  # Already absolute (legacy format)
    return os.path.join(BASE_CACHE_DIR, filepath)


def find_superset_record(bbox_3006, records):
    """
    Look for a record in 'records' whose bounding box is a superset of bbox_3006.
    Return the first match, or None if none found.
    """
    for rec in records:
        rbox = rec["bbox"]  # (xmin, ymin, xmax, ymax)
        if is_superset_bbox(rbox, bbox_3006):
            return rec
    return None


# ------------------------------------------------------------------------
# 4) Overpass logic
# ------------------------------------------------------------------------
def download_overpass_buildings(bbox_3006):
    """
    1) Transform bbox_3006 -> EPSG:4326 (lat/lon).
    2) Query Overpass for building footprints in that bounding box.
    3) Return a GeoDataFrame in EPSG:3006.
    """
    min_lon, min_lat, max_lon, max_lat = bounds_to_wgs84(bbox_3006, "EPSG:3006")

    query = f"""
    [out:json][timeout:{OVERPASS_SERVER_TIMEOUT_SECONDS}];
    (
      way["building"]({min_lat},{min_lon},{max_lat},{max_lon});
      relation["building"]["type"="multipolygon"]({min_lat},{min_lon},{max_lat},{max_lon});
    );
    out geom;
    """
    info(f"Querying Overpass for buildings in bbox={bbox_3006}")
    data = query_overpass_with_failover(query)
    source_endpoint = _overpass_source_endpoint(data)

    # Parse
    nodes = {}
    for elem in data.get("elements", []):
        if elem["type"] == "node":
            nid = elem["id"]
            lat = elem["lat"]
            lon = elem["lon"]
            nodes[nid] = (lat, lon)

    elements = data.get("elements", [])
    relations, relation_outlines = _building_relations(elements)

    building_rows = []
    polygons_4326 = []
    for elem in elements:
        if elem["type"] == "way":
            if elem["id"] in relation_outlines or "nodes" not in elem:
                continue
            coords = _way_coordinates_latlon(elem, nodes)
            if len(coords) <= 2:
                continue
            if coords[0] != coords[-1]:
                coords.append(coords[0])  # close ring
            polygon = Polygon([(lon, lat) for (lat, lon) in coords])
            tags = elem.get("tags", {})
        elif elem["type"] == "relation" and elem["id"] in relations:
            polygon, tags = relations[elem["id"]]
        else:
            continue

        row = {"osm_id": elem["id"], "osm_type": elem["type"]}
        for tag, column in BUILDING_TAG_COLUMNS.items():
            row[column] = tags.get(tag)
        building_rows.append(row)
        polygons_4326.append(polygon)

    building_columns = ["osm_id", "osm_type", *BUILDING_TAG_COLUMNS.values()]
    gdf_4326 = gpd.GeoDataFrame(
        building_rows,
        columns=building_columns,
        geometry=polygons_4326,
        crs="EPSG:4326",
    )
    _set_overpass_source_endpoint(gdf_4326, source_endpoint)
    gdf_3006 = gdf_4326.to_crs("EPSG:3006")
    _set_overpass_source_endpoint(gdf_3006, source_endpoint)
    return gdf_3006


def download_overpass_roads(bbox_3006):
    """
    1) Transform bbox_3006 -> EPSG:4326 (lat/lon).
    2) Query Overpass for roads (highways) in that bounding box.
    3) Return a GeoDataFrame in EPSG:3006.
    """
    min_lon, min_lat, max_lon, max_lat = bounds_to_wgs84(bbox_3006, "EPSG:3006")

    node_highways = "|".join(ROAD_NODE_HIGHWAYS)
    query = f"""
    [out:json][timeout:{OVERPASS_SERVER_TIMEOUT_SECONDS}];
    way["highway"]({min_lat},{min_lon},{max_lat},{max_lon})->.roads;
    .roads out geom;
    node(w.roads)["highway"~"^({node_highways})$"];
    out;
    """
    info(f"Querying Overpass for roads in bbox={bbox_3006}")
    data = query_overpass_with_failover(query)
    source_endpoint = _overpass_source_endpoint(data)

    # Parse
    nodes = {}
    for elem in data.get("elements", []):
        if elem["type"] == "node":
            nid = elem["id"]
            lat = elem["lat"]
            lon = elem["lon"]
            nodes[nid] = (lat, lon)

    node_tags = {
        elem["id"]: elem["tags"]
        for elem in data.get("elements", [])
        if elem["type"] == "node" and elem.get("tags")
    }

    road_rows = []
    road_geometries = []
    for elem in data.get("elements", []):
        if elem["type"] != "way" or "nodes" not in elem:
            continue

        tags = elem.get("tags", {})
        refs, coords = _way_node_ids_and_coordinates_latlon(elem, nodes)
        if len(coords) < 2:
            continue

        reverse = _is_reverse_oneway(tags.get("oneway"))
        if reverse:
            refs = list(reversed(refs))
            coords = list(reversed(coords))

        for segment_index, ((start_lat, start_lon), (end_lat, end_lon)) in enumerate(
            zip(coords[:-1], coords[1:])
        ):
            start_node = refs[segment_index] if segment_index < len(refs) else None
            end_node = (
                refs[segment_index + 1] if segment_index + 1 < len(refs) else None
            )
            road_geometries.append(
                LineString([(start_lon, start_lat), (end_lon, end_lat)])
            )
            road_rows.append(
                {
                    "osm_way_id": elem["id"],
                    "osm_start_node_id": start_node,
                    "osm_end_node_id": end_node,
                    "segment_index": segment_index,
                    "highway": tags.get("highway"),
                    "name": tags.get("name"),
                    "lanes": tags.get("lanes"),
                    "maxspeed": tags.get("maxspeed"),
                    "maxspeed_kmh": _parse_maxspeed_kmh(tags.get("maxspeed")),
                    "oneway": _is_oneway(tags),
                    "oneway_raw": tags.get("oneway"),
                    "bridge": tags.get("bridge"),
                    "tunnel": tags.get("tunnel"),
                    "junction": tags.get("junction"),
                    **{
                        column: tags.get(tag)
                        for tag, column in ROAD_EXTRA_TAG_COLUMNS.items()
                    },
                    **_road_node_columns("start", node_tags.get(start_node)),
                    **_road_node_columns("end", node_tags.get(end_node)),
                }
            )

    road_columns = [
        "osm_way_id",
        "osm_start_node_id",
        "osm_end_node_id",
        "segment_index",
        "highway",
        "name",
        "lanes",
        "maxspeed",
        "maxspeed_kmh",
        "oneway",
        "oneway_raw",
        "bridge",
        "tunnel",
        "junction",
        *ROAD_EXTRA_TAG_COLUMNS.values(),
        "start_node_highway",
        "start_node_crossing",
        "end_node_highway",
        "end_node_crossing",
    ]
    gdf_4326 = gpd.GeoDataFrame(
        road_rows, columns=road_columns, geometry=road_geometries, crs="EPSG:4326"
    )
    _set_overpass_source_endpoint(gdf_4326, source_endpoint)
    gdf_3006 = gdf_4326.to_crs("EPSG:3006")
    _set_overpass_source_endpoint(gdf_3006, source_endpoint)
    return gdf_3006


def _road_node_columns(end, tags):
    """Crossing and traffic-signal tags of a segment's start or end node."""
    tags = tags or {}
    highway = tags.get("highway")
    return {
        f"{end}_node_highway": highway if highway in ROAD_NODE_HIGHWAYS else None,
        f"{end}_node_crossing": tags.get("crossing"),
    }


def _way_coordinates_latlon(elem, nodes):
    geometry = elem.get("geometry")
    if geometry:
        return [(point["lat"], point["lon"]) for point in geometry]

    coords = []
    for ref in elem.get("nodes", []):
        if ref in nodes:
            coords.append(nodes[ref])  # (lat, lon)
    return coords


def _building_relations(elements):
    """Assemble the building relations in an Overpass response.

    Returns ``{relation id: (polygon, tags)}`` for the relations that form a
    closed shape, and the ids of the building ways that outline them. Those
    ways are described by their relation, so they are not added a second
    time; their tags fill in what the relation leaves out.
    """
    way_tags = {
        elem["id"]: elem.get("tags", {}) for elem in elements if elem["type"] == "way"
    }
    relations = {}
    outlines = set()
    for elem in elements:
        if elem["type"] != "relation":
            continue
        polygon = _relation_polygon_lonlat(elem)
        if polygon is None:
            continue
        refs = [
            member["ref"]
            for member in elem.get("members", [])
            if member.get("type") == "way"
            and member.get("role") != "inner"
            and member.get("ref") in way_tags
        ]
        outlines.update(refs)
        tags = _relation_building_tags(
            elem.get("tags", {}), [way_tags[ref] for ref in refs]
        )
        relations[elem["id"]] = (polygon, tags)
    return relations, outlines


def _relation_building_tags(relation_tags, outline_tags):
    """Return a building relation's tags, filled in from its outline ways.

    The relation's own value wins, except ``building=yes``, which a more
    specific building type on the outlines replaces. A tag the relation
    lacks is taken from the outline ways only when every one of them has
    it with the same value.
    """
    tags = dict(relation_tags)
    for tag in BUILDING_TAG_COLUMNS:
        current = tags.get(tag)
        if current is not None and not (tag == "building" and current == "yes"):
            continue
        values = {way.get(tag) for way in outline_tags}
        if len(values) == 1 and None not in values:
            tags[tag] = values.pop()
    return tags


def _member_rings(lines):
    """Polygons for the rings formed by relation member lines.

    A member that is closed on its own is used as drawn, so touching rings
    do not also enclose the space between them. Members that are pieces of
    a ring are joined with ``polygonize``; it turns nested rings into faces
    with holes, so each face's exterior is kept to recover the rings.
    """
    closed = [Polygon(line.coords) for line in lines if line.is_ring]
    pieces = [line for line in lines if not line.is_ring]
    joined = (
        [Polygon(face.exterior) for face in polygonize(unary_union(pieces))]
        if pieces
        else []
    )
    return closed + joined


def _relation_polygon_lonlat(elem):
    """Build a (Multi)Polygon in lon/lat from a multipolygon relation.

    Outer and inner rings may each be split over several member ways, which
    are joined first; a member that is closed on its own is used as drawn.
    Each inner ring is cut from every outer ring that contains it. An outer
    ring inside a courtyard does not contain the courtyard ring, so it stays
    part of the building. An inner ring that no outer ring contains is cut
    from every outer ring it overlaps. Members without the ``inner`` role
    count as outer rings. Returns ``None`` when no closed outer ring can be
    formed.
    """
    outer_lines = []
    inner_lines = []
    for member in elem.get("members", []):
        if member.get("type") != "way":
            continue
        points = [(p["lon"], p["lat"]) for p in member.get("geometry") or [] if p]
        if len(points) < 2:
            continue
        lines = inner_lines if member.get("role") == "inner" else outer_lines
        lines.append(LineString(points))

    outer = _member_rings(outer_lines)
    if not outer:
        return None
    inner = _member_rings(inner_lines)

    holes = [[] for _ in outer]
    for ring in inner:
        containing = [i for i, polygon in enumerate(outer) if polygon.covers(ring)]
        if containing:
            for i in containing:
                holes[i].append(ring)
        else:
            for i, polygon in enumerate(outer):
                if polygon.intersects(ring):
                    holes[i].append(ring)

    parts = [
        polygon.difference(unary_union(rings)) if rings else polygon
        for polygon, rings in zip(outer, holes)
    ]
    shape = unary_union(parts)
    if shape.is_empty:
        return None
    return shape


def _way_node_ids_and_coordinates_latlon(elem, nodes):
    refs = list(elem.get("nodes", []))
    coords = _way_coordinates_latlon(elem, nodes)

    if not elem.get("geometry"):
        refs = [ref for ref in refs if ref in nodes]

    if len(refs) != len(coords):
        refs = [None] * len(coords)

    return refs, coords


# ------------------------------------------------------------------------
# 5) Superset-based caching logic for Buildings
# ------------------------------------------------------------------------
def get_buildings_for_bbox(bbox_3006):
    """
    1) load metadata
    2) look for superset => filter local
    3) otherwise => Overpass => store => add to metadata
    4) return GDF in EPSG:3006
    """
    records = load_cache_metadata()
    building_records = [
        r
        for r in records
        if r["type"] == "buildings" and r.get("version") == BUILDING_CACHE_VERSION
    ]
    sup_rec = find_superset_record(bbox_3006, building_records)
    if sup_rec:
        debug("Found superset bounding box for buildings:", sup_rec["bbox"])
        gdf_all = gpd.read_file(
            to_absolute_path(sup_rec["filepath"]), layer=sup_rec["layer"]
        )
        _set_overpass_source_endpoint(gdf_all, sup_rec.get("source_endpoint"))
        subset_gdf = filter_gdf_to_bbox(gdf_all, bbox_3006)
        _set_overpass_source_endpoint(subset_gdf, sup_rec.get("source_endpoint"))
        saved_filename = to_absolute_path(sup_rec["filepath"])
        info(
            f"Subset size: {len(subset_gdf)} features for buildings in bbox={bbox_3006}"
        )
        return subset_gdf, saved_filename
    else:
        # Overpass
        new_gdf = download_overpass_buildings(bbox_3006)
        info(f"Downloaded {len(new_gdf)} building footprints from Overpass.")
        # store
        out_filename = os.path.join(
            CACHE_DIR,
            f"buildings_{bbox_3006[0]}_{bbox_3006[1]}_{bbox_3006[2]}_{bbox_3006[3]}.gpkg",
        )
        new_gdf.to_file(out_filename, layer="buildings", driver="GPKG")
        # update metadata
        record = {
            "type": "buildings",
            "version": BUILDING_CACHE_VERSION,
            "bbox": list(bbox_3006),
            "filepath": to_relative_path(out_filename),
            "layer": "buildings",
        }
        source_endpoint = _gdf_source_endpoint(new_gdf)
        if source_endpoint:
            record["source_endpoint"] = source_endpoint
        records.append(record)
        save_cache_metadata(records)
        return new_gdf, out_filename


# ------------------------------------------------------------------------
# 6) Superset-based caching logic for Roads
# ------------------------------------------------------------------------
def get_roads_for_bbox(bbox_3006):
    """
    1) load metadata
    2) look for superset => filter local
    3) otherwise => Overpass => store => add to metadata
    4) return GDF in EPSG:3006
    """
    records = load_cache_metadata()
    road_records = [
        r
        for r in records
        if r["type"] == "roads" and r.get("version") == ROAD_CACHE_VERSION
    ]
    sup_rec = find_superset_record(bbox_3006, road_records)
    if sup_rec:
        debug("Found superset bounding box for roads:", sup_rec["bbox"])
        gdf_all = gpd.read_file(
            to_absolute_path(sup_rec["filepath"]), layer=sup_rec["layer"]
        )
        _set_overpass_source_endpoint(gdf_all, sup_rec.get("source_endpoint"))
        saved_filename = to_absolute_path(sup_rec["filepath"])
        subset_gdf = filter_gdf_to_bbox(gdf_all, bbox_3006)
        _set_overpass_source_endpoint(subset_gdf, sup_rec.get("source_endpoint"))
        info(f"Subset size: {len(subset_gdf)} features for roads in bbox={bbox_3006}")
        return subset_gdf, saved_filename
    else:
        # Overpass
        new_gdf = download_overpass_roads(bbox_3006)
        info(f"Downloaded {len(new_gdf)} road features from Overpass.")
        # store
        out_filename = os.path.join(
            CACHE_DIR,
            f"roads_{bbox_3006[0]}_{bbox_3006[1]}_{bbox_3006[2]}_{bbox_3006[3]}.gpkg",
        )
        new_gdf.to_file(out_filename, layer="roads", driver="GPKG")
        # update metadata
        record = {
            "type": "roads",
            "version": ROAD_CACHE_VERSION,
            "bbox": list(bbox_3006),
            "filepath": to_relative_path(out_filename),
            "layer": "roads",
        }
        source_endpoint = _gdf_source_endpoint(new_gdf)
        if source_endpoint:
            record["source_endpoint"] = source_endpoint
        records.append(record)
        save_cache_metadata(records)
        return new_gdf, out_filename


# ------------------------------------------------------------------------
# Example usage
# ------------------------------------------------------------------------
"""
if __name__ == "__main__":
    # bigger bounding box in EPSG:3006
    bboxA = (267000, 6519000, 270000, 6521000)
    # smaller bounding box, fully inside A
    bboxB = (268000, 6519500, 269000, 6520000)

    print("=== Buildings: BBox A ===")
    bldgA = get_buildings_for_bbox(bboxA)
    print("Buildings A size:", len(bldgA))

    print("=== Buildings: BBox B (subset) ===")
    bldgB = get_buildings_for_bbox(bboxB)
    print("Buildings B size:", len(bldgB))

    print("=== Roads: BBox A ===")
    roadsA = get_roads_for_bbox(bboxA)
    print("Roads A size:", len(roadsA))

    print("=== Roads: BBox B (subset) ===")
    roadsB = get_roads_for_bbox(bboxB)
    print("Roads B size:", len(roadsB))
"""
