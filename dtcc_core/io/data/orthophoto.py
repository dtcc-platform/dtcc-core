"""Client for the DTCC LM orthophoto tile server.

The server lists original Lantmäteriet orthophoto GeoTIFFs for an EPSG:3006
bbox (``GET /items``) and streams each file on request (``GET /files/...``).
"""

import hashlib
import json
import math
import os
import re
import tempfile
from collections.abc import Sequence
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from urllib.parse import urlsplit

import rasterio
import requests

from dtcc_core.io.raster_tiles import without_tracebacks
from dtcc_core.model import Bounds
from dtcc_core.model.values.raster_tiles import RasterTile, RasterTileCollection

SERVER_URL_ENV = "DTCC_ORTHOPHOTO_URL"
CRS = "EPSG:3006"
MAX_DETAIL_CHARS = 1024
CREDENTIALS_MISSING_DETAIL = "LM credentials not configured"
BOUNDS_TOLERANCE = 0.001

_IDENTIFIER = re.compile(r"[a-z0-9_-]+")


class OrthophotoClientError(RuntimeError):
    """Failure talking to the tile server.

    ``failure_class`` is one of ``"timeout"``, ``"connection"``,
    ``"http_4xx"``, ``"http_5xx"``, ``"invalid_payload"`` or
    ``"configuration"`` (the server has no LM credentials). ``message`` and
    ``detail`` are truncated to ``MAX_DETAIL_CHARS`` characters.
    """

    def __init__(
        self,
        message: str,
        *,
        operation: str,
        target: str,
        failure_class: str,
        status_code: int | None = None,
        detail: str | None = None,
    ):
        message = message[:MAX_DETAIL_CHARS]
        super().__init__(message)
        self.message = message
        self.operation = operation
        self.target = target
        self.failure_class = failure_class
        self.status_code = status_code
        self.detail = None if detail is None else detail[:MAX_DETAIL_CHARS]


@dataclass(frozen=True)
class OrthophotoItem:
    """One manifest entry. ``bbox`` is the grid cell, which is the file extent."""

    id: str
    collection: str
    bbox: tuple[float, float, float, float]
    datetime: datetime
    spektraltyp: str
    href: str
    resolution: float | None
    size_bytes: int | None


def resolve_server_url(server_url: str | None) -> str:
    """Return the normalized service root.

    An explicit ``server_url`` takes precedence over ``DTCC_ORTHOPHOTO_URL``.
    The root may include a reverse-proxy path prefix; a trailing slash is
    dropped. Raises ValueError when neither is set or the URL is invalid.
    """
    if server_url is None:
        server_url = os.environ.get(SERVER_URL_ENV)
        if server_url is None:
            raise ValueError(
                f"No orthophoto server configured: pass server_url or set "
                f"{SERVER_URL_ENV}."
            )
    return _normalize_server_url(server_url)


def file_url(server_url: str, item: OrthophotoItem) -> str:
    """Return the download URL of ``item`` below the service root."""
    return _normalize_server_url(server_url) + item.href


def fetch_manifest(
    bounds: tuple[float, float, float, float],
    *,
    server_url: str,
    year: int | None,
    collection: str | None,
    spektraltyp: tuple[str, ...],
    timeout: tuple[float, float],
) -> list[OrthophotoItem]:
    """Fetch and validate the manifest of items covering ``bounds``.

    ``bounds`` is ``(xmin, ymin, xmax, ymax)`` in EPSG:3006. An empty
    ``spektraltyp`` leaves the server default. ``timeout`` is requests'
    ``(connect, read)`` pair. Items keep manifest order.
    """
    try:
        finite = len(bounds) == 4 and all(math.isfinite(v) for v in bounds)
    except OverflowError:
        finite = False
    if not finite:
        raise ValueError(f"bounds must be four finite numbers, got {bounds!r}")
    xmin, ymin, xmax, ymax = (float(v) for v in bounds)
    if xmin >= xmax or ymin >= ymax:
        raise ValueError(f"bounds must have min < max, got {bounds!r}")

    url = _normalize_server_url(server_url) + "/items"
    params: dict[str, object] = {
        "bbox": ",".join(repr(v) for v in (xmin, ymin, xmax, ymax))
    }
    if year is not None:
        params["year"] = year
    if collection is not None:
        params["collection"] = collection
    if spektraltyp:
        params["spektraltyp"] = ",".join(spektraltyp)

    def fail(message, failure_class, status_code=None, detail=None):
        return OrthophotoClientError(
            message,
            operation="manifest",
            target=url,
            failure_class=failure_class,
            status_code=status_code,
            detail=detail,
        )

    try:
        response = requests.get(
            url,
            params=params,
            timeout=timeout,
            headers={"Accept": "application/json"},
            allow_redirects=False,
            stream=True,
        )
    except requests.RequestException as error:
        _raise_transport_error(error, fail)
    try:
        _check_status(response, fail)
        try:
            body = response.content
        except requests.RequestException as error:
            _raise_transport_error(error, fail, status_code=response.status_code)
    finally:
        response.close()

    try:
        payload = json.loads(body)
    except ValueError:
        raise fail(
            "Manifest is not valid JSON",
            "invalid_payload",
            status_code=response.status_code,
            detail=_bounded_text(body),
        ) from None
    try:
        return _parse_manifest(payload)
    except ValueError as error:
        raise fail(
            f"Invalid manifest: {error}",
            "invalid_payload",
            status_code=response.status_code,
        ) from None


def download_tile(
    item: OrthophotoItem,
    *,
    server_url: str,
    timeout: tuple[float, float],
    cache_root: Path,
) -> Path:
    """Return the local path of ``item``'s original GeoTIFF, downloading it if needed.

    Files live at ``cache_root/orthophoto/<service>/<collection>/<id>.tif``,
    where ``<service>`` hashes the normalized service root. An existing file is
    returned as is. A missing file is streamed into a unique temporary file in
    the same directory, checked (status, encoding, length, GeoTIFF header, CRS
    and bounds against ``item.bbox``) and then atomically renamed into place.
    Local filesystem errors and unexpected exceptions propagate unchanged.
    """
    base = _normalize_server_url(server_url)
    service = hashlib.sha256(base.encode()).hexdigest()[:16]
    destination = (
        Path(cache_root) / "orthophoto" / service / item.collection / f"{item.id}.tif"
    )
    if destination.is_file():
        return destination

    url = file_url(base, item)

    def fail(message, failure_class, status_code=None, detail=None):
        return OrthophotoClientError(
            message,
            operation="download",
            target=url,
            failure_class=failure_class,
            status_code=status_code,
            detail=detail,
        )

    try:
        response = requests.get(
            url,
            timeout=timeout,
            headers={"Accept-Encoding": "identity"},
            allow_redirects=False,
            stream=True,
        )
    except requests.RequestException as error:
        _raise_transport_error(error, fail)
    try:
        _check_status(response, fail)
        if "Content-Encoding" in response.headers:
            raise fail(
                f"Unexpected Content-Encoding "
                f"{response.headers['Content-Encoding']!r}",
                "invalid_payload",
                status_code=200,
            )
        expected = _content_length(response, fail)
        destination.parent.mkdir(parents=True, exist_ok=True)
        fd, name = tempfile.mkstemp(
            dir=destination.parent, prefix=f".{item.id}.", suffix=".part"
        )
        temporary = Path(name)
        try:
            written = 0
            with open(fd, "wb") as output:
                try:
                    for chunk in response.iter_content(chunk_size=1 << 20):
                        output.write(chunk)
                        written += len(chunk)
                except requests.RequestException as error:
                    _raise_transport_error(error, fail, status_code=200)
            if expected is not None and written < expected:
                raise fail(
                    f"Body ended after {written} of {expected} bytes",
                    "connection",
                    status_code=200,
                )
            if expected is not None and written > expected:
                raise fail(
                    f"Body has {written} bytes, Content-Length is {expected}",
                    "invalid_payload",
                    status_code=200,
                )
            _check_geotiff(temporary, item, fail)
            os.replace(temporary, destination)
        except BaseException:
            try:
                temporary.unlink(missing_ok=True)
            except OSError:
                pass
            raise
    finally:
        response.close()
    return destination


def acquire_tiles(
    bounds: tuple[float, float, float, float],
    *,
    server_url: str,
    year: int | None,
    collection: str | None,
    spektraltyp: tuple[str, ...],
    timeout: tuple[float, float],
    cache_root: Path,
) -> RasterTileCollection:
    """Fetch the manifest for ``bounds`` and return its items as local tiles.

    The manifest is fetched on every call; each item is then downloaded in
    manifest order, or reused from ``cache_root``. The first error stops
    acquisition and propagates; files published before it stay cached.
    """
    items = fetch_manifest(
        bounds,
        server_url=server_url,
        year=year,
        collection=collection,
        spektraltyp=spektraltyp,
        timeout=timeout,
    )
    return download_items(
        items, server_url=server_url, timeout=timeout, cache_root=cache_root
    )


def download_items(
    items: Sequence[OrthophotoItem],
    *,
    server_url: str,
    timeout: tuple[float, float],
    cache_root: Path,
    failures: list | None = None,
) -> RasterTileCollection:
    """Download listed items in order, or reuse them from ``cache_root``.

    Without ``failures`` the first OrthophotoClientError propagates. With a list,
    each one is recorded as ``(item, error)`` without tracebacks, the item is left
    out and the next item is downloaded. Local errors always propagate; files
    published before them stay cached.
    """
    tiles = []
    for item in items:
        try:
            path = download_tile(
                item, server_url=server_url, timeout=timeout, cache_root=cache_root
            )
        except OrthophotoClientError as error:
            if failures is None:
                raise
            failures.append((item, without_tracebacks(error)))
            continue
        tiles.append(
            RasterTile(
                path=path,
                id=item.id,
                collection=item.collection,
                datetime=item.datetime,
                extent=Bounds(*item.bbox),
                crs=CRS,
                spektraltyp=item.spektraltyp,
                resolution=item.resolution,
                size_bytes=item.size_bytes,
            )
        )
    return RasterTileCollection(tiles=tiles)


def _normalize_server_url(server_url: str) -> str:
    # No message repeats the URL, since any part of it may hold a credential. A
    # parser error is replaced outside its handler, so nothing chains to it.
    try:
        parts = urlsplit(server_url)
        parts.port  # raises ValueError unless the port is a number in range
    except ValueError:
        parts = None
    if parts is not None and "@" in parts.netloc:
        raise ValueError(
            "Orthophoto server URL must not contain user information "
            "(user name or password)"
        )
    if parts is None or parts.scheme not in ("http", "https") or not parts.netloc:
        raise ValueError("Orthophoto server URL must be http(s)://host[:port][/prefix]")
    if "?" in server_url or "#" in server_url:
        raise ValueError("Orthophoto server URL must not have a query or fragment")
    return server_url.rstrip("/")


def _raise_transport_error(error, fail, status_code=None):
    """Raise a client error for a transport failure, else re-raise ``error``."""
    if isinstance(error, requests.Timeout):
        failure_class, message = "timeout", f"Request timed out: {error}"
    elif isinstance(
        error, (requests.ConnectionError, requests.exceptions.ChunkedEncodingError)
    ):
        failure_class, message = "connection", f"Connection failed: {error}"
    elif isinstance(error, requests.exceptions.ContentDecodingError):
        failure_class, message = "invalid_payload", f"Undecodable body: {error}"
    else:
        raise error
    raise fail(message, failure_class, status_code=status_code) from None


def _check_status(response, fail):
    status = response.status_code
    if status == 200:
        return
    detail = _bounded_text(_read_bounded(response))
    if 300 <= status < 400:
        location = response.headers.get("Location", "")
        raise fail(
            f"Unexpected redirect (HTTP {status}) to {location!r}",
            "invalid_payload",
            status_code=status,
            detail=detail,
        )
    if status == 503 and _json_detail(detail) == CREDENTIALS_MISSING_DETAIL:
        failure_class = "configuration"
    elif 500 <= status < 600:
        failure_class = "http_5xx"
    elif 400 <= status < 500:
        failure_class = "http_4xx"
    else:
        failure_class = "invalid_payload"
    raise fail(
        f"HTTP {status}: {detail}", failure_class, status_code=status, detail=detail
    )


def _content_length(response, fail) -> int | None:
    value = response.headers.get("Content-Length")
    if value is None:
        return None
    try:
        if re.fullmatch(r"[0-9]+", value):
            return int(value)
    except ValueError:
        pass
    raise fail(
        f"Malformed Content-Length {value!r}", "invalid_payload", status_code=200
    )


def _check_geotiff(path: Path, item: OrthophotoItem, fail) -> None:
    """Check the header of a downloaded file without reading pixels.

    Only the GeoTIFF driver may open the file, so formats that refer to other
    files or remote sources (such as VRT) are never opened, and GDAL sidecar
    files beside it (.msk, .ovr, .aux.xml) are not read.
    """
    try:
        with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR"), rasterio.open(
            path, driver="GTiff"
        ) as source:
            shape = (source.width, source.height, source.count)
            epsg = source.crs.to_epsg() if source.crs is not None else None
            bounds = tuple(source.bounds)
    except rasterio.errors.RasterioError:
        # GDAL's message names the temporary file in the cache; leave it out.
        raise fail(
            "Not a readable GeoTIFF", "invalid_payload", status_code=200
        ) from None
    if min(shape) <= 0:
        raise fail(
            f"Unexpected raster: width/height/bands {shape}",
            "invalid_payload",
            status_code=200,
        )
    if epsg != 3006:
        raise fail(
            f"Raster CRS is EPSG:{epsg}, expected {CRS}",
            "invalid_payload",
            status_code=200,
        )
    # Written as "all within" so that NaN bounds fail the check.
    if not all(abs(a - b) <= BOUNDS_TOLERANCE for a, b in zip(bounds, item.bbox)):
        raise fail(
            f"Raster bounds {bounds} do not match manifest bbox {item.bbox}",
            "invalid_payload",
            status_code=200,
        )


def _read_bounded(response) -> bytes:
    data = b""
    try:
        for chunk in response.iter_content(chunk_size=MAX_DETAIL_CHARS):
            data += chunk
            if len(data) >= MAX_DETAIL_CHARS:
                break
    except requests.RequestException:
        pass
    return data[:MAX_DETAIL_CHARS]


def _bounded_text(data: bytes) -> str:
    return data[:MAX_DETAIL_CHARS].decode("utf-8", errors="replace")


def _json_detail(text: str):
    try:
        payload = json.loads(text)
    except ValueError:
        return None
    return payload.get("detail") if isinstance(payload, dict) else None


def _parse_manifest(payload) -> list[OrthophotoItem]:
    if not isinstance(payload, dict):
        raise ValueError("manifest is not a JSON object")
    if payload.get("crs") != CRS:
        raise ValueError(f"manifest crs is {payload.get('crs')!r}, expected {CRS}")
    items = payload.get("items")
    if not isinstance(items, list):
        raise ValueError("manifest items is not a list")
    return [_parse_item(index, entry) for index, entry in enumerate(items)]


def _parse_item(index: int, entry) -> OrthophotoItem:
    if not isinstance(entry, dict):
        raise ValueError(f"item {index} is not a JSON object")

    def identifier(key):
        value = entry.get(key)
        if not isinstance(value, str) or not _IDENTIFIER.fullmatch(value):
            raise ValueError(f"item {index} has invalid {key} {value!r}")
        return value

    item_id = identifier("id")
    collection = identifier("collection")

    bbox = entry.get("bbox")
    if (
        not isinstance(bbox, list)
        or len(bbox) != 4
        or not all(_is_finite_number(v) for v in bbox)
        or not (bbox[0] < bbox[2] and bbox[1] < bbox[3])
    ):
        raise ValueError(f"item {index} has invalid bbox {bbox!r}")

    text = entry.get("datetime")
    try:
        timestamp = datetime.fromisoformat(text) if isinstance(text, str) else None
    except ValueError:
        timestamp = None
    if timestamp is None or timestamp.tzinfo is None:
        raise ValueError(f"item {index} has invalid datetime {text!r}")

    spektraltyp = entry.get("spektraltyp")
    if not isinstance(spektraltyp, str) or not spektraltyp:
        raise ValueError(f"item {index} has invalid spektraltyp {spektraltyp!r}")

    href = entry.get("href")
    if href != f"/files/{collection}/{item_id}.tif":
        raise ValueError(f"item {index} has unexpected href {href!r}")

    resolution = entry.get("resolution")
    if not (_is_finite_number(resolution) and resolution > 0):
        resolution = None
    size_bytes = entry.get("size_bytes")
    if type(size_bytes) is not int or size_bytes <= 0:
        size_bytes = None

    return OrthophotoItem(
        id=item_id,
        collection=collection,
        bbox=tuple(float(v) for v in bbox),
        datetime=timestamp,
        spektraltyp=spektraltyp,
        href=href,
        resolution=None if resolution is None else float(resolution),
        size_bytes=size_bytes,
    )


def _is_finite_number(value) -> bool:
    """True for a finite JSON int or float; bools and huge ints are excluded."""
    if type(value) not in (int, float):
        return False
    try:
        return math.isfinite(value)
    except OverflowError:
        return False
