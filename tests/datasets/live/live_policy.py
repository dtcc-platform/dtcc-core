"""Failure policy shared by the live dataset tests.

Transient upstream failures (connection errors, timeouts, HTTP 5xx and 429,
also when GDAL reports them for a remote raster)
skip the case: they do not show that a dataset contract broke. Every other
error fails it, including HTTP 4xx, TLS and URL errors, and malformed payloads.
Download helpers often wrap the original error with ``raise ... from``, so the
chain of causes is checked. See ``docs/datasets/qa.md``.
"""

from __future__ import annotations

import re

import aiohttp
import pytest
import requests
from rasterio.errors import RasterioError

from dtcc_core.datasets.dataset import DatasetUpstreamError, redact_url_secrets

OVERPASS_ALL_FAILED = "All Overpass endpoints failed"
_OVERPASS_ERRORS_HEADER = "Errors by endpoint:\n"

# Status codes reported in plain RuntimeError messages, e.g. the LiDAR tile
# lookup ("Request failed with status 503") and Overpass ("HTTP 504").
_STATUS_IN_MESSAGE = re.compile(r"(?:Request failed with status|\bHTTP) (\d{3})\b")

# requests error text for connection failures, timeouts, and retries that ran
# out on 5xx or 429 responses.
_TRANSIENT_REQUEST_TEXT = re.compile(
    r"timed out|Timeout|NewConnectionError|NameResolutionError"
    r"|Connection (?:aborted|refused|reset)|ConnectionResetError|RemoteDisconnected"
    r"|too many (?:5\d\d|429) error responses",
    re.IGNORECASE,
)


# GDAL /vsicurl/ errors, raised by rasterio when it reads a remote file such
# as the GHSL zips: "HTTP response code: 503" or "CURL error: <curl message>".
_GDAL_HTTP_STATUS = re.compile(r"HTTP response code: (\d{3})\b")
_GDAL_TRANSIENT_CURL = re.compile(
    r"CURL error: .*(?:timed out|Timeout was reached|Could not resolve host"
    r"|Couldn't connect|Failed to connect|Connection reset|Recv failure"
    r"|Empty reply)",
    re.IGNORECASE,
)

try:
    from rasterio._err import CPLE_BaseError
except ImportError:  # pragma: no cover - depends on the rasterio build
    _GDAL_ERRORS: tuple[type[BaseException], ...] = (RasterioError,)
else:
    _GDAL_ERRORS = (RasterioError, CPLE_BaseError)


def _is_transient_status(status: int | None) -> bool:
    return status is not None and (status == 429 or 500 <= status <= 599)


def _overpass_failure_is_transient(message: str) -> bool:
    """Classify an "All Overpass endpoints failed" message.

    When the message lists each endpoint's error, the failure is transient only
    if every endpoint failed transiently, so one mirror rejecting the query is
    not hidden by another timing out. Otherwise only the last error is known.
    """
    _, listed, by_endpoint = message.partition(_OVERPASS_ERRORS_HEADER)
    if listed:
        problems = [
            line.split(": ", 1)[-1]
            for line in by_endpoint.splitlines()
            if line.startswith("- ")
        ]
        if problems:
            return all(_overpass_error_is_transient(p) for p in problems)
    return _overpass_error_is_transient(message)


def _overpass_error_is_transient(message: str) -> bool:
    """Classify one Overpass endpoint error."""
    if "invalid JSON" in message:
        return False
    match = _STATUS_IN_MESSAGE.search(message)
    if match:
        return _is_transient_status(int(match.group(1)))
    # Otherwise the last error is a requests exception message.
    return bool(_TRANSIENT_REQUEST_TEXT.search(message))


def _is_transient_single(exc: BaseException) -> bool | None:
    """Classify one exception; None means it says nothing on its own."""
    if isinstance(exc, DatasetUpstreamError):
        return exc.is_transient or exc.status_code == 429
    if isinstance(exc, (requests.exceptions.SSLError, aiohttp.ClientSSLError)):
        return False
    if isinstance(exc, requests.HTTPError):
        return _is_transient_status(getattr(exc.response, "status_code", None))
    if isinstance(
        exc,
        (
            requests.ConnectionError,
            requests.Timeout,
            aiohttp.ClientConnectionError,
            TimeoutError,
        ),
    ):
        return True
    if isinstance(exc, aiohttp.ClientResponseError):
        return _is_transient_status(exc.status)
    if isinstance(exc, _GDAL_ERRORS):
        return _gdal_error_is_transient(str(exc))
    status = getattr(exc, "status", None)
    if isinstance(status, int):
        return _is_transient_status(status)
    if isinstance(exc, RuntimeError):
        message = str(exc)
        if message.startswith(OVERPASS_ALL_FAILED):
            return _overpass_failure_is_transient(message)
        match = _STATUS_IN_MESSAGE.search(message)
        if match:
            return _is_transient_status(int(match.group(1)))
    return None


def _gdal_error_is_transient(message: str) -> bool:
    """Classify a GDAL remote-read error by its HTTP status or curl message."""
    match = _GDAL_HTTP_STATUS.search(message)
    if match:
        return _is_transient_status(int(match.group(1)))
    return bool(_GDAL_TRANSIENT_CURL.search(message))


def is_transient_live_failure(exc: BaseException) -> bool:
    """Return whether a live-case exception is a transient upstream failure."""
    seen: set[int] = set()
    current: BaseException | None = exc
    while current is not None and id(current) not in seen:
        seen.add(id(current))
        verdict = _is_transient_single(current)
        if verdict is not None:
            return verdict
        current = current.__cause__
    return False


def skip_if_transient(exc: BaseException, case: str) -> None:
    """Skip the running test when ``exc`` is a transient upstream failure."""
    if is_transient_live_failure(exc):
        pytest.skip(transient_skip_reason(exc, case))


def transient_skip_reason(exc: BaseException, case: str) -> str:
    """Return a skip message for ``exc`` with credential query values removed."""
    return f"{case}: transient upstream failure: {redact_url_secrets(str(exc))}"
