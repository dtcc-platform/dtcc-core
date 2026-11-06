"""Failure policy shared by the live dataset tests.

Transient upstream failures (connection errors, timeouts, HTTP 5xx and 429)
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

from dtcc_core.datasets.dataset import DatasetUpstreamError, _redact_url_secrets

OVERPASS_ALL_FAILED = "All Overpass endpoints failed"

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


def _is_transient_status(status: int | None) -> bool:
    return status is not None and (status == 429 or 500 <= status <= 599)


def _overpass_last_error_is_transient(message: str) -> bool:
    """Classify the last mirror error in an "All Overpass endpoints failed" message."""
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
    status = getattr(exc, "status", None)
    if isinstance(status, int):
        return _is_transient_status(status)
    if isinstance(exc, RuntimeError):
        message = str(exc)
        if message.startswith(OVERPASS_ALL_FAILED):
            return _overpass_last_error_is_transient(message)
        match = _STATUS_IN_MESSAGE.search(message)
        if match:
            return _is_transient_status(int(match.group(1)))
    return None


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
    return f"{case}: transient upstream failure: {_redact_url_secrets(str(exc))}"
