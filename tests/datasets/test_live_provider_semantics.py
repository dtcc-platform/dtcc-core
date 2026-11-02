"""Regression tests for live provider test failure semantics."""

from __future__ import annotations

import pytest
import requests

from dtcc_core.datasets.dataset import DatasetDescriptor, DatasetUpstreamError
from tests.datasets.live.live_policy import (
    is_transient_live_failure,
    skip_if_transient,
)
from tests.datasets.live.test_live_sweden import _case_result_from_exception


def _upstream_error(failure_class: str) -> DatasetUpstreamError:
    return DatasetUpstreamError(
        dataset="weather",
        operation="fetch",
        target="fixture://live-case",
        failure_class=failure_class,
        message=f"weather fetch failed with {failure_class}",
    )


def test_transient_live_upstream_errors_may_skip() -> None:
    for failure_class in ("connection", "timeout", "http_5xx"):
        case = _case_result_from_exception(_upstream_error(failure_class))

        assert case.status == "skip"


def test_non_transient_live_upstream_errors_fail() -> None:
    for failure_class in ("http_4xx", "invalid_payload", "request"):
        case = _case_result_from_exception(_upstream_error(failure_class))

        assert case.status == "fail"


def test_untyped_live_exceptions_fail() -> None:
    case = _case_result_from_exception(ValueError("malformed live payload"))

    assert case.status == "fail"


def _wrapped(outer: BaseException, cause: BaseException) -> BaseException:
    outer.__cause__ = cause
    return outer


def _context_only(outer: BaseException, context: BaseException) -> BaseException:
    outer.__context__ = context
    return outer


def test_skip_message_redacts_credential_query_params() -> None:
    exc = requests.ConnectionError("Max retries exceeded with url: /feed?key=s3cr3t")

    with pytest.raises(pytest.skip.Exception) as skipped:
        skip_if_transient(exc, "transit")

    assert "s3cr3t" not in str(skipped.value)
    assert "key=REDACTED" in str(skipped.value)


def _http_error(
    status: int, url: str = "https://example.test/feed"
) -> requests.HTTPError:
    response = requests.Response()
    response.status_code = status
    response.url = url
    return requests.HTTPError(f"{status} Error: for url: {url}", response=response)


@pytest.mark.parametrize(
    "exc",
    [
        requests.ConnectionError("connection refused"),
        requests.Timeout("read timed out"),
        _http_error(503),
        _http_error(429),
        RuntimeError("All Overpass endpoints failed: HTTP 504"),
        RuntimeError(
            "All Overpass endpoints failed. Last error: HTTPSConnectionPool"
            "(host='overpass.example', port=443): Read timed out."
        ),
        RuntimeError("Request failed with status 503:\nService Unavailable"),
        _wrapped(
            RuntimeError("Failed to download DeSO 2025 polygons"), _http_error(502)
        ),
        _wrapped(
            RuntimeError("Failed to query SCB DeSO statistic"),
            requests.Timeout("read timed out"),
        ),
        TimeoutError(),
        DatasetUpstreamError(
            dataset="transit_vehicles",
            operation="fetch",
            target="fixture://live-case",
            failure_class="http_4xx",
            status_code=429,
            message="rate limited",
        ),
    ],
    ids=[
        "connection",
        "timeout",
        "http_503",
        "http_429",
        "overpass_all_http_504",
        "overpass_all_timed_out",
        "lidar_lookup_503",
        "wrapped_http_502",
        "wrapped_timeout",
        "asyncio_timeout",
        "typed_http_429",
    ],
)
def test_transient_untyped_failures_skip(exc: BaseException) -> None:
    assert is_transient_live_failure(exc)
    assert _case_result_from_exception(exc).status == "skip"


@pytest.mark.parametrize(
    "exc",
    [
        _http_error(404),
        RuntimeError("Lidar tile lookup failed: bad payload"),
        ValueError("malformed live payload"),
        RuntimeError("All Overpass endpoints failed: HTTP 400"),
        RuntimeError(
            "All Overpass endpoints failed: invalid JSON from https://o.example: x"
        ),
        RuntimeError("Request failed with status 404:\nNot Found"),
        _wrapped(
            RuntimeError("Failed to download DeSO 2025 polygons"), _http_error(400)
        ),
        requests.exceptions.SSLError("certificate verify failed"),
        RuntimeError(
            "All Overpass endpoints failed: HTTPSConnectionPool(host='o.example', "
            "port=443): Max retries exceeded with url: /api/interpreter (Caused by "
            "SSLError(SSLCertVerificationError('certificate has expired')))"
        ),
        RuntimeError("All Overpass endpoints failed: Exceeded 30 redirects."),
        _context_only(KeyError("stations"), requests.Timeout("read timed out")),
    ],
    ids=[
        "http_404",
        "other_runtime_error",
        "value_error",
        "overpass_all_http_400",
        "overpass_all_invalid_json",
        "lidar_lookup_404",
        "wrapped_http_400",
        "ssl_error",
        "overpass_all_ssl_error",
        "overpass_all_redirects",
        "fallback_error_while_handling_timeout",
    ],
)
def test_other_failures_fail(exc: BaseException) -> None:
    assert not is_transient_live_failure(exc)


def test_upstream_error_message_redacts_credential_query_params() -> None:
    url = "https://opendata.samtrafiken.se/gtfs-rt/sl/VehiclePositions.pb"
    errors = [
        _http_error(403, f"{url}?key=s3cr3t-value&format=pb"),
        requests.ConnectionError(
            "HTTPSConnectionPool(host='opendata.samtrafiken.se', port=443): Max "
            "retries exceeded with url: /gtfs-rt/sl/VehiclePositions.pb"
            "?key=s3cr3t-value (Caused by NameResolutionError)"
        ),
        requests.HTTPError(f"401 for url: {url}?format=pb&api_key=s3cr3t-value"),
        requests.HTTPError(f"401 for url: {url}?access_token=s3cr3t-value"),
        requests.HTTPError(f"401 for url: {url}?subscription-key=s3cr3t-value"),
        requests.HTTPError(f"403 for url: {url}?X-Amz-Signature=s3cr3t-value"),
        requests.HTTPError(f"401 for url: {url}?a=1&db_password=s3cr3t-value"),
        requests.HTTPError(f"401 for url: {url}?client_auth=s3cr3t-value"),
        requests.HTTPError(f"401 for url: {url}?a=1&x-api-key=s3cr3t-value"),
    ]
    for exc in errors:
        error = DatasetDescriptor.build_upstream_error(
            "transit_vehicles", "fetch_gtfs_realtime", url, exc
        )

        assert "s3cr3t-value" not in str(error)
        assert "REDACTED" in str(error)
        assert url in str(error)


def test_upstream_error_message_keeps_ordinary_query_params() -> None:
    url = "https://example.test/wfs?bbox=1,2,3,4&typeName=deso&authority=EPSG&author=x"
    error = DatasetDescriptor.build_upstream_error(
        "deso", "fetch", url, requests.HTTPError(f"500 for url: {url}")
    )

    assert "bbox=1,2,3,4&typeName=deso&authority=EPSG&author=x" in str(error)
