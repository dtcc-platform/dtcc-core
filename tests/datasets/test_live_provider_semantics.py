"""Regression tests for live provider test failure semantics."""

from __future__ import annotations

from dtcc_core.datasets.dataset import DatasetUpstreamError
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
