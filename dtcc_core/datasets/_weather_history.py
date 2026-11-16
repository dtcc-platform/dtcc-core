"""Shared transport and normalization for historical weather providers."""

from contextlib import contextmanager

import numpy as np

from ..io import HttpLimitError, HttpTransportError, ResponseCache, http_client
from ..model import Point, TimeSeries
from ..model.values.time_series import _to_times
from ..reproject.reproject import reproject_array
from .dataset import DatasetDescriptor, DatasetUpstreamError

_CACHE = ResponseCache()
_HOUR = 3_600_000_000_000


class NoMatchingStation(LookupError):
    """No station offers all requested measurements within the distance limit."""


class WeatherHistoryLimitError(ValueError):
    """A historical-weather request exceeds an explicit resource budget."""


def _invalid(target, message):
    return DatasetUpstreamError(
        dataset="weather_history",
        operation="parse",
        target=target,
        failure_class="invalid_payload",
        message=f"Invalid weather history from {target}: {message}",
    )


@contextmanager
def _guard(target):
    try:
        yield
    except (DatasetUpstreamError, WeatherHistoryLimitError):
        raise
    except HttpLimitError as exc:
        limit = (
            "max_requests"
            if exc.resource in {"logical_requests", "network_attempts"}
            else "max_download_bytes"
        )
        raise WeatherHistoryLimitError(f"Historical weather exceeded {limit}") from None
    except HttpTransportError as exc:
        if exc.kind == "decoding":
            raise _invalid(exc.target, "Cannot decode response") from None
        failure_class = exc.kind
        if exc.kind == "http_status":
            if 500 <= exc.status <= 599:
                failure_class = "http_5xx"
            elif 400 <= exc.status <= 499:
                failure_class = "http_4xx"
            else:
                failure_class = "http_error"
        elif exc.kind == "certificate":
            failure_class = "connection"
        raise DatasetUpstreamError(
            dataset="weather_history",
            operation="fetch",
            target=exc.target,
            failure_class=failure_class,
            status_code=exc.status,
            message=f"weather_history fetch failed for {exc.target}: {failure_class}",
        ) from None
    except (
        ValueError,
        TypeError,
        KeyError,
        IndexError,
        StopIteration,
        OverflowError,
    ) as exc:
        raise _invalid(target, type(exc).__name__) from exc


class _Client:
    """Weather-specific payload validation and observation counting."""

    def __init__(self, args):
        self.args = args
        self.observations = 0
        self.http = http_client(
            timeout=args.timeout_s,
            retries=2,
            respect_retry_after=False,
            max_requests=args.max_requests,
            max_attempts=args.max_requests,
            max_download_bytes=args.max_download_bytes,
            cache_ttl=args.cache_ttl_s,
            cache=_CACHE,
            cache_namespace="weather_history",
        )

    @property
    def requests(self):
        return self.http.counters.logical_requests

    @property
    def bytes(self):
        return self.http.counters.delivered_body_bytes

    def __enter__(self):
        self.http.__enter__()
        return self

    def __exit__(self, *exc):
        self.http.__exit__(*exc)

    def json(self, url, params=None):
        with _guard(url):
            query = dict(params or {})
            credentials = {"apikey": query.pop("apikey")} if "apikey" in query else None
            value = self.http.json(url, params=query, credentials=credentials)
            if not isinstance(value, dict) or value.get("error"):
                raise ValueError("Expected a successful JSON object")
            return value

    def text(self, url, cache_token=None):
        with _guard(url):
            return self.http.text(url, cache_token=cache_token, encoding="utf-8-sig")

    def count_observation(self):
        self.observations += 1
        if self.observations > self.args.max_observations:
            raise WeatherHistoryLimitError(
                "Historical weather exceeded max_observations"
            )


def _stamp(value):
    times, scalar = _to_times(value)
    if not scalar:
        raise ValueError("Expected a scalar time")
    return times[0]


def _ns(value):
    return int(_stamp(value).astype(np.int64))


def _link(resource, rel, media_type):
    from urllib.parse import urlparse

    href = next(
        item["href"]
        for item in resource["link"]
        if item["rel"] == rel and item["type"] == media_type
    )
    # SMHI publishes absolute links on its metobs download host.
    parsed = urlparse(href)
    if parsed.scheme != "https" or parsed.netloc != "opendata-download-metobs.smhi.se":
        raise ValueError("Unexpected SMHI resource host")
    return href


def _output_point(lon, lat, height, crs):
    coords = reproject_array(
        np.array([[lon, lat, height]], dtype=float), "EPSG:4326", crs
    )[0]
    point = Point(x=float(coords[0]), y=float(coords[1]), z=float(coords[2]))
    point.transform.srs = crs
    return point


def _store(result, key, times, values, unit, role, metadata, observations, args):
    series = TimeSeries(
        name=key[1],
        unit=unit,
        times=times,
        values=np.asarray(values, dtype=float),
        association="sample",
        role=role,
        description=metadata.get("description", ""),
    ).between(args.start, args.end)
    start, end = _ns(args.start), _ns(args.end)
    expected = max(0, -(-end // _HOUR) - -(-start // _HOUR))
    finite = int(np.isfinite(series.values).sum())
    metadata.update(
        location_id=key[0],
        requested_start=args.start,
        requested_end_exclusive=args.end,
        returned_start=str(series.start) if len(series) else None,
        returned_end=str(series.end) if len(series) else None,
        observation_count=len(series),
        missing_value_count=int(np.isnan(series.values).sum()),
        expected_hourly_count=expected,
        coverage_complete=(
            finite == expected
            and len(series) == expected
            and bool((series.times.astype(np.int64) % _HOUR == 0).all())
        ),
    )
    result.series[key] = series
    result.series_metadata[key] = metadata
    result.observations[key] = {
        str(stamp): observations[str(stamp)]
        for stamp in series.times
        if str(stamp) in observations
    }


def _finish_health(result, errors, parameters):
    DatasetDescriptor.apply_result_health_metadata(
        result.attributes,
        upstream_errors=errors,
        requested_parameters=parameters,
        fetched_parameters=[key[1] for key in result],
    )
    result.attributes["coverage_complete"] = (
        bool(result)
        and all(
            metadata["coverage_complete"]
            for metadata in result.series_metadata.values()
        )
        and len(result) == len(parameters)
    )
    if result and not result.attributes["coverage_complete"]:
        result.attributes["partial_result"] = True
