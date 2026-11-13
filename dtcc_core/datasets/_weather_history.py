"""Shared transport and normalization for historical weather providers."""

import json
import time
from collections import OrderedDict
from contextlib import contextmanager
from threading import Lock

import numpy as np

from ..model import Point, TimeSeries
from ..model.values.time_series import _to_times
from ..reproject.reproject import reproject_array
from .dataset import DatasetDescriptor, DatasetUpstreamError

_CACHE = OrderedDict()
_CACHE_LOCK = Lock()
_CACHE_BYTES = 32 * 1024 * 1024
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
    """Serial HTTP requests with retries, bounded reads and a bounded TTL cache."""

    def __init__(self, args):
        import requests
        from requests.adapters import HTTPAdapter
        from urllib3.util.retry import Retry

        self.args = args
        self.requests = 0
        self.bytes = 0
        self.observations = 0
        self.session = requests.Session()
        retry = Retry(
            total=2,
            backoff_factor=0.2,
            status_forcelist=(429, 500, 502, 503, 504),
            allowed_methods=("GET",),
            respect_retry_after_header=False,
        )
        self.session.mount("https://", HTTPAdapter(max_retries=retry))

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.session.close()

    def _read(self, url, params=None, cache_token=None):
        import requests

        self.requests += 1
        if self.requests > self.args.max_requests:
            raise WeatherHistoryLimitError("Historical weather exceeded max_requests")
        key = (url, tuple(sorted((params or {}).items())), cache_token)
        cacheable = self.args.cache_ttl_s > 0 and "apikey" not in (params or {})
        body = None
        if cacheable:
            with _CACHE_LOCK:
                cached = _CACHE.get(key)
                if cached and cached[0] > time.monotonic():
                    body = cached[1]
                    _CACHE.move_to_end(key)
        if body is None:
            try:
                with self.session.get(
                    url,
                    params=params,
                    timeout=self.args.timeout_s,
                    stream=True,
                ) as response:
                    response.raise_for_status()
                    chunks, size = [], 0
                    for chunk in response.iter_content(chunk_size=65536):
                        size += len(chunk)
                        if self.bytes + size > self.args.max_download_bytes:
                            raise WeatherHistoryLimitError(
                                "Historical weather exceeded max_download_bytes"
                            )
                        chunks.append(chunk)
                    body = b"".join(chunks)
                    cacheable = cacheable and "no-store" not in response.headers.get(
                        "Cache-Control", ""
                    )
            except requests.RequestException as exc:
                # Request exceptions may embed credential-bearing query strings.
                error = DatasetDescriptor.build_upstream_error(
                    "weather_history", "fetch", url, exc
                )
                error.message = (
                    f"weather_history fetch failed for {url}: {error.failure_class}"
                )
                error.args = (error.message,)
                raise error from None
            if cacheable and len(body) <= _CACHE_BYTES:
                with _CACHE_LOCK:
                    _CACHE[key] = (time.monotonic() + self.args.cache_ttl_s, body)
                    _CACHE.move_to_end(key)
                    while (
                        len(_CACHE) > 32
                        or sum(len(item[1]) for item in _CACHE.values()) > _CACHE_BYTES
                    ):
                        _CACHE.popitem(last=False)
        self.bytes += len(body)
        if self.bytes > self.args.max_download_bytes:
            raise WeatherHistoryLimitError(
                "Historical weather exceeded max_download_bytes"
            )
        return body

    def json(self, url, params=None):
        with _guard(url):
            value = json.loads(self._read(url, params))
            if not isinstance(value, dict) or value.get("error"):
                raise ValueError("Expected a successful JSON object")
            return value

    def text(self, url, cache_token=None):
        with _guard(url):
            return self._read(url, cache_token=cache_token).decode("utf-8-sig")

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
