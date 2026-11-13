"""Explicit, bounded public-response cache; no persistent or implicit global state."""

import hashlib
import json
import time
from collections import OrderedDict
from threading import Lock

from ._contracts import _UNSET
from ._policy import headers_lower, http_date_timestamp


def _directives(headers):
    value = headers_lower(headers).get("cache-control", "")
    return {
        part.strip().split("=", 1)[0].lower(): part.strip().partition("=")[2].strip('"')
        for part in value.split(",")
    }


def _eligible(request):
    headers = headers_lower(dict(request.headers))
    return (
        request.public
        and request.body is None
        and request.json is _UNSET
        and request.form is None
        and request.multipart is None
        and request.method == "GET"
        and request.cache_ttl > 0
        and not any(
            key in headers for key in ("authorization", "cookie", "proxy-authorization")
        )
        and not {"no-store", "no-cache", "max-age"} & _directives(headers).keys()
        and "pragma" not in headers
    )


def _key(request):
    identity = (
        request.method,
        request.url,
        request.params,
        sorted((key.lower(), value) for key, value in request.headers),
        request.cache_namespace,
        request.cache_token,
    )
    return hashlib.sha256(json.dumps(identity).encode()).digest()


def _lifetime(headers, directives, ttl, wall_clock):
    """Cap local TTL by origin freshness, including apparent age from Date."""
    try:
        now = wall_clock()
        date = http_date_timestamp(headers["date"]) if "date" in headers else now
        age = max(0, int(headers.get("age", "0")), now - date)
        if "max-age" in directives:
            return min(ttl, max(0, int(directives["max-age"]) - age))
        if "expires" in headers:
            fresh = http_date_timestamp(headers["expires"]) - date
            return min(ttl, max(0, fresh - age))
        return ttl
    except (ValueError, TypeError, OverflowError):
        return 0


class ResponseCache:
    """TTL/LRU cache for explicit public buffered GETs, safe to share across clients."""

    def __init__(
        self,
        *,
        max_entries=32,
        max_bytes=32 * 1024 * 1024,
        clock=time.monotonic,
        wall_clock=time.time,
    ):
        if (
            not isinstance(max_entries, int)
            or not isinstance(max_bytes, int)
            or min(max_entries, max_bytes) < 0
        ):
            raise ValueError("Cache limits must be nonnegative integers")
        self.max_entries, self.max_bytes = max_entries, max_bytes
        self._clock, self._wall_clock, self._lock = clock, wall_clock, Lock()
        self._entries = OrderedDict()
        self._bytes = 0

    def clear(self):
        with self._lock:
            self._entries.clear()
            self._bytes = 0

    def invalidate(self, request):
        with self._lock:
            self._remove(_key(request))

    def _remove(self, key):
        entry = self._entries.pop(key, None)
        if entry:
            self._bytes -= len(entry[2].body)

    def get(self, request):
        if not _eligible(request):
            return None
        key = _key(request)
        with self._lock:
            entry = self._entries.get(key)
            if entry is None:
                return None
            created, expiry, response = entry
            if self._clock() >= min(expiry, created + request.cache_ttl):
                self._remove(key)
                return None
            self._entries.move_to_end(key)
            # Headers are mutable convenience values: don't share the stored dict.
            from dataclasses import replace

            return replace(response, headers=dict(response.headers), _request=request)

    def put(self, request, response):
        headers = headers_lower(response.headers)
        directives = _directives(headers)
        if (
            not _eligible(request)
            or response.status != 200
            or {"no-store", "private", "no-cache"} & directives.keys()
            or "set-cookie" in headers
            or "vary" in headers
            or len(response.body) > self.max_bytes
            or not self.max_entries
        ):
            return
        ttl = _lifetime(headers, directives, request.cache_ttl, self._wall_clock)
        from dataclasses import replace

        now, key = self._clock(), _key(request)
        with self._lock:
            self._remove(key)
            # A response that is already stale replaces no fresh entry.
            if ttl <= 0:
                return
            self._entries[key] = (
                now,
                now + ttl,
                replace(response, headers=dict(response.headers)),
            )
            self._bytes += len(response.body)
            while len(self._entries) > self.max_entries or self._bytes > self.max_bytes:
                self._remove(next(iter(self._entries)))
