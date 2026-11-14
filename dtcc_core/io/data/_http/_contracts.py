"""Backend-neutral inputs, results and resource accounting for HTTP operations."""

from __future__ import annotations

import json
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from threading import Lock
from typing import Any, BinaryIO
from urllib.parse import urlsplit, urlunsplit

_UNSET = object()
DEFAULT_BODY_LIMIT = 32 * 1024 * 1024
_Scalar = str | int | float
_Fields = Mapping[str, _Scalar] | Sequence[tuple[str, _Scalar]]
_Headers = Mapping[str, str] | Sequence[tuple[str, str]]


def _pairs(value, name):
    """Snapshot fields as string pairs without guessing compound encodings."""
    items = value.items() if isinstance(value, Mapping) else value
    result = []
    for item in items:
        if not isinstance(item, (tuple, list)) or len(item) != 2:
            raise ValueError(f"{name} must contain name/value pairs")
        key, entry = item
        if not isinstance(key, str):
            raise ValueError(f"{name} names must be strings")
        if not isinstance(entry, str):
            if (
                name == "headers"
                or isinstance(entry, bool)
                or not isinstance(entry, (int, float))
                or (isinstance(entry, float) and not math.isfinite(entry))
            ):
                raise ValueError(f"Invalid {name} value type")
            entry = str(entry)
        result.append((key, entry))
    return tuple(result)


def _positive(value, name, *, zero=False):
    if not math.isfinite(value) or (value < 0 if zero else value <= 0):
        raise ValueError(
            f"{name} must be finite and {'nonnegative' if zero else 'positive'}"
        )


def _target(url):
    """Remove query values, fragments and user information from diagnostics."""
    parsed = urlsplit(url)
    host = parsed.hostname or ""
    if ":" in host:
        host = f"[{host}]"
    if parsed.port:
        host += f":{parsed.port}"
    return urlunsplit((parsed.scheme, host, parsed.path, "", ""))


def _origin(url):
    """Return scheme://host[:port] in lowercase, omitting the scheme's default port."""
    parsed = urlsplit(url)
    scheme, host, port = parsed.scheme.lower(), parsed.hostname or "", parsed.port
    if ":" in host:
        host = f"[{host}]"
    if port is not None and port != {"http": 80, "https": 443}.get(scheme):
        host += f":{port}"
    return f"{scheme}://{host}"


def _allowed_origin(value):
    try:
        parsed = urlsplit(value)
        valid = (
            parsed.scheme.lower() in {"http", "https"}
            and parsed.hostname
            and parsed.path in {"", "/"}
            and not (parsed.query or parsed.fragment or parsed.username)
        )
        origin = _origin(value)
    except (ValueError, TypeError, AttributeError):
        valid = False
    if not valid:
        raise ValueError("allowed_origins must contain scheme://host[:port] origins")
    return origin


def _rebuild_error(cls, message, state):
    error = cls.__new__(cls)
    RuntimeError.__init__(error, message)
    error.__dict__.update(state)
    return error


class HttpTransportError(RuntimeError):
    """Safe transport failure; raw error_body requires caller sanitization."""

    def __init__(
        self,
        request,
        kind,
        *,
        status=None,
        attempts=0,
        retry_allowed_by_policy=False,
        error_body=b"",
        retry_after=None,
        redirect_target=None,
        phase=None,
    ):
        self.method = request.method
        self.target = _target(request.url)
        self.kind = kind
        self.status = status
        self.attempts = attempts
        self.retry_allowed_by_policy = retry_allowed_by_policy
        self.error_body = error_body
        self.retry_after = retry_after
        self.redirect_target = redirect_target
        # "connect" before response headers arrive, "read" after, None without I/O.
        self.phase = phase
        super().__init__(
            f"{self.method} {self.target}: {kind}"
            + (f" (HTTP {status})" if status else "")
        )

    def __reduce__(self):
        # The constructor takes a Request, which the error does not keep.
        return _rebuild_error, (type(self), str(self), dict(self.__dict__))

    def metadata(self):
        """Return diagnostics suitable for recording without raw payloads."""
        return {
            key: getattr(self, key)
            for key in (
                "method",
                "target",
                "kind",
                "status",
                "attempts",
                "retry_allowed_by_policy",
                "phase",
            )
        }


class HttpLimitError(HttpTransportError):
    """A resource limit was crossed; counts and resource name are safe diagnostics."""

    def __init__(self, request, resource, limit, observed, **options):
        self.resource, self.limit, self.observed = resource, limit, observed
        super().__init__(request, "resource_limit", **options)

    def metadata(self):
        return dict(
            super().metadata(),
            resource=self.resource,
            limit=self.limit,
            observed=self.observed,
        )


@dataclass(frozen=True)
class Timeout:
    """Connect/read inactivity timeouts; total is supported only asynchronously."""

    connect: float = 10.0
    read: float = 30.0
    total: float | None = None

    def __post_init__(self):
        _positive(self.connect, "connect")
        _positive(self.read, "read")
        if self.total is not None:
            _positive(self.total, "total")


@dataclass(frozen=True)
class RetryPolicy:
    """Opt-in replay policy; retries counts additional attempts, excluding redirects."""

    retries: int = 0
    statuses: frozenset[int] = frozenset({429, 500, 502, 503, 504})
    connect: bool = True
    read: bool = True
    backoff: float = 0.2
    max_delay: float = 30.0
    jitter: float = 0.0
    respect_retry_after: bool = True

    def __post_init__(self):
        if not isinstance(self.retries, int) or self.retries < 0:
            raise ValueError("retries must be a nonnegative integer")
        for name in ("backoff", "max_delay", "jitter"):
            _positive(getattr(self, name), name, zero=True)
        object.__setattr__(self, "statuses", frozenset(self.statuses))


@dataclass(frozen=True)
class RedirectPolicy:
    """Redirect limits: hop count, cross-origin destinations and a host restriction.

    Same-origin hops need only max_hops. A cross-origin hop also needs its exact
    origin in allowed_origins. allowed_hosts, when set, restricts every request
    and destination to those hosts; it never permits a redirect by itself.
    """

    max_hops: int = 0
    allowed_origins: frozenset[str] = frozenset()
    allowed_hosts: frozenset[str] = frozenset()

    def __post_init__(self):
        if not isinstance(self.max_hops, int) or self.max_hops < 0:
            raise ValueError("max_hops must be a nonnegative integer")
        object.__setattr__(
            self,
            "allowed_origins",
            frozenset(_allowed_origin(value) for value in self.allowed_origins),
        )
        object.__setattr__(
            self,
            "allowed_hosts",
            frozenset(host.lower() for host in self.allowed_hosts),
        )


@dataclass(frozen=True)
class MultipartPart:
    """One named multipart value; binary streams remain owned by the caller."""

    name: str
    value: bytes | BinaryIO = field(repr=False)
    filename: str | None = None
    content_type: str = "application/octet-stream"


@dataclass(frozen=True)
class Request:
    """An explicit HTTP operation, with credentials and payload excluded from repr."""

    method: str
    url: str = field(repr=False)
    params: _Fields = field(default=(), repr=False)
    headers: _Headers = field(default=(), repr=False)
    body: bytes | None = field(default=None, repr=False)
    json: Any = field(default=_UNSET, repr=False)
    form: _Fields | None = field(default=None, repr=False)
    multipart: tuple[MultipartPart, ...] | None = field(default=None, repr=False)
    idempotent: bool = False
    public: bool = False
    cache_ttl: float = 0.0
    cache_namespace: str = field(default="", repr=False)
    cache_token: str = field(default="", repr=False)
    expected_status: frozenset[int] = frozenset(range(200, 300))

    def __post_init__(self):
        method = self.method.upper()
        if not method or not method.isalpha():
            raise ValueError("Invalid HTTP method")
        object.__setattr__(self, "method", method)
        try:
            parsed = urlsplit(self.url)
            valid = (
                parsed.scheme in {"http", "https"}
                and parsed.hostname
                and not parsed.username
                and not parsed.password
            )
            _target(self.url)
        except ValueError:
            valid = False
        if not valid:
            raise ValueError(
                "Expected an HTTP(S) URL without user information"
            ) from None
        bodies = (
            self.body is not None,
            self.json is not _UNSET,
            self.form is not None,
            self.multipart is not None,
        )
        if sum(bodies) > 1:
            raise ValueError("Choose one request body representation")
        if self.body is not None and not isinstance(self.body, bytes):
            raise ValueError("body must be bytes")
        for name in ("params", "headers"):
            object.__setattr__(self, name, _pairs(getattr(self, name) or (), name))
        if self.form is not None:
            object.__setattr__(self, "form", _pairs(self.form, "form"))
        if self.json is not _UNSET:
            # Encode now so every attempt sends the body as it was when built.
            object.__setattr__(
                self,
                "_json_body",
                json.dumps(self.json, allow_nan=False).encode("utf-8"),
            )
        if self.multipart is not None:
            object.__setattr__(self, "multipart", tuple(self.multipart))
        names = [name.lower() for name, _ in self.headers]
        if len(names) != len(set(names)):
            raise ValueError("Request header names must be unique, ignoring case")
        object.__setattr__(self, "expected_status", frozenset(self.expected_status))
        _positive(self.cache_ttl, "cache_ttl", zero=True)

    @property
    def target(self):
        return _target(self.url)

    @property
    def replayable(self):
        return not self.multipart or all(
            isinstance(part.value, bytes) for part in self.multipart
        )


@dataclass(frozen=True)
class Response:
    """A fully buffered response; headers and body are explicit, non-repr values."""

    status: int
    target: str
    headers: dict[str, str] = field(repr=False)
    body: bytes = field(repr=False)
    _request: Request = field(repr=False)

    def json(self):
        try:
            return json.loads(self.body)
        except (ValueError, UnicodeError):
            raise HttpTransportError(
                self._request, "decoding", status=self.status
            ) from None

    def text(self, encoding="utf-8"):
        try:
            return self.body.decode(encoding)
        except (UnicodeError, LookupError):
            raise HttpTransportError(
                self._request, "decoding", status=self.status
            ) from None


@dataclass
class Counters:
    """Client-local counters, including decoded bytes consumed on failed attempts."""

    logical_requests: int = 0
    network_attempts: int = 0
    cache_hits: int = 0
    received_body_bytes: int = 0
    delivered_body_bytes: int = 0


class Budget:
    """Thread-safe shared limits; reservations occur before sends or delivery."""

    def __init__(
        self,
        *,
        max_requests=None,
        max_attempts=None,
        max_received_bytes=None,
        max_delivered_bytes=None,
    ):
        self.limits = dict(
            logical_requests=max_requests,
            network_attempts=max_attempts,
            received_body_bytes=max_received_bytes,
            delivered_body_bytes=max_delivered_bytes,
        )
        for value in self.limits.values():
            if value is not None and (type(value) is not int or value < 0):
                raise ValueError("Budget limits must be nonnegative integers")
        self._used = dict.fromkeys(self.limits, 0)
        self._lock = Lock()
        self._exhausted = None

    def before_read(self, request):
        with self._lock:
            exhausted = self._exhausted
        if exhausted is not None:
            raise HttpLimitError(request, *exhausted)

    def consume(self, name, amount, request):
        with self._lock:
            exhausted = self._exhausted if name == "network_attempts" else None
            if name == "network_attempts" and exhausted is None:
                for key in ("received_body_bytes", "delivered_body_bytes"):
                    cap = self.limits[key]
                    if cap is not None and self._used[key] >= cap:
                        exhausted = (key, cap, self._used[key])
                        break
            if exhausted is None:
                limit = self.limits[name]
                observed = self._used[name] + amount
                if name == "received_body_bytes":
                    # Received bytes count even when the chunk crosses the limit.
                    self._used[name] = observed
                if limit is not None and observed > limit:
                    exhausted = (name, limit, observed)
                    if name in {"received_body_bytes", "delivered_body_bytes"}:
                        self._exhausted = exhausted
                elif name != "received_body_bytes":
                    self._used[name] = observed
        if exhausted is not None:
            raise HttpLimitError(request, *exhausted)
