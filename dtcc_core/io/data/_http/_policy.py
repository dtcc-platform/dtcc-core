"""Decisions shared by synchronous and asynchronous HTTP clients."""

import math
import random
import time
from dataclasses import replace
from datetime import UTC
from email.utils import parsedate_to_datetime
from urllib.parse import urljoin, urlsplit

from ._contracts import (
    _UNSET,
    Budget,
    Counters,
    HttpLimitError,
    HttpTransportError,
    RedirectPolicy,
    RetryPolicy,
    Timeout,
    _origin,
    _target,
)


def headers_lower(headers):
    result = {}
    for key, value in headers.items():
        key = key.lower()
        result[key] = result[key] + ", " + value if key in result else value
    return result


def redirect(request, location, policy, hops):
    allowed = (
        policy.max_hops > hops
        and request.method in {"GET", "HEAD"}
        and request.body is None
        and request.json is _UNSET
        and request.form is None
        and request.multipart is None
    )
    # Request validates schemes, user information and malformed ports.
    try:
        target = urljoin(request.url, location)
        updated = replace(request, url=target, params=())
        parsed = urlsplit(target)
    except ValueError:
        raise HttpTransportError(request, "redirect_policy") from None
    destination = _origin(target)
    if policy.allowed_hosts and parsed.hostname not in policy.allowed_hosts:
        allowed = False
    if parsed.scheme == "http" and urlsplit(request.url).scheme == "https":
        allowed = False
    headers = headers_lower(dict(request.headers))
    if _origin(request.url) != destination:
        explicit_origin = destination in policy.allowed_origins
        allowed = allowed and explicit_origin and request.public
        # Forward only public representation headers across origins.
        updated = replace(
            updated,
            headers=tuple(
                (key, value)
                for key, value in request.headers
                if key.lower() in {"accept", "accept-encoding", "user-agent"}
            ),
        )
        if "authorization" in headers or "cookie" in headers:
            allowed = False
    if not allowed:
        raise HttpTransportError(request, "redirect_policy")
    return updated


def replay_allowed(request):
    return request.replayable and (
        request.method in {"GET", "HEAD"} or request.idempotent
    )


def retry_allowed(request, policy, used, kind, phase, status=None):
    if used >= policy.retries or not replay_allowed(request):
        return False
    if kind == "http_status":
        return status in policy.statuses
    if kind in {"timeout", "connection"}:
        return policy.connect if phase == "connect" else policy.read
    return False


def http_date_timestamp(value):
    """Return the POSIX time of an HTTP date; a "-0000" zone is read as UTC."""
    date = parsedate_to_datetime(value)
    if date.tzinfo is None:
        date = date.replace(tzinfo=UTC)
    return date.timestamp()


def parse_retry_after(value, wall_clock):
    """Return uncapped seconds, or None for an absent/invalid Retry-After."""
    if not value:
        return None
    value = value.strip()
    try:
        if value.isascii() and value.isdecimal():
            delay = float(value)
        else:
            delay = http_date_timestamp(value) - wall_clock()
        return max(0.0, delay) if math.isfinite(delay) else None
    except (ValueError, TypeError, OverflowError):
        return None


def retry_delay(policy, used, headers, wall_clock, jitter_source):
    delay = min(policy.max_delay, policy.backoff * 2 ** min(used, 30))
    if policy.respect_retry_after:
        requested = parse_retry_after(
            headers_lower(headers).get("retry-after"), wall_clock
        )
        if requested is not None:
            delay = requested
    return min(policy.max_delay, delay + policy.jitter * jitter_source())


class Operation:
    """Per-operation accounting and decisions, independent of I/O scheduling."""

    def __init__(self, state, request, limit):
        if limit is not None and (not isinstance(limit, int) or limit < 0):
            raise ValueError("max_body_bytes must be a nonnegative integer")
        # Reject only requests the retry policy would otherwise replay.
        if (
            state.retry.retries
            and not request.replayable
            and (request.method in {"GET", "HEAD"} or request.idempotent)
        ):
            raise ValueError("Automatic retries cannot replay upload streams")
        if (
            state.redirects.allowed_hosts
            and urlsplit(request.url).hostname not in state.redirects.allowed_hosts
        ):
            raise ValueError("Request host is not in allowed_hosts")
        self.state, self.request, self.limit = state, request, limit
        self.attempts = self.retries = self.hops = self.size = 0
        self.local_errors = []
        state.counters.logical_requests += 1
        self._consume("logical_requests", 1, None)

    def _consume(self, name, amount, phase):
        """Charge the shared budget, reporting this operation's attempts on failure."""
        try:
            self.state.budget.consume(name, amount, self.request)
        except HttpLimitError as exc:
            raise HttpLimitError(
                self.request,
                exc.resource,
                exc.limit,
                exc.observed,
                attempts=self.attempts,
                phase=phase,
            ) from None

    def attempt(self):
        self._consume("network_attempts", 1, "connect")
        self.state.counters.network_attempts += 1
        self.attempts += 1
        self.size = 0

    def before_read(self):
        try:
            self.state.budget.before_read(self.request)
        except HttpLimitError as exc:
            raise HttpLimitError(
                self.request,
                exc.resource,
                exc.limit,
                exc.observed,
                attempts=self.attempts,
                phase="read",
            ) from None

    def received(self, chunk, *, limited=True):
        """Account for a chunk; limited=False skips max_body_bytes, not budgets."""
        self.size += len(chunk)
        self.state.counters.received_body_bytes += len(chunk)
        self._consume("received_body_bytes", len(chunk), "read")
        if limited and self.limit is not None and self.size > self.limit:
            raise HttpLimitError(
                self.request,
                "body_bytes",
                self.limit,
                self.size,
                attempts=self.attempts,
                phase="read",
            )

    def delivered(self, size, *, phase="read"):
        self._consume("delivered_body_bytes", size, phase)
        self.state.counters.delivered_body_bytes += size

    def error(self, kind, *, status=None, body=b"", phase=None, redirect_target=None):
        retry_allowed_by_policy = retry_allowed(
            self.request, self.state.retry, 0, kind, phase, status
        )
        return HttpTransportError(
            self.request,
            kind,
            status=status,
            attempts=self.attempts,
            retry_allowed_by_policy=retry_allowed_by_policy,
            error_body=body,
            redirect_target=redirect_target,
            phase=phase,
        )

    def retry(self, error, phase):
        return retry_allowed(
            self.request,
            self.state.retry,
            self.retries,
            error.kind,
            error.phase or phase,
            error.status,
        )

    def delay(self, headers):
        delay = retry_delay(
            self.state.retry,
            self.retries,
            headers,
            self.state.wall_clock,
            self.state.jitter_source,
        )
        self.retries += 1
        return delay

    def follow(self, location, *, status):
        try:
            self.request = redirect(
                self.request,
                location,
                self.state.redirects,
                self.hops,
            )
        except HttpTransportError:
            target = None
            try:
                url = urljoin(self.request.url, location)
                parsed = urlsplit(url)
                if parsed.scheme in {"http", "https"} and parsed.hostname:
                    target = _target(url)
            except ValueError:
                pass
            raise self.error(
                "redirect_policy", status=status, redirect_target=target, phase="read"
            ) from None
        self.hops += 1


class State:
    """Explicit client state; concurrency is scheduled by the caller."""

    def __init__(
        self,
        *,
        timeout=None,
        retry=None,
        redirects=None,
        budget=None,
        cache=None,
        wall_clock=time.time,
        jitter_source=random.random,
    ):
        self.timeout = timeout or Timeout()
        self.retry = retry or RetryPolicy()
        self.redirects = redirects or RedirectPolicy()
        self.budget = budget or Budget()
        self.cache = cache
        self.wall_clock = wall_clock
        self.jitter_source = jitter_source
        self.counters = Counters()

    def cached(self, operation):
        if self.cache is None:
            return None
        response = self.cache.get(operation.request)
        if response is not None:
            if response.status not in operation.request.expected_status:
                return None
            if operation.limit is not None and len(response.body) > operation.limit:
                raise HttpLimitError(
                    operation.request, "body_bytes", operation.limit, len(response.body)
                )
            # A cache hit performs no network I/O, so its errors have no phase.
            operation.delivered(len(response.body), phase=None)
            self.counters.cache_hits += 1
        return response

    def store(self, request, response):
        if self.cache is not None:
            self.cache.put(request, response)
