"""Synchronous transport over requests, with explicit retries and ownership."""

import json
import os
import random
import ssl
import time
from collections.abc import Callable
from contextlib import contextmanager

from ._cache import ResponseCache
from ._contracts import (
    _UNSET,
    DEFAULT_BODY_LIMIT,
    Budget,
    HttpTransportError,
    RedirectPolicy,
    Request,
    Response,
    RetryPolicy,
    Timeout,
)
from ._policy import Operation, State, headers_lower, parse_retry_after


def _caused_by(exc, error_type):
    """Search an exception's causes, contexts, arguments and urllib3 reasons."""
    pending, seen = [exc], set()
    while pending:
        current = pending.pop()
        if not isinstance(current, BaseException) or id(current) in seen:
            continue
        if isinstance(current, error_type):
            return True
        seen.add(id(current))
        pending += [current.__cause__, current.__context__, *current.args]
        pending.append(getattr(current, "reason", None))
    return False


class _Stream:
    def __init__(self, client, operation, response):
        self._client, self._operation, self._response = client, operation, response
        self.status = response.status_code
        self.headers = headers_lower(response.headers)
        self.target = operation.request.target
        self._started = False

    def iter_bytes(self, chunk_size=65536):
        """Deliver decoded chunks once; caller exceptions stay outside transport."""
        if chunk_size <= 0:
            raise ValueError("chunk_size must be positive")
        if self._started:
            raise ValueError("Response stream has already been consumed")
        self._started = True
        for chunk in self._client._chunks(self._response, self._operation, chunk_size):
            self._operation.delivered(len(chunk))
            yield chunk


class HttpClient:
    """Context-managed synchronous HTTP client; not shared between threads."""

    def __init__(
        self,
        *,
        timeout: Timeout | None = None,
        retry: RetryPolicy | None = None,
        redirects: RedirectPolicy | None = None,
        budget: Budget | None = None,
        cache: ResponseCache | None = None,
        ca_bundle: str | os.PathLike[str] | None = None,
        proxy: str | None = None,
        sleep: Callable[[float], None] = time.sleep,
        wall_clock: Callable[[], float] = time.time,
        jitter_source: Callable[[], float] = random.random,
    ):
        import requests

        self._requests = requests
        self._state = State(
            timeout=timeout,
            retry=retry,
            redirects=redirects,
            budget=budget,
            cache=cache,
            wall_clock=wall_clock,
            jitter_source=jitter_source,
        )
        if self._state.timeout.total is not None:
            raise ValueError("Synchronous total deadlines are not supported")
        if ca_bundle is not None:
            ca_bundle = os.fspath(ca_bundle)
        self._session = requests.Session()
        self._session.trust_env = False
        self._session.verify = (
            ca_bundle
            or os.environ.get("REQUESTS_CA_BUNDLE")
            or os.environ.get("CURL_CA_BUNDLE")
            or True
        )
        self._proxies = {"http": proxy, "https": proxy} if proxy else {}
        self._sleep, self._closed = sleep, False

    @property
    def counters(self):
        return self._state.counters

    def __enter__(self):
        if self._closed:
            raise ValueError("Client is closed")
        return self

    def __exit__(self, *exc):
        self.close()

    def close(self):
        if not self._closed:
            self._session.close()
        self._closed = True

    def _kwargs(self, request):
        headers = dict(request.headers)
        kwargs = dict(
            params=request.params,
            headers=headers,
        )
        if request.json is not _UNSET:
            kwargs["data"] = json.dumps(request.json).encode("utf-8")
            if "content-type" not in headers_lower(headers):
                headers["Content-Type"] = "application/json"
        elif request.multipart is not None:
            kwargs["files"] = [
                (p.name, (p.filename, p.value, p.content_type))
                for p in request.multipart
            ]
        elif request.form is not None:
            kwargs["data"] = list(request.form)
        else:
            kwargs["data"] = request.body
        return kwargs

    def _kind(self, exc):
        from urllib3.exceptions import ReadTimeoutError

        if isinstance(exc, self._requests.exceptions.SSLError) and _caused_by(
            exc, ssl.SSLCertVerificationError
        ):
            return "certificate"
        if exc.args and isinstance(exc.args[0], ReadTimeoutError):
            return "timeout"
        if isinstance(exc, self._requests.Timeout):
            return "timeout"
        if isinstance(exc, self._requests.exceptions.ContentDecodingError):
            return "decoding"
        return (
            "connection"
            if isinstance(
                exc,
                (
                    self._requests.ConnectionError,
                    self._requests.exceptions.ChunkedEncodingError,
                ),
            )
            else "request"
        )

    def _send(self, operation):
        request = operation.request
        kwargs = self._kwargs(request)
        timeout = (self._state.timeout.connect, self._state.timeout.read)
        try:
            prepared = self._session.prepare_request(
                self._requests.Request(request.method, request.url, **kwargs)
            )
            if any(prepared.hooks.values()):
                raise ValueError("Prepared response hooks bypass retry policy")
            settings = self._session.merge_environment_settings(
                prepared.url, proxies=self._proxies, stream=True, verify=None, cert=None
            )
            adapter = self._session.get_adapter(prepared.url)
            operation.attempt()
            # Session.send prepares its next redirect even with allow_redirects=False,
            # consuming the body. The adapter leaves it unread for our policy/limits.
            response = adapter.send(prepared, timeout=timeout, **settings)
            return response
        except self._requests.RequestException as exc:
            phase = "read" if isinstance(exc, self._requests.ReadTimeout) else "connect"
            raise operation.error(self._kind(exc), phase=phase) from None

    def _chunks(self, response, operation, size=65536, *, limited=True):
        iterator = response.iter_content(chunk_size=size)
        while True:
            operation.before_read()
            try:
                chunk = next(iterator, None)
            except self._requests.RequestException as exc:
                raise operation.error(
                    self._kind(exc), status=response.status_code, phase="read"
                ) from None
            if chunk is None:
                return
            operation.received(chunk, limited=limited)
            yield chunk

    def _status_error(self, response, operation):
        body = bytearray()
        # Error bodies have their own 8 KiB bound instead of max_body_bytes.
        for chunk in self._chunks(response, operation, 8192, limited=False):
            body.extend(chunk[: 8192 - len(body)])
            if len(body) == 8192:
                break
        return operation.error(
            "http_status", status=response.status_code, body=bytes(body)
        )

    def _open(self, operation):
        while True:
            response, error, headers = None, None, {}
            try:
                response = self._send(operation)
                headers = headers_lower(response.headers)
                status = response.status_code
                if status in {301, 302, 303, 307, 308}:
                    location = headers.get("location")
                    if not location:
                        raise operation.error("redirect_policy", status=status)
                    operation.follow(location, status=status)
                elif status in operation.request.expected_status:
                    return response
                else:
                    error = self._status_error(response, operation)
            except HttpTransportError as exc:
                error = exc
            except BaseException:
                if response is not None:
                    response.close()
                raise
            if response is not None:
                response.close()
            if error is None:
                continue
            error.retry_after = parse_retry_after(
                headers.get("retry-after"), self._state.wall_clock
            )
            if not operation.retry(error, "connect"):
                raise error from None
            self._sleep(operation.delay(headers))

    @contextmanager
    def stream(self, request, *, max_body_bytes=None):
        """Open a stream; no body replay occurs after handing it to the caller."""
        if self._closed:
            raise ValueError("Client is closed")
        operation = Operation(self._state, request, max_body_bytes)
        response = self._open(operation)
        try:
            yield _Stream(self, operation, response)
        finally:
            response.close()

    def get(
        self, url: str, *, max_body_bytes: int | None = DEFAULT_BODY_LIMIT, **options
    ) -> Response:
        """Buffer a GET response; options are the fields accepted by Request."""
        return self.request(
            Request("GET", url, **options), max_body_bytes=max_body_bytes
        )

    def post(
        self, url: str, *, max_body_bytes: int | None = DEFAULT_BODY_LIMIT, **options
    ) -> Response:
        """Buffer a POST response; retries still require explicit idempotency."""
        return self.request(
            Request("POST", url, **options), max_body_bytes=max_body_bytes
        )

    def request(self, request, *, max_body_bytes=DEFAULT_BODY_LIMIT):
        """Read a complete bounded body; eligible partial reads may be retried."""
        if self._closed:
            raise ValueError("Client is closed")
        operation = Operation(self._state, request, max_body_bytes)
        cached = self._state.cached(operation)
        if cached is not None:
            return cached
        while True:
            response = self._open(operation)
            try:
                body = b"".join(self._chunks(response, operation))
            except HttpTransportError as error:
                if not operation.retry(error, "read"):
                    raise
                response.close()
                self._sleep(operation.delay(headers_lower(response.headers)))
                continue
            finally:
                response.close()
            operation.delivered(len(body))
            result = Response(
                response.status_code,
                operation.request.target,
                headers_lower(response.headers),
                body,
                operation.request,
            )
            # A redirect's representation must not be stored under its original identity.
            if not operation.hops:
                self._state.store(request, result)
            return result
