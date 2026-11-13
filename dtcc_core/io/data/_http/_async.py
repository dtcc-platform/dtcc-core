"""Asynchronous transport over aiohttp with the shared transport contracts."""

import asyncio
import io
import os
import random
import ssl
import time
from collections.abc import Awaitable, Callable
from contextlib import asynccontextmanager

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


class _BorrowedStream(io.IOBase):
    """Let aiohttp close its payload wrapper without closing the caller's file."""

    def __init__(self, stream, errors):
        self._stream = stream
        self.errors = errors

    def read(self, size=-1):
        try:
            return self._stream.read(size)
        except Exception as exc:
            self.errors.append(exc)
            raise


class _AsyncStream:
    def __init__(self, client, operation, response):
        self._client, self._operation, self._response = client, operation, response
        self.status, self.headers = response.status, headers_lower(response.headers)
        self.target = operation.request.target
        self._started = self.closed = False

    async def iter_bytes(self, chunk_size=65536):
        """Deliver decoded chunks once without catching caller exceptions."""
        if chunk_size <= 0:
            raise ValueError("chunk_size must be positive")
        if self.closed:
            raise ValueError("Response stream is closed")
        if self._started:
            raise ValueError("Response stream has already been consumed")
        self._started = True
        async for chunk in self._client._chunks(
            self._response, self._operation, chunk_size
        ):
            self._operation.delivered(len(chunk))
            yield chunk


class AsyncHttpClient:
    """Client confined to one event loop; the caller schedules concurrent work."""

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
        sleep: Callable[[float], Awaitable[None]] = asyncio.sleep,
        wall_clock: Callable[[], float] = time.time,
        jitter_source: Callable[[], float] = random.random,
    ):
        import aiohttp

        self._aiohttp = aiohttp
        self._state = State(
            timeout=timeout,
            retry=retry,
            redirects=redirects,
            budget=budget,
            cache=cache,
            wall_clock=wall_clock,
            jitter_source=jitter_source,
        )
        self._ca_bundle = os.fspath(ca_bundle) if ca_bundle is not None else None
        self._proxy = proxy
        self._session, self._sleep, self._closed = None, sleep, False

    @property
    def counters(self):
        return self._state.counters

    async def __aenter__(self):
        await self._ensure_session()
        return self

    async def __aexit__(self, *exc):
        await self.close()

    async def _ensure_session(self):
        if self._closed:
            raise ValueError("Client is closed")
        if self._session is None:
            context = self._ssl_context()
            connector = self._aiohttp.TCPConnector(ssl=context)
            self._session = self._aiohttp.ClientSession(
                connector=connector,
                cookie_jar=self._aiohttp.DummyCookieJar(),
                trust_env=False,
                raise_for_status=False,
            )
            # aiohttp can replay a persistent-connection failure internally.
            # Disable that replay so budgets and retry policy own every attempt.
            if not hasattr(self._session, "_retry_connection"):
                session, self._session = self._session, None
                await session.close()
                raise RuntimeError(
                    "aiohttp no longer has ClientSession._retry_connection; "
                    "implicit connection replay cannot be disabled"
                )
            self._session._retry_connection = False

    def _ssl_context(self):
        """Trust ca_bundle alone when given, as requests does for the sync client."""
        if self._ca_bundle is None:
            from ..tls import download_ssl_context

            return download_ssl_context()
        if os.path.isdir(self._ca_bundle):
            return ssl.create_default_context(capath=self._ca_bundle)
        return ssl.create_default_context(cafile=self._ca_bundle)

    async def close(self):
        if self._session is not None:
            await self._session.close()
        self._closed = True

    def _kwargs(self, request, errors):
        headers = dict(request.headers)
        timeout = self._state.timeout
        kwargs = dict(
            params=request.params,
            proxy=self._proxy,
            headers=headers,
            allow_redirects=False,
            auto_decompress=True,
            raise_for_status=False,
            timeout=self._aiohttp.ClientTimeout(
                total=timeout.total,
                # ClientTimeout.connect would also bound waiting for a pool slot.
                sock_connect=timeout.connect,
                sock_read=timeout.read,
            ),
        )
        if request.json is not _UNSET:
            kwargs["data"] = request._json_body
            if "content-type" not in headers_lower(headers):
                headers["Content-Type"] = "application/json"
        elif request.multipart is not None:
            writer = self._aiohttp.MultipartWriter("form-data")
            for part in request.multipart:
                value = (
                    part.value
                    if isinstance(part.value, bytes)
                    else _BorrowedStream(part.value, errors)
                )
                payload = self._aiohttp.payload.get_payload(
                    value, content_type=part.content_type
                )
                disposition = {"name": part.name}
                if part.filename is not None:
                    disposition["filename"] = part.filename
                payload.set_content_disposition("form-data", **disposition)
                writer.append_payload(payload)
            kwargs["data"] = writer
        elif request.form is not None:
            kwargs["data"] = list(request.form)
        else:
            kwargs["data"] = request.body
        return kwargs

    def _kind(self, exc):
        from aiohttp.http_exceptions import ContentEncodingError

        if isinstance(exc.__cause__, ContentEncodingError):
            return "decoding"
        if isinstance(exc, TimeoutError):
            return "timeout"
        if isinstance(exc, self._aiohttp.ClientConnectorCertificateError):
            return "certificate"
        if isinstance(exc, self._aiohttp.ClientPayloadError):
            return "connection"
        return (
            "connection"
            if isinstance(exc, self._aiohttp.ClientConnectionError)
            else "request"
        )

    def _error(self, operation, exc, *, phase, status=None):
        if operation.local_errors:
            raise operation.local_errors[0] from None
        return operation.error(self._kind(exc), status=status, phase=phase)

    async def _send(self, operation):
        request = operation.request
        operation.local_errors.clear()
        kwargs = self._kwargs(request, operation.local_errors)
        operation.attempt()
        try:
            return await self._session.request(request.method, request.url, **kwargs)
        except (self._aiohttp.ClientError, TimeoutError) as exc:
            # Only a connect timeout is known to precede sending; the total
            # deadline can fire while waiting for response headers.
            phase = (
                "read"
                if isinstance(exc, TimeoutError)
                and not isinstance(exc, self._aiohttp.ConnectionTimeoutError)
                else "connect"
            )
            raise self._error(operation, exc, phase=phase) from None

    async def _chunks(self, response, operation, size=65536, *, limited=True):
        iterator = response.content.iter_chunked(size).__aiter__()
        while True:
            operation.before_read()
            try:
                chunk = await anext(iterator, None)
            except (self._aiohttp.ClientError, TimeoutError) as exc:
                raise self._error(
                    operation, exc, status=response.status, phase="read"
                ) from None
            if chunk is None:
                return
            operation.received(chunk, limited=limited)
            yield chunk

    async def _status_error(self, response, operation):
        body = bytearray()
        # Error bodies have their own 8 KiB bound instead of max_body_bytes.
        # A failed body read keeps the status error and the bytes read so far.
        try:
            async for chunk in self._chunks(response, operation, 8192, limited=False):
                body.extend(chunk[: 8192 - len(body)])
                if len(body) == 8192:
                    break
        except HttpTransportError:
            pass
        return operation.error(
            "http_status", status=response.status, body=bytes(body), phase="read"
        )

    async def _open(self, operation):
        while True:
            response, error, headers = None, None, {}
            try:
                response = await self._send(operation)
                headers = headers_lower(response.headers)
                status = response.status
                if status in {301, 302, 303, 307, 308}:
                    location = headers.get("location")
                    if not location:
                        raise operation.error(
                            "redirect_policy", status=status, phase="read"
                        )
                    operation.follow(location, status=status)
                elif status in operation.request.expected_status:
                    return response
                else:
                    error = await self._status_error(response, operation)
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
            await self._sleep(operation.delay(headers))

    @asynccontextmanager
    async def stream(self, request, *, max_body_bytes=None):
        """Open a stream; response cleanup also runs on cancellation."""
        await self._ensure_session()
        operation = Operation(self._state, request, max_body_bytes)
        response = await self._open(operation)
        stream = _AsyncStream(self, operation, response)
        try:
            yield stream
        finally:
            response.close()
            # Reads after the block would otherwise see an empty or broken body.
            stream.closed = True

    async def get(
        self, url: str, *, max_body_bytes: int | None = DEFAULT_BODY_LIMIT, **options
    ) -> Response:
        """Buffer a GET response; options are the fields accepted by Request."""
        return await self.request(
            Request("GET", url, **options), max_body_bytes=max_body_bytes
        )

    async def post(
        self, url: str, *, max_body_bytes: int | None = DEFAULT_BODY_LIMIT, **options
    ) -> Response:
        """Buffer a POST response; retries still require explicit idempotency."""
        return await self.request(
            Request("POST", url, **options), max_body_bytes=max_body_bytes
        )

    async def request(self, request, *, max_body_bytes=DEFAULT_BODY_LIMIT):
        """Read a bounded response; retry eligible failed reads before delivery."""
        await self._ensure_session()
        operation = Operation(self._state, request, max_body_bytes)
        cached = self._state.cached(operation)
        if cached is not None:
            return cached
        while True:
            response = await self._open(operation)
            try:
                chunks = [chunk async for chunk in self._chunks(response, operation)]
                body = b"".join(chunks)
            except HttpTransportError as error:
                if not operation.retry(error, "read"):
                    raise
                await self._sleep(operation.delay(headers_lower(response.headers)))
                continue
            finally:
                response.close()
            operation.delivered(len(body))
            result = Response(
                response.status,
                operation.request.target,
                headers_lower(response.headers),
                body,
                operation.request,
            )
            if not operation.hops:
                self._state.store(request, result)
            return result
