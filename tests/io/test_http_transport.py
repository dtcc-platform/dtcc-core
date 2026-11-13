"""Backend conformance against a loopback server, without public provider access."""

import asyncio
import gzip
import io
import json
import threading
import traceback
from collections import Counter
from contextlib import asynccontextmanager
from dataclasses import replace
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import parse_qs, urlsplit

import pytest

from dtcc_core.io.data._http import (
    AsyncHttpClient,
    Budget,
    HttpClient,
    HttpTransportError,
    MultipartPart,
    RedirectPolicy,
    Request,
    ResponseCache,
    RetryPolicy,
)


class Handler(BaseHTTPRequestHandler):
    def log_message(self, *args):
        pass

    def do_GET(self):
        parsed = urlsplit(self.path)
        self.server.calls[parsed.path] += 1
        query = parse_qs(parsed.query)
        status, headers, body = 200, {}, b'[1, null, "text"]'
        if parsed.path == "/echo":
            body = json.dumps({"query": query, "headers": dict(self.headers)}).encode()
        elif parsed.path == "/redirect":
            status, headers, body = (
                302,
                {"Location": query.get("to", ["/json"])[0]},
                b"",
            )
        elif parsed.path == "/retry":
            status = (
                503
                if self.server.calls[parsed.path] <= int(query.get("fail", [1])[0])
                else 200
            )
            headers = {"Retry-After": "0"}
            body = b"retry" if status == 503 else b"ok"
        elif parsed.path == "/status":
            status, body = int(query.get("code", [404])[0]), b"provider-secret"
            if "retry_after" in query:
                headers["Retry-After"] = query["retry_after"][0]
        elif parsed.path == "/bom":
            body = b"\xef\xbb\xbfhello"
        elif parsed.path == "/bad":
            body = b"secret is not JSON"
        elif parsed.path == "/gzip":
            body, headers = (
                b"not gzip" if "invalid" in query else gzip.compress(b"x" * 200000),
                {"Content-Encoding": "gzip"},
            )
        elif parsed.path == "/cache":
            headers = {key: values[0] for key, values in query.items()}
        elif parsed.path == "/broken":
            if self.server.calls[parsed.path] == 1:
                self.send_response(200)
                self.send_header("Content-Length", "200000")
                self.end_headers()
                self.wfile.write(b"x" * 70000)
                self.wfile.flush()
                self.close_connection = True
                return
            body = b"ok"
        self.send_response(status)
        self._write_headers(parsed.path, headers)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _write_headers(self, path, headers):
        if path == "/duplicate":
            self.send_header("Cache-Control", "max-age=60")
            self.send_header("Cache-Control", "no-store")
        for key, value in headers.items():
            self.send_header(key, value)

    def _read_body(self):
        if self.headers.get("Transfer-Encoding") != "chunked":
            return self.rfile.read(int(self.headers.get("Content-Length", 0)))
        chunks = []
        while True:
            line = self.rfile.readline().strip()
            if not line:
                return b"".join(chunks)
            length = int(line, 16)
            if length == 0:
                self.rfile.readline()
                return b"".join(chunks)
            chunks.append(self.rfile.read(length))
            assert self.rfile.read(2) == b"\r\n"

    def do_POST(self):
        body = self._read_body()
        self.server.received.append((dict(self.headers), body))
        self.send_response(200)
        result = json.dumps(
            {
                "body": body.decode(),
                "content_type": self.headers.get("Content-Type"),
            }
        ).encode()
        self.send_header("Content-Length", str(len(result)))
        self.end_headers()
        try:
            self.wfile.write(result)
        except (BrokenPipeError, ConnectionResetError):
            pass


@pytest.fixture
def server():
    calls = Counter()
    received = []

    httpd = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    httpd.calls, httpd.received = calls, received
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()
    try:
        yield f"http://127.0.0.1:{httpd.server_port}", calls, received
    finally:
        httpd.shutdown()
        httpd.server_close()
        thread.join()


@pytest.fixture(params=["sync", "async"])
def backend(request):
    return request.param


class Bridge:
    """Exercise the same assertions through the two scheduling interfaces."""

    def __init__(self, client, backend):
        self.client, self.backend = client, backend

    async def request(self, request, **kwargs):
        if self.backend == "sync":
            return self.client.request(request, **kwargs)
        return await self.client.request(request, **kwargs)

    @asynccontextmanager
    async def stream(self, request, **kwargs):
        if self.backend == "sync":
            with self.client.stream(request, **kwargs) as response:
                yield response
        else:
            async with self.client.stream(request, **kwargs) as response:
                yield response

    async def chunks(self, response, size=65536):
        if self.backend == "sync":
            for chunk in response.iter_bytes(size):
                yield chunk
        else:
            async for chunk in response.iter_bytes(size):
                yield chunk


def run(backend, scenario, **kwargs):
    async def main():
        if backend == "sync":
            with HttpClient(**kwargs) as client:
                await scenario(Bridge(client, backend))
        else:
            async with AsyncHttpClient(**kwargs) as client:
                await scenario(Bridge(client, backend))

    asyncio.run(main())


def test_buffered_decoding_and_repeated_query(backend, server):
    url, _, _ = server

    async def scenario(api):
        response = await api.request(Request("get", url + "/json"))
        assert response.json() == [1, None, "text"]
        response = await api.request(
            Request("GET", url + "/echo", params=(("a", "1"), ("a", "2")))
        )
        assert response.json()["query"] == {"a": ["1", "2"]}
        response = await api.request(Request("GET", url + "/bom"))
        assert response.text("utf-8-sig") == "hello"
        assert response.text().startswith("\ufeff")
        response = await api.request(Request("GET", url + "/bad?token=secret"))
        with pytest.raises(HttpTransportError) as error:
            response.json()
        assert error.value.kind == "decoding"
        assert "secret" not in repr(response)
        assert "secret" not in str(error.value)
        assert api.client.counters.logical_requests == 4

    run(backend, scenario)


def test_convenience_calls_mapping_inputs_and_shared_cache(backend, server):
    url, calls, _ = server

    async def scenario(api):
        async def call(method, *args, **kwargs):
            result = getattr(api.client, method)(*args, **kwargs)
            return await result if backend == "async" else result

        response = await call(
            "get",
            url + "/echo",
            params={"year": 2025, "resolution": 0.5},
            headers={"Accept": "application/json"},
            public=True,
            cache_ttl=60,
        )
        assert response.json()["query"] == {"year": ["2025"], "resolution": ["0.5"]}
        assert response.json()["headers"]["Accept"] == "application/json"
        cached = await api.request(
            Request(
                "GET",
                url + "/echo",
                params=(("year", "2025"), ("resolution", "0.5")),
                headers=(("Accept", "application/json"),),
                public=True,
                cache_ttl=60,
            )
        )
        assert cached.body == response.body and calls["/echo"] == 1
        assert api.client.counters.cache_hits == 1
        result = await call("post", url, json=None, expected_status={200})
        assert result.json()["body"] == "null"
        result = await call(
            "post", url, form={"data": "query with spaces", "year": 2025}
        )
        assert parse_qs(result.json()["body"]) == {
            "data": ["query with spaces"],
            "year": ["2025"],
        }
        result = await call("post", url, form=(("field", "one"), ("field", "two")))
        assert parse_qs(result.json()["body"]) == {"field": ["one", "two"]}
        with pytest.raises(HttpTransportError, match="resource_limit"):
            await call("get", url + "/echo", max_body_bytes=1)
        with pytest.raises(HttpTransportError, match="resource_limit"):
            await call("post", url, body=b"payload", max_body_bytes=1)
        with pytest.raises(HttpTransportError, match="http_status"):
            await call("get", url + "/status?code=201", expected_status={200})

    run(backend, scenario, cache=ResponseCache())


@pytest.mark.parametrize("retries", [0, 1])
@pytest.mark.parametrize("stream", [False, True])
def test_retry_after_exposed_to_outer_retry_loop(backend, server, retries, stream):
    url, calls, _ = server
    sleeps = []

    async def sleep_async(delay):
        sleeps.append(delay)

    async def scenario(api):
        request = Request(
            "GET", url + "/status", params={"code": 429, "retry_after": 90}
        )
        with pytest.raises(HttpTransportError) as caught:
            if stream:
                async with api.stream(request):
                    pytest.fail("429 must raise before handing out the stream")
            else:
                await api.request(request)
        error = caught.value
        assert error.status == 429 and error.retry_after == 90.0
        assert error.attempts == retries + 1 and calls["/status"] == retries + 1
        assert error.retry_allowed_by_policy is bool(retries)
        assert "retry_after" not in error.metadata()
        assert "retry_allowed_by_policy" in error.metadata()
        assert error.redirect_target is None
        assert sleeps == [3] * retries

    run(
        backend,
        scenario,
        retry=RetryPolicy(retries=retries, max_delay=3),
        sleep=sleeps.append if backend == "sync" else sleep_async,
    )


def test_retry_after_date_exposed_when_policy_ignores_header(backend, server):
    url, _, _ = server

    async def scenario(api):
        for value, expected in (("Thu, 01 Jan 1970 00:01:45 GMT", 5.0), ("bad", None)):
            with pytest.raises(HttpTransportError) as caught:
                await api.request(
                    Request(
                        "GET",
                        url + "/status",
                        params={"code": 503, "retry_after": value},
                    )
                )
            assert caught.value.retry_after == expected

    run(
        backend,
        scenario,
        retry=RetryPolicy(respect_retry_after=False),
        wall_clock=lambda: 100,
    )


@pytest.mark.parametrize(
    "destination, target",
    [
        ("/json?token=secret#secret", "/json"),
        (
            "https://user:secret@elsewhere.test/path?key=secret#secret",
            "https://elsewhere.test/path",
        ),
        ("https://host:bad/secret", None),
        ("file:///secret", None),
    ],
)
def test_rejected_redirect_keeps_status_and_sanitized_destination(
    backend, server, destination, target
):
    url, calls, _ = server

    async def scenario(api):
        with pytest.raises(HttpTransportError) as caught:
            await api.request(
                Request(
                    "GET",
                    url + "/redirect",
                    params={"to": destination},
                    expected_status={200, 302},
                )
            )
        error = caught.value
        assert error.status == 302 and error.kind == "redirect_policy"
        assert error.attempts == 1
        assert error.redirect_target == (
            url + target if target and target.startswith("/") else target
        )
        assert error.error_body == b""
        assert calls["/json"] == 0
        diagnostics = str(error) + repr(error.metadata()) + repr(error.redirect_target)
        diagnostics += "".join(traceback.format_exception(error))
        assert "secret" not in diagnostics

    run(backend, scenario)


def test_status_retries_and_counters(backend, server):
    url, calls, _ = server
    sleeps = []

    async def sleep_async(delay):
        sleeps.append(delay)

    async def scenario(api):
        result = await api.request(Request("GET", url + "/retry?fail=2"))
        assert result.body == b"ok"
        assert calls["/retry"] == 3
        counters = api.client.counters
        assert counters.network_attempts == 3
        assert counters.logical_requests == 1
        assert counters.received_body_bytes == 12
        assert counters.delivered_body_bytes == 2
        assert sleeps == [0, 0]

    run(
        backend,
        scenario,
        retry=RetryPolicy(retries=2),
        sleep=sleeps.append if backend == "sync" else sleep_async,
    )


def test_exhausted_status_and_safe_diagnostics(backend, server):
    url, calls, _ = server

    async def scenario(api):
        request = Request(
            "GET",
            url + "/status?code=503&apikey=private-key",
            headers=(("Authorization", "Bearer private-key"),),
        )
        with pytest.raises(HttpTransportError) as caught:
            await api.request(request)
        error = caught.value
        assert error.status == 503 and error.kind == "http_status"
        assert error.attempts == 2
        assert error.error_body == b"provider-secret"
        assert calls["/status"] == 2
        diagnostics = (
            str(error) + repr(error) + repr(request) + json.dumps(error.metadata())
        )
        diagnostics += "".join(traceback.format_exception(error))
        assert "private-key" not in diagnostics and "provider-secret" not in diagnostics

    run(backend, scenario, retry=RetryPolicy(retries=1, backoff=0))


def test_status_selection_and_default_no_retry(backend, server):
    url, calls, _ = server

    async def scenario(api):
        with pytest.raises(HttpTransportError):
            await api.request(Request("GET", url + "/status?code=503"))
        assert calls["/status"] == 1
        response = await api.request(Request("GET", url + "/status?code=201"))
        assert response.status == 201
        with pytest.raises(HttpTransportError):
            await api.request(
                Request("GET", url + "/status?code=201", expected_status={200})
            )

    run(backend, scenario)


def test_limits_and_compressed_body(backend, server):
    url, calls, _ = server

    async def scenario(api):
        request = Request("GET", url + "/json")
        response = await api.request(request, max_body_bytes=17)
        assert len(response.body) == 17
        with pytest.raises(HttpTransportError, match="resource_limit"):
            await api.request(request, max_body_bytes=16)
        with pytest.raises(HttpTransportError, match="resource_limit"):
            await api.request(Request("GET", url + "/gzip"), max_body_bytes=1000)
        assert calls["/gzip"] == 1
        assert api.client.counters.received_body_bytes > 1000

    run(backend, scenario, retry=RetryPolicy(retries=2))


def test_budget_attempts_and_received_bytes(backend, server):
    url, calls, _ = server

    async def scenario(api):
        with pytest.raises(HttpTransportError, match="resource_limit"):
            await api.request(Request("GET", url + "/retry?fail=10"))
        assert calls["/retry"] == 2
        assert api.client.counters.network_attempts == 2

    run(
        backend,
        scenario,
        retry=RetryPolicy(retries=10, backoff=0),
        budget=Budget(max_attempts=2),
    )

    async def bytes_scenario(api):
        with pytest.raises(HttpTransportError, match="resource_limit"):
            await api.request(Request("GET", url + "/json"))
        assert api.client.counters.received_body_bytes == 17

    run(backend, bytes_scenario, budget=Budget(max_received_bytes=16))


def test_stream_and_caller_failure(backend, server):
    url, _, _ = server
    failure = OSError("caller failed")

    async def scenario(api):
        async with api.stream(Request("GET", url + "/json")) as response:
            assert response.status == 200 and response.headers
            chunks = [part async for part in api.chunks(response, 3)]
            assert b"".join(chunks) == b'[1, null, "text"]'
            with pytest.raises(ValueError):
                _ = [part async for part in api.chunks(response)]
        with pytest.raises(OSError) as caught:
            async with api.stream(Request("GET", url + "/json")):
                raise failure
        assert caught.value is failure
        # Subsequent work still succeeds after early exit and caller failure.
        assert (await api.request(Request("GET", url + "/json"))).status == 200

    run(backend, scenario)


def test_partial_buffer_retry_and_stream_no_replay(backend, server):
    url, calls, _ = server

    async def scenario(api):
        result = await api.request(Request("GET", url + "/broken"))
        assert result.body == b"ok" and calls["/broken"] == 2
        assert api.client.counters.received_body_bytes >= 65538
        assert api.client.counters.delivered_body_bytes == 2
        calls["/broken"] = 0
        chunks = []
        with pytest.raises(HttpTransportError):
            async with api.stream(Request("GET", url + "/broken")) as response:
                async for chunk in api.chunks(response, 1024):
                    chunks.append(chunk)
        assert chunks and calls["/broken"] == 1

    run(backend, scenario, retry=RetryPolicy(retries=1, backoff=0))


def test_redirects(backend, server):
    url, calls, _ = server

    async def denied(api):
        with pytest.raises(HttpTransportError, match="redirect_policy"):
            await api.request(Request("GET", url + "/redirect"))
        assert calls["/json"] == 0

    run(backend, denied)

    async def allowed(api):
        result = await api.request(Request("GET", url + "/redirect"))
        assert result.json() == [1, None, "text"]
        assert api.client.counters.network_attempts == 2
        with pytest.raises(HttpTransportError, match="redirect_policy"):
            await api.request(Request("GET", url + "/redirect?to=/redirect"))

    run(backend, allowed, redirects=RedirectPolicy(max_hops=1))


def test_json_null_form_and_multipart(backend, server):
    url, _, received = server
    file = io.BytesIO(b"binary-data")

    async def scenario(api):
        response = await api.request(Request("POST", url, json=None))
        assert response.json()["body"] == "null"
        response = await api.request(
            Request("POST", url, form=(("a", "1"), ("a", "2")))
        )
        assert response.json()["body"] == "a=1&a=2"
        request = Request(
            "POST",
            url,
            multipart=(
                MultipartPart("files", b"first", "a.txt", "text/plain"),
                MultipartPart("files", file, "b.bin"),
            ),
        )
        await api.request(request)
        headers, body = received[-1]
        assert headers["Content-Type"].startswith("multipart/form-data")
        assert body.count(b'name="files"') == 2
        assert b'filename="b.bin"' in body and b"binary-data" in body
        assert not file.closed

    run(backend, scenario)


def test_cache_and_cache_hit_limits(backend, server):
    url, calls, _ = server
    cache = ResponseCache()

    async def scenario(api):
        request = Request("GET", url + "/json", public=True, cache_ttl=60)
        first = await api.request(request)
        first.headers["modified"] = "yes"
        second = await api.request(request)
        assert "modified" not in second.headers
        assert calls["/json"] == 1
        assert api.client.counters.cache_hits == 1
        assert api.client.counters.received_body_bytes == 17
        assert api.client.counters.delivered_body_bytes == 34
        with pytest.raises(HttpTransportError, match="resource_limit"):
            await api.request(request, max_body_bytes=16)
        assert calls["/json"] == 1
        cache.invalidate(request)
        await api.request(request)
        assert calls["/json"] == 2

    run(backend, scenario, cache=cache)


@pytest.mark.parametrize(
    "headers",
    [
        (("Authorization", "secret"),),
        (("Cookie", "secret"),),
        (("Cache-Control", "no-cache"),),
        (("Pragma", "no-cache"),),
    ],
)
def test_cache_request_bypass(backend, server, headers):
    url, calls, _ = server

    async def scenario(api):
        request = Request(
            "GET", url + "/json", headers=headers, public=True, cache_ttl=60
        )
        await api.request(request)
        await api.request(request)
        assert calls["/json"] == 2

    run(backend, scenario, cache=ResponseCache())


@pytest.mark.parametrize(
    "params",
    [
        (("Cache-Control", "No-StOrE"),),
        (("Cache-Control", "private"),),
        (("Cache-Control", "no-cache"),),
        (("Vary", "Accept"),),
        (("Set-Cookie", "session=secret"),),
        (("Cache-Control", "max-age=0"),),
    ],
)
def test_cache_response_bypass_and_no_cookie_retention(backend, server, params):
    url, calls, _ = server

    async def scenario(api):
        request = Request(
            "GET", url + "/cache", params=params, public=True, cache_ttl=60
        )
        await api.request(request)
        await api.request(request)
        assert calls["/cache"] == 2
        echo = await api.request(Request("GET", url + "/echo"))
        assert "Cookie" not in echo.json()["headers"]

    run(backend, scenario, cache=ResponseCache())


def test_cache_ttl_version_and_eviction(backend, server):
    url, calls, _ = server
    now = [0.0]
    cache = ResponseCache(max_entries=1, clock=lambda: now[0])

    async def scenario(api):
        request = Request("GET", url + "/json", public=True, cache_ttl=60)
        await api.request(request)
        now[0] = 2
        await api.request(replace(request, cache_ttl=1))
        assert calls["/json"] == 2
        await api.request(replace(request, cache_token="updated"))
        await api.request(request)
        assert calls["/json"] == 4
        cache.clear()
        await api.request(request)
        assert calls["/json"] == 5

    run(backend, scenario, cache=cache)


def test_cache_cannot_override_expected_status(backend, server):
    url, calls, _ = server

    async def scenario(api):
        request = Request("GET", url + "/json", public=True, cache_ttl=60)
        await api.request(request)
        with pytest.raises(HttpTransportError) as caught:
            await api.request(replace(request, expected_status={204}))
        assert caught.value.status == 200
        assert calls["/json"] == 2

    run(backend, scenario, cache=ResponseCache())


def test_duplicate_cache_control_headers_are_honored(backend, server):
    url, calls, _ = server

    async def scenario(api):
        request = Request("GET", url + "/duplicate", public=True, cache_ttl=60)
        response = await api.request(request)
        assert "no-store" in response.headers["cache-control"]
        await api.request(request)
        assert calls["/duplicate"] == 2

    run(backend, scenario, cache=ResponseCache())


@pytest.mark.parametrize(
    "headers",
    [
        (("Accept", "text/csv"), ("Accept", "application/json")),
        (("Accept", "application/json"), ("accept", "text/csv")),
    ],
)
def test_duplicate_request_headers_are_rejected(backend, server, headers):
    url, calls, _ = server

    async def scenario(api):
        for pairs in (headers, tuple(reversed(headers))):
            with pytest.raises(ValueError, match="unique"):
                await api.request(
                    Request("GET", url, headers=pairs, public=True, cache_ttl=60)
                )
        assert not calls

    run(backend, scenario, cache=ResponseCache())


@pytest.mark.parametrize("limit", [16, 17])
def test_received_budget_blocks_subsequent_sends(backend, server, limit):
    url, calls, _ = server

    async def scenario(api):
        request = Request("GET", url + "/json")
        if limit == 16:
            with pytest.raises(HttpTransportError, match="resource_limit"):
                await api.request(request)
        else:
            assert len((await api.request(request)).body) == 17
        for _ in range(3):
            with pytest.raises(HttpTransportError, match="resource_limit"):
                await api.request(request)
        assert calls["/json"] == 1
        assert api.client.counters.received_body_bytes == 17

    run(backend, scenario, budget=Budget(max_received_bytes=limit))


def test_exhausted_budget_stops_other_open_streams(backend, server):
    url, calls, _ = server

    async def scenario(api):
        request = Request("GET", url + "/json")
        async with api.stream(request) as first, api.stream(request) as second:
            with pytest.raises(HttpTransportError, match="resource_limit"):
                _ = [chunk async for chunk in api.chunks(first)]
            with pytest.raises(HttpTransportError, match="resource_limit"):
                _ = [chunk async for chunk in api.chunks(second)]
        assert calls["/json"] == 2
        assert api.client.counters.received_body_bytes == 17

    run(backend, scenario, budget=Budget(max_received_bytes=1))


def test_delivery_budget_prevents_more_downloads(backend, server):
    url, calls, _ = server

    async def scenario(api):
        for _ in range(2):
            with pytest.raises(HttpTransportError, match="resource_limit"):
                await api.request(Request("GET", url + "/json"))
        assert calls["/json"] == 1

    run(backend, scenario, budget=Budget(max_delivered_bytes=1))


def test_multipart_fields_have_no_invented_filename(backend, server):
    url, _, received = server

    async def scenario(api):
        from email.parser import BytesParser
        from email.policy import default

        stream = io.BytesIO(b"stream metadata")
        await api.request(
            Request(
                "POST",
                url,
                multipart=(
                    MultipartPart("metadata", b"{}", content_type="application/json"),
                    MultipartPart("stream", stream),
                    MultipartPart("file", b"contents", filename="file.bin"),
                ),
            )
        )
        headers, body = received[-1]
        message = BytesParser(policy=default).parsebytes(
            f"Content-Type: {headers['Content-Type']}\r\n\r\n".encode() + body
        )
        parts = list(message.iter_parts())
        assert [p.get_filename() for p in parts] == [None, None, "file.bin"]
        assert [p.get_payload(decode=True) for p in parts] == [
            b"{}",
            b"stream metadata",
            b"contents",
        ]
        assert not stream.closed

    run(backend, scenario)


def test_malformed_encoding_is_not_retried(backend, server):
    url, calls, _ = server

    async def scenario(api):
        with pytest.raises(HttpTransportError) as caught:
            await api.request(Request("GET", url + "/gzip?invalid=1"))
        assert caught.value.kind == "decoding"
        assert not caught.value.retry_allowed_by_policy
        assert caught.value.attempts == 1
        assert calls["/gzip"] == 1

    run(backend, scenario, retry=RetryPolicy(retries=1, backoff=0))


@pytest.mark.parametrize("exception_type", [OSError, ValueError, RuntimeError])
def test_upload_reader_exception_is_preserved(backend, server, exception_type):
    url, _, _ = server
    failure = exception_type("local upload error")

    class BrokenReader(io.BytesIO):
        def read(self, size=-1):
            raise failure

    async def scenario(api):
        reader = BrokenReader(b"file contents")
        with pytest.raises(exception_type) as caught:
            await api.request(
                Request(
                    "POST", url, multipart=(MultipartPart("file", reader, "file.bin"),)
                )
            )
        assert caught.value is failure
        assert not reader.closed

    run(backend, scenario)


@pytest.mark.parametrize("mode", ["buffered", "stream", "error_status"])
def test_async_reader_error_after_response_headers(mode, monkeypatch):
    from dtcc_core.io.data._http._async import _BorrowedStream

    release_reader = threading.Event()
    recorded_failure = threading.Event()
    failure = RuntimeError("late upload failure")
    original_read = _BorrowedStream.read

    def record_failure(self, size=-1):
        try:
            return original_read(self, size)
        except Exception:
            recorded_failure.set()
            raise

    monkeypatch.setattr(_BorrowedStream, "read", record_failure)

    class DelayedReader(io.BytesIO):
        def read(self, size=-1):
            assert release_reader.wait(5), "server did not release the upload reader"
            raise failure

    class EarlyResponse(BaseHTTPRequestHandler):
        def log_message(self, *args):
            pass

        def do_POST(self):
            self.send_response(503 if mode == "error_status" else 200)
            self.send_header("Content-Length", "2")
            self.end_headers()
            self.wfile.flush()
            # The test releases the reader once the client has received the headers.
            assert recorded_failure.wait(5), "upload failure was not recorded"

    httpd = ThreadingHTTPServer(("127.0.0.1", 0), EarlyResponse)
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()

    async def scenario():
        reader = DelayedReader(b"data")
        async with AsyncHttpClient() as client:
            original = client._session.request

            async def after_headers(*args, **kwargs):
                response = await original(*args, **kwargs)
                release_reader.set()
                return response

            client._session.request = after_headers
            request = Request(
                "POST",
                f"http://127.0.0.1:{httpd.server_port}/",
                multipart=(MultipartPart("file", reader, "file.bin"),),
            )
            with pytest.raises(RuntimeError) as caught:
                if mode == "stream":
                    async with client.stream(request) as response:
                        _ = [chunk async for chunk in response.iter_bytes()]
                else:
                    await client.request(request)
            assert caught.value is failure
            assert not reader.closed

    try:
        asyncio.run(scenario())
    finally:
        release_reader.set()
        httpd.shutdown()
        httpd.server_close()
        thread.join()


def test_status_error_body_ignores_caller_body_limit(backend, server):
    url, calls, _ = server

    async def scenario(api):
        with pytest.raises(HttpTransportError) as caught:
            await api.request(
                Request("GET", url + "/status?code=503"), max_body_bytes=4
            )
        error = caught.value
        assert (error.kind, error.status, error.attempts) == ("http_status", 503, 2)
        assert error.error_body == b"provider-secret"
        assert calls["/status"] == 2

    run(backend, scenario, retry=RetryPolicy(retries=1, backoff=0))


def test_unencodable_json_body_consumes_no_attempt(backend, server):
    url, _, _ = server

    async def scenario(api):
        with pytest.raises(TypeError):
            await api.request(Request("POST", url + "/echo", json={"x": object()}))
        assert api.client.counters.network_attempts == 0
        # The shared one-attempt budget is still available.
        assert (await api.request(Request("GET", url + "/json"))).status == 200

    run(backend, scenario, budget=Budget(max_attempts=1))


def test_non_idempotent_upload_stream_is_sent_with_retries_enabled(backend, server):
    url, _, received = server

    async def scenario(api):
        part = MultipartPart("file", io.BytesIO(b"payload"), "x.bin")
        request = Request("POST", url + "/upload", multipart=(part,))
        assert (await api.request(request)).status == 200
        assert len(received) == 1 and b"payload" in received[0][1]

    run(backend, scenario, retry=RetryPolicy(retries=2))
