"""Standalone author interface and native HTTPX2 failure-path contracts."""

import gzip
import io
import json
import pickle
import ssl
from dataclasses import asdict
from pathlib import Path

import httpx2
import pytest

from dtcc_core.io import HttpLimitError, HttpTransportError, ResponseCache, http_client
from dtcc_core.io._http import HttpClient, Request

from . import test_http_transport

server = test_http_transport.server

URL = "https://weather.test/archive?token=secret"


class Body(httpx2.SyncByteStream):
    def __init__(self, chunks):
        self.chunks, self.reads, self.closed = chunks, 0, False

    def __iter__(self):
        for chunk in self.chunks:
            if isinstance(chunk, BaseException):
                raise chunk
            self.reads += 1
            yield chunk

    def close(self):
        self.closed = True


@pytest.fixture
def native():
    clients = []

    def create(events, **options):
        remaining, calls, delays = iter(events), [], []

        def handler(request):
            calls.append(request)
            event = next(remaining)
            if isinstance(event, BaseException):
                raise event
            return event

        client = http_client(**options)
        backend = client._backend
        backend._client = httpx2.Client(
            transport=httpx2.MockTransport(handler), trust_env=False
        )
        backend._sleep = delays.append
        clients.append(client)
        return client, calls, delays

    yield create
    for client in clients:
        client.close()


def reply(body=b'{"temperature":12}', *, status=200, headers=None):
    stream = Body(body if isinstance(body, list) else [body])
    response = httpx2.Response(status, stream=stream, headers=headers)
    response.test_body = stream
    return response


@pytest.mark.parametrize(
    "value", [{"error": True}, [1, None], 42, 1.5, "text", True, None]
)
def test_plain_function_returns_json_unchanged(native, value):
    client, calls, _ = native([reply(json.dumps(value).encode())])

    def weather(day):
        with client as http:
            return http.json(URL, params={"date": day})

    assert weather("2026-01-01") == value
    assert calls[0].url.params["date"] == "2026-01-01"
    assert client._backend._client.is_closed


def test_defaults_and_new_clients_have_independent_budgets():
    with http_client() as first, http_client() as second:
        state = first._backend._state
        assert state.timeout.connect == state.timeout.read == 30
        assert state.retry.retries == 2
        assert state.retry.respect_retry_after
        assert state.redirects.max_hops == 5
        assert state.cache is None
        assert state.budget.limits == dict(
            logical_requests=2000,
            network_attempts=2000,
            received_body_bytes=100 * 1024 * 1024,
            delivered_body_bytes=100 * 1024 * 1024,
        )
        assert first.counters is not second.counters


@pytest.mark.parametrize(
    "options",
    [
        {"timeout": 0},
        {"timeout": float("nan")},
        {"timeout": True},
        {"retries": -1},
        {"retries": True},
        {"max_requests": -1},
        {"max_attempts": 1.5},
        {"max_download_bytes": True},
        {"cache_ttl": -1},
        {"cache_ttl": float("inf")},
        {"cache_ttl": "1"},
        {"cache": {}},
        {"cache_namespace": 1},
        {"trust_env": "yes"},
    ],
)
def test_invalid_configuration_before_io(options):
    with pytest.raises(ValueError):
        http_client(**options)


def test_overrides_and_helper_calls_share_budget(native):
    client, calls, _ = native(
        [reply(), reply()],
        timeout=60,
        max_requests=1,
        max_attempts=10,
        max_download_bytes=1000,
        retries=0,
    )

    def helper(http):
        return http.json(URL)

    with client:
        assert helper(client) == {"temperature": 12}
        with pytest.raises(HttpLimitError) as caught:
            helper(client)
        assert (caught.value.resource, caught.value.limit, caught.value.observed) == (
            "logical_requests",
            1,
            2,
        )
        assert client._backend._state.timeout.read == 60
    assert len(calls) == 1


def test_utf8_and_explicit_encoding(native):
    client, _, _ = native(
        [reply("å".encode()), reply(b"\xef\xbb\xbfhello"), reply(b"\xe5")]
    )
    with client:
        assert client.text(URL) == "å"
        assert client.text(URL, encoding="utf-8-sig") == "hello"
        assert client.text(URL, encoding="latin-1") == "å"


@pytest.mark.parametrize(
    "method,body,kwargs",
    [
        ("json", b"secret invalid json", {}),
        ("text", b"\xffsecret", {}),
        ("text", b"secret", {"encoding": "not-a-real-encoding"}),
    ],
)
def test_safe_decode_errors_are_not_retried(native, method, body, kwargs):
    response = reply(body)
    client, calls, delays = native([response])
    with client, pytest.raises(HttpTransportError) as caught:
        getattr(client, method)(URL, **kwargs)
    assert caught.value.kind == "decoding"
    assert "secret" not in str(caught.value)
    assert caught.value.__cause__ is None and caught.value.__suppress_context__
    assert response.is_closed and len(calls) == 1 and delays == []


@pytest.mark.parametrize("status", [429, 500, 502, 503, 504])
def test_status_retries_honor_retry_after(native, status):
    failed = reply(b"busy", status=status, headers={"Retry-After": "1"})
    client, calls, delays = native([failed, reply()])
    with client:
        assert client.json(URL) == {"temperature": 12}
        assert client.counters.network_attempts == 2
        assert client.counters.logical_requests == 1
    assert failed.is_closed and len(calls) == 2 and delays == [1]


def test_retry_exhaustion_and_override(native):
    client, calls, delays = native(
        [reply(status=503, headers={"Retry-After": "999"}) for _ in range(3)],
        respect_retry_after=False,
    )
    with client, pytest.raises(HttpTransportError) as caught:
        client.json(URL)
    assert (caught.value.kind, caught.value.status, caught.value.attempts) == (
        "http_status",
        503,
        3,
    )
    assert len(calls) == 3 and delays == [0.2, 0.4]


@pytest.mark.parametrize(
    "failure",
    [
        httpx2.ConnectError,
        httpx2.ConnectTimeout,
        httpx2.ReadTimeout,
        httpx2.ProxyError,
        httpx2.RemoteProtocolError,
    ],
)
def test_transport_retries(native, failure):
    client, calls, delays = native([failure("secret"), reply()])
    with client:
        assert client.json(URL) == {"temperature": 12}
    assert len(calls) == 2 and delays == [0.2]


def test_partial_reads_count_bytes_and_limit_stops_retry(native):
    responses = [reply([b"12345678", httpx2.ReadTimeout("secret")]) for _ in range(2)]
    client, calls, delays = native(responses, max_download_bytes=12)
    with client, pytest.raises(HttpLimitError) as caught:
        client.text(URL)
    assert (caught.value.resource, caught.value.limit, caught.value.observed) == (
        "received_body_bytes",
        12,
        16,
    )
    assert client.counters.received_body_bytes == 16
    assert client.counters.delivered_body_bytes == 0
    assert all(response.is_closed for response in responses)
    assert len(calls) == 2 and delays == [0.2]
    restored = pickle.loads(pickle.dumps(caught.value))
    assert restored.metadata() == caught.value.metadata()
    assert "secret" not in repr(restored.metadata())


def test_attempt_budget_blocks_retry(native):
    client, calls, delays = native(
        [reply(b"busy", status=503, headers={"Retry-After": "30"})],
        max_attempts=1,
    )
    with client:
        with pytest.raises(HttpTransportError) as caught:
            client.json(URL)
        assert (
            caught.value.kind,
            caught.value.status,
            caught.value.retry_after,
            caught.value.attempts,
        ) == ("http_status", 503, 30.0, 1)
        assert len(calls) == 1 and delays == []
        with pytest.raises(HttpLimitError) as limit_caught:
            client.json(URL)
        assert (
            limit_caught.value.resource == "network_attempts"
            and limit_caught.value.observed == 2
        )


def test_download_stops_at_crossing_chunk(native):
    response = reply([b"12345678"] * 100)
    client, _, _ = native([response], max_download_bytes=10)
    with client, pytest.raises(HttpLimitError):
        client.text(URL)
    assert response.test_body.reads == 2 and response.test_body.closed


@pytest.mark.parametrize(
    "headers",
    [
        {"Cache-Control": "no-store"},
        {"Cache-Control": "private"},
        {"Cache-Control": "no-cache"},
        {"Cache-Control": "max-age=0"},
        {"Cache-Control": "max-age=1", "Age": "2"},
        {"Set-Cookie": "session=secret"},
        {"Vary": "Accept"},
    ],
)
def test_response_cache_restrictions(native, headers):
    client, calls, _ = native([reply(headers=headers), reply()], cache_ttl=60)
    with client:
        client.json(URL)
        client.json(URL)
    assert len(calls) == 2


@pytest.mark.parametrize(
    "kwargs",
    [
        {"public": False},
        {"credentials": {"apikey": "secret"}},
        {"headers": {"Authorization": "secret"}},
        {"headers": {"Cookie": "secret"}},
    ],
)
def test_private_requests_never_cached(native, kwargs):
    client, calls, _ = native([reply(), reply()], cache_ttl=60)
    with client:
        client.json(URL, **kwargs)
        client.json(URL, **kwargs)
    assert len(calls) == 2
    if "credentials" in kwargs:
        assert calls[0].url.params["apikey"] == "secret"


def test_credential_conflicts_rejected_without_io(native):
    client, calls, _ = native([])
    with client, pytest.raises(ValueError, match="Credential names"):
        client.json(URL, params={"apikey": "one"}, credentials={"apikey": "two"})
    assert not calls


def test_cache_tokens_logical_limits_and_delivered_bytes(native):
    client, calls, _ = native(
        [reply(b"abc"), reply(b"def")],
        cache_ttl=60,
        max_download_bytes=10,
        max_requests=3,
    )
    with client:
        assert client.text(URL, cache_token=1) == "abc"
        assert client.text(URL, cache_token=1) == "abc"
        assert client.text(URL, cache_token=2) == "def"
        assert asdict(client.counters) == dict(
            logical_requests=3,
            network_attempts=2,
            cache_hits=1,
            received_body_bytes=6,
            delivered_body_bytes=9,
        )
        with pytest.raises(HttpLimitError, match="resource_limit") as caught:
            client.text(URL, cache_token=1)
        assert caught.value.resource == "logical_requests"
    assert len(calls) == 2


def test_cache_delivery_limit(native):
    client, calls, _ = native([reply(b"abc")], cache_ttl=60, max_download_bytes=5)
    with client:
        client.text(URL)
        with pytest.raises(HttpLimitError) as caught:
            client.text(URL)
        assert caught.value.resource == "delivered_body_bytes"
        assert caught.value.observed == 6
        assert client.counters.received_body_bytes == 3
    assert len(calls) == 1


def test_explicit_shared_cache_and_ttl_expiry(native):
    now = [0]
    cache = ResponseCache(clock=lambda: now[0])
    first, calls, _ = native([reply(b"first")], cache=cache, cache_ttl=10)
    with first:
        assert first.text(URL) == "first"
    second, second_calls, _ = native([reply(b"later")], cache=cache, cache_ttl=10)
    with second:
        assert second.text(URL) == "first"
        now[0] = 10
        assert second.text(URL) == "later"
    assert len(calls) == len(second_calls) == 1


def test_default_cache_is_off_and_cookie_state_is_isolated(native):
    client, calls, _ = native(
        [reply(headers={"Set-Cookie": "session=secret"}), reply()]
    )
    with client:
        client.json(URL)
        client.json(URL)
    assert len(calls) == 2 and all("cookie" not in call.headers for call in calls)


def test_redirect_counts_attempt_and_does_not_read_or_cache_original(native):
    redirect = reply(b"unread", status=302, headers={"Location": "/final"})
    client, calls, _ = native([redirect, reply()], cache_ttl=60)
    with client:
        assert client.json(URL) == {"temperature": 12}
        assert client.counters.network_attempts == 2
        assert client.counters.logical_requests == 1
    assert redirect.test_body.reads == 0 and redirect.is_closed
    assert calls[1].url == "https://weather.test/final"


@pytest.mark.parametrize(
    "location", ["https://elsewhere.test/?secret=key", "http://weather.test/final"]
)
def test_redirect_policy(native, location):
    response = reply(status=302, headers={"Location": location})
    client, calls, _ = native([response])
    with client, pytest.raises(HttpTransportError) as caught:
        client.json(URL)
    assert caught.value.kind == "redirect_policy" and len(calls) == 1
    assert "secret" not in str(caught.value) and response.is_closed


def test_redirect_hop_limit(native):
    responses = [reply(status=302, headers={"Location": "/again"}) for _ in range(6)]
    client, calls, _ = native(responses)
    with client, pytest.raises(HttpTransportError, match="redirect_policy"):
        client.text(URL)
    assert len(calls) == 6 and all(response.is_closed for response in responses)


def test_status_body_bound_and_failed_body_keeps_status(native):
    response = reply([b"secret" * 2000, b"unread"], status=404)
    client, calls, delays = native([response])
    with client, pytest.raises(HttpTransportError) as caught:
        client.text(URL)
    assert caught.value.status == 404 and len(caught.value.error_body) == 8192
    assert "secret" not in str(caught.value) and response.test_body.reads == 1
    assert len(calls) == 1 and not delays


def test_error_body_budget_stops_retry_and_preserves_status(native):
    response = reply(b"x" * (32 * 1024), status=503, headers={"Retry-After": "30"})
    client, calls, delays = native(
        [response],
        max_download_bytes=10,
    )
    with client:
        with pytest.raises(HttpTransportError) as caught:
            client.text(URL)
        assert (
            caught.value.kind,
            caught.value.status,
            caught.value.retry_after,
            caught.value.attempts,
        ) == ("http_status", 503, 30.0, 1)
        assert len(calls) == 1 and delays == []
        assert client.counters.received_body_bytes == 8192
        assert response.is_closed
        with pytest.raises(HttpLimitError) as limit_caught:
            client.text(URL)
        assert limit_caught.value.resource == "received_body_bytes"


@pytest.mark.parametrize("backend", ["native", "sync", "async"])
def test_large_status_body_leaves_budget_for_retry(backend, server):
    from dtcc_core.io import Budget, RetryPolicy

    url, calls, _ = server
    request = Request("GET", url + "/retry", params={"error_size": 32 * 1024})

    def check(response, counters):
        assert response.status == 200 and response.body == b"ok"
        assert calls["/retry"] == 2
        assert counters.logical_requests == 1
        assert counters.network_attempts == 2
        assert counters.received_body_bytes == 8192 + len(b"ok")
        assert counters.delivered_body_bytes == len(b"ok")

    if backend == "native":
        with http_client(max_download_bytes=12 * 1024, retries=1) as client:
            check(client.request(request.url, params=request.params), client.counters)
    else:

        async def scenario(api):
            check(await api.request(request), api.client.counters)

        test_http_transport.run(
            backend,
            scenario,
            budget=Budget(max_received_bytes=12 * 1024),
            retry=RetryPolicy(retries=1),
        )


def test_certificate_error_not_retried(native):
    failure = httpx2.ConnectError("secret")
    failure.__cause__ = ssl.SSLCertVerificationError(1, "secret")
    client, calls, delays = native([failure])
    with client, pytest.raises(HttpTransportError) as caught:
        client.text(URL)
    assert caught.value.kind == "certificate" and len(calls) == 1 and not delays


def test_local_exception_closes_everything_and_propagates(native):
    failure = RuntimeError("local failure")
    response = reply([b"partial", failure])
    client, _, _ = native([response])
    with pytest.raises(RuntimeError) as caught, client:
        client.text(URL)
    assert caught.value is failure
    assert response.is_closed and client._backend._client.is_closed
    client.close()
    with pytest.raises(ValueError, match="closed"):
        client.text(URL)
    with pytest.raises(ValueError, match="closed"), client:
        pass


@pytest.mark.parametrize(
    "env,explicit,directory,expected",
    [
        ({}, None, False, "default.pem"),
        (
            {"REQUESTS_CA_BUNDLE": "requests.pem", "CURL_CA_BUNDLE": "curl.pem"},
            None,
            False,
            "requests.pem",
        ),
        ({"CURL_CA_BUNDLE": "curl.pem"}, None, False, "curl.pem"),
        (
            {"REQUESTS_CA_BUNDLE": "requests.pem"},
            Path("explicit.pem"),
            False,
            Path("explicit.pem"),
        ),
        ({"REQUESTS_CA_BUNDLE": "certs"}, None, True, "certs"),
    ],
)
def test_native_ca_configuration(monkeypatch, env, explicit, directory, expected):
    from dtcc_core.io._http import _httpx

    for name in ("REQUESTS_CA_BUNDLE", "CURL_CA_BUNDLE"):
        monkeypatch.delenv(name, raising=False)
    for name, value in env.items():
        monkeypatch.setenv(name, value)
    captured = {}
    monkeypatch.setattr(_httpx, "where", lambda: "default.pem")
    monkeypatch.setattr(_httpx.os.path, "isdir", lambda path: directory)
    monkeypatch.setattr(
        _httpx.ssl, "create_default_context", lambda **kwargs: captured.update(kwargs)
    )
    monkeypatch.setattr(
        _httpx.httpx2, "Client", lambda **kwargs: captured.update(kwargs)
    )
    with http_client(
        ca_bundle=explicit, trust_env=True, proxy="http://proxy.test"
    ) as client:
        client._backend._build_client()
    assert captured["capath" if directory else "cafile"] == expected
    assert captured["trust_env"] and captured["proxy"] == "http://proxy.test"
    assert not captured["follow_redirects"]


def test_native_loopback_io_and_no_requests_delegation(server, monkeypatch):
    url, calls, _ = server

    def forbidden(*args, **kwargs):
        pytest.fail("Native backend delegated to Requests")

    monkeypatch.setattr(HttpClient, "__init__", forbidden)
    with http_client(max_download_bytes=1000, respect_retry_after=False) as http:
        assert http.json(url + "/json") == [1, None, "text"]
        assert http.text(url + "/retry") == "ok"
        with pytest.raises(HttpLimitError):
            http.text(url + "/gzip")
    assert calls["/retry"] == 2
    assert http.counters.received_body_bytes >= 200000


def test_native_real_partial_body_retry(server):
    url, calls, _ = server
    with http_client() as http:
        assert http.text(url + "/broken") == "ok"
        assert http.counters.received_body_bytes > 2
        assert http.counters.delivered_body_bytes == 2
    assert calls["/broken"] == 2


def test_native_decoded_gzip_limit(native):
    response = reply(gzip.compress(b"x" * 100000), headers={"Content-Encoding": "gzip"})
    client, _, _ = native([response], max_download_bytes=1000)
    with client, pytest.raises(HttpLimitError):
        client.text(URL)
    assert client.counters.received_body_bytes == 100000 and response.is_closed


def test_low_level_logical_budget_uses_structured_error():
    from dtcc_core.io._http import Budget
    from dtcc_core.io._http._policy import Operation, State

    state = State(budget=Budget(max_requests=0))
    with pytest.raises(HttpLimitError) as caught:
        Operation(state, Request("GET", URL), None)
    assert caught.value.resource == "logical_requests"
    assert caught.value.observed == 1 and state.counters.network_attempts == 0


def test_native_tls_trusted_bundle_and_untrusted_certificate(monkeypatch):
    import threading
    from collections import Counter
    from http.server import ThreadingHTTPServer

    from .test_http_transport import Handler

    fixtures = Path(__file__).with_name("fixtures")
    cert = fixtures / "http-test-cert.pem"
    trust = ssl.SSLContext(ssl.PROTOCOL_TLS_SERVER)
    trust.load_cert_chain(cert, fixtures / "http-test-key.pem")
    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    server.calls, server.received = Counter(), []
    server.socket = trust.wrap_socket(server.socket, server_side=True)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    for name in ("REQUESTS_CA_BUNDLE", "CURL_CA_BUNDLE"):
        monkeypatch.delenv(name, raising=False)
    try:
        url = f"https://localhost:{server.server_port}/json"
        with http_client() as http, pytest.raises(HttpTransportError) as caught:
            http.json(url)
        assert caught.value.kind == "certificate" and caught.value.attempts == 1
        with http_client(ca_bundle=cert) as http:
            assert http.json(url) == [1, None, "text"]
        monkeypatch.setenv("REQUESTS_CA_BUNDLE", str(cert))
        with http_client() as http:
            assert http.json(url) == [1, None, "text"]
    finally:
        server.shutdown()
        server.server_close()
        thread.join()


@pytest.mark.parametrize(
    "params,expected",
    [
        (None, "token=secret"),
        ({"day": "2026-01-01"}, "token=secret&day=2026-01-01"),
        ([("token", "second"), ("x", "a b")], "token=secret&token=second&x=a+b"),
    ],
)
def test_url_query_is_preserved_and_explicit_params_appended(native, params, expected):
    client, calls, _ = native([reply()])
    with client:
        client.json(URL, params=params)
    assert calls[0].url.query.decode() == expected


def test_redirect_location_retains_query(native):
    client, calls, _ = native(
        [reply(status=302, headers={"Location": "/next?version=123"}), reply()]
    )
    with client:
        client.json(URL, params={"day": "today"})
    assert str(calls[1].url) == "https://weather.test/next?version=123"


def test_shared_budget_counts_inflight_chunks_after_another_read_crosses_limit():
    from dtcc_core.io._http import Budget

    budget = Budget(max_received_bytes=10)
    request = Request("GET", URL)
    budget.consume("received_body_bytes", 8, request)
    for observed in (16, 24):
        with pytest.raises(HttpLimitError) as caught:
            budget.consume("received_body_bytes", 8, request)
        assert caught.value.resource == "received_body_bytes"
        assert caught.value.observed == observed


@pytest.mark.parametrize(
    "body, content_type, expected",
    [
        ({"json": {"region": "west"}}, "application/json", b'{"region": "west"}'),
        (
            {"form": [("mode", "bus"), ("mode", "tram")]},
            "application/x-www-form-urlencoded",
            b"mode=bus&mode=tram",
        ),
        ({"body": b"binary"}, None, b"binary"),
    ],
)
def test_buffered_post_serializes_body(native, body, content_type, expected):
    client, calls, _ = native([reply()])
    with client:
        assert client.json(URL, method="POST", **body) == {"temperature": 12}
    assert calls[0].method == "POST"
    assert calls[0].content == expected
    if content_type:
        assert calls[0].headers["content-type"] == content_type


@pytest.mark.parametrize("idempotent, attempts", [(False, 1), (True, 2)])
def test_post_replay_is_explicit(native, idempotent, attempts):
    client, calls, _ = native([reply(status=503), reply()], retries=1)
    with client:
        if idempotent:
            assert client.json(URL, method="POST", json={}, idempotent=True)
        else:
            with pytest.raises(HttpTransportError):
                client.json(URL, method="POST", json={})
    assert len(calls) == attempts


def test_binary_and_custom_timeout_and_redirect_policy(native):
    from dtcc_core.io import RedirectPolicy, Timeout

    client, calls, _ = native(
        [reply(b"\x00\xff")],
        timeout=Timeout(2, 7),
        redirects=RedirectPolicy(max_hops=0),
    )
    with client:
        assert client.bytes(URL) == b"\x00\xff"
        assert client._backend._state.timeout == Timeout(2, 7)
    assert len(calls) == 1
    with pytest.raises(ValueError, match="total"):
        http_client(timeout=Timeout(total=5))
    with pytest.raises(ValueError, match="redirects"):
        http_client(redirects=True)


def test_multipart_borrows_file_and_does_not_retry(native, tmp_path):
    from dtcc_core.io import MultipartPart

    path = tmp_path / "upload.bin"
    path.write_bytes(b"file payload")
    client, calls, _ = native([reply(status=503)], retries=2)
    with path.open("rb") as source, client:
        with pytest.raises(HttpTransportError):
            client.request(
                URL,
                method="POST",
                multipart=[MultipartPart("files", source, "upload.bin")],
            )
        assert not source.closed
    assert len(calls) == 1
    assert b'filename="upload.bin"' in calls[0].content
    assert b"file payload" in calls[0].content


@pytest.mark.parametrize("source_kind", ["memory", "file"])
@pytest.mark.parametrize("offset", [0, 15, 27])
def test_multipart_reads_from_current_position(server, tmp_path, source_kind, offset):
    from email.parser import BytesParser
    from email.policy import default

    from dtcc_core.io import MultipartPart

    url, _, received = server
    contents = b"private prefix:file payload"
    if source_kind == "file":
        path = tmp_path / "upload.bin"
        path.write_bytes(contents)
        source = path.open("rb")
    else:
        source = io.BytesIO(contents)
    with source, http_client() as client:
        source.seek(offset)
        client.request(
            url, method="POST", multipart=[MultipartPart("file", source, "upload.bin")]
        )
        assert not source.closed
        assert source.tell() == len(contents)
    headers, body = received[0]
    message = BytesParser(policy=default).parsebytes(
        f"Content-Type: {headers['Content-Type']}\r\n\r\n".encode() + body
    )
    part = next(message.iter_parts())
    assert part.get_filename() == "upload.bin"
    assert part.get_payload(decode=True) == contents[offset:]


@pytest.mark.parametrize("status", [301, 302, 303, 307, 308])
@pytest.mark.parametrize("method", ["GET", "POST"])
def test_mutating_json_during_send_preserves_redirect_error(native, status, method):
    payload = {"value": [1]}

    def events():
        payload["value"].append(object())
        yield reply(status=status, headers={"Location": "/next?token=private"})

    client, calls, _ = native(events())
    with client, pytest.raises(HttpTransportError) as caught:
        client.request(URL, method=method, json=payload)
    error = caught.value
    assert (error.kind, error.status, error.phase) == (
        "redirect_policy",
        status,
        "read",
    )
    assert error.redirect_target == "https://weather.test/next"
    assert error.attempts == len(calls) == 1
    assert calls[0].content == b'{"value": [1]}'
