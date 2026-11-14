"""Deterministic HTTP policy, ownership and failure-path checks."""

import asyncio
import io
import ssl
import time
from dataclasses import replace
from email.utils import formatdate
from pathlib import Path

import pytest
import requests
from requests.adapters import HTTPAdapter

from dtcc_core.io.data._http import (
    AsyncHttpClient,
    Budget,
    HttpClient,
    HttpTransportError,
    MultipartPart,
    RedirectPolicy,
    Request,
    Response,
    ResponseCache,
    RetryPolicy,
    Timeout,
)
from dtcc_core.io.data._http._policy import (
    Operation,
    State,
    parse_retry_after,
    redirect,
    retry_delay,
)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"url": "file:///etc/passwd"},
        {"url": "https://user:secret@host/"},
        {"url": "https://host:bad/"},
        {"method": "GET\n"},
        {"body": b"x", "json": {}},
        {"body": "text"},
        {"cache_ttl": float("nan")},
        {"params": {"token": ["secret"]}},
        {"params": {"token": {"secret": 1}}},
        {"params": {"flag": True}},
        {"params": {"missing": None}},
        {"params": {"n": float("inf")}},
        {"form": {"n": float("nan")}},
        {"headers": {"X-Number": 1}},
        {"headers": {"Accept": "one", "accept": "two"}},
        {"params": [(1, "secret")]},
        {"params": ["ab"]},
        {"form": [("secret",)]},
    ],
)
def test_request_validation(kwargs):
    options = dict(method="GET", url="https://example.test/")
    options.update(kwargs)
    with pytest.raises(ValueError) as error:
        Request(**options)
    assert "secret" not in str(error.value)


def test_request_snapshots_mapping_and_pair_inputs():
    params = {"year": 2025, "resolution": 0.5}
    headers = {"Accept": "application/json"}
    form = [["data", "first"], ["data", "second"]]
    request = Request(
        "POST", "https://example.test/", params=params, headers=headers, form=form
    )
    params["year"] = 2026
    headers["Accept"] = "text/plain"
    form[0][1] = "changed"
    assert request.params == (("year", "2025"), ("resolution", "0.5"))
    assert request.headers == (("Accept", "application/json"),)
    assert request.form == (("data", "first"), ("data", "second"))


@pytest.mark.parametrize(
    "value, expected",
    [
        (None, None),
        ("", None),
        (" 90 ", 90.0),
        ("0", 0.0),
        (formatdate(105, usegmt=True), 5.0),
        (formatdate(95, usegmt=True), 0.0),
        ("Thu Jan  1 00:01:45 1970", 5.0),
        ("invalid", None),
        ("-1", None),
        ("1.5", None),
        ("inf", None),
        ("9" * 400, None),
    ],
)
def test_parse_retry_after(value, expected):
    assert parse_retry_after(value, lambda: 100) == expected


@pytest.mark.parametrize(
    "factory",
    [
        lambda: Timeout(connect=0),
        lambda: Timeout(read=float("inf")),
        lambda: RetryPolicy(retries=-1),
        lambda: RetryPolicy(backoff=float("nan")),
        lambda: RedirectPolicy(max_hops=-1),
        lambda: Budget(max_attempts=-1),
        lambda: ResponseCache(max_bytes=-1),
        lambda: HttpClient(timeout=Timeout(total=1)),
    ],
)
def test_policy_validation(factory):
    with pytest.raises(ValueError):
        factory()


def test_delay_backoff_retry_after_and_caps():
    policy = RetryPolicy(retries=4, backoff=2, max_delay=10, jitter=1)
    assert retry_delay(policy, 0, {}, lambda: 100, lambda: 0.5) == 2.5
    assert retry_delay(policy, 2, {}, lambda: 100, lambda: 0) == 8
    assert retry_delay(policy, 3, {}, lambda: 100, lambda: 0) == 10
    assert retry_delay(policy, 0, {"retry-after": "99"}, lambda: 100, lambda: 0) == 10
    assert (
        retry_delay(
            policy,
            0,
            {"Retry-After": formatdate(105, usegmt=True)},
            lambda: 100,
            lambda: 0,
        )
        == 5
    )
    assert retry_delay(policy, 0, {"Retry-After": "bad"}, lambda: 100, lambda: 0) == 2
    assert (
        retry_delay(
            replace(policy, respect_retry_after=False),
            0,
            {"Retry-After": "9"},
            lambda: 0,
            lambda: 0,
        )
        == 2
    )


def test_redirect_credentials_hosts_and_downgrades():
    request = Request(
        "GET",
        "https://one.test/path",
        public=True,
        headers=(("Accept", "application/json"), ("X-Api-Key", "secret")),
    )
    policy = RedirectPolicy(max_hops=2, allowed_origins={"https://two.test"})
    result = redirect(request, "https://two.test/path", policy, 0)
    assert result.headers == (("Accept", "application/json"),)
    for target in (
        "http://one.test/path",
        "https://three.test/",
        "file:///etc/passwd",
        "https://user:secret@one.test/",
    ):
        with pytest.raises(HttpTransportError, match="redirect_policy"):
            redirect(request, target, policy, 0)
    for changed in (
        replace(request, public=False),
        replace(request, headers=(("Authorization", "secret"),)),
        replace(request, method="POST", idempotent=True),
        replace(request, json={}),
    ):
        with pytest.raises(HttpTransportError):
            redirect(changed, "https://two.test/path", policy, 0)
    with pytest.raises(HttpTransportError):
        redirect(request, "/path", replace(policy, allowed_hosts={"other.test"}), 0)
    with pytest.raises(HttpTransportError):
        redirect(request, "/path", policy, 2)


class FakeResponse:
    status_code = 200
    headers = {}

    def __init__(self, chunks=(), status=200):
        self.chunks, self.closed, self.status_code = chunks, False, status

    def iter_content(self, chunk_size):
        for chunk in self.chunks:
            if isinstance(chunk, BaseException):
                raise chunk
            yield chunk

    def close(self):
        self.closed = True


def test_sync_session_ownership_cleanup_and_failed_read_accounting(monkeypatch):
    session = requests.Session()
    closed = []
    monkeypatch.setattr(session, "close", lambda: closed.append(True))
    first = FakeResponse([b"partial", requests.ConnectionError("secret")])
    second = FakeResponse([b"ok"])
    responses = iter([first, second])
    monkeypatch.setattr(
        session.get_adapter("https://"), "send", lambda *args, **kwargs: next(responses)
    )
    monkeypatch.setattr(requests, "Session", lambda: session)
    with HttpClient(retry=RetryPolicy(retries=1, backoff=0)) as client:
        response = client.request(Request("GET", "https://example.test/?key=secret"))
        assert response.body == b"ok"
        assert client.counters.received_body_bytes == 9
        assert client.counters.delivered_body_bytes == 2
    assert first.closed and second.closed and closed
    with pytest.raises(ValueError, match="closed"):
        client.request(Request("GET", "https://example.test/"))


@pytest.mark.parametrize("retry_connect, retries", [(False, 1), (True, 2)])
def test_connection_retry_policy(monkeypatch, retry_connect, retries):
    session, calls = requests.Session(), []

    def fail(*args, **kwargs):
        calls.append(True)
        raise requests.ConnectTimeout("secret")

    monkeypatch.setattr(session.get_adapter("https://"), "send", fail)
    monkeypatch.setattr(requests, "Session", lambda: session)
    with HttpClient(
        retry=RetryPolicy(retries=1, connect=retry_connect, backoff=0)
    ) as client:
        with pytest.raises(HttpTransportError) as error:
            client.request(Request("GET", "https://example.test/"))
        assert error.value.kind == "timeout"
        assert len(calls) == retries
    session.close()


@pytest.mark.parametrize("idempotent, count", [(False, 1), (True, 2)])
def test_post_retries_require_idempotency(monkeypatch, idempotent, count):
    session, calls = requests.Session(), []

    def fail(*args, **kwargs):
        calls.append(True)
        raise requests.ConnectionError()

    monkeypatch.setattr(session.get_adapter("https://"), "send", fail)
    monkeypatch.setattr(requests, "Session", lambda: session)
    with HttpClient(retry=RetryPolicy(retries=1, backoff=0)) as client:
        with pytest.raises(HttpTransportError):
            client.request(
                Request(
                    "POST",
                    "https://example.test/",
                    form=(("query", "x"),),
                    idempotent=idempotent,
                )
            )
        assert len(calls) == count
    session.close()


def test_no_replay_of_file_uploads():
    part = MultipartPart("file", io.BytesIO(b"x"), "x.bin")
    with HttpClient(retry=RetryPolicy(retries=1)) as client:
        with pytest.raises(ValueError, match="upload streams"):
            client.request(
                Request(
                    "POST",
                    "https://example.test/",
                    multipart=(part,),
                    idempotent=True,
                )
            )
    assert not part.value.closed


def test_sync_caller_oserror_does_not_become_transport_failure(monkeypatch):
    session = requests.Session()
    response = FakeResponse([OSError("disk")])
    monkeypatch.setattr(
        session.get_adapter("https://"), "send", lambda *args, **kwargs: response
    )
    monkeypatch.setattr(requests, "Session", lambda: session)
    with HttpClient() as client, pytest.raises(OSError, match="disk"):
        client.request(Request("GET", "https://example.test/"))
    assert response.closed
    session.close()


def test_cache_age_short_ttl_byte_eviction_and_keys():
    now = [0]
    cache = ResponseCache(max_bytes=3, clock=lambda: now[0], wall_clock=lambda: 100)
    request = Request(
        "GET",
        "https://example.test/",
        public=True,
        cache_ttl=60,
        params=(("x", "1"), ("x", "2")),
    )
    response = Response(
        200,
        request.target,
        {"Cache-Control": "max-age=10", "Age": "8"},
        b"abc",
        request,
    )
    cache.put(request, response)
    assert cache.get(request).body == b"abc"
    assert cache.get(replace(request, params=tuple(reversed(request.params)))) is None
    assert cache.get(replace(request, headers=(("Accept", "text/plain"),))) is None
    assert cache.get(replace(request, cache_namespace="other")) is None
    now[0] = 2
    assert cache.get(request) is None
    cache.put(request, replace(response, headers={}))
    other = replace(request, cache_token="other")
    cache.put(other, response)
    assert cache.get(request) is None
    assert cache.get(other) is not None
    cache.put(replace(other, cache_token="oversized"), replace(response, body=b"1234"))
    assert cache.get(other) is not None


def test_budget_delivery_limit_on_cache_hit():
    request = Request("GET", "https://example.test/", public=True, cache_ttl=60)
    cache = ResponseCache()
    cache.put(request, Response(200, request.target, {}, b"abc", request))
    with HttpClient(cache=cache, budget=Budget(max_delivered_bytes=2)) as client:
        with pytest.raises(HttpTransportError, match="resource_limit"):
            client.request(request)
        assert client.counters.network_attempts == 0


def test_async_cancellation_closes_response_and_owned_session(monkeypatch):

    async def scenario():
        entered = asyncio.Event()

        class Content:
            def iter_chunked(self, size):
                return self

            def __aiter__(self):
                return self

            async def __anext__(self):
                entered.set()
                await asyncio.Future()

        class AsyncResponse:
            status, headers, content = 200, {}, Content()
            closed = False

            def close(self):
                self.closed = True

        response = AsyncResponse()
        async with AsyncHttpClient() as client:
            session = client._session

            async def send(*args, **kwargs):
                return response

            monkeypatch.setattr(session, "request", send)
            task = asyncio.create_task(
                client.request(Request("GET", "https://example.test/"))
            )
            await entered.wait()
            task.cancel()
            with pytest.raises(asyncio.CancelledError):
                await task
            assert response.closed
        assert session.closed

    asyncio.run(scenario())


def test_async_shared_attempt_budget_is_atomic(monkeypatch):
    import aiohttp

    async def scenario():
        sent = []

        async def fail(*args, **kwargs):
            sent.append(True)
            await asyncio.sleep(0)
            raise aiohttp.ClientConnectionError("secret")

        async with AsyncHttpClient(budget=Budget(max_attempts=3)) as client:
            monkeypatch.setattr(client._session, "request", fail)
            results = await asyncio.gather(
                *(
                    client.request(Request("GET", "https://example.test/"))
                    for _ in range(10)
                ),
                return_exceptions=True,
            )
            assert len(sent) == 3
            assert sum(error.kind == "resource_limit" for error in results) == 7
            assert client.counters.network_attempts == 3

    asyncio.run(scenario())


def test_read_timeout_uses_read_policy_and_closes_response(monkeypatch):
    from urllib3.exceptions import ReadTimeoutError

    with requests.Session() as session:
        response = FakeResponse(
            [requests.ConnectionError(ReadTimeoutError(None, None, "secret"))]
        )
        calls = []

        def send(*args, **kwargs):
            calls.append(True)
            return response

        monkeypatch.setattr(session.get_adapter("https://"), "send", send)
        monkeypatch.setattr(requests, "Session", lambda: session)
        with HttpClient(retry=RetryPolicy(retries=1, read=False)) as client:
            with pytest.raises(HttpTransportError) as caught:
                client.request(Request("GET", "https://example.test/"))
            assert caught.value.kind == "timeout" and len(calls) == 1
            assert response.closed


def test_async_owned_session_tls_and_retry_configuration(monkeypatch):
    from dtcc_core.io.data import tls

    contexts = []
    real_context = tls.download_ssl_context

    def ssl_context():
        context = real_context()
        contexts.append(context)
        return context

    monkeypatch.setattr(tls, "download_ssl_context", ssl_context)

    async def scenario():
        async with AsyncHttpClient() as client:
            session = client._session
            assert session.connector._ssl is contexts[0]
            assert session._retry_connection is False
            assert not session.trust_env
        assert session.closed
        with pytest.raises(ValueError, match="closed"):
            await client.request(Request("GET", "https://example.test/"))

    asyncio.run(scenario())


def test_cache_expires_and_get_body_bypass():
    request = Request("GET", "https://example.test/", public=True, cache_ttl=60)
    now = [0]
    cache = ResponseCache(clock=lambda: now[0], wall_clock=lambda: 100)
    response = Response(
        200,
        request.target,
        {
            "Date": formatdate(100, usegmt=True),
            "Expires": formatdate(110, usegmt=True),
            "Age": "8",
        },
        b"ok",
        request,
    )
    cache.put(request, response)
    assert cache.get(request) is not None
    now[0] = 2
    assert cache.get(request) is None
    with_body = replace(request, body=b"x")
    cache.put(with_body, replace(response, headers={}))
    assert cache.get(with_body) is None


def test_initial_host_allowlist_is_checked_without_network():
    with HttpClient(redirects=RedirectPolicy(allowed_hosts={"allowed.test"})) as client:
        with pytest.raises(ValueError, match="allowed_hosts"):
            client.request(Request("GET", "https://other.test/"))
        assert client.counters.network_attempts == 0


def test_sync_redirect_body_is_never_eagerly_buffered(monkeypatch):
    class UnreadBody:
        closed = False
        reads = 0

        def stream(self, *args, **kwargs):
            self.reads += 1
            yield b"x" * 1000000

        def close(self):
            self.closed = True

        def release_conn(self):
            pass

    class RedirectAdapter(HTTPAdapter):
        def send(self, request, **kwargs):
            response = requests.Response()
            response.status_code = 302
            response.headers["Location"] = "/destination"
            response.raw = raw
            response.request = request
            response.url = request.url
            return response

    raw = UnreadBody()
    with requests.Session() as session:
        session.mount("https://", RedirectAdapter(max_retries=0))
        monkeypatch.setattr(requests, "Session", lambda: session)
        with HttpClient(budget=Budget(max_received_bytes=1)) as client:
            with pytest.raises(HttpTransportError, match="redirect_policy"):
                client.request(
                    Request("GET", "https://example.test/"), max_body_bytes=1
                )
            assert raw.reads == 0 and raw.closed
            assert client.counters.received_body_bytes == 0


@pytest.mark.parametrize(
    "freshness",
    [
        {"Cache-Control": "max-age=60"},
        {"Expires": formatdate(940, usegmt=True)},
    ],
)
def test_cache_does_not_refresh_already_stale_responses(freshness):
    request = Request("GET", "https://example.test/", public=True, cache_ttl=60)
    cache = ResponseCache(wall_clock=lambda: 1000)
    headers = {"Date": formatdate(880, usegmt=True), **freshness}
    cache.put(request, Response(200, request.target, headers, b"ok", request))
    assert cache.get(request) is None


@pytest.mark.parametrize("client_type", [HttpClient, AsyncHttpClient])
def test_session_injection_is_not_part_of_the_api(client_type):
    with pytest.raises(TypeError, match="session"):
        client_type(session=object())


@pytest.mark.parametrize("client_type", [HttpClient, AsyncHttpClient])
def test_ca_bundle_cannot_disable_tls_verification(client_type):
    with pytest.raises(TypeError):
        client_type(ca_bundle=False)


def test_sync_explicit_connection_options(monkeypatch, tmp_path):
    bundle = tmp_path / "ca.pem"
    with HttpClient(ca_bundle=bundle, proxy="http://proxy.test:8080") as client:
        captured = []

        def send(request, **kwargs):
            captured.append(kwargs)
            return FakeResponse([b"ok"])

        monkeypatch.setattr(client._session.get_adapter("https://"), "send", send)
        client.request(Request("GET", "https://example.test/"))
        assert captured[0]["verify"] == str(bundle)
        assert captured[0]["proxies"] == {
            "http": "http://proxy.test:8080",
            "https": "http://proxy.test:8080",
        }
        assert client._session.trust_env is False


def _one_ca_certificate(path):
    certificates = Path(requests.certs.where()).read_text()
    first = certificates.split("-----END CERTIFICATE-----")[0]
    path.write_text(first[first.index("-----BEGIN") :] + "-----END CERTIFICATE-----\n")
    return path


@pytest.mark.parametrize("as_directory", [False, True])
def test_async_ca_bundle_replaces_default_trust(monkeypatch, tmp_path, as_directory):
    import aiohttp

    bundle = _one_ca_certificate(tmp_path / "ca.pem")
    contexts = []
    real_connector = aiohttp.TCPConnector

    def connector(**kwargs):
        contexts.append(kwargs["ssl"])
        return real_connector(**kwargs)

    monkeypatch.setattr(aiohttp, "TCPConnector", connector)

    async def scenario():
        async with AsyncHttpClient(
            ca_bundle=tmp_path if as_directory else bundle,
            proxy="http://proxy.test:8080",
        ) as client:
            assert (
                client._kwargs(Request("GET", "https://example.test/"), [])["proxy"]
                == "http://proxy.test:8080"
            )

    asyncio.run(scenario())
    (context,) = contexts
    assert context.verify_mode == ssl.CERT_REQUIRED and context.check_hostname
    # A directory is loaded lazily, so it contributes no certificates up front.
    assert len(context.get_ca_certs()) == (0 if as_directory else 1)


def test_redirect_allowlists_ignore_case_and_default_ports():
    request = Request("GET", "https://data.example.org/path", public=True)
    policy = RedirectPolicy(
        max_hops=1,
        allowed_hosts={"Data.Example.org", "MIRROR.example.org"},
        allowed_origins={"HTTPS://Mirror.Example.org:443"},
    )
    Operation(State(redirects=policy), request, None)
    for target in ("https://Mirror.example.org/x", "https://mirror.example.org:443/x"):
        assert redirect(request, target, policy, 0).url == target
    with pytest.raises(HttpTransportError, match="redirect_policy"):
        redirect(request, "https://mirror.example.org:8443/x", policy, 0)


@pytest.mark.parametrize(
    "value",
    ["mirror.example.org", "https://mirror.example.org/data", "ftp://mirror.test"],
)
def test_redirect_origins_must_be_origins(value):
    with pytest.raises(ValueError, match="allowed_origins"):
        RedirectPolicy(allowed_origins={value})


def _cached_response(request, headers):
    return Response(200, request.target, headers, b"ok", request)


def test_uncacheable_response_does_not_evict_fresh_entries():
    cache = ResponseCache(max_entries=1)
    fresh = Request("GET", "https://example.test/a", public=True, cache_ttl=60)
    stale = replace(fresh, url="https://example.test/b")
    cache.put(fresh, _cached_response(fresh, {}))
    cache.put(stale, _cached_response(stale, {"Cache-Control": "max-age=0"}))
    assert cache.get(fresh) is not None
    assert cache.get(stale) is None
    # A newer response that may not be reused also drops the older entry.
    cache.put(fresh, _cached_response(fresh, {"Cache-Control": "max-age=0"}))
    assert cache.get(fresh) is None


@pytest.mark.skipif(not hasattr(time, "tzset"), reason="needs time.tzset")
@pytest.mark.parametrize("freshness", ["max-age", "expires"])
def test_cache_reads_unknown_zone_http_dates_as_utc(monkeypatch, freshness):
    request = Request("GET", "https://example.test/", public=True, cache_ttl=60)
    # usegmt=False writes the "-0000" zone, which parses to a naive datetime.
    headers = {"Date": formatdate(1000)}
    if freshness == "max-age":
        headers["Cache-Control"] = "max-age=60"
    else:
        headers["Expires"] = formatdate(1060)
    monkeypatch.setenv("TZ", "Etc/GMT-2")
    time.tzset()
    try:
        cache = ResponseCache(wall_clock=lambda: 1000)
        cache.put(request, _cached_response(request, headers))
        assert cache.get(request) is not None
    finally:
        monkeypatch.undo()
        time.tzset()


def test_async_connect_timeout_excludes_waiting_for_a_pooled_connection():
    client = AsyncHttpClient(timeout=Timeout(connect=3, read=7))
    timeout = client._kwargs(Request("GET", "https://example.test/"), [])["timeout"]
    assert (timeout.connect, timeout.sock_connect, timeout.sock_read) == (None, 3, 7)


@pytest.mark.parametrize(
    ("exception", "phase", "attempts"),
    [("total", "read", 1), ("connect", "connect", 2)],
)
def test_async_timeouts_before_headers_use_their_phase_switch(
    monkeypatch, exception, phase, attempts
):
    import aiohttp

    async def timeout(*args, **kwargs):
        raise (
            TimeoutError() if exception == "total" else aiohttp.ConnectionTimeoutError()
        )

    monkeypatch.setattr(aiohttp.ClientSession, "request", timeout)

    async def scenario():
        async with AsyncHttpClient(
            retry=RetryPolicy(retries=1, read=False, backoff=0),
            timeout=Timeout(total=5),
        ) as client:
            with pytest.raises(HttpTransportError) as caught:
                await client.request(Request("GET", "https://example.test/"))
        error = caught.value
        assert (error.kind, error.phase, error.attempts) == ("timeout", phase, attempts)

    asyncio.run(scenario())


def test_budget_errors_report_the_attempts_made(monkeypatch):
    responses = [FakeResponse([b"retry"], status=503), FakeResponse([b"ok"])]
    session = requests.Session()
    monkeypatch.setattr(
        session.get_adapter("https://"),
        "send",
        lambda *args, **kwargs: responses.pop(0),
    )
    monkeypatch.setattr(requests, "Session", lambda: session)
    with HttpClient(
        retry=RetryPolicy(retries=1, backoff=0),
        budget=Budget(max_received_bytes=6),
        sleep=lambda delay: None,
    ) as client:
        with pytest.raises(HttpTransportError) as caught:
            client.request(Request("GET", "https://example.test/"))
    error = caught.value
    assert (error.kind, error.attempts, error.phase) == ("resource_limit", 2, "read")
    session.close()


def _certificate_errors():
    import aiohttp
    from urllib3.exceptions import MaxRetryError
    from urllib3.exceptions import SSLError as Urllib3SSLError

    def sync_error(cause):
        # requests wraps urllib3's MaxRetryError, whose reason wraps the ssl error.
        return requests.exceptions.SSLError(
            MaxRetryError(None, "https://example.test/", reason=Urllib3SSLError(cause))
        )

    verify = ssl.SSLCertVerificationError(1, "certificate verify failed")
    return {
        "sync": (sync_error(verify), sync_error(ssl.SSLEOFError(8, "EOF"))),
        "async": (
            aiohttp.ClientConnectorCertificateError(object(), verify),
            aiohttp.ClientConnectorSSLError(object(), ssl.SSLEOFError(8, "EOF")),
        ),
    }


@pytest.mark.parametrize("backend", ["sync", "async"])
def test_certificate_failures_are_not_retried(monkeypatch, backend):
    import aiohttp

    certificate, transient = _certificate_errors()[backend]
    calls = []

    def fail(error):
        def send(*args, **kwargs):
            calls.append(error)
            raise error

        return send

    def outcome(error):
        calls.clear()
        if backend == "sync":
            session = requests.Session()
            monkeypatch.setattr(session.get_adapter("https://"), "send", fail(error))
            monkeypatch.setattr(requests, "Session", lambda: session)
            with HttpClient(
                retry=RetryPolicy(retries=2), sleep=lambda delay: None
            ) as client:
                with pytest.raises(HttpTransportError) as caught:
                    client.request(Request("GET", "https://example.test/"))
            return caught.value

        send = fail(error)

        async def request(*args, **kwargs):
            send()

        monkeypatch.setattr(aiohttp.ClientSession, "request", request)

        async def scenario():
            async def sleep(delay):
                pass

            async with AsyncHttpClient(
                retry=RetryPolicy(retries=2), sleep=sleep
            ) as client:
                with pytest.raises(HttpTransportError) as caught:
                    await client.request(Request("GET", "https://example.test/"))
            return caught.value

        return asyncio.run(scenario())

    error = outcome(certificate)
    assert (error.kind, error.attempts, len(calls)) == ("certificate", 1, 1)
    assert not error.retry_allowed_by_policy
    # Other TLS failures, such as a dropped handshake, remain retryable.
    error = outcome(transient)
    assert (error.kind, error.attempts, len(calls)) == ("connection", 3, 3)


def test_async_client_requires_the_aiohttp_replay_switch(monkeypatch):
    import aiohttp

    closed = []

    class Session:
        def __init__(self, **kwargs):
            pass

        async def close(self):
            closed.append(True)

    monkeypatch.setattr(aiohttp, "TCPConnector", lambda **kwargs: None)
    monkeypatch.setattr(aiohttp, "ClientSession", Session)

    async def scenario():
        with pytest.raises(RuntimeError, match="_retry_connection"):
            async with AsyncHttpClient():
                pass

    asyncio.run(scenario())
    assert closed == [True]


@pytest.mark.parametrize("backend", ["sync", "async"])
def test_failed_error_body_read_keeps_the_status_error(monkeypatch, backend):
    import aiohttp

    sends = []
    if backend == "sync":
        session = requests.Session()

        def send(*args, **kwargs):
            sends.append(1)
            return FakeResponse([b"part", requests.ConnectionError("reset")], 404)

        monkeypatch.setattr(session.get_adapter("https://"), "send", send)
        monkeypatch.setattr(requests, "Session", lambda: session)
        with HttpClient(
            retry=RetryPolicy(retries=2), sleep=lambda delay: None
        ) as client:
            with pytest.raises(HttpTransportError) as caught:
                client.request(Request("GET", "https://example.test/"))
    else:

        class Content:
            async def iter_chunked(self, size):
                yield b"part"
                raise aiohttp.ClientPayloadError("reset")

        class AsyncResponse:
            status, headers, content = 404, {}, Content()

            def close(self):
                pass

        async def request(*args, **kwargs):
            sends.append(1)
            return AsyncResponse()

        monkeypatch.setattr(aiohttp.ClientSession, "request", request)

        async def scenario():
            async def sleep(delay):
                pass

            async with AsyncHttpClient(
                retry=RetryPolicy(retries=2), sleep=sleep
            ) as client:
                with pytest.raises(HttpTransportError) as caught:
                    await client.request(Request("GET", "https://example.test/"))
            return caught

        caught = asyncio.run(scenario())
    error = caught.value
    assert (error.kind, error.status, error.phase) == ("http_status", 404, "read")
    assert (error.attempts, len(sends), error.error_body) == (1, 1, b"part")


def test_errors_after_headers_report_the_read_phase():
    state = State(redirects=RedirectPolicy(max_hops=1))
    operation = Operation(state, Request("GET", "https://example.test/"), None)
    with pytest.raises(HttpTransportError) as caught:
        operation.follow("https://other.test/", status=302)
    assert caught.value.phase == "read"
    cache = ResponseCache()
    request = Request("GET", "https://example.test/", public=True, cache_ttl=60)
    cache.put(request, Response(200, request.target, {}, b"cached", request))
    operation = Operation(State(cache=cache), request, 1)
    with pytest.raises(HttpTransportError) as caught:
        operation.state.cached(operation)
    # A cache hit performs no network I/O.
    assert caught.value.phase is None


def test_transport_errors_survive_pickling_and_copying():
    import copy
    import pickle

    request = Request("GET", "https://example.test/data?key=secret")
    error = HttpTransportError(
        request, "timeout", status=504, attempts=3, error_body=b"x", phase="read"
    )
    for clone in (pickle.loads(pickle.dumps(error)), copy.copy(error)):
        assert type(clone) is HttpTransportError and str(clone) == str(error)
        assert clone.metadata() == error.metadata()
        assert clone.error_body == b"x"
    assert error.metadata()["phase"] == "read"


def test_json_body_is_snapshotted_when_the_request_is_built():
    payload = {"value": 1}
    request = Request("POST", "https://example.test/", json=payload, idempotent=True)
    payload["value"] = 2
    with HttpClient() as client:
        assert client._kwargs(request)["data"] == b'{"value": 1}'
    assert AsyncHttpClient()._kwargs(request, [])["data"] == b'{"value": 1}'


def test_absent_headers_and_params_mean_none():
    request = Request("GET", "https://example.test/", headers=None, params=None)
    assert (request.headers, request.params) == ((), ())


@pytest.mark.parametrize("value", [float("nan"), float("inf")])
def test_json_body_rejects_non_finite_numbers(value):
    with pytest.raises(ValueError):
        Request("POST", "https://example.test/", json={"value": value})


def test_cache_hit_over_delivery_budget_has_no_phase():
    cache = ResponseCache()
    request = Request("GET", "https://example.test/", public=True, cache_ttl=60)
    cache.put(request, Response(200, request.target, {}, b"cached", request))
    state = State(cache=cache, budget=Budget(max_delivered_bytes=1))
    operation = Operation(state, request, None)
    with pytest.raises(HttpTransportError) as caught:
        state.cached(operation)
    assert (caught.value.kind, caught.value.phase) == ("resource_limit", None)
