# Internal HTTP transport

`dtcc_core.io.data._http` provides synchronous `HttpClient` and asynchronous
`AsyncHttpClient` operations. It is an internal interface, with no re-export from
Core or I/O. Existing provider clients continue to use their own implementations.

The transport handles HTTP, response lifetime, retries, decoded-body budgets and
optional public-response caching. Providers retain discovery, schema validation,
units, dataset health, endpoint failover and error translation. Downloaders retain
concurrency scheduling, progress callbacks, temporary files, whole-file restarts,
content validation and atomic publication.

## Buffered requests

```python
from dtcc_core.io.data._http import HttpClient, Request, RetryPolicy, Timeout

with HttpClient(timeout=Timeout(connect=8, read=30),
                retry=RetryPolicy(retries=2)) as client:
    response = client.get(
        "https://service.example/observations",
        params={"parameter": "temperature", "year": 2025},
        expected_status={200},
    )
    data = response.json()
    print(client.counters.network_attempts)
```

`get(url, **options)` and `post(url, **options)` are thin buffered conveniences
over `request(Request(...))`, available on both clients (await the async methods).
Options are `Request` fields plus `max_body_bytes`. For example,
`client.post(url, json=payload, headers={"Accept": "application/json"})` returns a
`Response`; call `.json()` separately. `json=None` sends JSON null. POST retries
still require explicit `idempotent=True`. Use `Request` directly for other methods
and streams. There is one execution, validation, caching and accounting path.

The client constructors expose `timeout`, `retry`, `redirects`, `budget`, `cache`,
`ca_bundle` and `proxy` as named options. `sleep`, `wall_clock` and `jitter_source`
are also explicit for deterministic tests. No backend keyword arguments are
forwarded through the interface.

`request()` buffers a decoded body, with a default cap of 32 MiB. Override with
`max_body_bytes=...`; `None` removes the cap. `Response` contains `status`,
`headers` (lowercase names with repeated values combined), `target` (URL without credentials/query/fragment), and `body` bytes.
`json()` accepts any JSON value, including arrays and null; it does not interpret
provider error envelopes. `text(encoding="utf-8")` requires an explicit alternate
encoding when necessary, such as `utf-8-sig` for BOM removal.

Request query, form and header inputs accept mappings or sequences of pairs and
are snapshotted as immutable string pairs. Query/form values may be strings,
integers or finite floats; numbers use `str(value)`. Keys and header values must
be strings. Booleans, None, lists and nested objects are rejected: omit absent
parameters and encode provider-specific values explicitly. Query/form keys may
repeat when supplied as pairs, such as `params=(("parameter", "temperature"),
("parameter", "wind"))`, and retain their order. Request header names must be unique ignoring case;
duplicates are rejected before network or cache access. Set exactly one of `body` (bytes), `json` (including
None), `form`, or `multipart`. Multipart is a tuple of `MultipartPart(name, value,
filename, content_type)`, with bytes or caller-owned binary streams as values.
Repeated field names support multiple uploaded files. Requests may buffer multipart
uploads in memory; aiohttp streams file parts. This interface does not guarantee
bounded upload memory or enforce upload-byte budgets.

## Async streams

```python
from dtcc_core.io.data._http import AsyncHttpClient, Request, Budget

async def read_tile(url, consume_chunk):
    budget = Budget(max_attempts=10, max_received_bytes=100 * 1024 * 1024)
    async with AsyncHttpClient(budget=budget) as client:
        async with client.stream(Request("GET", url, expected_status={200})) as response:
            # Inspect response.status and response.headers before consuming the body.
            async for chunk in response.iter_bytes(1024 * 1024):
                consume_chunk(chunk)
```

Sync streaming uses `with client.stream(...)` and `for chunk in
response.iter_bytes(...)`. Streams have no default body cap; pass
`max_body_bytes` to set one. Iteration is single-use. Close the context even if
iteration stops early. Responses are also closed on exceptions and async
cancellation. Client-created sessions are reused and closed with their client;
upload files remain caller-owned. Both clients always own their sessions. Sync clients are
not thread-safe. Async clients belong to one event loop. Share a `Budget` or
`ResponseCache` explicitly when a batch needs aggregate state.

Reuse a client across a batch to reuse its connections. Keep scheduling and
provider validation outside the client:

```python
async def read_tiles(urls, consume_chunk):
    async with AsyncHttpClient() as client:
        for url in urls:
            async with client.stream(Request("GET", url, expected_status={200})) as response:
                async for chunk in response.iter_bytes(1024 * 1024):
                    consume_chunk(chunk)
```

## Retry, redirect and timeout policy

Retries default to zero. `RetryPolicy(retries=N)` allows N additional attempts
for GET/HEAD on connection failures, timeouts and statuses 429/500/502/503/504.
`connect` and `read` switches distinguish connection setup and response reads;
requests waiting for headers can report a read timeout, which uses the read
switch. Other connection failures before receiving headers use the connect
switch. Statuses, backoff, maximum delay, jitter and `respect_retry_after` are
configurable. Retry-After supports seconds and HTTP dates and is capped by
`max_delay`. Tests may inject the sleeper, wall clock and jitter source.

A body-bearing operation such as POST requires `Request(idempotent=True)` to
retry; that assertion belongs to the caller. Upload streams are never replayed:
with retries enabled, a GET/HEAD or idempotent request carrying one raises
`ValueError` before sending, and other requests carrying one are sent once.
Buffered body reads can retry before delivery and retain accounting for
discarded chunks. Stream read failures propagate without
replay, including failures before the first chunk. File restart belongs to the
caller. Adapter-level and aiohttp implicit connection retries are disabled so
there is one retry loop.

Redirects default to rejection. `RedirectPolicy(max_hops=N)` allows same-origin
GET/HEAD redirects without request bodies. `allowed_hosts`, when supplied, applies
to every destination; an initial request to another host raises `ValueError`
before sending. Hosts and origins compare case-insensitively, and an origin's
default port may be written or omitted. Cross-origin redirects additionally
require the exact destination origin in `allowed_origins` and a caller assertion
`Request(public=True)`. They are rejected for requests with Authorization/Cookie headers. Only Accept, Accept-Encoding and User-Agent headers
are forwarded across origins. Every hop consumes an attempt; HTTPS downgrades,
userinfo and non-HTTP schemes are rejected. Provider host restrictions are caller
configuration.

Rejected redirects retain their HTTP `status` and, when the destination is a
valid HTTP(S) URL, `redirect_target` with credentials/query/fragment removed.
Malformed or non-HTTP destinations have `redirect_target=None`. Rejected redirect
bodies remain unread and no raw Location is exposed. Including a redirect status
in `expected_status` does not bypass redirect policy. The transport therefore does
not reproduce provider diagnostics that included a redirect's raw Location/body.

`Timeout` requires finite positive connect and read inactivity timeouts; it is
not a complete operation deadline. Async `total` limits each backend attempt,
not retry sleeps or a whole batch. Sync clients reject `total`, because requests
cannot implement the same deadline reliably.

Owned sessions disable environment proxies and implicit netrc authentication.
Sync TLS uses requests' normal CA bundle, including REQUESTS_CA_BUNDLE and
CURL_CA_BUNDLE overrides. Async TLS uses the existing `download_ssl_context()`
helper, including SSL_CERT_FILE/SSL_CERT_DIR. TLS verification stays enabled.
Use `ca_bundle="/path/to/ca.pem"` (a file or a directory) to replace the default
trusted CAs in both clients, and `proxy="http://proxy.example:8080"` to route
requests through a proxy. Credentials
belong in explicit request headers; challenge-based authentication is unsupported.
Arbitrary backend session injection is not supported. The clients construct their
own sessions, so ownership, cookies, retries and verification have one policy.
Tests instrument those owned sessions without exposing an injection API.

The async adapter disables aiohttp's implicit connection replay through the
isolated `_retry_connection` switch. Backend conformance tests exercise the locked
dependency versions; dependency upgrades must retain exact attempt accounting.

## Limits and diagnostics

`Budget` optionally limits attempts, received bytes and delivered bytes. Counters
are client-local; shared budgets reserve attempts atomically across async tasks.
Counters distinguish:

| Counter | Meaning |
| --- | --- |
| logical_requests | Calls to request/stream after validation |
| network_attempts | Backend send invocations, including failed connects, retries and redirects |
| cache_hits | Buffered calls served by the explicit cache |
| received_body_bytes | Decoded chunks returned by the backend, including discarded/failed reads |
| delivered_body_bytes | Complete buffered bytes or streamed chunks delivered to the caller |

Cache hits consume delivery budgets but do not consume attempts or received-byte
budgets. Once a byte cap has been reached, new network sends are rejected. Crossing
or failing a byte budget also blocks further reads by operations sharing it.
An exact-boundary response may still finish its iteration to establish EOF;
concurrent reads already in flight can return data and are counted. Recovery uses
a new Budget; there is no implicit reset. Limits include the chunk that crosses a received-byte cap; this chunk is
not delivered. These are decoded-body limits, not exact wire-byte measurements or
peak-memory guarantees: backends can buffer/decompress before yielding. Error
bodies are read only up to 8 KiB. Unread redirect bodies are closed without being
counted. Observation counts, time intervals and provider-specific budgets remain
outside this package.

`HttpTransportError` exposes `kind`, `status`, `attempts`, `retry_allowed_by_policy`, `method`
and sanitized `target`. Kinds include timeout, connection, request, http_status,
redirect_policy, decoding and resource_limit. `.metadata()` returns only safe
diagnostics. `.error_body` is bounded raw provider data: callers must sanitize it
before exposing it. Body, headers, query values and upload content are excluded
from request/response repr; backend exception chains are suppressed. URLs may
still contain sensitive path segments, so callers must not place credentials in
paths. Local reader exceptions retain their types even when response headers arrive
before the reader fails. Callback exceptions and cancellation also retain their
own types. A complete early server response can stop an unfinished async upload,
following aiohttp behavior; the client does not force the server to consume it.

Dataset adapters can translate these neutral errors to their own error classes;
the transport never imports dataset types. No provider-specific status or payload
interpretation is built in.

`retry_allowed_by_policy` describes whether the configured transport policy
permits replay of this kind of failure for this request, before considering
already-used retries. It is false with retries disabled, even for 429, and can be
true after retry exhaustion. It does not decide an outer downloader's policy.

Errors while rejecting a response also expose `retry_after`: a nonnegative
duration in seconds parsed from Retry-After seconds or an HTTP date, or None for
missing/invalid values. Past dates yield zero. It is available even with retries
disabled or `respect_retry_after=False`, and is not capped by `RetryPolicy.max_delay`.
Automatic transport sleeps remain capped; an outer retry loop must apply its own
cap and fallback. `retry_after`, `redirect_target` and raw error bytes are explicit
attributes, excluded from `.metadata()` and the exception message.

For whole-file restarts, keep one retry loop in the downloader and disable
transport retries to avoid multiplying attempts. This example restarts into a
fresh temporary file, leaves local exceptions unchanged, and limits outer waits:

```python
import asyncio
import os
import tempfile
from pathlib import Path
from dtcc_core.io.data._http import AsyncHttpClient, HttpTransportError, Request

async def download_file(client, url, destination):
    destination = Path(destination)
    for attempt in range(3):
        temporary = None
        try:
            async with client.stream(Request("GET", url, expected_status={200})) as response:
                fd, name = tempfile.mkstemp(dir=destination.parent, suffix=".part")
                temporary = Path(name)
                with os.fdopen(fd, "wb") as output:
                    async for chunk in response.iter_bytes(1024 * 1024):
                        output.write(chunk)
            # Provider-specific file validation belongs before publication.
            os.replace(temporary, destination)
            return destination
        except HttpTransportError as error:
            transient = error.kind in {"connection", "timeout"} or (
                error.kind == "http_status" and error.status in {429, 500, 502, 503, 504}
            )
            if not transient or attempt == 2:
                raise
            delay = error.retry_after if error.retry_after is not None else 2 ** attempt
        finally:
            if temporary is not None:
                try:
                    temporary.unlink(missing_ok=True)
                except OSError:
                    pass
        await asyncio.sleep(min(30, delay))

async def download_batch(urls_and_paths):
    async with AsyncHttpClient() as client:  # retries default to zero
        for url, path in urls_and_paths:
            await download_file(client, url, path)
```

The destination directories must already exist. A shared byte budget, if supplied,
also counts discarded downloads; resource-limit errors are not retried here.

## Explicit response cache

Pass `cache=ResponseCache(...)` and opt individual requests in with
`public=True, cache_ttl=...`. The public assertion includes the URL and query:
callers must leave it false for API keys or other query credentials. Authorization,
Cookie and Proxy-Authorization requests are always bypassed. Owned sessions do not
retain response cookies.

The thread-safe TTL/LRU cache limits entries and retained body bytes, defaults to
32 entries/32 MiB, and has `clear()` and `invalidate(request)`. It caches only
complete buffered, bodyless GETs with status 200. Keys hash URL, ordered parameters,
request headers, namespace and optional `cache_token`. Tokens support provider
version/update markers. A shorter requested TTL constrains existing entries too.

Cache-Control is case-insensitive. Response no-store/private/no-cache, Vary and
Set-Cookie bypass storage. Request no-store/no-cache/max-age/Pragma bypass lookup
and storage. Server max-age or Expires freshness (with Age and apparent age from Date) caps the caller TTL.
No disk persistence, stale fallback, conditional revalidation, stream/POST caching
or concurrent-request coalescing is provided. Redirected responses are not cached.

`request()` caches raw successful bytes; later calls to `Response.json()` or
`.text()` do not invalidate them on decoding failure. Schema validation is outside
transport. Callers requiring validated cache entries must decode before explicit
`ResponseCache.put()` or invalidate a rejected response.

## Verification

The tests in `tests/io/test_http_transport.py` exercise both backends against a
loopback server, including multipart, truncation, compressed limits and caching.
`tests/io/test_http_policy.py` covers deterministic decisions, budgets, ownership
and cancellation. No real provider data or public-network access is required.
