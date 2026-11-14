"""Small synchronous JSON/text interface with owned connections and budgets."""

from collections.abc import Mapping

from .data._http import (
    Budget,
    RedirectPolicy,
    Request,
    ResponseCache,
    RetryPolicy,
    Timeout,
)
from .data._http._contracts import _pairs, _positive


class _HttpClient:
    """Serial JSON/text GETs; pass this client to helpers sharing an operation."""

    def __init__(self, backend, *, cache_ttl, cache_namespace):
        self._backend = backend
        self._cache_ttl, self._cache_namespace = cache_ttl, cache_namespace

    @property
    def counters(self):
        """Client-local requests, attempts, cache hits and decoded byte counts."""
        return self._backend.counters

    def __enter__(self):
        if self._backend._closed:
            raise ValueError("Client is closed")
        return self

    def __exit__(self, *exc):
        self.close()

    def close(self):
        """Close owned connections; repeated calls are harmless."""
        self._backend.close()

    def _get(
        self,
        url,
        *,
        params=None,
        headers=None,
        credentials=None,
        public=True,
        cache_token=None,
    ):
        query = _pairs(params or (), "params")
        if credentials is not None:
            if not isinstance(credentials, Mapping):
                raise ValueError("credentials must be a mapping of query credentials")
            secret = _pairs(credentials, "params")
            names = {name for name, _ in secret}
            if any(name in names for name, _ in query):
                raise ValueError("Credential names must not also occur in params")
            query += secret
            public = False
        return self._backend.request(
            Request(
                "GET",
                url,
                params=query,
                headers=headers,
                public=public,
                cache_ttl=self._cache_ttl,
                cache_namespace=self._cache_namespace,
                cache_token="" if cache_token is None else str(cache_token),
            )
        )

    def json(
        self,
        url,
        *,
        params=None,
        headers=None,
        credentials=None,
        public=True,
        cache_token=None,
    ):
        """Return a decoded JSON value; provider-specific validation is the caller's."""
        return self._get(
            url,
            params=params,
            headers=headers,
            credentials=credentials,
            public=public,
            cache_token=cache_token,
        ).json()

    def text(
        self,
        url,
        *,
        params=None,
        headers=None,
        credentials=None,
        public=True,
        cache_token=None,
        encoding="utf-8",
    ):
        """Return decoded text, using UTF-8 unless an encoding is supplied."""
        return self._get(
            url,
            params=params,
            headers=headers,
            credentials=credentials,
            public=public,
            cache_token=cache_token,
        ).text(encoding)


def http_client(
    *,
    timeout=30,
    retries=2,
    max_requests=2000,
    max_attempts=2000,
    max_download_bytes=100 * 1024 * 1024,
    cache_ttl=0,
    cache=None,
    cache_namespace="",
    respect_retry_after=True,
    ca_bundle=None,
    proxy=None,
    trust_env=False,
):
    """Create a serial context-managed JSON/text GET client using native HTTPX2.

    Requests, attempts and decoded byte budgets are shared across helper calls.
    Retries and redirects consume attempts; cache hits consume requests and
    delivered bytes. Received bytes include partial failed reads. Each byte
    budget is max_download_bytes. Cache is off until cache_ttl is positive;
    an explicitly supplied ResponseCache can share public responses across clients.

    Use credentials={"apikey": value} for query credentials to disable caching.
    Mark other private requests public=False. Errors are HttpTransportError;
    HttpLimitError additionally reports resource, limit and observed counts.
    Connections close on leaving the with block. Timeout is connect/read
    inactivity, not a total deadline. Environment proxies are opt-in (trust_env).
    """
    from .data._http._httpx import NativeHttpClient

    if isinstance(timeout, bool) or not isinstance(timeout, (int, float)):
        raise ValueError("timeout must be a positive number")
    _positive(timeout, "timeout")
    if type(retries) is not int:
        raise ValueError("retries must be a nonnegative integer")
    if isinstance(cache_ttl, bool) or not isinstance(cache_ttl, (int, float)):
        raise ValueError("cache_ttl must be a nonnegative number")
    _positive(cache_ttl, "cache_ttl", zero=True)
    if not all(type(value) is bool for value in (respect_retry_after, trust_env)):
        raise ValueError("respect_retry_after and trust_env must be booleans")
    if not isinstance(cache_namespace, str):
        raise ValueError("cache_namespace must be a string")
    if cache is not None and not isinstance(cache, ResponseCache):
        raise ValueError("cache must be a ResponseCache")
    budget = Budget(
        max_requests=max_requests,
        max_attempts=max_attempts,
        max_received_bytes=max_download_bytes,
        max_delivered_bytes=max_download_bytes,
    )
    return _HttpClient(
        NativeHttpClient(
            timeout=Timeout(connect=timeout, read=timeout),
            retry=RetryPolicy(retries=retries, respect_retry_after=respect_retry_after),
            redirects=RedirectPolicy(max_hops=5),
            budget=budget,
            cache=cache
            if cache is not None
            else (ResponseCache() if cache_ttl else None),
            ca_bundle=ca_bundle,
            proxy=proxy,
            trust_env=trust_env,
        ),
        cache_ttl=cache_ttl,
        cache_namespace=cache_namespace,
    )
