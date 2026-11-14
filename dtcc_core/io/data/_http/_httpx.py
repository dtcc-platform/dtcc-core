"""Native HTTPX2 buffered GETs using the shared HTTP policies and accounting."""

import os
import ssl
import time
from urllib.parse import urlencode, urlsplit, urlunsplit

import httpx2
from certifi import where

from ._contracts import HttpLimitError, HttpTransportError, Response
from ._policy import Operation, State, headers_lower, parse_retry_after


def _certificate_error(exc):
    pending, seen = [exc], set()
    while pending:
        current = pending.pop()
        if not isinstance(current, BaseException) or id(current) in seen:
            continue
        if isinstance(current, ssl.SSLCertVerificationError):
            return True
        seen.add(id(current))
        pending.extend((current.__cause__, current.__context__, *current.args))
    return False


def _kind(exc):
    if _certificate_error(exc):
        return "certificate"
    if isinstance(exc, httpx2.TimeoutException):
        return "timeout"
    if isinstance(exc, httpx2.DecodingError):
        return "decoding"
    if isinstance(exc, (httpx2.NetworkError, httpx2.RemoteProtocolError)):
        return "connection"
    return "request"


class NativeHttpClient:
    """Serial, owned HTTPX2 connection pool; no backend retries or cookie reuse."""

    def __init__(self, *, ca_bundle=None, proxy=None, trust_env=False, **policies):
        self._state = State(**policies)
        self._ca_bundle, self._proxy, self._trust_env = ca_bundle, proxy, trust_env
        self._client = None
        self._closed = False
        self._sleep = time.sleep

    @property
    def counters(self):
        return self._state.counters

    def _build_client(self):
        if self._client is None:
            bundle = (
                self._ca_bundle
                or os.environ.get("REQUESTS_CA_BUNDLE")
                or os.environ.get("CURL_CA_BUNDLE")
                or where()
            )
            trust = ssl.create_default_context(
                **({"capath": bundle} if os.path.isdir(bundle) else {"cafile": bundle})
            )
            timeout = self._state.timeout
            self._client = httpx2.Client(
                timeout=httpx2.Timeout(timeout.read, connect=timeout.connect),
                trust_env=self._trust_env,
                follow_redirects=False,
                verify=trust,
                proxy=self._proxy,
            )
        return self._client

    def close(self):
        if not self._closed and self._client is not None:
            self._client.close()
        self._closed = True

    def _chunks(self, response, operation, *, limited=True):
        iterator = response.iter_bytes()
        while True:
            operation.before_read()
            chunk = next(iterator, None)
            if chunk is None:
                return
            operation.received(chunk, limited=limited)
            yield chunk

    def _status_error(self, response, operation):
        body = bytearray()
        try:
            for chunk in self._chunks(response, operation, limited=False):
                body.extend(chunk[: 8192 - len(body)])
                if len(body) == 8192:
                    break
        except HttpLimitError:
            raise
        except (httpx2.RequestError, HttpTransportError):
            # Keep the status when its diagnostic body is truncated or broken.
            pass
        return operation.error(
            "http_status", status=response.status_code, body=bytes(body), phase="read"
        )

    def _attempt(self, operation):
        request = operation.request
        client = self._build_client()
        client.cookies.clear()
        response, headers = None, {}
        try:
            operation.attempt()
            url = request.url
            if request.params:
                parsed = urlsplit(url)
                query = "&".join(
                    filter(None, (parsed.query, urlencode(request.params)))
                )
                url = urlunsplit(parsed._replace(query=query))
            with client.stream("GET", url, headers=request.headers) as response:
                headers, status = headers_lower(response.headers), response.status_code
                if status in {301, 302, 303, 307, 308}:
                    location = headers.get("location")
                    if not location:
                        raise operation.error(
                            "redirect_policy", status=status, phase="read"
                        )
                    operation.follow(location, status=status)
                    return None, headers, None
                if status not in request.expected_status:
                    raise self._status_error(response, operation)
                body = b"".join(self._chunks(response, operation))
                return (
                    Response(status, request.target, headers, body, request),
                    headers,
                    None,
                )
        except httpx2.RequestError as exc:
            phase = (
                "read"
                if response is not None or isinstance(exc, httpx2.ReadTimeout)
                else "connect"
            )
            error = operation.error(
                _kind(exc),
                status=response.status_code if response is not None else None,
                phase=phase,
            )
        except (httpx2.InvalidURL, UnicodeError):
            error = operation.error("request")
        except HttpTransportError as exc:
            error = exc
        finally:
            client.cookies.clear()
        error.retry_after = parse_retry_after(
            headers.get("retry-after"), self._state.wall_clock
        )
        return None, headers, error

    def request(self, request):
        if self._closed:
            raise ValueError("Client is closed")
        operation = Operation(self._state, request, None)
        cached = self._state.cached(operation)
        if cached is not None:
            return cached
        while True:
            response, headers, error = self._attempt(operation)
            if error is not None:
                if not operation.retry(error, error.phase):
                    raise error from None
                self._sleep(operation.delay(headers))
                continue
            if response is None:
                continue
            operation.delivered(len(response.body))
            if not operation.hops:
                self._state.store(request, response)
            return response
