import json
import math
from datetime import datetime, timedelta, timezone

import pytest
import requests

from dtcc_core.io.data import orthophoto
from dtcc_core.io.data.orthophoto import (
    OrthophotoClientError,
    OrthophotoItem,
    fetch_manifest,
    file_url,
    resolve_server_url,
)

BASE = "http://tiles.test:8000"
BOUNDS = (673000.0, 6578000.0, 674000.0, 6581000.0)
TIMEOUT = (10.0, 150.0)


class FakeResponse:
    def __init__(self, status_code=200, body=b"", headers=None, exc=None):
        self.status_code = status_code
        self.body = body if isinstance(body, bytes) else body.encode()
        self.headers = {} if headers is None else headers
        self.exc = exc
        self.bytes_read = 0
        self.closed = False

    def iter_content(self, chunk_size=1):
        for start in range(0, len(self.body), chunk_size):
            chunk = self.body[start : start + chunk_size]
            self.bytes_read += len(chunk)
            yield chunk
        if self.exc is not None:
            raise self.exc

    @property
    def content(self):
        return b"".join(self.iter_content(65536))

    def close(self):
        self.closed = True


class FakeGet:
    def __init__(self, *responses):
        self.responses = list(responses)
        self.calls = []

    def __call__(self, url, **kwargs):
        self.calls.append({"url": url, **kwargs})
        response = self.responses.pop(0)
        if isinstance(response, BaseException):
            raise response
        return response


def item_json(**overrides):
    item = {
        "id": "o65775_6725_25_mr25",
        "collection": "orto-o2-2025",
        "bbox": [672500, 6577500, 675000, 6580000],
        "datetime": "2025-05-31T08:10:56Z",
        "spektraltyp": "rgbi",
        "resolution": 0.16,
        "size_bytes": 659037310,
        "href": "/files/orto-o2-2025/o65775_6725_25_mr25.tif",
    }
    for key, value in overrides.items():
        if value is _MISSING:
            item.pop(key, None)
        else:
            item[key] = value
    return item


_MISSING = object()


def manifest_response(items, crs="EPSG:3006"):
    payload = {"bbox": list(BOUNDS), "items": items}
    if crs is not _MISSING:
        payload["crs"] = crs
    return FakeResponse(200, json.dumps(payload))


def install(monkeypatch, *responses):
    fake = FakeGet(*responses)
    monkeypatch.setattr(orthophoto.requests, "get", fake)
    return fake


def fetch(**overrides):
    kwargs = {
        "server_url": BASE,
        "year": None,
        "collection": None,
        "spektraltyp": (),
        "timeout": TIMEOUT,
    }
    kwargs.update(overrides)
    bounds = kwargs.pop("bounds", BOUNDS)
    return fetch_manifest(bounds, **kwargs)


def fetch_error(monkeypatch, *responses):
    fake = install(monkeypatch, *responses)
    with pytest.raises(OrthophotoClientError) as info:
        fetch()
    assert len(fake.calls) == 1
    return info.value


# Server URL


def test_explicit_server_url_overrides_environment(monkeypatch):
    monkeypatch.setenv("DTCC_ORTHOPHOTO_URL", "http://env.test")
    assert resolve_server_url("http://explicit.test") == "http://explicit.test"


def test_environment_server_url_is_used_when_not_given(monkeypatch):
    monkeypatch.setenv("DTCC_ORTHOPHOTO_URL", "http://env.test/")
    assert resolve_server_url(None) == "http://env.test"


def test_missing_server_url_is_a_local_configuration_error(monkeypatch):
    monkeypatch.delenv("DTCC_ORTHOPHOTO_URL", raising=False)
    with pytest.raises(ValueError, match="DTCC_ORTHOPHOTO_URL"):
        resolve_server_url(None)


@pytest.mark.parametrize(
    "url, expected",
    [
        ("http://tiles.test:8000/", "http://tiles.test:8000"),
        ("https://host/proxy/lm", "https://host/proxy/lm"),
        ("https://host/proxy/lm/", "https://host/proxy/lm"),
    ],
)
def test_server_url_is_normalized(url, expected):
    assert resolve_server_url(url) == expected


@pytest.mark.parametrize(
    "url",
    [
        "",
        "ftp://tiles.test",
        "tiles.test:8000",
        "/relative",
        "http://",
        "http://tiles.test?key=1",
        "http://tiles.test/#frag",
    ],
)
def test_invalid_server_url_is_rejected(url):
    with pytest.raises(ValueError):
        resolve_server_url(url)


def test_prefixed_service_keeps_prefix_for_items_and_files(monkeypatch):
    fake = install(monkeypatch, manifest_response([item_json()]))
    (item,) = fetch(server_url="https://host/proxy/lm/")
    assert fake.calls[0]["url"] == "https://host/proxy/lm/items"
    assert (
        file_url("https://host/proxy/lm/", item)
        == "https://host/proxy/lm/files/orto-o2-2025/o65775_6725_25_mr25.tif"
    )


# Request


def test_manifest_request_sends_bbox_filters_and_timeout(monkeypatch):
    fake = install(monkeypatch, manifest_response([]))
    fetch(year=2025, collection="orto-o2-2025", spektraltyp=("rgb", "cir"))
    (call,) = fake.calls
    assert call["url"] == BASE + "/items"
    assert call["params"] == {
        "bbox": "673000.0,6578000.0,674000.0,6581000.0",
        "year": 2025,
        "collection": "orto-o2-2025",
        "spektraltyp": "rgb,cir",
    }
    assert call["timeout"] == TIMEOUT
    assert call["allow_redirects"] is False


def test_manifest_request_omits_unset_filters(monkeypatch):
    fake = install(monkeypatch, manifest_response([]))
    fetch()
    assert fake.calls[0]["params"] == {
        "bbox": "673000.0,6578000.0,674000.0,6581000.0"
    }


def test_each_call_fetches_a_fresh_manifest(monkeypatch):
    fake = install(
        monkeypatch, manifest_response([]), manifest_response([item_json()])
    )
    assert fetch() == []
    assert len(fetch()) == 1
    assert len(fake.calls) == 2


@pytest.mark.parametrize(
    "bounds",
    [
        (math.nan, 0.0, 1.0, 1.0),
        (0.0, 0.0, math.inf, 1.0),
        (1.0, 0.0, 1.0, 1.0),
        (0.0, 2.0, 1.0, 1.0),
        (0.0, 0.0, 10**400, 1.0),
    ],
)
def test_invalid_bounds_are_rejected_before_any_request(monkeypatch, bounds):
    fake = install(monkeypatch)
    with pytest.raises(ValueError):
        fetch(bounds=bounds)
    assert fake.calls == []


# Manifest parsing


def test_valid_item_is_parsed(monkeypatch):
    response = manifest_response([item_json()])
    install(monkeypatch, response)
    assert fetch() == [
        OrthophotoItem(
            id="o65775_6725_25_mr25",
            collection="orto-o2-2025",
            bbox=(672500.0, 6577500.0, 675000.0, 6580000.0),
            datetime=datetime(2025, 5, 31, 8, 10, 56, tzinfo=timezone.utc),
            spektraltyp="rgbi",
            href="/files/orto-o2-2025/o65775_6725_25_mr25.tif",
            resolution=0.16,
            size_bytes=659037310,
        )
    ]
    assert response.closed


def test_empty_manifest_returns_no_items(monkeypatch):
    install(monkeypatch, manifest_response([]))
    assert fetch() == []


def test_several_items_in_one_cell_are_distinct_and_keep_manifest_order(
    monkeypatch,
):
    newer = item_json()
    older = item_json(
        id="o65775_6725_25_im17",
        collection="orto-o2-2017",
        datetime="2017-06-01T00:00:00Z",
        spektraltyp="rgb",
        href="/files/orto-o2-2017/o65775_6725_25_im17.tif",
    )
    install(monkeypatch, manifest_response([newer, older]))
    items = fetch()
    assert [(i.collection, i.id) for i in items] == [
        ("orto-o2-2025", "o65775_6725_25_mr25"),
        ("orto-o2-2017", "o65775_6725_25_im17"),
    ]
    assert items[0].bbox == items[1].bbox


@pytest.mark.parametrize(
    "text, expected",
    [
        ("2025-05-31T08:10:56Z", datetime(2025, 5, 31, 8, 10, 56, tzinfo=timezone.utc)),
        (
            "2025-05-31T08:10:56.250Z",
            datetime(2025, 5, 31, 8, 10, 56, 250000, tzinfo=timezone.utc),
        ),
        (
            "2025-05-31T10:10:56+02:00",
            datetime(2025, 5, 31, 8, 10, 56, tzinfo=timezone.utc),
        ),
    ],
)
def test_timestamps_parse_to_comparable_aware_datetimes(monkeypatch, text, expected):
    install(monkeypatch, manifest_response([item_json(datetime=text)]))
    (item,) = fetch()
    assert item.datetime == expected
    assert item.datetime.tzinfo is not None


def test_offset_timestamps_compare_chronologically_not_as_text(monkeypatch):
    # 09:00+02:00 is 07:00Z: later as text, 30 minutes earlier in time.
    offset = item_json(datetime="2025-05-31T09:00:00+02:00")
    utc = item_json(
        id="o65775_6725_25_mr24",
        datetime="2025-05-31T07:30:00Z",
        href="/files/orto-o2-2025/o65775_6725_25_mr24.tif",
    )
    install(monkeypatch, manifest_response([offset, utc]))
    first, second = fetch()
    assert second.datetime - first.datetime == timedelta(minutes=30)


@pytest.mark.parametrize(
    "value, expected",
    [
        (0.16, 0.16),
        (1, 1.0),
        (None, None),
        (_MISSING, None),
        (True, None),
        ("0.16", None),
        (0, None),
        (-0.5, None),
        (math.nan, None),
        (math.inf, None),
        (10**400, None),
    ],
)
def test_optional_resolution_is_normalized_and_item_kept(monkeypatch, value, expected):
    install(monkeypatch, manifest_response([item_json(resolution=value)]))
    (item,) = fetch()
    assert item.resolution == expected


@pytest.mark.parametrize(
    "value, expected",
    [
        (659037310, 659037310),
        (None, None),
        (_MISSING, None),
        (True, None),
        ("659037310", None),
        (0, None),
        (-1, None),
        (1.5, None),
        (659037310.0, None),
        (10**400, 10**400),
    ],
)
def test_optional_size_bytes_is_normalized_and_item_kept(monkeypatch, value, expected):
    install(monkeypatch, manifest_response([item_json(size_bytes=value)]))
    (item,) = fetch()
    assert item.size_bytes == expected


@pytest.mark.parametrize("crs", ["EPSG:3021", None, _MISSING])
def test_manifest_crs_must_be_epsg_3006(monkeypatch, crs):
    error = fetch_error(monkeypatch, manifest_response([], crs=crs))
    assert error.failure_class == "invalid_payload"


@pytest.mark.parametrize(
    "payload",
    [
        [],
        {"crs": "EPSG:3006"},
        {"crs": "EPSG:3006", "items": {}},
        {"crs": "EPSG:3006", "items": ["not an item"]},
    ],
)
def test_manifest_structure_is_validated(monkeypatch, payload):
    error = fetch_error(monkeypatch, FakeResponse(200, json.dumps(payload)))
    assert error.failure_class == "invalid_payload"


@pytest.mark.parametrize(
    "overrides",
    [
        {"id": _MISSING},
        {"id": "O65775_UPPER", "href": "/files/orto-o2-2025/O65775_UPPER.tif"},
        {"id": "a.b", "href": "/files/orto-o2-2025/a.b.tif"},
        {"id": "a b", "href": "/files/orto-o2-2025/a b.tif"},
        {"id": "", "href": "/files/orto-o2-2025/.tif"},
        {"id": 7, "href": "/files/orto-o2-2025/7.tif"},
        {"collection": _MISSING},
        {
            "collection": "orto/o2",
            "href": "/files/orto/o2/o65775_6725_25_mr25.tif",
        },
        {
            "collection": "orto-o2-2025\n",
            "href": "/files/orto-o2-2025\n/o65775_6725_25_mr25.tif",
        },
        {"bbox": _MISSING},
        {"bbox": [672500, 6577500, 675000]},
        {"bbox": [672500, 6577500, 675000, "6580000"]},
        {"bbox": [True, 6577500, 675000, 6580000]},
        {"bbox": [672500, 6577500, 675000, math.nan]},
        {"bbox": [672500, 6577500, 675000, 10**400]},
        {"bbox": [675000, 6577500, 672500, 6580000]},
        {"bbox": [672500, 6580000, 675000, 6580000]},
        {"datetime": _MISSING},
        {"datetime": "not a date"},
        {"datetime": "2025-05-31T08:10:56"},
        {"datetime": 1748679056},
        {"spektraltyp": _MISSING},
        {"spektraltyp": ""},
        {"href": _MISSING},
        {"href": "/files/orto-o2-2017/o65775_6725_25_mr25.tif"},
        {"href": "/files/orto-o2-2025/other.tif"},
        {"href": "http://evil.test/files/orto-o2-2025/o65775_6725_25_mr25.tif"},
        {"href": "//evil.test/files/orto-o2-2025/o65775_6725_25_mr25.tif"},
        {"href": "/files/orto-o2-2025/../o65775_6725_25_mr25.tif"},
        {"href": "/files/orto-o2-2025%2Fo65775_6725_25_mr25.tif"},
        {"href": "/files/orto-o2-2025/o65775_6725_25_mr25.tif?x=1"},
        {"href": "files/orto-o2-2025/o65775_6725_25_mr25.tif"},
    ],
)
def test_invalid_required_fields_are_rejected(monkeypatch, overrides):
    error = fetch_error(monkeypatch, manifest_response([item_json(**overrides)]))
    assert error.failure_class == "invalid_payload"
    assert error.operation == "manifest"
    assert error.target == BASE + "/items"


# Failures


@pytest.mark.parametrize(
    "exc, failure_class",
    [
        (requests.ConnectTimeout("connect timed out"), "timeout"),
        (requests.ReadTimeout("read timed out"), "timeout"),
        (requests.ConnectionError("refused"), "connection"),
    ],
)
def test_transport_failures_are_classified(monkeypatch, exc, failure_class):
    error = fetch_error(monkeypatch, exc)
    assert error.failure_class == failure_class
    assert error.status_code is None


def test_body_interrupted_while_reading_is_a_connection_failure(monkeypatch):
    response = FakeResponse(
        200, b'{"crs": "EPSG:3006", "it', exc=requests.exceptions.ChunkedEncodingError()
    )
    error = fetch_error(monkeypatch, response)
    assert error.failure_class == "connection"
    assert error.status_code == 200
    assert response.closed


@pytest.mark.parametrize(
    "body",
    [b"{not json", b"<html><body>Bad gateway</body></html>", b""],
)
def test_malformed_success_body_is_invalid_payload(monkeypatch, body):
    error = fetch_error(monkeypatch, FakeResponse(200, body))
    assert error.failure_class == "invalid_payload"
    assert error.status_code == 200


@pytest.mark.parametrize(
    "status, body, failure_class",
    [
        (503, {"detail": "LM credentials not configured"}, "configuration"),
        (503, {"detail": "LM credentials not configured."}, "http_5xx"),
        (503, {"detail": "Service unavailable"}, "http_5xx"),
        (502, "<html>Bad gateway</html>", "http_5xx"),
        (504, {"detail": "LM timed out"}, "http_5xx"),
        (500, "", "http_5xx"),
        (400, {"detail": "bad bbox"}, "http_4xx"),
        (404, {"detail": "not found"}, "http_4xx"),
        (422, {"detail": []}, "http_4xx"),
        (301, "", "invalid_payload"),
        (302, "", "invalid_payload"),
        (307, "", "invalid_payload"),
        (201, {"crs": "EPSG:3006", "items": []}, "invalid_payload"),
    ],
)
def test_status_codes_are_classified_and_preserved(
    monkeypatch, status, body, failure_class
):
    if not isinstance(body, str):
        body = json.dumps(body)
    headers = {"Location": "http://elsewhere.test/"}
    response = FakeResponse(status, body, headers=headers)
    error = fetch_error(monkeypatch, response)
    assert error.failure_class == failure_class
    assert error.status_code == status
    assert error.operation == "manifest"
    assert response.closed
    if 300 <= status < 400:
        assert "redirect" in str(error).lower()


def test_error_detail_is_bounded_without_reading_whole_body(monkeypatch):
    response = FakeResponse(502, b"x" * 1_000_000)
    error = fetch_error(monkeypatch, response)
    assert error.detail == "x" * 1024
    assert len(str(error)) <= 1024
    assert response.bytes_read <= 1024 * 2


def test_detail_of_credentials_error_is_kept(monkeypatch):
    body = json.dumps({"detail": "LM credentials not configured"})
    error = fetch_error(monkeypatch, FakeResponse(503, body))
    assert "LM credentials not configured" in error.detail


def test_corrupt_encoded_body_is_invalid_payload(monkeypatch):
    response = FakeResponse(
        200, b"", exc=requests.exceptions.ContentDecodingError("bad gzip")
    )
    error = fetch_error(monkeypatch, response)
    assert error.failure_class == "invalid_payload"
    assert error.status_code == 200
    assert response.closed


@pytest.mark.parametrize(
    "response",
    [
        manifest_response([], crs="E" * 5000),
        manifest_response([item_json(id="A" * 5000)]),
        manifest_response([item_json(href="/" + "h" * 5000)]),
        FakeResponse(200, b"{" + b"x" * 5000),
        requests.ConnectionError("c" * 5000),
        requests.ReadTimeout("t" * 5000),
    ],
)
def test_error_messages_and_details_are_bounded(monkeypatch, response):
    error = fetch_error(monkeypatch, response)
    assert len(str(error)) <= 1024
    assert len(error.message) <= 1024
    assert error.detail is None or len(error.detail) <= 1024


def test_client_error_truncates_message_and_detail():
    error = OrthophotoClientError(
        "m" * 5000,
        operation="manifest",
        target=BASE + "/items",
        failure_class="http_5xx",
        detail="d" * 5000,
    )
    assert error.message == str(error) == "m" * 1024
    assert error.detail == "d" * 1024
