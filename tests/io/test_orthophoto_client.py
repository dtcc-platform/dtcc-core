import hashlib
import json
import math
import os
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np
import pytest
import rasterio
import requests
from affine import Affine
from rasterio.transform import from_bounds
from requests.structures import CaseInsensitiveDict

from dtcc_core.model.values.raster_tiles import RasterTileCollection

from dtcc_core.io.data import orthophoto
from dtcc_core.io.data.orthophoto import (
    OrthophotoClientError,
    OrthophotoItem,
    acquire_tiles,
    download_tile,
    fetch_manifest,
    file_url,
    resolve_server_url,
)

BASE = "http://tiles.test:8000"
BOUNDS = (673000.0, 6578000.0, 674000.0, 6581000.0)
TIMEOUT = (10.0, 150.0)


class FakeResponse:
    def __init__(
        self, status_code=200, body=b"", headers=None, exc=None, during=None
    ):
        self.status_code = status_code
        self.body = body if isinstance(body, bytes) else body.encode()
        self.headers = CaseInsensitiveDict(headers or {})
        self.exc = exc
        self.during = during
        self.bytes_read = 0
        self.closed = False

    def iter_content(self, chunk_size=1):
        for start in range(0, len(self.body), chunk_size):
            chunk = self.body[start : start + chunk_size]
            self.bytes_read += len(chunk)
            yield chunk
            if self.during is not None:
                during, self.during = self.during, None
                during()
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


# Original-file cache

CELL = (672500.0, 6577500.0, 672502.0, 6577502.0)
ITEM = OrthophotoItem(
    id="o65775_6725_25_mr25",
    collection="orto-o2-2025",
    bbox=CELL,
    datetime=datetime(2025, 5, 31, 8, 10, 56, tzinfo=timezone.utc),
    spektraltyp="rgbi",
    href="/files/orto-o2-2025/o65775_6725_25_mr25.tif",
    resolution=0.5,
    size_bytes=None,
)


def tiff_bytes(tmp_path, bounds=CELL, crs="EPSG:3006", transform=None):
    path = tmp_path / "src" / "tile.tif"
    path.parent.mkdir(exist_ok=True)
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        width=4,
        height=4,
        count=3,
        dtype="uint8",
        crs=crs,
        transform=transform or from_bounds(*bounds, 4, 4),
    ) as dst:
        dst.write(np.arange(48, dtype=np.uint8).reshape(3, 4, 4))
    data = path.read_bytes()
    path.unlink()
    return data


def tile_response(body, length=True, **kwargs):
    headers = kwargs.pop("headers", {})
    if length is True:
        headers["Content-Length"] = str(len(body))
    elif length is not None:
        headers["Content-Length"] = length
    return FakeResponse(200, body, headers=headers, **kwargs)


def download(cache_root, server_url=BASE, item=ITEM):
    return download_tile(
        item, server_url=server_url, timeout=TIMEOUT, cache_root=cache_root
    )


def files_under(root):
    return sorted(
        str(p.relative_to(root)) for p in Path(root).rglob("*") if p.is_file()
    )


def download_error(monkeypatch, cache_root, *responses):
    fake = install(monkeypatch, *responses)
    with pytest.raises(OrthophotoClientError) as info:
        download(cache_root)
    assert len(fake.calls) == 1
    assert not any(name.endswith(".tif") for name in files_under(cache_root))
    assert [n for n in files_under(cache_root) if not n.endswith(".tif")] == []
    assert info.value.operation == "download"
    return info.value


def test_download_publishes_original_bytes_under_service_namespace(
    monkeypatch, tmp_path
):
    body = tiff_bytes(tmp_path)
    cache_root = tmp_path / "cache"
    install(monkeypatch, tile_response(body))
    path = download(cache_root)
    service = path.parent.parent.name
    assert path == (
        cache_root / "orthophoto" / service / "orto-o2-2025" / "o65775_6725_25_mr25.tif"
    )
    assert len(service) == 16 and set(service) <= set("0123456789abcdef")
    assert path.read_bytes() == body
    assert files_under(tmp_path) == [str(path.relative_to(tmp_path))]


def test_download_request_uses_identity_encoding_and_no_redirects(
    monkeypatch, tmp_path
):
    fake = install(monkeypatch, tile_response(tiff_bytes(tmp_path)))
    download(tmp_path / "cache", server_url="https://host/proxy/lm/")
    (call,) = fake.calls
    assert call["url"] == (
        "https://host/proxy/lm/files/orto-o2-2025/o65775_6725_25_mr25.tif"
    )
    assert call["headers"]["Accept-Encoding"] == "identity"
    assert call["allow_redirects"] is False
    assert call["stream"] is True
    assert call["timeout"] == TIMEOUT


def test_service_namespace_includes_prefix_and_ignores_trailing_slash(
    monkeypatch, tmp_path
):
    body = tiff_bytes(tmp_path)
    cache_root = tmp_path / "cache"
    urls = ["http://h.test/a", "http://h.test/a/", "http://h.test/b", "http://g.test/a"]
    install(monkeypatch, *(tile_response(body) for _ in urls))
    paths = []
    for url in urls:
        paths.append(download(cache_root, server_url=url))
        if url == "http://h.test/a":
            paths[-1].unlink()
    assert paths[0] == paths[1]
    assert len({paths[1], paths[2], paths[3]}) == 3


def test_published_file_is_reused_without_request_or_header_read(
    monkeypatch, tmp_path
):
    body = tiff_bytes(tmp_path)
    cache_root = tmp_path / "cache"
    install(monkeypatch, tile_response(body))
    path = download(cache_root)
    path.write_bytes(b"cached bytes are trusted")

    def forbidden(*args, **kwargs):
        raise AssertionError("no request or header read on a cache hit")

    monkeypatch.setattr(orthophoto.requests, "get", forbidden)
    monkeypatch.setattr(orthophoto.rasterio, "open", forbidden)
    assert download(cache_root) == path
    assert path.read_bytes() == b"cached bytes are trusted"


def test_missing_file_is_downloaded_again(monkeypatch, tmp_path):
    body = tiff_bytes(tmp_path)
    cache_root = tmp_path / "cache"
    fake = install(monkeypatch, tile_response(body), tile_response(body))
    download(cache_root).unlink()
    assert download(cache_root).read_bytes() == body
    assert len(fake.calls) == 2


@pytest.mark.parametrize(
    "status, body, headers, failure_class",
    [
        (
            503,
            {"detail": "LM credentials not configured"},
            {"Content-Encoding": "gzip"},
            "configuration",
        ),
        (404, {"detail": "not found"}, {}, "http_4xx"),
        (502, {"detail": "LM error"}, {}, "http_5xx"),
        (504, {"detail": "LM timed out"}, {}, "http_5xx"),
        (302, "", {"Location": "http://elsewhere.test/x.tif"}, "invalid_payload"),
    ],
)
def test_download_status_is_classified_before_body_checks(
    monkeypatch, tmp_path, status, body, headers, failure_class
):
    if not isinstance(body, str):
        body = json.dumps(body)
    response = FakeResponse(status, body, headers=headers)
    error = download_error(monkeypatch, tmp_path / "cache", response)
    assert error.failure_class == failure_class
    assert error.status_code == status
    assert response.closed


@pytest.mark.parametrize(
    "exc, failure_class",
    [
        (requests.ConnectTimeout("connect"), "timeout"),
        (requests.ReadTimeout("read"), "timeout"),
        (requests.ConnectionError("refused"), "connection"),
    ],
)
def test_download_transport_failures_are_classified(
    monkeypatch, tmp_path, exc, failure_class
):
    error = download_error(monkeypatch, tmp_path / "cache", exc)
    assert error.failure_class == failure_class


def test_encoded_response_is_invalid_payload(monkeypatch, tmp_path):
    response = tile_response(
        tiff_bytes(tmp_path), headers={"Content-Encoding": "gzip"}
    )
    error = download_error(monkeypatch, tmp_path / "cache", response)
    assert error.failure_class == "invalid_payload"
    assert error.status_code == 200
    assert response.closed


def test_body_shorter_than_content_length_is_connection_failure(
    monkeypatch, tmp_path
):
    body = tiff_bytes(tmp_path)
    response = tile_response(body, length=str(len(body) + 100))
    error = download_error(monkeypatch, tmp_path / "cache", response)
    assert error.failure_class == "connection"
    assert error.status_code == 200


def test_body_longer_than_content_length_is_invalid_payload(monkeypatch, tmp_path):
    body = tiff_bytes(tmp_path)
    response = tile_response(body, length=str(len(body) - 1))
    error = download_error(monkeypatch, tmp_path / "cache", response)
    assert error.failure_class == "invalid_payload"


@pytest.mark.parametrize(
    "length",
    ["abc", "-5", "", "1e3", "10, 10", " 5000", "\u0665\u0660", "9" * 5000],
)
def test_malformed_content_length_is_rejected_before_reading_body(
    monkeypatch, tmp_path, length
):
    response = tile_response(tiff_bytes(tmp_path), length=length)
    error = download_error(monkeypatch, tmp_path / "cache", response)
    assert error.failure_class == "invalid_payload"
    assert error.status_code == 200
    assert response.bytes_read == 0
    assert response.closed


def test_body_interrupted_mid_stream_is_connection_failure(monkeypatch, tmp_path):
    body = tiff_bytes(tmp_path)
    response = tile_response(
        body[:200],
        length=str(len(body)),
        exc=requests.exceptions.ChunkedEncodingError("connection broken"),
    )
    error = download_error(monkeypatch, tmp_path / "cache", response)
    assert error.failure_class == "connection"
    assert error.status_code == 200
    assert response.closed


@pytest.mark.parametrize(
    "body",
    [b"<html><body>Bad gateway</body></html>", b"", "truncated"],
)
def test_unreadable_tiff_is_invalid_payload(monkeypatch, tmp_path, body):
    if body == "truncated":
        body = tiff_bytes(tmp_path)[:16]
    error = download_error(
        monkeypatch, tmp_path / "cache", tile_response(body, length=None)
    )
    assert error.failure_class == "invalid_payload"


@pytest.mark.parametrize("crs", ["EPSG:3021", None])
def test_tiff_in_other_crs_is_invalid_payload(monkeypatch, tmp_path, crs):
    response = tile_response(tiff_bytes(tmp_path, crs=crs))
    error = download_error(monkeypatch, tmp_path / "cache", response)
    assert error.failure_class == "invalid_payload"


@pytest.mark.parametrize(
    "transform",
    [
        Affine(math.nan, 0.0, math.nan, 0.0, math.nan, math.nan),
        Affine(0.5, 0.0, math.nan, 0.0, -0.5, CELL[3]),
    ],
)
def test_tiff_with_nan_georeferencing_is_invalid_payload(
    monkeypatch, tmp_path, transform
):
    response = tile_response(tiff_bytes(tmp_path, transform=transform))
    error = download_error(monkeypatch, tmp_path / "cache", response)
    assert error.failure_class == "invalid_payload"
    assert error.status_code == 200
    assert response.closed


def test_vrt_payload_is_never_opened_as_vrt(monkeypatch, tmp_path):
    source = tmp_path / "src" / "valid.tif"
    source.parent.mkdir()
    source.write_bytes(tiff_bytes(tmp_path))
    body = (
        '<VRTDataset rasterXSize="4" rasterYSize="4"><SRS>EPSG:3006</SRS>'
        f"<GeoTransform>{CELL[0]},0.5,0,{CELL[3]},0,-0.5</GeoTransform>"
        '<VRTRasterBand dataType="Byte" band="1"><SimpleSource>'
        f'<SourceFilename relativeToVRT="0">{source}</SourceFilename>'
        "<SourceBand>1</SourceBand></SimpleSource></VRTRasterBand></VRTDataset>"
    ).encode()
    with rasterio.open(source) as valid:
        assert valid.driver == "GTiff"
    opened = []
    real_open = rasterio.open

    def recording_open(*args, **kwargs):
        dataset = real_open(*args, **kwargs)
        opened.append(dataset.driver)
        return dataset

    monkeypatch.setattr(orthophoto.rasterio, "open", recording_open)
    error = download_error(monkeypatch, tmp_path / "cache", tile_response(body))
    assert error.failure_class == "invalid_payload"
    assert error.status_code == 200
    assert "VRT" not in opened


def test_tiff_for_wrong_cell_is_invalid_payload(monkeypatch, tmp_path):
    shifted = tuple(v + 2500.0 if i % 2 == 0 else v for i, v in enumerate(CELL))
    response = tile_response(tiff_bytes(tmp_path, bounds=shifted))
    error = download_error(monkeypatch, tmp_path / "cache", response)
    assert error.failure_class == "invalid_payload"


@pytest.mark.parametrize("edge", range(4))
@pytest.mark.parametrize("offset, accepted", [(0.0005, True), (0.002, False)])
def test_bounds_tolerance_is_absolute_per_edge(
    monkeypatch, tmp_path, edge, offset, accepted
):
    bounds = list(CELL)
    bounds[edge] += offset if edge >= 2 else -offset
    response = tile_response(tiff_bytes(tmp_path, bounds=tuple(bounds)))
    cache_root = tmp_path / "cache"
    if accepted:
        install(monkeypatch, response)
        assert download(cache_root).is_file()
    else:
        error = download_error(monkeypatch, cache_root, response)
        assert error.failure_class == "invalid_payload"


def test_failed_write_propagates_and_removes_temporary_file(monkeypatch, tmp_path):
    install(monkeypatch, tile_response(tiff_bytes(tmp_path)))

    class FullDisk:
        def __init__(self, fd, mode):
            os.close(fd)

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

        def write(self, data):
            raise OSError(28, "No space left on device")

    monkeypatch.setattr(orthophoto, "open", FullDisk, raising=False)
    with pytest.raises(OSError) as info:
        download(tmp_path / "cache")
    assert not isinstance(info.value, OrthophotoClientError)
    assert info.value.errno == 28
    assert files_under(tmp_path / "cache") == []


def test_failed_replace_propagates_and_removes_temporary_file(monkeypatch, tmp_path):
    install(monkeypatch, tile_response(tiff_bytes(tmp_path)))

    def deny(src, dst):
        raise PermissionError(13, "Permission denied")

    monkeypatch.setattr(orthophoto.os, "replace", deny)
    with pytest.raises(PermissionError):
        download(tmp_path / "cache")
    assert files_under(tmp_path / "cache") == []


@pytest.mark.parametrize("exc", [RuntimeError("boom"), KeyboardInterrupt()])
def test_unexpected_errors_propagate_unchanged_and_clean_up(
    monkeypatch, tmp_path, exc
):
    body = tiff_bytes(tmp_path)
    response = tile_response(body[:200], length=str(len(body)), exc=exc)
    install(monkeypatch, response)
    with pytest.raises(type(exc)):
        download(tmp_path / "cache")
    assert files_under(tmp_path / "cache") == []
    assert response.closed


def test_header_check_error_propagates_unchanged(monkeypatch, tmp_path):
    install(monkeypatch, tile_response(tiff_bytes(tmp_path)))

    def broken_open(*args, **kwargs):
        raise RuntimeError("unexpected")

    monkeypatch.setattr(orthophoto.rasterio, "open", broken_open)
    with pytest.raises(RuntimeError, match="unexpected"):
        download(tmp_path / "cache")
    assert files_under(tmp_path / "cache") == []


def test_concurrent_downloads_use_distinct_temporary_files(monkeypatch, tmp_path):
    body = tiff_bytes(tmp_path)
    cache_root = tmp_path / "cache"
    temporary_names = []
    real_mkstemp = orthophoto.tempfile.mkstemp

    def recording_mkstemp(**kwargs):
        fd, name = real_mkstemp(**kwargs)
        temporary_names.append(name)
        return fd, name

    monkeypatch.setattr(orthophoto.tempfile, "mkstemp", recording_mkstemp)
    inner = tile_response(body)
    outer = tile_response(body, during=lambda: download(cache_root))
    install(monkeypatch, outer, inner)
    path = download(cache_root)
    assert path.read_bytes() == body
    assert len(temporary_names) == 2
    assert len(set(temporary_names)) == 2
    assert Path(temporary_names[0]).parent == path.parent
    assert files_under(cache_root) == [str(path.relative_to(cache_root))]


def test_failed_download_leaves_competing_writer_files_intact(monkeypatch, tmp_path):
    body = tiff_bytes(tmp_path)
    cache_root = tmp_path / "cache"
    competitor = {}

    def competing_writer():
        directory = next((cache_root / "orthophoto").iterdir()) / "orto-o2-2025"
        competitor["published"] = directory / "o65775_6725_25_mr25.tif"
        competitor["temporary"] = directory / ".o65775_6725_25_mr25.other.part"
        competitor["published"].write_bytes(b"competitor original")
        competitor["temporary"].write_bytes(b"competitor in progress")

    response = tile_response(
        body[:200], length=str(len(body)), during=competing_writer
    )
    install(monkeypatch, response)
    with pytest.raises(OrthophotoClientError) as info:
        download(cache_root)
    assert info.value.failure_class == "connection"
    assert competitor["published"].read_bytes() == b"competitor original"
    assert competitor["temporary"].read_bytes() == b"competitor in progress"
    assert len(files_under(cache_root)) == 2


# Tile acquisition


def cell_item_json(item_id, collection, timestamp, **overrides):
    values = {
        "id": item_id,
        "collection": collection,
        "bbox": list(CELL),
        "datetime": timestamp,
        "href": f"/files/{collection}/{item_id}.tif",
        "resolution": 0.5,
        "size_bytes": None,
    }
    values.update(overrides)
    return item_json(**values)


NEWER = cell_item_json("o65775_6725_25_mr25", "orto-o2-2025", "2025-05-31T08:10:56Z")
OLDER = cell_item_json(
    "o65775_6725_25_im17",
    "orto-o2-2017",
    "2017-06-01T10:00:00+02:00",
    spektraltyp="rgb",
    resolution=None,
    size_bytes=1234,
)


def acquire(cache_root, **overrides):
    kwargs = {
        "server_url": BASE,
        "year": None,
        "collection": None,
        "spektraltyp": (),
        "timeout": TIMEOUT,
        "cache_root": cache_root,
    }
    kwargs.update(overrides)
    return acquire_tiles(BOUNDS, **kwargs)


def test_acquire_tiles_returns_records_in_manifest_order(monkeypatch, tmp_path):
    body = tiff_bytes(tmp_path)
    cache_root = tmp_path / "cache"
    fake = install(
        monkeypatch,
        manifest_response([NEWER, OLDER]),
        tile_response(body),
        tile_response(body),
    )
    tiles = acquire(
        cache_root, year=2025, collection="orto-o2-2025", spektraltyp=("rgb", "rgbi")
    )
    assert isinstance(tiles, RasterTileCollection)
    assert [call["url"] for call in fake.calls] == [
        BASE + "/items",
        BASE + "/files/orto-o2-2025/o65775_6725_25_mr25.tif",
        BASE + "/files/orto-o2-2017/o65775_6725_25_im17.tif",
    ]
    assert fake.calls[0]["params"]["year"] == 2025
    assert fake.calls[0]["params"]["collection"] == "orto-o2-2025"
    assert fake.calls[0]["params"]["spektraltyp"] == "rgb,rgbi"
    assert all(call["timeout"] == TIMEOUT for call in fake.calls)

    newer, older = tiles
    service_dir = newer.path.parent.parent
    assert service_dir.parent == cache_root / "orthophoto"
    assert newer.path == service_dir / "orto-o2-2025" / "o65775_6725_25_mr25.tif"
    assert older.path == service_dir / "orto-o2-2017" / "o65775_6725_25_im17.tif"
    assert newer.path.read_bytes() == older.path.read_bytes() == body
    assert (newer.id, newer.collection, newer.spektraltyp) == (
        "o65775_6725_25_mr25",
        "orto-o2-2025",
        "rgbi",
    )
    assert newer.datetime == datetime(2025, 5, 31, 8, 10, 56, tzinfo=timezone.utc)
    assert older.datetime == datetime(2017, 6, 1, 8, 0, tzinfo=timezone.utc)
    assert newer.extent.tuple == older.extent.tuple == CELL
    assert newer.crs == older.crs == "EPSG:3006"
    assert (newer.resolution, newer.size_bytes) == (0.5, None)
    assert (older.resolution, older.size_bytes) == (None, 1234)
    assert (older.id, older.collection, older.spektraltyp) == (
        "o65775_6725_25_im17",
        "orto-o2-2017",
        "rgb",
    )


def test_repeated_acquisition_refreshes_manifest_and_reuses_files(
    monkeypatch, tmp_path
):
    body = tiff_bytes(tmp_path)
    cache_root = tmp_path / "cache"
    install(monkeypatch, manifest_response([NEWER]), tile_response(body))
    (first,) = acquire(cache_root)
    fake = install(monkeypatch, manifest_response([NEWER, OLDER]), tile_response(body))
    second = acquire(cache_root)
    assert [call["url"] for call in fake.calls] == [
        BASE + "/items",
        BASE + "/files/orto-o2-2017/o65775_6725_25_im17.tif",
    ]
    assert second[0].path == first.path
    assert len(second) == 2


def test_empty_manifest_gives_empty_collection(monkeypatch, tmp_path):
    fake = install(monkeypatch, manifest_response([]))
    tiles = acquire(tmp_path / "cache")
    assert isinstance(tiles, RasterTileCollection)
    assert len(tiles) == 0
    assert len(fake.calls) == 1
    assert files_under(tmp_path) == []


def test_first_download_error_stops_acquisition(monkeypatch, tmp_path):
    body = tiff_bytes(tmp_path)
    cache_root = tmp_path / "cache"
    third = cell_item_json(
        "o65775_6725_25_mr24", "orto-o2-2024", "2024-06-01T00:00:00Z"
    )
    fake = install(
        monkeypatch,
        manifest_response([NEWER, OLDER, third]),
        tile_response(body),
        FakeResponse(502, json.dumps({"detail": "LM error"})),
    )
    with pytest.raises(OrthophotoClientError) as info:
        acquire(cache_root)
    assert info.value.failure_class == "http_5xx"
    assert info.value.target == BASE + "/files/orto-o2-2017/o65775_6725_25_im17.tif"
    assert len(fake.calls) == 3
    (published,) = files_under(cache_root)
    assert published.endswith("orto-o2-2025/o65775_6725_25_mr25.tif")
    assert (cache_root / published).read_bytes() == body


def test_manifest_error_propagates_before_downloads(monkeypatch, tmp_path):
    fake = install(monkeypatch, requests.ConnectionError("refused"))
    with pytest.raises(OrthophotoClientError) as info:
        acquire(tmp_path / "cache")
    assert info.value.operation == "manifest"
    assert len(fake.calls) == 1


def test_local_cache_error_propagates_unchanged(monkeypatch, tmp_path):
    cache_root = tmp_path / "cache"
    cache_root.write_text("not a directory")
    body = tiff_bytes(tmp_path)
    install(monkeypatch, manifest_response([NEWER]), tile_response(body))
    with pytest.raises(OSError) as info:
        acquire(cache_root)
    assert not isinstance(info.value, OrthophotoClientError)
