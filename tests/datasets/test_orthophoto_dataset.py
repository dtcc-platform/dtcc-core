import gc
import importlib
import json
import math
import dataclasses
import struct
from functools import partial
import traceback
import warnings
import weakref
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pytest
import rasterio
import requests
from affine import Affine
from pydantic import ValidationError
from rasterio.enums import ColorInterp
from rasterio.io import MemoryFile
from requests.structures import CaseInsensitiveDict

import dtcc_core.datasets as datasets
from dtcc_core.datasets import attach_dataset_context
import dtcc_core.io as dtcc_io
from dtcc_core.common import progress as progress_module
from dtcc_core.datasets.dataset import DatasetDescriptor, DatasetUpstreamError
from dtcc_core.datasets.package import export_model_package, load_model_package
from dtcc_core.datasets.publish import _validate_upload_package
from dtcc_core.datasets.schema import DatasetContext
from dtcc_core.datasets.orthophoto import OrthophotoArgs, OrthophotoDataset
from dtcc_core.io import raster as raster_io
from dtcc_core.io import raster_tiles
from dtcc_core.io.data import cache as data_cache
from dtcc_core.io.data import orthophoto as client
from dtcc_core.io.raster_tiles import READ_OVERHEAD_BYTES, WorkingMemoryError
from dtcc_core.model import Bounds, Raster, exchange
from dtcc_core.model.values.raster_tiles import RasterTileCollection
import tests.datasets.live.test_orthophoto_live as live
from tests.datasets.test_dataset_publish import RecordingUploader

# The package attribute is the registered dataset; the module is imported directly.
orthophoto_module = importlib.import_module("dtcc_core.datasets.orthophoto")

BASE = "http://tiles.test:8000"
X0, Y0 = 672500.0, 6577500.0


def cell(i, j=0, size=2.0):
    return (X0 + i * size, Y0 + j * size, X0 + (i + 1) * size, Y0 + (j + 1) * size)


# Fake tile server


class FakeResponse:
    def __init__(self, status_code=200, body=b"", headers=None, chunk=None):
        self.status_code = status_code
        self.body = body if isinstance(body, bytes) else body.encode()
        self.headers = CaseInsensitiveDict(headers or {})
        self.chunk = chunk

    def iter_content(self, chunk_size=1):
        # A fixed chunk, when set, replaces the size the client asks for.
        chunk_size = self.chunk or chunk_size
        for start in range(0, len(self.body), chunk_size):
            yield self.body[start : start + chunk_size]

    @property
    def content(self):
        return self.body

    def close(self):
        pass


class Server:
    """Answers /items with the listed items and /files/... with their bytes."""

    def __init__(self):
        self.items = []
        self.files = {}
        self.manifest_failure = None
        self.file_failures = {}
        self.calls = []
        self.before = None
        self.chunk = None
        self.without_length = set()

    def add(self, name, bounds, body, *, year=2020, spektraltyp="rgbi",
            resolution=0.5, collection=None):
        collection = collection or f"orto-{year}"
        href = f"/files/{collection}/{name}.tif"
        self.items.append(
            {
                "id": name,
                "collection": collection,
                "bbox": list(bounds),
                "datetime": f"{year}-06-01T10:00:00Z",
                "spektraltyp": spektraltyp,
                "resolution": resolution,
                "size_bytes": len(body),
                "href": href,
            }
        )
        self.files[href] = body
        return href

    def __call__(self, url, **kwargs):
        self.calls.append(url)
        if self.before is not None:
            self.before(url)
        path = url[len(BASE):]
        if path == "/items":
            if isinstance(self.manifest_failure, BaseException):
                raise self.manifest_failure
            if self.manifest_failure is not None:
                return self.manifest_failure
            payload = {"bbox": [0, 0, 1, 1], "crs": "EPSG:3006", "items": self.items}
            return FakeResponse(200, json.dumps(payload))
        if path in self.file_failures:
            return self.file_failures[path]
        body = self.files[path]
        headers = {} if path in self.without_length else {
            "Content-Length": str(len(body))
        }
        return FakeResponse(200, body, headers, chunk=self.chunk)

    def file_calls(self):
        return [url for url in self.calls if "/files/" in url]

    def manifest_calls(self):
        return [url for url in self.calls if url.endswith("/items")]


@pytest.fixture
def server(monkeypatch, tmp_path):
    fake = Server()
    monkeypatch.setattr(client.requests, "get", fake)
    monkeypatch.setattr(data_cache, "cache_dir", tmp_path / "cache")
    monkeypatch.delenv("DTCC_ORTHOPHOTO_URL", raising=False)
    return fake


@pytest.fixture
def no_http(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("no request may be made")

    monkeypatch.setattr(client.requests, "get", forbidden)


def tiff(tmp_path, bounds, *, bands=4, tag=1, pixel=0.5, nodata=None, data=None,
         damaged_block=None):
    """GeoTIFF bytes: band 1 a tag, bands 2-3 column and row, band 4 constant."""
    width = int(round((bounds[2] - bounds[0]) / pixel))
    height = int(round((bounds[3] - bounds[1]) / pixel))
    if data is None:
        cols, rows = np.meshgrid(np.arange(width), np.arange(height))
        layers = [np.full_like(cols, tag), cols % 250 + 1, rows % 250 + 1]
        layers += [np.full_like(cols, 77)] * (bands - 3)
        data = np.stack(layers[:bands]).astype(np.uint8)
    path = tmp_path / "src" / f"{tag}-{bands}-{width}.tif"
    path.parent.mkdir(exist_ok=True)
    profile = {}
    if damaged_block is not None:
        profile = dict(tiled=True, blockxsize=16, blockysize=16, compress="deflate")
    if bands > 3:
        profile["alpha"] = "unspecified"
    with rasterio.open(
        path, "w", driver="GTiff", width=width, height=height, count=bands,
        dtype="uint8", crs="EPSG:3006", nodata=nodata,
        transform=Affine(pixel, 0.0, bounds[0], 0.0, -pixel, bounds[3]), **profile,
    ) as dst:
        dst.write(data)
    body = bytearray(path.read_bytes())
    if damaged_block is not None:
        ifd = struct.unpack_from("<I", body, 4)[0]
        for entry in range(struct.unpack_from("<H", body, ifd)[0]):
            tag_id, kind, count, value = struct.unpack_from(
                "<HHII", body, ifd + 2 + 12 * entry
            )
            if tag_id == 324:
                offset = struct.unpack_from(f"<{count}I", body, value)[damaged_block]
        body[offset : offset + 32] = b"\xff" * 32
    path.unlink()
    return bytes(body)


def call(**kwargs):
    kwargs.setdefault("server_url", BASE)
    kwargs.setdefault("bounds", cell(0))
    return datasets.orthophoto(**kwargs)


def tags(raster):
    data = raster.data
    return set(np.unique(data[..., 0][data[..., 3] == 255]).tolist())


# Arguments


def test_defaults_validate_from_bounds_alone():
    args = OrthophotoArgs(bounds=[0.0, 0.0, 1.0, 1.0])
    assert args.product == "raster"
    assert args.spektraltyp is None
    assert (args.connect_timeout, args.read_timeout) == (10.0, 150.0)
    assert args.format is None and args.max_memory_bytes is None


INVALID = {
    "NaN bounds": {"bounds": [math.nan, 0.0, 1.0, 1.0]},
    "infinite bounds": {"bounds": [0.0, 0.0, math.inf, 1.0]},
    "inverted bounds": {"bounds": [1.0, 0.0, 0.0, 1.0]},
    "zero span": {"bounds": [0.0, 1.0, 1.0, 1.0]},
    "five values": {"bounds": [0.0, 0.0, 1.0, 1.0, 2.0]},
    "inverted Z": {"bounds": [0.0, 0.0, 2.0, 1.0, 1.0, 1.0]},
    "NaN Z": {"bounds": [0.0, 0.0, math.nan, 1.0, 1.0, 1.0]},
    "bool bounds": {"bounds": [True, 0.0, 2.0, 1.0]},
    "product": {"product": "png"},
    "bool year": {"year": True},
    "collection": {"collection": "Orto 2025!"},
    "spektraltyp": {"spektraltyp": ["xyz"]},
    "empty spektraltyp": {"spektraltyp": []},
    "zero resolution": {"resolution": 0},
    "negative resolution": {"resolution": -1.0},
    "NaN resolution": {"resolution": math.nan},
    "bool resolution": {"resolution": True},
    "zero budget": {"max_memory_bytes": 0},
    "bool budget": {"max_memory_bytes": True},
    "fractional budget": {"max_memory_bytes": 1.5},
    "zero connect timeout": {"connect_timeout": 0},
    "infinite read timeout": {"read_timeout": math.inf},
    "bool timeout": {"read_timeout": True},
    "format": {"format": "png"},
    "tiles with resolution": {"product": "tiles", "resolution": 1.0},
    "tiles with format": {"product": "tiles", "format": "tif"},
    "tiles with budget": {"product": "tiles", "max_memory_bytes": 10**9},
    "raster with cir": {"spektraltyp": ["cir"]},
    "raster with pan": {"spektraltyp": ["rgb", "pan"]},
    "unknown argument": {"colour": "red"},
}


@pytest.mark.parametrize("case", INVALID, ids=list(INVALID))
def test_invalid_arguments_are_rejected_before_any_request(no_http, case):
    kwargs = {"bounds": [0.0, 0.0, 1.0, 1.0], "server_url": BASE}
    kwargs.update(INVALID[case])
    with pytest.raises(ValidationError):
        datasets.orthophoto(**kwargs)


def test_valid_argument_forms():
    six = OrthophotoArgs(bounds=[0.0, 0.0, 5.0, 1.0, 1.0, 5.0])
    assert list(six.bounds) == [0.0, 0.0, 5.0, 1.0, 1.0, 5.0]
    tiles = OrthophotoArgs(bounds=[0, 0, 1, 1], product="tiles",
                           spektraltyp=["cir", "pan", "cir"])
    assert tiles.spektraltyp == ["cir", "pan"]
    assert OrthophotoArgs(bounds=[0, 0, 1, 1], collection="orto-o2-2025",
                          year=2025, resolution=0.5).collection == "orto-o2-2025"


def test_bounds_object_with_nan_z_is_rejected(no_http):
    bounds = Bounds(xmin=0.0, ymin=0.0, xmax=1.0, ymax=1.0, zmin=math.nan, zmax=1.0)
    with pytest.raises(ValidationError):
        datasets.orthophoto(bounds=bounds, server_url=BASE)


def test_bounds_object_with_zero_z_span_uses_the_horizontal_extent(server, tmp_path):
    server.add("a", cell(0), tiff(tmp_path, cell(0)))
    raster = datasets.orthophoto(bounds=Bounds(*cell(0)), server_url=BASE)
    assert raster.data.shape == (4, 4, 4)


# Server URL


def test_explicit_url_wins_over_the_environment(server, monkeypatch):
    monkeypatch.setenv("DTCC_ORTHOPHOTO_URL", "http://other.test")
    call()
    assert server.manifest_calls() == [BASE + "/items"]


def test_environment_url_is_used(server, monkeypatch):
    monkeypatch.setenv("DTCC_ORTHOPHOTO_URL", BASE + "/")
    datasets.orthophoto(bounds=cell(0))
    assert server.manifest_calls() == [BASE + "/items"]


@pytest.mark.parametrize("url", [None, "http://user:secret@tiles.test"])
def test_missing_or_credential_url_is_a_local_error(no_http, monkeypatch, url):
    monkeypatch.delenv("DTCC_ORTHOPHOTO_URL", raising=False)
    with pytest.raises(ValueError) as info:
        datasets.orthophoto(bounds=cell(0), server_url=url)
    assert not isinstance(info.value, ValidationError)
    assert "secret" not in str(info.value)


@pytest.mark.parametrize(
    "url",
    [
        "http://us：er:secret@tiles.test",
        "https:user:secret@tiles.test",
        "http://tiles.test?key=secret",
    ],
)
@pytest.mark.parametrize("source", ["explicit", "environment"])
def test_malformed_url_is_a_local_error_without_echo(no_http, monkeypatch, url, source):
    if source == "environment":
        monkeypatch.setenv("DTCC_ORTHOPHOTO_URL", url)
        url_argument = None
    else:
        monkeypatch.delenv("DTCC_ORTHOPHOTO_URL", raising=False)
        url_argument = url
    with pytest.raises(ValueError) as info:
        datasets.orthophoto(bounds=cell(0), server_url=url_argument)
    error = info.value
    assert not isinstance(error, ValidationError)
    assert error.__context__ is None and error.__cause__ is None
    assert "secret" not in "".join(traceback.format_exception(error))


def test_describe_is_deterministic_and_has_no_url(monkeypatch):
    monkeypatch.setenv("DTCC_ORTHOPHOTO_URL", BASE)
    first = json.dumps(datasets.orthophoto.describe(), sort_keys=True)
    monkeypatch.delenv("DTCC_ORTHOPHOTO_URL")
    assert json.dumps(datasets.orthophoto.describe(), sort_keys=True) == first
    assert "tiles.test" not in first and "http" not in first


# One manifest and early refusals


def test_one_manifest_per_call_drives_every_step(server, tmp_path):
    server.add("a", cell(0), tiff(tmp_path, cell(0), tag=1))
    server.add("b", cell(1), tiff(tmp_path, cell(1), tag=2))
    raster = call(bounds=(X0, Y0, X0 + 4.0, Y0 + 2.0))
    assert len(server.manifest_calls()) == 1
    assert len(server.file_calls()) == 2
    health = raster.dataset_context.health
    assert (health["items_listed"], health["items_used"]) == (2, 2)


def test_explicit_resolution_over_budget_is_refused_before_any_request(no_http):
    with pytest.raises(WorkingMemoryError):
        datasets.orthophoto(bounds=(0.0, 0.0, 1000.0, 1000.0), server_url=BASE,
                            resolution=0.1, max_memory_bytes=READ_OVERHEAD_BYTES)


def test_known_manifest_resolutions_over_budget_are_refused_before_downloads(
    server, tmp_path
):
    server.add("a", cell(0), tiff(tmp_path, cell(0)), resolution=0.001)
    # A coarser item alone would fit: the finest listed resolution decides.
    server.add("b", cell(0), tiff(tmp_path, cell(0), tag=2), resolution=1.0)
    with pytest.raises(WorkingMemoryError):
        call(max_memory_bytes=READ_OVERHEAD_BYTES + 1000)
    assert server.file_calls() == []


def test_unknown_manifest_resolution_never_proves_a_fit(server, tmp_path):
    server.add("a", cell(0), tiff(tmp_path, cell(0)), resolution=0.001)
    server.add("b", cell(0), tiff(tmp_path, cell(0), tag=2), resolution=None)
    raster = call()
    assert len(server.file_calls()) == 2
    assert raster.data.shape == (4, 4, 4)


# End to end


def test_default_product_is_an_rgba_mosaic(server, tmp_path):
    server.add("a", cell(0), tiff(tmp_path, cell(0), tag=9))
    raster = call()
    assert isinstance(raster, Raster)
    assert raster.data.shape == (4, 4, 4) and raster.data.dtype == np.uint8
    assert raster.crs == "EPSG:3006"
    assert raster.georef == Affine(0.5, 0.0, X0, 0.0, -0.5, Y0 + 2.0)
    assert (raster.data[..., 3] == 255).all() and tags(raster) == {9}


def test_tiles_product_returns_original_files(server, tmp_path):
    body = tiff(tmp_path, cell(0), bands=3)
    server.add("a", cell(0), body, spektraltyp="cir")
    tiles = call(product="tiles", spektraltyp=["cir"])
    assert isinstance(tiles, RasterTileCollection)
    (tile,) = tiles
    assert tile.path.read_bytes() == body
    assert tile.path.is_relative_to(data_cache.cache_dir)
    assert (tile.id, tile.spektraltyp) == ("a", "cir")
    assert tiles.dataset_context.health["items_used"] == 1


def test_tif_format_returns_georeferenced_rgba_bytes(server, tmp_path):
    server.add("a", cell(0), tiff(tmp_path, cell(0), tag=9))
    payload = call(format="tif")
    assert isinstance(payload, bytes)
    with MemoryFile(payload) as memory, memory.open() as src:
        assert src.count == 4 and src.crs.to_epsg() == 3006
        assert src.colorinterp[-1] == rasterio.enums.ColorInterp.alpha
        assert src.transform == Affine(0.5, 0.0, X0, 0.0, -0.5, Y0 + 2.0)
        assert (src.read(1) == 9).all() and (src.read(4) == 255).all()


def test_a_repeated_call_lists_again_and_reuses_the_cache(server, tmp_path):
    server.add("a", cell(0), tiff(tmp_path, cell(0)))
    call()
    call()
    assert len(server.manifest_calls()) == 2
    assert len(server.file_calls()) == 1


# TIFF budget


def test_tif_write_is_refused_when_output_and_strip_exceed_the_budget(
    server, tmp_path, monkeypatch
):
    # 32 x 32 output pixels with small mosaic blocks: writing needs more memory
    # than mosaicking. No manifest resolution, so no early refusal.
    big = cell(0, size=8.0)
    server.add("a", big, tiff(tmp_path, big, pixel=0.25), resolution=None)
    monkeypatch.setattr(raster_tiles, "STRIP_BYTES", 600)
    with pytest.raises(WorkingMemoryError) as info:
        call(bounds=big, max_memory_bytes=1)
    mosaic_estimate = info.value.estimate_bytes
    assert mosaic_estimate < 32 * 32 * 5 + READ_OVERHEAD_BYTES

    def forbidden(*args, **kwargs):
        raise AssertionError("the GeoTIFF must not be written")

    monkeypatch.setattr(raster_io, "_save_geotif", forbidden)
    with pytest.raises(WorkingMemoryError) as info:
        call(bounds=big, format="tif", max_memory_bytes=mosaic_estimate)
    # The output, one strip of all four bands, and GDAL's allowance.
    assert info.value.estimate_bytes == (
        32 * 32 * 4 + 32 * 32 * 4 + READ_OVERHEAD_BYTES
    )


def test_tif_file_size_does_not_depend_on_the_gdal_cache(
    server, tmp_path, monkeypatch
):
    big = cell(0, size=256.0)
    data = np.random.default_rng(0).integers(1, 255, (4, 512, 512), dtype=np.uint8)
    server.add("a", big, tiff(tmp_path, big, data=data))
    # The GDAL cache the GeoTIFF is written under, far below the 1 MiB output.
    monkeypatch.setattr(orthophoto_module, "READ_OVERHEAD_BYTES", 200_000)
    sizes = []
    real_size = orthophoto_module._written_size

    def recording_size(path):
        sizes.append(real_size(path))
        return sizes[-1]

    monkeypatch.setattr(orthophoto_module, "_written_size", recording_size)
    call(bounds=big, format="tif")
    # Random pixels barely compress; strips rewritten per band would double this.
    assert sizes and sizes[0] <= 512 * 512 * 4 * 1.05


def test_tif_read_is_refused_when_the_file_exceeds_the_budget(
    server, tmp_path, monkeypatch
):
    server.add("a", cell(0), tiff(tmp_path, cell(0)))
    monkeypatch.setattr(orthophoto_module, "_written_size", lambda path: 10**12)

    def forbidden(path):
        raise AssertionError("the GeoTIFF must not be read")

    monkeypatch.setattr(orthophoto_module, "_read_written", forbidden)
    with pytest.raises(WorkingMemoryError) as info:
        call(format="tif")
    assert info.value.estimate_bytes == 10**12


def test_tif_bytes_are_read_after_the_mosaic_is_released(
    server, tmp_path, monkeypatch
):
    server.add("a", cell(0), tiff(tmp_path, cell(0)))
    references = []
    real_build = orthophoto_module.build_mosaic
    real_read = orthophoto_module._read_written

    def recording_build(*args, **kwargs):
        mosaic = real_build(*args, **kwargs)
        references.append(weakref.ref(mosaic.raster.data))
        return mosaic

    def checking_read(path):
        assert references and references[0]() is None
        return real_read(path)

    monkeypatch.setattr(orthophoto_module, "build_mosaic", recording_build)
    monkeypatch.setattr(orthophoto_module, "_read_written", checking_read)
    assert isinstance(call(format="tif"), bytes)


def test_tif_bytes_are_read_after_a_partial_mosaic_is_released(
    server, tmp_path, monkeypatch
):
    big = cell(0, size=4.0)
    server.add("old", big, tiff(tmp_path, big, tag=1, pixel=0.125), year=2010)
    server.add("new", big, tiff(tmp_path, big, tag=2, pixel=0.125, damaged_block=-1),
               year=2020)
    # Small blocks: the damaged block fails after earlier blocks are composited.
    monkeypatch.setattr(raster_tiles, "STRIP_BYTES", 3000)
    references = []
    real_build = orthophoto_module.build_mosaic
    real_read = orthophoto_module._read_written
    real_raster_read = rasterio.io.DatasetReader.read

    def recording_raster_read(self, *args, **kwargs):
        arena = kwargs["out"]
        while arena.base is not None:
            arena = arena.base
        references.append(weakref.ref(arena))
        return real_raster_read(self, *args, **kwargs)

    def recording_build(*args, **kwargs):
        mosaic = real_build(*args, **kwargs)
        assert len(mosaic.failures) == 1 and mosaic.failures[0].pixels > 0
        references.append(weakref.ref(mosaic.raster.data))
        return mosaic

    def checking_read(path):
        assert [ref for ref in references if ref() is not None] == []
        return real_read(path)

    monkeypatch.setattr(rasterio.io.DatasetReader, "read", recording_raster_read)
    monkeypatch.setattr(orthophoto_module, "build_mosaic", recording_build)
    monkeypatch.setattr(orthophoto_module, "_read_written", checking_read)
    # With cyclic collection off, a buffer is released only when its last
    # reference goes.
    gc.disable()
    try:
        assert isinstance(call(bounds=big, format="tif"), bytes)
    finally:
        gc.enable()


def test_tif_without_valid_pixels_is_refused(server, tmp_path):
    empty = np.zeros((4, 4, 4), dtype=np.uint8)
    server.add("a", cell(0), tiff(tmp_path, cell(0), data=empty, nodata=0))
    with pytest.raises(ValueError, match="no valid pixel"):
        call(format="tif")


def test_tif_in_strict_mode_raises_the_upstream_failure_first(server):
    server.manifest_failure = FakeResponse(502, "bad gateway")
    with pytest.raises(DatasetUpstreamError):
        call(format="tif", strict_live=True)


# Strictness


MANIFEST_FAILURES = {
    "credentials": (
        FakeResponse(503, json.dumps({"detail": "LM credentials not configured"})),
        "configuration",
        503,
    ),
    "server error": (FakeResponse(500, "boom"), "http_5xx", 500),
    "timeout": (requests.Timeout("slow"), "timeout", None),
    "connection": (requests.ConnectionError("down"), "connection", None),
}


@pytest.mark.parametrize("case", MANIFEST_FAILURES, ids=list(MANIFEST_FAILURES))
def test_manifest_failures(server, case):
    failure, failure_class, status = MANIFEST_FAILURES[case]
    server.manifest_failure = failure
    with pytest.raises(DatasetUpstreamError) as info:
        call(strict_live=True)
    assert info.value.dataset == "orthophoto"
    assert info.value.failure_class == failure_class
    assert info.value.status_code == status
    assert isinstance(info.value.__cause__, client.OrthophotoClientError)
    result = call()
    assert result.data.shape == ()
    health = result.dataset_context.health
    assert health["status"] == "failed" and health["partial_result"] is True
    assert health["upstream_errors"][0]["failure_class"] == failure_class
    assert health["coverage_complete"] is False


def test_download_failure_strict_stops_and_default_continues(server, tmp_path):
    server.add("a", cell(0), tiff(tmp_path, cell(0), tag=1))
    href = server.add("b", cell(1), tiff(tmp_path, cell(1), tag=2))
    server.add("c", cell(2), tiff(tmp_path, cell(2), tag=3))
    server.file_failures[href] = FakeResponse(404, "missing")
    bounds = (X0, Y0, X0 + 6.0, Y0 + 2.0)
    with pytest.raises(DatasetUpstreamError) as info:
        call(bounds=bounds, strict_live=True)
    assert (info.value.failure_class, info.value.status_code) == ("http_4xx", 404)
    assert len(server.file_calls()) == 2
    server.calls.clear()
    raster = call(bounds=bounds)
    assert tags(raster) == {1, 3}
    health = raster.dataset_context.health
    assert health["status"] == "partial"
    assert (health["items_listed"], health["items_used"]) == (3, 2)
    assert health["upstream_errors"][0]["target"].endswith("/b.tif")


def test_mosaic_source_failure_strict_raises_after_mosaicking(server, tmp_path):
    server.add("a", cell(0), tiff(tmp_path, cell(0), tag=1), year=2010)
    server.add("b", cell(0), tiff(tmp_path, cell(0), bands=3, tag=2), year=2020)
    with pytest.raises(DatasetUpstreamError) as info:
        call(strict_live=True)
    assert info.value.failure_class == "invalid_payload"
    assert info.value.target == "orto-2020/b"
    raster = call()
    assert tags(raster) == {1}
    assert raster.dataset_context.health["status"] == "partial"


@pytest.mark.parametrize("strict", [False, True])
def test_local_errors_raise_in_both_modes(server, tmp_path, monkeypatch, strict):
    server.add("a", cell(0), tiff(tmp_path, cell(0)))
    blocked = tmp_path / "blocked"
    blocked.write_text("")
    monkeypatch.setattr(data_cache, "cache_dir", blocked)
    with pytest.raises(OSError) as info:
        call(strict_live=strict)
    assert not isinstance(info.value, DatasetUpstreamError)


# Partial contributors and coverage


def test_tile_failing_after_contributing_stays_a_source(server, tmp_path, monkeypatch):
    big = cell(0, size=4.0)
    server.add("old", big, tiff(tmp_path, big, tag=1, pixel=0.125), year=2010)
    server.add("new", big, tiff(tmp_path, big, tag=2, pixel=0.125, damaged_block=-1),
               year=2020)
    # Small blocks: the rows above the damaged bottom-right block come first.
    monkeypatch.setattr(raster_tiles, "STRIP_BYTES", 3000)
    raster = call(bounds=big)
    context = raster.dataset_context
    health = context.health
    new_pixels = int((raster.data[..., 0] == 2).sum())
    assert new_pixels > 0
    assert health["status"] == "partial" and health["coverage_complete"] is True
    assert health["items_used"] == 2
    items = {source["id"]: source for source in context.provenance.sources
             if isinstance(source, dict) and source.get("role") == "source_item"}
    assert items["new"]["pixels"] == new_pixels
    assert items["new"]["pixels"] + items["old"]["pixels"] == 32 * 32
    assert "Older imagery (old) overlaps new, which failed." in context.warnings


def test_gaps_without_errors_are_complete_with_incomplete_coverage(server, tmp_path):
    server.add("a", cell(0), tiff(tmp_path, cell(0)))
    raster = call(bounds=(X0, Y0, X0 + 4.0, Y0 + 2.0))
    health = raster.dataset_context.health
    assert health["status"] == "complete" and health["partial_result"] is False
    assert health["coverage_complete"] is False
    assert (health["valid_pixels"], health["total_pixels"]) == (16, 32)
    assert any("50.0%" in warning for warning in raster.dataset_context.warnings)


def test_empty_manifest_is_empty_without_upstream_errors(server):
    raster = call()
    health = raster.dataset_context.health
    assert raster.data.shape == ()
    assert health["status"] == "empty"
    assert health["upstream_error_count"] == 0
    assert health["coverage_complete"] is False


@pytest.mark.parametrize("listed", [False, True])
def test_no_imagery_warning_claims_no_missing_coverage(server, tmp_path, listed):
    if listed:
        # Listed and covering the bounds, but every pixel is declared nodata.
        empty = np.zeros((4, 4, 4), dtype=np.uint8)
        server.add("a", cell(0), tiff(tmp_path, cell(0), nodata=0, data=empty))
    context = call().dataset_context
    assert context.health["items_listed"] == int(listed)
    assert context.warnings == ["No orthophoto has imagery in the requested bounds."]


def test_tiles_fallback_warning_states_overlap_only(server, tmp_path):
    href = server.add("new", cell(0), tiff(tmp_path, cell(0), tag=2), year=2020)
    # The older tile holds only declared nodata: it overlaps but fills nothing.
    empty = np.zeros((4, 4, 4), dtype=np.uint8)
    server.add("old", cell(0), tiff(tmp_path, cell(0), nodata=0, data=empty),
               year=2010)
    server.file_failures[href] = FakeResponse(404, "missing")
    warnings = call(product="tiles").dataset_context.warnings
    assert "Older imagery (old) overlaps new, which failed." in warnings
    assert not any("fills" in warning for warning in warnings)


def test_raster_fallback_warning_states_overlap_without_filling(server, tmp_path):
    two = (X0, Y0, X0 + 4.0, Y0 + 2.0)
    href = server.add("new", cell(0), tiff(tmp_path, cell(0), tag=2), year=2020)
    # The older tile is a source on the right cell and nodata under the failed one.
    data = np.full((4, 4, 8), 9, dtype=np.uint8)
    data[:, :, :4] = 0
    server.add("old", two, tiff(tmp_path, two, nodata=0, data=data), year=2010)
    server.file_failures[href] = FakeResponse(404, "missing")
    raster = call(bounds=two)
    assert (raster.data[:, :4, 3] == 0).all()
    warnings = raster.dataset_context.warnings
    assert "Older imagery (old) overlaps new, which failed." in warnings
    assert not any("fills" in warning for warning in warnings)


# No paths


def test_context_carries_no_cache_paths(server, tmp_path):
    server.add("a", cell(0), tiff(tmp_path, cell(0), tag=1), year=2010)
    server.add("bad", cell(0), tiff(tmp_path, cell(0), bands=3), year=2020)
    server.add("broken", cell(0, size=4.0),
               # Block 2 is the bottom-left one, under the requested cell.
               tiff(tmp_path, cell(0, size=4.0), pixel=0.125, damaged_block=2),
               year=2021)
    href = server.add("html", cell(0), b"<html>bad gateway</html>", year=2022)
    raster = call()
    context = raster.dataset_context
    text = context.model_dump_json()
    assert context.health["upstream_error_count"] == 3
    assert str(data_cache.cache_dir) not in text
    assert str(tmp_path) not in text
    assert ".part" not in text


# Context


def item_sources(context):
    return [source for source in context.provenance.sources
            if isinstance(source, dict) and source.get("role") == "source_item"]


def test_context_records_sources_dates_grid_and_health(server, tmp_path):
    server.add("a", cell(0), tiff(tmp_path, cell(0), tag=1), year=2018)
    server.add("b", cell(1), tiff(tmp_path, cell(1), tag=2), year=2023)
    bounds = (X0, Y0, X0 + 4.0, Y0 + 2.0)
    context = call(bounds=bounds).dataset_context
    assert [s["id"] for s in item_sources(context)] == ["b", "a"]
    assert item_sources(context)[0]["datetime"].startswith("2023-06-01")
    assert context.metadata.collection_period == "2018-06-01/2023-06-01"
    health = context.health
    assert health["requested_bounds"] == list(bounds)
    assert health["bounds"] == list(bounds) and health["resolution"] == 0.5
    assert health["status"] == "complete" and health["coverage_complete"] is True
    assert any("dates" in warning for warning in context.warnings)
    json.dumps(health)


def test_successive_and_nested_calls_have_independent_contexts(server, tmp_path):
    server.add("a", cell(0), tiff(tmp_path, cell(0), tag=1))
    descriptor = datasets.orthophoto
    before = dict(vars(descriptor))
    first = call()
    nested = []

    def call_inside(url):
        if url.endswith("/items") and server.before is not None:
            server.before = None
            nested.append(call(bounds=(X0, Y0, X0 + 4.0, Y0 + 2.0)))

    server.before = call_inside
    second = call()
    server.before = None
    contexts = [first.dataset_context, second.dataset_context,
                nested[0].dataset_context]
    assert len({id(context) for context in contexts}) == 3
    assert contexts[0].health["coverage_complete"] is True
    assert contexts[2].health["coverage_complete"] is False
    assert contexts[1].health == contexts[0].health
    assert dict(vars(descriptor)) == before


# Export


def serve_partial(server, tmp_path):
    """Item a (2021) is served; b (2022, same cell) fails with 503."""
    server.add("a", cell(0), tiff(tmp_path, cell(0), tag=1), year=2021)
    href = server.add("b", cell(0), tiff(tmp_path, cell(0), tag=2), year=2022)
    server.file_failures[href] = FakeResponse(503, "busy")


def export(path, **kwargs):
    kwargs.setdefault("server_url", BASE)
    kwargs.setdefault("bounds", cell(0))
    return datasets.orthophoto.export(path, **kwargs)


def test_export_sidecar_carries_this_calls_context(server, tmp_path):
    serve_partial(server, tmp_path)
    direct = call().dataset_context
    server.calls.clear()
    path = tmp_path / "out" / "o.tif"
    result = export(path)
    assert len(server.manifest_calls()) == 1
    assert result.manifest_path == tmp_path / "out" / "o.manifest.json"
    assert json.loads(result.manifest_path.read_text()) == result.manifest
    context = DatasetContext.model_validate(result.manifest["dataset_context"])
    assert context.health["status"] == "partial"
    assert context.health["upstream_error_count"] == 1
    assert [source["id"] for source in item_sources(context)] == ["a"]
    assert "Older imagery (a) overlaps b, which failed." in context.warnings
    assert context.request.parameters["format"] == "tif"
    assert context.health == direct.health
    assert item_sources(context) == item_sources(direct)
    assert context.warnings == direct.warnings
    with rasterio.open(path) as src:
        assert src.count == 4 and src.colorinterp[-1] == ColorInterp.alpha


def test_export_sidecar_has_no_cache_paths(server, tmp_path):
    server.add("a", cell(0), tiff(tmp_path, cell(0), tag=1), year=2010)
    server.add("bad", cell(0), tiff(tmp_path, cell(0), bands=3), year=2020)
    server.add("broken", cell(0, size=4.0),
               tiff(tmp_path, cell(0, size=4.0), pixel=0.125, damaged_block=2),
               year=2021)
    server.add("html", cell(0), b"<html>bad gateway</html>", year=2022)
    result = export(tmp_path / "out" / "o.tif")
    text = result.manifest_path.read_text()
    assert result.manifest["dataset_context"]["health"]["upstream_error_count"] == 3
    assert str(data_cache.cache_dir) not in text
    assert str(tmp_path) not in text
    assert ".part" not in text


def test_nested_exports_have_independent_sidecars(server, tmp_path):
    server.add("a", cell(0), tiff(tmp_path, cell(0), tag=1))
    descriptor = datasets.orthophoto
    before = dict(vars(descriptor))
    nested = []

    def export_inside(url):
        if url.endswith("/items") and server.before is not None:
            server.before = None
            nested.append(export(tmp_path / "inner.tif",
                                 bounds=(X0, Y0, X0 + 4.0, Y0 + 2.0)))

    server.before = export_inside
    outer = export(tmp_path / "outer.tif")
    server.before = None
    health = [result.manifest["dataset_context"]["health"]
              for result in (outer, nested[0])]
    assert health[0]["coverage_complete"] is True
    assert health[1]["coverage_complete"] is False
    assert dict(vars(descriptor)) == before


def test_publish_uploads_a_sidecar_with_context(server, tmp_path):
    serve_partial(server, tmp_path)
    uploader = RecordingUploader()
    datasets.orthophoto.publish(
        dataset_key="orthophoto-test", format="tif", bounds=cell(0),
        server_url=BASE, output_dir=tmp_path / "package", keep_export=True,
        uploader=uploader,
    )
    (upload,) = uploader.calls
    assert len(server.manifest_calls()) == 1
    assert upload["manifest"]["dataset_context"]["health"]["status"] == "partial"
    _validate_upload_package(
        upload["manifest_path"], upload["files"], manifest=upload["manifest"]
    )


def test_export_without_manifest_writes_only_the_tiff(server, tmp_path):
    server.add("a", cell(0), tiff(tmp_path, cell(0)))
    out = tmp_path / "out"
    result = export(out / "o.tif", manifest=False)
    assert result.manifest is None and result.manifest_path is None
    assert [path.name for path in out.iterdir()] == ["o.tif"]


@pytest.mark.parametrize("strict", [False, True])
def test_failed_export_writes_nothing(server, tmp_path, strict):
    server.manifest_failure = FakeResponse(502, "bad gateway")
    out = tmp_path / "out"
    if strict:
        with pytest.raises(DatasetUpstreamError):
            export(out / "o.tif", strict_live=True)
    else:
        with pytest.raises(ValueError, match="could not be fetched"):
            export(out / "o.tif")
    assert not out.exists() or list(out.iterdir()) == []


class PackageRecorder:
    def __init__(self):
        self.calls = []

    def upload_package(self, *, dataset_key, manifest_path, files, manifest,
                       idempotency_key=None):
        self.calls.append(manifest)
        return {"published": dataset_key}


def test_legacy_package_of_a_partial_mosaic_warns_and_keeps_warnings(
    server, tmp_path
):
    serve_partial(server, tmp_path)
    raster = call()
    uploader = PackageRecorder()
    exports = [
        lambda: raster.export(tmp_path / "pkg"),
        lambda: raster.export(tmp_path / "pkg.dtccpkg"),
        lambda: raster.publish(dataset_key="orthophoto-test", uploader=uploader),
    ]
    for export_package in exports:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            export_package()
        messages = [str(warning.message) for warning in caught
                    if issubclass(warning.category, UserWarning)]
        assert len(messages) == 1
        assert "'partial'" in messages[0] and "canonical=True" in messages[0]
    manifest = json.loads(
        (tmp_path / "pkg" / "manifest.json").read_text(encoding="utf-8")
    )
    with zipfile.ZipFile(tmp_path / "pkg.dtccpkg") as archive:
        assert json.loads(archive.read("manifest.json")) == manifest
    assert "Older imagery (a) overlaps b, which failed." in (
        manifest["presentation"]["warnings"]
    )
    assert "health" not in manifest
    assert uploader.calls[0]["presentation"] == manifest["presentation"]


def test_complete_or_canonical_packages_do_not_warn(server, tmp_path):
    server.add("a", cell(0), tiff(tmp_path, cell(0)))
    complete = call()
    assert complete.dataset_context.health["status"] == "complete"
    serve_partial(server, tmp_path)
    partial = call(bounds=cell(0))
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        complete.export(tmp_path / "complete")
        partial.export(tmp_path / "canonical", canonical=True)


# Export paths: what each one carries

TWO_CELLS = (X0, Y0, X0 + 4.0, Y0 + 2.0)


def serve_black_and_gap(server, tmp_path):
    """Cell 0 holds valid black (RGB 0, NIR 5) and nodata 0; cell 1 is a gap."""
    data = np.zeros((4, 4, 4), dtype=np.uint8)
    data[3, :, :2] = 5
    data[0, 0, 0] = 200
    server.add("a", cell(0), tiff(tmp_path, cell(0), data=data, nodata=0))


def written_tiff(directory):
    (path,) = Path(directory).rglob("*.tif")
    return path


def extracted_tiff(archive_path, directory):
    with zipfile.ZipFile(archive_path) as archive:
        (name,) = [name for name in archive.namelist() if name.endswith(".tif")]
        archive.extract(name, directory)
    return Path(directory) / name


TIFF_PATHS = {
    "format tif": lambda raster, out: (
        out.mkdir(),
        (out / "o.tif").write_bytes(call(bounds=TWO_CELLS, format="tif")),
    )
    and out / "o.tif",
    "descriptor export": lambda raster, out: export(
        out / "o.tif", bounds=TWO_CELLS
    ).path,
    "descriptor publish": lambda raster, out: written_tiff(
        datasets.orthophoto.publish(
            dataset_key="orthophoto-test", format="tif", bounds=TWO_CELLS,
            server_url=BASE, output_dir=out, keep_export=True,
            uploader=RecordingUploader(),
        )
        and out
    ),
    "save with alpha": lambda raster, out: (
        out.mkdir(), raster.save(out / "o.tif", alpha=True)
    )
    and out / "o.tif",
    "save": lambda raster, out: (out.mkdir(), raster.save(out / "o.tif"))
    and out / "o.tif",
    "legacy package directory": lambda raster, out: written_tiff(
        raster.export(out / "pkg").path
    ),
    "legacy package archive": lambda raster, out: extracted_tiff(
        raster.export(out / "pkg.dtccpkg").path, out / "extracted"
    ),
    "canonical package supplement": lambda raster, out: written_tiff(
        raster.export(out / "pkg", canonical=True, format="tif").path
    ),
}
LABELLED_ALPHA = {"format tif", "descriptor export", "descriptor publish",
                  "save with alpha"}


@pytest.mark.parametrize("name", TIFF_PATHS, ids=list(TIFF_PATHS))
def test_every_tiff_path_keeps_the_rgba_pixels(server, tmp_path, name):
    # Generic paths write band 4 without an alpha label: a documented limitation.
    serve_black_and_gap(server, tmp_path)
    raster = call(bounds=TWO_CELLS)
    data = raster.data
    alpha = data[..., 3]
    assert ((data[..., :3] == 0).all(-1) & (alpha == 255)).sum() > 0
    assert (alpha == 0).sum() > 0
    path = TIFF_PATHS[name](raster, tmp_path / "out")
    with rasterio.open(path) as src:
        assert src.count == 4 and set(src.dtypes) == {"uint8"}
        assert src.crs.to_epsg() == 3006
        assert src.transform == raster.georef
        assert src.nodata is None
        np.testing.assert_array_equal(np.moveaxis(src.read(), 0, -1), data)
        assert (src.colorinterp[-1] == ColorInterp.alpha) == (name in LABELLED_ALPHA)
    np.testing.assert_array_equal(dtcc_io.load_raster(path).data, data)


def test_tiff_bytes_carry_no_context_or_paths(server, tmp_path):
    serve_black_and_gap(server, tmp_path)
    payload = call(bounds=TWO_CELLS, format="tif")
    assert type(payload) is bytes
    with MemoryFile(payload) as memory, memory.open() as src:
        assert src.tags() == {"AREA_OR_POINT": "Area"}
        assert src.descriptions == (None, None, None, None)
    assert str(data_cache.cache_dir).encode() not in payload
    assert str(tmp_path).encode() not in payload
    path = tmp_path / "bare.tif"
    path.write_bytes(payload)
    assert dtcc_io.load_raster(path).dataset_context is None


def forbidden_get(*args, **kwargs):
    raise AssertionError("no request may be made")


@pytest.mark.parametrize("target", ["pkg", "pkg.dtccpkg"])
def test_canonical_package_keeps_the_calls_context_without_requests(
    server, tmp_path, monkeypatch, target
):
    serve_partial(server, tmp_path)
    raster = call()
    monkeypatch.setattr(client.requests, "get", forbidden_get)
    package = raster.export(tmp_path / target, canonical=True)
    manifest = package.manifest.model_dump(mode="json")
    assert manifest["health"]["status"] == "partial"
    assert manifest["warnings"] == raster.dataset_context.warnings
    assert str(data_cache.cache_dir) not in json.dumps(manifest)
    loaded = load_model_package(tmp_path / target)
    np.testing.assert_array_equal(loaded.data, raster.data)
    assert loaded.crs == "EPSG:3006"
    assert loaded.dataset_context == raster.dataset_context


@pytest.mark.parametrize("status", ["failed", "empty"])
def test_failed_and_empty_results_export_canonically(server, tmp_path, status):
    if status == "failed":
        server.manifest_failure = FakeResponse(502, "bad gateway")
    raster = call()
    raster.export(tmp_path / "pkg", canonical=True)
    loaded = load_model_package(tmp_path / "pkg")
    assert loaded.dataset_context.health["status"] == status


def test_canonical_packages_of_two_results_keep_their_own_health(server, tmp_path):
    server.add("a", cell(0), tiff(tmp_path, cell(0)))
    whole = call()
    with_gap = call(bounds=TWO_CELLS)
    with_gap.export(tmp_path / "with_gap", canonical=True)
    whole.export(tmp_path / "whole", canonical=True)
    health = [load_model_package(tmp_path / name).dataset_context.health
              for name in ("whole", "with_gap")]
    assert [item["coverage_complete"] for item in health] == [True, False]


def test_canonical_publish_records_health(server, tmp_path):
    serve_partial(server, tmp_path)
    uploader = PackageRecorder()
    call().publish(dataset_key="orthophoto-test", uploader=uploader, canonical=True)
    (manifest,) = uploader.calls
    assert manifest["health"]["status"] == "partial"


class ForbiddenUploader:
    def upload_package(self, **kwargs):
        raise AssertionError("nothing may be uploaded")


EMPTY_TIFF_PATHS = {
    "save": lambda raster, out: raster.save(out / "o.tif"),
    "save with alpha": lambda raster, out: raster.save(out / "o.tif", alpha=True),
    "legacy package directory": lambda raster, out: raster.export(out / "pkg"),
    "legacy package archive": lambda raster, out: raster.export(out / "o.dtccpkg"),
    "canonical package supplement": lambda raster, out: raster.export(
        out / "can", canonical=True, format="tif"
    ),
    "publish": lambda raster, out: raster.publish(
        dataset_key="orthophoto-test", uploader=ForbiddenUploader()
    ),
}


# Legacy exports warn about the dropped "empty" health before they fail.
@pytest.mark.filterwarnings("ignore:A legacy v2 Dataset package:UserWarning")
@pytest.mark.parametrize("name", EMPTY_TIFF_PATHS, ids=list(EMPTY_TIFF_PATHS))
def test_an_empty_result_writes_no_tiff(server, tmp_path, name):
    raster = call()
    assert raster.data.shape == ()
    out = tmp_path / "out"
    out.mkdir()
    with pytest.raises((ValueError, OSError)):
        EMPTY_TIFF_PATHS[name](raster, out)
    assert list(out.rglob("*.tif")) == [] and list(out.rglob("*.dtccpkg")) == []


@pytest.mark.parametrize("path", ["call", "export", "publish"])
def test_tiff_after_a_manifest_failure_is_refused(server, tmp_path, path):
    server.manifest_failure = FakeResponse(502, "bad gateway")
    out = tmp_path / "out"
    with pytest.raises(ValueError, match="could not be fetched"):
        if path == "call":
            call(format="tif")
        elif path == "export":
            export(out / "o.tif")
        else:
            datasets.orthophoto.publish(
                dataset_key="orthophoto-test", format="tif", bounds=cell(0),
                server_url=BASE, uploader=ForbiddenUploader(),
            )
    assert not out.exists()


def test_tile_export_is_refused_before_any_request(no_http, tmp_path):
    kwargs = {"bounds": cell(0), "server_url": BASE, "product": "tiles"}
    with pytest.raises(ValidationError, match='product="tiles"'):
        datasets.orthophoto.export(tmp_path / "t.tif", **kwargs)
    with pytest.raises(ValidationError, match='product="tiles"'):
        datasets.orthophoto.publish(
            dataset_key="orthophoto-test", format="tif",
            uploader=ForbiddenUploader(), **kwargs,
        )
    assert list(tmp_path.iterdir()) == []


NO_SAVER = "Cannot save object: RasterTileCollection. No IO save method"
NOT_CANONICAL = "Canonical object RasterTileCollection is unsupported"

# Each path's own refusal: the exception class and a pattern of its message.
TILE_EXPORTS = {
    "package": (
        lambda tiles, out: export_model_package(tiles, out / "pkg"),
        ValueError, "Cannot infer a safe export format for RasterTileCollection",
    ),
    "package as tif": (
        lambda tiles, out: export_model_package(tiles, out / "pkg", format="tif"),
        ValueError, f"Could not export RasterTileCollection as 'tif': {NO_SAVER}",
    ),
    "package as json": (
        lambda tiles, out: export_model_package(tiles, out / "pkg", format="json"),
        ValueError, f"Could not export RasterTileCollection as 'json': {NO_SAVER}",
    ),
    "canonical package": (
        lambda tiles, out: export_model_package(tiles, out / "pkg", canonical=True),
        NotImplementedError, NOT_CANONICAL,
    ),
    "canonical encoding": (
        lambda tiles, out: exchange.dumps(tiles), NotImplementedError, NOT_CANONICAL,
    ),
    "bytes": (
        lambda tiles, out: DatasetDescriptor.export_to_bytes(tiles, "tif"),
        AttributeError, NO_SAVER,
    ),
    "save_raster": (
        lambda tiles, out: dtcc_io.save_raster(tiles, out / "t.tif"),
        RuntimeError, "Unable to save raster; type .*RasterTileCollection.* not "
        "supported",
    ),
}


@pytest.mark.parametrize("name", TILE_EXPORTS, ids=list(TILE_EXPORTS))
def test_tile_collection_export_is_refused_on_every_path(
    server, tmp_path, monkeypatch, name
):
    server.add("a", cell(0), tiff(tmp_path, cell(0)))
    tiles = call(product="tiles")
    monkeypatch.setattr(client.requests, "get", forbidden_get)
    out = tmp_path / "out"
    out.mkdir()
    export, expected, message = TILE_EXPORTS[name]
    with pytest.raises(expected, match=message) as info:
        export(tiles, out)
    assert type(info.value) is expected
    assert str(data_cache.cache_dir) not in str(info.value)
    assert [path for path in out.rglob("*") if path.is_file()] == []


# Progress


@pytest.fixture
def events(monkeypatch):
    """Every progress state the call reports, unthrottled."""
    recorded = []
    monkeypatch.setattr(
        orthophoto_module, "ProgressTracker",
        partial(progress_module.ProgressTracker, callback=recorded.append,
                mode="callback", min_update_interval=0),
    )
    return recorded


def subscribe(monkeypatch, subscriber):
    monkeypatch.setattr(
        orthophoto_module, "ProgressTracker",
        partial(progress_module.ProgressTracker, callback=subscriber,
                mode="callback", min_update_interval=0),
    )


def phases_seen(events):
    seen = []
    for state in events:
        if state["phase"] and (not seen or seen[-1] != state["phase"]):
            seen.append(state["phase"])
    return seen


def transfer_states(events):
    return [state for state in events if state["phase"] == "transfer"]


def random_tile(tmp_path, bounds, pixel=0.125):
    width = int(round((bounds[2] - bounds[0]) / pixel))
    height = int(round((bounds[3] - bounds[1]) / pixel))
    data = np.random.default_rng(0).integers(1, 255, (4, height, width),
                                             dtype=np.uint8)
    return tiff(tmp_path, bounds, data=data, pixel=pixel)


@pytest.mark.parametrize(
    "kwargs, expected",
    [
        ({}, ["discovery", "transfer", "headers", "mosaic"]),
        ({"format": "tif"}, ["discovery", "transfer", "headers", "mosaic", "export"]),
        ({"product": "tiles"}, ["discovery", "transfer"]),
    ],
    ids=["raster", "tif", "tiles"],
)
def test_progress_phases_run_in_order_to_100(server, tmp_path, events, kwargs,
                                             expected):
    server.add("a", cell(0), tiff(tmp_path, cell(0)))
    call(**kwargs)
    assert phases_seen(events) == expected
    percents = [state["percent"] for state in events]
    assert percents == sorted(percents)
    assert events[0]["percent"] == 0.0 and events[-1]["percent"] == 100.0


def test_transfer_progress_follows_content_length(server, tmp_path, events):
    big = cell(0, size=8.0)
    server.add("a", big, random_tile(tmp_path, big))
    server.chunk = 256
    call(bounds=big)
    within = [state["phases"]["transfer"]["progress"]
              for state in transfer_states(events)]
    assert len([value for value in within if 0 < value < 100]) > 5


def test_transfer_without_content_length_reports_bytes_only(
    server, tmp_path, events, monkeypatch
):
    first, second = cell(0, size=8.0), cell(1, size=8.0)
    href = server.add("a", first, random_tile(tmp_path, first))
    server.add("b", second, tiff(tmp_path, second, pixel=0.125))
    server.chunk = 256
    server.without_length.add(href)
    monkeypatch.setattr(orthophoto_module, "BYTES_MESSAGE_STEP", 4096)
    call(bounds=(first[0], first[1], second[2], second[3]))
    # One report per BYTES_MESSAGE_STEP received, each naming the item; the
    # phase holds at the item's start while the length is unknown.
    reports = [state for state in transfer_states(events)
               if "orto-2020/a" in (state["message"] or "")]
    assert len(reports) > 2 and all("MiB" in state["message"] for state in reports)
    assert {state["phases"]["transfer"]["progress"] for state in reports} == {0.0}


def test_mosaic_progress_moves_within_the_mosaic_phase(
    server, tmp_path, events, monkeypatch
):
    big = cell(0, size=4.0)
    server.add("a", big, tiff(tmp_path, big, pixel=0.125))
    # Small blocks: the 32 x 32 output is composited in many steps.
    monkeypatch.setattr(raster_tiles, "STRIP_BYTES", 3000)
    call(bounds=big)
    within = [state["phases"]["mosaic"]["progress"] for state in events
              if state["phase"] == "mosaic"]
    assert len([value for value in within if 0 < value < 100]) > 2


def test_a_long_transfer_reports_about_once_per_percent(server, tmp_path, events):
    big = cell(0, size=16.0)
    href = server.add("a", big, random_tile(tmp_path, big))
    server.chunk = 100
    call(bounds=big)
    assert len(server.files[href]) / server.chunk > 600
    assert len(transfer_states(events)) <= 105


RASTER_PHASES = {"discovery": 0.05, "transfer": 0.6, "headers": 0.05, "mosaic": 0.25}
TILE_BYTES = 659_037_310  # an original 2025 tile


def drive_transfer(sizes, step, with_length=True):
    """Transfer states for logical downloads of ``sizes`` bytes in ``step``s.

    Each state carries the call's unrounded percent as ``exact``.
    """
    recorded = []

    def record(state):
        exact = progress_module.get_progress().state.overall_percent
        recorded.append({**state, "exact": exact})

    items = [listed_item("orto-2025", f"t{index}") for index in range(len(sizes))]
    with progress_module.ProgressTracker(
        phases=RASTER_PHASES, callback=record, mode="callback",
        min_update_interval=0,
    ) as tracker:
        with tracker.phase("discovery"):
            pass
        progress = orthophoto_module._transfer_progress(items, tracker, [], [])
        with tracker.phase("transfer"):
            count = len(sizes)
            for index, size in enumerate(sizes):
                progress(index, count, 0, None)
                for written in [*range(0, size, step), size]:
                    progress(index, count, written, size if with_length else None)
            progress(count, count, 0, None)
    return transfer_states(recorded)


def whole_percents(states):
    return {int(state["percent"]) for state in states}


# Besides one report per whole percent of the call: the phase's entry and exit
# and its completion message.
BOUNDARY_STATES = 3


@pytest.mark.parametrize(
    "sizes, step",
    [([TILE_BYTES], 1024**2), ([1000] * 300, 100)],
    ids=["one full-size tile", "300 small tiles"],
)
def test_transfer_reports_once_per_whole_percent_of_the_call(sizes, step):
    states = drive_transfer(sizes, step)
    assert len(whole_percents(states)) > 60
    assert len(states) <= len(whole_percents(states)) + BOUNDARY_STATES
    count = len(sizes)
    assert states[-1]["message"] == f"{count} of {count} original orthophotos ready"
    # Between the phase's entry and its completion message, each report is the
    # first in a new whole percent of the call (rounded for float noise at exact
    # boundaries).
    reports = [int(round(state["exact"], 6)) for state in states[1:-2]]
    assert reports == sorted(set(reports))


def test_a_report_starts_each_whole_percent_of_the_call():
    recorded = []
    with progress_module.ProgressTracker(
        phases=RASTER_PHASES, callback=recorded.append, mode="callback",
        min_update_interval=0,
    ) as tracker:
        with tracker.phase("discovery"):
            pass
        with tracker.phase("transfer"):
            report = orthophoto_module._Reporter(tracker, "transfer")
            # Transfer starts at 5.26% of the call and spans 63.2%: 1.4% and
            # 1.9% of it are 6.15% and 6.46% of the call, one whole percent.
            for percent in (0.0, 1.4, 1.9, 3.0):
                report(percent, f"at {percent}")
    messages = [state["message"] for state in transfer_states(recorded)]
    assert messages[1:-1] == ["at 0.0", "at 1.4", "at 3.0"]


def test_a_full_size_transfer_without_length_reports_every_byte_step():
    states = drive_transfer([TILE_BYTES], 1024**2, with_length=False)
    steps = TILE_BYTES // orthophoto_module.BYTES_MESSAGE_STEP
    assert steps < len(states) <= steps + 1 + BOUNDARY_STATES


def test_a_many_block_mosaic_reports_once_per_whole_percent_of_the_call(
    server, tmp_path, events, monkeypatch
):
    big = cell(0, size=8.0)
    server.add("a", big, tiff(tmp_path, big, pixel=0.125))
    monkeypatch.setattr(raster_tiles, "STRIP_BYTES", 600)
    call(bounds=big)
    mosaic = [state for state in events if state["phase"] == "mosaic"]
    assert len(mosaic) > 5
    done, total = mosaic[-1]["message"].removeprefix("mosaic: ").split(" of ")
    assert done == total
    assert len(events) <= (
        len(whole_percents(events)) + BOUNDARY_STATES * len(phases_seen(events))
    )


def test_transfer_completion_counts_only_ready_originals(server, tmp_path, events):
    server.add("a", cell(0), tiff(tmp_path, cell(0)))
    href = server.add("b", cell(1), tiff(tmp_path, cell(1)))
    server.file_failures[href] = FakeResponse(404, "missing")
    call(bounds=TWO_CELLS)
    assert transfer_states(events)[-1]["message"] == (
        "1 of 2 original orthophotos ready"
    )


def test_cached_call_completes_the_transfer_without_requests(server, tmp_path,
                                                             events):
    server.add("a", cell(0), tiff(tmp_path, cell(0)))
    call()
    events.clear()
    server.calls.clear()
    call()
    assert server.file_calls() == []
    assert transfer_states(events)[-1]["phases"]["transfer"]["completed"] is True
    assert events[-1]["percent"] == 100.0


def test_download_failures_and_progress(server, tmp_path, events):
    server.add("a", cell(0), tiff(tmp_path, cell(0)))
    href = server.add("b", cell(1), tiff(tmp_path, cell(1)))
    server.file_failures[href] = FakeResponse(404, "missing")
    call(bounds=TWO_CELLS)
    assert phases_seen(events) == ["discovery", "transfer", "headers", "mosaic"]
    assert events[-1]["percent"] == 100.0
    events.clear()
    with pytest.raises(DatasetUpstreamError):
        call(bounds=TWO_CELLS, strict_live=True)
    assert phases_seen(events) == ["discovery", "transfer"]


def test_overlong_body_never_advances_past_its_item(server, tmp_path, events):
    href = server.add("a", cell(0), tiff(tmp_path, cell(0)))
    server.add("b", cell(1), tiff(tmp_path, cell(1)))
    body = server.files[href]
    server.file_failures[href] = FakeResponse(
        200, body + b"x" * 4096, {"Content-Length": str(len(body))}, chunk=64
    )
    call(bounds=TWO_CELLS)
    first = [state["phases"]["transfer"]["progress"]
             for state in transfer_states(events)
             if "orto-2020/a" in (state["message"] or "")]
    assert first and max(first) <= 50.0


def test_zero_content_length_is_a_recorded_failure(server, tmp_path, events):
    href = server.add("a", cell(0), tiff(tmp_path, cell(0)))
    server.file_failures[href] = FakeResponse(200, b"", {"Content-Length": "0"})
    health = call().dataset_context.health
    assert health["upstream_errors"][0]["failure_class"] == "invalid_payload"


def test_a_nested_call_leaves_the_outer_tracker_alone(server, tmp_path, events):
    server.add("a", cell(0), tiff(tmp_path, cell(0)))
    with progress_module.ProgressTracker(phases={"outer": 1.0}, mode="silent") as outer:
        with outer.phase("outer"):
            call()
            assert outer.state.phases["outer"].progress == 0.0
            assert progress_module.get_progress() is outer
    assert progress_module.get_progress() is None
    assert events[0]["percent"] == 0.0


def test_a_memory_refusal_closes_the_header_phase(server, tmp_path, events):
    server.add("a", cell(0), tiff(tmp_path, cell(0)), resolution=None)
    with pytest.raises(WorkingMemoryError):
        call(max_memory_bytes=1)
    assert events[-1]["phase"] == "headers"
    assert events[-1]["phases"]["headers"]["completed"] is True
    assert progress_module.get_progress() is None


def test_all_refused_headers_still_end_at_100(server, tmp_path, events):
    server.add("a", cell(0), tiff(tmp_path, cell(0), bands=3))
    raster = call()
    assert raster.dataset_context.health["status"] == "failed"
    assert phases_seen(events)[-1] == "mosaic"
    assert events[-1]["percent"] == 100.0


@pytest.mark.parametrize("phase", ["discovery", "transfer", "headers", "mosaic",
                                   "export"])
def test_the_first_subscriber_error_propagates_from_every_phase(
    server, tmp_path, monkeypatch, phase
):
    server.add("a", cell(0), tiff(tmp_path, cell(0)))
    seen = []
    errors = []

    def subscriber(state):
        seen.append(state["phase"])
        # The second state of the phase is reported inside it; every later call
        # (the phase's cleanup included) raises again.
        if errors or seen.count(phase) == 2:
            errors.append(RuntimeError(f"call {len(errors) + 1}"))
            raise errors[-1]

    subscribe(monkeypatch, subscriber)
    with pytest.raises(RuntimeError) as info:
        call(format="tif")
    assert info.value is errors[0]
    assert progress_module.get_progress() is None


def client_error():
    return client.OrthophotoClientError(
        "subscriber", operation="download", target="x", failure_class="connection"
    )


@pytest.mark.parametrize("strict", [False, True], ids=["default", "strict"])
def test_a_transfer_subscribers_client_error_escapes_unchanged(
    server, tmp_path, monkeypatch, strict
):
    server.add("a", cell(0), tiff(tmp_path, cell(0)))
    seen = []
    errors = []

    def subscriber(state):
        seen.append(state["phase"])
        if errors:
            errors.append(RuntimeError("cleanup"))
            raise errors[-1]
        if seen.count("transfer") == 2:
            errors.append(client_error())
            raise errors[-1]

    subscribe(monkeypatch, subscriber)
    with pytest.raises(client.OrthophotoClientError) as info:
        call(strict_live=strict)
    assert info.value is errors[0]


def test_a_discovery_subscribers_client_error_is_not_a_failed_result(
    server, tmp_path, monkeypatch
):
    server.add("a", cell(0), tiff(tmp_path, cell(0)))
    errors = []

    def subscriber(state):
        if state["phase"] == "discovery":
            errors.append(client_error())
            raise errors[-1]

    subscribe(monkeypatch, subscriber)
    with pytest.raises(client.OrthophotoClientError) as info:
        call()
    assert info.value is errors[0]


# Live test helpers (the live module is imported, so its test is not collected
# here)


def test_live_module_is_gated():
    assert live.pytestmark.name == "live"
    assert Path(live.__file__).parent.parts[-3:] == ("tests", "datasets", "live")


@pytest.mark.parametrize("value", [None, "", "   "])
def test_live_test_skips_without_a_server_url(monkeypatch, value):
    if value is None:
        monkeypatch.delenv("DTCC_ORTHOPHOTO_URL", raising=False)
    else:
        monkeypatch.setenv("DTCC_ORTHOPHOTO_URL", value)
    with pytest.raises(pytest.skip.Exception) as info:
        live._require_server_url()
    assert "DTCC_ORTHOPHOTO_URL" in str(info.value)
    assert "Chalmers" in str(info.value)


def test_live_test_runs_with_a_server_url(monkeypatch):
    monkeypatch.setenv("DTCC_ORTHOPHOTO_URL", BASE)
    assert live._require_server_url() is None


def upstream_error(failure_class):
    return DatasetUpstreamError(
        dataset="orthophoto", operation="list", target="x",
        failure_class=failure_class, status_code=None, message="m",
    )


@pytest.mark.parametrize("failure_class", ["connection", "timeout", "http_5xx"])
def test_live_test_skips_transient_failures(failure_class):
    with pytest.raises(pytest.skip.Exception, match=failure_class):
        live._skip_if_transient(upstream_error(failure_class))


@pytest.mark.parametrize("failure_class",
                         ["http_4xx", "invalid_payload", "configuration"])
def test_live_test_fails_on_other_failures(failure_class):
    error = upstream_error(failure_class)
    # pytest.raises does not catch a skip, so a skip would pass silently.
    try:
        live._skip_if_transient(error)
    except pytest.skip.Exception:
        pytest.fail(f"{failure_class} must fail the live test, not skip it")
    except DatasetUpstreamError as raised:
        assert raised is error
    else:
        pytest.fail("the failure must be re-raised")


def listed_item(collection, item_id):
    return client.OrthophotoItem(
        id=item_id, collection=collection, bbox=cell(0),
        datetime=datetime(2025, 6, 1, tzinfo=timezone.utc), spektraltyp="rgbi",
        href=f"/files/{collection}/{item_id}.tif", resolution=0.16, size_bytes=None,
    )


@pytest.mark.parametrize(
    "items, allowed",
    [
        ([("orto-o2-2025", "o65775_6725_25_mr25")], True),
        ([("orto-o2-2025", "o65775_6725_25_mr25"), ("orto-e-2008", "o65750_6700_50")],
         False),
        ([("orto-o2-2025", "o65775_6725_25_mr26")], False),
        ([("orto-o2-2024", "o65775_6725_25_mr25")], False),
    ],
    ids=["pinned", "extra item", "other id", "other collection"],
)
def test_live_download_guard_admits_only_the_pinned_item(items, allowed):
    calls = []
    guarded = live.guard_downloads(lambda items, **kwargs: calls.append(items))
    listed = [listed_item(*item) for item in items]
    if allowed:
        guarded(listed, server_url=BASE)
        assert calls == [listed]
    else:
        with pytest.raises(AssertionError):
            guarded(listed, server_url=BASE)
        assert calls == []


THREE_CELLS = (X0, Y0, X0 + 6.0, Y0 + 2.0)


@pytest.fixture
def live_like(server, tmp_path):
    """A strict mosaic with mixed 0.5 and 0.25 m sources and a gap, and tiles."""
    server.add("a", cell(0), tiff(tmp_path, cell(0), tag=1))
    server.add("b", cell(1), tiff(tmp_path, cell(1), tag=2, pixel=0.25),
               resolution=0.25)
    raster = call(bounds=THREE_CELLS, strict_live=True)
    server.calls.clear()
    tiles = call(bounds=THREE_CELLS, product="tiles", strict_live=True)
    return raster, tiles, list(server.calls)


def with_changes(raster, data=None, **health):
    context = raster.dataset_context
    changed = Raster(
        data=raster.data.copy() if data is None else data, georef=raster.georef,
        crs=raster.crs, nodata=raster.nodata,
    )
    return attach_dataset_context(
        changed,
        context.model_copy(update={"health": {**context.health, **health}}),
    )


def test_live_checks_pass_on_a_valid_result(live_like):
    raster, tiles, requests_made = live_like
    assert raster.dataset_context.health["coverage_complete"] is False
    live.check_mosaic(raster, THREE_CELLS, data_cache.cache_dir)
    live.check_tiles(tiles, THREE_CELLS, data_cache.cache_dir, requests_made)


def broken_mosaics(raster):
    health = raster.dataset_context.health
    alpha = raster.data.copy()
    # A transparent pixel, so neither the valid count nor the colour check sees it.
    row, col = np.argwhere(alpha[..., 3] == 0)[0]
    alpha[row, col, 3] = 128
    colour = raster.data.copy()
    row, col = np.argwhere(colour[..., 3] == 0)[0]
    colour[row, col, :3] = 7
    xmin, ymin, xmax, ymax = health["bounds"]
    step = health["resolution"]
    context = raster.dataset_context
    leaked = attach_dataset_context(
        Raster(data=raster.data.copy(), georef=raster.georef, crs=raster.crs,
               nodata=raster.nodata),
        context.model_copy(update={"warnings": context.warnings
                                   + [str(data_cache.cache_dir)]}),
    )
    return {
        "alpha 128": with_changes(raster, data=alpha),
        "colour under alpha 0": with_changes(raster, data=colour),
        "valid pixels off": with_changes(raster,
                                         valid_pixels=health["valid_pixels"] + 1),
        "total pixels off": with_changes(raster,
                                         total_pixels=health["total_pixels"] + 1),
        "grid shifted": with_changes(
            raster, bounds=[xmin + step, ymin, xmax + step, ymax]
        ),
        "partial": with_changes(raster, status="partial"),
        "cache path in warnings": leaked,
    }


@pytest.mark.parametrize("case", ["alpha 128", "colour under alpha 0",
                                  "valid pixels off", "total pixels off",
                                  "grid shifted", "partial",
                                  "cache path in warnings"])
def test_live_mosaic_check_rejects_broken_results(live_like, case):
    raster, _, _ = live_like
    with pytest.raises(AssertionError):
        live.check_mosaic(broken_mosaics(raster)[case], THREE_CELLS,
                          data_cache.cache_dir)


@pytest.mark.parametrize("case", ["missing file", "wrong size", "file request",
                                  "second manifest"])
def test_live_tiles_check_rejects_broken_results(live_like, case):
    _, tiles, requests_made = live_like
    if case == "missing file":
        tiles[0].path.unlink()
    elif case == "wrong size":
        changed = [dataclasses.replace(tiles[0], size_bytes=tiles[0].size_bytes + 1)]
        tiles = attach_dataset_context(
            RasterTileCollection(tiles=changed + list(tiles)[1:]),
            tiles.dataset_context,
        )
    elif case == "file request":
        requests_made = requests_made + [BASE + "/files/orto-2020/a.tif"]
    else:
        requests_made = requests_made + [BASE + "/items"]
    with pytest.raises(AssertionError):
        live.check_tiles(tiles, THREE_CELLS, data_cache.cache_dir, requests_made)


@pytest.fixture
def live_server(server, tmp_path, monkeypatch):
    """The fake server listing only the live test's pinned item."""
    monkeypatch.setenv("DTCC_ORTHOPHOTO_URL", BASE)
    ((collection, item_id),) = live.PINNED
    server.add(item_id, live.BOUNDS, tiff(tmp_path, live.BOUNDS, pixel=5.0),
               year=live.YEAR, collection=collection, resolution=5.0)
    return server


def test_live_test_body_passes_offline(live_server, monkeypatch):
    live.test_orthophoto_live_mosaic_and_tiles(monkeypatch)
    assert len(live_server.manifest_calls()) == 2
    assert len(live_server.file_calls()) == 1


LIVE_FAILURES = {
    "connection": lambda: requests.ConnectionError("down"),
    "timeout": lambda: requests.Timeout("slow"),
    "http_5xx": lambda: FakeResponse(503, "busy"),
    "http_4xx": lambda: FakeResponse(404, "missing"),
}


@pytest.mark.parametrize("call_number", [1, 2], ids=["mosaic call", "tiles call"])
@pytest.mark.parametrize("failure_class", LIVE_FAILURES)
def test_live_test_classifies_failures_of_both_calls(
    live_server, monkeypatch, call_number, failure_class
):
    failure = LIVE_FAILURES[failure_class]()

    def before(url):
        if url.endswith("/items") and len(live_server.manifest_calls()) == call_number:
            if isinstance(failure, BaseException):
                raise failure
            live_server.manifest_failure = failure

    live_server.before = before
    # pytest.raises does not catch a skip, so both outcomes are caught here.
    try:
        live.test_orthophoto_live_mosaic_and_tiles(monkeypatch)
    except pytest.skip.Exception as skip:
        assert failure_class in live.TRANSIENT, f"{failure_class} was skipped"
        assert failure_class in str(skip)
    except DatasetUpstreamError as error:
        assert failure_class not in live.TRANSIENT, f"{failure_class} failed"
        assert error.failure_class == failure_class
    else:
        pytest.fail("the failure must skip or fail the live test")
    assert len(live_server.manifest_calls()) == call_number


# Registration


def test_orthophoto_is_registered_with_tif_output():
    assert isinstance(datasets.orthophoto, OrthophotoDataset)
    assert datasets.orthophoto.list_supported_formats() == ["tif"]
    described = datasets.orthophoto.describe()
    assert described["data_category"] == "raw"
    assert described["result_kind"] == "raster"


def test_orthophoto_passes_the_dataset_qa_audit():
    from dtcc_core.datasets.qa import audit_dataset

    report = audit_dataset(datasets.orthophoto)
    assert report.failures() == ()
