import gc
import importlib
import json
import math
import struct
import traceback
import weakref
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
from dtcc_core.datasets.dataset import DatasetUpstreamError
from dtcc_core.datasets.publish import _validate_upload_package
from dtcc_core.datasets.schema import DatasetContext
from dtcc_core.datasets.orthophoto import OrthophotoArgs, OrthophotoDataset
from dtcc_core.io import raster as raster_io
from dtcc_core.io import raster_tiles
from dtcc_core.io.data import cache as data_cache
from dtcc_core.io.data import orthophoto as client
from dtcc_core.io.raster_tiles import READ_OVERHEAD_BYTES, WorkingMemoryError
from dtcc_core.model import Bounds, Raster
from dtcc_core.model.values.raster_tiles import RasterTileCollection
from tests.datasets.test_dataset_publish import RecordingUploader

# The package attribute is the registered dataset; the module is imported directly.
orthophoto_module = importlib.import_module("dtcc_core.datasets.orthophoto")

BASE = "http://tiles.test:8000"
X0, Y0 = 672500.0, 6577500.0


def cell(i, j=0, size=2.0):
    return (X0 + i * size, Y0 + j * size, X0 + (i + 1) * size, Y0 + (j + 1) * size)


# Fake tile server


class FakeResponse:
    def __init__(self, status_code=200, body=b"", headers=None):
        self.status_code = status_code
        self.body = body if isinstance(body, bytes) else body.encode()
        self.headers = CaseInsensitiveDict(headers or {})

    def iter_content(self, chunk_size=1):
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
        return FakeResponse(200, body, {"Content-Length": str(len(body))})

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
