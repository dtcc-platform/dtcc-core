import hashlib
import json
import re
import shutil
import sys
import zipfile
from dataclasses import replace
from pathlib import Path, PureWindowsPath

import numpy as np
import pytest
import rasterio

from dtcc_core.datasets.package import export_model_package, load_model_package
from dtcc_core.datasets.publish import DatasetPackageError, DatasetUploadClient
from dtcc_core.io.model import load_model, save_model
from dtcc_core.io.raster_tiles import WorkingMemoryError, build_mosaic, load_tile
from dtcc_core.model import Bounds
from dtcc_core.model.values.raster_tiles import RasterTileCollection
from tests.io.test_raster_tiles import CELL, band_values, write_tiff
from tests.model.test_raster_tiles import METADATA, REFERENCES, make_context, make_tile

FORMS = ["pkg", "pkg.dtccpkg"]
BUDGET = 2 * 1024**3


def cached_collection(cache):
    """Two tiles from two collections that share one tile id, in a local cache."""
    tiles = []
    for collection, spektraltyp, bands, nodata in (
        ("orto-o2-2025", "rgbi", 4, None),
        ("orto-o2-2017", "rgb", 3, 0),
    ):
        tile = make_tile(cache, collection=collection)
        write_tiff(tile.path, band_values(bands), nodata=nodata)
        tiles.append(replace(tile, extent=Bounds(*CELL), spektraltyp=spektraltyp))
    collection = RasterTileCollection(tiles=tiles)
    collection.dataset_context = make_context()
    return collection


def assert_same_rasters(collection, expected):
    for tile, raster in zip(collection, expected, strict=True):
        loaded = tile.load()
        np.testing.assert_array_equal(loaded.data, raster.data)
        assert loaded.georef == raster.georef
        assert loaded.crs == raster.crs
        np.testing.assert_equal(loaded.nodata, raster.nodata)


def metadata(collection):
    return [tuple(getattr(tile, name) for name in METADATA) for tile in collection]


@pytest.fixture
def opened(monkeypatch):
    """Record the path of every Rasterio open."""
    paths = []
    real_open = rasterio.open

    def recording_open(path, *args, **kwargs):
        paths.append(path)
        return real_open(path, *args, **kwargs)

    monkeypatch.setattr(rasterio, "open", recording_open)
    return paths


@pytest.mark.parametrize("form", FORMS)
def test_tiles_reload_and_reexport_after_the_cache_is_gone(
    tmp_path, monkeypatch, opened, form
):
    cache = tmp_path / "cache"
    original = cached_collection(cache)
    expected = [tile.load() for tile in original]
    digests = [hashlib.sha256(tile.path.read_bytes()).hexdigest() for tile in original]
    sizes = [tile.path.stat().st_size for tile in original]
    target = tmp_path / form
    package = original.export(target, canonical=True)
    shutil.rmtree(cache)

    tiles = [
        artifact for artifact in package.manifest.artifacts if artifact.role == "tile"
    ]
    assert [artifact.path for artifact in tiles] == REFERENCES
    assert [artifact.sha256 for artifact in tiles] == digests
    assert [artifact.size for artifact in tiles] == sizes
    assert str(tmp_path) not in package.manifest.model_dump_json()

    restored = load_model_package(target)
    assert type(restored) is RasterTileCollection
    assert metadata(restored) == metadata(original)
    assert not any(tile._unbound for tile in restored)
    if form == "pkg":
        assert all(tile.archive is None for tile in restored)
        assert [tile.path for tile in restored] == [
            (target / reference).resolve() for reference in REFERENCES
        ]
    else:
        with zipfile.ZipFile(target) as archive:
            compression = {i.filename: i.compress_type for i in archive.infolist()}
        assert compression == {
            "manifest.json": zipfile.ZIP_DEFLATED,
            "artifacts/model.dtcc": zipfile.ZIP_DEFLATED,
            **{reference: zipfile.ZIP_STORED for reference in REFERENCES},
        }
        assert [tile.path for tile in restored] == [Path(ref) for ref in REFERENCES]
        assert all(tile.archive == target.resolve() for tile in restored)
    opened.clear()
    assert_same_rasters(restored, expected)
    if form == "pkg":
        assert opened == [tile.path for tile in restored]
    else:
        assert opened == [
            f"/vsizip/{{{target.resolve()}}}/{reference}" for reference in REFERENCES
        ]
    assert restored.dataset_context == original.dataset_context

    if form == "pkg.dtccpkg":

        def forbidden(*args, **kwargs):
            raise AssertionError("archive members are read through zipfile")

        monkeypatch.setattr(shutil, "copyfile", forbidden)
    again = tmp_path / "again.dtccpkg"
    restored.export(again, canonical=True)
    reloaded = load_model_package(again)
    assert metadata(reloaded) == metadata(original)
    assert_same_rasters(reloaded, expected)
    assert reloaded.dataset_context == original.dataset_context


@pytest.fixture
def archive_tile(tmp_path):
    original = cached_collection(tmp_path / "cache")
    original.export(tmp_path / "pkg.dtccpkg", canonical=True)
    return load_model_package(tmp_path / "pkg.dtccpkg")[0]


def test_tiny_budget_is_refused_before_archive_pixels_are_read(
    archive_tile, monkeypatch
):
    def forbidden(self, *args, **kwargs):
        raise AssertionError("pixels must not be read")

    monkeypatch.setattr(rasterio.io.DatasetReader, "read", forbidden)
    with pytest.raises(WorkingMemoryError):
        archive_tile.load(max_memory_bytes=1)


@pytest.mark.parametrize(
    "member",
    [
        "/vsizip/pkg.dtccpkg/artifacts/tiles/orto-o2-2025/0-o65775_6725_25_mr25.tif",
        "artifacts\\tiles\\orto-o2-2025\\0-o65775_6725_25_mr25.tif",
        "artifacts/tiles/../tiles/orto-o2-2025/0-o65775_6725_25_mr25.tif",
        "/artifacts/tiles/orto-o2-2025/0-o65775_6725_25_mr25.tif",
    ],
)
def test_archive_member_must_be_a_portable_package_path(archive_tile, member):
    with pytest.raises(ValueError, match="not a portable package path"):
        load_tile(replace(archive_tile, path=member), max_memory_bytes=BUDGET)


def test_archive_must_be_an_existing_local_file(archive_tile, tmp_path):
    missing = replace(archive_tile, archive=tmp_path / "missing.dtccpkg")
    with pytest.raises(FileNotFoundError):
        load_tile(missing, max_memory_bytes=BUDGET)
    with pytest.raises(IsADirectoryError):
        load_tile(replace(archive_tile, archive=tmp_path), max_memory_bytes=BUDGET)
    for archive in ("/vsicurl/http://tiles.test/pkg.dtccpkg", "{pkg}.dtccpkg"):
        with pytest.raises(ValueError, match="local file path"):
            load_tile(replace(archive_tile, archive=archive), max_memory_bytes=BUDGET)


@pytest.mark.parametrize(
    "path", ["/vsizip/pkg.zip/tile.tif", "\\vsicurl\\http://tiles.test/tile.tif"]
)
def test_local_tile_path_must_not_be_a_gdal_virtual_path(tmp_path, path):
    tile = replace(make_tile(tmp_path), path=path)
    with pytest.raises(ValueError, match="GDAL virtual file path"):
        load_tile(tile, max_memory_bytes=BUDGET)


def test_unbound_tile_is_refused_by_load_tile_itself(tmp_path, monkeypatch):
    save_model(cached_collection(tmp_path / "cache"), tmp_path / "tiles.dtcc")
    tile = load_model(tmp_path / "tiles.dtcc")[0]
    monkeypatch.chdir(tmp_path)
    write_tiff(Path(REFERENCES[0]), band_values(4))
    with pytest.raises(ValueError, match="no package binding"):
        load_tile(tile, max_memory_bytes=BUDGET)


def test_relative_local_path_loads_from_its_directory(tmp_path, monkeypatch):
    write_tiff(tmp_path / "tile.tif", band_values(3))
    tile = replace(
        make_tile(tmp_path),
        path=Path("tile.tif"),
        extent=Bounds(*CELL),
        spektraltyp="rgb",
    )
    monkeypatch.chdir(tmp_path)
    assert tile.load().data.shape == (4, 4, 3)


@pytest.mark.skipif(sys.platform == "win32", reason="a colon is not a file name")
def test_local_path_that_looks_like_a_url_stays_local(
    tmp_path, monkeypatch, no_network
):
    write_tiff(tmp_path / "https:" / "tiles.test" / "tile.tif", band_values(3))
    tile = replace(
        make_tile(tmp_path),
        path=Path("https:/tiles.test/tile.tif"),
        extent=Bounds(*CELL),
        spektraltyp="rgb",
    )
    monkeypatch.chdir(tmp_path)
    assert tile.load().data.shape == (4, 4, 3)


@pytest.mark.skipif(sys.platform != "win32", reason="native Windows path spelling")
def test_native_windows_path_with_backslashes_loads(tmp_path):
    path = write_tiff(tmp_path / "tile.tif", band_values(3))
    native = str(PureWindowsPath(path.resolve()))
    assert "\\" in native and PureWindowsPath(native).is_absolute()
    tile = replace(
        make_tile(tmp_path), path=native, extent=Bounds(*CELL), spektraltyp="rgb"
    )
    assert tile.load().data.shape == (4, 4, 3)


# Rejections


@pytest.fixture
def package(tmp_path):
    cached_collection(tmp_path / "cache").export(tmp_path / "pkg", canonical=True)
    return tmp_path / "pkg"


def edit_manifest(package, edit):
    manifest_path = package / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    edit(manifest["artifacts"])
    manifest_path.write_text(json.dumps(manifest))


def test_tile_with_a_mismatched_sha256_is_refused(package):
    member = package / REFERENCES[0]
    data = bytearray(member.read_bytes())
    data[-1] ^= 0xFF
    member.write_bytes(bytes(data))
    with pytest.raises(ValueError, match="size/sha256 does not match"):
        load_model_package(package)


def test_reference_without_a_manifest_tile_is_refused(package):
    edit_manifest(package, lambda artifacts: artifacts.pop())
    with pytest.raises(ValueError, match="no manifest artifact with role 'tile'"):
        load_model_package(package)


def test_manifest_tile_without_a_reference_is_refused(package):
    extra = "artifacts/tiles/orto-o2-2025/9-extra.tif"
    shutil.copyfile(package / REFERENCES[1], package / extra)

    def add(artifacts):
        artifacts.append(dict(artifacts[-1], path=extra))

    edit_manifest(package, add)
    with pytest.raises(ValueError, match="not referenced by the canonical model"):
        load_model_package(package)


def test_reference_to_a_derived_artifact_is_refused(package):
    def derive(artifacts):
        artifacts[-1].update(role="derived", derived_from="artifacts/model.dtcc")

    edit_manifest(package, derive)
    with pytest.raises(ValueError, match="no manifest artifact with role 'tile'"):
        load_model_package(package)


@pytest.mark.parametrize("form", FORMS)
def test_format_is_refused_before_writing(tmp_path, form):
    collection = cached_collection(tmp_path / "cache")
    out = tmp_path / "out"
    out.mkdir()
    with pytest.raises(ValueError, match="format"):
        collection.export(out / form, canonical=True, format="tif")
    with pytest.raises(ValueError, match="format"):
        export_model_package(collection, out / form, canonical=True, format="png")
    assert list(out.iterdir()) == []


def tile_collection(tmp_path, empty):
    if not empty:
        return cached_collection(tmp_path / "cache")
    collection = RasterTileCollection()
    collection.dataset_context = make_context()
    return collection


class Forbidden:
    """Raises on every attribute access, call and truth test."""

    def __getattr__(self, name):
        raise AssertionError("nothing may be uploaded")

    def __call__(self, *args, **kwargs):
        raise AssertionError("nothing may be uploaded")

    def __bool__(self):
        raise AssertionError("nothing may be uploaded")


@pytest.mark.parametrize("empty", [False, True])
@pytest.mark.parametrize("form", FORMS)
def test_publish_is_refused_before_configuration_is_read(
    tmp_path, monkeypatch, form, empty
):
    def forbidden(*args, **kwargs):
        raise AssertionError("upload configuration must not be read")

    monkeypatch.setattr(DatasetUploadClient, "from_config", forbidden)
    collection = tile_collection(tmp_path, empty)
    with pytest.raises(NotImplementedError, match="RasterTileCollection"):
        collection.publish(dataset_key="orthophoto")
    package = collection.export(tmp_path / form, canonical=True)
    with pytest.raises(NotImplementedError, match="RasterTileCollection"):
        package.publish(dataset_key="orthophoto")


@pytest.mark.parametrize("empty", [False, True])
@pytest.mark.parametrize("form", FORMS)
def test_publish_never_touches_a_custom_uploader(tmp_path, form, empty):
    package = tile_collection(tmp_path, empty).export(tmp_path / form, canonical=True)
    with pytest.raises(NotImplementedError, match="RasterTileCollection"):
        package.publish(dataset_key="orthophoto", uploader=Forbidden())


@pytest.mark.parametrize("empty", [False, True])
@pytest.mark.parametrize("form", FORMS)
def test_upload_admission_refuses_a_tile_package(tmp_path, no_network, form, empty):
    target = tmp_path / form
    package = tile_collection(tmp_path, empty).export(target, canonical=True)
    root = target
    if form == "pkg.dtccpkg":
        root = tmp_path / "extracted"
        with zipfile.ZipFile(target) as archive:
            archive.extractall(root)
    client = DatasetUploadClient("http://upload.test", "token", session=Forbidden())
    with pytest.raises(DatasetPackageError, match="RasterTileCollection"):
        client.upload_package(
            dataset_key="orthophoto",
            manifest_path=root / "manifest.json",
            files=[root / artifact.path for artifact in package.artifacts],
        )


@pytest.mark.parametrize("intersects", [True, False])
@pytest.mark.parametrize("binding", ["archive", "unbound"])
def test_mosaic_refuses_package_tiles_before_opening_files(
    tmp_path, monkeypatch, binding, intersects
):
    original = cached_collection(tmp_path / "cache")
    if binding == "archive":
        original.export(tmp_path / "pkg.dtccpkg", canonical=True)
        other, reason = load_model_package(tmp_path / "pkg.dtccpkg")[1], "tile.load()"
    else:
        save_model(original, tmp_path / "tiles.dtcc")
        other, reason = load_model(tmp_path / "tiles.dtcc")[1], "no package binding"
        monkeypatch.chdir(tmp_path)
        write_tiff(Path(REFERENCES[1]), band_values(3))
    if not intersects:
        other = replace(other, extent=Bounds(0.0, 0.0, 1.0, 1.0))

    def forbidden(*args, **kwargs):
        raise AssertionError("no file may be opened")

    monkeypatch.setattr(rasterio, "open", forbidden)
    with pytest.raises(ValueError, match=re.escape(reason)):
        build_mosaic(
            [original[1], other], CELL, resolution=None, max_memory_bytes=BUDGET
        )


@pytest.mark.parametrize(
    "field",
    [
        {"derived_from": "artifacts/model.dtcc"},
        {"model_type": "City"},
        {"model_schema_version": 1},
        {"format": "png"},
        {"crs": "EPSG:3006"},
    ],
)
def test_tile_entry_carries_no_other_metadata(package, field):
    edit_manifest(package, lambda artifacts: artifacts[-1].update(field))
    with pytest.raises(ValueError, match="Tile artifact"):
        load_model_package(package)


def test_tile_entry_is_checked_before_any_artifact_is_read(package):
    edit_manifest(package, lambda artifacts: artifacts[-1].update(crs="EPSG:3006"))
    model = package / "artifacts" / "model.dtcc"
    data = bytearray(model.read_bytes())
    data[-1] ^= 0xFF
    model.write_bytes(bytes(data))
    with pytest.raises(ValueError, match="Tile artifact"):
        load_model_package(package)


@pytest.mark.skipif(sys.platform == "win32", reason="a colon is not a file name")
def test_mosaic_of_a_local_path_that_looks_like_a_url_stays_local(
    tmp_path, monkeypatch, no_network
):
    write_tiff(tmp_path / "https:" / "tiles.test" / "tile.tif", band_values(3))
    tile = replace(
        make_tile(tmp_path),
        path=Path("https:/tiles.test/tile.tif"),
        extent=Bounds(*CELL),
        spektraltyp="rgb",
    )
    monkeypatch.chdir(tmp_path)
    mosaic = build_mosaic([tile], CELL, resolution=None, max_memory_bytes=BUDGET)
    assert mosaic.failures == ()
    assert mosaic.valid_pixels == 16


def test_mosaic_of_a_generator_matches_the_list(tmp_path):
    tiles = list(cached_collection(tmp_path / "cache"))
    expected = build_mosaic(tiles, CELL, resolution=None, max_memory_bytes=BUDGET)
    mosaic = build_mosaic(
        (tile for tile in tiles), CELL, resolution=None, max_memory_bytes=BUDGET
    )
    assert expected.valid_pixels == 16
    assert mosaic.valid_pixels == expected.valid_pixels
    np.testing.assert_array_equal(mosaic.raster.data, expected.raster.data)


def test_archive_path_with_braces_is_refused_at_binding(tmp_path):
    target = tmp_path / "a{b}" / "pkg.dtccpkg"
    cached_collection(tmp_path / "cache").export(target, canonical=True)
    with pytest.raises(ValueError, match="without braces"):
        load_model_package(target)


@pytest.mark.parametrize("form", FORMS)
def test_failed_tile_copy_leaves_no_destination(tmp_path, form):
    collection = cached_collection(tmp_path / "cache")
    collection[1].path.unlink()
    out = tmp_path / "out"
    out.mkdir()
    with pytest.raises(FileNotFoundError):
        collection.export(out / form, canonical=True)
    assert list(out.iterdir()) == []


@pytest.mark.parametrize("form", FORMS)
def test_unbound_tiles_are_not_exported(tmp_path, monkeypatch, form):
    save_model(cached_collection(tmp_path / "cache"), tmp_path / "tiles.dtcc")
    restored = load_model(tmp_path / "tiles.dtcc")
    restored.dataset_context = make_context()
    monkeypatch.chdir(tmp_path)
    write_tiff(Path(REFERENCES[0]), band_values(4))
    out = tmp_path / "out"
    out.mkdir()
    with pytest.raises(ValueError, match="no package binding"):
        restored.export(out / form, canonical=True)
    assert list(out.iterdir()) == []
