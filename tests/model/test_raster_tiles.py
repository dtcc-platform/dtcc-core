from dataclasses import replace
from datetime import UTC, datetime
from pathlib import Path

import pytest
import rasterio

from dtcc_core.datasets.schema import (
    DatasetContext,
    DatasetIdentity,
    DatasetMetadata,
    DatasetPresentation,
    DatasetProvenance,
    DatasetRequest,
)
from dtcc_core.io.model import load_model, save_model
from dtcc_core.model import Bounds, City, Object, exchange
from dtcc_core.model import dtcc_pb2 as wire
from dtcc_core.model.values.raster_tiles import RasterTile, RasterTileCollection
from tests.io.test_raster_tiles import CELL, band_values, write_tiff


def make_tile(tmp_path, name="o65775_6725_25_mr25", collection="orto-o2-2025"):
    return RasterTile(
        path=tmp_path / "cache" / collection / f"{name}.tif",
        id=name,
        collection=collection,
        datetime=datetime(2025, 5, 31, 8, 10, 56, tzinfo=UTC),
        extent=Bounds(672500.0, 6577500.0, 675000.0, 6580000.0),
        crs="EPSG:3006",
        spektraltyp="rgbi",
        resolution=0.16,
        size_bytes=659037310,
    )


def make_context():
    return DatasetContext(
        identity=DatasetIdentity(name="orthophoto", title="Orthophoto"),
        metadata=DatasetMetadata(crs=["EPSG:3006"]),
        provenance=DatasetProvenance(sources=["LM"]),
        presentation=DatasetPresentation(),
        request=DatasetRequest(dataset_name="orthophoto", bounds=[0.0, 0.0, 1.0, 1.0]),
    )


@pytest.fixture
def no_raster_reads(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("tile records must not open raster files")

    monkeypatch.setattr(rasterio, "open", forbidden)


@pytest.fixture
def collection(tmp_path, no_raster_reads):
    tiles = [
        make_tile(tmp_path),
        make_tile(tmp_path, "o65775_6725_25_im17", "orto-o2-2017"),
    ]
    return RasterTileCollection(tiles=tiles)


def test_tile_keeps_metadata_without_touching_its_file(tmp_path, no_raster_reads):
    tile = make_tile(tmp_path)
    assert not tile.path.exists()
    assert tile.path == tmp_path / "cache" / "orto-o2-2025" / "o65775_6725_25_mr25.tif"
    assert tile.id == "o65775_6725_25_mr25"
    assert tile.collection == "orto-o2-2025"
    assert tile.datetime == datetime(2025, 5, 31, 8, 10, 56, tzinfo=UTC)
    assert tile.extent.tuple == (672500.0, 6577500.0, 675000.0, 6580000.0)
    assert tile.crs == "EPSG:3006"
    assert tile.spektraltyp == "rgbi"
    assert tile.resolution == 0.16
    assert tile.size_bytes == 659037310
    repr(tile)


def test_tile_records_are_immutable(tmp_path):
    tile = make_tile(tmp_path)
    with pytest.raises(AttributeError):
        tile.path = Path("elsewhere.tif")


def test_collection_is_a_sequence_of_tiles(collection):
    assert len(collection) == 2
    assert collection[0].id == "o65775_6725_25_mr25"
    assert collection[-1].collection == "orto-o2-2017"
    assert [tile.id for tile in collection] == [
        "o65775_6725_25_mr25",
        "o65775_6725_25_im17",
    ]
    assert collection[1:][0] is collection.tiles[1]


def test_empty_collection(no_raster_reads):
    collection = RasterTileCollection()
    assert len(collection) == 0
    assert list(collection) == []


def test_collection_summary_and_info_read_no_files(collection):
    assert "2" in repr(collection)
    collection.info(print=False)


def test_collection_carries_dataset_context(collection):
    context = make_context()
    collection.dataset_context = context
    assert collection.dataset_context is context
    assert collection.metadata.crs == ["EPSG:3006"]


@pytest.mark.parametrize(
    "kwargs",
    [{}, {"format": "json"}, {"format": "tif"}, {"canonical": True}],
)
def test_export_is_refused_before_writing_anything(collection, tmp_path, kwargs):
    collection.dataset_context = make_context()
    target = tmp_path / "package"
    with pytest.raises(NotImplementedError, match="RasterTileCollection"):
        collection.export(target, **kwargs)
    assert not target.exists()


@pytest.mark.parametrize("with_context", [True, False])
@pytest.mark.parametrize("with_uploader", [True, False])
def test_publish_is_refused_before_configuration_or_upload(
    collection, monkeypatch, with_context, with_uploader
):
    from dtcc_core.datasets.publish import DatasetUploadClient

    def forbidden(*args, **kwargs):
        raise AssertionError("upload configuration must not be read")

    monkeypatch.setattr(DatasetUploadClient, "from_config", forbidden)

    class Uploader:
        def upload(self, *args, **kwargs):
            raise AssertionError("nothing may be uploaded")

    if with_context:
        collection.dataset_context = make_context()
    kwargs = {"uploader": Uploader()} if with_uploader else {}
    with pytest.raises(NotImplementedError, match="RasterTileCollection"):
        collection.publish(dataset_key="orthophoto", **kwargs)


METADATA = (
    "id",
    "collection",
    "datetime",
    "extent",
    "crs",
    "spektraltyp",
    "resolution",
    "size_bytes",
)
REFERENCES = [
    "artifacts/tiles/orto-o2-2025/0-o65775_6725_25_mr25.tif",
    "artifacts/tiles/orto-o2-2017/1-o65775_6725_25_mr25.tif",
]


def shared_id_collection(base):
    """Two tiles from different collections that share one tile id."""
    return RasterTileCollection(
        tiles=[
            make_tile(base),
            replace(
                make_tile(base, collection="orto-o2-2017"),
                resolution=None,
                size_bytes=None,
            ),
        ]
    )


def metadata(collection):
    return [tuple(getattr(tile, name) for name in METADATA) for tile in collection]


def round_trip(source, route, tmp_path):
    if route == "bytes":
        return exchange.loads(exchange.dumps(source))
    if route == "proto":
        restored = RasterTileCollection()
        restored.from_proto(source.to_proto().SerializeToString())
        return restored
    save_model(source, tmp_path / "tiles.dtcc")
    return load_model(tmp_path / "tiles.dtcc", expected_type=RasterTileCollection)


@pytest.mark.parametrize("route", ["bytes", "proto", "dtcc"])
def test_collection_round_trips_metadata_and_package_references(
    tmp_path, no_raster_reads, route
):
    source = shared_id_collection(tmp_path / "cache")
    restored = round_trip(source, route, tmp_path)
    assert type(restored) is RasterTileCollection
    assert metadata(restored) == metadata(source)
    assert [tile.path for tile in restored] == [Path(ref) for ref in REFERENCES]
    assert all(tile.archive is None and tile._unbound for tile in restored)
    assert exchange.dumps(restored) == exchange.dumps(source)


def test_local_cache_paths_never_reach_the_wire(tmp_path, no_raster_reads):
    source = shared_id_collection(tmp_path / "cache")
    assert str(tmp_path) not in source.to_json()
    assert str(tmp_path).encode() not in exchange.dumps(source)


@pytest.mark.parametrize("route", ["bytes", "proto", "dtcc"])
def test_empty_collection_round_trips(tmp_path, route):
    restored = round_trip(RasterTileCollection(), route, tmp_path)
    assert type(restored) is RasterTileCollection
    assert len(restored) == 0
    assert exchange.dumps(restored) == exchange.dumps(RasterTileCollection())


def first_artifact(value):
    def mutate(records):
        records[0]["artifact"].text = value

    return mutate


def first_tile_in(collection, artifact):
    def mutate(records):
        records[0]["collection"].text = collection
        records[0]["artifact"].text = artifact

    return mutate


def duplicate_artifact(records):
    records[1]["artifact"].CopyFrom(records[0]["artifact"])


def missing_artifact(records):
    del records[0]["artifact"]


TILE = "0-o65775_6725_25_mr25.tif"
INVALID_RECORDS = {
    "absolute": (
        first_artifact(f"/artifacts/tiles/orto-o2-2025/{TILE}"),
        "not a portable package path",
    ),
    "traversal": (
        first_artifact(f"artifacts/tiles/orto-o2-2025/../orto-o2-2025/{TILE}"),
        "not a portable package path",
    ),
    "backslash": (
        first_artifact(f"artifacts\\tiles\\orto-o2-2025\\{TILE}"),
        "not a portable package path",
    ),
    "reserved": (
        first_tile_in("CON", f"artifacts/tiles/CON/{TILE}"),
        "not a portable package path",
    ),
    "outside": (
        first_artifact(f"artifacts/other/orto-o2-2025/{TILE}"),
        "must be 'artifacts/tiles/",
    ),
    "slash_in_collection": (
        first_tile_in("orto/2025", f"artifacts/tiles/orto/2025/{TILE}"),
        "must be 'artifacts/tiles/orto2025/",
    ),
    "duplicate": (duplicate_artifact, "must be 'artifacts/tiles/orto-o2-2017/1-"),
    "missing": (missing_artifact, "exactly the keys"),
}


@pytest.mark.parametrize("failure", INVALID_RECORDS)
def test_bytes_boundary_rejects_invalid_tile_records(failure):
    mutate, reason = INVALID_RECORDS[failure]
    pb = wire.ModelFile.FromString(exchange.dumps(shared_id_collection(Path("cache"))))
    mutate([value.object.values for value in pb.object.attributes["tiles"].list.values])
    with pytest.raises(ValueError, match=reason):
        exchange.loads(pb.SerializeToString())


def test_bytes_boundary_rejects_the_kind_under_a_city():
    city = City(id="city")
    city.add_child(Object(id="child"))
    pb = wire.ModelFile.FromString(exchange.dumps(city))
    pb.object.children[0].kind = wire.Object.RASTER_TILE_COLLECTION
    with pytest.raises(NotImplementedError, match="kind 19"):
        exchange.loads(pb.SerializeToString())


def test_collection_under_a_city_is_refused(tmp_path):
    city = City(id="city")
    city.children[RasterTileCollection] = [shared_id_collection(tmp_path)]
    with pytest.raises(NotImplementedError, match="RasterTileCollection"):
        exchange.dumps(city)


def test_restored_tile_never_loads_from_the_working_directory(tmp_path, monkeypatch):
    save_model(shared_id_collection(tmp_path / "cache"), tmp_path / "tiles.dtcc")
    restored = load_model(tmp_path / "tiles.dtcc")
    monkeypatch.chdir(tmp_path)
    write_tiff(Path(REFERENCES[0]), band_values(3))
    tile = restored[0]
    assert tile.path.is_file()
    with pytest.raises(ValueError, match="no package binding"):
        tile.load()
    local = make_tile(Path("unused"))
    local = replace(local, path=tile.path, extent=Bounds(*CELL), spektraltyp="rgb")
    assert not local._unbound
    assert local.load().data.shape == (4, 4, 3)


def test_save_is_refused_without_writing(collection, tmp_path):
    target = tmp_path / "tiles.json"
    with pytest.raises(AttributeError, match="RasterTileCollection"):
        collection.save(target)
    assert not target.exists()


@pytest.mark.parametrize("name", ["pkg", "pkg.dtccpkg"])
def test_canonical_package_export_is_refused_before_writing(collection, tmp_path, name):
    from dtcc_core.datasets.package import export_model_package

    collection.dataset_context = make_context()
    out = tmp_path / "out"
    out.mkdir()
    with pytest.raises(NotImplementedError, match="RasterTileCollection"):
        export_model_package(collection, out / name, canonical=True)
    assert list(out.iterdir()) == []


def test_load_delegates_to_load_tile_with_the_budget(tmp_path, monkeypatch):
    from dtcc_core.io import raster_tiles

    calls = []

    def recording_load_tile(tile, **kwargs):
        calls.append((tile, kwargs))
        return "raster"

    monkeypatch.setattr(raster_tiles, "load_tile", recording_load_tile)
    tile = make_tile(tmp_path)
    assert tile.load(max_memory_bytes=123) == "raster"
    assert tile.load() == "raster"
    assert calls == [
        (tile, {"max_memory_bytes": 123}),
        (tile, {"max_memory_bytes": raster_tiles.DEFAULT_MAX_MEMORY_BYTES}),
    ]
