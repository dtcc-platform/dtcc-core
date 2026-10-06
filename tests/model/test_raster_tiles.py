from datetime import datetime, timezone
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
from dtcc_core.model import Bounds
from dtcc_core.model.values.raster_tiles import RasterTile, RasterTileCollection


def make_tile(tmp_path, name="o65775_6725_25_mr25", collection="orto-o2-2025"):
    return RasterTile(
        path=tmp_path / "cache" / collection / f"{name}.tif",
        id=name,
        collection=collection,
        datetime=datetime(2025, 5, 31, 8, 10, 56, tzinfo=timezone.utc),
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
    assert tile.datetime == datetime(2025, 5, 31, 8, 10, 56, tzinfo=timezone.utc)
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


@pytest.mark.parametrize("method", ["to_proto", "to_json"])
def test_canonical_serialization_is_refused(collection, method):
    with pytest.raises(NotImplementedError, match="RasterTileCollection"):
        getattr(collection, method)()


def test_save_is_refused_without_writing(collection, tmp_path):
    target = tmp_path / "tiles.json"
    with pytest.raises(AttributeError, match="RasterTileCollection"):
        collection.save(target)
    assert not target.exists()
