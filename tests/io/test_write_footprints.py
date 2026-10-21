import pytest
from pathlib import Path
import json
from shapely.geometry import shape

import dtcc_core


@pytest.fixture
def data_dir():
    return Path(__file__).parent / ".." / "data"


@pytest.fixture
def building_shp_path(data_dir):
    return str((data_dir / "MinimalCase" / "PropertyMap.shp").resolve())


@pytest.fixture
def loaded_buildings(building_shp_path):
    return dtcc_core.io.load_footprints(building_shp_path)


@pytest.fixture
def simple_city(loaded_buildings):
    city = dtcc_core.model.City()
    city.add_buildings(loaded_buildings)
    return city


def test_write_geojson(simple_city, loaded_buildings):
    dtcc_core.io.save_footprints(simple_city, "footprints.geojson")
    assert Path("footprints.geojson").exists()
    with open("footprints.geojson") as f:
        data = json.load(f)
    assert data["type"] == "FeatureCollection"
    assert len(data["features"]) == 5

    Path("footprints.geojson").unlink()


def test_write_shp_zip(simple_city, loaded_buildings):
    dtcc_core.io.save_footprints(simple_city, "footprints.shp.zip")
    assert Path("footprints.shp.zip").exists()
    with open("footprints.shp.zip", "rb") as f:
        data = f.read()
        assert data[:2] == b"PK"  # zip file signature
    Path("footprints.shp.zip").unlink()


def test_source_bounds_vector_export_and_reload(simple_city, tmp_path):
    source_bounds = simple_city.buildings[0].attributes["source_bounds"]
    path = tmp_path / "footprints.geojson"
    dtcc_core.io.save_footprints(simple_city, path)
    data = json.loads(path.read_text())
    feature = data["features"][0]
    assert feature["properties"]["source_bounds"] == ",".join(map(str, source_bounds))

    restored = dtcc_core.io.load_footprints(path)
    assert restored[0].attributes["source_bounds"] == pytest.approx(
        shape(feature["geometry"]).bounds
    )


@pytest.mark.parametrize("extension", ["pb", "pb2"])
def test_footprint_protobuf_save_is_explicitly_unsupported(tmp_path, extension):
    path = tmp_path / f"footprints.{extension}"
    with pytest.raises(RuntimeError, match="format .pb.*not supported"):
        dtcc_core.io.save_footprints(dtcc_core.model.City(), path)
    assert not path.exists()


@pytest.mark.parametrize("extension", ["pb", "pb2"])
def test_footprint_protobuf_load_is_explicitly_unsupported(tmp_path, extension):
    path = tmp_path / f"footprints.{extension}"
    path.write_bytes(dtcc_core.model.City().to_proto().SerializeToString())
    with pytest.raises(RuntimeError, match="format .pb.*not supported"):
        dtcc_core.io.load_footprints(path)
