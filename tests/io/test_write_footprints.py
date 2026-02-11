import pytest
from pathlib import Path
import json

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
