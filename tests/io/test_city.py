import pytest
from pathlib import Path
from dtcc_core import io
from dtcc_core.model import City


@pytest.fixture
def test_data_dir():
    return Path(__file__).parent / ".." / "data"


@pytest.fixture
def cityjson_path(test_data_dir):
    return test_data_dir / "cityjson" / "DenHaag_01.city.json.zip"


@pytest.fixture
def mesh_path(test_data_dir):
    return test_data_dir / "meshes" / "9_cubes.obj"


@pytest.fixture
def cityjson_city(cityjson_path):
    return io.load_city(cityjson_path)


def test_load_cityjson(cityjson_city):
    assert isinstance(cityjson_city, City)
    assert len(cityjson_city.buildings) == 844


def test_load_city_mesh(mesh_path):
    city = io.load_city(mesh_path)
    assert isinstance(city, City)
    assert len(city.buildings) == 9


def test_buildings_to_df(cityjson_city):
    df = io.city.buildings_to_df(cityjson_city)
    assert len(df) == 844


def test_buildings_to_df_with_geometry(cityjson_city):
    df = io.city.buildings_to_df(cityjson_city)
    assert "geometry" in df.columns


if __name__ == "__main__":
    pytest.main()
