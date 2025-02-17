import pytest
import numpy as np
from pathlib import Path
from dtcc_core.model.object import RoadNetwork, RoadType
from dtcc_core.model.logging import info
from dtcc_core.io import roadnetwork
from dtcc_core.model import Bounds

# Try importing geopandas
try:
    import geopandas as gpd
except ImportError:
    gpd = None


@pytest.fixture
def data_path():
    return Path(__file__).parent / ".." / "data" / "road_network" / "test_road.shp"


@pytest.fixture
def test_bounds():
    return Bounds(xmin=58150, ymin=6398226, xmax=59150, ymax=6398984)


@pytest.fixture
def basic_roadnetwork(data_path):
    return roadnetwork.load(data_path)


def test_load_roadnetwork(basic_roadnetwork):
    assert isinstance(basic_roadnetwork, RoadNetwork)
    assert len(basic_roadnetwork.vertices) == 185
    assert len(basic_roadnetwork.edges) == 244
    assert len(basic_roadnetwork.length) == 244
    assert len(basic_roadnetwork.linestrings) == len(basic_roadnetwork.edges)


def test_load_roadnetwork_with_bounds(data_path, test_bounds):
    rn = roadnetwork.load(data_path, bounds=test_bounds)
    assert isinstance(rn, RoadNetwork)
    assert len(rn.edges) == 37

    for k, v in rn.attributes.items():
        assert len(v) == len(rn.edges)


def test_load_roadnetwork_no_geometry(data_path):
    rn = roadnetwork.load(data_path, load_geometry=False)
    assert isinstance(rn, RoadNetwork)
    assert len(rn.vertices) == 185
    assert len(rn.edges) == 244
    assert len(rn.length) == 244
    assert len(rn.linestrings) == 0


def test_roadnetwork_attributes(basic_roadnetwork):
    assert len(basic_roadnetwork.attributes) == 7
    for k, v in basic_roadnetwork.attributes.items():
        assert len(v) == 244


def test_roadnetwork_protobuf_conversion(basic_roadnetwork):
    pb = basic_roadnetwork.to_proto()
    rn2 = RoadNetwork()
    rn2.from_proto(pb)

    assert len(basic_roadnetwork.vertices) == len(rn2.vertices)
    assert np.allclose(basic_roadnetwork.vertices, rn2.vertices)
    assert len(basic_roadnetwork.edges) == len(rn2.edges)
    assert np.allclose(basic_roadnetwork.edges, rn2.edges)
    assert len(basic_roadnetwork.length) == len(rn2.length)
    assert np.allclose(basic_roadnetwork.length, rn2.length)


@pytest.mark.skipif(gpd is None, reason="Geopandas not installed")
def test_roadnetwork_to_dataframe(basic_roadnetwork):
    df = basic_roadnetwork.to_df()
    assert isinstance(df, gpd.GeoDataFrame)
    assert len(df) == len(basic_roadnetwork.edges)


if __name__ == "__main__":
    pytest.main()
