import pytest
from pathlib import Path
from dtcc_core import builder, io
from dtcc_core.model import RoadNetwork, GeometryType, Surface


@pytest.fixture
def test_data():
    data_path = (
        Path(__file__).parent / ".." / "data" / "road_network" / "test_road.shp"
    ).resolve()
    assert data_path.is_file()
    return data_path


@pytest.fixture
def road_network(test_data):
    return io.load_roadnetwork(test_data)


@pytest.fixture
def network_stats(road_network):
    """Common statistics used across tests."""
    edges = road_network.edges
    return {
        "largest_idx": max([e[0] for e in edges] + [e[1] for e in edges]),
        "num_roads": len(edges),
    }


def test_bidirectional_to_matrix(road_network, network_stats):
    """Test conversion to bidirectional matrix representation."""
    rn_matrix = road_network.to_matrix()

    # Check matrix dimensions
    assert rn_matrix.shape == (
        network_stats["largest_idx"] + 1,
        network_stats["largest_idx"] + 1,
    )

    # Check number of non-zero elements
    # bidirectional so each edge is counted twice, minus 1 because of the loop
    expected_nnz = (network_stats["num_roads"] * 2) - 1
    assert rn_matrix.nnz == expected_nnz

    # Check specific edge properties
    edge_12 = road_network.edges[12]
    length_12 = road_network.length[12]
    assert rn_matrix[edge_12[0], edge_12[1]] == length_12
    assert rn_matrix[edge_12[1], edge_12[0]] == length_12


def test_unidirectional_to_matrix(road_network, network_stats):
    """Test conversion to unidirectional matrix representation."""
    rn_matrix = road_network.to_matrix(bidirectional=False)
    assert rn_matrix.nnz == network_stats["num_roads"]


def test_to_polygon(road_network):
    """Test conversion to polygon representation."""
    road_width = 4
    rn_poly = road_network.to_surfaces(as_shapely=True, widths=road_width)

    assert len(rn_poly) == len(road_network.linestrings)
    assert rn_poly[0].geom_type == "Polygon"

    # Check area approximation
    expected_area = road_network.length[0] * road_width
    actual_area = rn_poly[0].area
    area_difference = abs(1 - (actual_area / expected_area))
    assert area_difference <= 0.02


def test_to_surfaces(road_network):
    """Test conversion to surface objects."""
    rn_surfaces = road_network.to_surfaces(widths=4, as_shapely=False)

    assert len(rn_surfaces) == len(road_network.linestrings)
    assert isinstance(rn_surfaces[0], Surface)
