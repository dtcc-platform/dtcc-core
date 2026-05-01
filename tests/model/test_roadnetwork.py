import pytest
import numpy as np
import dtcc_core
from dtcc_core.model.object import GeometryType, RoadNetwork, RoadType
from dtcc_core.model.geometry import LineString, MultiLineString


def test_create_roadnetwork():
    rn = RoadNetwork()
    assert len(rn.vertices) == 0
    assert len(rn.edges) == 0


def test_roadnetwork_bounds():
    rn = RoadNetwork()
    rn.vertices = np.array([(0, 0), (1, 0), (1, 1), (0, 2)])
    rn.edges = np.array(
        [
            (0, 1),
            (1, 2),
            (2, 3),
        ]
    )

    bounds = rn.bounds
    assert bounds.xmin == 0
    assert bounds.xmax == 1
    assert bounds.ymin == 0
    assert bounds.ymax == 2


def test_roadnetwork_to_arrays():
    rn = RoadNetwork()
    rn.vertices = np.array([(0, 0), (1, 0), (1, 1)])
    rn.edges = np.array([(0, 1), (1, 2)])
    rn.length = np.array([1.0, 2.0])
    rn.attributes = {"highway": ["residential", "primary"]}

    line = LineString()
    line.vertices = np.array([(0, 0), (1, 0), (1, 1)])
    multilinestring = MultiLineString()
    multilinestring.linestrings.append(line)
    rn.geometry[GeometryType.MULTILINESTRING] = multilinestring

    arrays = rn.to_arrays(include_geometry=True)

    assert arrays["vertices"].shape == (3, 2)
    assert arrays["edges"].shape == (2, 2)
    assert arrays["lengths"].shape == (2,)
    assert arrays["attributes"]["highway"].tolist() == ["residential", "primary"]
    assert arrays["line_vertices"].shape == (3, 2)
    assert arrays["line_offsets"].tolist() == [0, 3]


def test_roadnetwork_str_and_info(capsys):
    rn = RoadNetwork()
    rn.vertices = np.array([(0, 0), (1, 0), (1, 1)])
    rn.edges = np.array([(0, 1), (1, 2)])
    rn.length = np.array([1.0, 2.0])
    rn.attributes = {
        "highway": ["residential", "primary"],
        "oneway": [True, False],
    }
    rn.transform.srs = "EPSG:3006"

    summary = str(rn)
    details = rn.info(print=False)

    assert "DTCC RoadNetwork with 3 vertices, 2 edge(s)" in summary
    assert "DTCC RoadNetwork" in details
    assert "Vertices: 3" in details
    assert "Edges: 2" in details
    assert "CRS: EPSG:3006" in details
    assert "Total: 3.00" in details
    assert "residential: 1" in details
    assert "primary: 1" in details
    assert "One-way segments: 1" in details

    result = rn.info()
    captured = capsys.readouterr()
    assert result is None
    assert "DTCC RoadNetwork" in captured.out
    assert "Vertices: 3" in captured.out


def test_roadnetwork_to_matrix_without_self_loops():
    rn = RoadNetwork()
    rn.vertices = np.array([(0, 0), (1, 0), (1, 1)])
    rn.edges = np.array([(0, 1), (1, 2)])
    rn.length = np.array([1.0, 2.0])

    matrix = rn.to_matrix()

    assert matrix.nnz == 4
    assert matrix[0, 1] == 1.0
    assert matrix[1, 0] == 1.0
    assert matrix[1, 2] == 2.0
    assert matrix[2, 1] == 2.0


def test_roadnetwork_builder_methods_registered():
    assert dtcc_core.builder is not None
    assert hasattr(RoadNetwork, "to_matrix")
    assert hasattr(RoadNetwork, "to_surfaces")


if __name__ == "__main__":
    pytest.main()
