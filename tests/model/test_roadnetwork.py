import re
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
    rn.add_geometry(multilinestring, GeometryType.MULTILINESTRING)

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

    assert summary == "<RoadNetwork(num_vertices=3, num_edges=2, num_segments=2)>"
    assert "RoadNetwork" in details
    assert re.search(r"│Vertices\s*│3\s*│", details)
    assert re.search(r"│Edges\s*│2\s*│", details)
    assert re.search(r"│CRS\s*│EPSG:3006\s*│", details)
    assert re.search(r"│Total\s*│3.00\s*│", details)
    assert re.search(r"│residential\s*│1\s*│", details)
    assert re.search(r"│primary\s*│1\s*│", details)
    assert re.search(r"│One-way segments\s*│1\s*│", details)

    result = rn.info()
    captured = capsys.readouterr()
    assert result is None
    assert "RoadNetwork" in captured.out
    assert captured.out == details + "\n"


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
