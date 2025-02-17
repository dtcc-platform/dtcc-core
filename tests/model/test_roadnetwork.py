import pytest
import numpy as np
from dtcc_core.model.object import RoadNetwork, RoadType


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


if __name__ == "__main__":
    pytest.main()
