import pytest
import numpy as np
from dtcc_core.model import PointCloud, Bounds


@pytest.fixture
def pc():
    pc = PointCloud()
    pc.points = np.array([[0, 0, 0], [1, 1, 1], [2, 2, 2]])
    pc.calculate_bounds()
    pc.crs = "EPSG:3857"
    pc.origin = (0, 0)
    pc.classification = np.array([1, 2, 3])
    pc.intensity = np.array([1, 2, 3])
    pc.return_number = np.array([1, 1, 1])
    pc.num_returns = np.array([1, 1, 1])

    return pc


def test_empty():
    pc = PointCloud()
    assert pc.bounds.tuple == (0, 0, 0, 0)


def test_calc_bounds(pc):
    pc.calculate_bounds()
    assert pc.bounds.tuple == (0, 0, 2, 2)


def test_to_proto(pc):

    proto_pc = pc.to_proto()
    assert proto_pc.point_cloud.points == [0, 0, 0, 1, 1, 1, 2, 2, 2]
    assert proto_pc.point_cloud.classification == [1, 2, 3]


def test_to_protobuf_missing_fields():
    pc = PointCloud()
    pc.points = np.array([[0, 0, 0], [1, 1, 1], [2, 2, 2]])
    pc.calculate_bounds()
    proto_pc = pc.to_proto()
    assert proto_pc.point_cloud.points == [0, 0, 0, 1, 1, 1, 2, 2, 2]
    assert proto_pc.point_cloud.classification == []


def test_from_proto(pc):
    proto_pc = pc.to_proto()
    pc2 = PointCloud()
    pc2.from_proto(proto_pc)
    assert pc2.points.tolist() == pc.points.tolist()
    assert pc2.classification.tolist() == pc.classification.tolist()
    assert pc2.intensity.tolist() == pc.intensity.tolist()
    assert pc2.return_number.tolist() == pc.return_number.tolist()
    assert pc2.num_returns.tolist() == pc.num_returns.tolist()


if __name__ == "__main__":
    pytest.main()
