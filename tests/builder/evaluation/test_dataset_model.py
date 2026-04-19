import numpy as np

from dtcc_core.model.geometry.pointcloud import PointCloud
from dtcc_core.model.geometry.surface import Surface
from dtcc_core.builder.evaluation.dataset import EvalSample, GroundTruth


def test_ground_truth_supports_minimum_fields():
    gt = GroundTruth(roof_type="FLAT")
    assert gt.roof_type == "FLAT"
    assert gt.ridge_height is None
    assert gt.eave_height is None
    assert gt.geometry is None


def test_ground_truth_supports_optional_heights():
    gt = GroundTruth(roof_type="GABLED", ridge_height=8.0, eave_height=5.0)
    assert gt.ridge_height == 8.0
    assert gt.eave_height == 5.0


def test_eval_sample_bundles_inputs_and_truth():
    fp = Surface(vertices=np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], dtype=float))
    pc = PointCloud(points=np.zeros((10, 3)))
    gt = GroundTruth(roof_type="FLAT")

    sample = EvalSample(
        building_id="b001",
        footprint=fp,
        point_cloud=pc,
        ground_truth=gt,
    )

    assert sample.building_id == "b001"
    assert sample.footprint is fp
    assert sample.point_cloud is pc
    assert sample.ground_truth.roof_type == "FLAT"
