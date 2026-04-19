import pytest

from dtcc_core.builder.evaluation.dataset import GroundTruth
from dtcc_core.builder.evaluation.metrics import ridge_height_error, eave_height_error


class _StubResult:
    def __init__(self, ridge_height=None, eave_height=None):
        self.ridge_height = ridge_height
        self.eave_height = eave_height


def test_ridge_height_error_abs_difference():
    rec = {"classification": _StubResult(ridge_height=8.0)}
    gt = GroundTruth(roof_type="GABLED", ridge_height=7.5)
    assert ridge_height_error(rec, gt) == 0.5


def test_ridge_height_error_no_gt_returns_none():
    rec = {"classification": _StubResult(ridge_height=8.0)}
    gt = GroundTruth(roof_type="FLAT")
    assert ridge_height_error(rec, gt) is None


def test_ridge_height_error_no_classification_returns_none():
    gt = GroundTruth(roof_type="GABLED", ridge_height=7.5)
    assert ridge_height_error({}, gt) is None


def test_eave_height_error_behaves_symmetrically():
    rec = {"classification": _StubResult(eave_height=5.1)}
    gt = GroundTruth(roof_type="GABLED", eave_height=5.0)
    assert eave_height_error(rec, gt) == pytest.approx(0.1, abs=1e-9)
