from dtcc_core.builder.evaluation.dataset import GroundTruth
from dtcc_core.builder.evaluation.metrics import (
    detected_plane_count,
    expected_plane_count,
    plane_count_delta,
)


def test_detected_plane_count_from_diagnostics():
    rec = {"planes": [object(), object(), object()]}
    assert detected_plane_count(rec) == 3


def test_detected_plane_count_none_planes():
    assert detected_plane_count({"planes": None}) == 0
    assert detected_plane_count({}) == 0


def test_expected_plane_count_explicit_overrides_type_default():
    gt = GroundTruth(roof_type="GABLED", expected_plane_count=5)
    assert expected_plane_count(gt) == 5


def test_expected_plane_count_defaults_per_type():
    assert expected_plane_count(GroundTruth(roof_type="FLAT")) == 1
    assert expected_plane_count(GroundTruth(roof_type="GABLED")) == 2
    assert expected_plane_count(GroundTruth(roof_type="HIPPED")) == 4


def test_expected_plane_count_unknown_type_returns_none():
    assert expected_plane_count(GroundTruth(roof_type="MANSARD")) is None


def test_plane_count_delta_absolute():
    rec = {"planes": [1, 2, 3, 4]}
    assert plane_count_delta(rec, GroundTruth(roof_type="GABLED")) == 2


def test_plane_count_delta_none_when_expected_unknown():
    rec = {"planes": [1, 2]}
    assert plane_count_delta(rec, GroundTruth(roof_type="MANSARD")) is None
