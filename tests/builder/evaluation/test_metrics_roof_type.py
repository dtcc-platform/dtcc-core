from dtcc_core.model.object.building import Building
from dtcc_core.builder.evaluation.dataset import GroundTruth
from dtcc_core.builder.evaluation.metrics import roof_type_prediction, roof_type_correct


def test_roof_type_prediction_reads_from_attributes():
    b = Building(id="b1")
    b.attributes["roof_type"] = "GABLED"
    assert roof_type_prediction(b) == "GABLED"


def test_roof_type_prediction_missing_returns_none():
    b = Building(id="b2")
    assert roof_type_prediction(b) is None


def test_roof_type_correct_matches_ground_truth():
    b = Building(id="b3")
    b.attributes["roof_type"] = "FLAT"
    assert roof_type_correct(b, GroundTruth(roof_type="FLAT")) is True
    assert roof_type_correct(b, GroundTruth(roof_type="GABLED")) is False


def test_roof_type_correct_none_prediction_is_false():
    b = Building(id="b4")
    assert roof_type_correct(b, GroundTruth(roof_type="FLAT")) is False
