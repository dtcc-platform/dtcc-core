from dtcc_core.model.object.building import Building
from dtcc_core.model.object.object import GeometryType
from dtcc_core.model.geometry.surface import MultiSurface
from dtcc_core.builder.evaluation.metrics import stage_outcome


def test_stage_outcome_success_when_lod2_present_and_no_fallback():
    b = Building(id="b1")
    b.add_geometry(MultiSurface(), GeometryType.LOD2)
    assert stage_outcome(b) == "success"


def test_stage_outcome_reads_fallback_reason():
    b = Building(id="b2")
    b.attributes["fallback_reason"] = "insufficient_points"
    assert stage_outcome(b) == "insufficient_points"


def test_stage_outcome_unknown_when_no_lod2_and_no_reason():
    b = Building(id="b3")
    assert stage_outcome(b) == "unknown"
