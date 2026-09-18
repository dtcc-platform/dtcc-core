"""Stored building height is independent of geometry and modelling estimates."""

import numpy as np
import pytest

from dtcc_core import io
from dtcc_core.model import Building, BuildingPart, GeometryType, Surface


@pytest.mark.parametrize("cls", [Building, BuildingPart])
def test_height_has_one_stored_authority_and_no_fallback(cls, tmp_path):
    value = cls(id="b", attributes={"height": 99, "estimated_height": 8})
    value.add_geometry(
        Surface(vertices=np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 3.0]])),
        GeometryType.LOD2,
    )
    assert value.height is None and value.measured_height is None
    assert (
        value.lod2.bounds.depth == 3
    )  # An explicitly selected local-coordinate extent.
    value.height = 12.5
    assert value.measured_height == value.attributes["measured_height"] == 12.5
    assert value.estimated_height == 8 and value.attributes["height"] == 99
    value.measured_height = 0
    assert value.height == 0
    path = tmp_path / "building.dtcc"
    io.save_model(value, path)
    restored = io.load_model(path)
    assert restored.height == 0 and restored.estimated_height == 8
    restored.height = None
    assert "measured_height" not in restored.attributes
    assert restored.height is None and restored.estimated_height == 8


def test_estimate_validation_occurs_at_persistence_and_preserves_file(tmp_path):
    value = Building(id="b", attributes={"measured_height": 12.5})
    value.estimated_height = 8
    path = tmp_path / "building.dtcc"
    io.save_model(value, path)
    before = path.read_bytes()
    value.estimated_height = -1
    assert (
        value.estimated_height == -1 and value.height == 12.5
    )  # Cheap access while editing.
    with pytest.raises(ValueError, match="estimated_height"):
        io.save_model(value, path)
    assert path.read_bytes() == before
    io.save_model(value, path, validate_schema=False)
    with pytest.raises(ValueError, match="estimated_height"):
        io.load_model(path)
    assert io.load_model(path, validate_schema=False).estimated_height == -1
