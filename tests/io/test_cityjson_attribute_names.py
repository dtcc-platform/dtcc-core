"""DTCC property names and external spellings have one explicit boundary."""

import json
from copy import deepcopy
from pathlib import Path

import pytest

from dtcc_core import io
from dtcc_core.model import Building, Object
from dtcc_core.model._standard_schema import DEFAULT_VERSION, SCHEMA_ID

FIXTURES = Path(__file__).resolve().parents[2] / "sandbox/model_profiles/fixtures"


def test_mixed_city_names_survive_external_and_canonical_boundaries(tmp_path):
    source = json.loads((FIXTURES / "mixed-city.city.json").read_text())
    objects = source["CityObjects"]
    objects["building-1"]["attributes"].update(
        measuredHeight=12.5,
        estimated_height=8.0,
        roofType="gable",
        customLabel={"measuredHeight": "untouched", "storeysAboveGround": "untouched"},
    )
    plant = next(
        value
        for value in objects.values()
        if value["type"] == "SolitaryVegetationObject"
    )
    plant.setdefault("attributes", {}).update(crownDiameter=4.0, trunkDiameter=0.5)
    bench = next(
        value for value in objects.values() if value["type"] == "CityFurniture"
    )
    bench.setdefault("attributes", {}).update(
        measuredHeight="unrelated user metadata",
        storeysAboveGround="unrelated user metadata",
    )
    path = tmp_path / "source.city.json"
    path.write_text(json.dumps(source))
    city = io.load_city(path, strict=True)
    assert city.buildings[0].height == 12.5
    assert city.buildings[0].estimated_height == 8.0
    attrs = city.buildings[0].attributes
    assert attrs["measured_height"] == 12.5 and attrs["roof_type"] == "gable"
    assert "measuredHeight" not in attrs and "roofType" not in attrs
    assert attrs["customLabel"] == {
        "measuredHeight": "untouched",
        "storeysAboveGround": "untouched",
    }
    native_plant = next(
        obj
        for obj in city.get_children(Object)
        if obj.semantic_type.endswith("#SolitaryVegetationObject")
    )
    assert native_plant.attributes["crown_diameter"] == 4.0
    assert native_plant.attributes["trunk_diameter"] == 0.5
    city.save(tmp_path / "city.dtcc")
    restored = io.load_city(tmp_path / "city.dtcc")
    assert restored.schema_version == DEFAULT_VERSION == "0.9.0"
    assert restored.buildings[0].height == 12.5
    assert restored.buildings[0].estimated_height == 8.0
    restored.save(tmp_path / "result.city.json", strict=True)
    exported = json.loads((tmp_path / "result.city.json").read_text())
    assert {
        key: value.get("attributes", {})
        for key, value in exported["CityObjects"].items()
    } == {key: value.get("attributes", {}) for key, value in objects.items()}

    before = (tmp_path / "city.dtcc").read_bytes()
    restored.buildings[0].attributes["measured_height"] = -1
    with pytest.raises(ValueError, match="measured_height"):
        restored.save(tmp_path / "city.dtcc")
    assert (tmp_path / "city.dtcc").read_bytes() == before


@pytest.mark.parametrize("strict", [False, True])
def test_ordinary_and_strict_building_names_reject_lossy_collisions(tmp_path, strict):
    source = json.loads((FIXTURES / "buildings.city.json").read_text())
    path = tmp_path / "source.city.json"
    path.write_text(json.dumps(source))
    city = io.load_city(path, strict=strict)
    assert city.buildings[0].attributes["measured_height"] == 10.2
    target = tmp_path / "result.city.json"
    city.save(target, strict=strict)
    assert (
        json.loads(target.read_text())["CityObjects"]["building-1"]["attributes"][
            "measuredHeight"
        ]
        == 10.2
    )
    before = target.read_bytes()

    for include_external in (True, False):
        ambiguous = deepcopy(source)
        attrs = ambiguous["CityObjects"]["building-1"]["attributes"]
        attrs["measured_height"] = 99
        if not include_external:
            del attrs["measuredHeight"]
        path.write_text(json.dumps(ambiguous))
        with pytest.raises(ValueError, match="reserved.*mapping"):
            io.load_city(path, strict=strict)
    city.buildings[0].attributes["measuredHeight"] = 99
    with pytest.raises(ValueError, match="reserved.*mapping"):
        city.save(target, strict=strict)
    assert target.read_bytes() == before


def test_retired_standard_versions_and_generic_metadata_are_not_migrated(tmp_path):
    value = Building(
        id="b", attributes={"measuredHeight": "user metadata", "customLabel": 7}
    )
    path = tmp_path / "building.dtcc"
    io.save_model(value, path)
    assert (
        io.load_model(path).attributes == value.attributes
    )  # No global case conversion.
    value.schema_id, value.schema_version = SCHEMA_ID, "0.2.0"
    with pytest.raises(ValueError, match="Unsupported semantic schema"):
        io.save_model(value, path)
    io.save_model(value, path, validate_schema=False)
    with pytest.raises(ValueError, match="Unsupported semantic schema"):
        io.load_model(path)
    restored = io.load_model(path, validate_schema=False)
    assert (
        restored.schema_version == "0.2.0" and restored.attributes == value.attributes
    )
