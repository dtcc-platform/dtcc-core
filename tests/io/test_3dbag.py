"""A source mapping must preserve evidence and never invent height or dates."""

from copy import deepcopy
import json
from pathlib import Path

import pytest

from dtcc_core import io
from dtcc_core.datasets import load_model_package
from dtcc_core.model import exchange
from dtcc_core.io.three_d_bag import NAP_CRS


def source():
    # Synthetic coordinates and supplier facts, not an actual survey.
    fixture = (
        Path(__file__).parents[2]
        / "sandbox/model_profiles/fixtures/buildings.city.json"
    )
    data = json.loads(fixture.read_text())
    data["metadata"]["referenceSystem"] = "https://www.opengis.net/def/crs/EPSG/0/7415"
    data["CityObjects"]["building-1"]["attributes"] = {
        "b3_h_maaiveld": -0.41,
        "b3_h_dak_50p": 11.34,
        "b3_h_dak_min": None,
        "b3_dak_type": "multiple horizontal",
        "b3_pw_bron": "AHN3",
        "b3_pw_datum": 2014,
        "b3_val3dity_lod22": "[303]",
    }
    return data


def test_mapping_public_boundaries_and_source_preservation(tmp_path):
    data = source()
    before = deepcopy(data)
    path = tmp_path / "source.city.json"
    path.write_text(json.dumps(data))
    city = io.load_3dbag(path)
    assert json.loads(path.read_text()) == before
    b = city.buildings[0]
    attrs = b.attributes
    assert data == before
    assert all(
        attrs[k] == v for k, v in data["CityObjects"][b.id]["attributes"].items()
    )
    assert b.height is None
    assert "height_measurements" not in attrs
    assert b.building_parts[0].attributes.get("elevation_measurements") is None
    records = attrs["elevation_measurements"]
    assert [r["value"] for r in records] == [-0.41, 11.34]
    assert all(r["vertical_reference"] == NAP_CRS and r["unit"] == "m" for r in records)
    assert attrs["roof_type"]["value"] == "multiple horizontal"
    assert attrs["b3_pw_datum"] == 2014 and attrs["b3_val3dity_lod22"] == "[303]"
    payload = exchange.dumps(city)
    city.save(tmp_path / "mapped.dtcc")
    assert exchange.dumps(io.load_city(tmp_path / "mapped.dtcc")) == payload
    package = city.export(tmp_path / "mapped.dtccpkg", canonical=True)
    packaged = load_model_package(package.path)
    assert exchange.dumps(packaged) == payload
    assert packaged.dataset_context == city.dataset_context
    city.save(tmp_path / "mapped.city.json", strict=True)
    restored = io.load_city(tmp_path / "mapped.city.json", strict=True)
    restored.id = city.id
    assert exchange.dumps(restored) == payload
    # Use the normal loader for enriched data; never silently remap destinations.
    with pytest.raises(ValueError, match="destination already exists"):
        io.load_3dbag(tmp_path / "mapped.city.json")


def test_source_contract_failures_and_missing_values():
    for value in (True, "11.34", {"value": 11.34}):
        data = source()
        data["CityObjects"]["building-1"]["attributes"]["b3_h_dak_50p"] = value
        with pytest.raises(ValueError, match="b3_h_dak_50p"):
            io.load_3dbag(data, validate_schema=False)
    data = source()
    data["metadata"]["referenceSystem"] = "EPSG:3006"
    with pytest.raises(ValueError, match="EPSG:7415"):
        io.load_3dbag(data)
    data = source()
    data["CityObjects"]["building-1"]["attributes"]["b3_h_50p"] = 11.34
    with pytest.raises(ValueError, match="newer 3DBAG"):
        io.load_3dbag(data)
    for roof in ("unknown", "no points", "no planes", "future_supplier_code", None):
        data = source()
        attrs = data["CityObjects"]["building-1"]["attributes"]
        attrs.update(b3_dak_type=roof, b3_h_maaiveld=None, b3_h_dak_50p=None)
        b = io.load_3dbag(data).buildings[0]
        assert (
            "roof_type" not in b.attributes
            and "elevation_measurements" not in b.attributes
        )
        assert b.attributes["b3_dak_type"] == roof
    data["CityObjects"]["building-1"]["attributes"] = {}
    with pytest.raises(ValueError, match="No supported"):
        io.load_3dbag(data)


def test_elevation_schema_applies_to_ordinary_save_load_and_bypass(tmp_path):
    city = io.load_3dbag(source())
    path = tmp_path / "mapped.dtcc"
    city.save(path)
    before = path.read_bytes()
    record = city.buildings[0].attributes["elevation_measurements"][0]
    original = deepcopy(record)
    for key, value in [
        ("vertical_reference", ""),
        ("unit", "s"),
        ("value", True),
        ("reference", {"value": "ground"}),
        ("typo", "extra"),
    ]:
        record.clear()
        record.update(original)
        record[key] = value
        with pytest.raises(ValueError, match="elevation_measurements"):
            city.save(path)
        assert path.read_bytes() == before
    record.clear()
    record.update(original)
    del record["vertical_reference"]
    with pytest.raises(ValueError, match="vertical_reference"):
        city.save(path)
    city.save(path, validate_schema=False)
    with pytest.raises(ValueError, match="vertical_reference"):
        io.load_city(path)
    assert io.load_city(path, validate_schema=False).buildings[0].height is None
