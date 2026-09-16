"""The standard building vocabulary works through ordinary persistence boundaries."""

import json
from pathlib import Path

import numpy as np
import pytest

from dtcc_core import io
from dtcc_core.model import Building, BuildingPart, GeometryType, MultiSurface, Object, SemanticRegion, Surface
from dtcc_core.model.profiles import SemanticProfile
from dtcc_core.model._standard_schema import DEFAULT_VERSION

ROOT = Path(__file__).resolve().parents[2]
NS = 'https://github.com/dtcc-platform/dtcc-core/schemas/model#'
SURFACES = ['ClosureSurface', 'OuterCeilingSurface', 'OuterFloorSurface',
            'InteriorWallSurface', 'CeilingSurface', 'FloorSurface']


@pytest.fixture(scope='module')
def profile():
    return SemanticProfile(ROOT / 'dtcc_core/schemas/dtcc.yaml')


def test_building_attributes_and_surfaces_through_cityjson_and_native_io(tmp_path):
    source = json.loads((ROOT / 'sandbox/model_profiles/fixtures/buildings.city.json').read_text())
    source['CityObjects']['building-1']['attributes'].update({
        'description': 'Mixed-use building', 'class': '1000',
        'function': ['housing', 'commerce'], 'usage': ['residential', 'office'],
        'storeysAboveGround': 3, 'storeysBelowGround': 0,
        'code_context': {'class': 'urn:example:building-codes'},
    })
    part = source['CityObjects']['part-1']
    part['attributes'].update({'class': 'wing', 'function': [], 'storeysAboveGround': 2})
    # Synthetic assignments exercise the vocabulary/membership contract, not orientation.
    part['geometry'][0]['semantics'] = {
        'surfaces': [{'type': name} for name in SURFACES], 'values': [0, 1, 2, 3, 4, 5, None]}
    input_path = tmp_path / 'input.city.json'
    input_path.write_text(json.dumps(source))
    city = io.load_city(input_path, strict=True)
    building = city.buildings[0]
    assert building.attributes['storeys_above_ground'] == 3
    assert building.attributes['storeys_below_ground'] == 0
    native_part = building.building_parts[0]
    assert type(native_part) is BuildingPart
    assert 'storeys_below_ground' not in native_part.attributes
    assert native_part.attributes['function'] == []
    assert [r.semantic_type for r in native_part.lod2.regions] == [NS + s for s in SURFACES]
    path = tmp_path / 'building.dtcc'
    city.save(path)
    restored = io.load_city(path)
    assert restored.schema_version == DEFAULT_VERSION == '0.9.0'
    output_path = tmp_path / 'output.city.json'
    restored.save(output_path, strict=True)
    result = json.loads(output_path.read_text())
    for key, feature in source['CityObjects'].items():
        assert result['CityObjects'][key]['attributes'] == feature['attributes']
        assert result['CityObjects'][key]['geometry'][0]['semantics'] == feature['geometry'][0]['semantics']
    assert io.load_city(output_path, strict=True).buildings[0].attributes == building.attributes


def test_declared_attributes_reject_bad_values_without_coercion(profile, tmp_path):
    # Same declarations apply to native parts and to a generic Object classified in YAML.
    building = Object(id='generic-building', semantic_type=NS + 'Building')
    part = BuildingPart(id='part')
    building.add_child(part)
    path = tmp_path / 'building.dtcc'
    io.save_model(building, path)
    before = path.read_bytes()
    cases = [
        ('description', ['text']), ('class', 1000),
        ('function', 'housing'), ('usage', ['residential', 7]),
        ('storeys_above_ground', True), ('storeys_above_ground', 1.5),
        ('storeys_above_ground', '2'), ('storeys_below_ground', -1),
    ]
    for slot, value in cases:
        part.attributes = {slot: value}
        report = profile.validate(building)
        assert not report.valid and any(slot in issue.path for issue in report.issues)
        assert part.attributes == {slot: value}  # Validation does not normalize.
    with pytest.raises(ValueError, match='storeys_below_ground'):
        io.save_model(building, path)
    assert path.read_bytes() == before
    io.save_model(building, path, validate_schema=False)
    with pytest.raises(ValueError, match='storeys_below_ground'):
        io.load_model(path)
    assert io.load_model(path, validate_schema=False).get_children(BuildingPart)[0].attributes == part.attributes


def test_optional_counts_use_json_schema_integrality_and_preserve_native_values(tmp_path):
    building = Building(id='b', attributes={
        'storeys_above_ground': 2.0, 'storeys_below_ground': None,
        'function': [], 'usage': None, 'custom': {'storeys_above_ground': 'untouched'},
    })
    path = tmp_path / 'building.dtcc'
    io.save_model(building, path)
    restored = io.load_model(path)
    assert restored.attributes == building.attributes
    assert type(restored.attributes['storeys_above_ground']) is float
    # JSON Schema integer means an integral numeric value; native encoding stays exact.


def test_new_surfaces_receive_owner_rules_and_abstract_building_is_not_concrete(profile, tmp_path):
    other = Object(id='other')
    geometry = MultiSurface(surfaces=[
        Surface(vertices=np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]]))])
    geometry.regions = [SemanticRegion(NS + name, np.array([], dtype=np.int64)) for name in SURFACES]
    other.add_geometry(geometry, GeometryType.LOD2)
    report = profile.validate(other)
    assert len(report.issues) == len(SURFACES)
    assert all(issue.rule == 'target_type' for issue in report.issues)
    with pytest.raises(ValueError, match='Schema validation failed'):
        io.save_model(other, tmp_path / 'invalid.dtcc')
    # Geometry roots remain valid owners for these same regions.
    io.save_model(geometry, tmp_path / 'surface.dtcc')
    assert len(io.load_model(tmp_path / 'surface.dtcc').regions) == len(SURFACES)
    other.geometry.clear()
    other.semantic_type = NS + 'AbstractBuilding'
    with pytest.raises(ValueError, match='Unsupported concrete semantic type'):
        io.save_model(other, tmp_path / 'abstract.dtcc')


@pytest.mark.parametrize('strict', [False, True])
def test_storey_mapping_rejects_ambiguous_spelling(tmp_path, strict):
    source = json.loads((ROOT / 'sandbox/model_profiles/fixtures/buildings.city.json').read_text())
    attrs = source['CityObjects']['building-1']['attributes']
    attrs['storeysAboveGround'] = 4
    path = tmp_path / 'input.city.json'
    path.write_text(json.dumps(source))
    city = io.load_city(path, strict=strict)
    assert city.buildings[0].attributes['storeys_above_ground'] == 4
    attrs['storeys_above_ground'] = 5
    path.write_text(json.dumps(source))
    with pytest.raises(ValueError, match='reserved.*mapping'):
        io.load_city(path, strict=strict)
    city.buildings[0].attributes['storeysAboveGround'] = 5
    with pytest.raises(ValueError, match='reserved.*mapping'):
        city.save(tmp_path / 'ambiguous.city.json', strict=strict)
