"""Qualified metadata uses one schema validator at ordinary exchange boundaries."""

from copy import deepcopy
from pathlib import Path

import pytest
import yaml

from dtcc_core import io
from dtcc_core.datasets import load_model_package
from dtcc_core.model import Building, Object, exchange
from dtcc_core.model.profiles import SemanticProfile
from dtcc_core.model._standard_schema import DEFAULT_VERSION
from sandbox.model_profiles.qualified_values_example import example_city

ROOT = Path(__file__).resolve().parents[2]
SCHEMA = ROOT / f'dtcc_core/schemas/model/{DEFAULT_VERSION}/schema.yaml'


def test_qualified_values_survive_public_boundaries_without_inference(tmp_path):
    city = example_city()
    payload = exchange.dumps(city)
    city.save(tmp_path / 'qualified.dtcc')
    restored = io.load_city(tmp_path / 'qualified.dtcc')
    assert exchange.dumps(restored) == payload
    assert restored.schema_version == '0.9.0' and exchange.VERSION == 6
    building = restored.buildings[0]
    assert building.height is None and building.estimated_height == 9.0
    assert building.attributes['height_measurements'][1]['value'] == 1180
    package = city.export(tmp_path / 'qualified.dtccpkg', canonical=True)
    packaged = load_model_package(package.path)
    assert exchange.dumps(packaged) == payload
    assert packaged.dataset_context == city.dataset_context
    city.save(tmp_path / 'qualified.city.json', strict=True)
    result = io.load_city(tmp_path / 'qualified.city.json', strict=True)
    result.id = city.id
    assert exchange.dumps(result) == payload
    # Editing one qualified claim neither updates the scalar nor aliases source data.
    building.attributes['height_measurements'][0]['value'] = 13
    assert building.height is None
    assert city.buildings[0].attributes['height_measurements'][0]['value'] == 12.5


def test_nested_record_errors_have_native_paths_and_do_not_coerce():
    profile = SemanticProfile(SCHEMA)
    city = example_city()
    part = city.buildings[0].building_parts[0]
    measurement = city.buildings[0].attributes['height_measurements'][0]
    cases = [('value', True), ('value', -1), ('value', '12.5'), ('value', None),
             ('unit', 's'), ('unit', ''), ('status', 'calculated'),
             ('high_reference', ''), ('low_reference', {'value': 'ground'})]
    for slot, value in cases:
        record = {**measurement, slot: value}
        part.attributes = {'height_measurements': [record]}
        report = profile.validate(city)
        assert not report.valid
        assert any(f"objects['part-1'].attributes['height_measurements'][0]['{slot}']" in issue.path
                   for issue in report.issues)
        assert part.attributes['height_measurements'][0] == record
    for missing in ('value', 'unit', 'status', 'high_reference', 'low_reference'):
        record = {k: v for k, v in measurement.items() if k != missing}
        part.attributes = {'height_measurements': [record]}
        assert any(issue.path.endswith(f"[0]['{missing}']") and issue.rule == 'required'
                   for issue in profile.validate(city).issues)
    for value in ({'value': 'a'}, {'value': 10, 'code_space': 'codes'},
                  {'value': 'a', 'code_space': ' '},
                  {'value': 'a', 'code_space': 'codes', 'codeSpace': 'typo'}):
        part.attributes = {'class': value}
        assert not profile.validate(city).valid
    for value in ({'value': 1}, [None], ['not a record']):
        part.attributes = {'height_measurements': value}
        assert not profile.validate(city).valid
    part.attributes = {'height_measurements': [], 'class': '', 'usage': None}
    assert profile.validate(city).valid


def test_invalid_edit_save_load_bypass_and_file_preservation(tmp_path):
    city = example_city()
    path = tmp_path / 'qualified.dtcc'
    city.save(path)
    before = path.read_bytes()
    city.buildings[0].attributes['height_measurements'][0]['unit'] = 'unknown_unit'
    with pytest.raises(ValueError, match='unit'):
        city.save(path)
    assert path.read_bytes() == before
    city.save(path, validate_schema=False)
    with pytest.raises(ValueError, match='unit'):
        io.load_city(path)
    assert io.load_city(path, validate_schema=False).buildings[0].attributes == city.buildings[0].attributes
    json_path = tmp_path / 'invalid.city.json'
    city.save(json_path, strict=True, validate_schema=False)
    with pytest.raises(ValueError, match='unit'):
        io.load_city(json_path, strict=True)
    before = json_path.read_bytes()
    with pytest.raises(ValueError, match='unit'):
        city.save(json_path, strict=True)
    assert json_path.read_bytes() == before
    city.buildings[0].relations['serves'] = ['missing']
    with pytest.raises(ValueError, match='dangling ID'):
        city.save(path, validate_schema=False)


def test_record_support_is_schema_driven_and_does_not_create_entity_graphs(tmp_path):
    schema = yaml.safe_load(SCHEMA.read_text())
    schema['classes']['SurveyNote'] = {'attributes': {
        'confidence': {'range': 'float', 'required': True, 'minimum_value': 0, 'maximum_value': 1}}}
    schema['classes']['AbstractBuilding']['attributes']['survey_note'] = {'range': 'SurveyNote', 'inlined': True}
    path = tmp_path / 'schema.yaml'
    path.write_text(yaml.safe_dump(schema))
    profile = SemanticProfile(path)
    building = Building(id='b', attributes={'survey_note': {'confidence': .8}})
    assert profile.validate(building).valid
    building.attributes['survey_note']['confidence'] = 2
    assert profile.validate(building).issues[0].path == "objects['b'].attributes['survey_note']['confidence']"
    building = Object(id='b', semantic_type='https://github.com/dtcc-platform/dtcc-core/schemas/model#QualifiedCode',
                      attributes={'value': 'a', 'code_space': 'codes'})
    with pytest.raises(ValueError, match='Unsupported concrete semantic type'):
        io.save_model(building, tmp_path / 'record-as-object.dtcc')
    for extension, message in [
        ({'range': 'Building', 'inlined': True}, 'identifier-free'),
        ({'range': 'Building'}, 'cannot contain ID references'),
        ({'range': 'QualifiedCode'}, 'explicit inlining'),
    ]:
        changed = deepcopy(schema)
        changed['classes']['SurveyNote']['attributes']['target'] = extension
        path.write_text(yaml.safe_dump(changed))
        with pytest.raises(ValueError, match=message):
            SemanticProfile(path)
    schema['classes']['IdentifiedNote'] = {
        'is_a': 'SurveyNote', 'attributes': {'note_id': {'identifier': True}}}
    path.write_text(yaml.safe_dump(schema))
    with pytest.raises(ValueError, match='identifier-free'):
        SemanticProfile(path)
