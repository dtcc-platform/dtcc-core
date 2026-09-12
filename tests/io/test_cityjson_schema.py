"""Strict CityJSON evaluates the standard contract without weakening format admission."""

from copy import deepcopy
import json
from pathlib import Path
import zipfile

import pytest

from dtcc_core import io
from dtcc_core.io.cityjson.cityjson import load
from dtcc_core.io.cityjson.write_cityjson import to_cityjson
from dtcc_core.model._standard_schema import SCHEMA_ID, DEFAULT_VERSION

FIXTURES = Path(__file__).resolve().parents[2] / 'sandbox/model_profiles/fixtures'


@pytest.fixture
def source():
    return json.loads((FIXTURES / 'buildings.city.json').read_text())


def write_source(path, source):
    text = json.dumps(source)
    if path.suffix == '.zip':
        with zipfile.ZipFile(path, 'w') as archive:
            archive.writestr('city.json', text)
    else:
        path.write_text(text)


@pytest.mark.parametrize('suffix', ['.city.json', '.city.json.zip'])
def test_public_default_schema_and_bypass_preserve_files(source, tmp_path, suffix):
    path = tmp_path / ('city' + suffix)
    write_source(path, source)
    city = io.load_city(path, strict=True)
    assert (city.schema_id, city.schema_version) == (SCHEMA_ID, DEFAULT_VERSION)
    city.save_cityjson(path, strict=True)
    before = path.read_bytes()
    part = city.buildings[0].building_parts[0]
    part.attributes['storeys_above_ground'] = True
    with pytest.raises(ValueError, match=r"objects\['part-1'\].*storeys_above_ground"):
        city.save(path, strict=True)
    assert path.read_bytes() == before
    city.save_cityjson(path, strict=True, validate_schema=False)
    invalid = path.read_bytes()
    with pytest.raises(ValueError, match='storeys_above_ground'):
        io.load_city(path, strict=True)
    assert path.read_bytes() == invalid
    restored = io.load_city(path, strict=True, validate_schema=False)
    assert restored.buildings[0].building_parts[0].attributes['storeys_above_ground'] is True
    assert (restored.schema_id, restored.schema_version) == (SCHEMA_ID, DEFAULT_VERSION)
    with pytest.raises(ValueError, match='storeys_above_ground'):
        restored.save(tmp_path / 'invalid.dtcc')
    restored.save(tmp_path / 'bypassed.dtcc', validate_schema=False)
    assert io.load_city(tmp_path / 'bypassed.dtcc', validate_schema=False).schema_version == DEFAULT_VERSION


def test_direct_dictionary_boundaries_and_open_metadata(source):
    attrs = source['CityObjects']['building-1']['attributes']
    attrs['custom'] = {'usage': 123, 'arbitrary': ['kept']}
    attrs['measuredHeight'] = -1
    before = deepcopy(source)
    with pytest.raises(ValueError, match='measured_height'):
        load(source, strict=True)
    city = load(source, strict=True, validate_schema=False)
    assert source == before
    with pytest.raises(ValueError, match='measured_height'):
        to_cityjson(city, strict=True)
    exported = to_cityjson(city, strict=True, validate_schema=False)
    assert exported['CityObjects']['building-1']['attributes'] == attrs
    city.buildings[0].height = 12
    assert to_cityjson(city, strict=True)['CityObjects']['building-1']['attributes']['custom'] == attrs['custom']


def test_schema_bypass_does_not_bypass_format_or_native_integrity(source, tmp_path):
    malformed = deepcopy(source)
    malformed['CityObjects']['building-1']['geometry'][0]['boundaries'][0][0][0] = 999999
    with pytest.raises((ValueError, IndexError)):
        load(malformed, strict=True, validate_schema=False)
    unsupported = deepcopy(source)
    unsupported['CityObjects']['building-1']['address'] = []
    with pytest.raises(NotImplementedError, match='unsupported content'):
        load(unsupported, strict=True, validate_schema=False)
    duplicate = tmp_path / 'duplicate.city.json'
    duplicate.write_text('{"type":"CityJSON","type":"CityJSON"}')
    with pytest.raises(ValueError, match='Duplicate CityJSON key'):
        io.load_city(duplicate, strict=True, validate_schema=False)
    city = load(source, strict=True)
    target = tmp_path / 'city.city.json'
    city.save(target, strict=True)
    before = target.read_bytes()
    city.buildings[0].building_parts[0].lod2.regions[0].indices[0] = 999999
    with pytest.raises(ValueError):
        city.save(target, strict=True, validate_schema=False)
    assert target.read_bytes() == before


def test_opening_host_is_semantic_but_dangling_host_is_structural():
    source = json.loads((FIXTURES / 'openings.city.json').read_text())
    semantics = source['CityObjects']['part-1']['geometry'][0]['semantics']
    semantics['surfaces'][0]['type'] = 'GroundSurface'  # Children still refer to this host.
    with pytest.raises(ValueError, match='expected'):
        load(source, strict=True)
    city = load(source, strict=True, validate_schema=False)
    with pytest.raises(ValueError, match='expected'):
        to_cityjson(city, strict=True)
    to_cityjson(city, strict=True, validate_schema=False)
    city.buildings[0].building_parts[0].lod3.regions[1].parent = 999
    with pytest.raises(ValueError, match='Region parent'):
        to_cityjson(city, strict=True, validate_schema=False)


def test_export_honors_stored_schema_and_never_mutates_selection(source):
    city = load(source, strict=True)
    city.schema_id = city.schema_version = None
    to_cityjson(city, strict=True)
    assert city.schema_id is None and city.schema_version is None
    city.schema_id, city.schema_version = SCHEMA_ID, '0.4.0'
    with pytest.raises(ValueError, match='Unsupported semantic schema'):
        to_cityjson(city, strict=True)
    exported = to_cityjson(city, strict=True, validate_schema=False)
    assert city.schema_version == '0.4.0'
    # CityJSON carries no DTCC declaration; reimport selects the current default.
    assert load(exported, strict=True).schema_version == DEFAULT_VERSION
    city.schema_version = None
    with pytest.raises((ValueError, TypeError), match='schema_version'):
        to_cityjson(city, strict=True, validate_schema=False)


def test_validation_option_cannot_claim_permissive_validation(source, tmp_path):
    path = tmp_path / 'city.city.json'
    write_source(path, source)
    city = load(source)
    assert city.schema_id is None
    for call in [lambda flag: io.load_city(path, validate_schema=flag),
                 lambda flag: city.save(path, validate_schema=flag),
                 lambda flag: load(source, validate_schema=flag),
                 lambda flag: to_cityjson(city, validate_schema=flag)]:
        with pytest.raises(ValueError, match='requires strict=True'):
            call(True)
        with pytest.raises(TypeError, match='True or False'):
            call('false')
    with pytest.raises(TypeError, match='True or False'):
        load(source, strict=True, validate_schema=0)
    assert load(source, validate_schema=False).schema_id is None
