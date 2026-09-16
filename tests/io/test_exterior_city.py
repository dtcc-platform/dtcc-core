"""The exterior profile preserves physical features without inventing a graph."""

import copy
import json

import numpy as np
import pytest

from dtcc_core import io
from dtcc_core.model import MultiLineString, Object, SemanticRegion, exchange
from dtcc_core.io.cityjson.admission import COMPOSITE_SURFACE_ROLE
from dtcc_core.io.cityjson.semantics import SEMANTIC_NAMESPACE


from sandbox.model_profiles.exterior_city_example import exterior_source as source


def feature(city, id):
    return next(child for group in city.children.values() for child in group if child.id == id)


def test_exterior_public_roundtrip_and_direct_access(tmp_path):
    path = tmp_path / 'exterior.city.json'
    path.write_text(json.dumps(source()))
    city = io.load_city(path, strict=True)
    road = feature(city, 'road')
    assert type(road) is Object
    assert type(road.get_geometry(lod='0')) is MultiLineString
    np.testing.assert_array_equal(road.get_geometry(lod='0').linestrings[0].vertices,
                                  [[0, 0, 0], [2, 0, 0], [2, 2, 0]])
    assert road.geometry['cityjson-1'].role == COMPOSITE_SURFACE_ROLE
    assert road.get_geometry(id='cityjson-1').regions[0].attributes == {
        'function': 'carriageway', 'surface_material': 'asphalt'}
    assert feature(city, 'plants').attributes['average_height'] == 1.5
    assert feature(city, 'water').get_geometry(lod='2.2').regions[0].attributes['water_level'] == 'meanWaterLevel'
    city.save(tmp_path / 'exterior.dtcc')
    payload = exchange.dumps(city)
    restored = io.load_city(tmp_path / 'exterior.dtcc')
    assert exchange.dumps(restored) == payload
    restored.save(path, strict=True)
    exported = json.loads(path.read_text())
    for id, original in source()['CityObjects'].items():
        assert exported['CityObjects'][id]['attributes'] == original.get('attributes', {})
        assert [g['type'] for g in exported['CityObjects'][id]['geometry']] == [
            g['type'] for g in original['geometry']]
    reloaded = io.load_city(path, strict=True)
    # CityObjects is a JSON map: its order has no external meaning.
    for id in source()['CityObjects']:
        assert exchange.dumps(feature(reloaded, id)) == exchange.dumps(feature(city, id))


def test_owner_vocabulary_and_unmapped_values_rejected_even_with_bypass(tmp_path):
    for owner, region in [('road', 'WaterSurface'), ('water', 'TrafficArea'),
                          ('plants', 'RoofSurface'), ('canal', 'TrafficArea')]:
        data = source()
        geometry = copy.deepcopy(data['CityObjects']['road']['geometry'][1])
        geometry['semantics'] = {'surfaces': [{'type': region}], 'values': [0, 0, 0]}
        data['CityObjects'][owner]['geometry'] = [geometry]
        with pytest.raises(NotImplementedError, match='semantic regions'):
            io.load_cityjson(data, strict=True, validate_schema=False)
    data = source()
    data['CityObjects']['road']['geometry'][1]['semantics']['surfaces'][0]['function'] = ['carriageway']
    with pytest.raises(ValueError, match='scalar'):
        io.load_cityjson(data, strict=True, validate_schema=False)
    data = source()
    data['CityObjects']['plants']['attributes']['average_height'] = 1.5
    with pytest.raises(ValueError, match='reserved'):
        io.load_cityjson(data, strict=True, validate_schema=False)

    city = io.load_cityjson(source(), strict=True)
    path = tmp_path / 'exterior.city.json'
    city.save(path, strict=True)
    before = path.read_bytes()
    feature(city, 'road').get_geometry(id='cityjson-1').regions[0].semantic_type = SEMANTIC_NAMESPACE + 'WaterSurface'
    with pytest.raises(NotImplementedError, match='semantic regions'):
        city.save(path, strict=True, validate_schema=False)
    assert path.read_bytes() == before


def test_lines_reject_unrepresentable_state_and_quantization_without_overwrite(tmp_path):
    for indices in ([], [[]], [[0]], [[True, 1]], [[0, 999]], [[-1, 0]]):
        data = source()
        data['CityObjects']['road']['geometry'][0]['boundaries'] = indices
        with pytest.raises(ValueError, match='line'):
            io.load_cityjson(data, strict=True, validate_schema=False)
    city = io.load_cityjson(source(), strict=True)
    path = tmp_path / 'lines.city.json'
    city.save(path, strict=True)
    before = path.read_bytes()
    lines = feature(city, 'road').get_geometry(lod='0')
    lines.linestrings[0].vertices = lines.linestrings[0].vertices[:, :2]
    with pytest.raises(ValueError, match='3D vertices'):
        city.save(path, strict=True, validate_schema=False)
    assert path.read_bytes() == before
    lines.linestrings[0].vertices = np.array([[0., 0., 0.], [.00001, 0., 0.]])
    with pytest.raises(ValueError, match='collapses a line'):
        city.save(path, strict=True, validate_schema=False)
    assert path.read_bytes() == before
    lines.linestrings[0].vertices[1, 0] = 1
    lines.regions = [SemanticRegion(SEMANTIC_NAMESPACE + 'TrafficArea', np.array([0]))]
    with pytest.raises(NotImplementedError, match='Semantic regions'):
        city.save(path, strict=True, validate_schema=False)
    assert path.read_bytes() == before



def test_plant_height_schema_failure_and_explicit_bypass(tmp_path):
    city = io.load_cityjson(source(), strict=True)
    native, external = tmp_path / 'exterior.dtcc', tmp_path / 'exterior.city.json'
    city.save(native)
    city.save(external, strict=True)
    before = {path: path.read_bytes() for path in (native, external)}
    feature(city, 'plants').attributes['average_height'] = -1.
    for path in before:
        with pytest.raises(ValueError, match='average_height'):
            city.save(path, **({'strict': True} if path == external else {}))
        assert path.read_bytes() == before[path]
    payload = exchange.dumps(city, validate_schema=False)
    with pytest.raises(ValueError, match='average_height'):
        exchange.loads(payload)
    assert feature(exchange.loads(payload, validate_schema=False), 'plants').attributes['average_height'] == -1.
