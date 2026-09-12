"""Bounded exterior city, through the ordinary Python and persistence APIs.

Run: python -m sandbox.model_profiles.exterior_city_example OUTPUT
All coordinates and attributes in this example are synthetic tutorial data.
"""

import argparse
import copy
import json
from pathlib import Path

import numpy as np

from dtcc_core import io
from dtcc_core.datasets import load_model_package
from dtcc_core.datasets.schema import DatasetContext
from dtcc_core.model import Mesh, MultiLineString, Point, Solid, exchange

from .tin_relief_example import source as tin_source

SEMANTIC_NAMESPACE = 'https://github.com/dtcc-platform/dtcc-core/schemas/model#'


def exterior_source():
    def polygons(kind, names=()):
        result = {'type': kind, 'lod': '2.2',
                  'boundaries': [[[0, 1, 2]], [[0, 2, 3]], [[0, 3, 4]]]}
        if kind == 'Solid':
            result['boundaries'] = [[[0, 2, 1]], [[0, 1, 4]], [[1, 2, 4]], [[2, 0, 4]]]
        if names:
            result['semantics'] = {'surfaces': [{'type': name} for name in names],
                                   'values': [i % len(names) for i in range(len(result['boundaries']))]}
        if kind == 'Solid':
            result['boundaries'] = [result['boundaries']]
            if names:
                result['semantics']['values'] = [result['semantics']['values']]
        return result

    def lines():
        return {'type': 'MultiLineString', 'lod': '0', 'boundaries': [[0, 1, 2], [3, 4]]}

    road = polygons('CompositeSurface', ('TrafficArea', 'AuxiliaryTrafficArea', 'TransportationMarking'))
    road['semantics']['surfaces'][0].update(function='carriageway', surfaceMaterial='asphalt')
    hole = polygons('MultiSurface', ('TransportationHole',))
    water = polygons('Solid', ('WaterSurface', 'WaterGroundSurface', 'WaterClosureSurface'))
    water['semantics']['surfaces'][0]['waterLevel'] = 'meanWaterLevel'
    return {'type': 'CityJSON', 'version': '2.0',
            'transform': {'scale': [1, 1, 1], 'translate': [0, 0, 0]},
            'metadata': {}, 'vertices': [[0, 0, 0], [2, 0, 0], [2, 2, 0], [0, 2, 0], [0, 0, -1]],
            'CityObjects': {
                'road': {'type': 'Road', 'attributes': {'function': ['road'],
                    'class': {'value': 'primary', 'code_space': 'https://example.test/roads'}},
                         'geometry': [lines(), road, hole]},
                'rail': {'type': 'Railway', 'geometry': [polygons('MultiSurface', ('TrafficArea',))]},
                'canal': {'type': 'Waterway', 'geometry': [lines()]},
                'square': {'type': 'TransportSquare', 'geometry': [polygons('MultiSurface', ('AuxiliaryTrafficArea',))]},
                'water': {'type': 'WaterBody', 'geometry': [water]},
                'plants': {'type': 'PlantCover', 'attributes': {'averageHeight': 1.5},
                           'geometry': [polygons('MultiSurface'), polygons('Solid')]},
            }}



def source():
    """Combine the existing mixed and TIN fixtures with the exterior features."""
    data = json.loads((Path(__file__).parent / 'fixtures/mixed-city.city.json').read_text())

    def append(document, ids, offset=(0., 0., 0.)):
        start = len(data['vertices'])
        points = np.asarray(document['vertices'], dtype=np.float64) * document['transform']['scale']
        points += np.asarray(document['transform']['translate']) + offset
        target = (points - data['transform']['translate']) / data['transform']['scale']
        rounded = np.rint(target)
        # The authored examples share an exact millimetre grid.
        np.testing.assert_allclose(target, rounded, atol=1e-5, rtol=0)
        data['vertices'].extend(rounded.astype(np.int64).tolist())

        def shift(boundaries):
            return [shift(value) if isinstance(value, list) else value + start for value in boundaries]

        for id in ids:
            feature = copy.deepcopy(document['CityObjects'][id])
            for geometry in feature.get('geometry', []):
                geometry['boundaries'] = shift(geometry['boundaries'])
            data['CityObjects'][id] = feature

    tin = tin_source()
    append(tin, ['terrain-west', 'terrain-east'])
    exterior = exterior_source()
    origin = np.asarray(data['transform']['translate']) + [40., 0., 0.]
    append(exterior, exterior['CityObjects'], offset=origin)
    return data


def objects(city):
    result, stack = {}, [city]
    while stack:
        parent = stack.pop()
        for group in parent.children.values():
            for child in group:
                result[child.id] = child
                stack.append(child)
    return result


def verify(city, restored, tolerance=.00050001):
    """Check meanings and geometry, permitting only declared quantization error."""
    expected, actual = objects(city), objects(restored)
    assert expected.keys() == actual.keys()
    assert city.transform.srs == restored.transform.srs
    for id, value in expected.items():
        other = actual[id]
        assert type(value) is type(other)
        assert value.attributes == other.attributes and value.semantic_type == other.semantic_type
        assert value.transform.srs == other.transform.srs
        assert {child.id for group in value.children.values() for child in group} == {
            child.id for group in other.children.values() for child in group}
        assert list(value.geometry) == list(other.geometry)
        for key, record in value.geometry.items():
            target = other.geometry[key]
            assert (record.lod, record.role) == (target.lod, target.role)
            a, b = record.geometry, target.geometry
            assert type(a) is type(b) and a.transform.srs == b.transform.srs
            assert len(a.regions) == len(b.regions)
            for r, s in zip(a.regions, b.regions):
                assert (r.semantic_type, r.id, r.parent, r.attributes) == (
                    s.semantic_type, s.id, s.parent, s.attributes)
                np.testing.assert_array_equal(r.indices, s.indices)
            if type(a) is Point:
                pairs = [([a.x, a.y, a.z], [b.x, b.y, b.z])]
            elif type(a) is MultiLineString:
                assert len(a.linestrings) == len(b.linestrings)
                pairs = [(x.vertices, y.vertices) for x, y in zip(a.linestrings, b.linestrings)]
            elif type(a) is Mesh:
                pairs = [(a.vertices[a.faces], b.vertices[b.faces])]
            else:
                assert len(a.surfaces) == len(b.surfaces)
                pairs = []
                for x, y in zip(a.surfaces, b.surfaces):
                    assert len(x.holes) == len(y.holes)
                    pairs.extend(zip([x.vertices, *x.holes], [y.vertices, *y.holes]))
                if type(a) is Solid:
                    assert len(a.shells) == len(b.shells)
                    for x, y in zip(a.shells, b.shells):
                        np.testing.assert_array_equal(x, y)
            for a_coordinates, b_coordinates in pairs:
                np.testing.assert_allclose(a_coordinates, b_coordinates, atol=tolerance, rtol=0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    path = args.output / 'source.city.json'
    data = source()
    path.write_text(json.dumps(data, indent=2) + '\n')
    city = io.load_city(path, strict=True)
    features = objects(city)
    road = features['road']
    line_vertices = road.get_geometry(lod='0').linestrings[0].vertices
    surface = road.get_geometry(id='cityjson-1')
    assert line_vertices.shape == (3, 3)
    assert surface.regions[0].attributes['surface_material'] == 'asphalt'
    assert road.geometry['cityjson-1'].role == 'cityjson:CompositeSurface'
    assert features['bench-1'].id == 'bench-1'
    assert features['plants'].attributes['average_height'] == 1.5
    assert features['water'].get_geometry(lod='2.2').regions[0].attributes['water_level'] == 'meanWaterLevel'

    city.dataset_context = DatasetContext(
        identity={'name': 'exterior-city-example', 'title': 'Synthetic exterior city'},
        metadata={'crs': [city.transform.srs]},
        provenance={'sources': [{'description': 'Synthetic tutorial data, not a survey'}]},
        presentation={}, request={'dataset_name': 'exterior-city-example'},
    )
    payload = exchange.dumps(city)
    city.save(args.output / 'exterior.dtcc')
    assert exchange.dumps(io.load_city(args.output / 'exterior.dtcc')) == payload
    package = city.export(args.output / 'exterior.dtccpkg', canonical=True)
    restored = load_model_package(package.path)
    assert exchange.dumps(restored) == payload
    assert restored.dataset_context == city.dataset_context

    path = args.output / 'exterior.city.json'
    restored.save(path, strict=True)
    verify(city, io.load_city(path, strict=True))
    exported = json.loads(path.read_text())
    for id, original in data['CityObjects'].items():
        assert [g['type'] for g in exported['CityObjects'][id]['geometry']] == [
            g['type'] for g in original.get('geometry', [])]
    before = path.read_bytes()
    surface.regions[0].semantic_type = SEMANTIC_NAMESPACE + 'WaterSurface'
    try:
        city.save(path, strict=True, validate_schema=False)
    except NotImplementedError as error:
        assert 'semantic regions' in str(error)
    else:
        raise AssertionError('Cross-theme surface vocabulary was accepted')
    assert path.read_bytes() == before
    surface.regions[0].semantic_type = SEMANTIC_NAMESPACE + 'TrafficArea'
    features['plants'].attributes['average_height'] = -1.
    native_path = args.output / 'exterior.dtcc'
    native_before = native_path.read_bytes()
    for target in (native_path, path):
        try:
            city.save(target, **({'strict': True} if target == path else {}))
        except ValueError as error:
            assert 'average_height' in str(error)
        else:
            raise AssertionError('Negative average vegetation height was accepted')
    assert native_path.read_bytes() == native_before and path.read_bytes() == before
    unchecked = exchange.dumps(city, validate_schema=False)
    assert objects(exchange.loads(unchecked, validate_schema=False))['plants'].attributes['average_height'] == -1.
    features['plants'].attributes['average_height'] = 1.5

    report = {'schema_version': city.schema_version, 'wire_version': exchange.VERSION,
              'feature_types': sorted({value['type'] for value in data['CityObjects'].values()}),
              'features': len(features), 'representations': sum(len(v.geometry) for v in features.values()),
              'native_bytes': len(payload), 'native_and_package_exact': True,
              'strict_cityjson_verified': True, 'coordinate_tolerance': .00050001,
              'owner_failure_preserved_file_under_schema_bypass': True,
              'schema_failure_preserved_native_and_cityjson': True, 'semantic_bypass_checked': True}
    (args.output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
