"""Synthetic mixed building/terrain example using ordinary strict exchange APIs.

Run: python -m sandbox.model_profiles.tin_relief_example OUTPUT
"""

import argparse
import json
from pathlib import Path

import numpy as np

from dtcc_core import io
from dtcc_core.datasets import load_model_package
from dtcc_core.datasets.schema import (DatasetContext, DatasetIdentity, DatasetMetadata,
                                      DatasetProvenance, DatasetPresentation, DatasetRequest)
from dtcc_core.model import Terrain, Mesh, exchange


def source():
    data = json.loads((Path(__file__).parent / 'fixtures/buildings.city.json').read_text())
    start = len(data['vertices'])
    data['vertices'].extend([[0, 0, 0], [10000, 0, 1000], [10000, 10000, 2000],
                             [0, 10000, -1000], [20000, 0, 0],
                             [30000, 0, 1000], [30000, 10000, 2000]])
    data['CityObjects']['terrain-west'] = {
        'type': 'TINRelief', 'attributes': {'name': 'Synthetic west terrain', 'source': {'synthetic': True}},
        'geometry': [
            {'type': 'CompositeSurface', 'lod': '1.2',
             'boundaries': [[[start, start+1, start+2]], [[start, start+2, start+3]]]},
            {'type': 'CompositeSurface', 'lod': '2.2',
             'boundaries': [[[start, start+1, start+3]], [[start+1, start+2, start+3]]]},
        ],
    }
    data['CityObjects']['terrain-east'] = {
        'type': 'TINRelief', 'attributes': {'name': 'Synthetic east terrain'},
        'geometry': [{'type': 'CompositeSurface', 'lod': '1',
                      'boundaries': [[[start+4, start+5, start+6]]]}],
    }
    return data


def verify(city, restored):
    """Compare semantic geometry, allowing only CityJSON coordinate quantization."""
    originals = {t.id: t for t in city.get_children(Terrain)}
    results = {t.id: t for t in restored.get_children(Terrain)}
    assert originals.keys() == results.keys() == {'terrain-west', 'terrain-east'}
    for id, a in originals.items():
        b = results[id]
        assert a.attributes == b.attributes and a.semantic_type == b.semantic_type
        assert a.transform.srs == b.transform.srs
        assert list(a.geometry) == list(b.geometry)
        for key, record in a.geometry.items():
            other = b.geometry[key]
            assert (record.lod, record.role) == (other.lod, other.role)
            x, y = record.geometry, other.geometry
            assert type(x) is type(y) is Mesh and x.transform.srs == y.transform.srs
            np.testing.assert_allclose(x.vertices[x.faces], y.vertices[y.faces], atol=.00050001, rtol=0)
    # The existing building/part is also preserved through the mixed adapter path.
    for a, b in zip(city.buildings, restored.buildings):
        assert exchange.dumps(a) == exchange.dumps(b)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    path = args.output / 'source.city.json'
    path.write_text(json.dumps(source(), indent=2) + '\n')
    city = io.load_city(path, strict=True)
    west = next(t for t in city.get_children(Terrain) if t.id == 'terrain-west')
    mesh = west.get_geometry(id='cityjson-1')
    assert type(mesh) is Mesh and mesh.vertices.shape == (4, 3) and mesh.faces.shape == (2, 3)
    assert west.get_geometry(lod='2.2') is mesh
    city.dataset_context = DatasetContext(
        identity=DatasetIdentity(name='tin-example', title='Synthetic mixed terrain example'),
        metadata=DatasetMetadata(crs=[city.transform.srs]),
        provenance=DatasetProvenance(sources=[{'description': 'Synthetic tutorial data; not a terrain survey'}]),
        presentation=DatasetPresentation(), request=DatasetRequest(dataset_name='tin-example'),
    )
    payload = exchange.dumps(city)
    city.save(args.output / 'mixed.dtcc')
    assert exchange.dumps(io.load_city(args.output / 'mixed.dtcc')) == payload
    package = city.export(args.output / 'mixed.dtccpkg', canonical=True)
    restored = load_model_package(package.path)
    assert exchange.dumps(restored) == payload and restored.dataset_context == city.dataset_context
    target = args.output / 'mixed.city.json'
    city.save(target, strict=True)
    verify(city, io.load_city(target, strict=True))
    output = json.loads(target.read_text())
    for obj in output['CityObjects'].values():
        if obj['type'] == 'TINRelief':
            assert all(g['type'] == 'CompositeSurface' and 'semantics' not in g for g in obj['geometry'])
    before = target.read_bytes()
    mesh.markers = np.array([1, 2])
    try:
        city.save(target, strict=True)
    except NotImplementedError as error:
        assert 'markers' in str(error)
    else:
        raise AssertionError('Unmapped marker data was accepted')
    assert target.read_bytes() == before
    report = {'schema_version': city.schema_version, 'wire_version': exchange.VERSION,
              'terrain_features': 2, 'terrain_representations': 3, 'terrain_triangles': 5,
              'native_bytes': len(payload),
              'verification': 'Native/package and context exact; mixed strict CityJSON checked; marker rejection preserves existing output.'}
    (args.output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
