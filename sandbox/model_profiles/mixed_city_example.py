"""Mixed source/native city workflows, retaining the distinction between plant and Tree."""

import argparse
import json
from pathlib import Path

import numpy as np

from dtcc_core import io
from dtcc_core.datasets import load_model_package
from dtcc_core.datasets.schema import (DatasetContext, DatasetIdentity, DatasetMetadata,
                                      DatasetPresentation, DatasetProvenance, DatasetRequest)
from dtcc_core.model import Object, Landuse, Tree, exchange
from dtcc_core.model.object.landuse import LanduseClasses

from buildings_example import project

HERE = Path(__file__).resolve().parent
NS = 'https://github.com/dtcc-platform/dtcc-core/schemas/model#'


def save_package(city, path):
    city.profile_id, city.profile_version = 'https://example.org/dtcc/profiles/city', '0.2.0'
    city.dataset_context = DatasetContext(
        identity=DatasetIdentity(name=path.stem, title='Synthetic mixed city'),
        metadata=DatasetMetadata(crs=[city.transform.srs]),
        provenance=DatasetProvenance(sources=['local synthetic example']),
        presentation=DatasetPresentation(), request=DatasetRequest(dataset_name=path.stem),
    )
    package = city.export(path, canonical=True)
    restored = load_model_package(package.path)
    assert exchange.dumps(restored) == exchange.dumps(city)
    assert restored.dataset_context == city.dataset_context
    return restored


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    city = io.load_city(HERE / 'fixtures/mixed-city.city.json', strict=True)
    park = city.get_children(Landuse)[0]
    generic = {obj.id: obj for obj in city.get_children(Object)}
    plant, bench = generic['tree-1'], generic['bench-1']
    assert city.buildings[0].footprint().to_polygon().area == 96
    assert park.lod0.surfaces[0].to_polygon().area == 360
    assert bench.id == 'bench-1' and bench.lod1.x == 325024.123
    assert plant.attributes['species'] == 'Quercus robur' and not city.trees
    city.save(args.output / 'source.dtcc')
    restored = io.load_model(args.output / 'source.dtcc')
    assert exchange.dumps(restored) == exchange.dumps(city)
    restored.save(args.output / 'source.city.json', strict=True)
    reimported = io.load_city(args.output / 'source.city.json', strict=True)
    expected = {obj.id: obj for group in city.children.values() for obj in group}
    for group in reimported.children.values():
        for obj in group:
            assert exchange.dumps(obj) == exchange.dumps(expected[obj.id])
    source_records = project(save_package(city, args.output / 'source.dtccpkg'))

    # A separate, explicitly authored native example. This is not an implicit
    # source decoder: plant position, height and radius are supplied intentionally.
    native = city.copy()
    native.children[Object] = [obj for obj in native.get_children(Object) if obj.id != 'tree-1']
    tree = Tree(id='tree-1', semantic_type=NS + 'Tree',
                position=np.array([325016.123456789, 6400006.2356789, 12.125]),
                height=8.123456789, crown_radius=2.123456789,
                attributes={'species': 'Quercus robur'})
    tree.transform.srs = native.transform.srs
    native.add_child(tree)
    native.get_children(Landuse)[0].landuses = [LanduseClasses.GRASS]
    native.save(args.output / 'native.dtcc')
    restored = io.load_model(args.output / 'native.dtcc')
    assert exchange.dumps(restored) == exchange.dumps(native)
    assert restored.trees[0].height == tree.height
    np.testing.assert_array_equal(restored.trees[0].position, tree.position)
    native_records = project(save_package(native, args.output / 'native.dtccpkg'))

    # New semantics fit the same Object and v5 wire without editing core classes.
    extension = native.copy()
    charging_bench = next(obj for obj in extension.get_children(Object) if obj.id == 'bench-1')
    charging_bench.semantic_type = NS + 'ChargingBench'
    charging_bench.attributes['charging_power'] = 100.0
    extension.profile_version = '0.2.1'
    extension = exchange.loads(exchange.dumps(extension))
    data = {'source_objects': source_records, 'native_objects': native_records,
            'extension_objects': project(extension)}
    (args.output / 'profile.json').write_text(json.dumps(data, indent=2) + '\n')
    print('Passed mixed source/native file/package workflows, strict source CityJSON round trip, '
          '96 m2 footprint, 360 m2 park, tree measurements and bench point access')


if __name__ == '__main__':
    main()
