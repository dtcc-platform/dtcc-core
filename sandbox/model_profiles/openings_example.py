"""Load opening surfaces, access their host and geometry, save and triangulate."""

import argparse
import json
from pathlib import Path

import numpy as np

from dtcc_core import io
from dtcc_core.datasets import load_model_package
from dtcc_core.datasets.schema import (DatasetContext, DatasetIdentity, DatasetMetadata,
                                      DatasetPresentation, DatasetProvenance, DatasetRequest)
from dtcc_core.model import exchange

from buildings_example import project

HERE = Path(__file__).resolve().parent
NS = 'https://github.com/dtcc-platform/dtcc-core/schemas/model#'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    city = io.load_city(HERE / 'fixtures/openings.city.json', strict=True)
    geometry = city.buildings[0].building_parts[0].lod3
    window = geometry.regions_of(NS + 'Window')[0]
    host = geometry.regions[window.parent]
    polygons = [geometry.surfaces[index] for index in window.indices]
    assert host.attributes['name'] == 'front-wall' and len(polygons) == 1
    assert polygons[0].vertices.shape == (4, 3)
    city.save(args.output / 'openings.dtcc')
    restored = io.load_model(args.output / 'openings.dtcc')
    assert exchange.dumps(restored) == exchange.dumps(city)
    restored.save(args.output / 'openings.city.json', strict=True)
    reimported = io.load_city(args.output / 'openings.city.json', strict=True)
    assert exchange.dumps(reimported.buildings[0]) == exchange.dumps(city.buildings[0])
    city.profile_id, city.profile_version = 'https://example.org/dtcc/profiles/buildings', '0.3.0'
    city.dataset_context = DatasetContext(
        identity=DatasetIdentity(name='openings-example', title='Synthetic building openings'),
        metadata=DatasetMetadata(crs=[city.transform.srs]),
        provenance=DatasetProvenance(sources=['local synthetic CityJSON fixture']),
        presentation=DatasetPresentation(), request=DatasetRequest(dataset_name='openings-example'),
    )
    package = city.export(args.output / 'openings.dtccpkg', canonical=True)
    packaged = load_model_package(package.path)
    assert exchange.dumps(packaged) == exchange.dumps(city)
    assert packaged.dataset_context == city.dataset_context
    mesh = geometry.mesh()
    mesh.save(args.output / 'openings-mesh.dtcc')
    assert exchange.dumps(io.load_model(args.output / 'openings-mesh.dtcc')) == exchange.dumps(mesh)
    assert [r.parent for r in mesh.regions] == [r.parent for r in geometry.regions]
    triangles = mesh.vertices[mesh.faces]
    areas = np.linalg.norm(np.cross(triangles[:, 1] - triangles[:, 0],
                                    triangles[:, 2] - triangles[:, 0]), axis=1) / 2
    np.testing.assert_allclose([areas[r.indices].sum() for r in mesh.regions[:3]], [50, 4, 6])
    (args.output / 'profile.json').write_text(json.dumps(project(packaged), indent=2) + '\n')
    print(f'Passed file/package/CityJSON/mesh workflow: {window.attributes["name"]} on '
          f'{host.attributes["name"]}; wall/window/door areas 50/4/6 m2')


if __name__ == '__main__':
    main()
