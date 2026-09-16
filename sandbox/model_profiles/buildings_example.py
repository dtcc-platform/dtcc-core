"""Exercise the buildings workflow with DTCC's ordinary Python file/package APIs."""

import argparse
import json
from pathlib import Path

import numpy as np

from dtcc_core import io
from dtcc_core.datasets import load_model_package
from dtcc_core.datasets.schema import (DatasetContext, DatasetIdentity, DatasetMetadata,
                                      DatasetPresentation, DatasetProvenance, DatasetRequest)
from dtcc_core.model import exchange
from dtcc_core.model.profiles import _project

HERE = Path(__file__).resolve().parent
NAMESPACE = 'https://github.com/dtcc-platform/dtcc-core/schemas/model#'


def project(city):
    """Use the shared native projection for the isolated historical evaluators."""
    exchange.validate(city)
    records = _project(city)[0]
    # The optional evaluators deliberately exercise immutable historical schemas.
    # Reclassify only this disposable experiment projection, never native data.
    for record in records:
        if record['semantic_type'].startswith(NAMESPACE):
            record['semantic_type'] = 'https://example.org/dtcc/' + record['semantic_type'][len(NAMESPACE):]
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    city = io.load_city(HERE / 'fixtures/buildings.city.json', strict=True)
    building = city.buildings[0]
    footprint = building.footprint().to_polygon()
    assert footprint.area == 96.0 and len(footprint.interiors) == 1
    part = building.building_parts[0]
    roof = part.lod2.regions_of(NAMESPACE + 'RoofSurface')[0]
    assert part.id == 'part-1' and roof.indices.tolist() == [1, 2]
    assert roof.attributes['solar-potential'] == 42.5
    original = exchange.dumps(city)
    city.save(args.output.with_suffix('.dtcc'))
    restored = io.load_model(args.output.with_suffix('.dtcc'))
    assert exchange.dumps(restored) == original
    restored.save(args.output.with_suffix('.city.json'), strict=True)
    reimported = io.load_city(args.output.with_suffix('.city.json'), strict=True)
    # City is a native aggregate, not a CityJSON city object; its UUID is not a
    # source feature ID. Feature IDs and semantic assignments are preserved.
    assert reimported.buildings[0].id == building.id
    result_part = reimported.buildings[0].building_parts[0]
    for a, b in zip(part.lod2.surfaces, result_part.lod2.surfaces):
        np.testing.assert_allclose(a.vertices, b.vertices, atol=0.0005, rtol=0)
    assert result_part.lod2.regions[1].indices.tolist() == [1, 2]
    # Profile selection is explicit. Canonical exchange can carry its identity;
    # plain CityJSON has no mapping for it in this bounded adapter.
    city.profile_id = 'https://example.org/dtcc/profiles/buildings'
    city.profile_version = '0.3.0'
    city.dataset_context = DatasetContext(
        identity=DatasetIdentity(name='buildings-example', title='Synthetic buildings'),
        metadata=DatasetMetadata(crs=[city.transform.srs]),
        provenance=DatasetProvenance(sources=['local synthetic CityJSON fixture']),
        presentation=DatasetPresentation(), request=DatasetRequest(dataset_name='buildings-example'),
    )
    package = city.export(args.output.with_suffix('.dtccpkg'), canonical=True)
    packaged = load_model_package(package.path)
    assert exchange.dumps(packaged) == exchange.dumps(city)
    assert packaged.dataset_context == city.dataset_context
    # A schema-only subtype remains the same native SemanticRegion and wire shape.
    future = city.copy()
    future.profile_version = '0.3.1'
    future.buildings[0].building_parts[0].lod2.regions[1].semantic_type = NAMESPACE + 'SolarRoofSurface'
    future.buildings[0].building_parts[0].lod2.regions[1].attributes['efficiency'] = 0.22
    future = exchange.loads(exchange.dumps(future))
    args.output.write_text(json.dumps({'objects': project(packaged), 'extension_objects': project(future)}, indent=2)+'\n')
    print(f'Passed native file, package and CityJSON workflow; footprint={footprint.area} m2, roof polygons={len(roof.indices)}')


if __name__ == '__main__':
    main()
