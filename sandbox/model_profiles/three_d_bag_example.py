"""Verify the pinned real 3DBAG semantic mapping and ordinary exchange boundaries.

Run from the repository: python -m sandbox.model_profiles.three_d_bag_example SOURCE OUTPUT
The source tile is not downloaded, repaired or modified by this example.
"""

import argparse
import gc
import hashlib
import json
from pathlib import Path
from statistics import median
from time import perf_counter

from dtcc_core import io
from dtcc_core.datasets import load_model_package
from dtcc_core.model import exchange
from dtcc_core.model._standard_schema import validate_admitted
from .representations_example import compare

SOURCE_URL = 'https://3d.bk.tudelft.nl/opendata/cityjson/3dcities/v2.0/9-284-556.city.json'
SHA256 = '2bb5d22ae2cbfe2096041e3a79b3e826f43d3a8c9824eeb019feab4e0a2742ba'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    if hashlib.sha256(args.source.read_bytes()).hexdigest() != SHA256:
        raise ValueError('This example requires the documented, unchanged 3DBAG tile')
    args.output.mkdir(parents=True, exist_ok=True)
    baseline = io.load_city(args.source, strict=True, extent_policy='recompute')
    start = perf_counter()
    city = io.load_3dbag(args.source, extent_policy='recompute')
    loader_seconds = perf_counter() - start  # Schema is already initialized.
    city.id = baseline.id
    assert len(city.buildings) == 1110
    assert sum(len(b.attributes['elevation_measurements']) for b in city.buildings) == 5550
    assert all(b.height is None for b in city.buildings)
    city.dataset_context.provenance.sources.append({'url': SOURCE_URL, 'sha256': SHA256})
    payload = exchange.dumps(city)
    # Removing only enrichment must recover every native source fact, including
    # representation order, exact coordinates, shells, regions and source attrs.
    stripped = exchange.loads(payload)
    for b in stripped.buildings:
        del b.attributes['elevation_measurements']
        del b.attributes['roof_type']
    assert exchange.dumps(stripped) == exchange.dumps(baseline)
    del stripped
    gc.collect()
    timings = {'unmapped_schema': [], 'mapped_schema': [], 'encode': [], 'decode': []}
    for _ in range(5):
        for label, model in [('unmapped_schema', baseline), ('mapped_schema', city)]:
            start = perf_counter()
            validate_admitted(model, model.schema_id, model.schema_version)
            timings[label].append(perf_counter() - start)
    del baseline
    gc.collect()
    for _ in range(3):
        start = perf_counter()
        encoded = exchange.dumps(city)
        timings['encode'].append(perf_counter() - start)
        assert encoded == payload
        start = perf_counter()
        restored = exchange.loads(encoded)
        timings['decode'].append(perf_counter() - start)
        del restored
        gc.collect()
    path = args.output / 'mapped.dtcc'
    city.save(path)
    assert exchange.dumps(io.load_city(path)) == payload
    package = city.export(args.output / 'mapped.dtccpkg', canonical=True)
    restored = load_model_package(package.path)
    assert exchange.dumps(restored) == payload
    assert restored.dataset_context == city.dataset_context
    del restored
    city.save(args.output / 'mapped.city.json', strict=True)
    restored = io.load_city(args.output / 'mapped.city.json', strict=True)
    maximum_cityjson_error = compare(city, restored)
    del restored
    # Default validation rejects loss of the vertical reference atomically.
    record = city.buildings[0].attributes['elevation_measurements'][0]
    vertical_reference = record.pop('vertical_reference')
    try:
        city.save(path)
    except ValueError as error:
        assert 'vertical_reference' in str(error)
    else:
        raise AssertionError('Missing elevation reference was accepted')
    assert path.read_bytes() == payload
    record['vertical_reference'] = vertical_reference
    report = {
        'source': SOURCE_URL, 'sha256': SHA256, 'source_release': None,
        'schema_version': city.schema_version, 'wire_version': exchange.VERSION,
        'buildings': 1110, 'elevation_records': 5550, 'native_bytes': len(payload),
        'maximum_cityjson_coordinate_error_m': maximum_cityjson_error,
        'loader_seconds_schema_warm': loader_seconds,
        'warm_median_seconds': {key: median(values) for key, values in timings.items()},
        'timing_samples_seconds': timings,
        'example': city.buildings[0].attributes['elevation_measurements'][0],
        'verification': 'Unchanged source facts; exact native/package and package context; strict CityJSON metadata, shells and regions exact, coordinates within 0.00050001 m; rejected save preserves file.',
    }
    (args.output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
