"""Audit unchanged local CityJSON sources and exercise admitted building models.

No downloads, source filtering or automatic format upgrades. Unsupported imports
and unsupported/invalid geometry operations are reported separately. Failed
round-trip assertions fail the command; no source repairs are made.
"""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import platform
import tempfile
import time

import numpy as np

from dtcc_core import io
from dtcc_core.model import exchange
from dtcc_core.builder.meshing.backends import resolve_2d_mesher
from benchmark_buildings import measure


def audit(path, repeats):
    with path.open('rb') as stream:
        raw = stream.read(exchange.MAX_BYTES + 1)
    if len(raw) > exchange.MAX_BYTES:
        raise ValueError('Source exceeds the canonical byte limit')
    source = json.loads(raw)
    objects = source['CityObjects']
    geometries = [g for obj in objects.values() for g in obj.get('geometry', [])]
    report = {
        'file': path.name, 'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest(),
        'version': source['version'], 'vertices': len(source['vertices']),
        'object_types': dict(Counter(obj['type'] for obj in objects.values())),
        'geometry_types_and_lods': dict(Counter(f"{g['type']}@{g.get('lod')}" for g in geometries)),
        'top_level_keys': sorted(source), 'metadata_keys': sorted(source.get('metadata', {})),
        'object_keys': sorted({k for obj in objects.values() for k in obj}),
        'attribute_names': sorted({k for obj in objects.values() for k in obj.get('attributes', {})}),
        'geometry_keys': sorted({k for g in geometries for k in g}),
        'semantic_types': dict(Counter(s['type'] for g in geometries
                                      for s in g.get('semantics', {}).get('surfaces', []))),
        'objects_with_repeated_lods': sum(
            len(gs := obj.get('geometry', [])) != len({g.get('lod') for g in gs})
            for obj in objects.values()),
    }
    start = time.perf_counter()
    try:
        city = io.load_city(path, strict=True)
    except (ValueError, NotImplementedError) as error:
        report['strict_import'] = {'status': 'rejected', 'exception': type(error).__name__, 'reason': str(error)}
        return report
    report['strict_import'] = {'status': 'accepted', 'seconds': time.perf_counter() - start}
    report['canonical_codec'] = measure(city, repeats)
    original = exchange.dumps(city)
    with tempfile.TemporaryDirectory(prefix='dtcc-real-buildings-') as directory:
        directory = Path(directory)
        city.save(directory / 'source.dtcc')
        restored = io.load_model(directory / 'source.dtcc')
        assert exchange.dumps(restored) == original
        try:
            restored.save(directory / 'export.city.json', strict=True)
        except (ValueError, NotImplementedError) as error:
            cityjson_result = {'status': 'rejected', 'reason': str(error)}
        else:
            reimported = io.load_city(directory / 'export.city.json', strict=True)
            assert {b.id for b in restored.buildings} == {b.id for b in reimported.buildings}
            cityjson_result = {'status': 'exported and reimported'}
        meshes, mesh_seconds = [], 0.0
        # This acceptance path deliberately selects the admitted tutorial's LoD2.
        # Other datasets are audited; they are never coerced into this workflow.
        if any(b.building_parts or b.lod2 is None for b in city.buildings):
            raise NotImplementedError('Meshing audit requires top-level LoD2 buildings without parts')
        for index, building in enumerate(city.buildings):
            geometry = building.lod2
            before = exchange.dumps(geometry)
            start = time.perf_counter()
            try:
                mesh = geometry.mesh()
            except (ValueError, NotImplementedError) as error:
                assert exchange.dumps(geometry) == before
                meshes.append({'building_id': building.id, 'polygons': len(geometry.surfaces),
                               'status': 'rejected', 'reason': str(error)})
                continue
            mesh_seconds += time.perf_counter() - start
            assert exchange.dumps(geometry) == before
            assert len(mesh.regions) == len(geometry.regions)
            # All tutorial polygons are triangles: region area gives an
            # independent geometric check on the transferred memberships.
            triangles = mesh.vertices[mesh.faces]
            areas = np.linalg.norm(np.cross(triangles[:, 1] - triangles[:, 0],
                                            triangles[:, 2] - triangles[:, 0]), axis=1) / 2
            for region, original_region in zip(mesh.regions, geometry.regions):
                assert (region.semantic_type, region.id, region.attributes) == (
                    original_region.semantic_type, original_region.id, original_region.attributes)
                source_area = 0.0
                for i in original_region.indices:
                    s = geometry.surfaces[i]
                    assert len(s.vertices) == 3 and not s.holes, 'Area evidence requires source triangles'
                    source_area += np.linalg.norm(np.cross(s.vertices[1] - s.vertices[0],
                                                           s.vertices[2] - s.vertices[0])) / 2
                np.testing.assert_allclose(areas[region.indices].sum(), source_area, rtol=1e-8, atol=1e-6)
            mesh.save(directory / f'mesh-{index}.dtcc')
            assert exchange.dumps(io.load_model(directory / f'mesh-{index}.dtcc')) == exchange.dumps(mesh)
            meshes.append({'building_id': building.id, 'polygons': len(geometry.surfaces),
                           'status': 'passed', 'triangles': len(mesh.faces), 'region_triangle_counts': [len(r.indices) for r in mesh.regions]})
        report['workflow'] = {'canonical_roundtrip': 'exact', 'cityjson_export': cityjson_result,
                              'mesher': resolve_2d_mesher(), 'meshing_seconds_single_run': mesh_seconds,
                              'region_area_checks': 'passed for successfully meshed buildings', 'meshes': meshes}
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('paths', type=Path, nargs='+')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--repeats', type=int, default=3)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error('repeats must be positive')
    result = {'python': platform.python_version(), 'platform': platform.platform(),
              'codec_repeats': args.repeats, 'sources': [audit(p, args.repeats) for p in args.paths]}
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    for source in result['sources']:
        print(source['file'], source['strict_import'])


if __name__ == '__main__':
    main()
