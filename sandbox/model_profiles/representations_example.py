"""Checksum-pinned 3DBAG -> edit -> native file -> mesh -> canonical package.

Run in the ordinary DTCC environment. No network, source filtering or geometry
repair. Extent discrepancies are explicitly recomputed into Dataset Context.
Standard schema validation stays enabled. Strict CityJSON exchange is checked
before attaching the derived mesh, which that adapter cannot represent.
"""

import argparse
import hashlib
import json
from pathlib import Path
import platform
import resource
import sys
import time

import numpy as np
from pyproj import CRS

from dtcc_core import io
from dtcc_core.datasets import load_model_package
from dtcc_core.model import Field, Solid, exchange

SOURCE_SHA256 = '2bb5d22ae2cbfe2096041e3a79b3e826f43d3a8c9824eeb019feab4e0a2742ba'
SOURCE_URL = 'https://3d.bk.tudelft.nl/opendata/cityjson/3dcities/v2.0/9-284-556.city.json'
BUILDING_ID = 'NL.IMBAG.Pand.0503100000000030'
SECOND_PART_ID = 'NL.IMBAG.Pand.0503100000000010-0'


def objects(city):
    result, stack = {}, [city]
    while stack:
        obj = stack.pop()
        result[obj.id] = obj
        stack.extend(child for group in obj.children.values() for child in group)
    return result


def compare(source, restored):
    restored.id = source.id  # Native aggregate UUID has no CityJSON field.
    a, b = objects(source), objects(restored)
    assert set(a) == set(b)
    maximum_error = 0.0
    for id, old in a.items():
        new = b[id]
        assert type(old) is type(new) and old.attributes == new.attributes
        assert old.transform.srs == new.transform.srs
        assert {child.id for group in old.children.values() for child in group} == {
            child.id for group in new.children.values() for child in group}
        assert list(old.geometry) == list(new.geometry)
        for key, record in old.geometry.items():
            other = new.geometry[key]
            assert (record.lod, record.role) == (other.lod, other.role)
            x, y = record.geometry, other.geometry
            assert type(x) is type(y) and len(x.surfaces) == len(y.surfaces)
            assert x.transform.srs == y.transform.srs
            if isinstance(x, Solid):
                assert len(x.shells) == len(y.shells)
                for u, v in zip(x.shells, y.shells): np.testing.assert_array_equal(u, v)
            for u, v in zip(x.surfaces, y.surfaces):
                assert len(u.holes) == len(v.holes)
                for ur, vr in zip([u.vertices, *u.holes], [v.vertices, *v.holes]):
                    assert ur.shape == vr.shape
                    error = float(np.max(np.abs(ur-vr)))
                    maximum_error = max(maximum_error, error)
                    np.testing.assert_allclose(ur, vr, atol=0.00050001, rtol=0)
            assert len(x.regions) == len(y.regions)
            for u, v in zip(x.regions, y.regions):
                assert (u.id,u.semantic_type,u.attributes,u.parent) == (v.id,v.semantic_type,v.attributes,v.parent)
                np.testing.assert_array_equal(u.indices, v.indices)
    return maximum_error


def polygon_area(surface):
    def ring_area(ring):
        local = ring - ring[0]
        return np.linalg.norm(np.cross(local, np.roll(local, -1, axis=0)).sum(axis=0)) / 2
    return ring_area(surface.vertices) - sum(ring_area(hole) for hole in surface.holes)


def mesh_part(part):
    solid = part.get_geometry(lod='2.2')
    before = exchange.dumps(solid)
    start = time.perf_counter()
    mesh = solid.mesh(mesher='dtcc_mesher')
    mesh_seconds = time.perf_counter() - start
    assert exchange.dumps(solid) == before
    np.testing.assert_array_equal(mesh.transform.affine, solid.transform.affine)
    assert mesh.transform.srs == solid.transform.srs
    # This source uses metre axes and an identity affine: areas below are m².
    assert all(axis.unit_name == 'metre' for axis in CRS(solid.transform.srs).axis_info)
    np.testing.assert_array_equal(solid.transform.affine, np.eye(4))
    triangles = mesh.vertices[mesh.faces]
    areas = np.linalg.norm(np.cross(triangles[:, 1] - triangles[:, 0],
                                    triangles[:, 2] - triangles[:, 0]), axis=1) / 2
    assert np.all(np.isfinite(areas)) and np.all(areas > 0)
    source_areas = np.array([polygon_area(surface) for surface in solid.surfaces])
    np.testing.assert_allclose(areas.sum(), source_areas.sum(), rtol=1e-6, atol=1e-6)
    assert len(mesh.regions) == len(solid.regions)
    for old, new in zip(solid.regions, mesh.regions):
        assert (old.id, old.semantic_type, old.attributes, old.parent) == (
            new.id, new.semantic_type, new.attributes, new.parent)
        np.testing.assert_allclose(areas[new.indices].sum(), source_areas[old.indices].sum(),
                                   rtol=1e-6, atol=1e-6)
    mesh.fields = [Field(name='triangle_area', unit='m2', association='face', values=areas)]
    part.add_geometry(mesh, id='boundary_mesh', lod='2.2', role='boundary_mesh')
    return {'part_id': part.id, 'source_surfaces': len(solid.surfaces),
            'source_shells': len(solid.shells), 'triangles': len(mesh.faces),
            'region_triangle_counts': [len(region.indices) for region in mesh.regions],
            'area_m2': float(areas.sum()), 'field': 'triangle_area: float64, face, m2',
            'triangulation_seconds_single_run': mesh_seconds,
            'mesher': 'dtcc_mesher', 'region_area_relative_tolerance': 1e-6}


def main():
    from benchmark_buildings import measure

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    with args.source.open('rb') as stream:
        raw = stream.read()
    if hashlib.sha256(raw).hexdigest() != SOURCE_SHA256:
        raise ValueError('This evidence command requires the recorded unchanged 3DBAG source')
    args.output.mkdir(parents=True, exist_ok=True)
    seconds = {}

    def timed(name, operation):
        start = time.perf_counter()
        result = operation()
        seconds[name] = time.perf_counter() - start
        print(f'{name}: {seconds[name]:.3f} s', flush=True)
        return result

    city = timed('strict_import', lambda: io.load_city(args.source, strict=True, extent_policy='recompute'))
    city.dataset_context.provenance.sources = [{'url': SOURCE_URL, 'sha256': SOURCE_SHA256}]
    native = objects(city)
    assert len(city.buildings) == 1110 and len(native) == 2222
    records = [r for o in native.values() for r in o.geometry.values()]
    assert len(records) == 4443 and sum(type(r.geometry) is Solid for r in records) == 3333
    surface_count = sum(len(record.geometry.surfaces) for record in records)
    collapsed = sum(len(np.unique(ring, axis=0)) < 3 for record in records
                    for surface in record.geometry.surfaces for ring in [surface.vertices, *surface.holes])
    del records
    building = native[BUILDING_ID]
    part = building.building_parts[0]
    assert part.lod2 is None and type(part.get_geometry(lod='2.2')) is Solid
    edit = {'building_id': building.id, 'attribute': 'description',
            'was_present': 'description' in building.attributes,
            'previous_value': building.attributes.get('description'),
            'value': 'DTCC real-data integration example'}
    building.attributes['description'] = edit['value']
    # Check external interchange before adding content without a CityJSON mapping.
    timed('cityjson_save', lambda: city.save(args.output / 'export.city.json', strict=True))
    restored = timed('cityjson_load', lambda: io.load_city(args.output / 'export.city.json', strict=True))
    maximum_error = compare(city, restored)
    del restored
    payload = exchange.dumps(city)
    timed('native_save_before_mesh', lambda: city.save(args.output / 'imported.dtcc'))
    restored = timed('native_load_before_mesh', lambda: io.load_model(args.output / 'imported.dtcc'))
    assert exchange.dumps(restored) == payload
    assert (restored.schema_id, restored.schema_version) == (city.schema_id, city.schema_version)
    assert restored.dataset_context is None
    # Context is a package-manifest concern; carry the known source context into
    # the package explicitly after loading the standalone native model.
    restored.dataset_context = city.dataset_context
    city = restored
    native = objects(city)
    building = native[BUILDING_ID]
    part = building.building_parts[0]
    assert building.attributes['description'] == edit['value']
    mesh_report = timed('mesh_and_verify_selected_part', lambda: mesh_part(part))
    # This first-in-source part failed before the polygon-normal fix. Check it
    # explicitly as well; successful triangulation is not validity certification.
    second_solid = native[SECOND_PART_ID].get_geometry(lod='2.2')
    before = exchange.dumps(second_solid)
    result = second_solid.mesh(mesher='dtcc_mesher')
    source_probe = {'part_id': SECOND_PART_ID, 'outcome': 'meshed', 'triangles': len(result.faces)}
    del result
    assert exchange.dumps(second_solid) == before
    payload = exchange.dumps(city)
    path = args.output / 'model.dtcc'
    timed('native_save_with_mesh', lambda: city.save(path))
    restored = timed('native_load_with_mesh', lambda: io.load_model(path))
    assert exchange.dumps(restored) == payload
    del restored
    # Prove an invalid edit fails at the ordinary save boundary, before replacing
    # the last valid file. Avoid inventing or retaining any measured height.
    attributes = building.attributes.copy()
    building.attributes['storeys_above_ground'] = True
    try:
        city.save(path)
    except ValueError as exc:
        assert 'storeys_above_ground' in str(exc)
        schema_failure = str(exc)
    else:
        raise AssertionError('Invalid storey count unexpectedly saved')
    finally:
        building.attributes = attributes
    assert path.read_bytes() == payload
    package = timed('package_write', lambda: city.export(args.output / 'package.dtccpkg', canonical=True))
    restored = timed('package_read', lambda: load_model_package(package.path))
    assert exchange.dumps(restored) == payload
    assert restored.dataset_context == city.dataset_context
    del restored
    codec = measure(city, 3)
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform != 'darwin': peak *= 1024
    report = {'source_url':SOURCE_URL, 'sha256':SOURCE_SHA256, 'source_bytes':len(raw),
              'python':platform.python_version(), 'platform':platform.platform(), 'wire_version':exchange.VERSION,
              'schema_id': city.schema_id, 'schema_version': city.schema_version,
              'default_schema_validation': True,
              'buildings':1110, 'parts':1111, 'representations':4443, 'solids':3333,
              'surfaces':surface_count,
              'source_collapsed_rings_preserved':collapsed,
              'extent_discrepancies_retained':city.dataset_context.health['extent_discrepancy_count'],
              'maximum_cityjson_coordinate_error':maximum_error, 'public_operation_seconds_single_run':seconds,
              'edit': edit, 'mesh': mesh_report, 'additional_source_meshing_probe': source_probe,
              'invalid_edit_rejected_without_replacing_file': schema_failure,
              'final_representations': sum(len(obj.geometry) for obj in native.values()),
              'canonical_file_and_package':'exact; package also retains Dataset Context',
              'cityjson_roundtrip':'all feature IDs, attributes, representations, shells and regions checked',
              'codec_three_run_medians':codec, 'workflow_process_peak_MiB':peak/2**20,
              'limits':'Peak includes dependencies, multiple models and evidence comparisons; not isolated model size or codec overhead. Public timings are single observations; mesh timing includes verification. Codec medians are three warm runs with validation. No source repair, geometric validity certification or full-tile meshing. Derived mesh has no shell partition/volume cells; original Solid remains. CityJSON check precedes mesh attachment; Dataset Context is restored explicitly from the import for package export.'}
    (args.output / 'report.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
