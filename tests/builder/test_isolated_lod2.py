"""Independent LOD2 skins share reconstruction and meshing with joint surfaces."""
import importlib

import numpy as np
import pytest
from shapely.geometry import box

from dtcc_core import builder, datasets
from dtcc_core.builder.meshing.roof_validation import vertex_link_audit
from dtcc_core.builder.meshing.tetgen import has_volume_veto
from dtcc_core.model import GeometryType
from .test_roof_catalogue import building, city_of, model


def sample(*, observations=True):
    footprints = [box(-10, -5, 10, 5), box(7, -5, 27, 5)]
    return city_of([building(p, model('gable', [0., 8.+i, .5], p, key=f'b{i}'),
                            observations=observations) for i, p in enumerate(footprints)], slope=1.)


def field(mesh, name):
    return next(f.values for f in mesh.fields if f.name == name)


def check_skin(mesh, base_z):
    edges = np.concatenate([mesh.faces[:, [0, 1]], mesh.faces[:, [1, 2]], mesh.faces[:, [2, 0]]])
    keys, counts = np.unique(np.sort(edges, axis=1), axis=0, return_counts=True)
    assert counts.max() == 2
    assert np.any(counts == 1)
    np.testing.assert_allclose(mesh.vertices[keys[counts == 1], 2], base_z, atol=1e-8)
    assert vertex_link_audit(mesh)['nonmanifold_vertices'] == 0


@pytest.mark.parametrize('method', ['catalogue', 'ransac'])
def test_overlapping_buildings_keep_independent_fitted_skins(method):
    city = sample()
    original = [b.point_cloud.points.copy() for b in city.buildings]
    parts = city.build_surface_mesh(lod='lod2', lod2_method=method, isolated=True,
        merge_meshes=False, mesher='dtcc_mesher', report_mesh_quality=False)
    assert len(parts) == 3
    terrain, *skins = parts
    audit = terrain.stage_audit['roof_catalogue']
    assert audit['retained_nonflat_fits'] == 2
    assert not audit['block_fallbacks']
    assert np.all(terrain.markers == -2)
    for skin, local in zip(skins, audit['building_assemblies']):
        check_skin(skin, local['base_z'])
        assert local['roof_catalogue']['compact_surface']['topology']['valid']
        assert has_volume_veto(skin)
        assert np.all(field(skin, 'building') == local['building'])
        patches = [p for p in audit['patch_markers'] if p['block_id'] == local['building']]
        assert sum(p['area'] for p in patches) == pytest.approx(200.)
        for p in patches:
            xyz = skin.vertices[skin.faces[skin.markers == p['roof_marker']]].reshape(-1, 3)
            assert len(xyz)
            predicted = (xyz[:, :2]-p['origin']) @ np.asarray(p['coefficients'])[:2] + p['coefficients'][2]
            np.testing.assert_allclose(xyz[:, 2], predicted, atol=1e-8)
    assert has_volume_veto(terrain)
    for b, before in zip(city.buildings, original):
        np.testing.assert_array_equal(b.point_cloud.points, before)
        assert b.lod2 is None


def test_combining_skins_does_not_weld_terrain_or_overlapping_buildings():
    kwargs = dict(lod='lod2', isolated=True, mesher='dtcc_mesher', report_mesh_quality=False)
    parts = builder.build_city_surface_mesh(sample(), merge_meshes=False, **kwargs)
    merged = builder.build_city_surface_mesh(sample(), **kwargs)
    assert len(merged.vertices) == sum(len(p.vertices) for p in parts)
    np.testing.assert_array_equal(merged.vertices, np.vstack([p.vertices for p in parts]))
    owners = field(merged, 'building')
    vertex_sets = [set(merged.faces[owners == owner].ravel()) for owner in np.unique(owners)]
    assert all(not a & b for i, a in enumerate(vertex_sets) for b in vertex_sets[i+1:])
    assert has_volume_veto(merged)


def test_reuses_isolated_cleaning_and_original_observations(monkeypatch):
    city = sample().clean_footprints(isolated=True)
    orchestration = importlib.import_module('dtcc_core.builder.meshing.catalogue_surface')
    def unexpected(*args, **kwargs):
        raise AssertionError('Cleaning must not run twice')
    monkeypatch.setattr(orchestration, 'clean_regions', unexpected)
    mesh = city.build_surface_mesh(lod='lod2', isolated=True, mesher='dtcc_mesher', report_mesh_quality=False)
    audit = mesh.stage_audit['roof_catalogue']
    assert audit['reused_cleaning']
    assert audit['retained_nonflat_fits'] == 2


@pytest.mark.parametrize('method', ['catalogue', 'ransac'])
def test_missing_points_preserves_lod1_height_and_records_fallback(method):
    mesh = builder.build_city_surface_mesh(sample(observations=False), lod='lod2',
        lod2_method=method, isolated=True, mesher='dtcc_mesher', report_mesh_quality=False)
    audit = mesh.stage_audit['roof_catalogue']
    assert all(s['final_outcome'] == 'lod1_fallback' for s in audit['sources'])
    for patch in audit['patch_markers']:
        assert patch['coefficients'] == pytest.approx([0., 0., 10.])
        np.testing.assert_allclose(mesh.vertices[mesh.faces[mesh.markers == patch['roof_marker']], 2], 10.)


def test_dataset_build_from_city_supports_isolated_lod2():
    mesh = datasets.city_surface_mesh.build_from_city(sample(), bounds=(-25, -20, 45, 25), lod='lod2', isolated=True,
        stage_audit_enabled=True, mesher='dtcc_mesher', report_mesh_quality=False)
    assert mesh.stage_audit['isolated']
    assert mesh.stage_audit['roof_catalogue']['retained_nonflat_fits'] == 2


def test_empty_city_has_only_independent_terrain():
    mesh = city_of([]).build_surface_mesh(lod='lod2', isolated=True,
        mesher='dtcc_mesher', report_mesh_quality=False)
    assert len(mesh.faces) and np.all(mesh.markers == -2)
    assert np.all(field(mesh, 'building') == -1)
    assert mesh.stage_audit['roof_catalogue']['retained_fits'] == 0
    assert has_volume_veto(mesh)


def test_cleaned_group_keeps_source_fits_in_one_isolated_shell():
    city = sample().clean_footprints()
    assert len(city.buildings) == 1 and len(city.buildings[0].building_parts) == 2
    parts = city.build_surface_mesh(lod='lod2', isolated=True, merge_meshes=False,
        mesher='dtcc_mesher', report_mesh_quality=False)
    assert len(parts) == 2
    audit = parts[0].stage_audit['roof_catalogue']
    assert audit['retained_nonflat_fits'] == 2
    assert len(audit['building_assemblies']) == 1
    check_skin(parts[1], audit['building_assemblies'][0]['base_z'])


def test_outside_and_small_buildings_are_excluded_without_affecting_valid_skin():
    city = sample()
    for key, footprint in [('outside', box(44, 0, 54, 10)), ('tiny', box(30, 0, 31, 1))]:
        city.add_buildings([building(footprint, model('flat', [8.], footprint, key=key))])
    mesh = city.build_surface_mesh(lod='lod2', isolated=True, mesher='dtcc_mesher', report_mesh_quality=False)
    sources = mesh.stage_audit['roof_catalogue']['sources']
    assert {r['source_key'] for r in sources} == {'b0', 'b1'}


def test_coplanar_patch_merging_keeps_global_marker_mapping(monkeypatch):
    from dtcc_core.builder.roofs import RoofFitResult
    from dtcc_core.builder.roofs.ransac import RansacRoofModel
    adapter = importlib.import_module('dtcc_core.builder.roofs.ransac')
    def fit(footprint, points, ground, *, source_key):
        xmin, ymin, xmax, ymax = footprint.bounds
        cuts = np.linspace(xmin, xmax, 4)
        domains = tuple(box(a, ymin, b, ymax) for a, b in zip(cuts[:-1], cuts[1:]))
        planes = np.tile([0., 0., 8.], (3, 1))
        roof = RansacRoofModel(source_key, 'flat', planes.ravel(), np.zeros(2), planes, {}, domains)
        return RoofFitResult(roof, dict(accepted=True, method='ransac'))
    monkeypatch.setattr(adapter, 'fit_ransac', fit)
    mesh = sample().build_surface_mesh(lod='lod2', lod2_method='ransac', isolated=True,
        mesher='dtcc_mesher', report_mesh_quality=False)
    audit = mesh.stage_audit['roof_catalogue']
    total = len(audit['patch_markers'])
    assert total == 6
    xyz = mesh.vertices[mesh.faces]
    cross = np.cross(xyz[:, 1]-xyz[:, 0], xyz[:, 2]-xyz[:, 0])
    roofs = (field(mesh, 'building') >= 0) & (cross[:, 2] > 1e-8)
    assert set(mesh.markers[roofs]) == {6, 9}


def test_skin_validation_rejects_an_open_roof_edge():
    from dtcc_core.builder.meshing.roof_validation import validate_shell
    from dtcc_core.builder.geometry_builders.meshes import _shell_closure_audit
    parts = sample().build_surface_mesh(lod='lod2', isolated=True, merge_meshes=False,
        mesher='dtcc_mesher', report_mesh_quality=False)
    mesh = parts[1]
    audit = mesh.stage_audit['roof_catalogue']
    base = audit['building_assemblies'][0]['base_z']
    roof_marker = audit['patch_markers'][0]['roof_marker']
    remove = np.flatnonzero(mesh.markers == roof_marker)[0]
    mesh.faces = np.delete(mesh.faces, remove, axis=0)
    mesh.markers = np.delete(mesh.markers, remove)
    report = validate_shell(mesh, num_patches=len(audit['patch_markers']),
                            closure_audit=_shell_closure_audit, base_z=base)
    assert not report['valid']
    assert report['closure']['unexpected_open_edges'] > 0
