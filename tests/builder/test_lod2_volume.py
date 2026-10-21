"""Joint LOD2 shells survive domain closure and tetrahedralization."""
import numpy as np
import pytest
from shapely.geometry import box

from dtcc_core import builder
from dtcc_core.model import GeometryType
from dtcc_core.builder.geometry_builders import meshes as pipeline
from dtcc_core.builder.meshing.tetgen import is_tetgen_available, has_volume_veto
from dtcc_core.builder.meshing import tetgen_utils
from dtcc_core.datasets.city_volume_mesh import CityVolumeMeshDataset, CityVolumeMeshArgs
from .test_roof_catalogue import building, city_of, model


def sample(group=False):
    p = box(-10, -5, 10, 5)
    buildings = [building(p, model('hip', [0., 8., .6, .6], p, key='hip'))]
    if group:
        q = box(10.2, -5, 20, 5)
        buildings.append(building(q, model('gable', [0., 11., .4], q, key='gable')))
    return city_of(buildings)


def triangle_areas(vertices, faces):
    xyz = vertices[faces]
    return np.linalg.norm(np.cross(xyz[:, 1]-xyz[:, 0], xyz[:, 2]-xyz[:, 0]), axis=1)/2



def assert_roof_plane(points,patch):
    """Use physical normal distance, consistently for every roof orientation."""
    gradient=np.asarray(patch['coefficients'])[:2]
    expected=(points[:,:2]-patch['origin'])@gradient+patch['coefficients'][2]
    residual=(points[:,2]-expected)/np.sqrt(1+gradient@gradient)
    np.testing.assert_allclose(residual,0.,atol=1e-7,rtol=0)


def test_roof_plane_contract_uses_normal_distance_for_steep_planes():
    patch=dict(origin=np.zeros(2),coefficients=np.array([10.,0.,0.]))
    normal=np.array([-10.,0.,1.])/np.sqrt(101.)
    assert_roof_plane(np.array([normal*5e-8]),patch)
    with pytest.raises(AssertionError):assert_roof_plane(np.array([normal*2e-7]),patch)


def assert_volume_boundary(volume, shell, plc, catalogue):
    assert np.isfinite(volume.vertices).all()
    xyz = volume.vertices[volume.cells]
    det = np.linalg.det(xyz[:, 1:]-xyz[:, :1])
    assert np.all(det > 0)
    boundary = volume.boundary_faces
    cells_faces = np.concatenate([volume.cells[:, ids] for ids in
                                  ([0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3])])
    faces, counts = np.unique(np.sort(cells_faces, axis=1), axis=0, return_counts=True)
    assert counts.max() == 2
    assert set(map(tuple, faces[counts == 1])) == set(map(tuple, np.sort(boundary, axis=1)))
    edges = np.sort(boundary[:, [[0, 1], [1, 2], [2, 0]]].reshape(-1, 2), axis=1)
    assert np.all(np.unique(edges, axis=0, return_counts=True)[1] == 2)
    assert set(range(-6, 0)) <= set(volume.boundary_markers)
    original_area = triangle_areas(shell.vertices, shell.faces)
    final_area = triangle_areas(volume.vertices, boundary)
    for marker in np.unique(shell.markers):
        target = -1 if marker < 0 else marker
        assert final_area[volume.boundary_markers == target].sum() == pytest.approx(
            original_area[shell.markers == marker].sum(), rel=1e-7, abs=1e-7)
    for patch in catalogue['patch_markers']:
        points = volume.vertices[boundary[volume.boundary_markers == patch['roof_marker']]].reshape(-1, 3)
        if len(points):
            assert_roof_plane(points,patch)
    # The tetrahedra fill exactly the closed PLC outside the buildings.
    all_faces = np.vstack([plc.shell_faces, plc.boundary_facets])
    points = plc.vertices - plc.vertices.mean(axis=0)
    t = points[all_faces]
    expected_volume = abs(np.einsum('ij,ij->i', t[:, 0], np.cross(t[:, 1], t[:, 2])).sum()/6)
    assert det.sum()/6 == pytest.approx(expected_volume, rel=1e-7)


@pytest.mark.skipif(not is_tetgen_available(), reason='TetGen is not installed')
@pytest.mark.parametrize('group', [False, True])
@pytest.mark.parametrize('method', ['catalogue', 'ransac'])
def test_lod2_volume_reuses_joint_shell_and_preserves_boundary(monkeypatch, group, method):
    city = sample(group)
    surface = builder.build_city_surface_mesh(city, lod='lod2', lod2_method=method,
                                              max_mesh_size=8., report_mesh_quality=False)
    original = pipeline.tetgen_build_volume_mesh
    captured = {}
    def record(**kwargs):
        captured.update(kwargs)
        return original(**kwargs)
    monkeypatch.setattr(pipeline, 'tetgen_build_volume_mesh', record)
    audit = {}
    volume = builder.build_city_volume_mesh(city, lod='lod2', lod2_method=method,
                                            max_mesh_size=8., domain_height=40.,
                                            report_mesh_quality=False, stage_audit=audit)
    for key in ('vertices', 'faces', 'markers'):
        np.testing.assert_array_equal(getattr(surface, key), getattr(captured['mesh'], key))
    assert not has_volume_veto(captured['mesh'])
    cat = audit['roof_catalogue']
    assert cat['retained_nonflat_fits'] == surface.stage_audit['roof_catalogue']['retained_nonflat_fits']
    if method == 'catalogue' or group:
        assert cat['retained_nonflat_fits'] > 0
    assert not cat['block_fallbacks']
    assert len(cat['blocks']) == 1
    assert_volume_boundary(volume, surface, captured['prebuilt_plc'], cat)
    assert audit['attempts'][-1]['result']['status'] == 'success'


@pytest.mark.skipif(not is_tetgen_available(), reason='TetGen is not installed')
@pytest.mark.parametrize('entry', ['builder', 'city', 'dataset'])
def test_volume_lod1_remains_default(entry):
    city = sample()
    kwargs = dict(max_mesh_size=8., domain_height=40.)
    if entry == 'builder':
        call = lambda **kw: builder.build_city_volume_mesh(city, report_mesh_quality=False, **kw)
    elif entry == 'city':
        call = city.build_volume_mesh
    else:
        call = lambda **kw: CityVolumeMeshDataset().build_from_city(
            city, bounds=(-20, -20, 25, 20), report_mesh_quality=False, **kw)
    default, explicit = call(**kwargs), call(**kwargs, lod='lod1')
    for key in ('vertices', 'cells', 'boundary_faces', 'boundary_markers'):
        np.testing.assert_array_equal(getattr(default, key), getattr(explicit, key))


@pytest.mark.parametrize('kwargs', [dict(lod=None), dict(lod=GeometryType.LOD0),
    dict(lod='lod3'), dict(lod2_method='ransac'), dict(lod='lod2', lod2_method='bad'),
    dict(lod='lod2', isolated=True)])
def test_invalid_volume_requests_precede_terrain_work(kwargs):
    with pytest.raises(ValueError):
        builder.build_city_volume_mesh(None, **kwargs)


@pytest.mark.parametrize('kwargs', [dict(lod='lod0'), dict(lod=None), dict(lod2_method='ransac')])
def test_dataset_validates_output_lod(kwargs):
    with pytest.raises(ValueError):
        CityVolumeMeshArgs(bounds=(0, 0, 1, 1), **kwargs)


def test_top_cap_size_seeds_preserve_coarsely_spaced_boundary():
    # Interior target spacing is smaller than the prescribed outer segments.
    vertices = np.array([[0., 0., 0.], [16., 0., 0.], [16., 11., 0.], [0., 11., 0.]])
    loops = tetgen_utils._boundary_loops(vertices, 1e-8)
    top, faces, _ = tetgen_utils._remesh_top_cap_from_outer_boundary(
        vertices, boundary_loops=loops, z_top=30., top_cap_backend='dtcc_mesher',
        top_cap_max_mesh_size=3., top_cap_min_mesh_angle=25., tol=1e-8)
    boundary = tetgen_utils._triangle_mesh_boundary_vertex_indices(faces)
    assert set(map(tuple, top[boundary, :2])) == set(map(tuple, vertices[:, :2]))
    assert len(top) > len(vertices)
    assert triangle_areas(top, faces).sum() == pytest.approx(16.*11.)


def test_volume_preserves_surface_failure_diagnostics(monkeypatch):
    def fail_surface(city, **kwargs):
        kwargs['stage_audit']['roof_catalogue'] = {'block_fallbacks': [{'reason': 'invalid shell'}]}
        raise ValueError('invalid shell')
    monkeypatch.setattr(pipeline, '_require_tetgen', lambda: None)
    monkeypatch.setattr(pipeline, 'get_default_tetgen_switches', lambda: {})
    monkeypatch.setattr(pipeline, 'build_city_surface_mesh', fail_surface)
    audit = {}
    with pytest.raises(ValueError, match='invalid shell'):
        builder.build_city_volume_mesh(None, lod='lod2', stage_audit=audit)
    assert audit['roof_catalogue']['block_fallbacks'][0]['reason'] == 'invalid shell'
    assert audit['attempts'][0]['result']['status'] == 'failed'
