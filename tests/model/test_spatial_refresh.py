"""Explicit local-envelope refresh must not rewrite intrinsic grid domains."""
import numpy as np

from dtcc_core import io
from dtcc_core.model import Bounds, Grid, Mesh, Object, VolumeGrid


def test_mutation_refresh_and_exchange_preserve_domains_and_nonidentity_transforms(tmp_path):
    root, child = Object(id='root'), Object(id='child')
    root.add_child(child)
    mesh = Mesh(vertices=np.array([[0., 0., 1.], [2., 0., 1.], [0., 2., 3.]]),
                faces=np.array([[0, 1, 2]], dtype=np.int64))
    grid_domain = Bounds(-10., -20., -8., -16.)
    volume_domain = Bounds(-5., -5., -1., -1., -5., -1.)
    grid = Grid(width=2, height=4, _bounds=grid_domain.copy())
    volume = VolumeGrid(width=2, height=2, depth=2, _bounds=volume_domain.copy())
    zero_grid_domain = Bounds(-4., -4., -4., -3.)
    zero_volume_domain = Bounds(-3., -4., -3., -3., -4., -3.)
    zero_grid = Grid(width=2, height=2, _bounds=zero_grid_domain.copy())
    zero_volume = VolumeGrid(width=2, height=2, depth=2, _bounds=zero_volume_domain.copy())
    domains = {'grid': grid_domain, 'volume': volume_domain,
               'zero-grid': zero_grid_domain, 'zero-volume': zero_volume_domain}
    for id, geometry in [('mesh', mesh), ('grid', grid), ('volume', volume),
                         ('zero-grid', zero_grid), ('zero-volume', zero_volume)]:
        child.add_geometry(geometry, id=id)
    root.transform.set_translation(1000., 2000., 3000.)
    child.transform.set_translation(100., 200., 300.)
    mesh.transform.set_translation(10., 20., 30.)
    for value in (root, child, mesh):
        value.transform.srs = 'EPSG:3006'

    assert root.bounds == Bounds(-10., -20., 2., 2., -5., 3.)
    for id, domain in domains.items():
        assert child.get_geometry(id=id).bounds == domain
    mesh.vertices[1, 0] = 9.
    assert root.bounds.xmax == 2.  # Public NumPy mutation does not invalidate caches.
    assert root.calculate_bounds() == Bounds(-10., -20., 9., 2., -5., 3.)
    assert mesh.bounds.xmax == 9.
    for id, domain in domains.items():
        assert child.get_geometry(id=id).bounds == domain
    # A single transform is explicit and useful; aggregate bounds did not apply it.
    np.testing.assert_array_equal(mesh.transform(mesh.vertices)[1], [19., 20., 31.])

    path = tmp_path / 'spatial.dtcc'
    io.save_model(root, path)
    restored = io.load_model(path)
    restored_child = restored.get_children(Object)[0]
    restored_mesh = restored_child.get_geometry(id='mesh')
    assert restored.bounds == root.bounds
    for id, domain in domains.items():
        assert restored_child.get_geometry(id=id).bounds == domain
    np.testing.assert_array_equal(restored_mesh.vertices, mesh.vertices)
    for a, b in ((root, restored), (child, restored_child), (mesh, restored_mesh)):
        np.testing.assert_array_equal(a.transform.affine, b.transform.affine)
        assert a.transform.srs == b.transform.srs

    # Default initialization still uses unit cells; later resolution edits leave
    # that initialized domain intact, just as they do an explicitly supplied one.
    default_grid = Grid(width=2, height=3)
    default_volume = VolumeGrid(width=2, height=3, depth=4)
    assert default_grid.bounds == Bounds(0, 0, 2, 3)
    assert default_volume.bounds == Bounds(0, 0, 2, 3, 0, 4)
    default_grid.width = default_volume.width = 8
    assert default_grid.calculate_bounds() == Bounds(0, 0, 2, 3)
    assert default_volume.calculate_bounds() == Bounds(0, 0, 2, 3, 0, 4)
