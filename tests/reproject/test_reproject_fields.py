"""Coordinate reprojection preserves attached samples and their component convention."""

import warnings

import numpy as np
import pytest

from dtcc_core.model import Field, Mesh, MultiSurface, Object, Surface, VolumeMesh, exchange
from dtcc_core.reproject.reproject import reproject_mesh, reproject_object


@pytest.fixture
def vertices():
    return np.array([
        [500000., 6500000., 100.],
        [500010., 6500000., 100.],
        [500000., 6500010., 100.],
        [500000., 6500000., 110.],
    ])


@pytest.mark.parametrize('mesh_type', [Mesh, VolumeMesh])
def test_reproject_mesh_preserves_scalar_fields_and_topology(vertices, mesh_type):
    if mesh_type is Mesh:
        topology = 'faces'
        association = 'face'
        elements = np.array([[0, 1, 2]])
    else:
        topology = 'cells'
        association = 'cell'
        elements = np.array([[0, 1, 2, 3]])
    source = mesh_type(vertices=vertices, markers=np.array([7]), **{topology: elements})
    source.transform.srs = 'EPSG:3006'
    source.fields = [
        Field(name='temperature', unit='degC', description='Simulation result',
              association='vertex', values=np.array([[20.], [21.], [22.], [23.]])),
        Field(name='mass', unit='kg', association=association, values=np.array([5.])),
    ]
    before = exchange.dumps(source)
    _ = source.bounds

    with warnings.catch_warnings():
        warnings.simplefilter('error')
        result = reproject_mesh(source, 'EPSG:3006', 'EPSG:4326')

    assert type(result) is mesh_type
    assert result.transform.srs == 'EPSG:4326'
    assert not np.allclose(result.vertices[:, :2], vertices[:, :2])
    np.testing.assert_array_equal(result.vertices[:, 2], vertices[:, 2])
    np.testing.assert_array_equal(getattr(result, topology), elements)
    np.testing.assert_array_equal(result.markers, source.markers)
    assert result.bounds.xmax < 180 and result.bounds.ymax < 90
    for original, copied in zip(source.fields, result.fields, strict=True):
        assert exchange.dumps(copied) == exchange.dumps(original)
        assert not np.shares_memory(copied.values, original.values)
    assert exchange.dumps(source) == before
    assert exchange.dumps(exchange.loads(exchange.dumps(result))) == exchange.dumps(result)


def test_reproject_object_preserves_multisurface_and_surface_fields(vertices):
    surface = Surface(vertices=vertices[:3], fields=[
        Field(name='temperature', unit='degC', association='face', values=np.array([20.])),
    ])
    geometry = MultiSurface(surfaces=[surface], fields=[
        Field(name='label', association='face', values=np.array([7], dtype=np.uint8)),
    ])
    source = Object()
    source.add_geometry(geometry, id='survey', lod='2.2')
    before = exchange.dumps(source)

    result = reproject_object(source, 'EPSG:3006', 'EPSG:4326')

    copied = result.get_geometry(id='survey')
    assert not np.allclose(copied.surfaces[0].vertices[:, :2], vertices[:3, :2])
    for original, transformed in ((geometry, copied), (surface, copied.surfaces[0])):
        assert transformed.transform.srs == 'EPSG:4326'
        assert exchange.dumps(transformed.fields[0]) == exchange.dumps(original.fields[0])
        assert not np.shares_memory(transformed.fields[0].values, original.fields[0].values)
    assert exchange.dumps(source) == before


def test_reproject_volume_result_warns_and_preserves_vector_components(vertices):
    mesh = VolumeMesh(vertices=vertices, cells=np.array([[0, 1, 2, 3]]), fields=[
        Field(name='velocity', unit='m/s', dim=3, association='vertex',
              values=np.tile([1., 2., 3.], (4, 1))),
    ])
    source = Object()
    source.add_geometry(mesh, id='simulation')
    before = exchange.dumps(source)

    with pytest.warns(UserWarning, match='Vector field components are preserved unchanged') as caught:
        result = reproject_object(source, 'EPSG:3006', 'EPSG:3007')

    assert len(caught) == 1
    copied = result.get_geometry(id='simulation')
    assert isinstance(copied, VolumeMesh)
    assert copied.transform.srs == 'EPSG:3007'
    assert not np.allclose(copied.vertices[:, :2], vertices[:, :2])
    np.testing.assert_array_equal(copied.cells, mesh.cells)
    assert exchange.dumps(copied.fields[0]) == exchange.dumps(mesh.fields[0])
    copied.fields[0].values[0, 0] = 99.
    assert exchange.dumps(source) == before

    with warnings.catch_warnings():
        warnings.simplefilter('error')
        assert reproject_mesh(mesh, 'EPSG:3006', 'EPSG:3006') is mesh
