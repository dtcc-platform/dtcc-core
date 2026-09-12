"""One public Protobuf contract, including the previously legacy-only state."""

import numpy as np
import pytest
from affine import Affine

from dtcc_core import io, model
from dtcc_core.model import exchange


def test_numerical_and_object_state_share_one_message(tmp_path):
    city = model.City(id='city')
    terrain = model.Terrain(id='terrain')
    raster = model.Raster(data=np.array([[[2**64 - 1], [2**53 + 1]]], dtype=np.uint64),
                         georef=Affine(2, .25, 100, .5, -2, 200), crs='EPSG:3006')
    terrain.add_geometry(raster, id='elevation')
    grid = model.VolumeGrid(width=2, height=3, depth=4)
    grid.bounds = model.Bounds(xmin=10, ymin=20, zmin=30, xmax=14, ymax=26, zmax=38)
    grid.fields = [model.Field(name='velocity', association='cell', dim=3,
                              values=np.arange(72, dtype=np.float64).reshape(24, 3))]
    terrain.add_geometry(grid, id='simulation')
    city.add_child(terrain)
    roads = model.RoadNetwork(id='roads', vertices=np.array([[1., 2.], [3., 4.]]),
                              edges=np.array([[0, 1]], dtype=np.uint32), length=np.array([2.75]))
    city.add_child(roads)
    city.add_child(model.SensorCollection(id='sensors'))
    path = tmp_path / 'city.dtcc'
    city.save(path)
    # A generic generated Protobuf class decodes the whole file directly.
    message = model.proto.ModelFile.FromString(path.read_bytes())
    assert message.DESCRIPTOR.file.name == 'dtcc.proto'
    assert message.DESCRIPTOR.full_name == 'DTCC.ModelFile'
    restored = model.City()
    restored.from_proto(message)
    result = restored.get_children(model.Terrain)[0]
    result_raster = result.get_geometry(id='elevation')
    np.testing.assert_array_equal(result_raster.data, raster.data)
    assert result_raster.data.dtype == raster.data.dtype
    assert result_raster.georef == raster.georef and result_raster.crs == raster.crs
    result_grid = result.get_geometry(id='simulation')
    assert result_grid.bounds == grid.bounds
    np.testing.assert_array_equal(result_grid.fields[0].values, grid.fields[0].values)
    result_roads = restored.get_children(model.RoadNetwork)[0]
    np.testing.assert_array_equal(result_roads.vertices, roads.vertices)
    np.testing.assert_array_equal(result_roads.edges, roads.edges)
    np.testing.assert_array_equal(result_roads.length, roads.length)
    assert type(restored.get_children(model.SensorCollection)[0]) is model.SensorCollection


def test_new_file_wrappers_and_failed_read_preserve_receiver(tmp_path):
    cloud = model.PointCloud(points=np.array([[1.25, 2.5, 3.75]]), classification=np.array([2], dtype=np.uint8))
    mesh = model.VolumeMesh(vertices=np.eye(4, 3), cells=np.array([[0, 1, 2, 3]]),
                            markers=np.array([7], dtype=np.int32))
    raster = model.Raster(data=np.array([[1.00000000000001, np.nan]], dtype=np.float64))
    for value, reader in [(cloud, io.load_pointcloud), (mesh, io.load_volume_mesh), (raster, io.load_raster)]:
        path = tmp_path / f'{type(value).__name__}.dtcc'
        value.save(path)
        restored = reader(path, validate_schema=False)
        assert type(restored) is type(value)
        assert exchange.dumps(restored) == exchange.dumps(value)
    original = mesh.to_proto().SerializeToString()
    malformed = mesh.to_proto()
    malformed.geometry.volume_mesh.cells.data = b'bad'
    with pytest.raises(ValueError, match='data length'):
        mesh.from_proto(malformed, validate_schema=False)
    assert mesh.to_proto().SerializeToString() == original
    with pytest.raises(ValueError, match='Expected VolumeMesh'):
        mesh.from_proto(cloud.to_proto())
    assert mesh.to_proto().SerializeToString() == original
    mesh.cells[0, 0] = 99
    with pytest.raises(ValueError, match='out-of-range'):
        mesh.save(tmp_path / 'VolumeMesh.dtcc', validate_schema=False)
    assert (tmp_path / 'VolumeMesh.dtcc').read_bytes() == original


def test_retired_wire_versions_and_unknown_fields_fail_without_fallback():
    source = model.Point(x=1).to_proto()
    for version in range(1, 6):
        source.version = version
        with pytest.raises(ValueError, match='Unsupported model format/version'):
            model.Point().from_proto(source, validate_schema=False)
    source.version = exchange.VERSION
    with pytest.raises(NotImplementedError, match='unsupported fields'):
        model.Point().from_proto(source.SerializeToString() + b'\xa0\x06\x01')
    # An old root Object message cannot be mistaken for ModelFile.
    with pytest.raises(ValueError, match='Unsupported model format/version'):
        model.Object().from_proto(b'\x0a\x03old')


def test_both_protobuf_entry_points_enforce_the_semantic_boundary():
    value = model.Building(id='building', attributes={'measured_height': -1})
    with pytest.raises(ValueError, match='measured_height'):
        value.to_proto()
    message = value.to_proto(validate_schema=False)
    receiver = model.Building(id='unchanged')
    with pytest.raises(ValueError, match='measured_height'):
        receiver.from_proto(message)
    assert receiver.id == 'unchanged'
    receiver.from_proto(message, validate_schema=False)
    assert receiver.attributes['measured_height'] == -1


def test_malformed_value_metadata_fails_before_replacing_native_state():
    field = model.Field(name='original', association='sample', values=np.array([2.]))
    message = field.to_proto()
    message.field.dim = 2  # The array has one scalar component, not two.
    with pytest.raises(ValueError, match='Field.values'):
        field.from_proto(message, validate_schema=False)
    assert field.name == 'original' and field.dim == 1

    raster = model.Raster(data=np.array([[2.]]), crs='EPSG:3006')
    message = raster.to_proto()
    del message.raster.georef[-1]
    with pytest.raises(ValueError, match='six coefficients'):
        raster.from_proto(message, validate_schema=False)
    assert raster.crs == 'EPSG:3006'
    message = raster.to_proto()
    message.raster.data.dtype = 'object'
    with pytest.raises(ValueError, match='dtype'):
        raster.from_proto(message, validate_schema=False)
    assert raster.data.tolist() == [[2.]]
