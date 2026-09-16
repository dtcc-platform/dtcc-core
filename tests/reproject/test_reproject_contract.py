"""Reprojection cannot relabel its source or silently discard admitted state."""

import numpy as np
import pytest

from dtcc_core.model import Field, Object, Point, PointCloud, Surface, Tree, exchange
from dtcc_core.reproject import reproject_pointcloud
from dtcc_core.reproject.reproject import reproject_object


def cloud():
    value = PointCloud(points=np.array([[500000., 6500000., 100.]]),
                       classification=np.array([2], dtype=np.uint8))
    value.transform.srs = 'EPSG:3006'
    return value


def test_public_reprojection_preserves_source_and_native_arrays():
    source = cloud()
    before = exchange.dumps(source)
    _ = source.bounds
    result = reproject_pointcloud(source, None, 'EPSG:4326')
    assert exchange.dumps(source) == before
    assert result.transform is not source.transform and result.transform.srs == 'EPSG:4326'
    assert not np.shares_memory(source.classification, result.classification)
    assert result.points[0, 2] == source.points[0, 2]
    assert result.bounds.xmax < 180 and result.bounds.ymax < 90
    assert exchange.dumps(exchange.loads(exchange.dumps(result))) == exchange.dumps(result)
    source.transform.srs = ''
    assert reproject_pointcloud(source, 'EPSG:3006', 'EPSG:4326').transform.srs == 'EPSG:4326'


def test_reprojection_rejects_state_without_a_transformation_rule():
    source = cloud()
    with pytest.raises(TypeError, match='override_geometry_crs'):
        reproject_pointcloud(source, 'EPSG:32633', 'EPSG:4326', override_geometry_crs='false')
    source.fields = [Field(name='velocity', association='sample', values=np.array([2.]))]
    before = exchange.dumps(source)
    with pytest.raises(NotImplementedError, match='Fields'):
        reproject_pointcloud(source, 'EPSG:3006', 'EPSG:4326')
    assert exchange.dumps(source) == before
    source.fields.clear()
    source.transform.set_translation(1, 0, 0)
    with pytest.raises(NotImplementedError, match='affine'):
        reproject_pointcloud(source, 'EPSG:3006', 'EPSG:4326')
    source.transform.set_translation(0, 0, 0)
    with pytest.raises(NotImplementedError, match='two-axis'):
        reproject_pointcloud(source, 'EPSG:3006', 'EPSG:4978')
    integer_source = PointCloud(points=np.array([[500000, 6500000, 2**53 + 1]], dtype=np.uint64))
    with pytest.raises(ValueError, match='integer Z'):
        reproject_pointcloud(integer_source, 'EPSG:3006', 'EPSG:4326')
    surface = Surface(vertices=np.vstack((source.points, source.points + [1, 0, 0],
                                         source.points + [0, 1, 0])),
                      fields=[Field(name='temperature', association='face', values=np.array([2.]))])
    obj = Object(id='o')
    obj.add_geometry(surface, id='survey')
    before = exchange.dumps(obj)
    with pytest.raises(NotImplementedError, match='Fields'):
        reproject_object(obj, 'EPSG:3006', 'EPSG:4326')
    assert exchange.dumps(obj) == before
    obj.geometry.clear()
    obj.add_child(Object(id='child'))
    with pytest.raises(NotImplementedError, match='Nested object'):
        reproject_object(obj, 'EPSG:3006', 'EPSG:4326')
    obj.children.clear()
    obj.add_geometry(Point(), id='unsupported')
    with pytest.raises(NotImplementedError, match='Point reprojection'):
        reproject_object(obj, 'EPSG:3006', 'EPSG:4326')
    with pytest.raises(NotImplementedError, match='Tree object reprojection'):
        reproject_object(Tree(position=np.array([500000., 6500000., 100.])),
                         'EPSG:3006', 'EPSG:4326')
