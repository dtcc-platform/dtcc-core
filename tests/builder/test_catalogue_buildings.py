"""Catalogue city-model shells and CityJSON export preserve reconstructed roofs."""
import numpy as np
from collections import Counter
import pytest
from shapely import affinity
from shapely.geometry import box

from dtcc_core import builder, io
from dtcc_core.model import GeometryType
from dtcc_core.builder.geometry_builders.lod2 import is_watertight
from .test_roof_catalogue import SHAPES, building, model, city_of


@pytest.mark.parametrize('shape,parameters', SHAPES)
@pytest.mark.parametrize('offset', [(0., 0.), (320000., 6400000.)])
def test_catalogue_city_shell_uses_shared_fit_and_closed_geometry(shape, parameters, offset):
    polygon = affinity.translate(box(-12, -6, 12, 6), *offset)
    truth = model(shape, parameters, polygon)
    source = building(polygon, truth)
    footprint = source.lod0.vertices.copy()
    builder.build_lod2_buildings([source])
    assert source.attributes['roof_reconstruction']['shape'] == shape
    assert source.attributes['roof_reconstruction']['outcome'] == 'accepted_fit'
    assert is_watertight(source.lod2)
    directed = Counter()
    for surface in source.lod2.surfaces:
        vertices = np.round(surface.vertices, 6)
        directed.update((tuple(a), tuple(b)) for a, b in zip(vertices, np.roll(vertices, -1, axis=0)))
    assert all(count == 1 and directed[b, a] == 1 for (a, b), count in directed.items())
    np.testing.assert_array_equal(source.lod0.vertices, footprint)
    roofs = [s for s in source.lod2.surfaces if s.calculate_normal()[2] > .1]
    assert roofs
    for roof in roofs:
        np.testing.assert_allclose(roof.vertices[:, 2], truth.predict(roof.vertices[:, :2]), atol=.1)


def test_catalogue_city_export_round_trip_and_sparse_fallback(tmp_path):
    polygon = box(-12, -6, 12, 6)
    source = building(polygon, model('gable', [0., 8., .55], polygon))
    sparse_polygon = box(20, -5, 30, 5)
    sparse = building(sparse_polygon, model('flat', [8.], sparse_polygon, key='sparse'), observations=False)
    city = city_of([source, sparse])
    city.build_lod2_buildings(calculate_heights=False)
    assert source.lod2 is not None and sparse.lod1 is not None
    assert sparse.lod2 is None
    assert sparse.attributes['roof_reconstruction']['outcome'] == 'lod1_fallback'
    path = tmp_path/'city.city.json'
    city.save(path)
    loaded = io.load_city(path)
    restored = next(b for b in loaded.buildings if b.id == source.id)
    shell = restored.flatten_geometry(GeometryType.LOD2)
    assert shell is not None and is_watertight(shell)
    original_xyz = np.vstack([s.vertices for s in source.lod2.surfaces])
    loaded_xyz = np.vstack([s.vertices for s in shell.surfaces])
    np.testing.assert_allclose(loaded_xyz, original_xyz, atol=.002)


def test_switching_to_ransac_does_not_keep_catalogue_diagnostics():
    polygon = box(-12, -6, 12, 6)
    source = building(polygon, model('gable', [0., 8., .55], polygon))
    builder.build_lod2_buildings([source])
    assert source.attributes['roof_reconstruction']['method'] == 'catalogue'
    builder.build_lod2_buildings([source], lod2_method='ransac')
    assert source.lod2 is not None
    assert 'roof_reconstruction' not in source.attributes
