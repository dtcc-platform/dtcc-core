import numpy as np
import pytest

import matplotlib
matplotlib.use('Agg', force=True)
import matplotlib.pyplot as plt
from affine import Affine

from dtcc_core.model import (
    Bounds, City, Field, Grid, Mesh, Object, Point, PointCloud, Raster, Solid,
    Surface, VolumeGrid, exchange,
)


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close('all')


def test_object_preview_selects_detail_and_preserves_local_transform_and_holes():
    city = City(id='city')
    obj = Object(id='bench')
    obj.add_geometry(Point(x=1000, y=1000), id='location', lod='0')
    surface = Surface(vertices=np.array([[0.,0.,0.],[2.,0.,0.],[2.,2.,0.],[0.,2.,0.]]),
                      holes=[np.array([[.5,.5,0.],[.5,1.,0.],[1.,1.,0.],[1.,.5,0.]])])
    surface.transform.set_translation(10, 20, 3)
    obj.add_geometry(surface, id='detail', lod='2.2')
    city.add_child(obj)
    before = exchange.dumps(city)
    ax = city.plot(show=False)
    assert ax.name == '3d'
    assert 9 < ax.get_xlim()[0] < 10 and 12 < ax.get_xlim()[1] < 13
    assert any('ring outlines' in text.get_text() for text in ax.texts)
    assert exchange.dumps(city) == before
    coarse = city.plot(lod='0', show=False)
    assert coarse.get_xlim()[0] > 990
    with pytest.raises(ValueError, match='No geometry matches'):
        city.plot(representation='missing', show=False)


def test_vector_field_samples_use_association_locations_and_magnitude():
    city = City(id='city')
    obj = Object(id='flow')
    obj.add_geometry(Solid(surfaces=[Surface(vertices=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]]))],
                           shells=[np.array([0])]), id='solid', lod='3')
    cloud = PointCloud(points=np.array([[10.,20.,30.],[40.,50.,60.]]))
    cloud.fields = [Field(name='velocity', unit='m/s', association='vertex', dim=3,
                          values=np.array([[3.,4.,0.],[5.,12.,0.]]))]
    obj.add_geometry(cloud, id='samples')
    city.add_child(obj)
    ax = city.plot(field='velocity', show=False)
    np.testing.assert_allclose(ax.collections[0].get_array(), [5., 13.])
    assert ax.get_zlim()[0] < 30 and ax.get_zlim()[1] > 60
    assert 'magnitude' in ax.figure.axes[1].get_ylabel()


def test_preview_bounds_work_and_rejects_ambiguous_frames():
    cloud = PointCloud(points=np.column_stack([np.arange(1000), np.zeros(1000), np.ones(1000)]))
    ax = cloud.plot(max_elements=17, show=False)
    assert len(ax.collections[0].get_offsets()) == 17
    assert any('sampled' in text.get_text() for text in ax.texts)
    grid = VolumeGrid(width=2, height=3, depth=4)
    assert grid._bounds is None
    grid.plot(show=False)
    assert grid._bounds is None
    with pytest.raises(ValueError, match='positive integer'):
        grid.plot(max_elements=0, show=False)
    city = City(id='city')
    city.transform.set_translation(1, 2, 3)
    with pytest.raises(ValueError, match='Object transforms'):
        city.plot(show=False)
    city.transform.affine = np.eye(4)
    a, b = Object(id='a'), Object(id='b')
    a.transform.srs, b.transform.srs = 'EPSG:3006', 'EPSG:7415'
    a.add_geometry(Point(x=1), id='location')
    b.add_geometry(Point(x=2), id='location')
    city.add_children([a,b])
    with pytest.raises(ValueError, match='common CRS'):
        city.plot(show=False)


def test_raster_nodata_and_explicit_elevation_meaning():
    from affine import Affine
    from dtcc_core.model import Terrain
    terrain = Terrain(id='dem')
    raster = Raster(data=np.array([[2., np.nan], [-1., 4.]]), georef=Affine(10,0,100,0,-10,200))
    terrain.add_geometry(raster, id='dem')
    terrain.attributes['elevation_rasters'] = [{'geometry_id': 'dem', 'unit': 'm'}]
    ax = terrain.plot(representation='dem', show=False)
    np.testing.assert_array_equal(ax.collections[0].get_array(), [2., -1., 4.])
    assert ax.get_zlim()[0] < -1 and ax.get_zlim()[1] > 4
    flat = raster.plot(show=False)
    assert flat.get_zlim() == pytest.approx((-.5, .5))
    raster.data[:] = np.nan
    with pytest.raises(ValueError, match='No finite'):
        raster.plot(show=False)


@pytest.mark.parametrize('channels', [3, 4])
@pytest.mark.parametrize('floating', [False, True])
def test_image_raster_preserves_colours_alpha_and_model(channels, floating):
    data = np.array([[[255, 0, 0, 255], [0, 0, 0, 255]],
                     [[0, 255, 0, 128], [0, 0, 255, 0]]], dtype=np.uint8)[..., :channels]
    if floating:
        data = data.astype(float) / 255
    raster = Raster(data=data, georef=Affine(10, 0, 100, 0, -10, 200), crs='EPSG:3006')
    before = exchange.dumps(raster)
    ax = raster.plot(show=False)
    assert ax.name == 'rectilinear'
    assert len(ax.images) == 1 and len(ax.figure.axes) == 1
    np.testing.assert_array_equal(ax.images[0].get_array(), data)
    assert ax.get_xlim() == pytest.approx((100, 120))
    assert ax.get_ylim() == pytest.approx((180, 200))
    assert ax.get_aspect() == 1
    ax.figure.canvas.draw()
    assert exchange.dumps(raster) == before


@pytest.mark.parametrize('georef', [Affine(10, 0, 100, 0, -20, 200),
                                   Affine(10, 0, 100, 0, 20, 200),
                                   Affine(10, 3, 100, 2, -20, 200)])
def test_image_raster_places_pixel_centres_using_the_full_affine(georef):
    raster = Raster(data=np.full((2, 3, 3), 255, dtype=np.uint8), georef=georef)
    ax = raster.plot(show=False)
    artist = ax.images[0]
    assert artist.origin == 'upper'
    assert tuple(artist.get_extent()) == (0, 3, 2, 0)
    centres = np.array([[.5, .5], [2.5, 1.5]])
    plotted = ax.transData.inverted().transform(artist.get_transform().transform(centres))
    np.testing.assert_allclose(plotted, [georef * tuple(p) for p in centres])
    b = raster.bounds
    assert ax.get_xlim() == pytest.approx((b.xmin, b.xmax))
    assert ax.get_ylim() == pytest.approx((b.ymin, b.ymax))


@pytest.mark.parametrize('shape', [(21, 33), (1, 1000), (1000, 1)])
def test_large_image_preview_is_bounded_without_cropping(shape):
    raster = Raster(data=np.full((*shape, 4), 255, dtype=np.uint8),
                    georef=Affine(2, 0, 100, 0, -2, 200))
    ax = raster.plot(max_elements=17, show=False)
    sampled = ax.images[0].get_array()
    assert sampled.shape[0] * sampled.shape[1] <= 17
    assert sampled.shape[-1] == 4
    assert tuple(ax.images[0].get_extent()) == (0, shape[1], shape[0], 0)
    b = raster.bounds
    assert ax.get_xlim() == pytest.approx((b.xmin, b.xmax))
    assert ax.get_ylim() == pytest.approx((b.ymin, b.ymax))
    assert any('sampled' in text.get_text() for text in ax.texts)


def test_image_preview_uses_supplied_axes_and_show_option(monkeypatch):
    raster = Raster(data=np.full((2, 2, 3), 255, dtype=np.uint8))
    _, ax = plt.subplots()
    shown = []
    monkeypatch.setattr(plt, 'show', lambda: shown.append(True))
    assert raster.plot(ax=ax, show=False, theme='light') is ax
    assert shown == []
    raster.plot()
    assert shown == [True]
    ax3d = plt.figure().add_subplot(projection='3d')
    with pytest.raises(ValueError, match='2D axes'):
        raster.plot(ax=ax3d, show=False)
    with pytest.raises(ValueError, match='positive integer'):
        raster.plot(max_elements=0, show=False)
    with pytest.raises(ValueError, match='selectors require an Object'):
        raster.plot(lod='0', show=False)


def test_empty_image_and_unsupported_band_layouts_fail_clearly():
    with pytest.raises(ValueError, match='No image pixels'):
        Raster(data=np.empty((0, 2, 4), dtype=np.uint8)).plot(show=False)
    with pytest.raises(ValueError, match='scalar 2D raster'):
        Raster(data=np.zeros((2, 2, 2), dtype=np.uint8)).plot(show=False)
