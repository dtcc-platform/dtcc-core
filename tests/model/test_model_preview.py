import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt

from dtcc_core.model import (
    Bounds,
    City,
    Field,
    Grid,
    Mesh,
    Object,
    Point,
    PointCloud,
    Raster,
    Solid,
    Surface,
    VolumeGrid,
    exchange,
)


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close("all")


def test_object_preview_selects_detail_and_preserves_local_transform_and_holes():
    city = City(id="city")
    obj = Object(id="bench")
    obj.add_geometry(Point(x=1000, y=1000), id="location", lod="0")
    surface = Surface(
        vertices=np.array(
            [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 2.0, 0.0], [0.0, 2.0, 0.0]]
        ),
        holes=[
            np.array(
                [[0.5, 0.5, 0.0], [0.5, 1.0, 0.0], [1.0, 1.0, 0.0], [1.0, 0.5, 0.0]]
            )
        ],
    )
    surface.transform.set_translation(10, 20, 3)
    obj.add_geometry(surface, id="detail", lod="2.2")
    city.add_child(obj)
    before = exchange.dumps(city)
    ax = city.plot(show=False)
    assert ax.name == "3d"
    assert 9 < ax.get_xlim()[0] < 10 and 12 < ax.get_xlim()[1] < 13
    assert any("ring outlines" in text.get_text() for text in ax.texts)
    assert exchange.dumps(city) == before
    coarse = city.plot(lod="0", show=False)
    assert coarse.get_xlim()[0] > 990
    with pytest.raises(ValueError, match="No geometry matches"):
        city.plot(representation="missing", show=False)


def test_vector_field_samples_use_association_locations_and_magnitude():
    city = City(id="city")
    obj = Object(id="flow")
    obj.add_geometry(
        Solid(
            surfaces=[
                Surface(
                    vertices=np.array(
                        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
                    )
                )
            ],
            shells=[np.array([0])],
        ),
        id="solid",
        lod="3",
    )
    cloud = PointCloud(points=np.array([[10.0, 20.0, 30.0], [40.0, 50.0, 60.0]]))
    cloud.fields = [
        Field(
            name="velocity",
            unit="m/s",
            association="vertex",
            dim=3,
            values=np.array([[3.0, 4.0, 0.0], [5.0, 12.0, 0.0]]),
        )
    ]
    obj.add_geometry(cloud, id="samples")
    city.add_child(obj)
    ax = city.plot(field="velocity", show=False)
    np.testing.assert_allclose(ax.collections[0].get_array(), [5.0, 13.0])
    assert ax.get_zlim()[0] < 30 and ax.get_zlim()[1] > 60
    assert "magnitude" in ax.figure.axes[1].get_ylabel()


def test_preview_bounds_work_and_rejects_ambiguous_frames():
    cloud = PointCloud(
        points=np.column_stack([np.arange(1000), np.zeros(1000), np.ones(1000)])
    )
    ax = cloud.plot(max_elements=17, show=False)
    assert len(ax.collections[0].get_offsets()) == 17
    assert any("sampled" in text.get_text() for text in ax.texts)
    grid = VolumeGrid(width=2, height=3, depth=4)
    assert grid._bounds is None
    grid.plot(show=False)
    assert grid._bounds is None
    with pytest.raises(ValueError, match="positive integer"):
        grid.plot(max_elements=0, show=False)
    city = City(id="city")
    city.transform.set_translation(1, 2, 3)
    with pytest.raises(ValueError, match="Object transforms"):
        city.plot(show=False)
    city.transform.affine = np.eye(4)
    a, b = Object(id="a"), Object(id="b")
    a.transform.srs, b.transform.srs = "EPSG:3006", "EPSG:7415"
    a.add_geometry(Point(x=1), id="location")
    b.add_geometry(Point(x=2), id="location")
    city.add_children([a, b])
    with pytest.raises(ValueError, match="common CRS"):
        city.plot(show=False)


def test_raster_nodata_and_explicit_elevation_meaning():
    from affine import Affine

    from dtcc_core.model import Terrain

    terrain = Terrain(id="dem")
    raster = Raster(
        data=np.array([[2.0, np.nan], [-1.0, 4.0]]),
        georef=Affine(10, 0, 100, 0, -10, 200),
    )
    terrain.add_geometry(raster, id="dem")
    terrain.attributes["elevation_rasters"] = [{"geometry_id": "dem", "unit": "m"}]
    ax = terrain.plot(representation="dem", show=False)
    np.testing.assert_array_equal(ax.collections[0].get_array(), [2.0, -1.0, 4.0])
    assert ax.get_zlim()[0] < -1 and ax.get_zlim()[1] > 4
    flat = raster.plot(show=False)
    assert flat.get_zlim() == pytest.approx((-0.5, 0.5))
    raster.data[:] = np.nan
    with pytest.raises(ValueError, match="No finite"):
        raster.plot(show=False)
