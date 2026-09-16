import pytest
import numpy as np
from affine import Affine
from pathlib import Path
from dtcc_core import io
from dtcc_core.model import Raster
from dtcc_core.builder.model_conversion import raster_to_builder_gridfield
from dtcc_core.builder import _dtcc_builder


@pytest.fixture
def data_dir():
    return (Path(__file__).parent / ".." / "data" / "rasters").resolve()


@pytest.fixture
def dem_raster(data_dir):
    return io.load_raster(data_dir / "test_dem.tif")


@pytest.fixture
def dem_gridfield(dem_raster):
    return raster_to_builder_gridfield(dem_raster)


def test_gridfield_dimensions(dem_gridfield):
    assert dem_gridfield.grid.xsize == 20
    assert dem_gridfield.grid.ysize == 40
    assert dem_gridfield.grid.xstep == 2.0
    assert dem_gridfield.grid.ystep == 2.0


def test_gridfield_values(dem_gridfield):
    values = dem_gridfield.values
    assert len(values) == 20 * 40

    # Check specific points
    assert pytest.approx(values[0], rel=1e-5) == 76  # bottom left corner
    assert pytest.approx(values[1], rel=1e-5) == 76.2  # next to bottom left
    assert pytest.approx(values[(20 * 40) - 20], rel=1e-5) == 0  # top left corner
    assert pytest.approx(values[-1], rel=1e-5) == 3.8  # top right corner


def _sample_field(field, xs, ys):
    xx, yy = np.meshgrid(xs, ys)
    vertices = np.column_stack((xx.ravel(), yy.ravel(), np.zeros(xx.size)))
    faces = []
    for row in range(len(ys) - 1):
        for col in range(len(xs) - 1):
            i = row * len(xs) + col
            faces.extend([[i, i + 1, i + len(xs)],
                          [i + 1, i + len(xs) + 1, i + len(xs)]])
    ground = _dtcc_builder.create_mesh(vertices, np.array(faces), np.full(len(faces), -2))
    mesh = _dtcc_builder.build_terrain_surface_mesh_from_ground_mesh(ground, field, 0)
    return _dtcc_builder.mesh_as_arrays(mesh)[0].reshape((-1, 3))[:, 2].reshape(xx.shape)


@pytest.mark.parametrize("georef", [
    Affine(2, 0, 10, 0, -3, 29),  # north-up
    Affine(2, 0, 10, 0, 3, 20),   # south-up
    Affine(-2, 0, 16, 0, -3, 29), # reversed columns
])
def test_gridfield_samples_pixel_centers_and_clamps_outer_margins(georef):
    cols, rows = np.meshgrid(np.arange(3) + 0.5, np.arange(3) + 0.5)
    x, y = georef * (cols, rows)
    raster = Raster(data=x + 2 * y, georef=georef, nodata=-9999)
    before = raster.data.copy()
    field = raster_to_builder_gridfield(raster)
    xs, ys = [10, 11, 12, 15, 16], [20, 21.5, 23, 27.5, 29]
    expected = np.clip(xs, 11, 15)[None, :] + 2 * np.clip(ys, 21.5, 27.5)[:, None]
    np.testing.assert_allclose(_sample_field(field, xs, ys), expected)
    assert field.grid.xstep == 2
    assert field.grid.ystep == 3
    np.testing.assert_array_equal(raster.data, before)


def test_single_pixel_gridfield_is_constant_over_its_footprint():
    field = raster_to_builder_gridfield(Raster(data=np.array([[7.0]])))
    np.testing.assert_array_equal(_sample_field(field, [0, 0.5, 1], [0, 0.5, 1]), 7)


@pytest.mark.parametrize("raster, error, message", [
    (Raster(data=np.empty((0, 2))), ValueError, "nonempty 2D"),
    (Raster(data=np.zeros((2, 2, 3))), ValueError, "single-band"),
    (Raster(data=np.ones((2, 2)), georef=Affine(1, 0.1, 0, 0, -1, 2)),
     NotImplementedError, "rotated or skewed"),
    (Raster(data=np.ones((2, 2)), georef=Affine(0, 0, 0, 0, -1, 2)),
     ValueError, "spacing"),
    (Raster(data=np.array([[1, -9999]]), nodata=-9999), ValueError, "nodata"),
    (Raster(data=np.ma.array([[1, 2]], mask=[[False, True]])), ValueError, "masked"),
    (Raster(data=np.array([[1, np.nan]])), ValueError, "finite real"),
])
def test_raster_conversion_rejects_unsupported_input(raster, error, message):
    with pytest.raises(error, match=message):
        raster_to_builder_gridfield(raster)


@pytest.mark.parametrize("data, bounds, width, height, message", [
    (np.ones(2), (0, 0, 2, 2), 2, 2, "matching"),
    (np.ones(2), (0, 0, 2, 2), 0, 2, "positive"),
    (np.ones(2), (0, 0, 0, 2), 1, 2, "bounds"),
    (np.array([1, np.inf]), (0, 0, 1, 2), 1, 2, "finite real"),
    (np.array([1 + 1j, 2]), (0, 0, 1, 2), 1, 2, "finite real"),
])
def test_native_gridfield_validates_raw_inputs(data, bounds, width, height, message):
    with pytest.raises(ValueError, match=message):
        _dtcc_builder.create_gridfield(data, bounds, width, height)


def test_terrain_workflows_use_raster_pixel_locations():
    from dtcc_core.builder.geometry_builders.terrain import (
        adaptive_terrain_mesh, build_terrain_surface_mesh,
    )
    raster = Raster(data=np.array([[5., 6., 7.], [3., 4., 5.], [1., 2., 3.]]),
                    georef=Affine(1, 0, 0, 0, -1, 3))
    mesh = build_terrain_surface_mesh(raster, max_mesh_size=1, smoothing=0,
                                      report_mesh_quality=False)
    x, y, z = mesh.vertices.T
    np.testing.assert_allclose(z, np.clip(x, 0.5, 2.5) + 2 * np.clip(y, 0.5, 2.5) - 0.5)
    adaptive = adaptive_terrain_mesh(raster, max_error=0.1, smoothing=0)
    np.testing.assert_allclose(adaptive.vertices[:, :2].min(axis=0), [0.5, 0.5])
    np.testing.assert_allclose(adaptive.vertices[:, :2].max(axis=0), [2.5, 2.5])
    x, y, z = adaptive.vertices.T
    np.testing.assert_allclose(z, x + 2 * y - 0.5)
    raster.data[0, 0] = np.nan
    with pytest.raises(ValueError, match="fill missing data"):
        build_terrain_surface_mesh(raster, report_mesh_quality=False)


if __name__ == "__main__":
    pytest.main()
