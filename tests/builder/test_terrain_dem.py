"""Point-source truth, missing support, sampling and native DEM persistence."""

import numpy as np
import pytest

from dtcc_core import builder, io
from dtcc_core.datasets import load_model_package
from dtcc_core.model import Bounds, Raster, exchange
from dtcc_core.builder.raster.interpolation import fill_holes
from sandbox.model_profiles.dem_example import make_dem, source


def test_ground_raster_retains_zero_negative_and_crs_with_explicit_gaps():
    pc = source()
    raster = builder.build_terrain_raster(
        pc,
        1.0,
        Bounds(0, 0, 3.2, 3.2),
        window_size=0,
        radius=0.1,
        hole_fill="none",
    )
    assert raster.crs == pc.transform.srs and np.isnan(raster.nodata)
    assert raster.shape == (4, 4)
    # Cell centers remain anchored to xmin/ymin when dimensions are rounded up.
    assert raster.georef * (0.5, 3.5) == (0.5, 0.5)
    np.testing.assert_array_equal(raster.data[-1, :2], [0.0, -2.0])
    assert np.isnan(raster.data).sum() == 14  # Water sample did not become terrain.
    filled = fill_holes(raster)
    assert np.isfinite(filled.data).all() and set(np.unique(filled.data)) == {0.0, -2.0}
    assert np.isnan(raster.data).sum() == 14  # Source was not modified.
    all_points = builder.build_terrain_raster(
        pc,
        1.0,
        Bounds(0, 0, 3.2, 3.2),
        window_size=0,
        radius=0.1,
        hole_fill="none",
        ground_only=False,
    )
    assert 99.0 in all_points.data


def test_no_ground_or_support_fails_instead_of_substituting_points():
    pc = source()
    pc.classification[:] = 9
    with pytest.raises(ValueError, match="class 2"):
        builder.build_terrain_raster(pc, 1.0)
    pc.classification = np.array([], dtype=np.uint8)
    with pytest.raises(ValueError, match="classification"):
        builder.build_terrain_raster(pc, 1.0)
    with pytest.raises(ValueError, match="support"):
        builder.build_terrain_raster(
            source(), 1.0, Bounds(10, 10, 12, 12), window_size=0
        )
    with pytest.raises(ValueError, match="valid values"):
        fill_holes(Raster(data=np.full((2, 2), np.nan)))
    with pytest.raises(ValueError, match="cell-index limit"):
        builder.build_terrain_raster(source(), 1e-100, Bounds(0, 0, 10, 10))


def test_qualified_dem_rejects_ambiguous_coordinates():
    pc = source()
    kwargs = dict(unit="m", vertical_reference="EPSG:5613")
    pc.transform.srs = ""
    with pytest.raises(ValueError, match="CRS"):
        builder.build_terrain_dem(pc, 1.0, **kwargs)
    pc.transform.srs = "EPSG:4326"
    with pytest.raises(ValueError, match="projected"):
        builder.build_terrain_dem(pc, 1.0, **kwargs)
    pc.transform.srs = "EPSG:7415"
    with pytest.raises(ValueError, match="vertical axis unit"):
        builder.build_terrain_dem(pc, 1.0, unit="cm", vertical_reference="EPSG:5709")
    pc.transform.srs = "EPSG:3006"
    pc.transform.set_translation(1, 0, 0)
    with pytest.raises(ValueError, match="transform"):
        builder.build_terrain_dem(pc, 1.0, **kwargs)


def test_dem_metadata_and_arrays_survive_native_package_and_reject_invalid_write(
    tmp_path,
):
    terrain = make_dem()
    path = tmp_path / "terrain.dtcc"
    io.save_model(terrain, path)
    original = path.read_bytes()
    restored = io.load_model(path)
    assert exchange.dumps(restored) == exchange.dumps(terrain)
    package = terrain.export(tmp_path / "terrain.dtccpkg", canonical=True)
    restored = load_model_package(package.path)
    assert exchange.dumps(restored) == exchange.dumps(terrain)
    assert restored.dataset_context == terrain.dataset_context
    terrain.attributes["elevation_rasters"][0]["unit"] = "furlong"
    with pytest.raises(ValueError, match="unit"):
        io.save_model(terrain, path)
    assert path.read_bytes() == original
