from unittest.mock import ANY, patch

import pytest

from dtcc_core.io.data import wrapper
from dtcc_core.model import Bounds


def test_download_data_propagates_lidar_request_error():
    bounds = Bounds(xmin=0, ymin=0, xmax=10, ymax=10)

    with patch.object(
        wrapper,
        "download_lidar",
        side_effect=RuntimeError(
            'Request failed with status 404:\n{"detail":"No lidar tiles intersect the requested bounding box."}'
        ),
    ):
        with pytest.raises(RuntimeError, match="No lidar tiles intersect"):
            wrapper.download_data("lidar", "dtcc", bounds)


def test_download_data_rejects_empty_lidar_download_result():
    bounds = Bounds(xmin=0, ymin=0, xmax=10, ymax=10)

    with patch.object(wrapper, "download_lidar", return_value=None), patch.object(
        wrapper.io, "load_pointcloud"
    ) as load_pointcloud:
        with pytest.raises(
            RuntimeError,
            match="No lidar data available for the requested bounding box.",
        ):
            wrapper.download_data("lidar", "dtcc", bounds)

    load_pointcloud.assert_not_called()


def test_download_data_uses_cached_footprints_before_network():
    bounds = Bounds(xmin=0, ymin=0, xmax=10, ymax=10)
    cached_buildings = ["cached-building"]

    with patch.object(
        wrapper,
        "_load_cached_footprints",
        return_value=cached_buildings,
    ), patch.object(wrapper, "download_tiles") as download_tiles:
        result = wrapper.download_data("footprints", "dtcc", bounds)

    assert result == cached_buildings
    download_tiles.assert_not_called()


def test_download_data_falls_back_to_download_when_cache_misses():
    bounds = Bounds(xmin=0, ymin=0, xmax=10, ymax=10)

    with patch.object(
        wrapper,
        "_load_cached_footprints",
        return_value=None,
    ), patch.object(
        wrapper,
        "download_tiles",
        return_value=["/tmp/tile_0_0.gpkg"],
    ) as download_tiles, patch.object(
        wrapper.io,
        "load_footprints",
        return_value=["downloaded-building"],
    ) as load_footprints:
        result = wrapper.download_data("footprints", "dtcc", bounds)

    assert result == ["downloaded-building"]
    download_tiles.assert_called_once_with(
        bounds.tuple,
        ANY,
        server_url=wrapper.DTCC_GPKG_URL,
    )
    load_footprints.assert_called_once_with(["/tmp/tile_0_0.gpkg"], bounds=bounds)


def test_load_cached_footprints_accepts_empty_result_when_tile_covers_bounds(tmp_path):
    cache_dir = tmp_path / "downloaded-gpkg"
    cache_dir.mkdir()
    cached_tile = cache_dir / "tile_0_0.gpkg"
    cached_tile.touch()
    bounds = Bounds(xmin=100, ymin=100, xmax=200, ymax=200)

    with patch.object(wrapper, "GPKG_CACHE_DIR", str(tmp_path)), patch.object(
        wrapper.io,
        "load_footprints",
        return_value=[],
    ) as load_footprints:
        result = wrapper._load_cached_footprints(bounds)

    assert result == []
    load_footprints.assert_called_once_with([str(cached_tile)], bounds=bounds)


def test_load_cached_footprints_falls_back_when_empty_cache_does_not_cover_bounds(
    tmp_path,
):
    cache_dir = tmp_path / "downloaded-gpkg"
    cache_dir.mkdir()
    cached_tile = cache_dir / "tile_0_0.gpkg"
    cached_tile.touch()
    bounds = Bounds(xmin=9900, ymin=9900, xmax=10100, ymax=10100)

    with patch.object(wrapper, "GPKG_CACHE_DIR", str(tmp_path)), patch.object(
        wrapper.io,
        "load_footprints",
        return_value=[],
    ) as load_footprints:
        result = wrapper._load_cached_footprints(bounds)

    assert result is None
    load_footprints.assert_called_once_with([str(cached_tile)], bounds=bounds)


def test_load_cached_footprints_accepts_empty_result_when_adjacent_tiles_cover_bounds(
    tmp_path,
):
    cache_dir = tmp_path / "downloaded-gpkg"
    cache_dir.mkdir()
    left_tile = cache_dir / "tile_0_0.gpkg"
    right_tile = cache_dir / "tile_10000_0.gpkg"
    left_tile.touch()
    right_tile.touch()
    bounds = Bounds(xmin=9900, ymin=100, xmax=10100, ymax=200)

    with patch.object(wrapper, "GPKG_CACHE_DIR", str(tmp_path)), patch.object(
        wrapper.io,
        "load_footprints",
        return_value=[],
    ) as load_footprints:
        result = wrapper._load_cached_footprints(bounds)

    assert result == []
    load_footprints.assert_called_once_with(
        [str(left_tile), str(right_tile)],
        bounds=bounds,
    )
