# Copyright(C) 2023 Anders Logg and Dag Wästberg
# Licensed under the MIT License

"""Tests for vector I/O utilities."""

import pytest
from pathlib import Path
import tempfile
import shapely.geometry

from dtcc_core.io.vector_utils import (
    validate_vector_file,
    get_vector_driver,
    create_bounds_filter,
    determine_io_crs,
    safe_reproject_geometry,
    VECTOR_DRIVERS,
)
from dtcc_core.model import Bounds


class TestValidateVectorFile:
    """Tests for validate_vector_file function."""

    def test_validate_existing_file(self, tmp_path):
        """Test validation of an existing file."""
        test_file = tmp_path / "test.shp"
        test_file.touch()

        result = validate_vector_file(test_file)

        assert result == test_file
        assert isinstance(result, Path)

    def test_validate_file_not_found(self, tmp_path):
        """Test validation raises FileNotFoundError for missing file."""
        nonexistent = tmp_path / "nonexistent.shp"

        with pytest.raises(FileNotFoundError) as exc_info:
            validate_vector_file(nonexistent)

        assert "not found" in str(exc_info.value).lower()
        assert str(nonexistent) in str(exc_info.value)

    def test_validate_string_path(self, tmp_path):
        """Test validation accepts string paths."""
        test_file = tmp_path / "test.shp"
        test_file.touch()

        result = validate_vector_file(str(test_file))

        assert result == test_file
        assert isinstance(result, Path)

    def test_validate_multiple_files_allowed(self, tmp_path):
        """Test validation of multiple files when allowed."""
        file1 = tmp_path / "test1.shp"
        file2 = tmp_path / "test2.shp"
        file1.touch()
        file2.touch()

        result = validate_vector_file([file1, file2], allow_multiple=True)

        assert len(result) == 2
        assert result[0] == file1
        assert result[1] == file2
        assert all(isinstance(p, Path) for p in result)

    def test_validate_multiple_files_not_allowed(self, tmp_path):
        """Test validation raises ValueError for multiple files when not allowed."""
        file1 = tmp_path / "test1.shp"
        file2 = tmp_path / "test2.shp"
        file1.touch()
        file2.touch()

        with pytest.raises(ValueError) as exc_info:
            validate_vector_file([file1, file2], allow_multiple=False)

        assert "multiple" in str(exc_info.value).lower()

    def test_validate_multiple_files_one_missing(self, tmp_path):
        """Test validation raises error if one of multiple files is missing."""
        file1 = tmp_path / "test1.shp"
        file2 = tmp_path / "test2.shp"
        file1.touch()
        # file2 not created

        with pytest.raises(FileNotFoundError):
            validate_vector_file([file1, file2], allow_multiple=True)

    def test_validate_empty_list(self):
        """Test validation of empty list returns empty list."""
        result = validate_vector_file([], allow_multiple=True)
        assert result == []

    def test_validate_tuple_input(self, tmp_path):
        """Test validation accepts tuples as input."""
        file1 = tmp_path / "test1.shp"
        file2 = tmp_path / "test2.shp"
        file1.touch()
        file2.touch()

        result = validate_vector_file((file1, file2), allow_multiple=True)

        assert len(result) == 2
        assert all(isinstance(p, Path) for p in result)


class TestGetVectorDriver:
    """Tests for get_vector_driver function."""

    def test_get_driver_shp(self):
        """Test driver mapping for .shp files."""
        driver = get_vector_driver("test.shp")
        assert driver == "ESRI Shapefile"

    def test_get_driver_geojson(self):
        """Test driver mapping for .geojson files."""
        driver = get_vector_driver("test.geojson")
        assert driver == "GeoJSON"

    def test_get_driver_json(self):
        """Test driver mapping for .json files."""
        driver = get_vector_driver("test.json")
        assert driver == "GeoJSON"

    def test_get_driver_gpkg(self):
        """Test driver mapping for .gpkg files."""
        driver = get_vector_driver("test.gpkg")
        assert driver == "GPKG"

    def test_get_driver_case_insensitive(self):
        """Test driver mapping is case-insensitive."""
        assert get_vector_driver("test.SHP") == "ESRI Shapefile"
        assert get_vector_driver("test.GeoJSON") == "GeoJSON"
        assert get_vector_driver("test.GPKG") == "GPKG"

    def test_get_driver_unsupported_format(self):
        """Test error for unsupported format."""
        with pytest.raises(ValueError) as exc_info:
            get_vector_driver("test.xyz")

        assert "unsupported" in str(exc_info.value).lower()
        assert ".xyz" in str(exc_info.value)
        assert "supported formats" in str(exc_info.value).lower()

    def test_get_driver_with_path_object(self):
        """Test driver mapping accepts Path objects."""
        driver = get_vector_driver(Path("test.shp"))
        assert driver == "ESRI Shapefile"

    def test_get_driver_complex_path(self):
        """Test driver mapping with complex file paths."""
        driver = get_vector_driver("/path/to/some/file.geojson")
        assert driver == "GeoJSON"

    def test_all_supported_formats_mapped(self):
        """Test that all formats in VECTOR_DRIVERS are accessible."""
        for ext, expected_driver in VECTOR_DRIVERS.items():
            driver = get_vector_driver(f"test{ext}")
            assert driver == expected_driver


class TestCreateBoundsFilter:
    """Tests for create_bounds_filter function."""

    def test_create_filter_none_bounds(self):
        """Test that None bounds returns None."""
        result = create_bounds_filter(None)
        assert result is None

    def test_create_filter_contains_strategy(self):
        """Test creating filter with contains strategy."""
        bounds = Bounds(0, 0, 100, 100)
        result = create_bounds_filter(bounds, strategy='contains')

        assert result is not None
        assert 'geometry' in result
        assert 'strategy' in result
        assert isinstance(result['geometry'], shapely.geometry.Polygon)
        assert callable(result['strategy'])

    def test_create_filter_intersects_strategy(self):
        """Test creating filter with intersects strategy."""
        bounds = Bounds(0, 0, 100, 100)
        result = create_bounds_filter(bounds, strategy='intersects')

        assert result is not None
        assert 'geometry' in result
        assert 'strategy' in result
        assert callable(result['strategy'])

    def test_create_filter_invalid_strategy(self):
        """Test error for invalid strategy."""
        bounds = Bounds(0, 0, 100, 100)

        with pytest.raises(ValueError) as exc_info:
            create_bounds_filter(bounds, strategy='invalid')

        assert "unknown strategy" in str(exc_info.value).lower()
        assert "invalid" in str(exc_info.value)

    def test_filter_geometry_bounds(self):
        """Test that filter geometry matches input bounds."""
        bounds = Bounds(10, 20, 110, 120)
        result = create_bounds_filter(bounds)

        geom = result['geometry']
        bbox = geom.bounds  # (minx, miny, maxx, maxy)

        assert bbox[0] == pytest.approx(10)
        assert bbox[1] == pytest.approx(20)
        assert bbox[2] == pytest.approx(110)
        assert bbox[3] == pytest.approx(120)

    def test_filter_positive_buffer(self):
        """Test positive buffer expands bounds."""
        bounds = Bounds(0, 0, 100, 100)
        result = create_bounds_filter(bounds, buffer=10)

        geom = result['geometry']
        bbox = geom.bounds

        # Buffer expands by 10 on all sides
        assert bbox[0] == pytest.approx(-10)
        assert bbox[1] == pytest.approx(-10)
        assert bbox[2] == pytest.approx(110)
        assert bbox[3] == pytest.approx(110)

    def test_filter_negative_buffer(self):
        """Test negative buffer contracts bounds."""
        bounds = Bounds(0, 0, 100, 100)
        result = create_bounds_filter(bounds, buffer=-10)

        geom = result['geometry']
        bbox = geom.bounds

        # Buffer contracts by 10 on all sides
        assert bbox[0] == pytest.approx(10)
        assert bbox[1] == pytest.approx(10)
        assert bbox[2] == pytest.approx(90)
        assert bbox[3] == pytest.approx(90)

    def test_filter_zero_buffer(self):
        """Test zero buffer doesn't change bounds."""
        bounds = Bounds(0, 0, 100, 100)
        result = create_bounds_filter(bounds, buffer=0)

        geom = result['geometry']
        bbox = geom.bounds

        assert bbox[0] == pytest.approx(0)
        assert bbox[1] == pytest.approx(0)
        assert bbox[2] == pytest.approx(100)
        assert bbox[3] == pytest.approx(100)

    def test_contains_strategy_behavior(self):
        """Test contains strategy actually tests containment."""
        bounds = Bounds(0, 0, 100, 100)
        filter_dict = create_bounds_filter(bounds, strategy='contains')

        filter_geom = filter_dict['geometry']
        strategy = filter_dict['strategy']

        # Geometry fully inside bounds
        inside = shapely.geometry.Point(50, 50).buffer(10)
        assert strategy(filter_geom, inside)

        # Geometry fully outside bounds
        outside = shapely.geometry.Point(200, 200).buffer(10)
        assert not strategy(filter_geom, outside)

        # Geometry partially overlapping (intersects but not contained)
        overlap = shapely.geometry.Point(95, 95).buffer(10)
        assert not strategy(filter_geom, overlap)

    def test_intersects_strategy_behavior(self):
        """Test intersects strategy tests for any overlap."""
        bounds = Bounds(0, 0, 100, 100)
        filter_dict = create_bounds_filter(bounds, strategy='intersects')

        filter_geom = filter_dict['geometry']
        strategy = filter_dict['strategy']

        # Geometry fully inside bounds
        inside = shapely.geometry.Point(50, 50).buffer(10)
        assert strategy(filter_geom, inside)

        # Geometry fully outside bounds
        outside = shapely.geometry.Point(200, 200).buffer(10)
        assert not strategy(filter_geom, outside)

        # Geometry partially overlapping
        overlap = shapely.geometry.Point(95, 95).buffer(10)
        assert strategy(filter_geom, overlap)

    def test_filter_with_linestring(self):
        """Test filter works with LineString geometries."""
        bounds = Bounds(0, 0, 100, 100)

        # Contains strategy
        filter_contains = create_bounds_filter(bounds, strategy='contains')
        line_inside = shapely.geometry.LineString([(10, 10), (90, 90)])
        assert filter_contains['strategy'](filter_contains['geometry'], line_inside)

        # Intersects strategy
        filter_intersects = create_bounds_filter(bounds, strategy='intersects')
        line_crossing = shapely.geometry.LineString([(50, 50), (150, 150)])
        assert filter_intersects['strategy'](filter_intersects['geometry'], line_crossing)
        assert not filter_contains['strategy'](filter_contains['geometry'], line_crossing)

    def test_filter_usage_pattern(self):
        """Test recommended usage pattern with None checking."""
        bounds = Bounds(0, 0, 100, 100)
        filter_dict = create_bounds_filter(bounds)

        # Pattern: if filter_dict and filter_dict['strategy'](...)
        test_geom = shapely.geometry.Point(50, 50).buffer(5)

        if filter_dict and filter_dict['strategy'](filter_dict['geometry'], test_geom):
            result = "passes"
        else:
            result = "fails"

        assert result == "passes"

    def test_filter_none_bounds_usage_pattern(self):
        """Test usage pattern with None bounds."""
        filter_dict = create_bounds_filter(None)
        test_geom = shapely.geometry.Point(50, 50).buffer(5)

        # Pattern should handle None gracefully
        if filter_dict and filter_dict['strategy'](filter_dict['geometry'], test_geom):
            result = "passes"
        else:
            result = "no filter"

        assert result == "no filter"


class TestDetermineIOCrs:
    """Tests for determine_io_crs function."""

    def test_load_without_target_crs(self):
        """Test load scenario without target CRS returns source CRS."""
        result = determine_io_crs("EPSG:3006", None, context="footprints")
        assert result == "EPSG:3006"

    def test_load_with_target_crs(self):
        """Test load scenario with target CRS validates and returns it."""
        result = determine_io_crs("EPSG:3006", "EPSG:4326", context="footprints")
        assert result == "EPSG:4326"

    def test_save_to_geojson_auto_wgs84(self):
        """Test save to GeoJSON automatically uses WGS84."""
        result = determine_io_crs(
            "EPSG:3006", None,
            output_filepath="data.geojson",
            context="trees"
        )
        assert result == "EPSG:4326"

    def test_save_to_geojson_json_extension(self):
        """Test save to .json file also uses WGS84."""
        result = determine_io_crs(
            "EPSG:3006", None,
            output_filepath="data.json",
            context="landuse"
        )
        assert result == "EPSG:4326"

    def test_save_to_shapefile_preserve_crs(self):
        """Test save to shapefile preserves source CRS."""
        result = determine_io_crs(
            "EPSG:3006", None,
            output_filepath="data.shp",
            context="footprints"
        )
        assert result == "EPSG:3006"

    def test_save_to_gpkg_preserve_crs(self):
        """Test save to GeoPackage preserves source CRS."""
        result = determine_io_crs(
            "EPSG:32633", None,
            output_filepath="data.gpkg",
            context="roadnetwork"
        )
        assert result == "EPSG:32633"

    def test_explicit_output_crs_overrides_format(self):
        """Test explicit output_crs overrides format requirements."""
        result = determine_io_crs(
            "EPSG:3006", "EPSG:32633",
            output_filepath="data.geojson",
            context="trees"
        )
        assert result == "EPSG:32633"

    def test_no_log_when_same_crs(self):
        """Test no reprojection when source equals target."""
        result = determine_io_crs(
            "EPSG:3006", "EPSG:3006",
            context="footprints"
        )
        assert result == "EPSG:3006"

    def test_log_reprojection_disabled(self):
        """Test log_reprojection=False parameter accepted."""
        result = determine_io_crs(
            "EPSG:3006", "EPSG:4326",
            context="trees",
            log_reprojection=False
        )
        assert result == "EPSG:4326"

    def test_context_parameter(self):
        """Test context parameter is accepted."""
        result = determine_io_crs("EPSG:3006", "EPSG:4326", context="road network")
        assert result == "EPSG:4326"

    def test_empty_context(self):
        """Test empty context works."""
        result = determine_io_crs("EPSG:3006", "EPSG:4326", context="")
        assert result == "EPSG:4326"

    def test_invalid_target_crs_raises_error(self):
        """Test invalid target CRS raises ValueError."""
        with pytest.raises(ValueError) as exc_info:
            determine_io_crs("EPSG:3006", "INVALID:CRS")
        assert "Invalid CRS" in str(exc_info.value)

    def test_path_object_accepted(self):
        """Test Path objects are accepted for output_filepath."""
        result = determine_io_crs(
            "EPSG:3006", None,
            output_filepath=Path("data.geojson")
        )
        assert result == "EPSG:4326"


class TestSafeReprojectGeometry:
    """Tests for safe_reproject_geometry function."""

    def test_single_geometry_reprojection(self):
        """Test reprojection of a single geometry."""
        point = shapely.geometry.Point(10, 60)  # Lon/lat in WGS84
        result = safe_reproject_geometry(point, "EPSG:4326", "EPSG:3857")

        assert isinstance(result, shapely.geometry.Point)
        # Coordinates should be different after reprojection
        assert result.x != point.x or result.y != point.y

    def test_list_of_geometries_reprojection(self):
        """Test reprojection of a list of geometries."""
        points = [
            shapely.geometry.Point(0, 0),
            shapely.geometry.Point(1, 1),
            shapely.geometry.Point(2, 2),
        ]
        result = safe_reproject_geometry(points, "EPSG:4326", "EPSG:3857")

        assert isinstance(result, list)
        assert len(result) == 3
        assert all(isinstance(p, shapely.geometry.Point) for p in result)

    def test_same_crs_returns_original(self):
        """Test same CRS returns original geometry unchanged."""
        point = shapely.geometry.Point(0, 0)
        result = safe_reproject_geometry(point, "EPSG:3006", "EPSG:3006")

        assert result is point  # Should be the exact same object

    def test_same_crs_list_returns_original(self):
        """Test same CRS returns original list unchanged."""
        points = [shapely.geometry.Point(0, 0), shapely.geometry.Point(1, 1)]
        result = safe_reproject_geometry(points, "EPSG:3006", "EPSG:3006")

        assert result is points  # Should be the exact same object

    def test_error_context_in_message(self):
        """Test error context parameter is accepted."""
        point = shapely.geometry.Point(0, 0)

        with pytest.raises(RuntimeError):
            safe_reproject_geometry(
                point, "EPSG:3006", "INVALID:CRS",
                error_context="building 123"
            )

    def test_raise_on_error_true(self):
        """Test raise_on_error=True raises RuntimeError on failure."""
        point = shapely.geometry.Point(0, 0)

        with pytest.raises(RuntimeError):
            safe_reproject_geometry(point, "EPSG:3006", "INVALID:CRS")

    def test_raise_on_error_false(self):
        """Test raise_on_error=False returns original on failure."""
        point = shapely.geometry.Point(0, 0)

        result = safe_reproject_geometry(
            point, "EPSG:3006", "INVALID:CRS",
            raise_on_error=False
        )

        assert result is point  # Should return original

    def test_raise_on_error_false_list(self):
        """Test raise_on_error=False returns original list on failure."""
        points = [shapely.geometry.Point(0, 0), shapely.geometry.Point(1, 1)]

        result = safe_reproject_geometry(
            points, "EPSG:3006", "INVALID:CRS",
            raise_on_error=False
        )

        assert result is points  # Should return original

    def test_linestring_reprojection(self):
        """Test reprojection works with LineString."""
        line = shapely.geometry.LineString([(0, 0), (1, 1)])
        result = safe_reproject_geometry(line, "EPSG:4326", "EPSG:3857")

        assert isinstance(result, shapely.geometry.LineString)

    def test_polygon_reprojection(self):
        """Test reprojection works with Polygon."""
        polygon = shapely.geometry.Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
        result = safe_reproject_geometry(polygon, "EPSG:4326", "EPSG:3857")

        assert isinstance(result, shapely.geometry.Polygon)

    def test_multipolygon_reprojection(self):
        """Test reprojection works with MultiPolygon."""
        poly1 = shapely.geometry.Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
        poly2 = shapely.geometry.Polygon([(2, 2), (3, 2), (3, 3), (2, 3)])
        multipoly = shapely.geometry.MultiPolygon([poly1, poly2])

        result = safe_reproject_geometry(multipoly, "EPSG:4326", "EPSG:3857")

        assert isinstance(result, shapely.geometry.MultiPolygon)

    def test_error_context_road_network(self):
        """Test error context for road network use case."""
        lines = [shapely.geometry.LineString([(0, 0), (1, 1)])]

        with pytest.raises(RuntimeError):
            safe_reproject_geometry(
                lines, "EPSG:3006", "INVALID:CRS",
                error_context="road network"
            )
