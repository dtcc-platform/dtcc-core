# Copyright(C) 2024 Dag Wästberg
# Licensed under the MIT License

"""Tests for CRS handling in vector I/O."""

import pytest
import numpy as np
import tempfile
import fiona
import shapely.geometry
from pathlib import Path

from dtcc_core.io.vector_utils import (
    validate_crs,
    get_format_required_crs,
    reproject_shapely_geometry,
    get_geometry_crs,
    set_geometry_crs,
)
from dtcc_core.model import Building, City, Tree, GeometryType
from dtcc_core.model.geometry import Surface


class TestValidateCrs:
    """Tests for validate_crs function."""

    def test_validate_valid_epsg(self):
        """Test validation of valid EPSG codes."""
        assert validate_crs("EPSG:4326") == "EPSG:4326"
        assert validate_crs("EPSG:3006") == "EPSG:3006"
        assert validate_crs("EPSG:3857") == "EPSG:3857"

    def test_validate_none(self):
        """Test validation of None returns None."""
        assert validate_crs(None) is None

    def test_validate_empty_string(self):
        """Test validation of empty string returns None."""
        assert validate_crs("") is None

    def test_validate_invalid_crs(self):
        """Test validation of invalid CRS raises ValueError."""
        with pytest.raises(ValueError, match="Invalid CRS"):
            validate_crs("invalid_crs")

    def test_validate_invalid_epsg(self):
        """Test validation of invalid EPSG code raises ValueError."""
        with pytest.raises(ValueError, match="Invalid CRS"):
            validate_crs("EPSG:99999999")


class TestGetFormatRequiredCrs:
    """Tests for get_format_required_crs function."""

    def test_geojson_extension(self):
        """Test GeoJSON requires EPSG:4326."""
        assert get_format_required_crs("data.geojson") == "EPSG:4326"

    def test_json_extension(self):
        """Test .json extension requires EPSG:4326."""
        assert get_format_required_crs("data.json") == "EPSG:4326"

    def test_shapefile_no_requirement(self):
        """Test shapefile has no CRS requirement."""
        assert get_format_required_crs("data.shp") is None

    def test_gpkg_no_requirement(self):
        """Test GeoPackage has no CRS requirement."""
        assert get_format_required_crs("data.gpkg") is None

    def test_case_insensitive(self):
        """Test extension matching is case-insensitive."""
        assert get_format_required_crs("data.GEOJSON") == "EPSG:4326"
        assert get_format_required_crs("data.GeoJson") == "EPSG:4326"

    def test_path_object(self):
        """Test works with Path objects."""
        assert get_format_required_crs(Path("data.geojson")) == "EPSG:4326"
        assert get_format_required_crs(Path("data.shp")) is None


class TestReprojectShapelyGeometry:
    """Tests for reproject_shapely_geometry function."""

    def test_reproject_point_wgs84_to_sweref99(self):
        """Test reprojecting a point from WGS84 to SWEREF99 TM."""
        # Gothenburg approximate location
        point_wgs84 = shapely.geometry.Point(11.9746, 57.7089)
        point_sweref = reproject_shapely_geometry(
            point_wgs84, "EPSG:4326", "EPSG:3006"
        )

        # Expected approximate coordinates in SWEREF99 TM
        # (these are approximate, so we use larger tolerance)
        assert abs(point_sweref.x - 319180) < 1000  # ~319180
        assert abs(point_sweref.y - 6399862) < 1000  # ~6399862

    def test_reproject_point_sweref99_to_wgs84(self):
        """Test reprojecting a point from SWEREF99 TM to WGS84."""
        point_sweref = shapely.geometry.Point(319180, 6399862)
        point_wgs84 = reproject_shapely_geometry(
            point_sweref, "EPSG:3006", "EPSG:4326"
        )

        # Expected approximate coordinates in WGS84
        assert abs(point_wgs84.x - 11.9746) < 0.01
        assert abs(point_wgs84.y - 57.7089) < 0.01

    def test_reproject_polygon(self):
        """Test reprojecting a polygon geometry."""
        polygon_wgs84 = shapely.geometry.box(11.97, 57.70, 11.98, 57.71)
        polygon_sweref = reproject_shapely_geometry(
            polygon_wgs84, "EPSG:4326", "EPSG:3006"
        )

        assert polygon_sweref.geom_type == "Polygon"
        assert polygon_sweref.is_valid
        # Check that coordinates are in SWEREF99 range (6-7 digits)
        coords = list(polygon_sweref.exterior.coords)
        assert all(300000 < x < 400000 for x, y in coords)
        assert all(6300000 < y < 6500000 for x, y in coords)

    def test_reproject_linestring(self):
        """Test reprojecting a linestring geometry."""
        line_wgs84 = shapely.geometry.LineString([
            (11.97, 57.70),
            (11.98, 57.71),
            (11.99, 57.72)
        ])
        line_sweref = reproject_shapely_geometry(
            line_wgs84, "EPSG:4326", "EPSG:3006"
        )

        assert line_sweref.geom_type == "LineString"
        assert len(list(line_sweref.coords)) == 3

    def test_reproject_same_crs_no_change(self):
        """Test reprojecting with same source and target CRS returns unchanged."""
        point = shapely.geometry.Point(11.9746, 57.7089)
        result = reproject_shapely_geometry(point, "EPSG:4326", "EPSG:4326")

        assert result.x == point.x
        assert result.y == point.y

    def test_reproject_multipolygon(self):
        """Test reprojecting a multipolygon geometry."""
        poly1 = shapely.geometry.box(11.97, 57.70, 11.975, 57.705)
        poly2 = shapely.geometry.box(11.98, 57.71, 11.985, 57.715)
        multipoly_wgs84 = shapely.geometry.MultiPolygon([poly1, poly2])

        multipoly_sweref = reproject_shapely_geometry(
            multipoly_wgs84, "EPSG:4326", "EPSG:3006"
        )

        assert multipoly_sweref.geom_type == "MultiPolygon"
        assert len(multipoly_sweref.geoms) == 2


class TestGetGeometryCrs:
    """Tests for get_geometry_crs function."""

    def test_get_crs_from_building(self):
        """Test getting CRS from Building object."""
        building = Building()
        surface = Surface()
        surface.transform.srs = "EPSG:4326"
        building.add_geometry(surface, GeometryType.LOD0)

        assert get_geometry_crs(building) == "EPSG:4326"

    def test_get_crs_from_tree(self):
        """Test getting CRS from Tree object with direct transform."""
        tree = Tree(id=1, position=[100, 200, 0], height=10)
        tree.transform.srs = "EPSG:3006"

        assert get_geometry_crs(tree) == "EPSG:3006"

    def test_get_crs_uses_fallback(self):
        """Test getting CRS uses fallback when not set."""
        building = Building()
        surface = Surface()
        building.add_geometry(surface, GeometryType.LOD0)

        assert get_geometry_crs(building) == "EPSG:3006"  # Default fallback

    def test_get_crs_custom_fallback(self):
        """Test getting CRS with custom fallback."""
        building = Building()
        surface = Surface()
        building.add_geometry(surface, GeometryType.LOD0)

        assert get_geometry_crs(building, fallback="EPSG:4326") == "EPSG:4326"

    def test_get_crs_from_first_geometry(self):
        """Test getting CRS from first geometry in dict."""
        building = Building()
        surface1 = Surface()
        surface1.transform.srs = "EPSG:4326"
        surface2 = Surface()
        surface2.transform.srs = "EPSG:3006"

        building.add_geometry(surface1, GeometryType.LOD0)
        building.add_geometry(surface2, GeometryType.LOD1)

        # Should return from first geometry with CRS set
        assert get_geometry_crs(building) == "EPSG:4326"


class TestSetGeometryCrs:
    """Tests for set_geometry_crs function."""

    def test_set_crs_on_building(self):
        """Test setting CRS on Building object."""
        building = Building()
        surface = Surface()
        building.add_geometry(surface, GeometryType.LOD0)

        set_geometry_crs(building, "EPSG:4326")

        assert building.geometry[GeometryType.LOD0].transform.srs == "EPSG:4326"

    def test_set_crs_on_tree(self):
        """Test setting CRS on Tree object."""
        tree = Tree(id=1, position=[100, 200, 0], height=10)

        set_geometry_crs(tree, "EPSG:3006")

        assert tree.transform.srs == "EPSG:3006"

    def test_set_crs_on_multiple_geometries(self):
        """Test setting CRS on object with multiple geometry types."""
        building = Building()
        surface1 = Surface()
        surface2 = Surface()
        building.add_geometry(surface1, GeometryType.LOD0)
        building.add_geometry(surface2, GeometryType.LOD1)

        set_geometry_crs(building, "EPSG:4326")

        assert building.geometry[GeometryType.LOD0].transform.srs == "EPSG:4326"
        assert building.geometry[GeometryType.LOD1].transform.srs == "EPSG:4326"

    def test_set_crs_overwrites_existing(self):
        """Test setting CRS overwrites existing value."""
        building = Building()
        surface = Surface()
        surface.transform.srs = "EPSG:3006"
        building.add_geometry(surface, GeometryType.LOD0)

        set_geometry_crs(building, "EPSG:4326")

        assert building.geometry[GeometryType.LOD0].transform.srs == "EPSG:4326"


class TestLoadWithTargetCrs:
    """Integration tests for loading vector data with target CRS."""

    def create_test_shapefile(self, filepath, crs="EPSG:3006"):
        """Helper to create a test shapefile."""
        schema = {
            "geometry": "Polygon",
            "properties": {"id": "str", "name": "str"}
        }

        # Create a simple polygon in the specified CRS
        if crs == "EPSG:3006":
            # SWEREF99 TM coordinates (Gothenburg area)
            polygon = shapely.geometry.box(319000, 6399000, 319100, 6399100)
        else:  # WGS84
            # WGS84 coordinates (Gothenburg area)
            polygon = shapely.geometry.box(11.97, 57.70, 11.98, 57.71)

        with fiona.open(
            filepath, "w", driver="ESRI Shapefile", crs=crs, schema=schema
        ) as dst:
            dst.write({
                "geometry": shapely.geometry.mapping(polygon),
                "properties": {"id": "test1", "name": "Test Building"}
            })

    def test_load_footprints_native_crs(self):
        """Test loading footprints without CRS conversion."""
        from dtcc_core.io.footprints import load

        with tempfile.TemporaryDirectory() as tmpdir:
            test_file = Path(tmpdir) / "buildings.shp"
            self.create_test_shapefile(test_file, crs="EPSG:3006")

            buildings = load(test_file)

            assert len(buildings) == 1
            assert get_geometry_crs(buildings[0]).lower() == "epsg:3006"

            # Verify coordinates are in SWEREF99 range
            surface = buildings[0].geometry[GeometryType.LOD0]
            coords = surface.vertices
            assert np.all(coords[:, 0] > 300000)
            assert np.all(coords[:, 1] > 6000000)

    def test_load_footprints_with_reprojection(self):
        """Test loading footprints with CRS reprojection."""
        from dtcc_core.io.footprints import load

        with tempfile.TemporaryDirectory() as tmpdir:
            test_file = Path(tmpdir) / "buildings.shp"
            self.create_test_shapefile(test_file, crs="EPSG:3006")

            buildings = load(test_file, target_crs="EPSG:4326")

            assert len(buildings) == 1
            assert get_geometry_crs(buildings[0]).upper() == "EPSG:4326"

            # Verify coordinates are in WGS84 range
            surface = buildings[0].geometry[GeometryType.LOD0]
            coords = surface.vertices
            assert np.all(np.abs(coords[:, 0]) < 180)  # Longitude
            assert np.all(np.abs(coords[:, 1]) < 90)   # Latitude


class TestSaveWithOutputCrs:
    """Integration tests for saving vector data with output CRS."""

    def test_save_footprints_geojson_auto_wgs84(self):
        """Test saving footprints to GeoJSON auto-converts to WGS84."""
        from dtcc_core.io.footprints import save

        # Create a building in SWEREF99 TM
        building = Building()
        surface = Surface()
        surface.vertices = np.array([
            [319000, 6399000, 0],
            [319100, 6399000, 0],
            [319100, 6399100, 0],
            [319000, 6399100, 0],
        ])
        surface.transform.srs = "EPSG:3006"
        building.add_geometry(surface, GeometryType.LOD0)
        building.id = "test1"

        city = City()
        city.add_buildings([building])

        with tempfile.TemporaryDirectory() as tmpdir:
            out_file = Path(tmpdir) / "buildings.geojson"
            save(city, out_file)

            # Verify file was created and has WGS84 CRS
            assert out_file.exists()

            with fiona.open(out_file) as src:
                # Check CRS is WGS84
                assert src.crs.to_epsg() == 4326

                # Check coordinates were transformed
                feature = next(iter(src))
                coords = feature["geometry"]["coordinates"][0]
                # All coordinates should be in WGS84 range
                for lon, lat in coords:
                    assert -180 <= lon <= 180
                    assert -90 <= lat <= 90

    def test_save_footprints_shapefile_preserves_crs(self):
        """Test saving footprints to shapefile preserves original CRS."""
        from dtcc_core.io.footprints import save

        # Create a building in SWEREF99 TM
        building = Building()
        surface = Surface()
        surface.vertices = np.array([
            [319000, 6399000, 0],
            [319100, 6399000, 0],
            [319100, 6399100, 0],
            [319000, 6399100, 0],
        ])
        surface.transform.srs = "EPSG:3006"
        building.add_geometry(surface, GeometryType.LOD0)
        building.id = "test1"

        city = City()
        city.add_buildings([building])

        with tempfile.TemporaryDirectory() as tmpdir:
            out_file = Path(tmpdir) / "buildings.shp"
            save(city, out_file)

            # Verify file was created and has correct CRS
            assert out_file.exists()

            with fiona.open(out_file) as src:
                # Check CRS is SWEREF99 TM
                assert src.crs.to_epsg() == 3006

    def test_save_footprints_with_explicit_output_crs(self):
        """Test saving footprints with explicit output CRS."""
        from dtcc_core.io.footprints import save

        # Create a building in SWEREF99 TM
        building = Building()
        surface = Surface()
        surface.vertices = np.array([
            [319000, 6399000, 0],
            [319100, 6399000, 0],
            [319100, 6399100, 0],
            [319000, 6399100, 0],
        ])
        surface.transform.srs = "EPSG:3006"
        building.add_geometry(surface, GeometryType.LOD0)
        building.id = "test1"

        city = City()
        city.add_buildings([building])

        with tempfile.TemporaryDirectory() as tmpdir:
            out_file = Path(tmpdir) / "buildings.shp"
            save(city, out_file, output_crs="EPSG:4326")

            # Verify file has WGS84 CRS
            with fiona.open(out_file) as src:
                assert src.crs.to_epsg() == 4326

                # Verify coordinates were transformed
                feature = next(iter(src))
                coords = feature["geometry"]["coordinates"][0]
                for lon, lat in coords:
                    assert -180 <= lon <= 180
                    assert -90 <= lat <= 90


class TestRoundtrip:
    """Tests for roundtrip load-save-load cycles."""

    def test_roundtrip_footprints_with_reprojection(self):
        """Test roundtrip preserves data through reprojection."""
        from dtcc_core.io.footprints import load, save

        with tempfile.TemporaryDirectory() as tmpdir:
            # Create original file in SWEREF99
            original_file = Path(tmpdir) / "original.shp"
            schema = {
                "geometry": "Polygon",
                "properties": {"id": "str"}
            }
            polygon = shapely.geometry.box(319000, 6399000, 319100, 6399100)

            with fiona.open(
                original_file, "w", driver="ESRI Shapefile",
                crs="EPSG:3006", schema=schema
            ) as dst:
                dst.write({
                    "geometry": shapely.geometry.mapping(polygon),
                    "properties": {"id": "test1"}
                })

            # Load original in native CRS for comparison
            buildings_original = load(original_file)

            # Load and save with reprojection
            buildings_wgs84 = load(original_file, target_crs="EPSG:4326")

            city = City()
            city.add_buildings(buildings_wgs84)

            output_file = Path(tmpdir) / "output.shp"
            save(city, output_file, output_crs="EPSG:3006")

            # Load again
            buildings_reloaded = load(output_file)

            # Verify CRS is correct
            assert get_geometry_crs(buildings_reloaded[0]).lower() == "epsg:3006"

            # Verify coordinates are approximately the same
            # Note: Comparing bounds since vertex count might differ (closed vs open rings)
            original_surface = buildings_original[0].geometry[GeometryType.LOD0]
            reloaded_surface = buildings_reloaded[0].geometry[GeometryType.LOD0]

            original_bounds = (
                np.min(original_surface.vertices[:, 0]),
                np.min(original_surface.vertices[:, 1]),
                np.max(original_surface.vertices[:, 0]),
                np.max(original_surface.vertices[:, 1]),
            )
            reloaded_bounds = (
                np.min(reloaded_surface.vertices[:, 0]),
                np.min(reloaded_surface.vertices[:, 1]),
                np.max(reloaded_surface.vertices[:, 0]),
                np.max(reloaded_surface.vertices[:, 1]),
            )

            # After roundtrip through WGS84, bounds should be close
            # but not exact due to floating point and projection
            assert np.allclose(original_bounds, reloaded_bounds, rtol=1e-5, atol=100)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
