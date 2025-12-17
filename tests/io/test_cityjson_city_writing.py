"""
Unit tests for writing City objects to CityJSON format.

Tests the high-level functionality of converting DTCC City objects
to CityJSON format, including buildings, terrain, and file I/O.
"""

import pytest
import json
import numpy as np
import tempfile
from pathlib import Path

from dtcc_core.model import City, Surface, MultiSurface, Mesh, Building, BuildingPart
from dtcc_core.model.object.object import GeometryType
from dtcc_core.model.object.terrain import Terrain
from dtcc_core.model.geometry import Bounds

from dtcc_core.io.cityjson.write_cityjson import (
    to_cityjson,
    save,
    CityJSONConfig,
)


class TestCityToCityJSON:
    """Test conversion of City objects to CityJSON format."""

    def test_empty_city_conversion(self):
        """Test conversion of an empty city."""
        city = City()
        result = to_cityjson(city)

        assert result["type"] == "CityJSON"
        assert result["version"] == "2.0"
        assert "transform" in result
        assert result["CityObjects"] == {}
        assert result["vertices"] == []
        assert result["transform"]["scale"] == [0.001, 0.001, 0.001]

    def test_city_with_bounds(self):
        """Test conversion of city with geographical bounds."""
        city = City()
        bounds = Bounds(
            xmin=10.0, ymin=20.0, zmin=0.0, xmax=110.0, ymax=120.0, zmax=50.0
        )
        city.bounds = bounds

        result = to_cityjson(city)

        assert "metadata" in result
        assert "geographicalExtent" in result["metadata"]
        expected_extent = [10.0, 20.0, 0.0, 110.0, 120.0, 50.0]
        assert result["metadata"]["geographicalExtent"] == expected_extent

    def test_city_with_single_building(self):
        """Test conversion of city with a single building."""
        city = City()

        # Create a building with a simple mesh
        building = Building()
        building.id = "building_001"
        building.attributes = {"height": 25.0, "function": "residential"}

        # Create a simple box mesh (6 faces)
        vertices = np.array(
            [
                [0, 0, 0],
                [10, 0, 0],
                [10, 10, 0],
                [0, 10, 0],  # Base
                [0, 0, 25],
                [10, 0, 25],
                [10, 10, 25],
                [0, 10, 25],  # Top
            ]
        )

        faces = np.array(
            [
                [0, 1, 2],
                [0, 2, 3],  # Bottom
                [4, 7, 6],
                [4, 6, 5],  # Top
                [0, 4, 5],
                [0, 5, 1],  # Front
                [2, 6, 7],
                [2, 7, 3],  # Back
                [0, 3, 7],
                [0, 7, 4],  # Left
                [1, 5, 6],
                [1, 6, 2],  # Right
            ]
        )

        mesh = Mesh(vertices=vertices, faces=faces)
        building.add_geometry(mesh, GeometryType.LOD2)
        city.add_building(building)

        result = to_cityjson(city)

        # Verify structure
        assert "building_001" in result["CityObjects"]
        building_data = result["CityObjects"]["building_001"]

        assert building_data["type"] == "Building"
        assert building_data["attributes"]["height"] == 25.0
        assert building_data["attributes"]["function"] == "residential"
        assert len(building_data["geometry"]) == 1

        # Check geometry
        geometry = building_data["geometry"][0]
        assert geometry["type"] == "MultiSurface"
        assert geometry["lod"] == 2.0
        assert len(geometry["boundaries"]) == 12  # 12 faces from the box
        assert len(result["vertices"]) == 8  # 8 vertices from the box

    def test_city_with_building_parts(self):
        """Test conversion of building with building parts."""
        city = City()

        # Create main building
        building = Building()
        building.id = "building_001"

        # Create building part
        building_part = BuildingPart()
        building_part.id = "part_001"
        building_part.attributes = {"usage": "residential"}

        # Add geometry to building part
        surface = Surface(
            vertices=np.array([[0, 0, 0], [5, 0, 0], [5, 5, 0], [0, 5, 0]])
        )
        building_part.add_geometry(surface, GeometryType.LOD1)

        building.add_child(building_part)
        city.add_building(building)

        result = to_cityjson(city)

        # Verify structure
        assert "building_001" in result["CityObjects"]
        assert "part_001" in result["CityObjects"]

        building_data = result["CityObjects"]["building_001"]
        part_data = result["CityObjects"]["part_001"]

        assert "children" in building_data
        assert "part_001" in building_data["children"]
        assert part_data["parents"] == ["building_001"]
        assert len(part_data["geometry"]) == 1
        assert part_data["geometry"][0]["lod"] == 1.0

    def test_city_with_terrain(self):
        """Test conversion of city with terrain."""
        city = City()

        # Create terrain
        terrain = Terrain()
        terrain.id = "terrain_001"
        terrain.attributes = {"type": "ground"}

        # Create terrain mesh
        vertices = np.array([[0, 0, 0], [100, 0, 2], [50, 100, 1], [0, 100, 0.5]])
        faces = np.array([[0, 1, 2], [0, 2, 3]])
        mesh = Mesh(vertices=vertices, faces=faces)
        terrain.add_geometry(mesh, GeometryType.MESH)

        city.add_child(terrain)

        result = to_cityjson(city)

        # Verify terrain conversion
        assert "terrain_001" in result["CityObjects"]
        terrain_data = result["CityObjects"]["terrain_001"]

        assert terrain_data["type"] == "TINRelief"
        assert terrain_data["attributes"]["type"] == "ground"
        assert len(terrain_data["geometry"]) == 1

        # Check terrain uses CompositeSurface
        geometry = terrain_data["geometry"][0]
        assert geometry["type"] == "CompositeSurface"
        assert len(geometry["boundaries"]) == 2  # 2 faces

        # Check terrain uses GroundSurface semantics
        assert "semantics" in geometry
        semantics = geometry["semantics"]
        assert len(semantics["surfaces"]) == 2
        for surface in semantics["surfaces"]:
            assert surface["type"] == "GroundSurface"

    def test_city_with_multiple_buildings(self):
        """Test conversion of city with multiple buildings."""
        city = City()

        # Create two buildings
        for i in range(2):
            building = Building()
            building.id = f"building_{i:03d}"
            building.attributes = {"building_id": i}

            # Simple surface geometry
            surface = Surface(
                vertices=np.array(
                    [
                        [i * 20, 0, 0],
                        [i * 20 + 10, 0, 0],
                        [i * 20 + 10, 10, 0],
                        [i * 20, 10, 0],
                    ]
                )
            )
            building.add_geometry(surface, GeometryType.LOD1)
            city.add_building(building)

        result = to_cityjson(city)

        # Verify both buildings are present
        assert "building_000" in result["CityObjects"]
        assert "building_001" in result["CityObjects"]
        assert len(result["CityObjects"]) == 2
        assert len(result["vertices"]) == 8  # 4 vertices per building

    def test_coordinate_scaling(self):
        """Test coordinate scaling functionality."""
        city = City()

        building = Building()
        # Use coordinates that will show scaling clearly
        surface = Surface(
            vertices=np.array([[1.5, 2.5, 3.5], [2.5, 2.5, 3.5], [2.5, 3.5, 3.5]])
        )
        building.add_geometry(surface, GeometryType.LOD1)
        city.add_building(building)

        # Test with default scale (0.001)
        result = to_cityjson(city, scale=0.001)

        # Vertices should be scaled by 1/0.001 = 1000 and rounded to integers
        expected_vertices = [
            [1500, 2500, 3500],  # 1.5*1000, 2.5*1000, 3.5*1000
            [2500, 2500, 3500],  # 2.5*1000, 2.5*1000, 3.5*1000
            [2500, 3500, 3500],  # 2.5*1000, 3.5*1000, 3.5*1000
        ]
        assert result["vertices"] == expected_vertices

        # Check transform
        assert result["transform"]["scale"] == [0.001, 0.001, 0.001]

    def test_custom_config(self):
        """Test conversion with custom configuration."""
        city = City()

        building = Building()
        surface = Surface(
            vertices=np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]])
        )
        building.add_geometry(surface, GeometryType.LOD1)
        city.add_building(building)

        # Custom configuration
        config = CityJSONConfig(
            semantic_types={"surface": "CustomWallSurface"},
            lod_mapping={GeometryType.LOD1: 1.5},
        )

        result = to_cityjson(city, config=config)

        # Check custom semantic type
        building_data = list(result["CityObjects"].values())[0]
        geometry = building_data["geometry"][0]
        assert geometry["semantics"]["surfaces"][0]["type"] == "CustomWallSurface"

        # Check custom LOD mapping
        assert geometry["lod"] == 1.5


class TestFileOperations:
    """Test file saving and loading operations."""

    def test_save_to_json_file(self):
        """Test saving city to JSON file."""
        city = City()
        building = Building()
        building.attributes = {"test": "value"}
        surface = Surface(vertices=np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0]]))
        building.add_geometry(surface, GeometryType.LOD1)
        city.add_building(building)

        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as tmp:
            tmp_path = Path(tmp.name)

        try:
            # Save the city
            save(city, tmp_path)

            # Verify file was created and contains valid JSON
            assert tmp_path.exists()

            with open(tmp_path, "r") as f:
                loaded_data = json.load(f)

            assert loaded_data["type"] == "CityJSON"
            assert loaded_data["version"] == "2.0"
            assert len(loaded_data["CityObjects"]) == 1
        except Exception as e:
            pytest.fail(f"Saving to JSON file failed: {e}")

        finally:
            # Clean up
            if tmp_path.exists():
                tmp_path.unlink()

    def test_save_unsupported_format(self):
        """Test saving with unsupported file extension."""
        city = City()

        with tempfile.NamedTemporaryFile(suffix=".txt") as tmp:
            tmp_path = Path(tmp.name)

            with pytest.raises(ValueError, match="Unsupported file format"):
                save(city, tmp_path)

    def test_save_with_custom_config(self):
        """Test saving with custom configuration."""
        city = City()
        building = Building()
        surface = Surface(vertices=np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0]]))
        building.add_geometry(surface, GeometryType.LOD1)
        city.add_building(building)

        config = CityJSONConfig(semantic_types={"surface": "CustomType"})

        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as tmp:
            tmp_path = Path(tmp.name)

        try:
            save(city, tmp_path, config=config)

            with open(tmp_path, "r") as f:
                loaded_data = json.load(f)

            # Check that custom config was applied
            building_data = list(loaded_data["CityObjects"].values())[0]
            geometry = building_data["geometry"][0]
            assert geometry["semantics"]["surfaces"][0]["type"] == "CustomType"

        finally:
            if tmp_path.exists():
                tmp_path.unlink()


class TestEdgeCases:
    """Test edge cases and error conditions."""

    def test_city_with_empty_geometry(self):
        """Test handling of buildings with empty geometry."""
        city = City()
        building = Building()
        building.id = "empty_building"

        # Add empty mesh
        empty_mesh = Mesh()
        building.add_geometry(empty_mesh, GeometryType.MESH)
        city.add_building(building)

        result = to_cityjson(city)

        # Building should exist but have no geometry
        assert "empty_building" in result["CityObjects"]
        building_data = result["CityObjects"]["empty_building"]
        assert len(building_data["geometry"]) == 0

    def test_city_with_mixed_geometry_types(self):
        """Test building with multiple geometry types."""
        city = City()
        building = Building()

        # Add different geometry types
        mesh = Mesh(
            vertices=np.array([[0, 0, 0], [1, 0, 0], [0.5, 1, 0]]),
            faces=np.array([[0, 1, 2]]),
        )
        surface = Surface(vertices=np.array([[0, 0, 1], [1, 0, 1], [0.5, 1, 1]]))

        building.add_geometry(mesh, GeometryType.LOD1)
        building.add_geometry(surface, GeometryType.LOD2)
        city.add_building(building)

        result = to_cityjson(city)

        # Should have both geometries
        building_data = list(result["CityObjects"].values())[0]
        assert len(building_data["geometry"]) == 2

        # Check LOD values
        lods = [geom["lod"] for geom in building_data["geometry"]]
        assert 1.0 in lods
        assert 2.0 in lods

    def test_large_coordinate_values(self):
        """Test handling of large coordinate values."""
        city = City()
        building = Building()

        # Use large coordinates
        surface = Surface(
            vertices=np.array(
                [
                    [1000000, 2000000, 100],
                    [1000010, 2000000, 100],
                    [1000010, 2000010, 100],
                ]
            )
        )
        building.add_geometry(surface, GeometryType.LOD1)
        city.add_building(building)

        result = to_cityjson(city, scale=0.001)

        # Should handle large coordinates properly
        assert len(result["vertices"]) == 3
        # Coordinates should be scaled and rounded to integers
        for vertex in result["vertices"]:
            assert all(isinstance(coord, int) for coord in vertex)

    def test_city_with_null_attributes(self):
        """Test handling of None/null attribute values."""
        city = City()
        building = Building()
        building.id = "test_building"
        building.attributes = {"name": None, "height": 10.0, "description": None}
        city.add_building(building)

        result = to_cityjson(city)

        # Should handle None values gracefully
        building_data = result["CityObjects"]["test_building"]
        assert building_data["attributes"]["name"] is None
        assert building_data["attributes"]["height"] == 10.0
        assert building_data["attributes"]["description"] is None
