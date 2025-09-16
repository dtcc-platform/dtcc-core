import pytest
import numpy as np
import json
from pathlib import Path
import tempfile

from dtcc_core.model import City
from dtcc_core.model import Surface, MultiSurface, Mesh, Building, BuildingPart
from dtcc_core.model.object.object import GeometryType
from dtcc_core.model.object.terrain import Terrain
from dtcc_core.model.geometry import Bounds

from dtcc_core.io.cityjson.write_cityjson import (
    to_cityjson,
    to_cityjson_surface,
    to_cityjson_multisurface, 
    to_cityjson_mesh,
    to_cityjson_terrain_mesh,
    save,
    _geometry_type_to_lod,
    _add_object_geometries
)


class TestGeometryConverters:
    """Test individual geometry conversion functions."""
    
    def test_to_cityjson_surface(self):
        """Test Surface to CityJSON conversion."""
        # Create a simple triangular surface
        vertices = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0], 
            [0.5, 1.0, 0.0]
        ])
        surface = Surface(vertices=vertices)
        
        vertices_list = []
        result = to_cityjson_surface(surface, vertices_list, 1000.0)
        
        # Check structure
        assert result["type"] == "MultiSurface"
        assert "boundaries" in result
        assert "semantics" in result
        assert len(result["boundaries"]) == 1
        assert len(vertices_list) == 3
        
        # Check vertex indexing
        boundary = result["boundaries"][0][0]  # First surface, exterior ring
        assert boundary == [0, 1, 2]
        
        # Check scaled vertices
        expected_vertices = [[0, 0, 0], [1000, 0, 0], [500, 1000, 0]]
        assert vertices_list == expected_vertices
    
    def test_to_cityjson_surface_with_holes(self):
        """Test Surface with holes conversion."""
        # Exterior ring
        vertices = np.array([
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 2.0, 0.0],
            [0.0, 2.0, 0.0]
        ])
        
        # Interior hole
        hole = np.array([
            [0.5, 0.5, 0.0],
            [1.5, 0.5, 0.0],
            [1.5, 1.5, 0.0],
            [0.5, 1.5, 0.0]
        ])
        
        surface = Surface(vertices=vertices, holes=[hole])
        
        vertices_list = []
        result = to_cityjson_surface(surface, vertices_list, 1000.0)
        
        # Check that we have exterior + interior rings
        assert len(result["boundaries"][0]) == 2  # exterior + hole
        assert len(vertices_list) == 8  # 4 exterior + 4 hole vertices
        
        # Check boundary indices
        exterior_boundary = result["boundaries"][0][0]
        hole_boundary = result["boundaries"][0][1]
        assert exterior_boundary == [0, 1, 2, 3]
        assert hole_boundary == [4, 5, 6, 7]
    
    def test_to_cityjson_multisurface(self):
        """Test MultiSurface to CityJSON conversion."""
        # Create two surfaces
        surface1 = Surface(vertices=np.array([
            [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]
        ]))
        surface2 = Surface(vertices=np.array([
            [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.5, 1.0, 0.0]
        ]))
        
        multisurface = MultiSurface(surfaces=[surface1, surface2])
        
        vertices_list = []
        result = to_cityjson_multisurface(multisurface, vertices_list, 1000.0)
        
        # Check structure
        assert result["type"] == "MultiSurface"
        assert len(result["boundaries"]) == 2
        assert len(vertices_list) == 6  # 3 + 3 vertices
        
        # Check semantics
        assert len(result["semantics"]["surfaces"]) == 2
        assert result["semantics"]["values"] == [0, 1]
    
    def test_to_cityjson_mesh(self):
        """Test Mesh to CityJSON conversion."""
        # Create a simple triangular mesh
        vertices = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 1.0, 0.0],
            [0.5, 0.5, 1.0]
        ])
        faces = np.array([
            [0, 1, 2],  # Bottom triangle
            [0, 1, 3],  # Side triangle 1
            [1, 2, 3],  # Side triangle 2
            [0, 2, 3]   # Side triangle 3
        ])
        mesh = Mesh(vertices=vertices, faces=faces)
        
        vertices_list = []
        result = to_cityjson_mesh(mesh, vertices_list, 1000.0)
        
        # Check structure
        assert result["type"] == "MultiSurface"
        assert len(result["boundaries"]) == 4  # 4 faces
        assert len(vertices_list) == 4  # 4 unique vertices
        
        # Check face boundaries
        expected_boundaries = [
            [[0, 1, 2]],  # Face 0
            [[0, 1, 3]],  # Face 1  
            [[1, 2, 3]],  # Face 2
            [[0, 2, 3]]   # Face 3
        ]
        assert result["boundaries"] == expected_boundaries

    def test_to_cityjson_terrain_mesh(self):
        """Test terrain Mesh to CityJSON CompositeSurface conversion."""
        # Create a simple triangular terrain mesh
        vertices = np.array([
            [0.0, 0.0, 0.0],
            [100.0, 0.0, 1.0],
            [50.0, 100.0, 0.5],
            [50.0, 50.0, 2.0]
        ])
        faces = np.array([
            [0, 1, 3],  # Triangle 1
            [1, 2, 3],  # Triangle 2
            [2, 0, 3]   # Triangle 3
        ])
        terrain_mesh = Mesh(vertices=vertices, faces=faces)
        
        vertices_list = []
        result = to_cityjson_terrain_mesh(terrain_mesh, vertices_list, 1000.0)
        
        # Check structure - should be CompositeSurface for terrain
        assert result["type"] == "CompositeSurface"
        assert len(result["boundaries"]) == 3  # 3 faces
        assert len(vertices_list) == 4  # 4 unique vertices
        
        # Check semantics - should be GroundSurface for terrain
        assert "semantics" in result
        semantics = result["semantics"]
        assert len(semantics["surfaces"]) == 3
        for surface in semantics["surfaces"]:
            assert surface["type"] == "GroundSurface"
        assert semantics["values"] == [0, 1, 2]
        
        # Check face boundaries
        expected_boundaries = [
            [[0, 1, 3]],  # Face 0
            [[1, 2, 3]],  # Face 1
            [[2, 0, 3]]   # Face 2
        ]
        assert result["boundaries"] == expected_boundaries


class TestObjectConversion:
    """Test conversion of DTCC objects to CityJSON."""
    
    def test_geometry_type_to_lod(self):
        """Test LOD mapping from GeometryType."""
        assert _geometry_type_to_lod(GeometryType.LOD0) == 0.0
        assert _geometry_type_to_lod(GeometryType.LOD1) == 1.0
        assert _geometry_type_to_lod(GeometryType.LOD2) == 2.0
        assert _geometry_type_to_lod(GeometryType.LOD3) == 3.0
        assert _geometry_type_to_lod(GeometryType.MESH) == 2.0
        assert _geometry_type_to_lod(GeometryType.SURFACE) == 2.0
    
    def test_add_object_geometries(self):
        """Test adding geometries from DTCC object."""
        # Create building with a mesh
        building = Building()
        mesh = Mesh(
            vertices=np.array([[0, 0, 0], [1, 0, 0], [0.5, 1, 0]]),
            faces=np.array([[0, 1, 2]])
        )
        building.add_geometry(mesh, GeometryType.MESH)
        
        obj_data = {"type": "Building", "geometry": [], "attributes": {}}
        vertices = []
        
        _add_object_geometries(building, obj_data, vertices, 1000.0)
        
        # Check that geometry was added
        assert len(obj_data["geometry"]) == 1
        assert obj_data["geometry"][0]["type"] == "MultiSurface"
        assert obj_data["geometry"][0]["lod"] == 2.0
        assert len(vertices) == 3


class TestCityConversion:
    """Test full City to CityJSON conversion."""
    
    def test_to_cityjson_empty_city(self):
        """Test conversion of empty city."""
        city = City()
        result = to_cityjson(city)
        
        # Check basic structure
        assert result["type"] == "CityJSON"
        assert result["version"] == "2.0"
        assert "transform" in result
        assert result["CityObjects"] == {}
        assert result["vertices"] == []
    
    def test_to_cityjson_city_with_bounds(self):
        """Test conversion of city with bounds."""
        city = City()
        city.bounds = Bounds(xmin=0.0, ymin=0.0, zmin=0.0, 
                           xmax=100.0, ymax=100.0, zmax=10.0)
        
        result = to_cityjson(city)
        
        # Check metadata
        assert "metadata" in result
        assert "geographicalExtent" in result["metadata"]
        assert result["metadata"]["geographicalExtent"] == [0.0, 0.0, 0.0, 100.0, 100.0, 10.0]
    
    def test_to_cityjson_city_with_building(self):
        """Test conversion of city with buildings."""
        city = City()
        
        # Create building with mesh geometry
        building = Building()
        building.id = "building_001"
        building.attributes = {"height": 10.0, "function": "residential"}
        
        mesh = Mesh(
            vertices=np.array([
                [0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [10.0, 10.0, 0.0], [0.0, 10.0, 0.0],  # Base
                [0.0, 0.0, 10.0], [10.0, 0.0, 10.0], [10.0, 10.0, 10.0], [0.0, 10.0, 10.0]  # Top
            ]),
            faces=np.array([
                [0, 1, 2], [0, 2, 3],  # Bottom
                [4, 7, 6], [4, 6, 5],  # Top
                [0, 4, 5], [0, 5, 1],  # Front
                [2, 6, 7], [2, 7, 3],  # Back
                [0, 3, 7], [0, 7, 4],  # Left
                [1, 5, 6], [1, 6, 2]   # Right
            ])
        )
        building.add_geometry(mesh, GeometryType.LOD2)
        city.add_building(building)
        
        result = to_cityjson(city)
        
        # Check building conversion
        assert "building_001" in result["CityObjects"]
        building_data = result["CityObjects"]["building_001"]
        assert building_data["type"] == "Building"
        assert building_data["attributes"]["height"] == 10.0
        assert building_data["attributes"]["function"] == "residential"
        assert len(building_data["geometry"]) == 1
        assert building_data["geometry"][0]["lod"] == 2.0
        
        # Check vertices were added
        assert len(result["vertices"]) == 8  # 8 unique vertices
    
    def test_to_cityjson_city_with_building_parts(self):
        """Test conversion of building with building parts."""
        city = City()
        
        # Create building with building part
        building = Building()
        building.id = "building_001"
        
        building_part = BuildingPart()
        building_part.id = "part_001"
        building_part.attributes = {"usage": "residential"}
        
        # Add geometry to building part
        surface = Surface(vertices=np.array([
            [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]
        ]))
        building_part.add_geometry(surface, GeometryType.LOD1)
        
        building.add_child(building_part)
        city.add_building(building)
        
        result = to_cityjson(city)
        
        # Check building and part structure
        assert "building_001" in result["CityObjects"]
        assert "part_001" in result["CityObjects"]
        
        building_data = result["CityObjects"]["building_001"]
        part_data = result["CityObjects"]["part_001"]
        
        assert "children" in building_data
        assert "part_001" in building_data["children"]
        assert part_data["parents"] == ["building_001"]
        assert len(part_data["geometry"]) == 1
    
    def test_to_cityjson_city_with_terrain(self):
        """Test conversion of city with terrain."""
        city = City()
        
        # Create terrain with mesh
        terrain = Terrain()
        terrain.id = "terrain_001"
        
        mesh = Mesh(
            vertices=np.array([
                [0.0, 0.0, 0.0], [100.0, 0.0, 2.0], 
                [50.0, 100.0, 1.0], [0.0, 100.0, 0.5]
            ]),
            faces=np.array([[0, 1, 2], [0, 2, 3]])
        )
        terrain.add_geometry(mesh, GeometryType.MESH)
        city.add_child(terrain)
        
        result = to_cityjson(city)
        
        # Check terrain conversion
        assert "terrain_001" in result["CityObjects"]
        terrain_data = result["CityObjects"]["terrain_001"]
        assert terrain_data["type"] == "TINRelief"
        assert len(terrain_data["geometry"]) == 1
        
        # Check terrain uses CompositeSurface geometry
        terrain_geometry = terrain_data["geometry"][0]
        assert terrain_geometry["type"] == "CompositeSurface"
        assert len(terrain_geometry["boundaries"]) == 2  # 2 faces
        
        # Check terrain uses GroundSurface semantics
        assert "semantics" in terrain_geometry
        semantics = terrain_geometry["semantics"]
        assert len(semantics["surfaces"]) == 2
        for surface in semantics["surfaces"]:
            assert surface["type"] == "GroundSurface"
    
    def test_to_cityjson_coordinate_scaling(self):
        """Test coordinate scaling and transformation."""
        city = City()
        building = Building()
        
        # Create geometry with decimal coordinates - must have faces to not be filtered out
        mesh = Mesh(
            vertices=np.array([
                [0.123, 0.456, 0.789],
                [1.234, 2.345, 3.456],
                [0.5, 1.0, 2.0]
            ]),
            faces=np.array([[0, 1, 2]])  # Add a face so it's not filtered as empty
        )
        building.add_geometry(mesh, GeometryType.MESH)
        city.add_building(building)
        
        # Test with scale factor 0.001 (default)
        result = to_cityjson(city, scale=0.001)
        
        # Vertices should be scaled by 1/0.001 = 1000 and rounded to integers
        expected_vertices = [
            [123, 456, 789],    # 0.123*1000, 0.456*1000, 0.789*1000
            [1234, 2345, 3456], # 1.234*1000, 2.345*1000, 3.456*1000
            [500, 1000, 2000]   # 0.5*1000, 1.0*1000, 2.0*1000
        ]
        assert result["vertices"] == expected_vertices
        
        # Check transform is correct
        assert result["transform"]["scale"] == [0.001, 0.001, 0.001]


class TestFileSaving:
    """Test saving CityJSON to files."""
    
    def test_save_to_json(self):
        """Test saving city to JSON file."""
        city = City()
        building = Building()
        building.attributes = {"test": "value"}
        city.add_building(building)
        
        with tempfile.NamedTemporaryFile(suffix='.json', delete=False) as tmp:
            tmp_path = Path(tmp.name)
        
        try:
            # Save the city
            save(city, tmp_path)
            
            # Verify file was created and contains valid JSON
            assert tmp_path.exists()
            
            with open(tmp_path, 'r') as f:
                loaded_data = json.load(f)
            
            assert loaded_data["type"] == "CityJSON"
            assert loaded_data["version"] == "2.0"
            
        finally:
            # Clean up
            if tmp_path.exists():
                tmp_path.unlink()
    
    def test_save_unsupported_format(self):
        """Test saving with unsupported file extension."""
        city = City()
        
        with tempfile.NamedTemporaryFile(suffix='.txt') as tmp:
            tmp_path = Path(tmp.name)
            
            with pytest.raises(ValueError, match="Unsupported file format"):
                save(city, tmp_path)


class TestEdgeCases:
    """Test edge cases and error conditions."""
    
    def test_empty_geometries(self):
        """Test handling of empty geometries."""
        city = City()
        building = Building()
        
        # Add empty mesh
        empty_mesh = Mesh()
        building.add_geometry(empty_mesh, GeometryType.MESH)
        city.add_building(building)
        
        # Should not crash
        result = to_cityjson(city)
        
        # Building should exist but have no geometry
        building_data = list(result["CityObjects"].values())[0]
        assert len(building_data["geometry"]) == 0
    
    def test_multiple_geometry_types(self):
        """Test object with multiple geometry types."""
        city = City()
        building = Building()
        
        # Add different geometry types
        mesh = Mesh(vertices=np.array([[0, 0, 0], [1, 0, 0], [0.5, 1, 0]]),
                   faces=np.array([[0, 1, 2]]))
        surface = Surface(vertices=np.array([[0, 0, 1], [1, 0, 1], [0.5, 1, 1]]))
        
        building.add_geometry(mesh, GeometryType.LOD1)
        building.add_geometry(surface, GeometryType.LOD2)
        city.add_building(building)
        
        result = to_cityjson(city)
        
        # Should have both geometries
        building_data = list(result["CityObjects"].values())[0]
        assert len(building_data["geometry"]) == 2
        
        lods = [geom["lod"] for geom in building_data["geometry"]]
        assert 1.0 in lods
        assert 2.0 in lods