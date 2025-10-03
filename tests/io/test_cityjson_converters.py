"""
Unit tests for CityJSON geometry converters.

Tests the refactored converter architecture including:
- BaseGeometryConverter functionality
- Specialized converters (Surface, MultiSurface, Mesh, TerrainMesh)
- GeometryConverterFactory
- CityJSONConfig
- Error handling and edge cases
"""

import pytest
import numpy as np
from unittest.mock import Mock

from dtcc_core.model import Surface, MultiSurface, Mesh
from dtcc_core.model.object.object import GeometryType

from dtcc_core.io.cityjson.converters import (
    convert_surface,
    convert_multisurface,
    convert_mesh,
    convert_terrain_mesh,
    get_converter,
    get_terrain_converter,
    CityJSONConfig,
    geometry_type_to_lod,
    scale_vertices,
    create_boundary,
)


# CityJSONConfig tests
def test_default_config():
    """Test default configuration values."""
    config = CityJSONConfig()

    assert config.semantic_types['surface'] == 'WallSurface'
    assert config.semantic_types['terrain'] == 'GroundSurface'
    assert config.semantic_types['building'] == 'WallSurface'

    assert config.lod_mapping[GeometryType.LOD0] == 0.0
    assert config.lod_mapping[GeometryType.LOD1] == 1.0
    assert config.lod_mapping[GeometryType.LOD2] == 2.0
    assert config.lod_mapping[GeometryType.LOD3] == 3.0
    assert config.lod_mapping[GeometryType.MESH] == 2.0


def test_custom_config():
    """Test custom configuration values."""
    custom_semantics = {
        'surface': 'CustomWall',
        'terrain': 'CustomGround',
        'building': 'CustomBuilding'
    }
    custom_lod_mapping = {GeometryType.LOD0: 0.5, GeometryType.LOD1: 1.5}

    config = CityJSONConfig(
        semantic_types=custom_semantics,
        lod_mapping=custom_lod_mapping
    )

    assert config.semantic_types['surface'] == 'CustomWall'
    assert config.semantic_types['terrain'] == 'CustomGround'
    assert config.lod_mapping[GeometryType.LOD0] == 0.5
    assert config.lod_mapping[GeometryType.LOD1] == 1.5
    # Non-customized values should remain default
    assert config.lod_mapping[GeometryType.LOD2] == 2.0


# Utility function tests
def test_scale_vertices():
    """Test vertex scaling and addition to global list."""
    test_vertices = np.array([
        [1.5, 2.3, 3.7],
        [4.1, 5.9, 6.2],
        [7.8, 9.4, 10.1]
    ])
    global_vertices = []

    offset = scale_vertices(test_vertices, scale=2.0, vertices_list=global_vertices)

    assert offset == 0  # First addition, so offset should be 0
    assert len(global_vertices) == 3
    # Check scaled and rounded values: 1.5*2=3, 2.3*2=4.6→5, 3.7*2=7.4→7, etc.
    expected = [[3, 5, 7], [8, 12, 12], [16, 19, 20]]
    assert global_vertices == expected


def test_scale_vertices_empty_array():
    """Test vertex scaling with empty array."""
    with pytest.raises(ValueError, match="Invalid vertices"):
        scale_vertices(np.array([]), 1.0, [])


def test_scale_vertices_invalid_input():
    """Test vertex scaling with invalid input."""
    with pytest.raises(ValueError, match="Invalid vertices"):
        scale_vertices(None, 1.0, [])


def test_create_boundary():
    """Test boundary creation with offset."""
    indices = np.array([0, 1, 2, 3])
    boundary = create_boundary(indices, offset=5)

    assert boundary == [5, 6, 7, 8]


def test_create_boundary_empty_array():
    """Test boundary creation with empty array."""
    with pytest.raises(ValueError, match="Invalid indices"):
        create_boundary(np.array([]), 0)


# Surface converter tests
def test_convert_surface_simple():
    """Test conversion of simple surface without holes."""
    config = CityJSONConfig()

    # Create mock surface
    mock_surface = Mock()
    mock_surface.vertices = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0]
    ])
    mock_surface.holes = []

    vertices = []
    result = convert_surface(mock_surface, vertices, scale=1000.0, config=config)

    assert result['type'] == 'MultiSurface'
    assert 'boundaries' in result
    assert 'semantics' in result
    assert len(result['boundaries']) == 1
    assert len(result['boundaries'][0]) == 1  # One exterior ring
    assert result['boundaries'][0][0] == [0, 1, 2, 3]  # Vertex indices
    assert len(vertices) == 4

    # Check semantics
    assert result['semantics']['surfaces'][0]['type'] == 'WallSurface'
    assert result['semantics']['values'] == [[0]]


def test_convert_surface_with_holes():
    """Test conversion of surface with holes."""
    config = CityJSONConfig()

    # Create mock surface with hole
    mock_surface = Mock()
    mock_surface.vertices = np.array([
        [0.0, 0.0, 0.0],   # Exterior
        [2.0, 0.0, 0.0],
        [2.0, 2.0, 0.0],
        [0.0, 2.0, 0.0]
    ])
    mock_surface.holes = [
        np.array([
            [0.5, 0.5, 0.0],   # Hole
            [1.5, 0.5, 0.0],
            [1.5, 1.5, 0.0],
            [0.5, 1.5, 0.0]
        ])
    ]

    vertices = []
    result = convert_surface(mock_surface, vertices, scale=1000.0, config=config)

    assert len(result['boundaries'][0]) == 2  # Exterior + hole
    assert result['boundaries'][0][0] == [0, 1, 2, 3]  # Exterior boundary
    assert result['boundaries'][0][1] == [4, 5, 6, 7]  # Hole boundary
    assert len(vertices) == 8  # 4 exterior + 4 hole vertices


# MultiSurface converter tests
def test_convert_multisurface():
    """Test conversion of MultiSurface."""
    config = CityJSONConfig()

    # Create mock MultiSurface
    mock_surface1 = Mock()
    mock_surface1.vertices = np.array([[0, 0, 0], [1, 0, 0], [0.5, 1, 0]])
    mock_surface1.holes = []

    mock_surface2 = Mock()
    mock_surface2.vertices = np.array([[1, 0, 0], [2, 0, 0], [1.5, 1, 0]])
    mock_surface2.holes = []

    mock_multisurface = Mock()
    mock_multisurface.surfaces = [mock_surface1, mock_surface2]

    vertices = []
    result = convert_multisurface(mock_multisurface, vertices, scale=1000.0, config=config)

    assert result['type'] == 'MultiSurface'
    assert len(result['boundaries']) == 2
    assert len(vertices) == 6  # 3 + 3 vertices

    # Check semantics
    assert len(result['semantics']['surfaces']) == 2
    assert result['semantics']['values'] == [0, 1]


# Mesh converter tests
def test_convert_mesh():
    """Test conversion of Mesh."""
    config = CityJSONConfig()

    # Create mock mesh
    mock_mesh = Mock()
    mock_mesh.vertices = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.5, 1.0, 0.0],
        [0.5, 0.5, 1.0]
    ])
    mock_mesh.faces = np.array([
        [0, 1, 2],  # Face 1
        [0, 1, 3],  # Face 2
    ])

    vertices = []
    result = convert_mesh(mock_mesh, vertices, scale=1000.0, config=config)

    assert result['type'] == 'MultiSurface'
    assert len(result['boundaries']) == 2  # 2 faces
    assert len(vertices) == 4  # 4 vertices

    # Check face boundaries
    assert result['boundaries'][0] == [[0, 1, 2]]
    assert result['boundaries'][1] == [[0, 1, 3]]

    # Check semantics
    assert len(result['semantics']['surfaces']) == 2
    assert result['semantics']['values'] == [0, 1]


# Terrain mesh converter tests
def test_convert_terrain_mesh():
    """Test conversion of terrain Mesh."""
    config = CityJSONConfig()

    # Create mock terrain mesh
    mock_mesh = Mock()
    mock_mesh.vertices = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.5, 1.0, 0.0]
    ])
    mock_mesh.faces = np.array([[0, 1, 2]])

    vertices = []
    result = convert_terrain_mesh(mock_mesh, vertices, scale=1000.0, config=config)

    assert result['type'] == 'CompositeSurface'
    assert len(result['boundaries']) == 1
    assert len(vertices) == 3

    # Check terrain semantics
    assert result['semantics']['surfaces'][0]['type'] == 'GroundSurface'


# Factory/registry tests
def test_get_converter():
    """Test that factory creates appropriate converter functions."""
    config = CityJSONConfig()

    # Test Surface converter creation
    surface_converter = get_converter(Surface, config)
    assert callable(surface_converter)

    # Test MultiSurface converter creation
    multisurface_converter = get_converter(MultiSurface, config)
    assert callable(multisurface_converter)

    # Test Mesh converter creation
    mesh_converter = get_converter(Mesh, config)
    assert callable(mesh_converter)


def test_get_terrain_converter():
    """Test terrain converter creation."""
    config = CityJSONConfig()
    terrain_converter = get_terrain_converter(config)
    assert callable(terrain_converter)


def test_get_converter_error_handling():
    """Test factory error handling for unknown types."""
    config = CityJSONConfig()

    with pytest.raises(ValueError, match="No converter available"):
        get_converter(str, config)


def test_get_converter_uses_config():
    """Test that factory passes config to converter functions."""
    custom_config = CityJSONConfig(semantic_types={'surface': 'CustomType'})

    converter = get_converter(Surface, custom_config)

    # Create mock surface to test
    mock_surface = Mock()
    mock_surface.vertices = np.array([[0, 0, 0], [1, 0, 0], [0.5, 1, 0]])
    mock_surface.holes = []

    vertices = []
    result = converter(mock_surface, vertices, scale=1.0)

    # Check that custom semantic type is used
    assert result['semantics']['surfaces'][0]['type'] == 'CustomType'


# GeometryTypeToLod tests
def test_standard_geometry_types():
    """Test LOD mapping for standard geometry types."""
    config = CityJSONConfig()

    assert geometry_type_to_lod(GeometryType.LOD0, config) == 0.0
    assert geometry_type_to_lod(GeometryType.LOD1, config) == 1.0
    assert geometry_type_to_lod(GeometryType.LOD2, config) == 2.0
    assert geometry_type_to_lod(GeometryType.LOD3, config) == 3.0
    assert geometry_type_to_lod(GeometryType.MESH, config) == 2.0
    assert geometry_type_to_lod(GeometryType.SURFACE, config) == 2.0


def test_unknown_geometry_type():
    """Test handling of unknown geometry types."""
    config = CityJSONConfig()
    assert geometry_type_to_lod(None, config) == 2.0  # Should return default


def test_custom_lod_mapping():
    """Test custom LOD mapping."""
    custom_config = CityJSONConfig(lod_mapping={GeometryType.LOD0: 0.5})
    assert geometry_type_to_lod(GeometryType.LOD0, custom_config) == 0.5


# Integration tests
def test_full_conversion_workflow():
    """Test a complete conversion workflow."""
    # Create a surface
    surface = Surface(vertices=np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0]
    ]))

    # Get converter and convert
    config = CityJSONConfig()
    converter = get_converter(Surface, config)
    vertices = []
    result = converter(surface, vertices, scale=1000.0)

    # Verify the result
    assert result['type'] == 'MultiSurface'
    assert len(result['boundaries']) == 1
    assert len(vertices) == 4
    assert all(isinstance(v, int) for v in vertices[0])  # Vertices should be integers


def test_config_persistence():
    """Test that configuration is properly maintained across conversions."""
    custom_config = CityJSONConfig(semantic_types={'surface': 'CustomWall'})

    converter = get_converter(Surface, custom_config)

    # Create mock surface
    mock_surface = Mock()
    mock_surface.vertices = np.array([[0, 0, 0], [1, 0, 0], [0.5, 1, 0]])
    mock_surface.holes = []

    vertices = []
    result = converter(mock_surface, vertices, scale=1.0)

    # Check that custom semantic type is used
    assert result['semantics']['surfaces'][0]['type'] == 'CustomWall'