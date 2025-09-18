"""
Simplified geometry converters for CityJSON format.

A simpler alternative to the ABC-based design - just functions and utilities.
"""

import numpy as np
from typing import Dict, List
from dataclasses import dataclass

from ...model import Surface, MultiSurface, Mesh
from ...model.object.object import GeometryType


@dataclass
class CityJSONConfig:
    """Configuration for CityJSON conversion."""

    semantic_types: Dict[str, str] = None
    lod_mapping: Dict[GeometryType, float] = None

    def __post_init__(self):
        if self.semantic_types is None:
            self.semantic_types = {
                "surface": "WallSurface",
                "terrain": "GroundSurface",
                "building": "WallSurface",
            }

        default_lod_mapping = {
            GeometryType.LOD0: 0.0,
            GeometryType.LOD1: 1.0,
            GeometryType.LOD2: 2.0,
            GeometryType.LOD3: 3.0,
            GeometryType.MESH: 2.0,
            GeometryType.SURFACE: 2.0,
            GeometryType.MULTISURFACE: 2.0,
        }

        if self.lod_mapping is None:
            self.lod_mapping = default_lod_mapping
        else:
            self.lod_mapping = {**default_lod_mapping, **self.lod_mapping}


# Utility functions (shared by all converters)
def scale_vertices(vertices: np.ndarray, scale: float, vertices_list: List) -> int:
    """Scale vertices and add to global vertices list. Returns offset."""
    if not isinstance(vertices, np.ndarray) or vertices.size == 0:
        raise ValueError("Invalid vertices: must be non-empty numpy array")

    vert_offset = len(vertices_list)
    try:
        scaled_vertices = np.round(vertices * scale).astype(int)
        # Convert NumPy int64 to Python int for JSON serialization
        python_int_vertices = [[int(x), int(y), int(z)] for x, y, z in scaled_vertices.tolist()]
        vertices_list.extend(python_int_vertices)
    except (ValueError, TypeError) as e:
        raise ValueError(f"Error scaling vertices: {e}")

    return vert_offset


def create_boundary(indices: np.ndarray, offset: int) -> List:
    """Create boundary from indices with offset applied."""
    if not isinstance(indices, np.ndarray) or indices.size == 0:
        raise ValueError("Invalid indices: must be non-empty numpy array")

    try:
        # Convert NumPy int64 to Python int for JSON serialization
        result = (indices + offset).tolist()
        return [int(x) for x in result]
    except (ValueError, TypeError) as e:
        raise ValueError(f"Error creating boundary: {e}")


# Converter functions (one per geometry type)
def convert_surface(
    surface: Surface, vertices: List, scale: float, config: CityJSONConfig = None
) -> Dict:
    """Convert Surface to CityJSON format."""
    config = config or CityJSONConfig()
    vert_offset = scale_vertices(surface.vertices, scale, vertices)

    # Create boundary array for the surface (exterior ring)
    boundary = [create_boundary(np.arange(len(surface.vertices)), vert_offset)]

    # Add holes if they exist
    for hole in surface.holes:
        hole_offset = scale_vertices(hole, scale, vertices)
        boundary.append(create_boundary(np.arange(len(hole)), hole_offset))

    return {
        "type": "MultiSurface",
        "boundaries": [boundary],
        "semantics": {
            "surfaces": [{"type": config.semantic_types["surface"]}],
            "values": [[0]],
        },
    }


def convert_multisurface(
    multisurface: MultiSurface,
    vertices: List,
    scale: float,
    config: CityJSONConfig = None,
) -> Dict:
    """Convert MultiSurface to CityJSON format."""
    config = config or CityJSONConfig()
    boundaries = []
    semantic_surfaces = []
    semantic_values = []

    for i, surface in enumerate(multisurface.surfaces):
        vert_offset = scale_vertices(surface.vertices, scale, vertices)

        # Create boundary for this surface (exterior ring)
        boundary = [create_boundary(np.arange(len(surface.vertices)), vert_offset)]

        # Add holes if they exist
        for hole in surface.holes:
            hole_offset = scale_vertices(hole, scale, vertices)
            boundary.append(create_boundary(np.arange(len(hole)), hole_offset))

        boundaries.append(boundary)
        semantic_surfaces.append({"type": config.semantic_types["surface"]})
        semantic_values.append(i)

    return {
        "type": "MultiSurface",
        "boundaries": boundaries,
        "semantics": {"surfaces": semantic_surfaces, "values": semantic_values},
    }


def convert_mesh(
    mesh: Mesh, vertices: List, scale: float, config: CityJSONConfig = None
) -> Dict:
    """Convert Mesh to CityJSON format."""
    config = config or CityJSONConfig()
    vert_offset = scale_vertices(mesh.vertices, scale, vertices)

    # Convert faces to boundaries (each face becomes a surface)
    boundaries = []
    for face in mesh.faces:
        # Adjust indices to account for vertex offset
        face_indices = create_boundary(face, vert_offset)
        # CityJSON expects boundaries as [[exterior_ring]]
        boundaries.append([face_indices])

    # Create semantic information for mesh faces
    semantic_surfaces = [{"type": config.semantic_types["surface"]}] * len(mesh.faces)
    semantic_values = list(range(len(mesh.faces)))

    return {
        "type": "MultiSurface",
        "boundaries": boundaries,
        "semantics": {"surfaces": semantic_surfaces, "values": semantic_values},
    }


def convert_terrain_mesh(
    mesh: Mesh, vertices: List, scale: float, config: CityJSONConfig = None
) -> Dict:
    """Convert terrain Mesh to CityJSON format."""
    config = config or CityJSONConfig()
    vert_offset = scale_vertices(mesh.vertices, scale, vertices)

    # Convert faces to boundaries (each face becomes a surface in the composite)
    boundaries = []
    for face in mesh.faces:
        # Adjust indices to account for vertex offset
        face_indices = create_boundary(face, vert_offset)
        # For CompositeSurface, each face is represented as [exterior_ring]
        boundaries.append([face_indices])

    # Create semantic information for terrain faces
    # All faces in terrain should be GroundSurface
    semantic_surfaces = [{"type": config.semantic_types["terrain"]}] * len(mesh.faces)
    semantic_values = list(range(len(mesh.faces)))

    return {
        "type": "CompositeSurface",
        "boundaries": boundaries,
        "semantics": {"surfaces": semantic_surfaces, "values": semantic_values},
    }


# Factory/registry pattern (simpler than inheritance)
_geometry_converters = {
    Surface: convert_surface,
    MultiSurface: convert_multisurface,
    Mesh: convert_mesh,
}


def get_converter(geometry_type: type, config: CityJSONConfig = None):
    """Get appropriate converter function for geometry type."""
    converter_func = _geometry_converters.get(geometry_type)
    if converter_func is None:
        raise ValueError(f"No converter available for geometry type: {geometry_type}")

    # Return a partial function with config baked in
    def converter(geometry, vertices, scale):
        return converter_func(geometry, vertices, scale, config)

    return converter


def get_terrain_converter(config: CityJSONConfig = None):
    """Get terrain-specific mesh converter."""

    def converter(geometry, vertices, scale):
        return convert_terrain_mesh(geometry, vertices, scale, config)

    return converter


def geometry_type_to_lod(
    geom_type: GeometryType, config: CityJSONConfig = None
) -> float:
    """Convert DTCC GeometryType to CityJSON LOD value."""
    config = config or CityJSONConfig()
    return config.lod_mapping.get(geom_type, 2.0)
