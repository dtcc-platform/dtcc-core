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
    # Rounding mode for quantization: "round" (default), "floor", or "ceil"
    rounding_mode: str = "round"

    def __post_init__(self):
        """
        Populate default mappings and validate rounding behavior.

        Raises
        ------
        ValueError
            If ``rounding_mode`` is not one of ``\"round\"``, ``\"floor\"``, or ``\"ceil\"``.
        """
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

        # Validate rounding_mode
        if self.rounding_mode not in ("round", "floor", "ceil"):
            raise ValueError(f"Invalid rounding_mode: {self.rounding_mode}")


# Utility functions (shared by all converters)
def scale_vertices(vertices: np.ndarray, scale: float, vertices_list: List, rounding_mode: str = "round") -> int:
    """Scale vertices and add to global vertices list. Returns offset.

    Parameters
    ----------
    vertices : np.ndarray
        Float coordinates (N x 3)
    scale : float
        Quantization multiplier (usually 1/transform_scale)
    vertices_list : List
        Global list of quantized integer vertices (mutated)
    rounding_mode : str
        "round" (default), "floor", or "ceil"
    """
    if not isinstance(vertices, np.ndarray) or vertices.size == 0:
        raise ValueError("Invalid vertices: must be non-empty numpy array")

    vert_offset = len(vertices_list)
    try:
        q = vertices * scale
        if rounding_mode == "round":
            q = np.round(q)
        elif rounding_mode == "floor":
            q = np.floor(q)
        elif rounding_mode == "ceil":
            q = np.ceil(q)
        else:
            raise ValueError(f"Unknown rounding_mode: {rounding_mode}")

        scaled_vertices = q.astype(int)
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


class VertexIndexer:
    """
    Global vertex indexer that quantizes coordinates and deduplicates vertices.

    Maintains:
    - vertices: List[List[int]] of unique quantized vertices
    - _index: Dict[(int,int,int) -> int] mapping vertex to global index
    """

    def __init__(self):
        """Initialize empty vertex storage and lookup index."""
        self.vertices: List[List[int]] = []
        self._index = {}

    def _quantize(self, points: np.ndarray, factor: float, rounding_mode: str = "round") -> np.ndarray:
        if not isinstance(points, np.ndarray) or points.size == 0:
            raise ValueError("Invalid vertices: must be non-empty numpy array")
        q = points * factor
        if rounding_mode == "round":
            q = np.round(q)
        elif rounding_mode == "floor":
            q = np.floor(q)
        elif rounding_mode == "ceil":
            q = np.ceil(q)
        else:
            raise ValueError(f"Unknown rounding_mode: {rounding_mode}")
        return q.astype(int)

    def add_points(self, points: np.ndarray, factor: float, rounding_mode: str = "round") -> List[int]:
        """
        Quantize and add points to the global vertex list.

        Parameters
        ----------
        points : np.ndarray
            Array of shape (N, 3) with floating point coordinates.
        factor : float
            Quantization multiplier applied before rounding.
        rounding_mode : str, default "round"
            Rounding strategy: ``\"round\"``, ``\"floor\"``, or ``\"ceil\"``.

        Returns
        -------
        list[int]
            Global indices corresponding to the (de-duplicated) quantized vertices.

        Raises
        ------
        ValueError
            If ``rounding_mode`` is invalid.
        """
        q = self._quantize(points, factor, rounding_mode)
        indices: List[int] = []
        for x, y, z in q.tolist():
            key = (int(x), int(y), int(z))
            idx = self._index.get(key)
            if idx is None:
                idx = len(self.vertices)
                self.vertices.append([int(x), int(y), int(z)])
                self._index[key] = idx
            indices.append(idx)
        return indices


# Converter functions (one per geometry type)
def convert_surface(
    surface: Surface,
    vertices: List,
    scale: float,
    config: CityJSONConfig = None,
    indexer=None,
) -> Dict:
    """Convert Surface to CityJSON format.

    Parameters
    ----------
    surface : Surface
    vertices : list
        Global vertices list (mutated). When using indexer, this should reference
        indexer.vertices for consistency.
    scale : float
        Quantization multiplier (1/transform_scale)
    config : CityJSONConfig
    indexer : VertexIndexer, optional
        If provided, global deduplication and quantization are handled by the indexer.
    """
    config = config or CityJSONConfig()

    if indexer is not None:
        # Exterior ring
        outer_idx = indexer.add_points(surface.vertices, scale, config.rounding_mode)
        boundary = [outer_idx]
        # Holes
        for hole in surface.holes:
            hole_idx = indexer.add_points(hole, scale, config.rounding_mode)
            boundary.append(hole_idx)
    else:
        vert_offset = scale_vertices(surface.vertices, scale, vertices, config.rounding_mode)
        boundary = [create_boundary(np.arange(len(surface.vertices)), vert_offset)]
        for hole in surface.holes:
            hole_offset = scale_vertices(hole, scale, vertices, config.rounding_mode)
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
    indexer=None,
) -> Dict:
    """Convert MultiSurface to CityJSON format."""
    config = config or CityJSONConfig()
    boundaries = []
    semantic_surfaces = []
    semantic_values = []

    for i, surface in enumerate(multisurface.surfaces):
        if indexer is not None:
            # Exterior ring
            outer_idx = indexer.add_points(surface.vertices, scale, config.rounding_mode)
            boundary = [outer_idx]
            # Holes
            for hole in surface.holes:
                hole_idx = indexer.add_points(hole, scale, config.rounding_mode)
                boundary.append(hole_idx)
        else:
            vert_offset = scale_vertices(surface.vertices, scale, vertices, config.rounding_mode)
            boundary = [create_boundary(np.arange(len(surface.vertices)), vert_offset)]
            for hole in surface.holes:
                hole_offset = scale_vertices(hole, scale, vertices, config.rounding_mode)
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
    mesh: Mesh, vertices: List, scale: float, config: CityJSONConfig = None, indexer=None
) -> Dict:
    """Convert Mesh to CityJSON format."""
    config = config or CityJSONConfig()

    boundaries = []
    if indexer is not None:
        # Add all mesh vertices and get their global indices
        global_indices = indexer.add_points(mesh.vertices, scale, config.rounding_mode)
        # Convert faces to boundaries using mapped indices
        for face in mesh.faces:
            face_list = [int(global_indices[int(i)]) for i in np.asarray(face).tolist()]
            boundaries.append([face_list])
    else:
        vert_offset = scale_vertices(mesh.vertices, scale, vertices, config.rounding_mode)
        for face in mesh.faces:
            face_indices = create_boundary(face, vert_offset)
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
    mesh: Mesh, vertices: List, scale: float, config: CityJSONConfig = None, indexer=None
) -> Dict:
    """Convert terrain Mesh to CityJSON format."""
    config = config or CityJSONConfig()

    boundaries = []
    if indexer is not None:
        global_indices = indexer.add_points(mesh.vertices, scale, config.rounding_mode)
        for face in mesh.faces:
            face_list = [int(global_indices[int(i)]) for i in np.asarray(face).tolist()]
            boundaries.append([face_list])
    else:
        vert_offset = scale_vertices(mesh.vertices, scale, vertices, config.rounding_mode)
        for face in mesh.faces:
            face_indices = create_boundary(face, vert_offset)
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
        """
        Convert a geometry to CityJSON using the registered converter.

        Parameters
        ----------
        geometry : Surface or MultiSurface or Mesh
            Geometry instance to convert.
        vertices : list
            Shared vertex accumulator mutated during conversion.
        scale : float
            Quantization multiplier (1/transform_scale).

        Returns
        -------
        dict
            CityJSON geometry representation.
        """
        return converter_func(geometry, vertices, scale, config)

    return converter


def get_terrain_converter(config: CityJSONConfig = None):
    """Get terrain-specific mesh converter."""

    def converter(geometry, vertices, scale):
        """
        Convert a terrain mesh to CityJSON.

        Parameters
        ----------
        geometry : Mesh
            Terrain mesh to convert.
        vertices : list
            Shared vertex accumulator mutated during conversion.
        scale : float
            Quantization multiplier (1/transform_scale).

        Returns
        -------
        dict
            CityJSON geometry representation tailored for terrain meshes.
        """
        return convert_terrain_mesh(geometry, vertices, scale, config)

    return converter


def geometry_type_to_lod(
    geom_type: GeometryType, config: CityJSONConfig = None
) -> float:
    """Convert DTCC GeometryType to CityJSON LOD value."""
    config = config or CityJSONConfig()
    return config.lod_mapping.get(geom_type, 2.0)
