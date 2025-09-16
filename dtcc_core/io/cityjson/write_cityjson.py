import json
from pathlib import Path
from typing import Dict, List

from ...model import City
from ...model import Surface, MultiSurface, Mesh, Building, BuildingPart
from ...model.object.object import GeometryType
from ...model.object.terrain import Terrain
import numpy as np


def to_cityjson_surface(surface: Surface, vertices: List, scale: float) -> Dict:
    """Convert a DTCC Surface to CityJSON geometry format."""
    vert_offset = len(vertices)
    # Scale and convert vertices to integers
    scaled_vertices = np.round(surface.vertices * scale).astype(int)
    vertices.extend(scaled_vertices.tolist())

    # Create boundary array for the surface (exterior ring)
    boundary = [list(range(vert_offset, len(vertices)))]

    # Add holes if they exist
    for hole in surface.holes:
        hole_offset = len(vertices)
        scaled_hole = np.round(hole * scale).astype(int)
        vertices.extend(scaled_hole.tolist())
        boundary.append(list(range(hole_offset, len(vertices))))

    return {
        "type": "MultiSurface",
        "boundaries": [boundary],
        "semantics": {"surfaces": [{"type": "WallSurface"}], "values": [[0]]},
    }


def to_cityjson_multisurface(
    multisurface: MultiSurface, vertices: List, scale: float
) -> Dict:
    """Convert a DTCC MultiSurface to CityJSON geometry format."""
    boundaries = []
    semantic_surfaces = []
    semantic_values = []

    for i, surface in enumerate(multisurface.surfaces):
        vert_offset = len(vertices)
        # Scale and convert vertices to integers
        scaled_vertices = np.round(surface.vertices * scale).astype(int)
        vertices.extend(scaled_vertices.tolist())

        # Create boundary for this surface (exterior ring)
        boundary = [list(range(vert_offset, len(vertices)))]

        # Add holes if they exist
        for hole in surface.holes:
            hole_offset = len(vertices)
            scaled_hole = np.round(hole * scale).astype(int)
            vertices.extend(scaled_hole.tolist())
            boundary.append(list(range(hole_offset, len(vertices))))

        boundaries.append(boundary)
        semantic_surfaces.append({"type": "WallSurface"})
        semantic_values.append(i)

    return {
        "type": "MultiSurface",
        "boundaries": boundaries,
        "semantics": {"surfaces": semantic_surfaces, "values": semantic_values},
    }


def to_cityjson_mesh(mesh: Mesh, vertices: List, scale: float) -> Dict:
    """Convert a DTCC Mesh to CityJSON geometry format."""
    vert_offset = len(vertices)
    # Scale and convert vertices to integers
    scaled_vertices = np.round(mesh.vertices * scale).astype(int)
    vertices.extend(scaled_vertices.tolist())

    # Convert faces to boundaries (each face becomes a surface)
    boundaries = []
    for face in mesh.faces:
        # Adjust indices to account for vertex offset
        face_indices = (face + vert_offset).tolist()
        # CityJSON expects boundaries as [[exterior_ring]]
        boundaries.append([face_indices])

    # Create semantic information for mesh faces
    semantic_surfaces = [{"type": "WallSurface"}] * len(mesh.faces)
    semantic_values = list(range(len(mesh.faces)))

    return {
        "type": "MultiSurface",
        "boundaries": boundaries,
        "semantics": {"surfaces": semantic_surfaces, "values": semantic_values},
    }


def to_cityjson_terrain_mesh(mesh: Mesh, vertices: List, scale: float) -> Dict:
    """Convert a DTCC terrain Mesh to CityJSON CompositeSurface format.
    
    According to CityJSON specification, TINRelief objects should use
    CompositeSurface geometry type with GroundSurface semantics.
    """
    vert_offset = len(vertices)
    # Scale and convert vertices to integers
    scaled_vertices = np.round(mesh.vertices * scale).astype(int)
    vertices.extend(scaled_vertices.tolist())

    # Convert faces to boundaries (each face becomes a surface in the composite)
    boundaries = []
    for face in mesh.faces:
        # Adjust indices to account for vertex offset
        face_indices = (face + vert_offset).tolist()
        # For CompositeSurface, each face is represented as [exterior_ring]
        boundaries.append([face_indices])

    # Create semantic information for terrain faces
    # All faces in terrain should be GroundSurface
    semantic_surfaces = [{"type": "GroundSurface"}] * len(mesh.faces)
    semantic_values = list(range(len(mesh.faces)))

    return {
        "type": "CompositeSurface",
        "boundaries": boundaries,
        "semantics": {"surfaces": semantic_surfaces, "values": semantic_values},
    }


def to_cityjson(city: City, scale: float = 0.001) -> Dict:
    """Convert a DTCC City to CityJSON format.

    Parameters
    ----------
    city : City
        The DTCC City object to convert
    scale : float, optional
        Scale factor for coordinate compression (default 0.001)

    Returns
    -------
    Dict
        CityJSON formatted dictionary
    """
    # Initialize CityJSON structure
    cityjson = {
        "type": "CityJSON",
        "version": "2.0",
        "transform": {"scale": [scale, scale, scale], "translate": [0.0, 0.0, 0.0]},
        "CityObjects": {},
        "vertices": [],
    }

    # Add metadata if city has bounds
    if city.bounds is not None:
        cityjson["metadata"] = {
            "referenceSystem": "https://www.opengis.net/def/crs/EPSG/0/7415",
            "geographicalExtent": [
                city.bounds.xmin,
                city.bounds.ymin,
                city.bounds.zmin,
                city.bounds.xmax,
                city.bounds.ymax,
                city.bounds.zmax,
            ],
        }

    # Convert scale to multiplication factor for vertex transformation
    scale_factor = 1.0 / scale
    vertices = []

    # Process buildings
    if Building in city.children:
        for building in city.children[Building]:
            building_data = {
                "type": "Building",
                "attributes": building.attributes.copy(),
                "geometry": [],
            }

            # Add building children if they exist
            if building.children:
                building_data["children"] = []

                # Process BuildingParts
                if BuildingPart in building.children:
                    for part in building.children[BuildingPart]:
                        part_data = {
                            "type": "BuildingPart",
                            "attributes": part.attributes.copy(),
                            "geometry": [],
                            "parents": [building.id],
                        }

                        # Add geometries from the part
                        _add_object_geometries(part, part_data, vertices, scale_factor)

                        cityjson["CityObjects"][part.id] = part_data
                        building_data["children"].append(part.id)

            # Add geometries directly attached to the building
            _add_object_geometries(building, building_data, vertices, scale_factor)

            cityjson["CityObjects"][building.id] = building_data

    # Process terrain
    if Terrain in city.children:
        for terrain in city.children[Terrain]:
            terrain_data = {
                "type": "TINRelief",
                "attributes": terrain.attributes.copy(),
                "geometry": [],
            }

            # Add terrain geometries
            _add_object_geometries(terrain, terrain_data, vertices, scale_factor)

            cityjson["CityObjects"][terrain.id] = terrain_data

    # Add all processed vertices to the final structure
    cityjson["vertices"] = vertices

    return cityjson


def _add_object_geometries(obj, obj_data: Dict, vertices: List, scale_factor: float):
    """Add geometries from a DTCC object to CityJSON object data."""
    # Check if this is a terrain object to use appropriate converter
    is_terrain = obj_data.get("type") == "TINRelief"
    
    for geom_type, geometry in obj.geometry.items():
        if geometry is None:
            continue

        # Skip empty geometries
        if isinstance(geometry, Mesh) and (
            len(geometry.vertices) == 0 or len(geometry.faces) == 0
        ):
            continue
        elif isinstance(geometry, Surface) and len(geometry.vertices) == 0:
            continue
        elif isinstance(geometry, MultiSurface) and len(geometry.surfaces) == 0:
            continue

        lod = _geometry_type_to_lod(geom_type)

        if isinstance(geometry, Mesh):
            # Use terrain-specific converter for terrain objects
            if is_terrain:
                geom_data = to_cityjson_terrain_mesh(geometry, vertices, scale_factor)
            else:
                geom_data = to_cityjson_mesh(geometry, vertices, scale_factor)
            geom_data["lod"] = lod
            obj_data["geometry"].append(geom_data)

        elif isinstance(geometry, Surface):
            geom_data = to_cityjson_surface(geometry, vertices, scale_factor)
            geom_data["lod"] = lod
            obj_data["geometry"].append(geom_data)

        elif isinstance(geometry, MultiSurface):
            geom_data = to_cityjson_multisurface(geometry, vertices, scale_factor)
            geom_data["lod"] = lod
            obj_data["geometry"].append(geom_data)


def _geometry_type_to_lod(geom_type: GeometryType) -> float:
    """Convert DTCC GeometryType to CityJSON LOD value."""
    lod_mapping = {
        GeometryType.LOD0: 0.0,
        GeometryType.LOD1: 1.0,
        GeometryType.LOD2: 2.0,
        GeometryType.LOD3: 3.0,
        GeometryType.MESH: 2.0,  # Default mesh to LOD2
        GeometryType.SURFACE: 2.0,  # Default surface to LOD2
        GeometryType.MULTISURFACE: 2.0,  # Default multisurface to LOD2
    }
    return lod_mapping.get(geom_type, 2.0)


def save(city: City, path: Path, scale: float = 0.001):
    """Save a city to a CityJSON file.

    Parameters
    ----------
    city : City
        The city to save.
    path : str or Path
        Path to the file.
    scale : float, optional
        Scale factor for coordinate compression (default 0.001)

    Raises
    ------
    ValueError
        If the file format is not supported
    """
    cj = to_cityjson(city, scale=scale)
    path = Path(path)

    if path.suffix.lower() == ".json":
        with open(path, "w") as file:
            json.dump(cj, file, indent=2)
    else:
        raise ValueError(
            f"Unsupported file format: {path.suffix}. Only .json is supported."
        )
