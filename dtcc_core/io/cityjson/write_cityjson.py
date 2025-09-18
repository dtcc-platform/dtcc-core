import json
from pathlib import Path
from typing import Dict, List

from dtcc_core.model import (
    City,
    Surface,
    MultiSurface,
    Mesh,
    Building,
    BuildingPart,
    Terrain,
    GeometryType,
)

import numpy as np

from .converters import (
    get_converter,
    get_terrain_converter,
    geometry_type_to_lod,
    CityJSONConfig,
    convert_surface,
    convert_multisurface,
    convert_mesh,
    convert_terrain_mesh,
)


def to_cityjson_surface(
    surface: Surface, vertices: List, scale: float, config: CityJSONConfig = None
) -> Dict:
    """Convert a DTCC Surface to CityJSON geometry format."""
    return convert_surface(surface, vertices, scale, config)


def to_cityjson_multisurface(
    multisurface: MultiSurface,
    vertices: List,
    scale: float,
    config: CityJSONConfig = None,
) -> Dict:
    """Convert a DTCC MultiSurface to CityJSON geometry format."""
    return convert_multisurface(multisurface, vertices, scale, config)


def to_cityjson_mesh(
    mesh: Mesh, vertices: List, scale: float, config: CityJSONConfig = None
) -> Dict:
    """Convert a DTCC Mesh to CityJSON geometry format."""
    return convert_mesh(mesh, vertices, scale, config)


def to_cityjson_terrain_mesh(
    mesh: Mesh, vertices: List, scale: float, config: CityJSONConfig = None
) -> Dict:
    """Convert a DTCC terrain Mesh to CityJSON CompositeSurface format.

    According to CityJSON specification, TINRelief objects should use
    CompositeSurface geometry type with GroundSurface semantics.
    """
    return convert_terrain_mesh(mesh, vertices, scale, config)


def to_cityjson(
    city: City, scale: float = 0.001, config: CityJSONConfig = None
) -> Dict:
    """Convert a DTCC City to CityJSON format.

    Parameters
    ----------
    city : City
        The DTCC City object to convert
    scale : float, optional
        Scale factor for coordinate compression (default 0.001)
    config : CityJSONConfig, optional
        Configuration for conversion settings

    Returns
    -------
    Dict
        CityJSON formatted dictionary
    """
    config = config or CityJSONConfig()

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
                float(city.bounds.xmin),
                float(city.bounds.ymin),
                float(city.bounds.zmin),
                float(city.bounds.xmax),
                float(city.bounds.ymax),
                float(city.bounds.zmax),
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
                        _add_object_geometries(
                            part, part_data, vertices, scale_factor, config
                        )

                        cityjson["CityObjects"][part.id] = part_data
                        building_data["children"].append(part.id)

            # Add geometries directly attached to the building
            _add_object_geometries(
                building, building_data, vertices, scale_factor, config
            )

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
            _add_object_geometries(
                terrain, terrain_data, vertices, scale_factor, config
            )

            cityjson["CityObjects"][terrain.id] = terrain_data

    # Add all processed vertices to the final structure
    cityjson["vertices"] = vertices

    return cityjson


def _add_object_geometries(
    obj,
    obj_data: Dict,
    vertices: List,
    scale_factor: float,
    config: CityJSONConfig = None,
):
    """Add geometries from a DTCC object to CityJSON object data."""
    config = config or CityJSONConfig()

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

        lod = geometry_type_to_lod(geom_type, config)

        if isinstance(geometry, Mesh):
            # Use terrain-specific converter for terrain objects
            if obj_data.get("type") == "TINRelief":
                geom_data = to_cityjson_terrain_mesh(
                    geometry, vertices, scale_factor, config
                )
            else:
                geom_data = to_cityjson_mesh(geometry, vertices, scale_factor, config)
            geom_data["lod"] = lod
            obj_data["geometry"].append(geom_data)

        elif isinstance(geometry, Surface):
            geom_data = to_cityjson_surface(geometry, vertices, scale_factor, config)
            geom_data["lod"] = lod
            obj_data["geometry"].append(geom_data)

        elif isinstance(geometry, MultiSurface):
            geom_data = to_cityjson_multisurface(
                geometry, vertices, scale_factor, config
            )
            geom_data["lod"] = lod
            obj_data["geometry"].append(geom_data)


def save(city: City, path: Path, scale: float = 0.001, config: CityJSONConfig = None):
    """Save a city to a CityJSON file.

    Parameters
    ----------
    city : City
        The city to save.
    path : str or Path
        Path to the file.
    scale : float, optional
        Scale factor for coordinate compression (default 0.001)
    config : CityJSONConfig, optional
        Configuration for conversion settings

    Raises
    ------
    ValueError
        If the file format is not supported
    """
    cj = to_cityjson(city, scale=scale, config=config)
    path = Path(path)

    if path.suffix.lower() == ".json":
        with open(path, "w") as file:
            json.dump(cj, file, indent=2)
    else:
        raise ValueError(
            f"Unsupported file format: {path.suffix}. Only .json is supported."
        )
