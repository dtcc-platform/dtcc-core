import json
import zipfile
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

from .converters import (
    get_converter,
    get_terrain_converter,
    geometry_type_to_lod,
    CityJSONConfig,
    convert_surface,
    convert_multisurface,
    convert_mesh,
    convert_terrain_mesh,
    VertexIndexer,
)


def to_cityjson_surface(
    surface: Surface,
    vertices: List,
    scale: float,
    config: CityJSONConfig = None,
    indexer=None,
) -> Dict:
    """Convert a DTCC Surface to CityJSON geometry format."""
    return convert_surface(surface, vertices, scale, config, indexer=indexer)


def to_cityjson_multisurface(
    multisurface: MultiSurface,
    vertices: List,
    scale: float,
    config: CityJSONConfig = None,
    indexer=None,
) -> Dict:
    """Convert a DTCC MultiSurface to CityJSON geometry format."""
    return convert_multisurface(multisurface, vertices, scale, config, indexer=indexer)


def to_cityjson_mesh(
    mesh: Mesh,
    vertices: List,
    scale: float,
    config: CityJSONConfig = None,
    indexer=None,
) -> Dict:
    """Convert a DTCC Mesh to CityJSON geometry format."""
    return convert_mesh(mesh, vertices, scale, config, indexer=indexer)


def to_cityjson_terrain_mesh(
    mesh: Mesh,
    vertices: List,
    scale: float,
    config: CityJSONConfig = None,
    indexer=None,
) -> Dict:
    """Convert a DTCC terrain Mesh to CityJSON CompositeSurface format.

    According to CityJSON specification, TINRelief objects should use
    CompositeSurface geometry type with GroundSurface semantics.
    """
    return convert_terrain_mesh(mesh, vertices, scale, config, indexer=indexer)


def to_cityjson(
    city: City, scale: float = 0.001, config: CityJSONConfig = None
) -> Dict:
    """Convert a DTCC City to CityJSON format.

    Parameters
    ----------
    city : City
        The DTCC City object to convert
    scale : float, optional
        CityJSON transform scale (default 0.001). Quantization factor is 1/scale.
    config : CityJSONConfig, optional
        Configuration for conversion settings

    Returns
    -------
    Dict
        CityJSON formatted dictionary
    """
    config = config or CityJSONConfig()

    # Validate scale
    if scale <= 0 or not (scale == scale):
        raise ValueError("Transform scale must be a positive finite number")

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
            "geographicalExtent": [
                float(city.bounds.xmin),
                float(city.bounds.ymin),
                float(city.bounds.zmin),
                float(city.bounds.xmax),
                float(city.bounds.ymax),
                float(city.bounds.zmax),
            ],
        }
        # Optional CRS if available on city
        try:
            crs = getattr(city, "crs", None)
        except Exception:
            crs = None
        if not crs and hasattr(city, "attributes"):
            crs = city.attributes.get("crs")
        if crs:
            cityjson["metadata"]["referenceSystem"] = crs

    # Quantization factor for integer vertices
    quantize_factor = 1.0 / scale
    indexer = VertexIndexer()
    vertices = indexer.vertices  # keep alias for converters that expect a list

    # Process buildings deterministically
    if Building in city.children:
        buildings = list(city.children[Building])
        buildings.sort(key=lambda b: getattr(b, "id", ""))

        for building in buildings:
            building_data = {
                "type": "Building",
                "attributes": building.attributes.copy(),
                "geometry": [],
            }

            # Process BuildingParts deterministically
            if building.children and BuildingPart in building.children:
                building_data["children"] = []
                parts = list(building.children[BuildingPart])
                parts.sort(key=lambda p: getattr(p, "id", ""))
                for part in parts:
                    part_data = {
                        "type": "BuildingPart",
                        "attributes": part.attributes.copy(),
                        "geometry": [],
                        "parents": [building.id],
                    }

                    # Add geometries from the part
                    _add_object_geometries(
                        part, part_data, vertices, quantize_factor, config, indexer=indexer
                    )

                    cityjson["CityObjects"][part.id] = part_data
                    building_data["children"].append(part.id)

            # Add geometries directly attached to the building
            _add_object_geometries(
                building, building_data, vertices, quantize_factor, config, indexer=indexer
            )

            cityjson["CityObjects"][building.id] = building_data

    # Process terrain deterministically
    if Terrain in city.children:
        terrains = list(city.children[Terrain])
        terrains.sort(key=lambda t: getattr(t, "id", ""))
        for terrain in terrains:
            terrain_data = {
                "type": "TINRelief",
                "attributes": terrain.attributes.copy(),
                "geometry": [],
            }

            # Add terrain geometries
            _add_object_geometries(
                terrain, terrain_data, vertices, quantize_factor, config, indexer=indexer
            )

            cityjson["CityObjects"][terrain.id] = terrain_data

    # Add all processed vertices to the final structure
    cityjson["vertices"] = indexer.vertices

    return cityjson


def _add_object_geometries(
    obj,
    obj_data: Dict,
    vertices: List,
    scale_factor: float,
    config: CityJSONConfig = None,
    indexer=None,
):
    """Add geometries from a DTCC object to CityJSON object data."""
    config = config or CityJSONConfig()

    # Iterate deterministically by geometry type
    for geom_type in sorted(list(obj.geometry.keys()), key=lambda gt: getattr(gt, "name", str(gt))):
        geometry = obj.geometry.get(geom_type)
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
                    geometry, vertices, scale_factor, config, indexer=indexer
                )
            else:
                geom_data = to_cityjson_mesh(
                    geometry, vertices, scale_factor, config, indexer=indexer
                )
            geom_data["lod"] = lod
            obj_data["geometry"].append(geom_data)

        elif isinstance(geometry, Surface):
            geom_data = to_cityjson_surface(
                geometry, vertices, scale_factor, config, indexer=indexer
            )
            geom_data["lod"] = lod
            obj_data["geometry"].append(geom_data)

        elif isinstance(geometry, MultiSurface):
            geom_data = to_cityjson_multisurface(
                geometry, vertices, scale_factor, config, indexer=indexer
            )
            geom_data["lod"] = lod
            obj_data["geometry"].append(geom_data)


def save(
    city: City,
    path: Path,
    scale: float = 0.001,
    config: CityJSONConfig = None,
    indent: int | None = 2,
    ensure_ascii: bool = False,
):
    """Save a city to a CityJSON file.

    Parameters
    ----------
    city : City
        The city to save.
    path : str or Path
        Path to the file.
    scale : float, optional
        CityJSON transform scale (default 0.001). Quantization factor is 1/scale.
    config : CityJSONConfig, optional
        Configuration for conversion settings
    indent : int | None, optional
        Indentation for pretty JSON. Use None for compact JSON.
    ensure_ascii : bool, optional
        If False (default), write UTF-8 characters directly. If True, escape non-ASCII.

    Raises
    ------
    ValueError
        If the file format is not supported
    """
    cj = to_cityjson(city, scale=scale, config=config)
    path = Path(path)

    suffix = path.suffix.lower()
    two_level = "".join([s.lower() for s in path.suffixes[-2:]])

    if suffix == ".json":
        with open(path, "w", encoding="utf-8") as file:
            json.dump(cj, file, indent=indent, ensure_ascii=ensure_ascii)
    elif suffix == ".zip" or two_level == ".json.zip":
        # Determine inner JSON name
        if two_level == ".json.zip":
            # Drop the trailing .json.zip and add .json
            base = path.name
            inner_name = base[: -len(".json.zip")] + ".json"
        else:
            inner_name = path.with_suffix("").name + ".json"

        json_bytes = json.dumps(cj, indent=indent, ensure_ascii=ensure_ascii).encode("utf-8")
        with zipfile.ZipFile(path, mode="w", compression=zipfile.ZIP_DEFLATED) as zf:
            zf.writestr(inner_name, json_bytes)
    else:
        raise ValueError(
            f"Unsupported file format: {path.suffix}. Only .json and .json.zip are supported."
        )
