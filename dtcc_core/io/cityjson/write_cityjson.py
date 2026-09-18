import json
import math
import zipfile
from pathlib import Path
from typing import Dict, List

import numpy as np

from dtcc_core.model import (
    Building,
    BuildingPart,
    City,
    GeometryType,
    Mesh,
    MultiSurface,
    Surface,
    Terrain,
)

from .converters import (
    CityJSONConfig,
    VertexIndexer,
    convert_mesh,
    convert_multisurface,
    convert_surface,
    convert_terrain_mesh,
    geometry_type_to_lod,
    get_converter,
    get_terrain_converter,
)


def to_cityjson_surface(
    surface: Surface,
    vertices: list,
    scale: float,
    config: CityJSONConfig = None,
    indexer=None,
) -> dict:
    """Convert a DTCC Surface to CityJSON geometry format."""
    return convert_surface(surface, vertices, scale, config, indexer=indexer)


def to_cityjson_multisurface(
    multisurface: MultiSurface,
    vertices: list,
    scale: float,
    config: CityJSONConfig = None,
    indexer=None,
) -> dict:
    """Convert a DTCC MultiSurface to CityJSON geometry format."""
    return convert_multisurface(multisurface, vertices, scale, config, indexer=indexer)


def to_cityjson_mesh(
    mesh: Mesh,
    vertices: list,
    scale: float,
    config: CityJSONConfig = None,
    indexer=None,
) -> dict:
    """Convert a DTCC Mesh to CityJSON geometry format."""
    return convert_mesh(mesh, vertices, scale, config, indexer=indexer)


def to_cityjson_terrain_mesh(
    mesh: Mesh,
    vertices: list,
    scale: float,
    config: CityJSONConfig = None,
    indexer=None,
) -> dict:
    """Convert a DTCC terrain Mesh to CityJSON CompositeSurface format.

    According to CityJSON specification, TINRelief objects should use
    CompositeSurface geometry type with GroundSurface semantics.
    """
    return convert_terrain_mesh(mesh, vertices, scale, config, indexer=indexer)


def to_cityjson(
    city: City,
    scale: float = 0.001,
    config: CityJSONConfig = None,
    *,
    strict=False,
    validate_schema=None,
) -> dict:
    """Convert a DTCC City to CityJSON format.

    Strict mode evaluates the root's selected standard schema by default.
    validate_schema=False skips only semantics; explicit True requires strict=True.
    The native schema declaration is consulted but not encoded in CityJSON.

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
    from .admission import schema_validation

    validate_schema = schema_validation(strict, validate_schema)
    config = config or CityJSONConfig()
    if strict:
        from ...model._standard_schema import validate_admitted
        from ...model.exchange import _schema_for_model
        from .admission import validate_export

        validate_export(city)
        schema_id, schema_version = _schema_for_model(city)
        if validate_schema:
            validate_admitted(city, schema_id, schema_version)
        city.calculate_bounds()  # Public arrays may have changed since the cache was computed.

    # Validate scale
    if scale <= 0 or not math.isfinite(scale):
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
        if strict:
            crs = city.transform.srs
        if crs:
            cityjson["metadata"]["referenceSystem"] = crs

    if strict and city.transform.srs:
        cityjson.setdefault("metadata", {})["referenceSystem"] = city.transform.srs

    # Quantization factor for integer vertices
    quantize_factor = 1.0 / scale
    indexer = VertexIndexer()
    vertices = indexer.vertices  # keep alias for converters that expect a list

    # Use the same recursive traversal for legacy and strict exports.
    def add_feature(feature, parent=None):
        if strict:
            from .admission import feature_type

            name = feature_type(feature)
        else:
            name = type(feature).__name__
        from .attributes import write_attributes

        data = {
            "type": name,
            "attributes": write_attributes(feature.attributes, name),
            "geometry": [],
        }
        if parent is not None:
            data["parents"] = [parent.id]
        parts = sorted(feature.children.get(BuildingPart, []), key=lambda p: p.id)
        if parts:
            data["children"] = [part.id for part in parts]
        if strict:
            from ...model import Mesh, MultiLineString, Point, Solid
            from .converters import convert_mesh, convert_multisurface

            for record in feature.geometry.values():
                geometry = record.geometry
                if type(geometry) is Point:
                    indices = indexer.add_points(
                        np.array([[geometry.x, geometry.y, geometry.z]]),
                        quantize_factor,
                        config.rounding_mode,
                    )
                    data["geometry"].append(
                        {"type": "MultiPoint", "lod": record.lod, "boundaries": indices}
                    )
                    continue
                if type(geometry) is MultiLineString:
                    boundaries = []
                    for line in geometry.linestrings:
                        indices = indexer.add_points(
                            line.vertices, quantize_factor, config.rounding_mode
                        )
                        if (
                            len(set(indices)) < 2
                            and len(np.unique(line.vertices, axis=0)) >= 2
                        ):
                            raise ValueError(
                                "Requested CityJSON quantization collapses a line"
                            )
                        boundaries.append(indices)
                    data["geometry"].append(
                        {
                            "type": "MultiLineString",
                            "lod": record.lod,
                            "boundaries": boundaries,
                        }
                    )
                    continue
                if type(geometry) is Mesh:
                    encoded = convert_mesh(
                        geometry, vertices, quantize_factor, config, indexer=indexer
                    )
                    encoded["type"] = "CompositeSurface"
                    encoded["lod"] = record.lod
                    # No relief surface vocabulary is asserted in this subset.
                    del encoded["semantics"]
                    if any(
                        len(set(polygon[0])) != 3 for polygon in encoded["boundaries"]
                    ):
                        raise ValueError(
                            "Requested CityJSON quantization collapses a TINRelief triangle"
                        )
                    data["geometry"].append(encoded)
                    continue
                # The converter uses the shared polygon/region representation;
                # Solid shells are restored explicitly below.
                encoded = convert_multisurface(
                    geometry, vertices, quantize_factor, config, indexer=indexer
                )
                encoded["lod"] = record.lod
                from .admission import COMPOSITE_SURFACE_ROLE

                if record.role == COMPOSITE_SURFACE_ROLE:
                    encoded["type"] = "CompositeSurface"
                for surface, polygon in zip(geometry.surfaces, encoded["boundaries"]):
                    for source_ring, ring in zip(
                        [surface.vertices, *surface.holes], polygon
                    ):
                        if (
                            len(set(ring)) < 3
                            and len(np.unique(source_ring, axis=0)) >= 3
                        ):
                            raise ValueError(
                                "Requested CityJSON quantization collapses a polygon ring"
                            )
                if not geometry.regions:
                    encoded["semantics"] = {
                        "surfaces": [],
                        "values": [None] * len(geometry.surfaces),
                    }
                if type(geometry) is Solid:
                    encoded["type"] = "Solid"
                    flat = encoded["boundaries"]
                    encoded["boundaries"] = [
                        [flat[i] for i in shell] for shell in geometry.shells
                    ]
                    values = encoded["semantics"]["values"]
                    encoded["semantics"]["values"] = [
                        [values[i] for i in shell] for shell in geometry.shells
                    ]
                data["geometry"].append(encoded)
        else:
            _add_object_geometries(
                feature, data, vertices, quantize_factor, config, indexer=indexer
            )
        cityjson["CityObjects"][feature.id] = data
        for part in parts:
            add_feature(part, feature)

    roots = (
        [child for group in city.children.values() for child in group]
        if strict
        else city.children.get(Building, [])
    )
    for feature in sorted(roots, key=lambda b: b.id):
        add_feature(feature)

    # Process terrain deterministically
    if not strict and Terrain in city.children:
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
                terrain,
                terrain_data,
                vertices,
                quantize_factor,
                config,
                indexer=indexer,
            )

            cityjson["CityObjects"][terrain.id] = terrain_data

    # Add all processed vertices to the final structure
    cityjson["vertices"] = indexer.vertices

    return cityjson


def _add_object_geometries(
    obj,
    obj_data: dict,
    vertices: list,
    scale_factor: float,
    config: CityJSONConfig = None,
    indexer=None,
):
    """Add geometries from a DTCC object to CityJSON object data."""
    config = config or CityJSONConfig()

    # Iterate deterministically by geometry type
    for key, record in sorted(obj.geometry.items()):
        geometry = record.geometry
        geom_type = (
            GeometryType.from_str("lod" + record.lod)
            if record.lod in ("0", "1", "2", "3")
            else record.role
        )
        if record.lod is not None and record.lod not in ("0", "1", "2", "3"):
            raise NotImplementedError("Fractional LoD requires strict CityJSON export")
        if record.role is not None:
            try:
                geom_type = GeometryType.from_str(record.role)
            except ValueError:
                raise NotImplementedError(
                    f"No legacy CityJSON mapping for role {record.role!r}"
                )
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
    strict: bool = False,
    validate_schema=None,
):
    """Save a city to a CityJSON file.

    Strict mode validates the selected standard schema before opening the output.
    Pass validate_schema=False to skip semantics, retaining strict format checks.

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
    cj = to_cityjson(
        city, scale=scale, config=config, strict=strict, validate_schema=validate_schema
    )
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

        json_bytes = json.dumps(cj, indent=indent, ensure_ascii=ensure_ascii).encode(
            "utf-8"
        )
        with zipfile.ZipFile(path, mode="w", compression=zipfile.ZIP_DEFLATED) as zf:
            zf.writestr(inner_name, json_bytes)
    else:
        raise ValueError(
            f"Unsupported file format: {path.suffix}. Only .json and .json.zip are supported."
        )
