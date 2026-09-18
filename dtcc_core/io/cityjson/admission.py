"""Explicit faithful subset of CityJSON 2.0, shared by its load/export boundary.

This is representation admission, not a substitute for a semantic profile or a
complete CityJSON conformance validator. Geometry conversion stays in utils and
converters. Unsupported content is rejected before it can be silently discarded.
"""

import re

import numpy as np

from ...model import (
    Building,
    BuildingPart,
    City,
    Landuse,
    Mesh,
    MultiLineString,
    MultiSurface,
    Object,
    Point,
    Solid,
    Terrain,
    exchange,
)
from .attributes import read_attributes
from .semantics import SEMANTIC_NAMESPACE, validate_region_owner
from .utils import build_geometry, build_tin

# CityJSON 2.0.2 JSON Schema allows integer 0–3 and decimal sublevels 0–3.
LOD = re.compile(r"[0-3](?:\.[0-3])?\Z")

# External-format mappings only: no profile inference or source code-list guesses.
FEATURE_TYPES = {
    "Building": Building,
    "BuildingPart": BuildingPart,
    "LandUse": Landuse,
    "SolitaryVegetationObject": Object,
    "CityFurniture": Object,
    "TINRelief": Terrain,
    **{
        name: Object
        for name in (
            "Road",
            "Railway",
            "Waterway",
            "TransportSquare",
            "WaterBody",
            "PlantCover",
        )
    },
}
FEATURE_GEOMETRIES = {
    "Building": (MultiSurface, Solid),
    "BuildingPart": (MultiSurface, Solid),
    "LandUse": (MultiSurface,),
    "SolitaryVegetationObject": (MultiSurface, Solid, Point),
    "CityFurniture": (MultiSurface, Solid, Point),
    "TINRelief": (Mesh,),
    **{
        name: (MultiLineString, MultiSurface)
        for name in ("Road", "Railway", "Waterway", "TransportSquare")
    },
    "WaterBody": (MultiLineString, MultiSurface, Solid),
    "PlantCover": (MultiSurface, Solid),
}


# Preserve the source polygon aggregate's distinct meaning on the generic carrier.
COMPOSITE_SURFACE_ROLE = "cityjson:CompositeSurface"


def schema_validation(strict, validate_schema):
    """An omitted option follows mode; explicit validation requires strict admission."""
    if validate_schema is None:
        return bool(strict)
    exchange._schema_flag(validate_schema)
    if validate_schema and not strict:
        raise ValueError("validate_schema=True requires strict=True for CityJSON")
    return validate_schema


def feature_type(value):
    """Resolve the bounded external spelling without downcasting native state."""
    for name, cls in FEATURE_TYPES.items():
        if type(value) is cls and (
            value.semantic_type == SEMANTIC_NAMESPACE + name
            or (cls not in (Object, Terrain) and value.semantic_type is None)
        ):
            return name
    raise NotImplementedError(
        "Object type/semantic URI needs an explicit CityJSON mapping"
    )


def _keys(value, allowed, where):
    if not isinstance(value, dict):
        raise ValueError(f"{where} must be an object")
    extra = set(value) - allowed
    if extra:
        raise NotImplementedError(f"{where}: unsupported content {sorted(extra)}")


def load_city(cj, *, extent_policy="validate"):
    if extent_policy not in ("validate", "recompute"):
        raise ValueError("extent_policy must be validate or recompute")
    _keys(
        cj,
        {"type", "version", "transform", "vertices", "CityObjects", "metadata"},
        "CityJSON",
    )
    if cj.get("type") != "CityJSON" or cj.get("version") != "2.0":
        raise ValueError("Strict import requires CityJSON version 2.0")
    exchange._validate_attributes(cj, max_depth=exchange.MAX_DEPTH)
    transform = cj.get("transform")
    _keys(transform, {"scale", "translate"}, "CityJSON transform")
    for key in ("scale", "translate"):
        vector = transform.get(key)
        if (
            not isinstance(vector, list)
            or len(vector) != 3
            or any(type(x) not in (int, float) for x in vector)
        ):
            raise ValueError("CityJSON scale/translate must have three finite numbers")
    scale, translate = (
        np.asarray(transform["scale"]),
        np.asarray(transform["translate"]),
    )
    for vector in (scale, translate):
        if (
            vector.shape != (3,)
            or vector.dtype.kind not in "ifu"
            or not np.isfinite(vector).all()
        ):
            raise ValueError("CityJSON scale/translate must have three finite numbers")
    if (scale <= 0).any():
        raise ValueError("CityJSON scale must be positive")
    source = cj.get("vertices")
    if not isinstance(source, list):
        raise ValueError("CityJSON vertices must be an array")
    for vertex in source:
        if (
            not isinstance(vertex, list)
            or len(vertex) != 3
            or any(type(x) is not int for x in vertex)
        ):
            raise ValueError("CityJSON vertices must be integer triples")
    vertices = np.asarray(source)
    if vertices.size == 0 and source == []:
        vertices = np.empty((0, 3), dtype=np.int64)
    if vertices.ndim != 2 or vertices.shape[1] != 3 or vertices.dtype.kind not in "iu":
        raise ValueError("CityJSON vertices must be integer triples")
    # Avoid silently rounding integer source coordinates before dequantization.
    for coordinate in vertices.flat:
        exchange._require_exact_double(coordinate)
    vertices = vertices.astype(np.float64) * scale + translate
    metadata = cj.get("metadata", {})
    _keys(metadata, {"referenceSystem", "geographicalExtent"}, "CityJSON metadata")
    city = City(semantic_type=SEMANTIC_NAMESPACE + "City")
    city.transform.srs = metadata.get("referenceSystem", "")
    extent = metadata.get("geographicalExtent")
    if extent is not None and (
        not isinstance(extent, list)
        or len(extent) != 6
        or any(type(x) not in (int, float) for x in extent)
    ):
        raise ValueError("CityJSON geographicalExtent must contain six finite numbers")
    objects = cj.get("CityObjects")
    if not isinstance(objects, dict):
        raise ValueError("CityObjects must be an object")
    native = {}
    parents = {}
    for id, obj in objects.items():
        _keys(
            obj,
            {
                "type",
                "attributes",
                "geometry",
                "children",
                "parents",
                "geographicalExtent",
            },
            id,
        )
        if not isinstance(obj.get("type"), str):
            raise ValueError(f"{id}: CityJSON feature type must be a string")
        if obj.get("type") not in FEATURE_TYPES:
            raise NotImplementedError(f"{id}: unsupported strict CityJSON feature type")
        value = FEATURE_TYPES[obj["type"]](
            id=id,
            attributes=read_attributes(obj.get("attributes", {}), obj["type"]),
            semantic_type=SEMANTIC_NAMESPACE + obj["type"],
        )
        native[id] = value
        value.transform.srs = city.transform.srs
        geometries = obj.get("geometry", [])
        if not isinstance(geometries, list):
            raise ValueError(f"{id}: geometry must be an array")
        for index, geometry in enumerate(geometries):
            _keys(
                geometry, {"type", "lod", "boundaries", "semantics"}, f"{id} geometry"
            )
            lod = geometry.get("lod")
            if not isinstance(lod, str) or not LOD.fullmatch(lod):
                raise ValueError("Strict import requires a CityJSON 2.0 LoD string")
            native_geometry = (
                build_tin(geometry, vertices)
                if obj["type"] == "TINRelief"
                else build_geometry(geometry, vertices)
            )
            if type(native_geometry) not in FEATURE_GEOMETRIES[obj["type"]]:
                raise NotImplementedError(
                    f"{id}: unsupported geometry for {obj['type']}"
                )
            validate_region_owner(native_geometry.regions, obj["type"])
            native_geometry.transform.srs = city.transform.srs
            role = (
                COMPOSITE_SURFACE_ROLE
                if type(native_geometry) is MultiSurface
                and geometry["type"] == "CompositeSurface"
                else None
            )
            value.add_geometry(
                native_geometry, id=f"cityjson-{index}", lod=lod, role=role
            )
        children = obj.get("children", [])
        if (
            not isinstance(children, list)
            or any(not isinstance(x, str) for x in children)
            or len(set(children)) != len(children)
        ):
            raise ValueError(f"{id}: children must be unique IDs")
        if children and obj["type"] not in ("Building", "BuildingPart"):
            raise NotImplementedError(
                "Strict containment currently supports building parts only"
            )
        if obj["type"] != "BuildingPart":
            if "parents" in obj:
                raise ValueError("Top-level CityJSON feature cannot have parents")
        else:
            parent = obj.get("parents")
            if (
                not isinstance(parent, list)
                or len(parent) != 1
                or not isinstance(parent[0], str)
            ):
                raise ValueError("BuildingPart requires exactly one parent ID")
            parents[id] = parent[0]
    for id, obj in objects.items():
        for child in obj.get("children", []):
            if child not in native or parents.get(child) != id:
                raise ValueError(f"{id}: dangling or inconsistent child {child!r}")
            native[id].add_child(native[child])
    for id, parent in parents.items():
        if parent not in native or id not in objects[parent].get("children", []):
            raise ValueError(f"{id}: dangling or inconsistent parent {parent!r}")
        if objects[parent]["type"] not in ("Building", "BuildingPart"):
            raise ValueError("BuildingPart parent must be a Building or BuildingPart")
    for id, value in native.items():
        if id not in parents:
            city.add_child(value)
    # A disconnected cycle would be invisible when walking the root alone.
    reachable, stack = set(), [city]
    while stack:
        value = stack.pop()
        if value.id in reachable:
            raise ValueError("CityJSON containment is cyclic or shared")
        reachable.add(value.id)
        stack.extend(child for group in value.children.values() for child in group)
    if len(reachable) != len(native) + 1:
        raise ValueError("CityJSON contains unreachable/cyclic objects")
    exchange.validate(city)
    discrepancies = []
    for source_id, extent, value in [
        (None, metadata.get("geographicalExtent"), city)
    ] + [
        (id, objects[id].get("geographicalExtent"), value)
        for id, value in native.items()
    ]:
        if _check_extent(extent, value, scale):
            if extent_policy == "validate":
                raise ValueError(
                    "CityJSON extent does not match represented geometry; use extent_policy='recompute' to retain the discrepancy in Dataset Context"
                )
            discrepancies.append(
                {"object_id": source_id, "source_extent": list(extent)}
            )
    if discrepancies:
        from ...datasets.schema import (
            DatasetContext,
            DatasetIdentity,
            DatasetMetadata,
            DatasetPresentation,
            DatasetProvenance,
            DatasetRequest,
        )

        city.dataset_context = DatasetContext(
            identity=DatasetIdentity(name="cityjson-import", title="CityJSON import"),
            metadata=DatasetMetadata(
                crs=[city.transform.srs] if city.transform.srs else []
            ),
            provenance=DatasetProvenance(
                processing_steps=[
                    {
                        "operation": "recompute source extent summaries from unchanged geometry",
                        "discrepancies": discrepancies,
                    }
                ]
            ),
            presentation=DatasetPresentation(),
            request=DatasetRequest(
                dataset_name="cityjson-import",
                parameters={"extent_policy": extent_policy},
            ),
            health={
                "status": "warning",
                "extent_discrepancy_count": len(discrepancies),
            },
            warnings=[
                f"Recomputed {len(discrepancies)} inconsistent CityJSON extent summaries; source summaries retained in provenance."
            ],
        )
    return city


def _check_extent(extent, value, scale):
    if extent is None:
        return
    if (
        not isinstance(extent, list)
        or len(extent) != 6
        or any(type(x) not in (int, float) for x in extent)
    ):
        raise ValueError("CityJSON geographicalExtent must contain six finite numbers")
    bounds = value.bounds
    if bounds is None:
        raise ValueError("CityJSON extent has no represented geometry")
    actual = [
        bounds.xmin,
        bounds.ymin,
        bounds.zmin,
        bounds.xmax,
        bounds.ymax,
        bounds.zmax,
    ]
    # Source summaries can predate coordinate quantization. Permit at most half
    # a source grid cell per axis; coordinates themselves are never changed.
    if not np.allclose(extent, actual, atol=np.tile(scale * 0.5, 2) + 1e-8, rtol=0):
        return True


def validate_export(city):
    exchange.validate(city)
    if type(city) is not City:
        raise ValueError("CityJSON export requires City")
    stack = [city]
    while stack:
        value = stack.pop()
        name = "City" if value is city else feature_type(value)
        if value is city and value.semantic_type not in (
            None,
            SEMANTIC_NAMESPACE + "City",
        ):
            raise NotImplementedError(
                "City semantic URI needs an explicit CityJSON extension mapping"
            )
        if type(value) is Landuse and value.landuses:
            raise NotImplementedError(
                "Native Landuse codes need an explicit CityJSON mapping"
            )
        if value.relations or value.profile_id is not None:
            raise NotImplementedError(
                "CityJSON export cannot preserve named relations/profile identity without an extension"
            )
        if value is city and (value.attributes or value.geometry):
            raise NotImplementedError(
                "City root attributes/geometry need a CityJSON mapping"
            )
        if value.transform.srs not in ("", city.transform.srs):
            raise NotImplementedError("CityJSON requires one CRS throughout the model")
        if not np.array_equal(value.transform.affine, np.eye(4)):
            raise NotImplementedError(
                "CityJSON export requires geometry in the global coordinate frame"
            )
        for group in value.children.values():
            for child in group:
                child_name = feature_type(child)
                if (value is city and child_name == "BuildingPart") or (
                    value is not city
                    and (
                        name not in ("Building", "BuildingPart")
                        or child_name != "BuildingPart"
                    )
                ):
                    raise NotImplementedError(
                        "Strict CityJSON containment supports top-level features and building parts"
                    )
                stack.append(child)
        for index, (key, record) in enumerate(value.geometry.items()):
            if key != f"cityjson-{index}" or (
                record.role is not None
                and not (
                    record.role == COMPOSITE_SURFACE_ROLE
                    and type(record.geometry) is MultiSurface
                )
            ):
                raise NotImplementedError(
                    "CityJSON representation IDs/roles require an extension mapping"
                )
            if record.lod is None or not LOD.fullmatch(record.lod):
                raise ValueError("Strict export requires a CityJSON 2.0 LoD string")
            geometry = record.geometry
            if type(geometry) not in FEATURE_GEOMETRIES[name]:
                raise NotImplementedError(
                    f"Unsupported strict CityJSON geometry for {name}"
                )
            validate_region_owner(geometry.regions, name)
            if type(geometry) is Mesh:
                if not len(geometry.faces):
                    raise ValueError("TINRelief requires nonempty triangle boundaries")
                if geometry.markers.size or geometry.normals.size:
                    raise NotImplementedError(
                        "Mesh markers/normals require an explicit CityJSON mapping"
                    )
                if len(np.unique(geometry.faces)) != len(geometry.vertices):
                    raise NotImplementedError(
                        "Unused mesh vertices have no strict TINRelief mapping"
                    )
            lines = geometry.linestrings if type(geometry) is MultiLineString else []
            if type(geometry) is MultiLineString and (
                not lines
                or any(
                    len(line.vertices) < 2 or line.vertices.shape[1] != 3
                    for line in lines
                )
            ):
                raise ValueError(
                    "CityJSON MultiLineString requires nonempty lines with at least two 3D vertices"
                )
            surfaces = (
                []
                if type(geometry) in (Point, Mesh, MultiLineString)
                else geometry.surfaces
            )
            if type(geometry) not in (Point, Mesh, MultiLineString) and (
                not surfaces or any(not len(s.vertices) for s in surfaces)
            ):
                raise ValueError("CityJSON MultiSurface boundaries must be nonempty")
            for item in (geometry, *surfaces, *lines):
                if (
                    item is not geometry
                    and type(geometry) is not MultiLineString
                    and item.normal.size
                ):
                    raise NotImplementedError(
                        "Stored Surface normals have no explicit CityJSON mapping"
                    )
                if item.fields or (item is not geometry and item.regions):
                    raise NotImplementedError(
                        "Nested fields/regions require an explicit CityJSON mapping"
                    )
                if not np.array_equal(item.transform.affine, np.eye(4)):
                    raise NotImplementedError(
                        "CityJSON export requires geometry in the global coordinate frame"
                    )
                if item.transform.srs not in ("", city.transform.srs):
                    raise NotImplementedError(
                        "CityJSON requires one CRS throughout the model"
                    )
