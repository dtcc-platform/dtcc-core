from .utils import get_terrain_mesh, get_root_objects, set_buildings
from ...model.object.city import City
from ...model.object.object import GeometryType
from ...model.object.terrain import Terrain
from ...model.geometry import Bounds
import numpy as np
from pathlib import Path
import json
import zipfile
from .attributes import _AttributeNameError


def _read_json(stream, strict):
    if not strict:
        return json.load(stream)
    from ...model.exchange import MAX_BYTES

    data = stream.read(MAX_BYTES + 1)
    if len(data) > MAX_BYTES:
        raise ValueError("CityJSON input exceeds 256 MiB limit")

    def unique_pairs(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"Duplicate CityJSON key {key!r}")
            result[key] = value
        return result

    return json.loads(data, object_pairs_hook=unique_pairs)


def setup_city(cj_obj: dict):
    """Set up city and transform vertices."""
    if not isinstance(cj_obj, dict):
        raise ValueError("CityJSON object must be a dictionary")

    city = City()

    # Handle coordinate transformation
    if "transform" in cj_obj:
        transform = cj_obj["transform"]
        if "scale" in transform and "translate" in transform:
            scale = np.array(transform["scale"])
            translate = np.array(transform["translate"])
        else:
            raise ValueError("Invalid transform: missing scale or translate")
    else:
        scale = np.array([1.0, 1.0, 1.0])
        translate = np.array([0.0, 0.0, 0.0])

    # Handle metadata and bounds
    if "metadata" in cj_obj:
        metadata = cj_obj["metadata"]
        if "geographicalExtent" in metadata:
            extent = metadata["geographicalExtent"]
            if len(extent) >= 6:
                try:
                    city.bounds = Bounds(
                        xmin=float(extent[0]),
                        ymin=float(extent[1]),
                        zmin=float(extent[2]),
                        xmax=float(extent[3]),
                        ymax=float(extent[4]),
                        zmax=float(extent[5]),
                    )
                except (ValueError, IndexError) as e:
                    print(f"Warning: Invalid geographical extent: {e}")

    # Transform vertices
    if "vertices" not in cj_obj:
        raise ValueError("Missing vertices in CityJSON")

    try:
        vertices = np.array(cj_obj["vertices"])
        if vertices.size == 0:
            raise ValueError("Empty vertices array")
        verts = vertices * scale + translate
    except (ValueError, TypeError) as e:
        raise ValueError(f"Invalid vertices format: {e}")

    return city, verts


def load(
    cityjson_path: str | dict,
    *,
    strict=False,
    extent_policy="validate",
    validate_schema=None,
) -> City:
    """Load CityJSON; strict mode evaluates the default standard schema.

    validate_schema=False bypasses semantics only. Omission follows strict mode;
    explicit True requires strict=True. Strict imports record the selected schema
    on the City; the source CityJSON does not declare a DTCC schema version.
    """
    from .admission import schema_validation

    validate_schema = schema_validation(strict, validate_schema)
    if not strict and extent_policy != "validate":
        raise ValueError("extent_policy requires strict CityJSON import")
    try:
        if isinstance(cityjson_path, dict):
            cj = cityjson_path
        else:
            cityjson_path = Path(cityjson_path)
            if not cityjson_path.exists():
                raise FileNotFoundError(f"File {cityjson_path} not found")

            if cityjson_path.suffix == ".zip":
                with zipfile.ZipFile(cityjson_path, "r") as z:
                    files = z.namelist()
                    if len(files) != 1 or not files[0].endswith(".json"):
                        raise ValueError(
                            "Invalid CityJSON zip file: must contain exactly one .json file"
                        )
                    with z.open(files[0]) as f:
                        cj = _read_json(f, strict)
            else:
                with open(cityjson_path, "rb") as f:
                    cj = _read_json(f, strict)

        if strict:
            from .admission import load_city
            from ...model.exchange import _schema_for_model
            from ...model._standard_schema import validate_admitted

            city = load_city(cj, extent_policy=extent_policy)
            schema_id, schema_version = _schema_for_model(city)
            if validate_schema:
                validate_admitted(city, schema_id, schema_version)
            city.schema_id, city.schema_version = schema_id, schema_version
            return city

        # Validate CityJSON structure
        if not isinstance(cj, dict):
            raise ValueError("Invalid CityJSON: root must be a dictionary")
        if cj.get("type") != "CityJSON":
            raise ValueError(
                "Not a valid CityJSON file: missing or incorrect 'type' field"
            )
        if "CityObjects" not in cj:
            raise ValueError("Invalid CityJSON: missing 'CityObjects' field")
        if "vertices" not in cj:
            raise ValueError("Invalid CityJSON: missing 'vertices' field")

        city, verts = setup_city(cj)
        cj_obj = cj["CityObjects"]

        root_objects = get_root_objects(cj_obj)

        # Process buildings
        if "Building" in root_objects:
            try:
                set_buildings(cj_obj, root_objects["Building"], verts, city)
            except _AttributeNameError:
                raise
            except Exception as e:
                print(f"Warning: Error processing buildings: {e}")

        # Process terrain
        if "TINRelief" in root_objects:
            try:
                tin = get_terrain_mesh(root_objects["TINRelief"], verts)
                if len(tin.vertices) > 0:
                    terrain = Terrain()
                    terrain.add_geometry(tin, GeometryType.MESH)
                    city.children[Terrain].append(terrain)
            except Exception as e:
                print(f"Warning: Error processing terrain: {e}")

        # Handle unsupported object types
        supported_types = {"Building", "TINRelief"}
        for obj_type, objects in root_objects.items():
            if obj_type not in supported_types:
                print(
                    f"Warning: Object type '{obj_type}' not yet supported ({len(objects)} objects)"
                )

        return city

    except json.JSONDecodeError as e:
        raise ValueError(f"Invalid JSON format: {e}")
    except (KeyError, TypeError) as e:
        raise ValueError(f"Invalid CityJSON structure: {e}")
