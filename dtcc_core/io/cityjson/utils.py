import numpy as np
from collections import defaultdict
from tqdm import tqdm

from ...model import Mesh, MultiSurface, Surface, Building, City
from ...model.object.city import CityObject
from ...model.object.building import BuildingPart
from ...model.object.object import GeometryType


def tin_geom_to_mesh(tin_geom_pair: tuple, verts) -> Mesh:
    """Convert TIN geometry to DTCC Mesh."""
    mesh = Mesh()
    faces = []
    uuid, tin_geom = tin_geom_pair

    # Collect all faces from geometry
    for geom in tin_geom.get("geometry", []):
        boundaries = geom.get("boundaries", [])
        faces.extend([face[0] for face in boundaries if face])

    if not faces:
        return mesh

    # Build vertex dictionary and remap indices
    vert_dict = {}
    for face in faces:
        for vertex_idx in face:
            if vertex_idx < len(verts):
                vert_dict[vertex_idx] = verts[vertex_idx]

    if not vert_dict:
        return mesh

    # Create sorted vertex mapping
    fv_map = sorted(vert_dict.items())
    vert_map = {old_idx: new_idx for new_idx, (old_idx, _) in enumerate(fv_map)}

    # Build mesh - ensure JSON serializable types
    mesh.vertices = np.array([vertex for _, vertex in fv_map])
    mesh.faces = np.array([[int(vert_map[vertex_idx]) for vertex_idx in face] for face in faces])

    return mesh


def get_terrain_mesh(tin_obj: dict, verts: np.ndarray) -> Mesh:
    """Extract terrain mesh from TIN objects."""
    if not tin_obj:
        return Mesh()

    if len(tin_obj) == 1:
        return tin_geom_to_mesh(tin_obj[0], verts)

    # For multiple TIN objects, combine them
    # TODO: Implement proper mesh merging
    # For now, return the first one
    return tin_geom_to_mesh(tin_obj[0], verts)


def build_multisurface(geom, verts):
    """Build a MultiSurface from a CityJSON geometry."""
    if not isinstance(geom, dict):
        raise ValueError("Geometry must be a dictionary")

    ms = MultiSurface()

    # Extract boundaries based on geometry type
    if geom["type"] in ("MultiSurface", "CompositeSurface"):
        boundaries = geom["boundaries"]
    elif geom["type"] == "Solid":
        boundaries = geom["boundaries"][0] if geom["boundaries"] else []
    else:
        raise ValueError(f"Unsupported geometry type: {geom['type']}")

    for surface_boundaries in boundaries:
        if not surface_boundaries:
            continue

        s = Surface()
        # First boundary is the outer ring
        outer_ring = surface_boundaries[0]
        s.vertices = np.array([verts[idx] for idx in outer_ring])

        # Remaining boundaries are holes
        for hole_ring in surface_boundaries[1:]:
            hole_vertices = np.array([verts[idx] for idx in hole_ring])
            s.holes.append(hole_vertices)

        ms.surfaces.append(s)

    return ms


def get_geom_to_use(geom_list, lod=2, prefer_surface=True) -> dict:
    """
    Get the geometry to use from a list of possible geometries.

    Parameters
    ----------
    geom_list : list
        List of possible geometries
    lod : int
        Preferred level of detail (2 by default)
    prefer_surface : bool
        If True, prefer surface geometries over solid

    Returns
    -------
    dict
        The selected geometry
    """
    if not geom_list:
        return {}

    if len(geom_list) == 1:
        return geom_list[0]

    # Find candidates by LOD preference (exact match, then decreasing)
    candidates = []
    current_lod = lod
    while not candidates and current_lod >= 0:
        candidates = [g for g in geom_list if g.get("lod") == current_lod]
        current_lod -= 1

    # If no LOD-specific candidates, use all
    if not candidates:
        candidates = geom_list

    if len(candidates) == 1:
        return candidates[0]

    # Apply type preference
    for geom in candidates:
        geom_type = geom.get("type", "")
        if "Surface" in geom_type and prefer_surface:
            return geom
        elif "Solid" in geom_type and not prefer_surface:
            return geom

    # Return first candidate as fallback
    return candidates[0]


def get_building_geometry(cj_obj, building, verts, lod=2):
    """Extract building geometry and children geometries."""
    if not isinstance(cj_obj, dict) or not isinstance(building, dict):
        raise ValueError("Invalid input: cj_obj and building must be dictionaries")

    building_root_geom = None
    building_children = []

    # Process root geometry
    if "geometry" in building:
        building_geometry = building["geometry"]
        if not isinstance(building_geometry, list):
            raise ValueError("Building geometry must be a list")

        geom = get_geom_to_use(building_geometry, lod=lod)
        if geom:  # Only process if we got a valid geometry
            ms = build_multisurface(geom, verts)
            building_root_geom = (geom, ms)

    # Process children geometries
    for child_id in building.get("children", []):
        if child_id not in cj_obj:
            print(f"Warning: Child {child_id} not found in CityObjects")
            continue

        child = cj_obj[child_id]
        if "geometry" not in child:
            continue

        child_geometry = child["geometry"]
        if not isinstance(child_geometry, list):
            continue

        geom = get_geom_to_use(child_geometry, lod=lod)
        if geom:
            ms = build_multisurface(geom, verts)
            building_children.append((geom, ms))

    return building_root_geom, building_children


def _convert_lod_to_geometry_type(lod_value):
    """Convert LOD value to GeometryType enum."""
    if isinstance(lod_value, str):
        lod_value = lod_value.split(".")[0]
    else:
        lod_value = str(lod_value).split(".")[0]
    return GeometryType.from_str(f"lod{lod_value}")


def build_dtcc_building(cj_obj, uuid, cj_building, verts, parent_city, lod=2):
    """Build a DTCC Building from CityJSON building data."""
    building = Building()
    if parent_city is not None:
        parent_city.children[Building].append(building)

    building.id = uuid
    building.attributes = cj_building.get("attributes", {})

    building_root_geom, building_children = get_building_geometry(
        cj_obj, cj_building, verts, lod=lod
    )

    # Add root geometry
    if building_root_geom is not None:
        geom, ms = building_root_geom
        lod_type = _convert_lod_to_geometry_type(geom.get("lod", 1))
        building.geometry[lod_type] = ms

    # Add building parts
    for geom, ms in building_children:
        lod_type = _convert_lod_to_geometry_type(geom.get("lod", 1))
        building_part = BuildingPart()
        building.children[BuildingPart].append(building_part)
        building_part.geometry[lod_type] = ms

    return building


def get_root_buildings(cj_obj: dict):
    """
    Extract top-level building entries from a CityJSON CityObject map.

    Parameters
    ----------
    cj_obj : dict
        Mapping of CityJSON object ids to object dictionaries.

    Returns
    -------
    list[tuple[str, dict]]
        List of ``(uuid, building_obj)`` pairs where objects have no parents and type ``\"Building\"``.
    """
    cj_root_buildings = []
    for k, v in cj_obj.items():
        if "parents" not in v:
            if v["type"] == "Building":
                cj_root_buildings.append((k, v))
    return cj_root_buildings


def get_root_objects(cj_obj: dict):
    """
    Group top-level CityJSON objects by type.

    Parameters
    ----------
    cj_obj : dict
        Mapping of CityJSON object ids to object dictionaries.

    Returns
    -------
    dict[str, list[tuple[str, dict]]]
        Dictionary keyed by CityJSON type with lists of ``(uuid, object)`` pairs that lack parents.
    """
    root_objects = defaultdict(list)
    for k, v in cj_obj.items():
        if "parents" not in v:
            root_objects[v["type"]].append((k, v))
    return root_objects


def set_buildings(
    cj_obj: dict,
    root_buildings: tuple[str, dict],
    verts: np.ndarray,
    parent_city: City,
    lod=2,
) -> [Building]:
    """
    Build DTCC buildings from root CityJSON entries and attach to a city.

    Parameters
    ----------
    cj_obj : dict
        CityJSON CityObject dictionary.
    root_buildings : iterable[tuple[str, dict]]
        Root building entries as ``(uuid, obj)`` pairs (e.g., output of ``get_root_buildings``).
    verts : np.ndarray
        Global vertex array used to construct geometries.
    parent_city : City
        City instance to which constructed buildings are added as children.
    lod : int, default 2
        Level of detail to extract when building geometries.

    Returns
    -------
    None
        Buildings are instantiated and appended to ``parent_city``; the function does not return a value.
    """
    buildings = []
    for uuid, v in tqdm(
        root_buildings, desc="Loading CityJson", unit=" building", ncols=120
    ):
        buildings.append(
            build_dtcc_building(cj_obj, uuid, v, verts, parent_city, lod=lod)
        )


def set_cityobject(
    cj_obj: dict,
    root_object: tuple[str, dict],
    verts: np.ndarray,
    parent_city: City,
    lod=2,
) -> list:
    """Set city objects from root objects (generic function for future extension)."""
    city_objects = []
    for uuid, v in root_object:
        # For now, treat all objects as buildings
        # TODO: Add support for other object types
        city_objects.append(
            build_dtcc_building(cj_obj, uuid, v, verts, parent_city, lod=lod)
        )
    return city_objects
