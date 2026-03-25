from ...model import Building, GeometryType, MultiSurface, Surface
from ..polygons.polygons import split_polygon_sides
from ..cleaning import ConditioningOptions, condition_building_footprints

from ..polygons.surface import clean_multisurface, clean_surface

from ..register import register_model_method
import shapely
from shapely.geometry import Polygon
from shapely.ops import unary_union
from shapely.validation import make_valid
from ..logging import debug, info, warning, error

from typing import List, Tuple, Union


@register_model_method
def get_footprint(building: Building, geom_type: GeometryType = None) -> Surface:
    """
    Extract footprint surface from building geometry.

    Parameters
    ----------
    building : Building
        The building to extract footprint from.
    geom_type : GeometryType, optional
        Specific geometry type to use. If None, uses highest available LOD.

    Returns
    -------
    Surface
        The building footprint as a surface, or None if no geometry found.
    """
    lod_levels = [
        GeometryType.LOD0,
        GeometryType.LOD1,
        GeometryType.LOD2,
        GeometryType.LOD3,
    ]

    geom = None

    if geom_type is not None:
        geom = building.flatten_geometry(geom_type)

    if geom is None:
        for lod in lod_levels:
            geom = building.flatten_geometry(lod)
            if geom is not None:
                break

    if geom is None:
        warning(f"Building {building.id} has no LOD geometry.")
        return None
    height = geom.bounds.zmax

    footprint = geom.to_polygon()
    if footprint.geom_type == "MultiPolygon":
        # Merge all polygons and return the largest by area
        if not footprint.geoms:
            warning(f"Building {building.id} MultiPolygon is empty.")
            return None
        merged = unary_union(footprint.geoms)
        if merged.geom_type == "MultiPolygon":
            largest = max(merged.geoms, key=lambda p: p.area)
            footprint = largest
        else:
            footprint = merged

    s = Surface()
    s.from_polygon(footprint, height)
    return s


def _extract_polygon_parts(geom) -> list[Polygon]:
    if geom is None or geom.is_empty:
        return []
    polygons: list[Polygon] = []
    stack = [geom]
    while stack:
        current = stack.pop()
        if current.is_empty:
            continue
        if isinstance(current, Polygon):
            polygons.append(current)
            continue
        if hasattr(current, "geoms"):
            stack.extend(reversed(list(current.geoms)))
    return [polygon for polygon in polygons if not polygon.is_empty and polygon.area > 0]


def _extract_building_metadata(
    building: Building,
    lod: GeometryType,
) -> tuple[float, float | None, float | None]:
    geometry = building.flatten_geometry(lod)
    roof_z = 0.0
    if geometry is not None:
        try:
            roof_z = float(geometry.bounds.zmax)
        except (AttributeError, TypeError):
            roof_z = 0.0

    try:
        height = float(building.height)
    except (AttributeError, TypeError):
        height = None
    if height is not None and height <= 0:
        height = None

    ground_height = building.attributes.get("ground_height")
    try:
        ground_height = float(ground_height) if ground_height is not None else None
    except (TypeError, ValueError):
        ground_height = None

    return roof_z, height, ground_height


def _area_weighted_value(
    source_indices: list[int],
    source_areas: list[float],
    values: list[float | None],
    *,
    default: float | None,
) -> float | None:
    weighted_values: list[float] = []
    weights: list[float] = []
    for index in source_indices:
        if index < 0 or index >= len(values):
            continue
        value = values[index]
        if value is None:
            continue
        weight = source_areas[index] if source_areas[index] > 0 else 1.0
        weighted_values.append(float(value))
        weights.append(float(weight))
    if not weighted_values:
        return default
    return float(
        sum(value * weight for value, weight in zip(weighted_values, weights))
        / sum(weights)
    )


def _build_conditioned_buildings(
    source_buildings: List[Building],
    polygons: list[Polygon],
    source_map: list[list[int]],
    *,
    lod: GeometryType,
) -> List[Building]:
    source_areas = [0.0] * len(source_buildings)
    source_roof_z: list[float | None] = [None] * len(source_buildings)
    source_height: list[float | None] = [None] * len(source_buildings)
    source_ground: list[float | None] = [None] * len(source_buildings)

    for index, building in enumerate(source_buildings):
        geometry = building.flatten_geometry(lod)
        if geometry is not None:
            polygon = geometry.to_polygon(simplify=0.0)
            if polygon is not None and not polygon.is_empty:
                source_areas[index] = float(max(polygon.area, 0.0))
        roof_z, height, ground = _extract_building_metadata(building, lod)
        source_roof_z[index] = roof_z
        source_height[index] = height
        source_ground[index] = ground

    conditioned_buildings: List[Building] = []
    for polygon, indices in zip(polygons, source_map):
        originals = [
            source_buildings[index]
            for index in indices
            if 0 <= index < len(source_buildings)
        ]
        roof_z = _area_weighted_value(
            indices,
            source_areas,
            source_roof_z,
            default=0.0,
        )
        height = _area_weighted_value(
            indices,
            source_areas,
            source_height,
            default=None,
        )
        ground_height = _area_weighted_value(
            indices,
            source_areas,
            source_ground,
            default=None,
        )

        surface = Surface()
        surface.from_polygon(polygon, roof_z or 0.0)

        building = Building()
        building.add_geometry(surface, GeometryType.LOD0)
        building.attributes = merge_building_attributes(originals)
        if height is not None:
            building.attributes["height"] = height
        if ground_height is not None:
            building.attributes["ground_height"] = ground_height
        conditioned_buildings.append(building)

    return conditioned_buildings


def _condition_buildings_with_shared_cleaner(
    buildings: List[Building],
    *,
    lod: GeometryType,
    options: ConditioningOptions,
    operation_name: str,
    return_index_map: bool = False,
) -> Union[List[Building], Tuple[List[Building], List[List[int]]]]:
    info(
        f"{operation_name}: running shared footprint conditioner on {len(buildings)} buildings."
    )
    result = condition_building_footprints(
        buildings,
        lod=lod,
        options=options,
    )
    conditioned_buildings = _build_conditioned_buildings(
        buildings,
        result.polygons,
        result.source_map,
        lod=lod,
    )
    info(
        f"{operation_name}: produced {len(conditioned_buildings)} conditioned buildings."
    )
    if return_index_map:
        return conditioned_buildings, result.source_map
    return conditioned_buildings


def merge_building_footprints(
    buildings: List[Building],
    lod: GeometryType = GeometryType.LOD0,
    max_distance: float = 0.5,
    min_area: float = 10.0,
    return_index_map: bool = False,
) -> List[Building] | Tuple[List[Building], List[List[int]]]:
    """
    Merge nearby building footprints into single buildings.

    Parameters
    ----------
    buildings : List[Building]
        List of buildings to merge.
    lod : GeometryType, default GeometryType.LOD0
        Level of detail to use for footprint extraction.
    max_distance : float, default 0.5
        Maximum distance between buildings to merge.
    min_area : float, default 10
        Minimum area threshold for merged footprints.
    return_index_map : bool, default False
        When True, also return the mapping from each merged building to the
        original building indices it represents.

    Returns
    -------
    Union[List[Building], Tuple[List[Building], List[List[int]]]]
        Merged buildings, optionally paired with the index map.
    """
    return _condition_buildings_with_shared_cleaner(
        buildings,
        lod=lod,
        options=ConditioningOptions(
            precision_grid=None,
            min_feature_size=0.0,
            merge_distance=max_distance,
            min_area=min_area,
            min_hole_area=0.0,
        ),
        operation_name="merge_building_footprints",
        return_index_map=return_index_map,
    )


def merge_building_attributes(buildings: List[Building]) -> dict:
    """
    Merge attributes from multiple buildings into a single dictionary.

    Parameters
    ----------
    buildings : List[Building]
        List of buildings whose attributes to merge.

    Returns
    -------
    dict
        Merged attributes dictionary.
    """
    attributes = {}
    for building in buildings:
        for k, v in building.attributes.items():
            if v:
                attributes[k] = v
    return attributes


def simplify_building_footprints(
    buildings: List[Building],
    tolerance: float = 0.5,
    method: str = "vwp",
    lod: GeometryType = GeometryType.LOD0,
    return_index_map: bool = False,
) -> Union[List[Building], Tuple[List[Building], List[List[int]]]]:
    """
    Simplify the building footprints by reducing the number of vertices while maintaining the overall shape.

    Parameters
    ----------
    buildings : List[Building]
        A list of `Building` objects whose footprints need to be simplified.
    tolerance : float, optional
        The tolerance for simplification. A higher value results in a more simplified footprint (default is 0.5).
    method : str, optional
        The simplification method to use. Options are 'rdp' for Ramer-Douglas-Peucker algorithm, 'vw' for
        Visvalingam-Whyatt algorithm or 'vwp' for topology preserving Visvalingam-Whyatt (default is 'vwp').
    lod : GeometryType, optional
        The level of detail of the geometry to simplify. Typically set to `GeometryType.LOD0` (default).
    return_index_map : bool, optional
        When True, also return a mapping from each simplified building to its
        original index in `buildings`.

    Returns
    -------
    Union[List[Building], Tuple[List[Building], List[List[int]]]]
        Simplified buildings, optionally paired with the index map.
    """

    if method not in ["vwp", "vw", "rdp"]:
        warning(
            f"Unknown polygon simplification method: {method}. "
            "The shared footprint conditioner ignores this argument."
        )
    if tolerance < 0:
        raise ValueError("tolerance must be non-negative.")

    info(
        "simplify_building_footprints() delegates to the shared conditioner; "
        "the 'method' argument is retained for API compatibility only."
    )
    return _condition_buildings_with_shared_cleaner(
        buildings,
        lod=lod,
        options=ConditioningOptions(
            precision_grid=None,
            min_feature_size=tolerance,
            merge_distance=0.0,
            min_area=0.0,
            min_hole_area=0.0,
        ),
        operation_name="simplify_building_footprints",
        return_index_map=return_index_map,
    )


def clean_building_footprints(
    buildings: List[Building],
    clearance: float = 0.5,
    smallest_hole_area: float = 1.0,
    return_index_map: bool = False,
) -> Union[List[Building], Tuple[List[Building], List[List[int]]]]:
    """
    Clean building footprints by removing overlaps, small holes, and ensuring clearance.

    Parameters
    ----------
    buildings : List[Building]
        List of buildings to clean.
    clearance : float, default 0.5
        Minimum clearance distance in meters.
    remove_overlaps : bool, default True
        Whether to remove overlapping footprints.
    smallest_hole_area : float, default 1.0
        Minimum area of holes to keep.
    Returns : List[Building]
        List of cleaned buildings.
    """

    return _condition_buildings_with_shared_cleaner(
        buildings,
        lod=GeometryType.LOD0,
        options=ConditioningOptions(
            precision_grid=None,
            min_feature_size=clearance,
            merge_distance=0.0,
            min_area=0.0,
            min_hole_area=smallest_hole_area,
        ),
        operation_name="clean_building_footprints",
        return_index_map=return_index_map,
    )


def fix_building_footprint_clearance(
    buildings: List[Building],
    clearance: float = 0.5,
    lod: GeometryType = GeometryType.LOD0,
    return_index_map: bool = False,
) -> Union[List[Building], Tuple[List[Building], List[List[int]]]]:
    """
    Fix clearance issues in building footprints.

    Parameters
    ----------
    buildings : List[Building]
        List of buildings to fix.
    clearance : float, default 0.5
        Minimum clearance distance in meters.
    lod : GeometryType, default GeometryType.LOD0
        Level of detail to fix.
    return_index_map : bool, optional
        When True, also return a mapping from each fixed building to its
        original index in `buildings`.

    Returns
    -------
    Union[List[Building], Tuple[List[Building], List[List[int]]]]
        Buildings with fixed clearances, optionally paired with the index map.
    """
    return _condition_buildings_with_shared_cleaner(
        buildings,
        lod=lod,
        options=ConditioningOptions(
            precision_grid=None,
            min_feature_size=clearance,
            merge_distance=0.0,
            min_area=0.0,
            min_hole_area=0.0,
        ),
        operation_name="fix_building_footprint_clearance",
        return_index_map=return_index_map,
    )


def split_footprint_walls(
    buildings: List[Building], max_wall_length: Union[float, List[float]] = 10.0
) -> List[Building]:
    """
    Split long walls in building footprints into shorter segments.

    Parameters
    ----------
    buildings : List[Building]
        List of buildings to process.
    max_wall_length : Union[float, List[float]], default 10
        Maximum wall length in meters. Can be single value or list per building.

    Returns
    -------
    List[Building]
        List of buildings with split walls.
    """
    split_buildings = []
    if isinstance(max_wall_length, (int, float)):
        max_wall_length = [max_wall_length] * len(buildings)
    elif len(max_wall_length) != len(buildings):
        error(
            "max_wall_length must be a single value or a list of values for each building."
        )
        return
    for building, wall_length in zip(buildings, max_wall_length):
        lod0 = building.lod0
        if lod0 is None:
            continue

        footprint = lod0.to_polygon()
        if footprint is None or footprint.is_empty:
            continue
        footprint = split_polygon_sides(footprint, wall_length)
        building_surface = Surface()
        building_surface.from_polygon(footprint, lod0.zmax)
        split_building = building.copy()
        split_building.add_geometry(building_surface, GeometryType.LOD0)
        split_building.calculate_bounds()
        split_buildings.append(split_building)

    return split_buildings


def clean_building_geometry(
    building: Building, lod=GeometryType.LOD2, tol=1e-2
) -> Building:
    """
    Clean building geometry by removing degenerate elements.

    Parameters
    ----------
    building : Building
        Building to clean.
    lod : GeometryType, default GeometryType.LOD2
        Level of detail to clean.
    tol : float, default 1e-2
        Tolerance for cleaning operations.

    Returns
    -------
    Building
        Building with cleaned geometry.
    """
    building_geom = building.geometry.get(lod, None)
    if building_geom is None:
        return building
    if isinstance(building_geom, MultiSurface):
        cleaned_geom = clean_multisurface(building_geom, tol)
    elif isinstance(building_geom, Surface):
        cleaned_geom = clean_surface(building_geom, tol)
    else:
        warning(f"Unsupported geometry type: {type(building_geom)}")
        return building
    building.add_geometry(cleaned_geom, lod)
    return building
