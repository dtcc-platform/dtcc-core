from ...model import Building, GeometryType, MultiSurface, Surface
from ..polygons.polygons import (
    polygon_merger,
    simplify_polygon,
    remove_slivers,
    split_polygon_sides,
)

from polyforge import (
    merge_close_polygons,
    fix_clearance,
    simplify_rdp,
    simplify_vwp,
    robust_fix_geometry,
)
from polyforge import MergeStrategy, GeometryConstraints

from ..polygons.surface import clean_multisurface, clean_surface

from ..register import register_model_method
from shapely.geometry import Polygon, MultiPolygon
from shapely.ops import unary_union
from ..logging import debug, info, warning, error

from typing import List, Tuple, Union
import warnings


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
    if len(buildings) <= 1:
        if return_index_map:
            return buildings, [[i] for i in range(len(buildings))]
        return buildings

    source_indices: List[int] = []
    footprints: List[Polygon] = []
    building_heights: List[float] = []

    for idx, building in enumerate(buildings):
        flattened_geom = building.get_footprint(lod)
        if flattened_geom is None:
            warning(f"Building {building.id} has no geometry at LOD {lod}. Skipping.")
            continue
        footprint = flattened_geom.to_polygon()
        if footprint is None or footprint.is_empty:
            warning(f"Building {building.id} produced an empty footprint. Skipping.")
            continue
        source_indices.append(idx)
        building_heights.append(flattened_geom.zmax)
        footprints.append(footprint)

    merged_footprint, merged_indices = merge_close_polygons(
        footprints,
        max_distance,
        merge_strategy=MergeStrategy.SELECTIVE_BUFFER,
        preserve_holes=True,
        insert_vertices=True,
        return_mapping=True,
    )

    merged_buildings: List[Building] = []
    merged_indices_global: List[List[int]] = []

    for idx, footprint in enumerate(merged_footprint):
        if footprint.geom_type == "MultiPolygon":
            ValueError("Merged footprint is a MultiPolygon")
        if footprint.is_empty or footprint.area < min_area:
            warning(f"Empty or too small footprint: {footprint.area}")
            continue

        local_indices = merged_indices[idx]
        global_indices = [source_indices[i] for i in local_indices]

        num = sum(building_heights[i] * footprints[i].area for i in local_indices)
        den = sum(footprints[i].area for i in local_indices)
        height: float = num / den if den > 0 else 0.0

        building_surface = Surface()
        building_surface.from_polygon(footprint, height)

        building = Building()
        building.add_geometry(building_surface, GeometryType.LOD0)

        original_buildings = [buildings[i] for i in global_indices]
        building.attributes = merge_building_attributes(original_buildings)

        merged_buildings.append(building)
        merged_indices_global.append(global_indices)

    if return_index_map:
        return merged_buildings, merged_indices_global
    return merged_buildings


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

    simplified_buildings: List[Building] = []
    index_map: List[List[int]] = []
    if method not in ["vwp", "rdp"]:
        warning(
            f"Unknown polygon simplification method: {method}. Falling back to 'vwp"
        )
        method = "vwp"

    for idx, building in enumerate(buildings):
        lod_geom = building.flatten_geometry(lod)
        if lod_geom is None:
            continue
        footprint = lod_geom.to_polygon()
        if footprint is None or footprint.is_empty:
            continue

        if method == "vwp":
            footprint = simplify_vwp(footprint, tolerance)
        elif method == "rdp":
            footprint = simplify_vwp(footprint, tolerance)

        building_surface = Surface()
        building_surface.from_polygon(footprint, lod_geom.zmax)
        simplified_building = building.copy()
        simplified_building.add_geometry(building_surface, GeometryType.LOD0)
        simplified_building.calculate_bounds()
        simplified_buildings.append(simplified_building)
        if return_index_map:
            index_map.append([idx])

    if return_index_map:
        return simplified_buildings, index_map
    return simplified_buildings


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

    fixed_buildings = []
    index_map: List[List[int]] = []
    constraints = GeometryConstraints(
        must_be_valid=True,
        min_clearance=clearance,
        min_hole_area=smallest_hole_area,
        allow_multipolygon=False,
    )

    for idx, building in enumerate(buildings):
        lod0 = building.lod0
        if lod0 is None:
            continue
        footprint = lod0.to_polygon()
        if footprint is None or footprint.is_empty:
            continue
        with warnings.catch_warnings(record=True) as caught_warnings:
            warnings.simplefilter("always", UserWarning)
            fixed_footprint, _ = robust_fix_geometry(footprint, constraints=constraints)
        for caught_warning in caught_warnings:
            if issubclass(caught_warning.category, UserWarning):
                warning(
                    f"Building {building.id}: {str(caught_warning.message).strip()}"
                )
        building_surface = Surface()
        building_surface.from_polygon(fixed_footprint, lod0.zmax)
        fixed_building = building.copy()
        fixed_building.add_geometry(building_surface, GeometryType.LOD0)
        fixed_building.calculate_bounds()
        fixed_buildings.append(fixed_building)
        if return_index_map:
            index_map.append([idx])

    if return_index_map:
        return fixed_buildings, index_map
    return fixed_buildings


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
    fixed_buildings: List[Building] = []
    index_map: List[List[int]] = []

    for idx, building in enumerate(buildings):
        lod_geom = building.flatten_geometry(lod)
        if lod_geom is None:
            continue
        footprint = lod_geom.to_polygon()
        if footprint is None or footprint.is_empty:
            continue
        # clean_surface
        footprint = fix_clearance(footprint, clearance)
        building_surface = Surface()
        building_surface.from_polygon(footprint, lod_geom.zmax)
        fixed_building = building.copy()
        fixed_building.add_geometry(building_surface, GeometryType.LOD0)
        fixed_building.calculate_bounds()
        fixed_buildings.append(fixed_building)
        if return_index_map:
            index_map.append([idx])

    if return_index_map:
        return fixed_buildings, index_map
    return fixed_buildings


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
