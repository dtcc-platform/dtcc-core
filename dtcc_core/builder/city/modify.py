from ..polygons.polygons import (
    polygon_merger,
    simplify_polygon,
    remove_slivers,
    fix_clearance,
)

from ..building.modify import (
    _condition_buildings_with_shared_cleaner,
    clean_building_geometry,
    fix_building_footprint_clearance,
    merge_building_footprints,
    simplify_building_footprints,
)
from ..cleaning import ConditioningOptions

import numpy as np

from ..register import register_model_method
from dtcc_core.model import City, Bounds, Terrain, Building, GeometryType, Raster
from statistics import mean
import shapely
import dataclasses
from copy import deepcopy
from collections import defaultdict
from shapely.geometry import MultiPolygon, Polygon, JOIN_STYLE, CAP_STYLE
from ..logging import info, warning, error


@register_model_method
def simplify_buildings(city: City, tolerance=0.1) -> City:
    """
    Simplify the footprint of buildings in a `City` object.

    Args:
        city (City): The `City` object to simplify the buildings of.
        tolerance (float): The tolerance for simplification (default 0.1).

    Returns:
        City: A new `City` object with the simplified buildings.
    """
    simplified_city = deepcopy(city)
    simplified_city.buildings = []
    for b in city.buildings:
        b = dataclasses.replace(b)
        b.footprint = simplify_polygon(b.footprint, tolerance)
        simplified_city.buildings.append(b)
    return simplified_city


@register_model_method
def remove_small_buildings(city: City, min_area=10) -> City:
    """
    Remove small buildings from a `City` object.

    Args:
        city (City): The `City` object to remove small buildings from.
        min_area (float): The minimum area in square meters for a building to be kept (default 10).

    Returns:
        City: A new `City` object with the small buildings removed.
    """
    filtered_city = deepcopy(city)
    filtered_city.buildings = []
    for b in city.buildings:
        if b.footprint.area > min_area:
            filtered_city.buildings.append(b)
    return filtered_city


@register_model_method
def merge_buildings(
    city: City,
    max_distance=0.15,
    min_area=10,
    simplify=True,
    properties_merge_strategy="list",
    height_merge_strategy="mean",
) -> City:
    """
    Merge buildings that are close together.

    Parameters
    ----------
    max_distance : float, optional
        The maximum distance in meters between buildings to consider them close enough to merge (default 0.15).
    min_area : float, optional
        The minimum area in square meters for a building to be kept (default 10).
    simplify : bool, optional
        Whether to simplify the merged buildings (default True).
    properties_merge_strategy : str, optional
        The strategy for merging properties. Options are 'list' and 'sample'. 'list' will create a list of all properties for the merged building. 'sample' will pick a property value from a random building (default "list").
    height_merge_strategy : str, optional
        The strategy for merging heights. Options are 'mean', 'area_weighted' and 'max' .
        'mean' will take the mean height of the merged buildings.
        'area_weighted' will take the area weighted mean height of the merged buildings.
        'max' will take the maximum height of the merged buildings (default "mean").

    Returns
    -------
    City
        A new `City` object with the merged buildings.
    """
    warning(
        "merge_buildings() now delegates to the shared footprint conditioner; "
        "legacy property and height merge strategies are no longer applied separately."
    )
    merged_city = deepcopy(city)
    merged = _condition_buildings_with_shared_cleaner(
        city.buildings,
        lod=GeometryType.LOD0,
        options=ConditioningOptions(
            precision_grid=None,
            min_feature_size=max(max_distance / 2.0, 1e-3) if simplify else 0.0,
            merge_distance=max_distance,
            min_area=min_area,
            min_hole_area=0.0,
        ),
        operation_name="merge_buildings",
        return_index_map=False,
    )
    merged_city.replace_buildings(merged)
    return merged_city


@register_model_method
def fix_building_clearance(
    city: City, target_clearance: float, min_angle: float, accepted_tol_fraction=0.9
) -> City:
    """
    Fix the clearance of the footprints in the building models. After running
    each building should have a minimum_clearance of `tol` meters and a minimum
    angle between each edge of `min_angle` degrees.
    """
    warning(
        "fix_building_clearance() now delegates to the shared footprint conditioner; "
        "min_angle and accepted_tol_fraction are ignored."
    )
    fixed_city = deepcopy(city)
    fixed_city.replace_buildings(
        fix_building_footprint_clearance(
            city.buildings,
            clearance=target_clearance,
            lod=GeometryType.LOD0,
            return_index_map=False,
        )
    )
    return fixed_city


def clean_building_surfaces(city: City, lod: GeometryType, tol=1e-2) -> City:
    """
    Clean building surfaces in a city by removing degenerate geometry elements.

    This function applies geometry cleaning operations to all buildings in the city
    to remove invalid or degenerate surfaces, vertices, and other problematic elements.

    Parameters
    ----------
    city : City
        The city object containing buildings to clean.
    lod : GeometryType
        The level of detail geometry to clean (e.g., LOD0, LOD1, LOD2).
    tol : float, default 1e-2
        Tolerance value for cleaning operations in meters.

    Returns
    -------
    City
        The city object with cleaned building surfaces.
    """
    for building in city.buildings:
        building = clean_building_geometry(building, lod, tol)
    return city


@register_model_method
def add_flat_terrain(city: City, ground_level: float = None, buffer: float = 0) -> City:
    """
    Add flat terrain to a city model.

    ground_level (float): The ground level of the terrain (default None). If None, the ground level will be set to the
    zmin of the city bounds
    buffer (float): The buffer to add to the bounds (default 0).
    """

    terrain = Terrain()
    raster = Raster()
    if ground_level is None:
        ground_level = city.bounds.zmin
    raster.data = np.ones((1, 1)) * ground_level
    if buffer != 0:
        city.bounds.buffer(buffer)
    raster.set_bounds(city.bounds)
    terrain.add_geometry(raster, GeometryType.RASTER)
    city.add_terrain(terrain)
    return city


# @register_model_method
# def calculate_bounds(city: City, buffer: float = 0) -> City:
#     """
#     Calculate the bounds of a `City` object.

#     Args:
#         city (City): The `City` object to calculate the bounds of.
#         buffer (float): The buffer to add to the bounds (default 0).

#     Returns:
#         City: A new `City` object with the bounds calculated.
#     """
#     footprints = [b.footprint for b in city.buildings]
#     bounds = MultiPolygon(footprints).bounds
#     city.bounds = Bounds(bounds[0], bounds[1], bounds[2], bounds[3])
#     if buffer != 0:
#         city.bounds.buffer(buffer)
#     return city
