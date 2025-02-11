from ...model import Building, GeometryType, MultiSurface, Surface, GeometryType
from ..polygons.polygons import (
    polygon_merger,
    simplify_polygon,
    remove_slivers,
    fix_clearance,
    split_polygon_sides,
)

from ..polygons.surface import clean_multisurface, clean_surface

from ..register import register_model_method
from shapely.geometry import Polygon
from ..logging import debug, info, warning, error

from typing import List, Tuple, Union
from statistics import mean
import numpy as np


@register_model_method
def get_footprint(building: Building, geom_type: GeometryType = None) -> Surface:
    geom = None
    if geom_type is not None:
        geom = building.flatten_geometry(geom_type)
    else:
        lod_levels = [
            GeometryType.LOD0,
            GeometryType.LOD1,
            GeometryType.LOD2,
            GeometryType.LOD3,
        ]
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
        geom.to_polygon()
    # print(f"get footprint has {len(footprint.exterior.coords)} vertices")
    s = Surface()
    s.from_polygon(footprint, height)
    return s


def merge_building_footprints(
    buildings: List[Building],
    lod: GeometryType = GeometryType.LOD0,
    max_distance: float = 0.5,
    min_area=10,
) -> List[Building]:
    if len(buildings) <= 1:
        return buildings

    buildings_geom = [building.flatten_geometry(lod) for building in buildings]
    # print(buildings_geom)
    buildings_geom = [geom for geom in buildings_geom if geom is not None]
    building_heights = [geom.zmax for geom in buildings_geom]
    footprints = [geom.to_polygon(max_distance / 4) for geom in buildings_geom]
    merged_footprint, merged_indices = polygon_merger(
        footprints, max_distance, min_area=min_area
    )
    merged_buildings = []
    for idx, footprint in enumerate(merged_footprint):
        if footprint.geom_type == "MultiPolygon":
            ValueError("Merged footprint is a MultiPolygon")
        footprint = footprint.simplify(1e-2, True)
        if footprint.geom_type == "MultiPolygon":
            ValueError("simplified footprint is a MultiPolygon")
        footprint = remove_slivers(footprint, max_distance / 2)
        if footprint.geom_type == "MultiPolygon":
            ValueError("de-slivered footprint is a MultiPolygon")

        if footprint.is_empty or footprint.area < min_area:
            warning(f"Empty or too small footprint: {footprint.area}")
            continue
        indices = merged_indices[idx]

        original_buildings = [buildings[i] for i in indices]
        original_footprints = [footprints[i] for i in indices]
        area_argsort = np.argsort([footprint.area for footprint in original_footprints])

        original_buildings = [original_buildings[i] for i in area_argsort]
        building_attributes = merge_building_attributes(original_buildings)

        height = sum([building_heights[i] * footprints[i].area for i in indices]) / sum(
            [footprints[i].area for i in indices]
        )

        building_surface = Surface()
        building_surface.from_polygon(footprint, height)

        building = Building()
        building.add_geometry(building_surface, GeometryType.LOD0)
        building.attributes = building_attributes
        merged_buildings.append(building)
    return merged_buildings


def merge_building_attributes(buildings: List[Building]) -> dict:
    attributes = {}
    for building in buildings:
        for k, v in building.attributes.items():
            if v:
                attributes[k] = v
    return attributes


def simplify_building_footprints(
    buildings: List[Building],
    tolerance: float = 0.5,
    lod: GeometryType = GeometryType.LOD0,
) -> List[Building]:
    simplified_buildings = []
    for building in buildings:
        lod0 = building.lod0
        if lod0 is None:
            continue
        footprint = lod0.to_polygon()
        if footprint is None or footprint.is_empty:
            continue
        footprint = footprint.simplify(tolerance, True)
        building_surface = Surface()
        building_surface.from_polygon(footprint, lod0.zmax)
        simplified_building = building.copy()
        simplified_building.add_geometry(building_surface, GeometryType.LOD0)
        simplified_building.calculate_bounds()
        simplified_buildings.append(simplified_building)
    return simplified_buildings


def fix_building_footprint_clearance(
    buildings: List[Building],
    clearance: float = 0.5,
    lod: GeometryType = GeometryType.LOD0,
) -> List[Building]:
    fixed_buildings = []
    for building in buildings:
        lod0 = building.lod0
        if lod0 is None:
            continue
        footprint = lod0.to_polygon()
        if footprint is None or footprint.is_empty:
            continue
        footprint = fix_clearance(footprint, clearance)
        building_surface = Surface()
        building_surface.from_polygon(footprint, lod0.zmax)
        fixed_building = building.copy()
        fixed_building.add_geometry(building_surface, GeometryType.LOD0)
        fixed_building.calculate_bounds()
        fixed_buildings.append(fixed_building)
    return fixed_buildings


def split_footprint_walls(
    buildings: List[Building], max_wall_length: Union[float, List[float]] = 10
) -> List[Building]:
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
