from __future__ import annotations

import importlib

from . import model_conversion
from . import meshing
from . import geometry
from . import city
from . import building
from . import polygons
from . import roadnetwork

from .building.modify import (
    merge_building_footprints,
    simplify_building_footprints,
    fix_building_footprint_clearance,
    split_footprint_walls,
)

from .city.modify import (
    merge_buildings,
    fix_building_clearance,
    clean_building_surfaces,
)

from .register import register_model_method

from .geometry_builders.meshes import (
    build_conditioned_footprints,
    build_city_surface_mesh,
    build_city_flat_mesh,
    build_city_volume_mesh,
)
from .meshing import (
    available_2d_meshers,
    get_default_2d_mesher,
    set_default_2d_mesher,
)

_LAZY_IMPORTS = {
    "pointcloud": "dtcc_core.builder.pointcloud",
    "raster": "dtcc_core.builder.raster",
    "trees": "dtcc_core.builder.trees",
    "build_terrain_surface_mesh": "dtcc_core.builder.geometry_builders.terrain",
    "build_terrain_raster": "dtcc_core.builder.geometry_builders.terrain",
    "flat_terrain": "dtcc_core.builder.geometry_builders.terrain",
    "flatten_terrain_raster": "dtcc_core.builder.geometry_builders.terrain",
    "adaptive_terrain_mesh": "dtcc_core.builder.geometry_builders.terrain",
    "extract_roof_points": "dtcc_core.builder.geometry_builders.buildings",
    "compute_building_heights": "dtcc_core.builder.geometry_builders.buildings",
    "set_building_heights_from_attribute": "dtcc_core.builder.geometry_builders.buildings",
    "build_lod1_buildings": "dtcc_core.builder.geometry_builders.buildings",
    "build_lod2_buildings": "dtcc_core.builder.geometry_builders.lod2",
    "extrude_building": "dtcc_core.builder.geometry_builders.buildings",
    "building_heights_from_pointcloud": "dtcc_core.builder.geometry_builders.buildings",
    "tree_raster_from_pointcloud": "dtcc_core.builder.trees.create",
    "find_tree_tops": "dtcc_core.builder.trees.create",
    "trees_from_pointcloud": "dtcc_core.builder.trees.create",
    "tree_crown_polygons": "dtcc_core.builder.trees.create",
}


def __getattr__(name: str):
    if name in _LAZY_IMPORTS:
        module = importlib.import_module(_LAZY_IMPORTS[name])
        if name in {"pointcloud", "raster", "trees"}:
            value = module
        else:
            value = getattr(module, name)
        globals()[name] = value
        return value
    raise AttributeError(f"module '{__name__}' has no attribute '{name}'")


from .meshing.boundary_conformance import conform_boundary

__all__ = [
    "extract_roof_points",
    "compute_building_heights",
    "build_lod1_buildings",
    "build_lod2_buildings",
    "build_conditioned_footprints",
    "build_city_surface_mesh",
    "build_city_flat_mesh",
    "build_terrain_surface_mesh",
    "build_terrain_raster",
    "flat_terrain",
    "flatten_terrain_raster",
    "merge_building_footprints",
    "simplify_building_footprints",
    "fix_building_footprint_clearance",
    "split_footprint_walls",
    "merge_buildings",
    "fix_building_clearance",
    "clean_building_surfaces",
    "building_heights_from_pointcloud",
    "set_building_heights_from_attribute",
    "tree_raster_from_pointcloud",
    "find_tree_tops",
    "trees_from_pointcloud",
    "tree_crown_polygons",
    "build_city_volume_mesh",
    "adaptive_terrain_mesh",
    "available_2d_meshers",
    "get_default_2d_mesher",
    "set_default_2d_mesher",
    "conform_boundary",
]
