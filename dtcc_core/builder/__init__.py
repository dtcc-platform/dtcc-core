from . import model_conversion
from . import builders
from . import parameters
from . import city_methods
from . import meshing
from . import pointcloud
from . import raster
from . import geometry
from . import city
from . import building
from . import roadnetwork

# from .builders import (
#     build,
#     build_city,
#     calculate_bounds,
#     build_building_meshes,
#     build_city_surface_mesh,
#     build_volume_mesh,
# )

from .geometry_builders.terrain import (
    build_terrain_mesh,
    build_terrain_raster,
    flat_terrain,
)
from .geometry_builders.buildings import (
    extract_roof_points,
    compute_building_heights,
    build_lod1_buildings,
    building_heights_from_pointcloud,

)


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

from .geometry_builders.meshes import build_surface_mesh

# # Add model extensions
# from ..model import City, PointCloud

# City.add_methods(city_methods.compute_building_points, "compute_building_points")
# City.add_methods(city_methods.compute_building_heights, "compute_building_heights")


__all__ = [
    "extract_roof_points",
    "compute_building_heights",
    "build_lod1_buildings",
    "build_surface_mesh",
    "build_terrain_mesh",
    "build_terrain_raster",
    "flat_terrain",
    "merge_building_footprints",
    "simplify_building_footprints",
    "fix_building_footprint_clearance",
    "split_footprint_walls",
]
