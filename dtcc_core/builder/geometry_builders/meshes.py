from ...model import (
    Mesh,
    VolumeMesh,
    Building,
    Terrain,
    City,
    Surface,
    MultiSurface,
    GeometryType,
)

from ..model_conversion import (
    create_builder_polygon,
    create_builder_surface,
    create_builder_multisurface,
    builder_mesh_to_mesh,
    mesh_to_builder_mesh,
    create_builder_city,
    raster_to_builder_gridfield,
)

from .. import _dtcc_builder

from ..polygons.polygons import (
    polygon_merger,
    simplify_polygon,
    remove_slivers,
    fix_clearance,
)

from ..building.modify import (
    merge_building_footprints,
    simplify_building_footprints,
    fix_building_footprint_clearance,
)

from ..meshing.convert import mesh_to_raster

from ..logging import debug, info, warning, error


def build_city_mesh(
    city: City,
    lod: GeometryType = GeometryType.LOD1,
    min_building_detail: float = 0.5,
    min_building_area: float = 15.0,
    merge_buildings: bool = True,
    merge_tolerance: float = 0.5,
    building_mesh_triangle_size: float = 5.0,
    max_mesh_size: float = 10.0,
    min_mesh_angle: float = 25.0,
    merge_meshes: bool = True,
    smoothing: int = 0,
    sort_triangles: bool = False,
) -> Mesh:
    """
    Build a mesh from the surfaces of the buildings in the city.

    Parameters
    ----------
    `city` : model.City
        The city to build the mesh from.
    `min_building_detail` : float, optional
        The minimum detail of the buildin to resolve, by default 0.5.
    `min_building_area` : float, optional
        The smallest building to include, by default 15.0.
    `merge_buildings` : bool, optional
        merge building footprints, by default True.
    `max_mesh_size` : float, optional
        The maximum size of the mesh, by default 1.0.
    `min_mesh_angle` : float, optional
        The minimum angle of the mesh, by default 30.0.
    `merge_meshes` : bool, optional
        Whether to merge the meshes to a single mesh, by default True.

    `smoothing` : float, optional
        The smoothing of the mesh, by default 0.0.

    Returns
    -------
    `model.Mesh`
    """
    buildings = city.buildings
    if merge_buildings:
        info(f"Merging {len(buildings)} buildings...")
        merged_buildings = merge_building_footprints(
            buildings, lod, min_area=min_building_area
        )
        simplifed_footprints = simplify_building_footprints(
            merged_buildings, min_building_detail / 2, lod=GeometryType.LOD0
        )
        clearance_fix = fix_building_footprint_clearance(
            simplifed_footprints, min_building_detail
        )
        building_footprints = [
            b.get_footprint(GeometryType.LOD0) for b in clearance_fix
        ]
        info(f"After merging: {len(building_footprints)} buildings.")
    else:
        building_footprints = [b.get_footprint(lod) for b in buildings]

    subdomain_resolution = [building_mesh_triangle_size] * len(building_footprints)

    terrain = city.terrain
    if terrain is None:
        raise ValueError("City has no terrain data. Please compute terrain first.")
    terrain_raster = terrain.raster
    terrain_mesh = terrain.mesh
    if terrain_raster is None and terrain_mesh is None:
        raise ValueError("City terrain has no data. Please compute terrain first.")
    if terrain_raster is None and terrain_mesh is not None:
        terrain_raster = mesh_to_raster(terrain_mesh, cell_size=max_mesh_size)
    builder_dem = raster_to_builder_gridfield(terrain_raster)

    builder_surfaces = [
        create_builder_surface(p) for p in building_footprints if p is not None
    ]
    builder_mesh = _dtcc_builder.build_city_surface_mesh(
        builder_surfaces,
        subdomain_resolution,
        builder_dem,
        max_mesh_size,
        min_mesh_angle,
        smoothing,
        merge_meshes,
        sort_triangles,
    )
    if merge_meshes:
        result_mesh = builder_mesh_to_mesh(builder_mesh[0])
    else:
        result_mesh = [builder_mesh_to_mesh(bm) for bm in builder_mesh]
    return result_mesh
