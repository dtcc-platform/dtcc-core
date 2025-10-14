import dtcc
from dtcc_core.builder import _dtcc_builder

from dtcc_core.model import (
    Mesh,
    VolumeMesh,
    Building,
    Terrain,
    City,
    Surface,
    MultiSurface,
    GeometryType,
    PointCloud,
)

from dtcc_core.builder.model_conversion import (
    create_builder_polygon,
    create_builder_surface,
    create_builder_multisurface,
    builder_mesh_to_mesh,
    builder_volume_mesh_to_volume_mesh,
    mesh_to_builder_mesh,
    create_builder_city,
    raster_to_builder_gridfield,
)





from dtcc_core.builder.meshing.convert import mesh_to_raster

from dtcc.logging import info, warning

def build_mesh(
    city: City,
    lod: GeometryType = None,
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
    if lod is None:
        lods = [GeometryType.LOD2, GeometryType.LOD1, GeometryType.LOD0]
        for test_lod in lods:
            if buildings and buildings[0].get_footprint(test_lod) is not None:
                lod = test_lod
                info(f"Using LOD {lod.name} for building footprints")
                break

    if merge_buildings:
        info(f"Merging {len(buildings)} buildings...")
        merged_buildings = dtcc.merge_building_footprints(
            buildings, lod, min_area=min_building_area
        )
        simplifed_footprints = dtcc.simplify_building_footprints(
            merged_buildings, min_building_detail / 2, lod=GeometryType.LOD0
        )
        clearance_fix = dtcc.fix_building_footprint_clearance(
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
        warning("City terrain has no data. Attempting to build terrain")
        city.build_terrain()
        terrain_raster = city.terrain.raster
        terrain_mesh = city.terrain.mesh
        if terrain_raster is None and terrain_mesh is None:
            raise ValueError("Failed to build city terrain")
    if terrain_raster is None and terrain_mesh is not None:
        terrain_raster = mesh_to_raster(terrain_mesh, cell_size=max_mesh_size)
    # builder_dem = raster_to_builder_gridfield(terrain_raster)

    xmin, ymin, xmax, ymax = terrain_raster.bounds.tuple
    builder_polygons = [
        create_builder_polygon(p.to_polygon()) for p in building_footprints if p is not None
    ]
    builder_mesh = _dtcc_builder.build_ground_mesh_spade(
        builder_polygons,
        subdomain_resolution,
        xmin,
        ymin,
        xmax,
        ymax,
        max_mesh_size,
        min_mesh_angle,
        sort_triangles,
    )
    
    result_mesh = builder_mesh.from_cpp()
    
    return result_mesh
# Poseidon (57.6971779, 11.9795910)
x0 = 319995.962899
y0 = 6399009.716755
L = 500.0

bounds = dtcc.Bounds(x0 - 0.5 * L, y0 - 0.5 * L, x0 + 0.5 * L, y0 + 0.5 * L)

# Download pointcloud and building footprints
pointcloud = dtcc.download_pointcloud(bounds=bounds)
buildings = dtcc.download_footprints(bounds=bounds)

# Remove global outliers
pointcloud = pointcloud.remove_global_outliers(3.0)

# Build terrain raster
raster = dtcc.build_terrain_raster(pointcloud, cell_size=2, radius=3, ground_only=True)

# Extract roof points and compute building heights
buildings = dtcc.extract_roof_points(buildings, pointcloud)
buildings = dtcc.compute_building_heights(buildings, raster, overwrite=True)

# Create city and add geometries
city = dtcc.City()
city.add_terrain(raster)
city.add_buildings(buildings, remove_outside_terrain=True)

# Build surface mesh
mesh = build_mesh(city=city)
# Offset to origin
# mesh.offset_to_origin()

# Save mesh to file
mesh.save("spade_surface_mesh.vtu")

