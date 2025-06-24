from dtcc_core.builder.register import register_model_method

from dtcc_core.model import Terrain, City, Building, Bounds, PointCloud


@register_model_method
def build_terrain(
    city: City,
    pc: PointCloud = None,
    cell_size: float = 2.0,
    build_mesh=True,
    max_triangle_size=5.0,
    smoothing=3,
) -> City:
    """
    Build terrain for a city using a point cloud.

    Args:
        city (City): The city object to build terrain for.
        pc (PointCloud): The point cloud to use for building the terrain.
        cell_size (float): The size of the cells in the raster (default is 2).

    Returns:
        City: The city object with the terrain added.
    """
    from dtcc_core.builder import build_terrain_raster, build_terrain_mesh

    if pc is not None and not isinstance(pc, PointCloud):
        raise ValueError("pc must be a PointCloud object")

    if pc is None:
        pc = city.pointcloud
        if pc is None:
            raise ValueError(
                "No point cloud provided and city has no point cloud geometry\n"
            )

    if len(pc.points) == 0:
        raise ValueError("Point cloud has no points")

    raster = build_terrain_raster(pc, cell_size=cell_size, ground_only=True)

    terrain = Terrain()
    terrain.add_raster(raster)
    if build_mesh:
        mesh = build_terrain_mesh(
            raster,
            max_mesh_size=max_triangle_size,
            smoothing=smoothing,
        )
        terrain.add_mesh(mesh)

    city.add_terrain(terrain)
    return city
