from dtcc_core.model import Mesh, Bounds
from typing import Tuple, Union, List

from tqdm import tqdm


def tile_surface_mesh(
    mesh: Mesh,
    tile_size: Union[float, Tuple[float, float]] = 100.0,
    progress: bool = False,
) -> list[Mesh]:
    """
    Tile a surface mesh into smaller meshes of a specified size.

    Args:
        mesh (Mesh): The input surface mesh to be tiled.
        tile_size (float|tuple): The size of each tile (in the same units as the mesh coordinates). if a single float
        is provided, it will be used for both width and height.
        bounds (Bounds): Optional bounds to limit the tiling area.

    Returns:
        list[Mesh]: A list of tiled mesh objects.
    """
    from dtcc_core.builder.meshing import SurfaceMeshClipper

    surface_bounds = mesh.bounds
    surface_tile_bounds = surface_bounds.tiles(tile_size)
    tiled_meshes = []
    clipper = SurfaceMeshClipper(mesh)

    for idx, tile in tqdm(
        enumerate(surface_tile_bounds),
        total=len(surface_tile_bounds),
        desc="Tiling mesh",
        disable=not progress,
    ):
        tiled_mesh = clipper.clip_to_bounds(tile)
        tiled_meshes.append(tiled_mesh)
    return tiled_meshes


def extrude_mesh_base(mesh: Mesh, extrude_to=0) -> Mesh:
    """
    Extrude a surface mesh downwards to create a closed mesh

    Args:
        mesh (Mesh): The input surface mesh to be extruded.
        z_base (float): The Z value to extrude the mesh down to.

    Returns:
        Mesh: The resulting volume mesh.
    """
