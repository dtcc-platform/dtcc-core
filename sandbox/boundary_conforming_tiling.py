# Copyright(C) 2026 DTCC
# Licensed under the MIT License

"""Demo: boundary-conforming tiled city mesh.

Builds a city surface mesh, demonstrates boundary conformance
by treating some meshes as "already built" (old meshes) and conforming a new
tile's boundary to match them before final tiling and export.
"""
import dtcc
from pathlib import Path
from time import time

from dtcc_core.model import Bounds, City
from dtcc_core.builder.meshing import conform_boundary
from dtcc_core.io.meshes import load_mesh, save

# ── Parameters ───────────────────────────────────────────────────────────
tile_size = 250.0

outdir = Path("conforming_tiling_output")


def build_surface_mesh(bounds: dtcc.Bounds, smoothing_iterations: int = 0) -> dtcc.Mesh:
    city = dtcc.City()
    city.bounds = bounds
    city.download_pointcloud(bounds=bounds, filter_on_z_bounds=True)
    city.download_footprints(bounds=bounds)
    city.building_heights_from_pointcloud()
    lod = [dtcc.GeometryType.LOD1 for _ in city.buildings]
    return city.build_surface_mesh(
        lod=lod, treat_lod0_as_holes=False, smoothing=smoothing_iterations
    )

if __name__ == "__main__":
    
    outdir.mkdir(exist_ok=True, parents=True)
    h = 1000.0

    # Tile 1 (new) — center tile
    bounds_tile_1 = Bounds(319891, 6399790, 319891 + h, 6399790 + h, 0, 200)

    # Tile 2 (old) — south tile
    bounds_tile_2 = Bounds(319891, 6399790 - h, 319891 + h, 6399790, 0, 200)

    print("Building tile 1 (new)...")
    tile_1_mesh = build_surface_mesh(bounds_tile_1, smoothing_iterations=0)

    print("Building tile 2 (old)...")
    tile_2_mesh = build_surface_mesh(bounds_tile_2, smoothing_iterations=5)

    new_tile = tile_1_mesh
    new_tile.save(outdir / "tile_original.vtk")
    old_meshes = [tile_2_mesh]
    for i, m in enumerate(old_meshes):
        m.save(outdir / f"old_tile_{i}.vtk")
    
    start = time()
    conformed_mesh = conform_boundary(new_mesh = new_tile, 
                                      old_meshes= old_meshes,
                                      smoothing_iterations=5,
                                      tolerance=10.0)
    print(f"Boundary conformance took {time() - start:.2f}s")

    conformed_mesh.save(outdir / "conformed_tile.vtk")

    conformed_mesh.offset_to_origin()
    conformed_mesh.offset([0, 0, 1])
    conformed_mesh.save(outdir / "conformed_tile_offset.vtk")


    conformed_tiles = conformed_mesh.tile(tile_size=tile_size)

    print_tiles = [t.create_printable_solid() for t in conformed_tiles]

    tiles_outdir = outdir / "tiles"
    tiles_outdir.mkdir(exist_ok=True, parents=True) 

    for i, t in enumerate(print_tiles):
        t.save(tiles_outdir / f"tiled_city_{i}.stl")