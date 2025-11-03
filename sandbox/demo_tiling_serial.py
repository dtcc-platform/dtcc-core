# Copyright (C) 2025 DTCC Contributors
# Licensed under the MIT License
# Demo: Serial tiled mesh building for large cities

import dtcc
import time
import psutil
import os
from pathlib import Path


def get_memory_usage():
    """Get current memory usage in MB."""
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / 1024 / 1024


def main():
    print("=" * 80)
    print("DEMO: SERIAL TILED MESH BUILDING")
    print("=" * 80)
    print()

    # City parameters (same for both serial and parallel)
    x0 = 319995.962899
    y0 = 6399009.716755
    L = 2000.0  # 2000m x 2000m area (larger than serial default)
    bounds = dtcc.Bounds(x0 - 0.5 * L, y0 - 0.5 * L, x0 + 0.5 * L, y0 + 0.5 * L)

    print(f"City bounds: {bounds.bndstr}")
    print(f"Area: {bounds.width:.0f}m x {bounds.height:.0f}m = {bounds.area/1e6:.2f} km²")
    print()

    # Download data
    print("Step 1: Downloading data...")
    print("-" * 80)
    mem_start = get_memory_usage()
    time_start = time.perf_counter()

    pointcloud = dtcc.download_pointcloud(bounds=bounds)
    buildings = dtcc.download_footprints(bounds=bounds)

    time_download = time.perf_counter() - time_start
    mem_after_download = get_memory_usage()

    print(f"  Pointcloud: {len(pointcloud.points)} points")
    print(f"  Buildings: {len(buildings)} footprints")
    print(f"  Download time: {time_download:.2f}s")
    print(f"  Memory increase: {mem_after_download - mem_start:.1f}MB")
    print()

    # Remove outliers
    print("Step 2: Preprocessing data...")
    print("-" * 80)
    time_start = time.perf_counter()

    pointcloud = pointcloud.remove_global_outliers(3.0)

    time_preprocess = time.perf_counter() - time_start
    print(f"  Removed outliers, {len(pointcloud.points)} points remain")
    print(f"  Preprocessing time: {time_preprocess:.2f}s")
    print()

    # Build terrain
    print("Step 3: Building terrain...")
    print("-" * 80)
    time_start = time.perf_counter()

    raster = dtcc.build_terrain_raster(pointcloud, cell_size=2, radius=3, ground_only=True)

    time_terrain = time.perf_counter() - time_start
    print(f"  Terrain raster: {raster.data.shape}")
    print(f"  Terrain time: {time_terrain:.2f}s")
    print()

    # Extract roof points and building heights
    print("Step 4: Computing building heights...")
    print("-" * 80)
    time_start = time.perf_counter()

    buildings = dtcc.extract_roof_points(buildings, pointcloud)
    buildings = dtcc.compute_building_heights(buildings, raster, overwrite=True)

    time_heights = time.perf_counter() - time_start
    print(f"  Buildings with heights: {len([b for b in buildings if b.height > 0])}")
    print(f"  Height computation time: {time_heights:.2f}s")
    print()

    # Create city
    print("Step 5: Creating city object...")
    print("-" * 80)
    time_start = time.perf_counter()

    city = dtcc.City()
    city.add_terrain(raster)
    city.add_buildings(buildings, remove_outside_terrain=True)

    time_city = time.perf_counter() - time_start
    print(f"  City buildings: {len(city.buildings)}")
    print(f"  City creation time: {time_city:.2f}s")
    print()

    # Build mesh WITH SERIAL TILING
    print("Step 6: Building surface mesh (SERIAL TILING - n_workers=1)...")
    print("-" * 80)
    print("Parameters:")
    print(f"  Tile size: 1000.0m x 1000.0m")
    print(f"  Tile overlap: 50.0m")
    print(f"  Workers: 1 (SERIAL PROCESSING)")
    print(f"  LOD: LOD0")
    print(f"  Building merge: True")
    print(f"  Triangle size: 5.0m")
    print(f"  Max mesh size: 10.0m")
    print(f"  Min mesh angle: 25.0°")
    print()

    mem_start_mesh = get_memory_usage()
    time_mesh_start = time.perf_counter()

    mesh = dtcc.build_city_mesh(
        city,
        lod=dtcc.GeometryType.LOD0,
        use_tiling=True,
        tile_size=1000.0,
        tile_overlap=50.0,
        n_workers=1,  # SERIAL PROCESSING
        merge_buildings=True,
        building_mesh_triangle_size=5.0,
        max_mesh_size=10.0,
        min_mesh_angle=25.0,
        merge_meshes=True,
        sort_triangles=False,
    )

    time_mesh = time.perf_counter() - time_mesh_start
    mem_after_mesh = get_memory_usage()
    peak_memory_mesh = mem_after_mesh - mem_start

    print(f"  Mesh vertices: {len(mesh.vertices)}")
    print(f"  Mesh faces: {len(mesh.faces)}")
    print(f"  Mesh time: {time_mesh:.2f}s")
    print(f"  Peak memory usage: {peak_memory_mesh:.1f}MB")
    print(f"  Memory per vertex: {peak_memory_mesh / len(mesh.vertices) * 1000:.2f}KB" if len(mesh.vertices) > 0 else "")
    print()

    # Save mesh
    print("Step 7: Saving mesh...")
    print("-" * 80)
    output_path = "demo_mesh_serial_tiling.vtu"
    time_start = time.perf_counter()

    mesh.save(output_path)

    time_save = time.perf_counter() - time_start
    file_size = os.path.getsize(output_path) / 1024 / 1024

    print(f"  Output file: {output_path}")
    print(f"  File size: {file_size:.1f}MB")
    print(f"  Save time: {time_save:.2f}s")
    print()

    # Summary
    print("=" * 80)
    print("SUMMARY - SERIAL TILING")
    print("=" * 80)
    total_time = time_download + time_preprocess + time_terrain + time_heights + time_city + time_mesh + time_save
    print(f"Total processing time: {total_time:.2f}s")
    print(f"  - Download: {time_download:.2f}s ({time_download/total_time*100:.1f}%)")
    print(f"  - Preprocessing: {time_preprocess:.2f}s ({time_preprocess/total_time*100:.1f}%)")
    print(f"  - Terrain: {time_terrain:.2f}s ({time_terrain/total_time*100:.1f}%)")
    print(f"  - Heights: {time_heights:.2f}s ({time_heights/total_time*100:.1f}%)")
    print(f"  - City: {time_city:.2f}s ({time_city/total_time*100:.1f}%)")
    print(f"  - Mesh: {time_mesh:.2f}s ({time_mesh/total_time*100:.1f}%)")
    print(f"  - Save: {time_save:.2f}s ({time_save/total_time*100:.1f}%)")
    print()
    print(f"Peak memory usage: {mem_after_mesh:.1f}MB")
    print(f"Output mesh: {len(mesh.vertices)} vertices, {len(mesh.faces)} faces")
    print()
    print("Note: Run demo_tiling_parallel.py to compare with parallel processing")
    print("=" * 80)


if __name__ == "__main__":
    main()
