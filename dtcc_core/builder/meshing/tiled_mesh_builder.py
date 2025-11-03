# Copyright (C) 2025 DTCC Contributors
# Licensed under the MIT License
# Tiled mesh builder for processing large cities with limited memory

import numpy as np
from typing import List, Tuple, Optional, Union
from pathlib import Path
from multiprocessing import Pool
import psutil
from tqdm import tqdm

from ...model import City, Building, Mesh, Bounds, GeometryType
from ..logging import info, warning, error, debug
from .shared_memory_backend import SharedMeshStore, TileMetadata


def calculate_optimal_tile_size(available_memory_fraction: float = 0.3) -> float:
    """
    Calculate optimal tile size based on available system memory.

    Args:
        available_memory_fraction: Fraction of available memory to use (default 0.3 = 30%)

    Returns:
        Recommended tile size in coordinate units (approximation based on memory)
    """
    # Get available memory in bytes
    mem = psutil.virtual_memory()
    available_bytes = mem.available * available_memory_fraction

    # Rough estimate: 8 bytes per vertex (3 floats = 24 bytes)
    # Each vertex appears in ~6 triangles, so ~144 bytes per unique vertex
    # For a square tile of size S x S with element size ~10m:
    # Num vertices ~ (S/10)^2, so memory ~ (S/10)^2 * 144
    # S^2 / 100 * 144 = available_bytes
    # S^2 = available_bytes * 100 / 144
    # S = sqrt(available_bytes * 100 / 144)

    tile_size = np.sqrt(available_bytes * 100 / 144)
    return float(tile_size)


def _process_mesh_tile(args) -> Tuple[int, np.ndarray, np.ndarray, Optional[np.ndarray], Tuple[float, float, float, float]]:
    """
    Worker function to process a single mesh tile.

    Args:
        args: Tuple of (tile_id, tile_bounds, buildings, terrain_raster, mesh_params)

    Returns:
        Tuple of (tile_id, vertices, faces, markers, bounds)
    """
    from ..geometry_builders import meshes as mesh_module

    tile_id, tile_bounds, buildings, terrain_raster, mesh_params = args

    try:
        # Create a temporary city object for this tile
        temp_city = City()
        temp_city.add_buildings(buildings)

        # Add terrain if available
        if terrain_raster is not None:
            from ...model import Terrain
            terrain = Terrain()
            terrain.add_geometry(terrain_raster, GeometryType.RASTER)
            temp_city.add_terrain(terrain)

        # Build mesh for this tile
        mesh = mesh_module.build_city_mesh(
            city=temp_city,
            **mesh_params
        )

        if mesh is None or len(mesh.vertices) == 0:
            # Return empty mesh
            return (
                tile_id,
                np.array([], dtype=np.float64).reshape(0, 3),
                np.array([], dtype=np.int32).reshape(0, 3),
                None,
                tile_bounds
            )

        return (
            tile_id,
            mesh.vertices,
            mesh.faces,
            mesh.markers if hasattr(mesh, 'markers') else None,
            tile_bounds
        )
    except Exception as e:
        import traceback
        error(f"Error processing tile {tile_id}: {str(e)}")
        traceback.print_exc()
        # Return empty mesh on error
        return (
            tile_id,
            np.array([], dtype=np.float64).reshape(0, 3),
            np.array([], dtype=np.int32).reshape(0, 3),
            None,
            tile_bounds
        )


def build_city_mesh_tiled(
    city: City,
    lod: GeometryType = None,
    tile_size: float = 1000.0,
    overlap: float = 50.0,
    n_workers: Optional[int] = None,
    min_building_detail: float = 0.5,
    min_building_area: float = 15.0,
    merge_buildings: bool = True,
    merge_tolerance: float = 0.5,
    building_mesh_triangle_size: float = 5.0,
    max_mesh_size: float = 10.0,
    min_mesh_angle: float = 25.0,
    smoothing: int = 0,
    sort_triangles: bool = False,
    progress: bool = True,
    use_disk_store: bool = True,
    temp_dir: Optional[str] = None,
    keep_intermediate: bool = False,
    snap_distance: float = 0.01,
) -> Mesh:
    """
    Build a city surface mesh using spatial tiling for memory efficiency.

    This function divides the city into spatial tiles, processes each tile independently,
    and merges the results incrementally to avoid memory spikes. Uses disk-backed storage
    by default to handle arbitrarily large cities with limited RAM.

    Args:
        city: City object to mesh
        lod: Level of detail for building footprints
        tile_size: Size of each tile in coordinate units (default 1000.0)
        overlap: Overlap between tiles to handle boundary objects (default 50.0)
        n_workers: Number of worker processes (default: CPU count)
        min_building_detail: Minimum building detail threshold
        min_building_area: Minimum building area to include
        merge_buildings: Whether to merge building footprints
        merge_tolerance: Tolerance for building merging
        building_mesh_triangle_size: Triangle size for building meshes
        max_mesh_size: Maximum mesh element size
        min_mesh_angle: Minimum mesh angle
        smoothing: Smoothing iterations
        sort_triangles: Whether to sort triangles
        progress: Whether to show progress bar
        use_disk_store: Use disk-backed HDF5 storage for memory efficiency (default True)
        temp_dir: Directory for temporary HDF5 files (default: system temp)
        keep_intermediate: Keep intermediate HDF5 files after completion (default False)
        snap_distance: Vertex snap distance for boundary stitching (default 0.01).
                      Vertices within this distance are merged together at tile boundaries.
                      Use smaller values for precise geometry (0.001),
                      larger for coarse meshes (0.1).

    Returns:
        Merged mesh containing the entire city

    Memory Behavior:
        - With use_disk_store=True: Memory bounded by tile_size, scales to arbitrarily large cities
        - With use_disk_store=False: Memory = tile_size * num_workers (can cause OOM on large cities)
    """
    if n_workers is None:
        n_workers = max(1, len(psutil.Process().cpu_affinity()))

    city_bounds = city.calculate_bounds()
    info(f"Building tiled mesh for city with bounds: {city_bounds}")

    # Create tiles
    tiles = city_bounds.tiles(tile_size)
    info(f"Created {len(tiles)} tiles of size {tile_size}x{tile_size}")

    # Prepare mesh parameters
    mesh_params = {
        'lod': lod,
        'min_building_detail': min_building_detail,
        'min_building_area': min_building_area,
        'merge_buildings': merge_buildings,
        'merge_tolerance': merge_tolerance,
        'building_mesh_triangle_size': building_mesh_triangle_size,
        'max_mesh_size': max_mesh_size,
        'min_mesh_angle': min_mesh_angle,
        'smoothing': smoothing,
        'sort_triangles': sort_triangles,
    }

    # Prepare tile arguments
    tile_args = []
    for tile_id, tile_bounds in enumerate(tiles):
        # Filter buildings for this tile (with overlap)
        from copy import deepcopy
        expanded_bounds = deepcopy(tile_bounds)
        expanded_bounds.buffer(overlap)

        buildings_in_tile = []
        for b in city.buildings:
            b_bounds = b.calculate_bounds()
            # Check if building bounds overlap with expanded tile bounds
            if not (b_bounds.xmax < expanded_bounds.xmin or
                    b_bounds.xmin > expanded_bounds.xmax or
                    b_bounds.ymax < expanded_bounds.ymin or
                    b_bounds.ymin > expanded_bounds.ymax):
                buildings_in_tile.append(b)

        if len(buildings_in_tile) == 0:
            continue  # Skip empty tiles

        tile_args.append((
            tile_id,
            tile_bounds.tuple,
            buildings_in_tile,
            city.terrain.raster if city.has_terrain() else None,
            mesh_params
        ))

    if len(tile_args) == 0:
        warning("No tiles have buildings")
        return Mesh(vertices=np.array([]), faces=np.array([], dtype=np.int32))

    # Use disk-backed storage if requested
    if use_disk_store:
        store = SharedMeshStore(temp_dir=temp_dir, cleanup_on_delete=not keep_intermediate)
        info(f"Using disk-backed HDF5 storage for memory efficiency")
        info(f"Temporary file: {store.store_file}")
    else:
        store = None

    # Process tiles in parallel
    valid_tile_count = 0
    with Pool(n_workers) as pool:
        for tile_id, vertices, faces, markers, bounds in tqdm(
            pool.imap_unordered(_process_mesh_tile, tile_args),
            total=len(tile_args),
            desc="Processing tiles",
            disable=not progress
        ):
            if len(vertices) > 0:
                if use_disk_store:
                    # Write to disk immediately, keep out of memory
                    store.add_tile(tile_id, bounds, vertices, faces, markers)
                    valid_tile_count += 1
                    # Explicitly free memory
                    del vertices, faces, markers
                else:
                    # Keep in memory (for backward compatibility)
                    pass

    if use_disk_store:
        if valid_tile_count == 0:
            warning("No valid meshes generated from tiles")
            store.clear()
            return Mesh(vertices=np.array([]), faces=np.array([], dtype=np.int32))

        # Merge from disk incrementally
        info(f"Merging {valid_tile_count} tile meshes from disk")
        result_mesh = _merge_meshes_from_disk(store, progress=progress, snap_distance=snap_distance)
        return result_mesh
    else:
        # Original in-memory approach (for backward compatibility)
        meshes = []
        with Pool(n_workers) as pool:
            results = list(tqdm(
                pool.imap_unordered(_process_mesh_tile, tile_args),
                total=len(tile_args),
                desc="Processing tiles",
                disable=not progress
            ))

        # Collect non-empty meshes
        for tile_id, vertices, faces, markers, bounds in results:
            if len(vertices) > 0:
                mesh = Mesh(vertices=vertices, faces=faces)
                if markers is not None:
                    mesh.markers = markers
                meshes.append(mesh)

        if len(meshes) == 0:
            warning("No valid meshes generated from tiles")
            return Mesh(vertices=np.array([]), faces=np.array([], dtype=np.int32))

        # Merge meshes incrementally
        info(f"Merging {len(meshes)} tile meshes")
        return _merge_meshes_incremental(meshes, progress=progress)


def _merge_meshes_from_disk(
    store: SharedMeshStore,
    batch_size: int = 2,
    progress: bool = True,
    snap_distance: float = 0.01,
) -> Mesh:
    """
    Merge mesh tiles from disk-backed HDF5 storage with boundary snapping.

    Loads tiles in batches from disk, merges them with vertex snapping at boundaries,
    and writes results back to disk. This approach keeps memory usage bounded
    regardless of number of tiles while properly handling mesh stitching.

    Args:
        store: SharedMeshStore instance with tile data
        batch_size: Number of tiles to merge at once (default 2, keep memory low)
        progress: Whether to show progress bar
        snap_distance: Vertex snap distance for boundary stitching (default 0.01).
                      Vertices within this distance are merged together.
                      Use smaller values for precise geometry (0.001),
                      larger for coarse meshes (0.1).

    Returns:
        Final merged mesh from disk

    Note:
        Boundary handling: When tile boundaries cut through buildings, the overlap
        region ensures buildings are included in adjacent tiles. The snap_distance
        parameter controls how aggressively vertices at boundaries are merged.
    """
    from .meshing import (
        merge_meshes as merge_meshes_func,
        remove_degenerate_faces,
        remove_duplicate_faces,
        remove_internal_faces,
    )

    tiles = store.list_tiles()
    merge_iteration = 0

    while len(tiles) > 1:
        merge_iteration += 1
        info(f"Merge iteration {merge_iteration}: {len(tiles)} tiles remaining")

        new_tiles = []

        for i in tqdm(
            range(0, len(tiles), batch_size),
            desc=f"Merging batch (iteration {merge_iteration})",
            disable=not progress
        ):
            batch_ids = tiles[i:i + batch_size]

            # Load only this batch from disk
            batch_meshes = []
            for tile_id in batch_ids:
                vertices, faces, markers = store.get_tile(tile_id)
                if len(vertices) > 0:
                    mesh = Mesh(vertices=vertices, faces=faces)
                    if markers is not None:
                        mesh.markers = markers
                    batch_meshes.append(mesh)

            if len(batch_meshes) == 0:
                continue

            # Merge this batch with boundary snapping
            if len(batch_meshes) == 1:
                merged = batch_meshes[0]
            else:
                # Use weld=True to remove duplicate vertices
                # Use snap to merge vertices within snap_distance (for boundary stitching)
                merged = merge_meshes_func(batch_meshes, weld=True, snap=snap_distance)

            # Clean up duplicate/internal faces created by the merge
            merged = remove_degenerate_faces(merged, min_area=1e-8)
            merged = remove_duplicate_faces(merged)
            merged = remove_internal_faces(merged, angle_threshold=170.0)

            # Write merged result back to disk with new ID
            # Use high IDs to avoid collisions
            new_tile_id = max(store.tile_metadata.keys()) + 1000 if store.tile_metadata else 1000
            bounds = (0.0, 0.0, 1.0, 1.0)  # Bounds no longer meaningful after merge
            store.add_tile(new_tile_id, bounds, merged.vertices, merged.faces,
                          merged.markers if hasattr(merged, 'markers') else None)
            new_tiles.append(new_tile_id)

            # Explicitly free memory
            del batch_meshes
            del merged

        tiles = new_tiles

    # Load final result
    if len(tiles) > 0:
        final_id = tiles[0]
        vertices, faces, markers = store.get_tile(final_id)
        return Mesh(vertices=vertices, faces=faces)
    else:
        return Mesh(vertices=np.array([]), faces=np.array([], dtype=np.int32))


def _merge_meshes_incremental(meshes: List[Mesh], batch_size: int = 4, progress: bool = True) -> Mesh:
    """
    Merge meshes incrementally to avoid memory spikes.

    Args:
        meshes: List of meshes to merge
        batch_size: Number of meshes to merge at once
        progress: Whether to show progress bar

    Returns:
        Single merged mesh
    """
    from .meshing import merge_meshes as merge_meshes_func

    if len(meshes) == 0:
        return Mesh(vertices=np.array([]), faces=np.array([], dtype=np.int32))

    if len(meshes) == 1:
        return meshes[0]

    current_meshes = meshes.copy()

    while len(current_meshes) > 1:
        new_meshes = []

        # Merge in batches
        for i in tqdm(
            range(0, len(current_meshes), batch_size),
            desc="Merging mesh batches",
            disable=not progress
        ):
            batch = current_meshes[i:i + batch_size]
            if len(batch) == 1:
                new_meshes.append(batch[0])
            else:
                merged = merge_meshes_func(batch, weld=True)
                new_meshes.append(merged)

        current_meshes = new_meshes

    return current_meshes[0]


def should_use_tiling(city: City, tile_size: float, available_memory_fraction: float = 0.5) -> bool:
    """
    Determine whether to use tiling based on city size and available memory.

    Args:
        city: City object to analyze
        tile_size: Proposed tile size
        available_memory_fraction: Fraction of available memory to consider

    Returns:
        True if tiling is recommended, False otherwise
    """
    if len(city.buildings) == 0:
        return False

    # Estimate memory needed
    num_buildings = len(city.buildings)
    avg_faces_per_building = 50  # Conservative estimate
    estimated_faces = num_buildings * avg_faces_per_building

    # Each face = ~40 bytes, plus overhead
    estimated_memory = estimated_faces * 40

    # Get available memory
    mem = psutil.virtual_memory()
    available = mem.available * available_memory_fraction

    # Use tiling if estimated memory exceeds available
    return estimated_memory > available
