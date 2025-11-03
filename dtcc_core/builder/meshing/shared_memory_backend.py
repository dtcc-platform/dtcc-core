# Copyright (C) 2025 DTCC Contributors
# Licensed under the MIT License
# Shared memory backend for large-scale mesh processing

import h5py
import numpy as np
import tempfile
from pathlib import Path
from typing import Optional, Tuple, Dict, List
from dataclasses import dataclass
import threading


@dataclass
class TileMetadata:
    """Metadata about a mesh tile."""
    tile_id: int
    bounds: Tuple[float, float, float, float]
    num_vertices: int
    num_faces: int
    vertex_offset: int
    face_offset: int


class SharedMeshStore:
    """
    Manages shared memory storage for mesh data using HDF5 and memory-mapped arrays.

    This class allows concurrent access to large mesh data without keeping everything
    in RAM. It uses HDF5 for structured storage and memmap for efficient array access.
    """

    def __init__(self, temp_dir: Optional[str] = None, cleanup_on_delete: bool = True):
        """
        Initialize the shared mesh store.

        Args:
            temp_dir: Directory for temporary files. If None, uses system temp dir.
            cleanup_on_delete: Whether to clean up files when store is deleted.
        """
        if temp_dir is None:
            temp_dir = tempfile.gettempdir()

        self.temp_dir = Path(temp_dir)
        self.cleanup_on_delete = cleanup_on_delete

        # Create unique HDF5 file for this store
        self.store_file = self.temp_dir / f"mesh_store_{id(self)}.h5"
        self.tile_metadata: Dict[int, TileMetadata] = {}
        self.lock = threading.Lock()

        # Initialize HDF5 file structure
        with h5py.File(self.store_file, 'w') as f:
            f.create_group('tiles')
            f.create_group('metadata')
            f.attrs['version'] = '1.0'

    def add_tile(self,
                 tile_id: int,
                 bounds: Tuple[float, float, float, float],
                 vertices: np.ndarray,
                 faces: np.ndarray,
                 markers: Optional[np.ndarray] = None) -> None:
        """
        Add a mesh tile to the store.

        Args:
            tile_id: Unique identifier for the tile
            bounds: Bounding box (xmin, ymin, xmax, ymax)
            vertices: Vertex array (N, 3) with float64 dtype
            faces: Face array (M, 3) with int32 dtype
            markers: Optional marker array (M,) with int32 dtype
        """
        with self.lock:
            with h5py.File(self.store_file, 'a') as f:
                tile_group = f['tiles'].create_group(f'tile_{tile_id}')

                # Store mesh data
                tile_group.create_dataset('vertices', data=vertices, dtype='float64')
                tile_group.create_dataset('faces', data=faces, dtype='int32')
                if markers is not None:
                    tile_group.create_dataset('markers', data=markers, dtype='int32')

                # Store bounds
                tile_group.attrs['bounds'] = bounds

                # Store metadata
                metadata = TileMetadata(
                    tile_id=tile_id,
                    bounds=bounds,
                    num_vertices=len(vertices),
                    num_faces=len(faces),
                    vertex_offset=0,  # Will be computed during merge
                    face_offset=0
                )
                self.tile_metadata[tile_id] = metadata

    def get_tile(self, tile_id: int) -> Tuple[np.ndarray, np.ndarray, Optional[np.ndarray]]:
        """
        Retrieve a tile from the store.

        Args:
            tile_id: Tile identifier

        Returns:
            Tuple of (vertices, faces, markers)
        """
        with h5py.File(self.store_file, 'r') as f:
            tile_group = f['tiles'][f'tile_{tile_id}']
            vertices = np.array(tile_group['vertices'])
            faces = np.array(tile_group['faces'])
            markers = None
            if 'markers' in tile_group:
                markers = np.array(tile_group['markers'])

        return vertices, faces, markers

    def list_tiles(self) -> List[int]:
        """Return list of all tile IDs in the store."""
        return sorted(list(self.tile_metadata.keys()))

    def get_metadata(self, tile_id: int) -> TileMetadata:
        """Get metadata for a tile."""
        return self.tile_metadata[tile_id]

    def get_total_stats(self) -> Dict[str, int]:
        """Get total statistics across all tiles."""
        total_vertices = sum(m.num_vertices for m in self.tile_metadata.values())
        total_faces = sum(m.num_faces for m in self.tile_metadata.values())
        return {
            'num_tiles': len(self.tile_metadata),
            'total_vertices': total_vertices,
            'total_faces': total_faces,
        }

    def clear(self) -> None:
        """Clear all data from the store."""
        with self.lock:
            if self.store_file.exists():
                self.store_file.unlink()
            self.tile_metadata.clear()

            # Reinitialize
            with h5py.File(self.store_file, 'w') as f:
                f.create_group('tiles')
                f.create_group('metadata')

    def __del__(self):
        """Clean up temporary files."""
        if self.cleanup_on_delete and self.store_file.exists():
            try:
                self.store_file.unlink()
            except Exception:
                pass


class MemoryMappedMeshArray:
    """
    Wrapper for memory-mapped mesh arrays to handle large vertex/face arrays.
    """

    def __init__(self, filepath: str, shape: Tuple[int, int], dtype: str = 'float64'):
        """
        Create a memory-mapped array.

        Args:
            filepath: Path to the memmap file
            shape: Shape of the array
            dtype: Data type
        """
        self.filepath = filepath
        self.array = np.memmap(filepath, dtype=dtype, mode='w+', shape=shape)

    def flush(self) -> None:
        """Flush changes to disk."""
        if isinstance(self.array, np.memmap):
            self.array.flush()

    def __getitem__(self, key):
        return self.array[key]

    def __setitem__(self, key, value):
        self.array[key] = value

    def __len__(self):
        return len(self.array)

    @property
    def shape(self):
        return self.array.shape
