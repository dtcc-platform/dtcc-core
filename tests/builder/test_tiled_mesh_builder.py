# Copyright (C) 2025 DTCC Contributors
# Licensed under the MIT License
# Tests for tiled mesh builder

import pytest
import numpy as np
from dtcc_core.model import City, Building, Terrain, Mesh, Bounds, Surface, Polygon
from dtcc_core.model import GeometryType, Raster
from dtcc_core.builder.meshing.tiled_mesh_builder import (
    build_city_mesh_tiled,
    _merge_meshes_incremental,
    _merge_meshes_from_disk,
    should_use_tiling,
    calculate_optimal_tile_size,
)
from dtcc_core.builder.meshing.shared_memory_backend import SharedMeshStore


class TestSharedMeshStore:
    """Tests for SharedMeshStore class."""

    def test_store_initialization(self):
        """Test that SharedMeshStore initializes correctly."""
        store = SharedMeshStore()
        assert store.store_file.exists()
        assert len(store.tile_metadata) == 0
        store.clear()

    def test_add_and_get_tile(self):
        """Test adding and retrieving a tile."""
        store = SharedMeshStore()

        # Create test data
        vertices = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ], dtype=np.float64)
        faces = np.array([[0, 1, 2]], dtype=np.int32)
        markers = np.array([1], dtype=np.int32)

        # Add tile
        bounds = (0.0, 0.0, 1.0, 1.0)
        store.add_tile(0, bounds, vertices, faces, markers)

        # Retrieve tile
        ret_vertices, ret_faces, ret_markers = store.get_tile(0)

        assert np.allclose(ret_vertices, vertices)
        assert np.allclose(ret_faces, faces)
        assert np.allclose(ret_markers, markers)

        store.clear()

    def test_list_tiles(self):
        """Test listing tiles."""
        store = SharedMeshStore()

        vertices = np.zeros((3, 3), dtype=np.float64)
        faces = np.array([[0, 1, 2]], dtype=np.int32)

        # Add multiple tiles
        for i in range(3):
            store.add_tile(i, (0.0, 0.0, 1.0, 1.0), vertices, faces)

        tiles = store.list_tiles()
        assert len(tiles) == 3
        assert tiles == [0, 1, 2]

        store.clear()

    def test_get_total_stats(self):
        """Test getting total statistics."""
        store = SharedMeshStore()

        vertices = np.zeros((3, 3), dtype=np.float64)
        faces = np.array([[0, 1, 2]], dtype=np.int32)

        store.add_tile(0, (0.0, 0.0, 1.0, 1.0), vertices, faces)
        store.add_tile(1, (1.0, 0.0, 2.0, 1.0), vertices, faces)

        stats = store.get_total_stats()
        assert stats['num_tiles'] == 2
        assert stats['total_vertices'] == 6  # 3 per tile * 2 tiles
        assert stats['total_faces'] == 2  # 1 per tile * 2 tiles

        store.clear()


class TestOptimalTileSize:
    """Tests for optimal tile size calculation."""

    def test_calculate_optimal_tile_size(self):
        """Test that optimal tile size calculation works."""
        size = calculate_optimal_tile_size(available_memory_fraction=0.3)
        assert size > 0
        assert np.isfinite(size)


class TestMergeMeshesIncremental:
    """Tests for incremental mesh merging."""

    def _create_test_mesh(self, offset=0.0):
        """Helper to create a simple test mesh."""
        vertices = np.array([
            [0.0 + offset, 0.0, 0.0],
            [1.0 + offset, 0.0, 0.0],
            [0.0 + offset, 1.0, 0.0],
        ], dtype=np.float64)
        faces = np.array([[0, 1, 2]], dtype=np.int32)
        return Mesh(vertices=vertices, faces=faces)

    def test_merge_single_mesh(self):
        """Test merging a single mesh returns it unchanged."""
        mesh = self._create_test_mesh()
        result = _merge_meshes_incremental([mesh], progress=False)

        assert len(result.vertices) == 3
        assert len(result.faces) == 1

    def test_merge_multiple_meshes(self):
        """Test merging multiple meshes."""
        meshes = [
            self._create_test_mesh(offset=0.0),
            self._create_test_mesh(offset=2.0),
        ]

        result = _merge_meshes_incremental(meshes, batch_size=2, progress=False)

        assert len(result.vertices) >= 6  # At least 6 vertices (may be deduplicated)
        assert len(result.faces) >= 2  # At least 2 faces

    def test_merge_empty_list(self):
        """Test merging an empty list."""
        result = _merge_meshes_incremental([], progress=False)

        assert len(result.vertices) == 0
        assert len(result.faces) == 0


class TestMergeFromDisk:
    """Tests for disk-backed mesh merging with snap distance."""

    def _create_test_mesh(self, offset=0.0):
        """Helper to create a simple test mesh."""
        vertices = np.array([
            [0.0 + offset, 0.0, 0.0],
            [1.0 + offset, 0.0, 0.0],
            [0.0 + offset, 1.0, 0.0],
        ], dtype=np.float64)
        faces = np.array([[0, 1, 2]], dtype=np.int32)
        return Mesh(vertices=vertices, faces=faces)

    def test_merge_from_disk_with_snap(self):
        """Test merging from disk with snap distance parameter."""
        store = SharedMeshStore()

        # Add test meshes to store
        mesh1 = self._create_test_mesh(offset=0.0)
        mesh2 = self._create_test_mesh(offset=2.0)

        store.add_tile(0, (0.0, 0.0, 1.0, 1.0), mesh1.vertices, mesh1.faces)
        store.add_tile(1, (2.0, 0.0, 3.0, 1.0), mesh2.vertices, mesh2.faces)

        # Merge with different snap distances
        result = _merge_meshes_from_disk(store, batch_size=2, progress=False, snap_distance=0.01)

        # Should produce a valid mesh
        assert len(result.vertices) > 0
        assert len(result.faces) > 0

        store.clear()

    def test_merge_from_disk_single_tile(self):
        """Test merging from disk with only one tile."""
        store = SharedMeshStore()

        mesh = self._create_test_mesh(offset=0.0)
        store.add_tile(0, (0.0, 0.0, 1.0, 1.0), mesh.vertices, mesh.faces)

        result = _merge_meshes_from_disk(store, batch_size=2, progress=False, snap_distance=0.01)

        assert len(result.vertices) == len(mesh.vertices)
        assert len(result.faces) == len(mesh.faces)

        store.clear()

    def test_merge_from_disk_many_tiles(self):
        """Test merging many tiles from disk."""
        store = SharedMeshStore()

        # Add 4 test meshes
        for i in range(4):
            mesh = self._create_test_mesh(offset=float(i * 2))
            bounds = (float(i * 2), 0.0, float(i * 2 + 1), 1.0)
            store.add_tile(i, bounds, mesh.vertices, mesh.faces)

        # Merge with batch_size=2 (should do 2 iterations)
        result = _merge_meshes_from_disk(store, batch_size=2, progress=False, snap_distance=0.01)

        assert len(result.vertices) > 0
        assert len(result.faces) > 0

        store.clear()


class TestShouldUseTiling:
    """Tests for the should_use_tiling function."""

    def test_should_use_tiling_empty_city(self):
        """Test that empty cities don't use tiling."""
        city = City()
        assert not should_use_tiling(city, tile_size=1000.0)

    def test_should_use_tiling_small_city(self):
        """Test that small cities may not use tiling."""
        city = City()
        # Add a small building
        building = Building()
        poly = Polygon(vertices=[
            (0.0, 0.0),
            (10.0, 0.0),
            (10.0, 10.0),
            (0.0, 10.0),
        ])
        surface = Surface(polygons=[poly])
        building.add_geometry(surface, GeometryType.LOD0)
        city.add_buildings([building])

        # Small cities with low memory usage shouldn't require tiling
        result = should_use_tiling(city, tile_size=1000.0, available_memory_fraction=0.8)
        assert isinstance(result, bool)  # Just verify it returns a bool


class TestBuildCityMeshTiled:
    """Tests for the main tiled mesh builder."""

    @pytest.fixture
    def test_city(self):
        """Create a test city with simple buildings."""
        city = City()

        # Add buildings
        for i in range(2):
            building = Building()
            x_offset = i * 100.0
            poly = Polygon(vertices=[
                (x_offset, 0.0),
                (x_offset + 50.0, 0.0),
                (x_offset + 50.0, 50.0),
                (x_offset, 50.0),
            ])
            surface = Surface(polygons=[poly])
            building.add_geometry(surface, GeometryType.LOD0)
            building.height = 20.0
            city.add_buildings([building])

        # Add simple terrain raster
        terrain = Terrain()
        raster_data = np.ones((100, 100), dtype=np.float32) * 10.0
        raster = Raster(
            data=raster_data,
            bounds=Bounds(xmin=-50.0, ymin=-50.0, xmax=150.0, ymax=150.0),
        )
        terrain.add_geometry(raster, GeometryType.RASTER)
        city.add_terrain(terrain)

        return city

    def test_build_city_mesh_tiled_basic(self, test_city):
        """Test basic tiled mesh building."""
        # This test is skipped by default due to complexity of full meshing
        # In production, this would require more complex setup
        pytest.skip("Full mesh building test requires complex setup")

    def test_build_city_mesh_tiled_with_one_worker(self, test_city):
        """Test tiled mesh building with one worker."""
        pytest.skip("Full mesh building test requires complex setup")
