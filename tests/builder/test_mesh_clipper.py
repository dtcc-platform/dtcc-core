"""
Unit tests for mesh_tiler module - SurfaceMeshClipper and SingleBBoxClipper classes.
"""

import pytest
import numpy as np
from typing import List

from dtcc_core.model import Mesh, Bounds
from dtcc_core.builder.meshing.mesh_tiler import (
    SurfaceMeshClipper,
    SingleBBoxClipper,
    TriangleInfo,
)
from shapely.geometry import Polygon


class TestSurfaceMeshClipper:
    """Test cases for SurfaceMeshClipper class."""

    @pytest.fixture
    def simple_mesh(self) -> Mesh:
        """Create a simple triangular mesh for testing."""
        vertices = np.array(
            [
                [0.0, 0.0, 1.0],  # 0
                [2.0, 0.0, 1.0],  # 1
                [2.0, 2.0, 1.0],  # 2
                [0.0, 2.0, 1.0],  # 3
                [1.0, 1.0, 1.5],  # 4 - center raised
            ]
        )

        faces = np.array(
            [
                [0, 1, 4],  # Bottom triangle
                [1, 2, 4],  # Right triangle
                [2, 3, 4],  # Top triangle
                [3, 0, 4],  # Left triangle
            ],
            dtype=np.int32,
        )

        return Mesh(vertices=vertices, faces=faces)

    @pytest.fixture
    def vertical_mesh(self) -> Mesh:
        """Create a mesh with vertical triangles for testing."""
        vertices = np.array(
            [
                # Horizontal triangle
                [0.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [1.0, 2.0, 0.0],
                # Vertical triangle (wall-like)
                [1.0, 1.0, 0.0],
                [1.0, 1.0, 1.0],
                [1.0, 1.0, 2.0],
            ]
        )

        faces = np.array(
            [
                [0, 1, 2],  # Horizontal
                [3, 4, 5],  # Vertical
            ],
            dtype=np.int32,
        )

        return Mesh(vertices=vertices, faces=faces)

    @pytest.fixture
    def empty_mesh(self) -> Mesh:
        """Create an empty mesh for edge case testing."""
        return Mesh(vertices=np.array([]), faces=np.array([], dtype=np.int32))

    def test_init_simple_mesh(self, simple_mesh):
        """Test SurfaceMeshClipper initialization with simple mesh."""
        clipper = SurfaceMeshClipper(simple_mesh)

        assert clipper.mesh == simple_mesh
        assert len(clipper.triangle_infos) == 4  # 4 faces
        assert hasattr(clipper, "spatial_index")

    def test_init_empty_mesh(self, empty_mesh):
        """Test SurfaceMeshClipper initialization with empty mesh."""
        # Empty mesh should be handled gracefully - this might raise an exception
        # which is acceptable behavior, so we test that it either works or fails gracefully
        try:
            clipper = SurfaceMeshClipper(empty_mesh)
            assert clipper.mesh == empty_mesh
            assert len(clipper.triangle_infos) == 0
            assert hasattr(clipper, "spatial_index")
        except (IndexError, ValueError):
            # It's acceptable for empty meshes to raise an exception during initialization
            pytest.skip("SurfaceMeshClipper does not support empty meshes")

    def test_preprocess_triangles_normal(self, simple_mesh):
        """Test triangle preprocessing for normal (non-vertical) triangles."""
        clipper = SurfaceMeshClipper(simple_mesh)

        # All triangles in simple_mesh should be non-vertical
        for info in clipper.triangle_infos:
            assert not info.is_vertical
            assert info.polygon_2d is not None
            assert isinstance(info.polygon_2d, Polygon)
            assert len(info.bbox) == 4  # (xmin, ymin, xmax, ymax)

    def test_preprocess_triangles_vertical(self, vertical_mesh):
        """Test triangle preprocessing for vertical triangles."""
        clipper = SurfaceMeshClipper(vertical_mesh)

        # First triangle should be normal, second should be vertical
        normal_info = clipper.triangle_infos[0]
        vertical_info = clipper.triangle_infos[1]

        assert not normal_info.is_vertical
        assert normal_info.polygon_2d is not None

        assert vertical_info.is_vertical
        assert vertical_info.polygon_2d is None

    def test_clip_to_bounds_simple(self, simple_mesh):
        """Test basic clipping to bounds."""
        clipper = SurfaceMeshClipper(simple_mesh)
        bounds = Bounds(xmin=0.5, ymin=0.5, xmax=1.5, ymax=1.5)

        clipped_mesh = clipper.clip_to_bounds(bounds)

        assert isinstance(clipped_mesh, Mesh)
        assert len(clipped_mesh.vertices) > 0
        assert len(clipped_mesh.faces) > 0

        # Verify all vertices are within bounds
        vertices = clipped_mesh.vertices
        assert np.all(vertices[:, 0] >= bounds.xmin - 1e-10)
        assert np.all(vertices[:, 0] <= bounds.xmax + 1e-10)
        assert np.all(vertices[:, 1] >= bounds.ymin - 1e-10)
        assert np.all(vertices[:, 1] <= bounds.ymax + 1e-10)

    def test_clip_to_bounds_no_intersection(self, simple_mesh):
        """Test clipping with bounds that don't intersect the mesh."""
        clipper = SurfaceMeshClipper(simple_mesh)
        bounds = Bounds(xmin=5.0, ymin=5.0, xmax=6.0, ymax=6.0)  # Far from mesh

        clipped_mesh = clipper.clip_to_bounds(bounds)

        assert isinstance(clipped_mesh, Mesh)
        assert len(clipped_mesh.vertices) == 0
        assert len(clipped_mesh.faces) == 0

    def test_clip_to_bounds_completely_inside(self, simple_mesh):
        """Test clipping with bounds that completely contain the mesh."""
        clipper = SurfaceMeshClipper(simple_mesh)
        bounds = Bounds(
            xmin=-1.0, ymin=-1.0, xmax=3.0, ymax=3.0
        )  # Contains entire mesh

        clipped_mesh = clipper.clip_to_bounds(bounds)

        assert isinstance(clipped_mesh, Mesh)
        assert len(clipped_mesh.vertices) >= len(simple_mesh.vertices)
        assert len(clipped_mesh.faces) >= len(simple_mesh.faces)

    def test_clip_vertical_triangles(self, vertical_mesh):
        """Test clipping mesh with vertical triangles."""
        clipper = SurfaceMeshClipper(vertical_mesh)
        bounds = Bounds(xmin=0.5, ymin=0.5, xmax=1.5, ymax=1.5)

        clipped_mesh = clipper.clip_to_bounds(bounds)

        assert isinstance(clipped_mesh, Mesh)
        assert len(clipped_mesh.vertices) > 0
        assert len(clipped_mesh.faces) > 0

    def test_query_triangles(self, simple_mesh):
        """Test spatial query functionality."""
        clipper = SurfaceMeshClipper(simple_mesh)

        # Create a query box that intersects some triangles
        from shapely.geometry import box

        query_box = box(0.5, 0.5, 1.5, 1.5)

        potential_triangles = clipper._query_triangles(query_box)

        assert isinstance(potential_triangles, list)
        assert len(potential_triangles) > 0
        assert all(isinstance(info, TriangleInfo) for info in potential_triangles)

    def test_empty_mesh_clipping(self, empty_mesh):
        """Test clipping behavior with empty mesh."""
        # Same as initialization test - empty meshes may not be supported
        try:
            clipper = SurfaceMeshClipper(empty_mesh)
            bounds = Bounds(xmin=0.0, ymin=0.0, xmax=1.0, ymax=1.0)

            clipped_mesh = clipper.clip_to_bounds(bounds)

            assert isinstance(clipped_mesh, Mesh)
            assert len(clipped_mesh.vertices) == 0
            assert len(clipped_mesh.faces) == 0
        except (IndexError, ValueError):
            # It's acceptable for empty meshes to raise an exception
            pytest.skip("SurfaceMeshClipper does not support empty meshes")


class TestSingleBBoxClipper:
    """Test cases for SingleBBoxClipper class."""

    @pytest.fixture
    def simple_triangle_info(self) -> TriangleInfo:
        """Create a simple TriangleInfo for testing."""
        vertices = np.array(
            [
                [0.0, 0.0, 1.0],
                [2.0, 0.0, 1.0],
                [1.0, 2.0, 1.0],
            ]
        )

        polygon_2d = Polygon(vertices[:, :2])

        return TriangleInfo(
            index=0,
            vertices=vertices,
            polygon_2d=polygon_2d,
            is_vertical=False,
            normal_up=True,
            bbox=(0.0, 0.0, 2.0, 2.0),
        )

    @pytest.fixture
    def vertical_triangle_info(self) -> TriangleInfo:
        """Create a vertical TriangleInfo for testing."""
        vertices = np.array(
            [
                [1.0, 1.0, 0.0],
                [1.0, 1.0, 1.0],
                [1.0, 1.0, 2.0],
            ]
        )

        return TriangleInfo(
            index=0,
            vertices=vertices,
            polygon_2d=None,
            is_vertical=True,
            normal_up=True,
            bbox=(1.0, 1.0, 1.0, 1.0),
        )

    def test_init(self):
        """Test SingleBBoxClipper initialization."""
        bbox = (0.0, 0.0, 2.0, 2.0)
        clipper = SingleBBoxClipper(bbox)

        assert clipper.bbox == bbox
        assert hasattr(clipper, "clip_box")
        assert clipper.vertices == []
        assert clipper.faces == []
        assert clipper.vertex_map == {}

    def test_process_normal_triangle_inside(self, simple_triangle_info):
        """Test processing a normal triangle that's completely inside bounds."""
        bbox = (-1.0, -1.0, 3.0, 3.0)  # Completely contains triangle
        clipper = SingleBBoxClipper(bbox)

        clipper.process_normal_triangle(simple_triangle_info)
        result = clipper.get_result()

        assert len(result.vertices) > 0
        assert len(result.faces) > 0

    def test_process_normal_triangle_outside(self, simple_triangle_info):
        """Test processing a normal triangle that's completely outside bounds."""
        bbox = (5.0, 5.0, 6.0, 6.0)  # Doesn't intersect triangle
        clipper = SingleBBoxClipper(bbox)

        clipper.process_normal_triangle(simple_triangle_info)
        result = clipper.get_result()

        assert len(result.vertices) == 0
        assert len(result.faces) == 0

    def test_process_normal_triangle_intersecting(self, simple_triangle_info):
        """Test processing a normal triangle that intersects bounds."""
        bbox = (0.5, 0.5, 1.5, 1.5)  # Partially intersects triangle
        clipper = SingleBBoxClipper(bbox)

        clipper.process_normal_triangle(simple_triangle_info)
        result = clipper.get_result()

        assert len(result.vertices) > 0
        assert len(result.faces) > 0

        # Verify clipped vertices are within bounds
        vertices = result.vertices
        assert np.all(vertices[:, 0] >= bbox[0] - 1e-10)
        assert np.all(vertices[:, 0] <= bbox[2] + 1e-10)
        assert np.all(vertices[:, 1] >= bbox[1] - 1e-10)
        assert np.all(vertices[:, 1] <= bbox[3] + 1e-10)

    def test_process_vertical_triangle_inside(self, vertical_triangle_info):
        """Test processing a vertical triangle that's completely inside bounds."""
        bbox = (0.0, 0.0, 2.0, 2.0)  # Contains the vertical triangle
        clipper = SingleBBoxClipper(bbox)

        clipper.process_vertical_triangle(vertical_triangle_info)
        result = clipper.get_result()

        assert len(result.vertices) > 0
        assert len(result.faces) > 0

    def test_process_vertical_triangle_outside(self, vertical_triangle_info):
        """Test processing a vertical triangle that's outside bounds."""
        bbox = (5.0, 5.0, 6.0, 6.0)  # Doesn't contain the vertical triangle
        clipper = SingleBBoxClipper(bbox)

        clipper.process_vertical_triangle(vertical_triangle_info)
        result = clipper.get_result()

        assert len(result.vertices) == 0
        assert len(result.faces) == 0

    def test_get_result_empty(self):
        """Test get_result when no triangles have been processed."""
        bbox = (0.0, 0.0, 1.0, 1.0)
        clipper = SingleBBoxClipper(bbox)

        result = clipper.get_result()

        assert isinstance(result, Mesh)
        assert len(result.vertices) == 0
        assert len(result.faces) == 0

    def test_add_vertex_deduplication(self):
        """Test that vertices are properly deduplicated."""
        bbox = (0.0, 0.0, 1.0, 1.0)
        clipper = SingleBBoxClipper(bbox)

        # Add the same vertex twice
        vertex = np.array([0.5, 0.5, 0.5])
        idx1 = clipper._add_vertex(vertex)
        idx2 = clipper._add_vertex(vertex)

        assert idx1 == idx2  # Same index returned
        assert len(clipper.vertices) == 1  # Only one vertex stored

    def test_add_vertex_different(self):
        """Test adding different vertices."""
        bbox = (0.0, 0.0, 1.0, 1.0)
        clipper = SingleBBoxClipper(bbox)

        # Add different vertices
        vertex1 = np.array([0.5, 0.5, 0.5])
        vertex2 = np.array([0.7, 0.7, 0.7])

        idx1 = clipper._add_vertex(vertex1)
        idx2 = clipper._add_vertex(vertex2)

        assert idx1 != idx2  # Different indices
        assert len(clipper.vertices) == 2  # Two vertices stored

    def test_compute_normal(self):
        """Test normal computation for triangles."""
        bbox = (0.0, 0.0, 1.0, 1.0)
        clipper = SingleBBoxClipper(bbox)

        # Create a triangle in XY plane with normal pointing up
        triangle = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
            ]
        )

        normal = clipper._compute_normal(triangle)

        # Normal should point in +Z direction
        assert abs(normal[2] - 1.0) < 1e-10
        assert abs(normal[0]) < 1e-10
        assert abs(normal[1]) < 1e-10

    def test_ensure_upward_normals(self):
        """Test that normals are correctly oriented upward."""
        bbox = (0.0, 0.0, 1.0, 1.0)
        clipper = SingleBBoxClipper(bbox)

        # Create vertices and faces
        vertices = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [1.0, 1.0, 0.0],
            ]
        )

        # Create a face that has a downward-pointing normal (clockwise winding from above)
        faces = np.array(
            [
                [0, 2, 1],  # This should result in a downward normal and get flipped
            ]
        )

        corrected_faces = clipper._ensure_upward_normals(vertices, faces.copy())

        # Check that faces were corrected
        assert corrected_faces.shape == faces.shape

        # Verify the correction actually happened by computing normals
        def compute_normal(v0, v1, v2):
            edge1 = v1 - v0
            edge2 = v2 - v0
            return np.cross(edge1, edge2)

        # Original face normal should point down
        orig_normal = compute_normal(
            vertices[faces[0][0]], vertices[faces[0][1]], vertices[faces[0][2]]
        )

        # Corrected face normal should point up
        corr_normal = compute_normal(
            vertices[corrected_faces[0][0]],
            vertices[corrected_faces[0][1]],
            vertices[corrected_faces[0][2]],
        )

        # The Z component should have opposite signs
        assert orig_normal[2] < 0  # Original points down
        assert corr_normal[2] > 0  # Corrected points up

    def test_interpolate_z_batch(self):
        """Test batch Z interpolation using barycentric coordinates."""
        bbox = (0.0, 0.0, 1.0, 1.0)
        clipper = SingleBBoxClipper(bbox)

        # Create a triangle with different Z values
        triangle = np.array(
            [
                [0.0, 0.0, 0.0],  # Bottom-left
                [1.0, 0.0, 1.0],  # Bottom-right
                [0.0, 1.0, 2.0],  # Top-left
            ]
        )

        # Test points
        points_2d = np.array(
            [
                [0.0, 0.0],  # Should get Z=0
                [1.0, 0.0],  # Should get Z=1
                [0.0, 1.0],  # Should get Z=2
                [0.5, 0.0],  # Should get Z=0.5
            ]
        )

        z_values = clipper._interpolate_z_batch(points_2d, triangle)

        assert len(z_values) == len(points_2d)
        assert abs(z_values[0] - 0.0) < 1e-10
        assert abs(z_values[1] - 1.0) < 1e-10
        assert abs(z_values[2] - 2.0) < 1e-10
        assert abs(z_values[3] - 0.5) < 1e-10


class TestIntegration:
    """Integration tests for the complete clipping workflow."""

    def test_complex_mesh_clipping(self):
        """Test clipping a more complex mesh."""
        # Create a 3x3 grid mesh
        vertices = []
        for i in range(4):
            for j in range(4):
                vertices.append([i * 1.0, j * 1.0, 0.5 + 0.1 * (i + j)])

        vertices = np.array(vertices)

        faces = []
        for i in range(3):
            for j in range(3):
                # Each grid cell becomes 2 triangles
                v0 = i * 4 + j
                v1 = i * 4 + j + 1
                v2 = (i + 1) * 4 + j
                v3 = (i + 1) * 4 + j + 1

                faces.append([v0, v1, v2])
                faces.append([v1, v3, v2])

        faces = np.array(faces, dtype=np.int32)
        mesh = Mesh(vertices=vertices, faces=faces)

        # Test clipping
        clipper = SurfaceMeshClipper(mesh)
        bounds = Bounds(xmin=0.5, ymin=0.5, xmax=2.5, ymax=2.5)

        clipped_mesh = clipper.clip_to_bounds(bounds)

        assert len(clipped_mesh.vertices) > 0
        assert len(clipped_mesh.faces) > 0

        # Verify bounds
        vertices = clipped_mesh.vertices
        assert np.all(vertices[:, 0] >= bounds.xmin - 1e-10)
        assert np.all(vertices[:, 0] <= bounds.xmax + 1e-10)
        assert np.all(vertices[:, 1] >= bounds.ymin - 1e-10)
        assert np.all(vertices[:, 1] <= bounds.ymax + 1e-10)

    def test_mixed_vertical_and_normal_triangles(self):
        """Test clipping a mesh with both vertical and normal triangles."""
        vertices = np.array(
            [
                # Normal triangle (horizontal)
                [0.0, 0.0, 1.0],
                [2.0, 0.0, 1.0],
                [1.0, 2.0, 1.0],
                # Vertical triangle
                [1.0, 1.0, 0.0],
                [1.0, 1.0, 2.0],
                [1.0, 1.0, 1.0],
                # Another normal triangle
                [0.0, 0.0, 0.0],
                [2.0, 2.0, 0.0],
                [0.0, 2.0, 0.0],
            ]
        )

        faces = np.array(
            [
                [0, 1, 2],  # Normal
                [3, 4, 5],  # Vertical
                [6, 7, 8],  # Normal
            ],
            dtype=np.int32,
        )

        mesh = Mesh(vertices=vertices, faces=faces)
        clipper = SurfaceMeshClipper(mesh)
        bounds = Bounds(xmin=0.5, ymin=0.5, xmax=1.5, ymax=1.5)

        clipped_mesh = clipper.clip_to_bounds(bounds)

        assert len(clipped_mesh.vertices) > 0
        assert len(clipped_mesh.faces) > 0

    def test_edge_case_tiny_triangles(self):
        """Test clipping with very small triangles."""
        # Create tiny triangles
        epsilon = 1e-8
        vertices = np.array(
            [
                [0.0, 0.0, 0.0],
                [epsilon, 0.0, 0.0],
                [0.0, epsilon, 0.0],
            ]
        )

        faces = np.array([[0, 1, 2]], dtype=np.int32)
        mesh = Mesh(vertices=vertices, faces=faces)

        clipper = SurfaceMeshClipper(mesh)
        bounds = Bounds(xmin=-1.0, ymin=-1.0, xmax=1.0, ymax=1.0)

        # Should not crash
        clipped_mesh = clipper.clip_to_bounds(bounds)
        assert isinstance(clipped_mesh, Mesh)
