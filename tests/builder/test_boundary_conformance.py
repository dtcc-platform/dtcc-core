# Copyright(C) 2026 DTCC
# Licensed under the MIT License

"""Tests for boundary conformance module."""

import numpy as np
import pytest

from dtcc_core.model.geometry.mesh import Mesh
from dtcc_core.builder.meshing.boundary_conformance import (
    conform_boundary,
    _classify_vertices,
    _find_boundary_vertices,
    _detect_contact_vertices,
    _find_shared_boundary_edge,
    _sample_old_mesh_heights,
    _build_adjacency,
    _laplacian_smooth_z,
    _compute_mean_edge_length,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_grid_mesh(nx, ny, x0, y0, dx, dy, z=0.0, marker=-2):
    """Create a regular grid mesh.

    Returns a Mesh with (nx*ny) vertices and 2*(nx-1)*(ny-1) triangular faces.
    """
    xs = np.arange(nx) * dx + x0
    ys = np.arange(ny) * dy + y0
    xx, yy = np.meshgrid(xs, ys)
    vertices = np.column_stack([xx.ravel(), yy.ravel(), np.full(nx * ny, z)])

    faces = []
    for j in range(ny - 1):
        for i in range(nx - 1):
            idx = j * nx + i
            faces.append([idx, idx + 1, idx + nx])
            faces.append([idx + 1, idx + nx + 1, idx + nx])
    faces = np.array(faces, dtype=np.int64)
    markers = np.full(len(faces), marker, dtype=np.int64)

    mesh = Mesh()
    mesh.vertices = vertices.astype(np.float64)
    mesh.faces = faces
    mesh.markers = markers
    return mesh


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def grid_mesh_a():
    """5x5 grid, x in [0,4], y in [0,4], z=0, markers -2 (ground)."""
    return _make_grid_mesh(5, 5, 0, 0, 1, 1, z=0.0, marker=-2)


@pytest.fixture
def grid_mesh_b():
    """5x5 grid, x in [4,8], y in [0,4], z=1.0, markers -2 (ground).
    Shares boundary at x=4 with grid_mesh_a."""
    return _make_grid_mesh(5, 5, 4, 0, 1, 1, z=1.0, marker=-2)


@pytest.fixture
def building_mesh():
    """Grid mesh with 2 faces marked as building (marker=0)."""
    mesh = _make_grid_mesh(5, 5, 0, 0, 1, 1, z=0.0, marker=-2)
    mesh.markers[0] = 0
    mesh.markers[1] = 0
    return mesh


@pytest.fixture
def corner_mesh_c():
    """5x5 grid, x in [4,8], y in [4,8], z=2.0.
    Shares corner with mesh_a at (4,4) and edge with mesh_b."""
    return _make_grid_mesh(5, 5, 4, 4, 1, 1, z=2.0, marker=-2)


@pytest.fixture
def single_triangle():
    """A single triangle mesh."""
    mesh = Mesh()
    mesh.vertices = np.array([[0, 0, 0], [1, 0, 0], [0.5, 1, 0]], dtype=np.float64)
    mesh.faces = np.array([[0, 1, 2]], dtype=np.int64)
    mesh.markers = np.array([-2], dtype=np.int64)
    return mesh


@pytest.fixture
def tetrahedron_surface():
    """Closed tetrahedron surface — no boundary edges."""
    mesh = Mesh()
    mesh.vertices = np.array(
        [[0, 0, 0], [1, 0, 0], [0.5, 1, 0], [0.5, 0.5, 1]], dtype=np.float64
    )
    mesh.faces = np.array(
        [[0, 1, 2], [0, 1, 3], [1, 2, 3], [0, 2, 3]], dtype=np.int64
    )
    mesh.markers = np.array([-2, -2, -2, -2], dtype=np.int64)
    return mesh


# ---------------------------------------------------------------------------
# TestClassifyVertices
# ---------------------------------------------------------------------------

class TestClassifyVertices:
    def test_all_terrain(self, grid_mesh_a):
        result = _classify_vertices(
            grid_mesh_a.faces, grid_mesh_a.markers, grid_mesh_a.num_vertices
        )
        assert result.shape == (grid_mesh_a.num_vertices,)
        assert not result.any()

    def test_with_buildings(self, building_mesh):
        result = _classify_vertices(
            building_mesh.faces, building_mesh.markers, building_mesh.num_vertices
        )
        # Faces 0 and 1 are building → their vertices should be True
        building_verts = np.unique(building_mesh.faces[:2].ravel())
        for vi in building_verts:
            assert result[vi], f"Vertex {vi} should be classified as building"

    def test_empty_markers(self, grid_mesh_a):
        result = _classify_vertices(
            grid_mesh_a.faces, np.empty(0), grid_mesh_a.num_vertices
        )
        assert not result.any()

    def test_marker_length_mismatch(self, grid_mesh_a):
        bad_markers = np.array([-2, -2], dtype=np.int64)
        result = _classify_vertices(
            grid_mesh_a.faces, bad_markers, grid_mesh_a.num_vertices
        )
        assert not result.any()


# ---------------------------------------------------------------------------
# TestFindSharedBoundaryEdge
# ---------------------------------------------------------------------------

class TestFindSharedBoundaryEdge:
    def test_x_max_shared(self, grid_mesh_a, grid_mesh_b):
        """Detect x_max edge of mesh_a shared with x_min of mesh_b."""
        edge_type, edge_value = _find_shared_boundary_edge(
            grid_mesh_a.vertices, grid_mesh_b, tolerance=0.1
        )
        assert edge_type == 'x_max'
        np.testing.assert_allclose(edge_value, 4.0, atol=1e-6)

    def test_y_max_shared(self):
        """Detect y_max edge shared."""
        # Mesh at y=[0,4], another at y=[4,8]
        old_mesh = _make_grid_mesh(5, 5, 0, 0, 1, 1, z=0.0, marker=-2)
        new_mesh = _make_grid_mesh(5, 5, 0, 4, 1, 1, z=1.0, marker=-2)
        edge_type, edge_value = _find_shared_boundary_edge(
            new_mesh.vertices, old_mesh, tolerance=0.1
        )
        assert edge_type == 'y_min'
        np.testing.assert_allclose(edge_value, 4.0, atol=1e-6)

    def test_no_shared_edge(self):
        """No shared edge when meshes are far apart."""
        old_mesh = _make_grid_mesh(5, 5, 0, 0, 1, 1, z=0.0, marker=-2)
        far_mesh = _make_grid_mesh(5, 5, 100, 100, 1, 1, z=1.0, marker=-2)
        edge_type, edge_value = _find_shared_boundary_edge(
            far_mesh.vertices, old_mesh, tolerance=1.0
        )
        assert edge_type is None
        assert edge_value is None


# ---------------------------------------------------------------------------
# TestFindBoundaryVertices
# ---------------------------------------------------------------------------

class TestFindBoundaryVertices:
    def test_open_grid(self, grid_mesh_a):
        result = _find_boundary_vertices(grid_mesh_a.faces, grid_mesh_a.num_vertices)
        # Perimeter vertices of a 5x5 grid
        xs = grid_mesh_a.vertices[:, 0]
        ys = grid_mesh_a.vertices[:, 1]
        on_border = (xs == 0) | (xs == 4) | (ys == 0) | (ys == 4)
        np.testing.assert_array_equal(result, on_border)

    def test_single_triangle(self, single_triangle):
        result = _find_boundary_vertices(
            single_triangle.faces, single_triangle.num_vertices
        )
        assert result.all(), "All vertices of a single triangle are boundary"

    def test_closed_tetrahedron(self, tetrahedron_surface):
        result = _find_boundary_vertices(
            tetrahedron_surface.faces, tetrahedron_surface.num_vertices
        )
        assert not result.any(), "Closed tetrahedron has no boundary vertices"


# ---------------------------------------------------------------------------
# TestBuildAdjacency
# ---------------------------------------------------------------------------

class TestBuildAdjacency:
    def test_single_triangle(self, single_triangle):
        adj = _build_adjacency(single_triangle.faces, single_triangle.num_vertices)
        assert adj.shape == (3, 3)
        # Every pair is connected
        dense = adj.toarray()
        assert (dense == dense.T).all(), "Adjacency must be symmetric"
        for i in range(3):
            for j in range(3):
                if i != j:
                    assert dense[i, j] > 0

    def test_grid_symmetry(self, grid_mesh_a):
        adj = _build_adjacency(grid_mesh_a.faces, grid_mesh_a.num_vertices)
        dense = adj.toarray()
        np.testing.assert_array_equal(dense, dense.T)


# ---------------------------------------------------------------------------
# TestSampleOldMeshHeights
# ---------------------------------------------------------------------------

class TestSampleOldMeshHeights:
    def test_point_inside_triangle(self):
        mesh = Mesh()
        mesh.vertices = np.array(
            [[0, 0, 0], [2, 0, 1], [1, 2, 2]], dtype=np.float64
        )
        mesh.faces = np.array([[0, 1, 2]], dtype=np.int64)
        mesh.markers = np.array([-2], dtype=np.int64)
        # Centroid is at (1, 2/3)
        query = np.array([[1.0, 2.0 / 3.0]])
        z, valid = _sample_old_mesh_heights(query, mesh)
        assert valid[0]
        expected_z = (0 + 1 + 2) / 3.0
        np.testing.assert_allclose(z[0], expected_z, atol=1e-10)

    def test_point_on_edge(self):
        mesh = Mesh()
        mesh.vertices = np.array(
            [[0, 0, 0], [2, 0, 4], [1, 2, 2]], dtype=np.float64
        )
        mesh.faces = np.array([[0, 1, 2]], dtype=np.int64)
        mesh.markers = np.array([-2], dtype=np.int64)
        # Midpoint of edge v0-v1 = (1, 0)
        query = np.array([[1.0, 0.0]])
        z, valid = _sample_old_mesh_heights(query, mesh)
        np.testing.assert_allclose(z[0], 2.0, atol=1e-10)

    def test_point_outside_fallback(self):
        mesh = Mesh()
        mesh.vertices = np.array(
            [[0, 0, 5], [2, 0, 5], [1, 2, 5]], dtype=np.float64
        )
        mesh.faces = np.array([[0, 1, 2]], dtype=np.int64)
        mesh.markers = np.array([-2], dtype=np.int64)
        # Point far outside
        query = np.array([[10.0, 10.0]])
        z, valid = _sample_old_mesh_heights(query, mesh)
        # Should fall back to nearest vertex z = 5.0
        assert not valid[0]
        np.testing.assert_allclose(z[0], 5.0, atol=1e-10)


# ---------------------------------------------------------------------------
# TestLaplacianSmoothZ
# ---------------------------------------------------------------------------

class TestLaplacianSmoothZ:
    def test_frozen_unchanged(self, grid_mesh_a):
        adj = _build_adjacency(grid_mesh_a.faces, grid_mesh_a.num_vertices)
        # Set varying z
        verts = grid_mesh_a.vertices.copy()
        verts[:, 2] = np.random.RandomState(42).rand(len(verts))
        freeze = np.ones(len(verts), dtype=bool)  # freeze all
        z = _laplacian_smooth_z(verts, adj, freeze, 20)
        np.testing.assert_array_equal(z, verts[:, 2])

    def test_converges_to_flat(self):
        """Free vertices between two frozen boundaries converge to linear."""
        # 5x2 grid gives a strip with 5 columns of vertices (indices 0-4, 5-9)
        mesh = _make_grid_mesh(5, 2, 0, 0, 1, 1, z=0.0)
        # Pin left column (x=0) at z=0 and right column (x=4) at z=4
        for i in range(mesh.num_vertices):
            mesh.vertices[i, 2] = mesh.vertices[i, 0]  # z = x initially for pinned
        # Reset interior to 0 so smoothing has work to do
        left_col = mesh.vertices[:, 0] == 0
        right_col = mesh.vertices[:, 0] == 4
        interior = ~left_col & ~right_col
        mesh.vertices[interior, 2] = 0.0

        adj = _build_adjacency(mesh.faces, mesh.num_vertices)
        freeze = left_col | right_col
        z = _laplacian_smooth_z(mesh.vertices, adj, freeze, 200)
        # Frozen endpoints preserved
        np.testing.assert_array_equal(z[left_col], 0.0)
        np.testing.assert_array_equal(z[right_col], 4.0)
        # Interior should approach linear interpolation (z ≈ x)
        for i in np.where(interior)[0]:
            expected_z = mesh.vertices[i, 0]
            np.testing.assert_allclose(z[i], expected_z, atol=0.3)

    def test_zero_iterations(self, grid_mesh_a):
        adj = _build_adjacency(grid_mesh_a.faces, grid_mesh_a.num_vertices)
        verts = grid_mesh_a.vertices.copy()
        verts[:, 2] = np.arange(len(verts), dtype=float)
        freeze = np.zeros(len(verts), dtype=bool)
        z = _laplacian_smooth_z(verts, adj, freeze, 0)
        np.testing.assert_array_equal(z, verts[:, 2])


# ---------------------------------------------------------------------------
# TestComputeMeanEdgeLength
# ---------------------------------------------------------------------------

class TestComputeMeanEdgeLength:
    def test_unit_square(self):
        mesh = _make_grid_mesh(2, 2, 0, 0, 1, 1, z=0.0)
        mel = _compute_mean_edge_length(mesh.faces, mesh.vertices)
        # Edges: 4 unit-length + 1 diagonal (sqrt(2)) ... but counted per face
        assert mel > 0
        assert mel < 2.0

    def test_equilateral(self):
        verts = np.array([[0, 0, 0], [1, 0, 0], [0.5, np.sqrt(3)/2, 0]])
        faces = np.array([[0, 1, 2]])
        mel = _compute_mean_edge_length(faces, verts)
        np.testing.assert_allclose(mel, 1.0, atol=1e-10)


# ---------------------------------------------------------------------------
# TestConformBoundary (integration)
# ---------------------------------------------------------------------------

class TestConformBoundary:
    def test_basic_conformance(self, grid_mesh_b, grid_mesh_a):
        """Boundary z matches old mesh after conformance."""
        result = conform_boundary(grid_mesh_b, [grid_mesh_a])
        # Vertices at x=4 in result should have z close to 0 (from grid_mesh_a)
        at_boundary = np.isclose(result.vertices[:, 0], 4.0)
        if at_boundary.any():
            np.testing.assert_allclose(
                result.vertices[at_boundary, 2], 0.0, atol=0.1
            )

    def test_building_vertices_preserved(self, building_mesh, grid_mesh_b):
        """Building vertex z unchanged."""
        building_verts_before = building_mesh.vertices.copy()
        is_building = _classify_vertices(
            building_mesh.faces, building_mesh.markers, building_mesh.num_vertices
        )
        result = conform_boundary(building_mesh, [grid_mesh_b])
        np.testing.assert_array_equal(
            result.vertices[is_building, 2],
            building_verts_before[is_building, 2],
        )

    def test_no_old_meshes(self, grid_mesh_a):
        """Returns copy unchanged."""
        result = conform_boundary(grid_mesh_a, [])
        np.testing.assert_array_equal(result.vertices, grid_mesh_a.vertices)

    def test_no_contact(self, grid_mesh_a):
        """Far-away old mesh, no pinning, copy unchanged."""
        far_mesh = _make_grid_mesh(5, 5, 100, 100, 1, 1, z=5.0)
        result = conform_boundary(grid_mesh_a, [far_mesh])
        np.testing.assert_array_equal(result.vertices, grid_mesh_a.vertices)

    def test_corner_adjacency(self, grid_mesh_b, grid_mesh_a, corner_mesh_c):
        """Two old meshes at corner — averaged z."""
        result = conform_boundary(grid_mesh_b, [grid_mesh_a, corner_mesh_c])
        assert result is not grid_mesh_b

    def test_input_not_mutated(self, grid_mesh_b, grid_mesh_a):
        """Original mesh untouched."""
        original_verts = grid_mesh_b.vertices.copy()
        conform_boundary(grid_mesh_b, [grid_mesh_a])
        np.testing.assert_array_equal(grid_mesh_b.vertices, original_verts)

    def test_auto_tolerance(self, grid_mesh_b, grid_mesh_a):
        """Works when tolerance=None."""
        result = conform_boundary(grid_mesh_b, [grid_mesh_a], tolerance=None)
        assert result.num_vertices == grid_mesh_b.num_vertices

    def test_idempotent(self, grid_mesh_b, grid_mesh_a):
        """Second call converges — pinned vertices remain pinned."""
        r1 = conform_boundary(grid_mesh_b, [grid_mesh_a], smoothing_iterations=200)
        r2 = conform_boundary(r1, [grid_mesh_a], smoothing_iterations=200)
        # Pinned boundary vertices must be identical
        at_boundary = np.isclose(r1.vertices[:, 0], 4.0)
        np.testing.assert_allclose(
            r1.vertices[at_boundary, 2],
            r2.vertices[at_boundary, 2],
            atol=1e-10,
        )
        # Interior converges with many iterations
        np.testing.assert_allclose(r1.vertices, r2.vertices, atol=0.05)

    def test_sparse_boundary_contact(self):
        """New boundary vertices fall between old mesh boundary vertices.

        This is the actual bug case - with vertex proximity matching, new
        boundary vertices were not detected because they didn't fall near
        old boundary vertices. Line projection solves this.
        """
        # Old mesh: coarse 3x3 grid at z=0.0
        # Boundary vertices only at positions 0, 2, 4 in both x and y
        old_mesh = _make_grid_mesh(3, 3, 0, 0, 2, 2, z=0.0, marker=-2)

        # New mesh: fine 5x5 grid at z=1.0, shares x=4 boundary with old mesh
        # New boundary vertices at y=0, 1, 2, 3, 4
        # With line projection, ALL should be pinned regardless of old mesh sampling
        new_mesh = _make_grid_mesh(5, 5, 4, 0, 1, 1, z=1.0, marker=-2)

        result = conform_boundary(new_mesh, [old_mesh], tolerance=1.0)

        at_boundary = np.isclose(result.vertices[:, 0], 4.0)
        assert at_boundary.sum() == 5, "Should have 5 boundary vertices at x=4"

        # All boundary vertices should conform to old mesh height (z≈0)
        np.testing.assert_allclose(
            result.vertices[at_boundary, 2], 0.0, atol=0.15,
            err_msg="All boundary vertices should conform to old mesh height"
        )

    def test_edge_projection_accuracy(self):
        """Verify vertices are projected onto the exact edge line."""
        # Old mesh at x=[0,4]
        old_mesh = _make_grid_mesh(5, 5, 0, 0, 1, 1, z=0.0, marker=-2)

        # New mesh at x=[4,8], but with slight XY noise on boundary vertices
        new_mesh = _make_grid_mesh(5, 5, 4, 0, 1, 1, z=1.0, marker=-2)
        # Add noise to boundary vertices at x≈4
        at_boundary = np.isclose(new_mesh.vertices[:, 0], 4.0)
        np.random.seed(42)
        new_mesh.vertices[at_boundary, 0] += np.random.uniform(-0.05, 0.05, at_boundary.sum())
        new_mesh.vertices[at_boundary, 1] += np.random.uniform(-0.05, 0.05, at_boundary.sum())

        result = conform_boundary(new_mesh, [old_mesh], tolerance=0.2)

        # After conformance, boundary vertices should be projected onto EXACT x=4.0 line
        # (this happens internally during z-sampling, then Laplacian smoothing preserves it)
        at_result_boundary = np.abs(result.vertices[:, 0] - 4.0) < 0.001
        # Z-values should be pinned (smoothing iterations=0 to test pure projection)
        result_no_smooth = conform_boundary(new_mesh, [old_mesh], tolerance=0.2, smoothing_iterations=0)
        at_boundary_after = np.abs(result_no_smooth.vertices[:, 0] - 4.0) < 0.1
        if at_boundary_after.sum() > 0:
            np.testing.assert_allclose(
                result_no_smooth.vertices[at_boundary_after, 2], 0.0, atol=0.15
            )

    def test_y_axis_boundary(self):
        """Test conformance on a y-axis boundary (not just x-axis)."""
        # Old mesh at y=[0,4]
        old_mesh = _make_grid_mesh(5, 5, 0, 0, 1, 1, z=0.0, marker=-2)

        # New mesh at y=[4,8], shares y=4 boundary
        new_mesh = _make_grid_mesh(5, 5, 0, 4, 1, 1, z=2.0, marker=-2)

        result = conform_boundary(new_mesh, [old_mesh], tolerance=0.2)

        # Vertices at y≈4 should be pinned to z≈0
        at_boundary = np.isclose(result.vertices[:, 1], 4.0)
        if at_boundary.sum() > 0:
            np.testing.assert_allclose(
                result.vertices[at_boundary, 2], 0.0, atol=0.15,
                err_msg="Y-axis boundary should conform correctly"
            )

    def test_large_mesh_performance(self):
        """Verify line projection approach maintains performance."""
        # Large old mesh (10000 vertices = 100x100 grid)
        old_mesh = _make_grid_mesh(100, 100, 0, 0, 1, 1, z=0.0, marker=-2)

        # New mesh at boundary (25 vertices = 5x5 grid)
        new_mesh = _make_grid_mesh(5, 5, 99, 0, 1, 1, z=1.0, marker=-2)

        import time
        start = time.time()
        result = conform_boundary(new_mesh, [old_mesh], tolerance=2.0)
        elapsed = time.time() - start

        # Should complete quickly (< 1 second)
        # Line projection is faster than KDTree since it only processes detected edge vertices
        assert elapsed < 1.0, f"Took {elapsed:.2f}s - should be fast"

        # Verify correctness
        at_boundary = np.isclose(result.vertices[:, 0], 99.0)
        if at_boundary.sum() > 0:
            np.testing.assert_allclose(
                result.vertices[at_boundary, 2], 0.0, atol=0.15
            )
