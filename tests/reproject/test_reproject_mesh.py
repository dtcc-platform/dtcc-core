import pytest
import numpy as np
from dtcc_core.reproject.reproject import reproject_mesh
from dtcc_core.model import Mesh


@pytest.fixture
def sample_vertices():
    """Sample vertices in SWEREF99 TM (EPSG:3006)"""
    return np.array([
        [500000.0, 6500000.0, 100.0],
        [500100.0, 6500000.0, 100.0],
        [500050.0, 6500100.0, 100.0],
        [500050.0, 6500050.0, 150.0],
    ])


@pytest.fixture
def sample_faces():
    """Sample triangular faces"""
    return np.array([
        [0, 1, 2],
        [0, 1, 3],
        [0, 2, 3],
        [1, 2, 3],
    ], dtype=np.int64)


@pytest.fixture
def sample_markers():
    """Sample face markers"""
    return np.array([1, 1, 2, 2])


@pytest.fixture
def basic_mesh(sample_vertices, sample_faces):
    """Create a basic Mesh without markers"""
    return Mesh(vertices=sample_vertices, faces=sample_faces)


@pytest.fixture
def mesh_with_markers(sample_vertices, sample_faces, sample_markers):
    """Create a Mesh with markers"""
    return Mesh(vertices=sample_vertices, faces=sample_faces, markers=sample_markers)


@pytest.fixture
def mesh_with_normals(sample_vertices, sample_faces):
    """Create a Mesh with normal vectors"""
    normals = np.array([
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, -1.0],
    ])
    return Mesh(vertices=sample_vertices, faces=sample_faces, normals=normals)


class TestReprojectionBasics:
    """Test basic reprojection functionality"""

    def test_basic_reprojection_shape_preserved(self, basic_mesh):
        """Reprojection should preserve shape of mesh vertices"""
        result = reproject_mesh(basic_mesh, "EPSG:3006", "EPSG:4326")
        assert result.vertices.shape == basic_mesh.vertices.shape

    def test_faces_preserved(self, basic_mesh):
        """Faces should be preserved during reprojection"""
        result = reproject_mesh(basic_mesh, "EPSG:3006", "EPSG:4326")
        np.testing.assert_array_equal(result.faces, basic_mesh.faces)

    def test_reprojection_coordinates_change(self, basic_mesh):
        """Coordinates should change after reprojection"""
        result = reproject_mesh(basic_mesh, "EPSG:3006", "EPSG:4326")
        # Coordinates should be different (SWEREF99 vs WGS84)
        assert not np.allclose(result.vertices[:, :2], basic_mesh.vertices[:, :2])

    def test_z_coordinates_preserved(self, basic_mesh):
        """Z coordinates should be preserved during reprojection"""
        result = reproject_mesh(basic_mesh, "EPSG:3006", "EPSG:4326")
        np.testing.assert_array_equal(result.vertices[:, 2], basic_mesh.vertices[:, 2])

    def test_markers_preserved(self, mesh_with_markers):
        """Markers should be preserved during reprojection"""
        result = reproject_mesh(mesh_with_markers, "EPSG:3006", "EPSG:4326")
        np.testing.assert_array_equal(result.markers, mesh_with_markers.markers)

    def test_normals_preserved(self, mesh_with_normals):
        """Normals should be preserved during reprojection"""
        result = reproject_mesh(mesh_with_normals, "EPSG:3006", "EPSG:4326")
        np.testing.assert_array_equal(result.normals, mesh_with_normals.normals)

    def test_mesh_is_copied(self, basic_mesh):
        """Result should be a new mesh, not modifying the original"""
        original_vertices = basic_mesh.vertices.copy()
        result = reproject_mesh(basic_mesh, "EPSG:3006", "EPSG:4326")

        # Original mesh vertices should not change
        np.testing.assert_array_equal(basic_mesh.vertices, original_vertices)

        # Result should have different vertices
        assert not np.allclose(result.vertices, original_vertices)


class TestEmptyMesh:
    """Test edge cases with empty or minimal meshes"""

    def test_empty_mesh(self):
        """Test reprojection with empty mesh"""
        mesh = Mesh(vertices=np.empty((0, 3)), faces=np.empty((0, 3), dtype=np.int64))
        result = reproject_mesh(mesh, "EPSG:3006", "EPSG:4326")
        assert result.vertices.shape == (0, 3)
        assert result.faces.shape == (0, 3)

    def test_single_triangle_mesh(self):
        """Test reprojection with minimal (single triangle) mesh"""
        vertices = np.array([
            [500000.0, 6500000.0, 100.0],
            [500100.0, 6500000.0, 100.0],
            [500050.0, 6500100.0, 100.0],
        ])
        faces = np.array([[0, 1, 2]], dtype=np.int64)
        mesh = Mesh(vertices=vertices, faces=faces)

        result = reproject_mesh(mesh, "EPSG:3006", "EPSG:4326")
        assert result.vertices.shape == (3, 3)
        assert result.faces.shape == (1, 3)

    def test_vertices_without_faces(self):
        """Test mesh with vertices but no faces"""
        vertices = np.array([
            [500000.0, 6500000.0, 100.0],
            [500100.0, 6500000.0, 100.0],
        ])
        mesh = Mesh(vertices=vertices, faces=np.empty((0, 3), dtype=np.int64))

        result = reproject_mesh(mesh, "EPSG:3006", "EPSG:4326")
        assert result.vertices.shape == (2, 3)
        assert result.faces.shape == (0, 3)


class TestCoordinateAccuracy:
    """Test coordinate transformation accuracy"""

    def test_known_coordinate_transformation(self):
        """Test transformation with known coordinates"""
        # Known points in SWEREF99 TM (EPSG:3006): Stockholm area
        vertices = np.array([
            [674032.0, 6580822.0, 50.0],
            [674132.0, 6580822.0, 50.0],
            [674082.0, 6580922.0, 50.0],
        ])
        faces = np.array([[0, 1, 2]], dtype=np.int64)
        mesh = Mesh(vertices=vertices, faces=faces)

        result = reproject_mesh(mesh, "EPSG:3006", "EPSG:4326")

        # All points should be roughly in Stockholm area (59.33°N, 18.06°E)
        lats = result.vertices[:, 1]
        lons = result.vertices[:, 0]
        assert np.all((59.0 < lats) & (lats < 60.0))  # Reasonable latitudes
        assert np.all((17.0 < lons) & (lons < 19.0))  # Reasonable longitudes

    def test_round_trip_transformation(self, basic_mesh):
        """Test that round-trip transformation is approximately reversible"""
        # Forward transformation
        intermediate = reproject_mesh(basic_mesh, "EPSG:3006", "EPSG:4326")

        # Backward transformation
        result = reproject_mesh(intermediate, "EPSG:4326", "EPSG:3006")

        # Should be very close to original (within reasonable tolerance)
        np.testing.assert_allclose(result.vertices, basic_mesh.vertices, rtol=1e-5, atol=1e-3)

    def test_round_trip_preserves_faces(self, basic_mesh):
        """Test round-trip transformation preserves face connectivity"""
        # Forward transformation
        intermediate = reproject_mesh(basic_mesh, "EPSG:3006", "EPSG:4326")

        # Backward transformation
        result = reproject_mesh(intermediate, "EPSG:4326", "EPSG:3006")

        # Faces should be exactly preserved
        np.testing.assert_array_equal(result.faces, basic_mesh.faces)


class TestDifferentCRSTypes:
    """Test reprojection between different types of CRS"""

    def test_projected_to_geographic(self, basic_mesh):
        """Test reprojection from projected (EPSG:3006) to geographic (EPSG:4326)"""
        result = reproject_mesh(basic_mesh, "EPSG:3006", "EPSG:4326")
        # Geographic coordinates should be in reasonable ranges
        assert np.all((-180 <= result.vertices[:, 0]) & (result.vertices[:, 0] <= 180))
        assert np.all((-90 <= result.vertices[:, 1]) & (result.vertices[:, 1] <= 90))

    def test_projected_to_projected(self, basic_mesh):
        """Test reprojection between projected CRS (EPSG:3006 to EPSG:32633 - UTM Zone 33N)"""
        result = reproject_mesh(basic_mesh, "EPSG:3006", "EPSG:32633")
        # Coordinates should change (even if slightly) but remain in projected ranges
        # Check that arrays are not exactly equal
        assert not np.array_equal(result.vertices[:, :2], basic_mesh.vertices[:, :2])
        assert result.vertices.shape == basic_mesh.vertices.shape

    def test_utm_zones(self):
        """Test reprojection between different UTM zones"""
        vertices = np.array([
            [500000.0, 6500000.0, 100.0],
            [500100.0, 6500000.0, 100.0],
            [500050.0, 6500100.0, 100.0],
        ])
        faces = np.array([[0, 1, 2]], dtype=np.int64)
        mesh = Mesh(vertices=vertices, faces=faces)

        # Reproject to UTM Zone 32N
        result = reproject_mesh(mesh, "EPSG:32633", "EPSG:32632")
        assert result.vertices.shape == vertices.shape
        assert not np.allclose(result.vertices[:, :2], vertices[:, :2])


class TestLargeMeshes:
    """Test with large meshes"""

    def test_large_mesh_with_many_vertices(self):
        """Test reprojection with a mesh having many vertices"""
        # Create a mesh with 1000 vertices
        n = 1000
        t = np.linspace(0, 2 * np.pi, n)
        x = 500000 + 100 * np.cos(t)
        y = 6500000 + 100 * np.sin(t)
        z = 100 + 10 * np.sin(5 * t)  # Varying z
        vertices = np.column_stack((x, y, z))

        # Create faces (simplified triangulation)
        faces = []
        for i in range(n - 2):
            faces.append([0, i + 1, i + 2])
        faces = np.array(faces, dtype=np.int64)

        mesh = Mesh(vertices=vertices, faces=faces)
        result = reproject_mesh(mesh, "EPSG:3006", "EPSG:4326")

        assert result.vertices.shape == (n, 3)
        assert result.faces.shape == faces.shape
        assert not np.allclose(result.vertices[:, :2], vertices[:, :2])

    def test_complex_mesh_with_many_faces(self):
        """Test mesh with many faces"""
        # Create a grid of vertices
        grid_size = 10
        vertices = []
        for i in range(grid_size):
            for j in range(grid_size):
                x = 500000.0 + i * 10.0
                y = 6500000.0 + j * 10.0
                z = 100.0 + np.sin(i * 0.5) * np.cos(j * 0.5) * 10.0
                vertices.append([x, y, z])
        vertices = np.array(vertices)

        # Create faces from grid
        faces = []
        for i in range(grid_size - 1):
            for j in range(grid_size - 1):
                # Each grid cell becomes 2 triangles
                v0 = i * grid_size + j
                v1 = i * grid_size + (j + 1)
                v2 = (i + 1) * grid_size + j
                v3 = (i + 1) * grid_size + (j + 1)

                faces.append([v0, v1, v2])
                faces.append([v1, v3, v2])
        faces = np.array(faces, dtype=np.int64)

        mesh = Mesh(vertices=vertices, faces=faces)
        result = reproject_mesh(mesh, "EPSG:3006", "EPSG:4326")

        assert result.vertices.shape == vertices.shape
        assert result.faces.shape == faces.shape
        np.testing.assert_array_equal(result.faces, mesh.faces)


class TestTopologyPreservation:
    """Test that mesh topology is preserved during reprojection"""

    def test_face_orientation_preserved(self):
        """Test that face vertex order (orientation) is preserved"""
        vertices = np.array([
            [500000.0, 6500000.0, 100.0],
            [500100.0, 6500000.0, 100.0],
            [500100.0, 6500100.0, 100.0],
            [500000.0, 6500100.0, 100.0],
        ])
        # Explicitly ordered faces
        faces = np.array([
            [0, 1, 2],  # Counter-clockwise
            [0, 2, 3],  # Counter-clockwise
        ], dtype=np.int64)

        mesh = Mesh(vertices=vertices, faces=faces)
        result = reproject_mesh(mesh, "EPSG:3006", "EPSG:4326")

        # Face connectivity and order should be exactly preserved
        np.testing.assert_array_equal(result.faces, faces)

    def test_vertex_indices_remain_valid(self):
        """Test that all face indices remain valid after reprojection"""
        vertices = np.random.rand(50, 3) * 1000 + np.array([500000, 6500000, 100])
        faces = []
        for i in range(0, 45, 3):
            faces.append([i, i + 1, i + 2])
        faces = np.array(faces, dtype=np.int64)

        mesh = Mesh(vertices=vertices, faces=faces)
        result = reproject_mesh(mesh, "EPSG:3006", "EPSG:4326")

        # All face indices should be valid
        assert np.all(result.faces >= 0)
        assert np.all(result.faces < len(result.vertices))


class TestAttributePreservation:
    """Test that all mesh attributes are properly preserved"""

    def test_all_attributes_preserved(self):
        """Test that all optional mesh attributes are preserved"""
        vertices = np.array([
            [500000.0, 6500000.0, 100.0],
            [500100.0, 6500000.0, 100.0],
            [500050.0, 6500100.0, 100.0],
        ])
        faces = np.array([[0, 1, 2]], dtype=np.int64)
        markers = np.array([5])
        normals = np.array([[0.0, 0.0, 1.0]])

        mesh = Mesh(vertices=vertices, faces=faces, markers=markers, normals=normals)
        result = reproject_mesh(mesh, "EPSG:3006", "EPSG:4326")

        np.testing.assert_array_equal(result.faces, faces)
        np.testing.assert_array_equal(result.markers, markers)
        np.testing.assert_array_equal(result.normals, normals)

    def test_empty_attributes_preserved(self):
        """Test that empty attributes are preserved correctly"""
        vertices = np.array([
            [500000.0, 6500000.0, 100.0],
            [500100.0, 6500000.0, 100.0],
            [500050.0, 6500100.0, 100.0],
        ])
        faces = np.array([[0, 1, 2]], dtype=np.int64)

        mesh = Mesh(
            vertices=vertices,
            faces=faces,
            markers=np.empty(0),
            normals=np.empty(0)
        )
        result = reproject_mesh(mesh, "EPSG:3006", "EPSG:4326")

        assert len(result.markers) == 0
        assert len(result.normals) == 0


if __name__ == "__main__":
    pytest.main()
