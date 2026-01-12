import pytest
import numpy as np
from dtcc_core.reproject.reproject import reproject_surface
from dtcc_core.model import Surface


@pytest.fixture
def sample_vertices():
    """Sample vertices in SWEREF99 TM (EPSG:3006)"""
    return np.array([
        [500000.0, 6500000.0, 100.0],
        [500100.0, 6500000.0, 100.0],
        [500100.0, 6500100.0, 100.0],
        [500000.0, 6500100.0, 100.0],
    ])


@pytest.fixture
def sample_hole():
    """Sample hole vertices in SWEREF99 TM (EPSG:3006)"""
    return np.array([
        [500025.0, 6500025.0, 100.0],
        [500075.0, 6500025.0, 100.0],
        [500075.0, 6500075.0, 100.0],
        [500025.0, 6500075.0, 100.0],
    ])


@pytest.fixture
def basic_surface(sample_vertices):
    """Create a basic Surface without holes"""
    return Surface(vertices=sample_vertices)


@pytest.fixture
def surface_with_holes(sample_vertices, sample_hole):
    """Create a Surface with holes"""
    return Surface(vertices=sample_vertices, holes=[sample_hole])


@pytest.fixture
def surface_with_multiple_holes(sample_vertices):
    """Create a Surface with multiple holes"""
    hole1 = np.array([
        [500020.0, 6500020.0, 100.0],
        [500040.0, 6500020.0, 100.0],
        [500040.0, 6500040.0, 100.0],
        [500020.0, 6500040.0, 100.0],
    ])
    hole2 = np.array([
        [500060.0, 6500060.0, 100.0],
        [500080.0, 6500060.0, 100.0],
        [500080.0, 6500080.0, 100.0],
        [500060.0, 6500080.0, 100.0],
    ])
    return Surface(vertices=sample_vertices, holes=[hole1, hole2])


class TestReprojectionBasics:
    """Test basic reprojection functionality"""

    def test_same_crs_returns_original(self, basic_surface):
        """When source and target CRS are the same, return original surface"""
        result = reproject_surface(basic_surface, "EPSG:3006", "EPSG:3006")
        assert result is basic_surface

    def test_basic_reprojection_shape_preserved(self, basic_surface):
        """Reprojection should preserve shape of surface"""
        result = reproject_surface(basic_surface, "EPSG:3006", "EPSG:4326")
        assert result.vertices.shape == basic_surface.vertices.shape

    def test_reprojection_coordinates_change(self, basic_surface):
        """Coordinates should change after reprojection"""
        result = reproject_surface(basic_surface, "EPSG:3006", "EPSG:4326")
        # Coordinates should be different (SWEREF99 vs WGS84)
        assert not np.allclose(result.vertices[:, :2], basic_surface.vertices[:, :2])

    def test_z_coordinates_preserved(self, basic_surface):
        """Z coordinates should be preserved during reprojection"""
        result = reproject_surface(basic_surface, "EPSG:3006", "EPSG:4326")
        np.testing.assert_array_equal(result.vertices[:, 2], basic_surface.vertices[:, 2])

    def test_holes_preserved(self, surface_with_holes):
        """Holes should be preserved during reprojection"""
        result = reproject_surface(surface_with_holes, "EPSG:3006", "EPSG:4326")
        assert len(result.holes) == len(surface_with_holes.holes)
        assert result.holes[0].shape == surface_with_holes.holes[0].shape

    def test_multiple_holes_preserved(self, surface_with_multiple_holes):
        """Multiple holes should be preserved during reprojection"""
        result = reproject_surface(surface_with_multiple_holes, "EPSG:3006", "EPSG:4326")
        assert len(result.holes) == 2
        for i, hole in enumerate(result.holes):
            assert hole.shape == surface_with_multiple_holes.holes[i].shape

    def test_holes_coordinates_change(self, surface_with_holes):
        """Hole coordinates should also be reprojected"""
        result = reproject_surface(surface_with_holes, "EPSG:3006", "EPSG:4326")
        # Hole coordinates should be different
        assert not np.allclose(result.holes[0][:, :2], surface_with_holes.holes[0][:, :2])

    def test_holes_z_coordinates_preserved(self, surface_with_holes):
        """Z coordinates in holes should be preserved"""
        result = reproject_surface(surface_with_holes, "EPSG:3006", "EPSG:4326")
        np.testing.assert_array_equal(result.holes[0][:, 2], surface_with_holes.holes[0][:, 2])


class TestEmptySurface:
    """Test edge cases with empty or minimal surfaces"""

    def test_empty_surface(self):
        """Test reprojection with empty surface"""
        surface = Surface(vertices=np.empty((0, 3)))
        result = reproject_surface(surface, "EPSG:3006", "EPSG:4326")
        assert result.vertices.shape == (0, 3)

    def test_surface_with_empty_holes_list(self, basic_surface):
        """Test surface with empty holes list"""
        basic_surface.holes = []
        result = reproject_surface(basic_surface, "EPSG:3006", "EPSG:4326")
        assert len(result.holes) == 0

    def test_triangle_surface(self):
        """Test reprojection with minimal (triangular) surface"""
        vertices = np.array([
            [500000.0, 6500000.0, 100.0],
            [500100.0, 6500000.0, 100.0],
            [500050.0, 6500100.0, 100.0],
        ])
        surface = Surface(vertices=vertices)
        result = reproject_surface(surface, "EPSG:3006", "EPSG:4326")
        assert result.vertices.shape == (3, 3)


class TestCoordinateAccuracy:
    """Test coordinate transformation accuracy"""

    def test_known_coordinate_transformation(self):
        """Test transformation with known coordinates"""
        # Known point in SWEREF99 TM (EPSG:3006): Stockholm area
        vertices = np.array([
            [674032.0, 6580822.0, 50.0],
            [674132.0, 6580822.0, 50.0],
            [674132.0, 6580922.0, 50.0],
            [674032.0, 6580922.0, 50.0],
        ])
        surface = Surface(vertices=vertices)

        result = reproject_surface(surface, "EPSG:3006", "EPSG:4326")

        # All points should be roughly in Stockholm area (59.33°N, 18.06°E)
        lats = result.vertices[:, 1]
        lons = result.vertices[:, 0]
        assert np.all((59.0 < lats) & (lats < 60.0))  # Reasonable latitudes
        assert np.all((17.0 < lons) & (lons < 19.0))  # Reasonable longitudes

    def test_round_trip_transformation(self, basic_surface):
        """Test that round-trip transformation is approximately reversible"""
        # Forward transformation
        intermediate = reproject_surface(basic_surface, "EPSG:3006", "EPSG:4326")

        # Backward transformation
        result = reproject_surface(intermediate, "EPSG:4326", "EPSG:3006")

        # Should be very close to original (within reasonable tolerance)
        np.testing.assert_allclose(result.vertices, basic_surface.vertices, rtol=1e-5, atol=1e-3)

    def test_round_trip_with_holes(self, surface_with_holes):
        """Test round-trip transformation preserves holes accurately"""
        # Forward transformation
        intermediate = reproject_surface(surface_with_holes, "EPSG:3006", "EPSG:4326")

        # Backward transformation
        result = reproject_surface(intermediate, "EPSG:4326", "EPSG:3006")

        # Vertices should match
        np.testing.assert_allclose(result.vertices, surface_with_holes.vertices, rtol=1e-5, atol=1e-3)

        # Holes should match
        assert len(result.holes) == len(surface_with_holes.holes)
        for i in range(len(result.holes)):
            np.testing.assert_allclose(result.holes[i], surface_with_holes.holes[i], rtol=1e-5, atol=1e-3)


class TestDifferentCRSTypes:
    """Test reprojection between different types of CRS"""

    def test_projected_to_geographic(self, basic_surface):
        """Test reprojection from projected (EPSG:3006) to geographic (EPSG:4326)"""
        result = reproject_surface(basic_surface, "EPSG:3006", "EPSG:4326")
        # Geographic coordinates should be in reasonable ranges
        assert np.all((-180 <= result.vertices[:, 0]) & (result.vertices[:, 0] <= 180))
        assert np.all((-90 <= result.vertices[:, 1]) & (result.vertices[:, 1] <= 90))

    def test_projected_to_projected(self, basic_surface):
        """Test reprojection between projected CRS (EPSG:3006 to EPSG:32633 - UTM Zone 33N)"""
        result = reproject_surface(basic_surface, "EPSG:3006", "EPSG:32633")
        # Coordinates should change (even if slightly) but remain in projected ranges
        # Check that arrays are not exactly equal
        assert not np.array_equal(result.vertices[:, :2], basic_surface.vertices[:, :2])
        assert result.vertices.shape == basic_surface.vertices.shape

    def test_utm_zones(self):
        """Test reprojection between different UTM zones"""
        # Create surface in UTM Zone 33N
        vertices = np.array([
            [500000.0, 6500000.0, 100.0],
            [500100.0, 6500000.0, 100.0],
            [500100.0, 6500100.0, 100.0],
            [500000.0, 6500100.0, 100.0],
        ])
        surface = Surface(vertices=vertices)

        # Reproject to UTM Zone 32N
        result = reproject_surface(surface, "EPSG:32633", "EPSG:32632")
        assert result.vertices.shape == vertices.shape
        assert not np.allclose(result.vertices[:, :2], vertices[:, :2])


class TestNormalPreservation:
    """Test that normal vectors are handled correctly"""

    def test_normal_is_not_modified(self, basic_surface):
        """Normal vector should not be automatically updated during reprojection"""
        # Set a normal vector
        basic_surface.normal = np.array([0.0, 0.0, 1.0])

        result = reproject_surface(basic_surface, "EPSG:3006", "EPSG:4326")

        # The normal vector is not reprojected (it's a direction, not a position)
        # In the current implementation, it's simply not copied
        # This test just ensures no errors occur
        assert result is not None


class TestLargeSurfaces:
    """Test with large surfaces"""

    def test_large_surface_with_many_vertices(self):
        """Test reprojection with a surface having many vertices"""
        # Create a surface with 1000 vertices
        n = 1000
        t = np.linspace(0, 2 * np.pi, n)
        x = 500000 + 100 * np.cos(t)
        y = 6500000 + 100 * np.sin(t)
        z = np.full(n, 100.0)
        vertices = np.column_stack((x, y, z))

        surface = Surface(vertices=vertices)
        result = reproject_surface(surface, "EPSG:3006", "EPSG:4326")

        assert result.vertices.shape == (n, 3)
        assert not np.allclose(result.vertices[:, :2], vertices[:, :2])

    def test_surface_with_many_holes(self):
        """Test surface with many holes"""
        vertices = np.array([
            [500000.0, 6500000.0, 100.0],
            [501000.0, 6500000.0, 100.0],
            [501000.0, 6501000.0, 100.0],
            [500000.0, 6501000.0, 100.0],
        ])

        # Create 10 small holes
        holes = []
        for i in range(10):
            hole = np.array([
                [500100.0 + i * 80, 6500100.0, 100.0],
                [500150.0 + i * 80, 6500100.0, 100.0],
                [500150.0 + i * 80, 6500150.0, 100.0],
                [500100.0 + i * 80, 6500150.0, 100.0],
            ])
            holes.append(hole)

        surface = Surface(vertices=vertices, holes=holes)
        result = reproject_surface(surface, "EPSG:3006", "EPSG:4326")

        assert len(result.holes) == 10
        for hole in result.holes:
            assert hole.shape == (4, 3)


if __name__ == "__main__":
    pytest.main()
