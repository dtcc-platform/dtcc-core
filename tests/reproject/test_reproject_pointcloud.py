import pytest
import numpy as np
from unittest.mock import Mock, patch, MagicMock
from dtcc_core.reproject.reproject import reproject_pointcloud
from dtcc_core.model import PointCloud
from dtcc_core.model.geometry.transform import Transform


@pytest.fixture
def sample_points():
    """Sample point cloud data in SWEREF99 TM (EPSG:3006)"""
    return np.array(
        [
            [500000.0, 6500000.0, 100.0],
            [500100.0, 6500100.0, 105.0],
            [500200.0, 6500200.0, 110.0],
        ]
    )


@pytest.fixture
def sample_classification():
    """Sample classification data"""
    return np.array([2, 2, 6])


@pytest.fixture
def sample_intensity():
    """Sample intensity data"""
    return np.array([100, 150, 200])


@pytest.fixture
def basic_pointcloud(sample_points, sample_classification, sample_intensity):
    """Create a basic PointCloud with identity transform"""
    pc = PointCloud(
        points=sample_points,
        classification=sample_classification,
        intensity=sample_intensity,
    )
    pc.transform = Transform()
    pc.transform.srs = "EPSG:3006"
    return pc


@pytest.fixture
def pointcloud_no_crs(sample_points, sample_classification, sample_intensity):
    """Create a PointCloud without CRS information"""
    pc = PointCloud(
        points=sample_points,
        classification=sample_classification,
        intensity=sample_intensity,
    )
    pc.transform = Transform()
    pc.transform.srs = None
    return pc


@pytest.fixture
def pointcloud_with_transform(sample_points, sample_classification, sample_intensity):
    """Create a PointCloud with non-identity transform"""
    pc = PointCloud(
        points=sample_points,
        classification=sample_classification,
        intensity=sample_intensity,
    )
    transform = Transform()
    transform.srs = "EPSG:3006"
    # Set up a translation transform
    transform.translation = np.array([100.0, 200.0, 0.0])
    pc.transform = transform
    return pc


class TestReprojectionBasics:
    """Test basic reprojection functionality"""

    def test_same_crs_returns_original(self, basic_pointcloud):
        """When source and target CRS are the same, return original pointcloud"""
        result = reproject_pointcloud(basic_pointcloud, "EPSG:3006", "EPSG:3006")
        assert result is basic_pointcloud

    def test_missing_target_crs_raises_error(self, basic_pointcloud):
        """Missing target CRS should raise ValueError"""
        with pytest.raises(ValueError, match="Target CRS not provided"):
            reproject_pointcloud(basic_pointcloud, "EPSG:3006", None)

    def test_empty_target_crs_raises_error(self, basic_pointcloud):
        """Empty target CRS should raise ValueError"""
        with pytest.raises(ValueError, match="Target CRS not provided"):
            reproject_pointcloud(basic_pointcloud, "EPSG:3006", "")

    def test_basic_reprojection_shape_preserved(self, basic_pointcloud):
        """Reprojection should preserve shape of point cloud"""
        result = reproject_pointcloud(basic_pointcloud, "EPSG:3006", "EPSG:4326")
        assert result.points.shape == basic_pointcloud.points.shape

    def test_reprojection_preserves_classification(self, basic_pointcloud):
        """Classification data should be preserved"""
        result = reproject_pointcloud(basic_pointcloud, "EPSG:3006", "EPSG:4326")
        np.testing.assert_array_equal(
            result.classification, basic_pointcloud.classification
        )

    def test_reprojection_preserves_intensity(self, basic_pointcloud):
        """Intensity data should be preserved"""
        result = reproject_pointcloud(basic_pointcloud, "EPSG:3006", "EPSG:4326")
        np.testing.assert_array_equal(result.intensity, basic_pointcloud.intensity)

    def test_reprojection_updates_target_crs(self, basic_pointcloud):
        """Target CRS should be set on output"""
        result = reproject_pointcloud(basic_pointcloud, "EPSG:3006", "EPSG:4326")
        assert result.transform.srs == "EPSG:4326"

    def test_reprojection_coordinates_change(self, basic_pointcloud):
        """Coordinates should change after reprojection"""
        result = reproject_pointcloud(basic_pointcloud, "EPSG:3006", "EPSG:4326")
        # Coordinates should be different (SWEREF99 vs WGS84)
        assert not np.allclose(result.points[:, :2], basic_pointcloud.points[:, :2])

    def test_z_coordinates_preserved(self, basic_pointcloud):
        """Z coordinates should be preserved during reprojection"""
        result = reproject_pointcloud(basic_pointcloud, "EPSG:3006", "EPSG:4326")
        np.testing.assert_array_equal(
            result.points[:, 2], basic_pointcloud.points[:, 2]
        )


class TestCRSDetection:
    """Test CRS detection and default behavior"""

    def test_no_src_crs_uses_pointcloud_crs(self, basic_pointcloud):
        """When src_crs is None, use pointcloud's CRS"""
        result = reproject_pointcloud(basic_pointcloud, None, "EPSG:4326")
        # Should successfully reproject using the pointcloud's EPSG:3006
        assert result.transform.srs == "EPSG:4326"

    def test_empty_src_crs_uses_pointcloud_crs(self, basic_pointcloud):
        """When src_crs is empty string, use pointcloud's CRS"""
        result = reproject_pointcloud(basic_pointcloud, "", "EPSG:4326")
        assert result.transform.srs == "EPSG:4326"

    @patch("dtcc_core.reproject.reproject.warning")
    def test_no_crs_anywhere_defaults_to_3006(self, mock_warning, pointcloud_no_crs):
        """When no CRS is available, default to EPSG:3006 with warning"""
        result = reproject_pointcloud(pointcloud_no_crs, None, "EPSG:4326")
        mock_warning.assert_called_once_with(
            "Source CRS not provided, defaulting to EPSG:3006"
        )
        assert result.transform.srs == "EPSG:4326"

    @patch("dtcc_core.reproject.reproject.warning")
    def test_override_geometry_crs_flag(self, mock_warning, basic_pointcloud):
        """Test override_geometry_crs flag behavior"""
        # When override is True, provided src_crs should be used
        result = reproject_pointcloud(
            basic_pointcloud,
            "EPSG:4326",  # Different from pointcloud's EPSG:3006
            "EPSG:32633",
            override_geometry_crs=True,
        )
        assert result.transform.srs == "EPSG:32633"
        # Should not warn when override is explicit
        assert not mock_warning.called

    def test_provided_src_crs_uses_geometry_crs_by_default(self, basic_pointcloud):
        """When src_crs is provided but override is False, use geometry CRS"""
        result = reproject_pointcloud(
            basic_pointcloud,
            "EPSG:4326",  # This should be ignored
            "EPSG:32633",
            override_geometry_crs=False,
        )
        # Should use the pointcloud's EPSG:3006, not the provided EPSG:4326
        assert result.transform.srs == "EPSG:32633"


class TestTransformHandling:
    """Test handling of pointcloud transforms"""

    def test_identity_transform_uses_points_directly(self, basic_pointcloud):
        """With identity transform, use points directly"""
        result = reproject_pointcloud(basic_pointcloud, "EPSG:3006", "EPSG:4326")
        # Should work without errors
        assert result is not None


class TestEdgeCases:
    """Test edge cases and boundary conditions"""

    def test_single_point(self):
        """Test reprojection with single point"""
        pc = PointCloud(points=np.array([[500000.0, 6500000.0, 100.0]]))
        pc.transform = Transform()
        pc.transform.srs = "EPSG:3006"

        result = reproject_pointcloud(pc, "EPSG:3006", "EPSG:4326")
        assert result.points.shape == (1, 3)

    def test_large_pointcloud(self):
        """Test reprojection with large point cloud"""
        large_points = np.random.rand(10000, 3) * 1000 + np.array(
            [500000, 6500000, 100]
        )
        pc = PointCloud(points=large_points)
        pc.transform = Transform()
        pc.transform.srs = "EPSG:3006"

        result = reproject_pointcloud(pc, "EPSG:3006", "EPSG:4326")
        assert result.points.shape == large_points.shape

    def test_preserves_all_optional_attributes(self):
        """Test that all optional pointcloud attributes are preserved"""
        points = np.array([[500000.0, 6500000.0, 100.0]])
        pc = PointCloud(
            points=points,
            classification=np.array([2]),
            intensity=np.array([100]),
            return_number=np.array([1]),
            num_returns=np.array([1]),
        )
        pc.transform = Transform()
        pc.transform.srs = "EPSG:3006"

        result = reproject_pointcloud(pc, "EPSG:3006", "EPSG:4326")

        np.testing.assert_array_equal(result.classification, pc.classification)
        np.testing.assert_array_equal(result.intensity, pc.intensity)
        np.testing.assert_array_equal(result.return_number, pc.return_number)
        np.testing.assert_array_equal(result.num_returns, pc.num_returns)


class TestCoordinateAccuracy:
    """Test coordinate transformation accuracy"""

    def test_known_coordinate_transformation(self):
        """Test transformation with known coordinates"""
        # Known point in SWEREF99 TM (EPSG:3006): Stockholm roughly
        # X: 674032, Y: 6580822
        points = np.array([[674032.0, 6580822.0, 0.0]])
        pc = PointCloud(points=points)
        pc.transform = Transform()
        pc.transform.srs = "EPSG:3006"

        result = reproject_pointcloud(pc, "EPSG:3006", "EPSG:4326")

        # Should be roughly 59.33°N, 18.06°E (Stockholm)
        lat, lon = result.points[0, 1], result.points[0, 0]
        assert 59.0 < lat < 60.0  # Reasonable latitude
        assert 17.0 < lon < 19.0  # Reasonable longitude

    def test_round_trip_transformation(self):
        """Test that round-trip transformation is approximately reversible"""
        points = np.array([[674032.0, 6580822.0, 50.0]])
        pc = PointCloud(points=points)
        pc.transform = Transform()
        pc.transform.srs = "EPSG:3006"

        # Forward transformation
        intermediate = reproject_pointcloud(pc, "EPSG:3006", "EPSG:4326")

        # Backward transformation
        result = reproject_pointcloud(intermediate, "EPSG:4326", "EPSG:3006")

        # Should be very close to original (within reasonable tolerance)
        np.testing.assert_allclose(result.points, points, rtol=1e-5, atol=1e-3)


if __name__ == "__main__":
    pytest.main()
