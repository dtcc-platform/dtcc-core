from dtcc_core.logging import info, warning, error

from typing import TypeVar, TYPE_CHECKING, Union, List

from ....model.geometry import Bounds


if TYPE_CHECKING:
    from ....model.geometry import PointCloud
    from ....model.values import Raster

    T_Pointcloud = TypeVar("T_Pointcloud", bound=PointCloud)


class PointCloudBuilderMixin:

    def rasterize(
        self: "T_Pointcloud",
        cell_size: float,
        bounds: Bounds = None,
        window_size: int = 3,
        radius: float = 0,
        ground_only: bool = True,
        fill_holes: bool = True,
    ) -> "Raster":
        """
        Rasterize the point cloud into a raster.

        Args:
            bounds (Bounds): The bounds to rasterize within.
            cell_size (float): The size of the cells in the raster.
            raster_type (str): The type of raster to create (e.g., 'elevation').
            smoothing (int): The number of smoothing iterations to apply.

        Returns:
            Raster: The resulting raster object.
        """
        from dtcc_core.builder.pointcloud.convert import rasterize

        return rasterize(
            self,
            cell_size=cell_size,
            bounds=bounds,
            window_size=window_size,
            radius=radius,
            ground_only=ground_only,
            fill_holes=fill_holes,
        )


class PointcloudFilterMixin:
    """
    Mixin for filtering point clouds based on various criteria.
    """

    def remove_global_outliers(
        self: "T_Pointcloud", margin: float = 3.0
    ) -> "T_Pointcloud":
        """
        Remove outliers from the point cloud based on Z-value deviations. This method is non-mutating
        (does not modify data in-place).

        Args:
            margin (float): The margin in standard deviations to consider a point an outlier.

        Returns:
            T_Pointcloud: A new pointcloud with the outliers removed.
        """
        from dtcc_core.builder.pointcloud.filter import (
            find_global_outliers,
            remove_global_outliers,
        )

        return remove_global_outliers(self, margin)

    def statistical_outlier_filter(
        self: "T_Pointcloud", neighbours: int = 5, outlier_margin: float = 1.0
    ) -> "T_Pointcloud":
        """
        Remove statistical outliers from the point cloud. This method is non-mutating
        (does not modify data in-place).

        Args:
            neighbours (int): The number of neighbours to consider for the outlier detection.
            outlier_margin (float): The margin in standard deviations to consider a point an outlier.

        Returns:
            T_Pointcloud: A new point cloud with the outliers removed.
        """
        from dtcc_core.builder.pointcloud.filter import (
            find_statistical_outliers,
            statistical_outlier_filter,
        )

        return statistical_outlier_filter(self, neighbours, outlier_margin)

    def classification_filter(
        self: "T_Pointcloud", classes: Union[int, List[int]], keep: bool = False
    ) -> "T_Pointcloud":
        """
        Filter the point cloud based on classification. This method is non-mutating
        (does not modify data in-place).

        Args:
            classes (Union[int, List[int]]): The classification(s) to filter by.
            keep (bool): If True, keep points with the specified classification(s); otherwise, remove them.

        Returns:
            T_Pointcloud: A new point cloud with the specified classification(s) filtered out or kept.
        """
        from dtcc_core.builder.pointcloud.filter import (
            find_classification,
            classification_filter,
        )

        return classification_filter(self, classes, keep)

    def crop(self: "T_Pointcloud", bounds: Bounds) -> "T_Pointcloud":
        """
        Crop the point cloud to the specified bounds. This method is non-mutating
        (does not modify data in-place).

        Args:
            bounds (Bounds): The bounds to crop the point cloud to.

        Returns:
            T_Pointcloud: A new point cloud cropped to the specified bounds.
        """
        from dtcc_core.builder.pointcloud.filter import pts_in_bounds, crop

        return crop(self, bounds)

    def z_filter(self: "T_Pointcloud", zmin, zmax) -> "T_Pointcloud":
        """
        Filter the point cloud based on Z-values. This method is non-mutating
        (does not modify data in-place).

        Args:
            zmin (float): The minimum Z-value to keep.
            zmax (float): The maximum Z-value to keep.

        Returns:
            Pointcloud: A new point cloud with points outside the specified Z-range removed.
        """
        from dtcc_core.builder.pointcloud.filter import z_range_filter

        return z_range_filter(self, zmin, zmax)

    def remove_vegetation(self: "T_Pointcloud") -> "T_Pointcloud":
        """
        Remove vegetation points from a point cloud. This method is non-mutating
        (does not modify data in-place).

        Returns:
            T_Pointcloud: A new point cloud with vegetation points removed.
        """
        from dtcc_core.builder.pointcloud.filter import remove_vegetation

        return remove_vegetation(self)

    def get_vegetation(self: "T_Pointcloud") -> "T_Pointcloud":
        """
        Get the vegetation points from the point cloud and return as a new pointcloud. This method is non-mutating
        (does not modify data in-place).

        Returns:
            T_Pointcloud: A new point cloud containing only the vegetation points.
        """
        from dtcc_core.builder.pointcloud.filter import get_vegetation

        return get_vegetation(self)
