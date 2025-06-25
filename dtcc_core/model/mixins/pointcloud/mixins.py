from dtcc_core.logging import info, warning, error

from typing import TypeVar, TYPE_CHECKING

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
