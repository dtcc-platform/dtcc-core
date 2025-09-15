from dtcc_core.logging import info, warning, error

from typing import TypeVar, TYPE_CHECKING


if TYPE_CHECKING:
    from ....model.object import City

    T_City = TypeVar("T_City", bound=City)


class CityModifyingMixin:

    def building_heights_from_pointcloud(
        self: "T_City",
        statistical_outlier_removal: bool = True,
        roof_outlier_neighbors=5,
        roof_outlier_margin=1.5,
        overwrite=True,
        keep_roof_points=False,
    ) -> "T_City":
        """
        Calculate building heights from point cloud and set as building attribute.

        Args:
            self (City): The city object to modify.
            statistical_outlier_removal : bool, default True
            Whether to apply statistical outlier removal to roof points.
        roof_outlier_neighbors : int, default 5
            Number of neighbors for outlier detection.
        roof_outlier_margin : float, default 1.5
            Margin for statistical outlier removal.
        overwrite : bool, default False
            Whether to overwrite existing height values.
        keep_roof_points : bool, default False
            Whether to keep extracted roof points as building geometry.

        Returns:
            City: The modified city object with building heights set as an attribute.
        """
        from dtcc_core.builder import (
            building_heights_from_pointcloud as _building_heights_from_pointcloud,
        )

        pc = self.pointcloud
        if pc is None or len(pc.points) == 0:
            raise ValueError("City has no point cloud geometry")

        if self.num_buildings == 0:
            raise ValueError("City has no buildings")

        terrain = self.terrain
        if terrain is None or terrain.raster is None:
            info("building city terrain raster from point cloud")
            self.build_terrain(cell_size=2.0, build_mesh=False)
        raster = self.terrain.raster
        if raster is None:
            raise ValueError("Failed to find or build terrain raster")

        buildings = self.buildings
        buildings_with_heights = _building_heights_from_pointcloud(
            buildings,
            pc,
            raster,
            statistical_outlier_remover=statistical_outlier_removal,
            roof_outlier_neighbors=roof_outlier_neighbors,
            roof_outlier_margin=roof_outlier_margin,
            overwrite=overwrite,
            keep_roof_points=keep_roof_points,
        )

        self.replace_buildings(buildings_with_heights)
        return self
