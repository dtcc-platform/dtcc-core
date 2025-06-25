from dtcc_core.logging import info, warning, error


from pathlib import Path
from typing import Union
from typing import TypeVar, TYPE_CHECKING

from ....model.object import GeometryType
from ....model.geometry import PointCloud, Bounds
from ....model.values import Raster

if TYPE_CHECKING:
    from ....model.object import City

    T_City = TypeVar("T_City", bound=City)


class CityLoaderMixin:
    def load_footprints(
        self: "City", path: Union[str, Path], user_city_bounds=False
    ) -> "City":
        """
        Load building footprints from a shapefile or geopackage into a City object.

        Args:
            self (City): The City object to load the footprints into.
            path (str): The path to the shapefile containing the building footprints.
            user_city_bounds (bool, optional): only load footprints within the set bounds

        Returns:
            City: The updated City object with the loaded footprints.
        """

        # non-model imports must go here to avoid circular imports
        import dtcc_core.io as io

        if user_city_bounds:
            if self.bounds is None or self.bounds.area == 0:
                raise ValueError(
                    "City bounds must be set before loading footprints with user_city_bounds=True"
                )
            footprints = io.footprints.load(path, bounds=self.bounds)
        else:
            bounds = io.footprints.building_bounds(path, 2)
            footprints = io.footprints.load(path)
            if self.bounds is None:
                self.bounds = bounds
            else:
                self.bounds.union(bounds)
        self.add_buildings(footprints)
        return self

    def load_pointcloud(
        self: "City",
        path: Union[Path, str],
        remove_global_outliers: float = 3.0,
        user_city_bounds=False,
    ) -> "City":
        """

        Load pointcloud from las or csv file into city object
        """
        # non-model imports must go here to avoid circular imports
        import dtcc_core.io as io

        if user_city_bounds:
            if self.bounds is None or self.bounds.area == 0:
                raise ValueError(
                    "City bounds must be set before loading footprints with user_city_bounds=True"
                )
            pc = io.pointcloud.load(path, bounds=self.bounds)
        else:
            pc = io.load_pointcloud(path)
            bounds = pc.bounds
            if self.bounds is None:
                self.bounds = bounds
            else:
                self.bounds.union(bounds)
        if remove_global_outliers > 0:
            pc = pc.remove_global_outliers(remove_global_outliers)
        self.add_geometry(pc, GeometryType.POINT_CLOUD)

        return self

    def load_terrain_raster(self: "City", path: Union[Path, str]) -> "City":
        # non-model imports must go here to avoid circular imports
        import dtcc_core.io as io

        raster = io.load_raster(path)
        self.add_terrain(raster)
        return self


def _get_download_bounds(city_bounds, user_bounds):
    """
    Determine the bounds to use for downloading data.
    """

    if user_bounds is not None and not isinstance(user_bounds, Bounds):
        raise TypeError("bounds must be a Bounds object")

    if city_bounds is None:
        if user_bounds is None:
            error("No bounds found, either set city.bounds or pass in bounds")
        else:
            return user_bounds
    else:
        if user_bounds is not None:
            bounds = city_bounds.intersect(user_bounds)
        else:
            bounds = city_bounds
    return bounds


class CityDownloadMixin:
    def download_footprints(
        self: "City",
        bounds: Union[Bounds, None] = None,
    ) -> "City":
        """
        Download building footprints from a URL and load them into a City object.

        Args:
            self (City): The City object to load the footprints into.
            url (str): The URL to download the building footprints from.
            user_city_bounds (bool, optional): only load footprints within the set bounds

        Returns:
            City: The updated City object with the loaded footprints.
        """

        # non-model imports must go here to avoid circular imports
        import dtcc_core.io as io

        download_bounds = _get_download_bounds(self.bounds, bounds)

        footprints = io.data.download_footprints(bounds=download_bounds)
        self.add_buildings(footprints)
        self.calculate_bounds()
        return self

    def download_pointcloud(self: "City", bounds: Union[Bounds, None] = None) -> "City":
        """
        Download pointcloud from a URL and load it into a City object.

        Args:
            self (City): The City object to load the pointcloud into.
            bounds (Bounds, optional): The bounds to filter the pointcloud.

        Returns:
            City: The updated City object with the loaded pointcloud.
        """

        # non-model imports must go here to avoid circular imports
        import dtcc_core.io as io

        download_bounds = _get_download_bounds(self.bounds, bounds)

        pc = io.data.download_pointcloud(bounds=download_bounds)
        self.add_point_cloud(pc)
        self.calculate_bounds()
        return self
