# Copyright(C) 2023 Dag Wästberg
# Licensed under the MIT License
from email.headerregistry import Address

import numpy as np
from typing import Union
from dataclasses import dataclass, field
from affine import Affine
from copy import deepcopy

from ..geometry.bounds import Bounds
from ..logging import info, warning, error
from ..model import Model


# FIXME: Make Raster fit the UML diagram
# FIXME: Make Raster own a Grid that holds Transform and Bounds


@dataclass(repr=False)
class Raster(Model):
    """
    A georeferenced n-dimensional raster of values.

    This class represents a georeferenced n-dimensional raster of values, where `data` is a
    NumPy array of shape (height, width, channels) or (height, width) if `channels` is 1.

    Attributes
    ----------
    data : np.ndarray
        The data of the raster as a NumPy array.
    georef : Affine
        The georeference of the raster.
    nodata : float
        The value representing nodata or missing data.
    crs : str
        The coordinate reference system (CRS) information.

    """

    # Shape () denotes an empty raster; initialize its sentinel deterministically.
    data: np.ndarray = field(default_factory=lambda: np.array(np.nan))
    georef: Affine = field(default_factory=Affine.identity)
    nodata: float = np.nan
    crs: str = ""

    def _info_sections(self):
        sections = super()._info_sections()
        sections[0][2].extend(
            [("No data", self.nodata), ("Georeference", str(self.georef))]
        )
        if self.data.ndim >= 2:
            sections[0][2].append(("Bounds", self.bounds.bndstr))
        return sections

    def _summary_items(self):
        return [
            ("shape", self.data.shape),
            ("dtype", str(self.data.dtype)),
            ("crs", self.crs),
        ]

    @property
    def shape(self):
        """
        Get the shape of the raster data.

        Returns
        -------
        tuple
            A tuple representing the shape of the raster data.

        """
        return self.data.shape

    @property
    def height(self):
        """
        Get the height (number of rows) of the raster data.

        Returns
        -------
        int
            The height of the raster data.

        """
        if len(self.data.shape) < 2:
            return 0
        return self.data.shape[0]

    @property
    def width(self):
        """
        Get the width (number of columns) of the raster data.

        Returns
        -------
        int
            The width of the raster data.

        """
        if len(self.data.shape) < 2:
            return 0
        return self.data.shape[1]

    @property
    def channels(self):
        """
        Get the number of channels in the raster data.

        Returns
        -------
        int
            The number of channels in the raster data.

        """
        if len(self.data.shape) < 2:
            return 0
        if len(self.data.shape) == 2:
            return 1
        else:
            return self.data.shape[2]

    @property
    def bounds(self):
        """
        Get the spatial bounds of the raster.

        Returns
        -------
        Bounds
            The spatial bounds of the raster.

        """

        corners = [
            self.georef * (x, y)
            for x, y in (
                (0, 0),
                (self.width, 0),
                (0, self.height),
                (self.width, self.height),
            )
        ]
        x, y = zip(*corners)
        return Bounds(min(x), min(y), max(x), max(y), 0, 0)

    def set_bounds(self, bounds: Bounds):
        """
        Set the spatial bounds of the raster.

        Parameters
        ----------
        bounds : Bounds
            The spatial bounds to set.

        Returns
        -------
        None

        """
        self.georef = Affine(
            bounds.width / self.width,
            0,
            bounds.xmin,
            0,
            -bounds.height / self.height,
            bounds.ymax,
        )

    def calculate_bounds(self):
        """
        Calculate the spatial bounds of the raster.

        Returns
        -------
        Bounds
            The spatial bounds of the raster.

        """
        return self.bounds

    @property
    def cell_size(self):
        """
        Get the cell size (pixel size) of the raster.

        Returns
        -------
        tuple
            A tuple containing the horizontal and vertical cell sizes.

        """
        return (self.georef.a, self.georef.e)

    @property
    def min(self):
        """
        Get the minimum value of the raster.

        Returns
        -------
        float
            The minimum value of the raster.

        """
        return self.data.min()

    @property
    def max(self):
        """
        Get the maximum value of the raster.

        Returns
        -------
        float
            The maximum value of the raster.

        """
        return self.data.max()

    def pixel_to_georef(self, x: float, y: float):
        """get the georeferenced coordinate of a given pixel"""

        return self.georef * (y, x)

    def get_value(self, x: float, y: float):
        """
        Get the value at the given coordinate.

        Parameters
        ----------
        x : float
            The x-coordinate.
        y : float
            The y-coordinate.

        Returns
        -------
        data
            The value at the given coordinate.

        """
        col, row = ~self.georef * (x, y)
        try:
            data = self.data[int(row), int(col)]
        except IndexError:
            error_str = f"IndexError in get_value at ({x}, {y})"
            error_str += f"\ncol: {col}, row: {row}"
            error_str += f"\ngeoref: {self.georef}"
            error_str += f"\nshape: {self.data.shape}"
            error(error_str)
            raise
            # data = self.nodata
        return data

    def copy(self, no_data=False):
        """
        Create a copy of the raster.

        Parameters
        ----------
        no_data : bool, default False
            When ``True``, copy only metadata (georef, nodata, crs) and omit data array.

        Returns
        -------
        Raster
            Deep copy of the raster or metadata-only copy when ``no_data`` is set.
        """
        if not no_data:
            return deepcopy(self)
        else:
            copy_raster = Raster()
            copy_raster.georef = self.georef
            copy_raster.nodata = self.nodata
            copy_raster.crs = self.crs
            return copy_raster
