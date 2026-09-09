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
from .. import dtcc_pb2 as proto


def _validate_integer_pixels(values, dtype):
    if dtype.kind not in "biu" or values.size == 0:
        return
    if not np.isfinite(values).all() or np.any(values != np.floor(values)):
        raise ValueError("Raster protobuf integer pixels must be finite integers.")
    # Exclusive power-of-two upper bounds remain exact as floats, unlike int64.max.
    if dtype.kind == "b":
        lower, upper = 0, 2
    else:
        bits = np.iinfo(dtype).bits
        upper = 2 ** (bits - 1) if dtype.kind == "i" else 2 ** bits
        lower = -upper if dtype.kind == "i" else 0
    if np.any(values < lower) or np.any(values >= upper):
        raise ValueError(f"Raster protobuf pixels are outside the range of {dtype}.")

# FIXME: Make Raster fit the UML diagram
# FIXME: Make Raster own a Grid that holds Transform and Bounds


@dataclass
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

    def __str__(self):
        """
        Return a string representation of the Raster.

        Returns
        -------
        str
            A string representation of the Raster.

        """
        return f"DTCC Raster with {self.data.shape} values"

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
            for x, y in ((0, 0), (self.width, 0), (0, self.height), (self.width, self.height))
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

    def to_proto(self) -> proto.Raster:
        """
        Convert the Raster object to a protobuf representation.

        Returns
        -------
        proto.Raster
            A protobuf representation of the Raster.

        """
        if not isinstance(self.data, np.ndarray) or self.data.dtype.kind not in "biuf":
            raise ValueError("Raster.data must be a real numeric NumPy array.")
        if self.data.ndim not in (0, 2, 3):
            raise ValueError("Raster.data must have shape (H, W) or (H, W, C).")
        if self.data.ndim == 0 and not np.isnan(self.data):
            raise ValueError(
                "Raster.data scalar values are unsupported; use shape (1, 1) "
                "for one pixel or Raster() for an empty raster."
            )
        if self.data.ndim == 3 and self.channels < 1:
            raise ValueError("Raster.data must have at least one channel.")
        if not isinstance(self.georef, Affine) or not np.isfinite(self.georef).all():
            raise ValueError("Raster.georef must be a finite Affine transform.")
        if self.data.dtype.kind in "biu":
            encoded = self.data.astype(np.float32)
            _validate_integer_pixels(encoded, self.data.dtype)
            if not np.array_equal(encoded.astype(self.data.dtype), self.data):
                raise ValueError(
                    "Raster protobuf float32 values cannot preserve these integer pixels; "
                    "use a lossless raster format such as GeoTIFF."
                )

        pb = proto.Raster()
        pb.height = self.height
        pb.width = self.width
        pb.channels = self.channels
        # The scalar sentinel has no pixels and must never emit a value.
        if self.data.ndim != 0:
            pb.values.extend(self.data.flatten())
        pb.nodata = self.nodata
        pb.dtype = self.data.dtype.name
        pb.georef.extend(
            [
                self.georef.a,
                self.georef.b,
                self.georef.c,
                self.georef.d,
                self.georef.e,
                self.georef.f,
            ]
        )
        pb.crs = self.crs

        return pb

    def from_proto(self, pb: Union[proto.Raster, bytes]):
        """
        Initialize the Raster object from a protobuf representation.

        Parameters
        ----------
        pb : Union[proto.Raster, bytes]
            A protobuf representation of the Raster or a bytes object.

        Returns
        -------
        None

        """
        if isinstance(pb, bytes):
            pb = proto.Raster.FromString(pb)

        height, width, channels = pb.height, pb.width, pb.channels
        # Older payloads stored dimensions only in grid and omitted channels.
        if height == 0 and width == 0 and pb.HasField("grid"):
            height, width = pb.grid.height, pb.grid.width
        if min(height, width, channels) < 0:
            raise ValueError("Raster protobuf dimensions must be nonnegative.")
        if channels == 0 and (height or width):
            channels = 1
        # The previous empty-raster writer emitted one uninitialized scalar,
        # despite all dimensions being zero. It never represented a pixel.
        legacy_empty = height == width == channels == 0 and len(pb.values) == 1
        if not legacy_empty and len(pb.values) != height * width * channels:
            raise ValueError(
                "Raster protobuf value count does not match height * width * channels."
            )
        if len(pb.georef) not in (0, 6):
            raise ValueError("Raster protobuf georef must contain exactly six coefficients.")
        # An absent georef and dtype mean identity and float64 in legacy files.
        georef = Affine.identity() if not pb.georef else Affine(*pb.georef)
        if not np.isfinite(georef).all():
            raise ValueError("Raster protobuf georef coefficients must be finite.")
        try:
            dtype = np.dtype(pb.dtype) if pb.dtype else np.dtype(float)
        except TypeError as exc:
            raise ValueError(f"Invalid Raster protobuf dtype: {pb.dtype!r}.") from exc
        if dtype.kind not in "biuf":
            raise ValueError("Raster protobuf dtype must be real numeric.")

        values = np.empty(0) if legacy_empty else np.array(pb.values)
        _validate_integer_pixels(values, dtype)
        if channels == 0:
            data = np.array(np.nan)
        elif channels == 1:
            data = values.astype(dtype).reshape((height, width))
        else:
            data = values.astype(dtype).reshape((height, width, channels))
        self.data = data
        self.nodata = pb.nodata
        self.georef = georef
        self.crs = pb.crs
