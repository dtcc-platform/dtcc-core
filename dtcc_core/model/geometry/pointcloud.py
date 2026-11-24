# Copyright(C) 2023 Anders Logg
# Licensed under the MIT License

from dataclasses import dataclass, field

import numpy as np

from ..mixins.pointcloud.mixins import PointCloudBuilderMixin, PointcloudFilterMixin
from .bounds import Bounds
from .geometry import Geometry


@dataclass(repr=False)
class PointCloud(PointCloudBuilderMixin, PointcloudFilterMixin, Geometry):
    """Represents a set of points in 3D.

    Attributes
    ----------
    points : np.ndarray
        The points of the point cloud as (n,3) dimensional numpy array.
    classification : np.ndarray
        The classification of the points as (n,) dimensional numpy array.
    intensity : np.ndarray
        The intensity of the points as (n,) dimensional numpy array.
    return_number : np.ndarray
        The return number of the points as (n,) dimensional numpy array.
    num_returns : np.ndarray
        The number of returns of the points as (n,) dimensional numpy array.
    """

    points: np.ndarray = field(default_factory=lambda: np.empty((0, 3)))
    classification: np.ndarray = field(default_factory=lambda: np.empty(0))
    intensity: np.ndarray = field(default_factory=lambda: np.empty(0))
    return_number: np.ndarray = field(default_factory=lambda: np.empty(0))
    num_returns: np.ndarray = field(default_factory=lambda: np.empty(0))

    def _summary_items(self):
        return [("num_points", len(self.points))] + super()._summary_items()

    def __len__(self):
        """
        Get the number of points in the PointCloud.

        Returns
        -------
        int
            The number of points in the PointCloud.

        """
        return self.points.shape[0]

    def used_classifications(self) -> set:
        """
        Get the set of unique classifications used in the PointCloud.

        Returns
        -------
        set
            A set containing unique classifications.

        """
        return set(np.unique(self.classification))

    def calculate_bounds(self):
        """
        Calculate the bounds of the point cloud and update the bounds attribute.

        Returns
        -------
        None

        """
        """Calculate the bounds of the point cloud and update the bounds attribute."""
        if len(self.points) == 0:
            self._bounds = Bounds()
        else:
            self._bounds = Bounds(
                xmin=self.points[:, 0].min(),
                xmax=self.points[:, 0].max(),
                ymin=self.points[:, 1].min(),
                ymax=self.points[:, 1].max(),
                zmin=self.points[:, 2].min(),
                zmax=self.points[:, 2].max(),
            )
        return self._bounds

    def remove_points(self, indices: np.ndarray) -> "PointCloud":
        """
        Remove points and their vertex- or sample-associated time series values.

        Parameters
        ----------
        indices : np.ndarray
            An array of indices specifying which points to remove.

        Returns
        -------
        PointCloud
            This point cloud, with the points removed.

        """
        indices = np.array(indices).astype(int)
        if len(indices) == 0:
            return self
        keep = np.ones(len(self.points), dtype=bool)
        keep[indices] = False
        series_list = []
        for series in self.time_series:
            if series.association in ("vertex", "sample"):
                if series.num_elements != len(self.points):
                    raise ValueError(
                        f"TimeSeries {series.name!r} has {series.num_elements} "
                        f"elements for {len(self.points)} points."
                    )
                series = series._copy_with_values(series.values[:, keep])
            series_list.append(series)

        selected = {"points": self.points[keep]}
        for name in ("classification", "intensity", "return_number", "num_returns"):
            values = getattr(self, name)
            if len(values):
                selected[name] = np.asarray(values)[keep]
        for name, values in selected.items():
            setattr(self, name, values)
        self.time_series = series_list
        self.calculate_bounds()
        return self

    def keep_points(self, indices: np.ndarray) -> "PointCloud":
        """
        Keep only the points specified by the given indices.

        Parameters
        ----------
        indices : np.ndarray
            Indices of the points to keep.

        Returns
        -------
        PointCloud
            This point cloud, with all other points removed.
        """

        removed_indices = np.setdiff1d(np.arange(len(self.points)), indices)
        return self.remove_points(removed_indices)

    def merge(self, other):
        """
        Merge another point cloud into this point cloud.

        Time series require an explicit mapping and raise NotImplementedError.

        Parameters
        ----------
        other : PointCloud
            Another PointCloud object to merge into this one.

        Returns
        -------
        PointCloud
            This point cloud, with ``other`` merged in.

        """

        if other.time_series or (self.time_series and len(other.points)):
            raise NotImplementedError(
                "Merging point clouds cannot transfer time series"
            )
        if len(other.points) == 0:
            return self

        if len(self.points) == 0:
            self.points = other.points
        else:
            self.points = np.concatenate((self.points, other.points))

        for name in ("classification", "intensity", "return_number", "num_returns"):
            incoming = getattr(other, name)
            if len(incoming) == len(other.points):
                current = getattr(self, name)
                # Empty defaults are float arrays; do not promote LAS integer
                # attributes when the first populated array is merged.
                merged = (
                    incoming.copy()
                    if len(current) == 0
                    else np.concatenate((current, incoming))
                )
                setattr(self, name, merged)
        self.calculate_bounds()
        return self

    def offset(self, offset: list | np.ndarray) -> "PointCloud":
        """
        Offset the point cloud by the given offset.

        Parameters
        ----------
        offset : np.ndarray
            The offset to apply to the point cloud.

        Returns
        -------
        None

        """
        self.points += offset
        self.calculate_bounds()
        return self
