# Copyright(C) 2023 Anders Logg
# Licensed under the MIT License

from dataclasses import dataclass, field
from typing import Union
import numpy as np


from .geometry import Geometry
from .bounds import Bounds
from .. import dtcc_pb2 as proto


@dataclass
class PointCloud(Geometry):
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

    def __str__(self):
        """
        Return a string representation of the PointCloud, containing its boundaries
        and number of points.

        Returns
        -------
        str
            A string representation of the PointCloud.

        """
        return f"DTCC PointCloud on {self.bounds} with {len(self.points)} points"

    def __repr__(self):
        result = f"DTCC PointCloud on {self.bounds} with {len(self.points)} points"
        return result

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
        Remove points from the point cloud using the given indices.

        Parameters
        ----------
        indices : np.ndarray
            An array of indices specifying which points to remove.

        Returns
        -------
        None

        """
        if len(indices) == 0:
            return self
        self.points = np.delete(self.points, indices, axis=0)
        if len(self.classification) > 0:
            self.classification = np.delete(self.classification, indices, axis=0)
        if len(self.intensity) > 0:
            self.intensity = np.delete(self.intensity, indices, axis=0)
        if len(self.return_number) > 0:
            self.return_number = np.delete(self.return_number, indices, axis=0)
        if len(self.num_returns) > 0:
            self.num_returns = np.delete(self.num_returns, indices, axis=0)
        self.calculate_bounds()
        return self

    def keep_points(self, indices: np.ndarray) -> "PointCloud":
        """
        Keep only the points specified by the given indices.
        :param indices: indices of points to keep
        :return: PointCloud
        """

        removed_indices = np.setdiff1d(np.arange(len(self.points)), indices)
        return self.remove_points(removed_indices)

    def merge(self, other):
        """
        Merge another point cloud into this point cloud.

        Parameters
        ----------
        other : PointCloud
            Another PointCloud object to merge into this one.

        Returns
        -------
        None

        """

        if len(other.points) == 0:
            return self

        if len(self.points) == 0:
            self.points = other.points
        else:
            self.points = np.concatenate((self.points, other.points))

        if len(other.classification) == len(other.points):
            self.classification = np.concatenate(
                (self.classification, other.classification)
            )
        if len(other.intensity) == len(other.points):
            self.intensity = np.concatenate((self.intensity, other.intensity))
        if len(other.return_number) == len(other.points):
            self.return_number = np.concatenate(
                (self.return_number, other.return_number)
            )
        if len(other.num_returns) == len(other.points):
            self.num_returns = np.concatenate((self.num_returns, other.num_returns))
        self.calculate_bounds()
        return self

    def offset(self, offset: Union[list, np.ndarray]) -> "PointCloud":
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

    def to_proto(self) -> proto.Geometry:
        """Return a protobuf representation of the PointCloud.

        Returns
        -------
        proto.Geometry
            A protobuf representation of the PointCloud and as a Geometry.
        """

        # Handle Geometry fields
        pb = Geometry.to_proto(self)

        # Handle specific fields
        _pb = proto.PointCloud()
        _pb.points.extend(self.points.flatten())
        _pb.classification.extend(self.classification)
        _pb.intensity.extend(self.intensity)
        _pb.return_number.extend(self.return_number)
        _pb.num_returns.extend(self.num_returns)
        pb.point_cloud.CopyFrom(_pb)

        return pb

    def from_proto(self, pb: Union[proto.Geometry, bytes]):
        """Initialize PointCloud from a protobuf representation.

        Parameters
        ----------
        pb: Union[proto.Geometry, bytes]
            The protobuf message or its serialized bytes representation.
        """

        # Handle byte representation
        if isinstance(pb, bytes):
            pb = proto.Geometry.FromString(pb)

        # Handle Geometry fields
        Geometry.from_proto(self, pb)

        # Handle specific fields
        _pb = pb.point_cloud
        self.points = np.array(_pb.points).reshape(-1, 3)
        self.classification = np.array(_pb.classification).astype(np.uint8)
        self.intensity = np.array(_pb.intensity).astype(np.uint16)
        self.return_number = np.array(_pb.return_number).astype(np.uint8)
        self.num_returns = np.array(_pb.num_returns).astype(np.uint8)
