# Copyright(C) 2025 Dag Wästberg
# Licensed under the MIT License

from __future__ import annotations


from pathlib import Path
from typing import Union, List

from dtcc_core.model.model import Model

from dtcc_core.model.geometry import Bounds, PointCloud
from ..pointcloud import load as load_pointcloud

from enum import Enum, auto
from shapely.geometry import box
from shapely.strtree import STRtree


class PointCloudContainerType(Enum):
    FILESYSTEM = auto()


class PointCloudDirectory(Model):

    def __init__(self, files: List[Path] = None, bounds: List[Bounds] = None):
        self.container_type = PointCloudContainerType.FILESYSTEM
        self.file_list = files if files is not None else []
        self.bounds_list = bounds if bounds is not None else []
        self._bounds = None
        for b in self.bounds_list:
            if self._bounds is None:
                self._bounds = b
            else:
                self._bounds.union(b)
        self._rtree = STRtree(
            [box(b.xmin, b.ymin, b.xmax, b.ymax) for b in self.bounds_list]
        )

    @property
    def bounds(self) -> Bounds:
        return self._bounds

    def pointcloud(
        self, bounds: Bounds, points_only=False, points_classification_only=False
    ) -> PointCloud:

        needed_files = self._rtree.query(
            box(bounds.xmin, bounds.ymin, bounds.xmax, bounds.ymax)
        )
        files_to_load = [self.file_list[i] for i in needed_files]
        pc = load_pointcloud(
            files_to_load,
            bounds=bounds,
            points_only=points_only,
            points_classification_only=points_classification_only,
        )
        return pc

    def __str__(self):
        return f"PointCloudDirectory with {len(self.file_list)} files"

    def __len__(self):
        return len(self.file_list)

    def to_proto(self):
        pass

    def from_proto(self, pb):
        pass
