import dtcc_core
from dtcc_core import datasets
from dtcc_core.builder import tree_raster_from_pointcloud
from dtcc_core.model import PointCloud, Building, Raster
from pydantic import BaseModel, Field
from typing import Optional, Literal

from dtcc_core.datasets import DatasetDescriptor, DatasetBaseArgs


class TreeArgs(DatasetBaseArgs):
    min_height: float = Field(
        2.0, description="Minimum height of trees to include (in meters)"
    )
    cell_size: float = Field(
        0.5, description="Cell size of the output tree raster (in meters)"
    )
    smallest_cluster: float = Field(
        4, description="Minimum size of tree clusters to retain (in m^2)"
    )
    fill_hole_size: float = Field(
        1.5, description="Maximum size of holes to fill in the raster (in m^2)"
    )

    sigma: float = Field(
        0.5, description="Standard deviation for Gaussian filter to smooth the raster"
    )

    format: Optional[Literal["tif"]] = Field("tif", description="Output file format")


class TreeDataset(DatasetDescriptor):
    name = "trees"
    description = "Generate a raster of tree heights from point cloud data."
    ArgsModel = TreeArgs

    def build(self, args: TreeArgs):
        bounds = self.parse_bounds(args.bounds)
        pc: PointCloud = dtcc_core.io.data.download_pointcloud(bounds=bounds)

        tree_raster: Raster = tree_raster_from_pointcloud(
            pc,
            cell_size=args.cell_size,
            shortest_tree=args.min_height,
            smallest_cluster=args.smallest_cluster,
            fill_hole_size=args.fill_hole_size,
            sigma=args.sigma,
        )

        if args.format is not None:
            return self.export_to_bytes(tree_raster, args.format)
        return tree_raster
