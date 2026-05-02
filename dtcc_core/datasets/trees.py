import dtcc_core
from dtcc_core import datasets
from dtcc_core.builder import tree_raster_from_pointcloud
from dtcc_core.model import PointCloud, Building, Raster, City
from dtcc_core.io.trees import save_trees
from pydantic import BaseModel, Field
from typing import Optional, Literal

from dtcc_core.datasets import DatasetDescriptor, DatasetBaseArgs

import tempfile


class TreeArgs(DatasetBaseArgs):
    tree_type: Literal["urban", "mixed", "dense", "arid"] = Field(
        "urban", description="Type of trees to detect"
    )
    cell_size: float = Field(
        0.5, description="Cell size of the output tree raster (in meters)"
    )

    vector_geometry: Literal["point", "circle"] = Field(
        "point",
        description="Save trees as points or circles when exporting to vector format",
    )

    format: Optional[Literal["tif", "gpkg", "geojson"]] = Field(
        None, description="Output file format"
    )


class TreesDataset(DatasetDescriptor):
    name = "trees"
    description = "Raster of tree heights or vector points representing trees from point cloud data."
    ArgsModel = TreeArgs
    data_category = "derived"
    result_kind = "tree_collection"
    python_return_type = "list[dtcc_core.model.Tree]"

    def build(self, args: TreeArgs):
        bounds = self.parse_bounds(args.bounds)
        pc: PointCloud = dtcc_core.io.data.download_pointcloud(bounds=bounds)
        if args.format == "tif":
            raster = tree_raster_from_pointcloud(
                pc,
                None,
                tree_type=args.tree_type,
                cell_size=args.cell_size,
            )
            return self.export_to_bytes(raster, "tif")
        else:
            city = City()
            city.add_point_cloud(pc)
            city.bounds = bounds
            trees = city.build_trees_from_pointcloud(tree_type=args.tree_type)
            if args.format is None:
                return trees
            elif args.format in ("gpkg", "geojson", "json"):
                if args.format in ("geojson", "json"):
                    args.format = "json"
                as_circles = args.vector_geometry == "circle"
                return self.export_to_bytes(
                    trees,
                    args.format,
                    save_callable=save_trees,
                    as_circles=as_circles,
                )
