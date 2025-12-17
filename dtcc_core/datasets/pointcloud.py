import os

import dtcc_core
from dtcc_core.model import PointCloud, Bounds
from typing import Literal, Optional, List, Tuple, Sequence, Union
from pydantic import BaseModel, Field
import tempfile

from io import BytesIO
from .dataset import DatasetDescriptor


class PointCloudArgs(BaseModel):
    bounds: Sequence[float] = Field(
        ...,
        description="Bounding box [minx, miny, [optional_minz], maxx, maxy, [optional_maxz]]",
    )
    classifications: Union[
        int, list[int], Literal["all", "terrain", "buildings", "vegetation"]
    ] = Field(
        None,
        description="Which classifications to include (e.g., [2, 5] for ground and vegetation)",
    )
    source: Literal["LM"] = Field("LM", description="Data source")
    format: Optional[Literal["copc", "las", "laz"]] = Field(
        None, description="Output file format"
    )
    remove_outliers: bool = Field(
        False, description="Whether to remove global outliers from the point cloud"
    )
    remove_outlier_threshold: float = Field(
        3.0, description="Threshold for outlier removal"
    )

    crs: Optional[str] = Field(
        None,
        description="Coordinate reference system (e.g., 'EPSG:3006'). Defaults to source CRS.",
    )


class PointCloudDataset(DatasetDescriptor):
    name = "pointcloud"
    ArgsModel = PointCloudArgs

    def validate(self, kwargs) -> PointCloudArgs:
        args = self.ArgsModel(**kwargs)
        if len(args.bounds) not in (4, 6):
            raise ValueError("Bounds must be a list of 4 or 6 floats.")

        return args

    def build(self, args: PointCloudArgs):
        if len(args.bounds) == 4:
            bounds = Bounds(
                xmin=args.bounds[0],
                ymin=args.bounds[1],
                xmax=args.bounds[2],
                ymax=args.bounds[3],
            )
        else:
            bounds = Bounds(
                xmin=args.bounds[0],
                ymin=args.bounds[1],
                zmin=args.bounds[2],
                xmax=args.bounds[3],
                ymax=args.bounds[4],
                zmax=args.bounds[5],
            )
        pc: PointCloud = dtcc_core.io.data.download_pointcloud(bounds=bounds)
        if args.classifications is not None and args.classifications not in (
            "all",
            "vegetation",
        ):
            classifications: List[int] = []
            if isinstance(args.classifications, int):
                classifications = [args.classifications]
            elif isinstance(args.classifications, list):
                classifications = args.classifications
            if isinstance(args.classifications, str):
                if args.classifications == "terrain":
                    classifications = [2, 8]
                elif args.classifications == "buildings":
                    classifications = [6, 9]
                elif args.classifications == "vegetation":
                    classifications = [3, 4, 5, 7]
            if len(classifications) > 0:
                pc = pc.classification_filter(classifications)
        if args.classifications == "vegetation":
            pc = pc.get_vegetation()
        if args.remove_outliers:
            pc = pc.remove_global_outliers(args.remove_outlier_threshold)

        if args.format is not None:
            with tempfile.NamedTemporaryFile(
                suffix="." + args.format, delete=True
            ) as buffer_file:
                pc.save(buffer_file.name, format=args.format)
                buffer_file.seek(0)
                pc_data = buffer_file.read()
            return pc_data
        return pc

    def show_options(self):
        return self.ArgsModel.model_json_schema()
