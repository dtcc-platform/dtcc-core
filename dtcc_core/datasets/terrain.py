import dtcc_core
from dtcc_core.model import City, PointCloud, Bounds, Terrain, Raster
from typing import Literal, Optional, List, Tuple, Sequence, Union
from pydantic import BaseModel, Field
import tempfile

from .dataset import DatasetDescriptor, DatasetBaseArgs


class TerrainArgs(DatasetBaseArgs):
    raster_resolution: float = Field(
        2, description="Resolution of the terrain raster in meters"
    )
    mesh_resolution: float = Field(
        5, description="Resolution of the terrain mesh in meters"
    )
    smoothing: int = Field(
        3, description="Number of smoothing iterations to apply to the terrain mesh"
    )

    remove_outliers: bool = Field(
        True, description="Whether to remove global outliers from the terrain raster"
    )
    remove_outlier_threshold: float = Field(
        3.0, description="Threshold for outlier removal"
    )
    format: Optional[Literal["tif", "obj", "stl"]] = Field(
        None, description="Output file format"
    )


class TerrainDataset(DatasetDescriptor):
    name = "terrain"
    ArgsModel = TerrainArgs

    def build(self, args: TerrainArgs):
        bounds = self.parse_bounds(args.bounds)
        pc = dtcc_core.io.data.download_pointcloud(bounds=bounds)
        if args.remove_outliers:
            pc = pc.remove_global_outliers(args.remove_outlier_threshold)
        if args.format == "tif":
            raster = dtcc_core.builder.build_terrain_raster(
                pc, cell_size=args.raster_resolution
            )
            self.export_to_bytes(raster, "tif")
        else:
            terrain_mesh = dtcc_core.builder.build_terrain_mesh(
                pc,
                max_mesh_size=args.mesh_resolution,
                smoothing=args.smoothing,
            )
            if args.format is None:
                return terrain_mesh
            elif args.format in ("obj", "stl"):
                return self.export_to_bytes(terrain_mesh, args.format)
