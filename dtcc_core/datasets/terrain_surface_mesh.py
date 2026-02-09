import dtcc_core
from dtcc_core.model import City, PointCloud, Bounds, Terrain, Raster
from typing import Literal, Optional, List, Tuple, Sequence, Union
from pydantic import BaseModel, Field
import tempfile

from .dataset import DatasetDescriptor, DatasetBaseArgs
from dtcc_core.common.progress import ProgressTracker, report_progress


class TerrainSurfaceMeshArgs(DatasetBaseArgs):
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


class TerrainSurfaceMeshDataset(DatasetDescriptor):
    name = "terrain_surface_mesh"
    description = "Generate a terrain surface mesh from point cloud data."
    ArgsModel = TerrainSurfaceMeshArgs

    def build(self, args: TerrainSurfaceMeshArgs):
        progress_phases = {
            "download_pointcloud": 0.40,
            "remove_outliers": 0.10,
            "build_terrain": 0.40,
            "export": 0.10,
        }
        with ProgressTracker(total=1.0, phases=progress_phases) as progress:
            bounds = self.parse_bounds(args.bounds)

            with progress.phase(
                "download_pointcloud", "Downloading point cloud data..."
            ):
                pc = dtcc_core.io.data.download_pointcloud(bounds=bounds)

            with progress.phase(
                "remove_outliers",
                (
                    f"Removing outliers (threshold={args.remove_outlier_threshold})..."
                    if args.remove_outliers
                    else "Skipping outlier removal"
                ),
            ):
                if args.remove_outliers:
                    pc = pc.remove_global_outliers(args.remove_outlier_threshold)

            with progress.phase(
                "build_terrain",
                (
                    "Building terrain raster..."
                    if args.format == "tif"
                    else "Building terrain surface mesh..."
                ),
            ):
                if args.format == "tif":
                    result = dtcc_core.builder.build_terrain_raster(
                        pc, cell_size=args.raster_resolution
                    )
                else:
                    result = dtcc_core.builder.build_terrain_surface_mesh(
                        pc,
                        max_mesh_size=args.mesh_resolution,
                        smoothing=args.smoothing,
                    )

            with progress.phase(
                "export",
                (
                    f"Exporting terrain to {args.format}..."
                    if args.format
                    else "Preparing terrain result..."
                ),
            ):
                if args.format == "tif":
                    return self.export_to_bytes(result, "tif")
                elif args.format is None:
                    return result
                elif args.format in ("obj", "stl"):
                    return self.export_to_bytes(result, args.format)
