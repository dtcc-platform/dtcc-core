import dtcc_core
from dtcc_core.model import PointCloud, Bounds
from typing import Literal, Optional, List, Tuple, Sequence, Union
from pydantic import BaseModel, Field
import tempfile

from .dataset import DatasetDescriptor, DatasetBaseArgs
from dtcc_core.common.progress import ProgressTracker, report_progress

class PointCloudArgs(DatasetBaseArgs):
    classifications: Union[
        int, list[int], Literal["all", "terrain", "buildings", "vegetation"]
    ] = Field(
        "all",
        description="Which classifications to include (e.g., [2, 9] for ground and water, or 'vegetation' for all vegetation classes). 'all' includes all points.",
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
    description = "Download point cloud data with optional classification filtering and outlier removal."
    ArgsModel = PointCloudArgs

    def build(self, args: PointCloudArgs):
        
        progress_phases = {
            "download_pointcloud": 0.30,
            "preprocess": 0.20,
            "extract_vegetation": 0.10,
            "remove_outliers": 0.20,
            "extract": 0.20,  # The big one
        }
        with ProgressTracker(total=1.0, phases=progress_phases) as progress:
            bounds = self.parse_bounds(args.bounds)
            
            with progress.phase("download_pointcloud", "Downloading point cloud..."):
                pc: PointCloud = dtcc_core.io.data.download_pointcloud(bounds=bounds)

            with progress.phase("preprocess", "Preprocessing point cloud..."):
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

            with progress.phase("extract_vegetation", "Extracting vegetation points..."):
                if args.classifications == "vegetation":
                    pc = pc.get_vegetation()

            with progress.phase("remove_outliers", "Removing outliers..."):
                if args.remove_outliers:
                    pc = pc.remove_global_outliers(args.remove_outlier_threshold)

            with progress.phase("extract", "Extracting point cloud data..."):
                if args.format is not None:
                    return self.export_to_bytes(pc, args.format)
                return pc