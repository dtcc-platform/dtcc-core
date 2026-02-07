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

    @staticmethod
    def _resolve_classifications(classifications) -> List[int]:
        """Resolve classification parameter to a list of class IDs."""
        if isinstance(classifications, int):
            return [classifications]
        elif isinstance(classifications, list):
            return classifications
        elif isinstance(classifications, str):
            mapping = {
                "terrain": [2, 8],
                "buildings": [6, 9],
                "vegetation": [3, 4, 5, 7],
            }
            return mapping.get(classifications, [])
        return []

    def build(self, args: PointCloudArgs):
        progress_phases = {
            "download_pointcloud": 0.50,
            "filter_classifications": 0.15,
            "remove_outliers": 0.15,
            "export": 0.20,
        }
        with ProgressTracker(total=1.0, phases=progress_phases) as progress:
            bounds = self.parse_bounds(args.bounds)

            with progress.phase("download_pointcloud", "Downloading point cloud data..."):
                pc: PointCloud = dtcc_core.io.data.download_pointcloud(bounds=bounds)

            with progress.phase(
                "filter_classifications",
                "Filtering point cloud classifications..."
                if args.classifications not in ("all", None)
                else "Using all classifications",
            ):
                if args.classifications == "vegetation":
                    report_progress(percent=30, message="Extracting vegetation points...")
                    pc = pc.get_vegetation()
                elif args.classifications not in ("all", None):
                    classifications = self._resolve_classifications(args.classifications)
                    if classifications:
                        report_progress(
                            percent=30,
                            message=f"Filtering to classes {classifications}...",
                        )
                        pc = pc.classification_filter(classifications)

            with progress.phase(
                "remove_outliers",
                f"Removing outliers (threshold={args.remove_outlier_threshold})..."
                if args.remove_outliers
                else "Skipping outlier removal",
            ):
                if args.remove_outliers:
                    pc = pc.remove_global_outliers(args.remove_outlier_threshold)

            with progress.phase(
                "export",
                f"Exporting to {args.format}..." if args.format
                else "Preparing point cloud result...",
            ):
                if args.format is not None:
                    return self.export_to_bytes(pc, args.format)
                return pc