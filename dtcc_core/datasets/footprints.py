import dtcc_core
from dtcc_core.model import Building, City
from typing import List, Literal, Optional
from pydantic import Field

from .dataset import DatasetDescriptor, DatasetBaseArgs
from dtcc_core.common.progress import ProgressTracker


class FootprintsArgs(DatasetBaseArgs):
    source: Literal["OSM", "LM"] = Field(
        "LM", description="Data source for building footprints"
    )
    smallest_building_size: float = Field(
        0.0, description="Smallest building size to include (in square meters)"
    )
    calculate_heights: bool = Field(
        False, description="Whether to calculate building heights from point cloud"
    )
    format: Optional[Literal["geojson", "gpkg", "shp.zip"]] = Field(
        None, description="Output file format"
    )


class FootprintsDataset(DatasetDescriptor):
    name = "building_footprints"
    description = "Building footprints"
    ArgsModel = FootprintsArgs
    data_category = "raw"
    result_kind = "building_footprints"
    python_return_type = "list[dtcc_core.model.Building]"

    def build(self, args: FootprintsArgs):
        progress_phases = {
            "download_footprints": 0.20,
            "download_pointcloud": 0.25,
            "compute_building_heights": 0.40,
            "export": 0.15,
        }

        with ProgressTracker(total=1.0, phases=progress_phases) as progress:
            bounds = self.parse_bounds(args.bounds)
            city = City()
            city.bounds = bounds
            with progress.phase(
                "download_footprints", "Downloading building footprints..."
            ):
                city.download_footprints()

            if args.calculate_heights:
                with progress.phase(
                    "download_pointcloud", "Downloading point cloud data..."
                ):
                    city.download_pointcloud()
                with progress.phase(
                    "compute_building_heights", "Computing building heights..."
                ):
                    city.building_heights_from_pointcloud(keep_roof_points=False)

            else:
                progress.update(increment=65)

            with progress.phase("export", "Exporting footprints..."):
                if args.format is None:
                    return city.buildings
                else:
                    return self.export_to_bytes(
                        city, args.format, save_callable=dtcc_core.io.footprints.save
                    )
