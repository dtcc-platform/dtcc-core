import dtcc_core
from dtcc_core.model import City, Bounds
from typing import Literal, Optional, List, Tuple, Sequence, Union
from pydantic import BaseModel, Field
import tempfile

from .dataset import DatasetDescriptor, DatasetBaseArgs
from dtcc_core.common.progress import ProgressTracker, report_progress


class CityArgs(DatasetBaseArgs):
    source: Literal["OSM", "LM"] = Field(
        "LM", description="Data source for building footprints"
    )
    smallest_building_size: float = Field(
        15.0, description="Smallest building size to include (in square meters)"
    )
    format: Optional[Literal["cityjson", "json"]] = Field(
        None, description="Output file format"
    )


class CityDataset(DatasetDescriptor):
    name = "city"
    description = "City model from point cloud data."
    ArgsModel = CityArgs

    def build(self, args: CityArgs):
        progress_phases = {
            "download_pointcloud": 0.30,
            "download_footprints": 0.30,
            "build_lod1": 0.15,
            "build_terrain": 0.10,
            "export": 0.15,
        }
        with ProgressTracker(total=1.0, phases=progress_phases) as progress:
            bounds = self.parse_bounds(args.bounds)
            city = City()
            city.bounds = bounds

            with progress.phase(
                "download_pointcloud", "Downloading point cloud data..."
            ):
                city.download_pointcloud()

            with progress.phase(
                "download_footprints", "Downloading building footprints..."
            ):
                city.download_footprints()

            with progress.phase("build_terrain", "Building terrain mesh..."):
                city.build_terrain(build_mesh=True)

            with progress.phase("build_lod1", "Building LOD1 city model..."):
                city.build_lod1_buildings()

            with progress.phase("export", "Exporting city model..."):
                if args.format is None:
                    return city
                else:
                    return self.export_to_bytes(city, "json")
