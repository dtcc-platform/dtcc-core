import dtcc_core
from dtcc_core.model import City, Bounds
from typing import Literal, Optional, List, Tuple, Sequence, Union
from pydantic import BaseModel, Field
import tempfile

from .dataset import DatasetDescriptor, DatasetBaseArgs
from dtcc_core.common.progress import ProgressTracker, report_progress


class BuildingArgs(DatasetBaseArgs):
    source: Literal["OSM", "LM"] = Field(
        "LM", description="Data source for building footprints"
    )
    smallest_building_size: float = Field(
        15.0, description="Smallest building size to include (in square meters)"
    )
    place_on_zero: bool = Field(
        False, description="Whether to place buildings on Z=0 plane"
    )
    format: Optional[Literal["obj", "stl", "cityjson", "json"]] = Field(
        None, description="Output file format"
    )


class BuildingDataset(DatasetDescriptor):
    name = "buildings"
    description = "3D buildings (LoD1) from point cloud and building footprints."
    ArgsModel = BuildingArgs

    def build(self, args: BuildingArgs):
        progress_phases = {
            "download_pointcloud": 0.25,
            "download_footprints": 0.20,
            "remove_outliers": 0.05,
            "compute_building_heights": 0.15,
            "build_lod1": 0.25,
            "export": 0.10,
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

            with progress.phase("remove_outliers", "Removing point cloud outliers..."):
                city.add_point_cloud(city.pointcloud.remove_global_outliers(3.0))

            with progress.phase(
                "compute_building_heights", "Computing building heights..."
            ):
                report_progress(percent=10, message="Building terrain raster...")
                terrain_raster = dtcc_core.builder.build_terrain_raster(
                    city.pointcloud,
                    cell_size=2.0,
                    ground_only=True,
                    _report_progress=False,
                )
                report_progress(percent=40, message="Extracting roof points...")
                buildings = dtcc_core.builder.extract_roof_points(
                    city.buildings, city.pointcloud
                )
                report_progress(
                    percent=70, message="Computing heights from roof points..."
                )
                buildings = dtcc_core.builder.compute_building_heights(
                    buildings, terrain_raster, overwrite=True
                )

            with progress.phase("build_lod1", "Extruding LOD1 geometry..."):
                buildings = dtcc_core.builder.build_lod1_buildings(buildings)

            with progress.phase(
                "export",
                (
                    f"Exporting buildings to {args.format}..."
                    if args.format
                    else "Preparing building result..."
                ),
            ):
                if args.format is None:
                    return buildings
                elif args.format in ("cityjson", "json"):
                    return self.export_to_bytes(city, "json", as_text=True)
                elif args.format in ("obj", "stl"):
                    report_progress(percent=20, message="Extracting building meshes...")
                    building_meshes = [
                        b.lod1.mesh(weld=True, snap=0.005) for b in buildings
                    ]
                    if args.place_on_zero:
                        for b in building_meshes:
                            b.offset([0, 0, -b.bounds.zmin])
                    report_progress(percent=60, message="Merging meshes...")
                    merged_mesh = dtcc_core.builder.meshing.merge_meshes(
                        building_meshes
                    )
                    report_progress(percent=80, message=f"Writing {args.format}...")
                    return self.export_to_bytes(merged_mesh, args.format)
