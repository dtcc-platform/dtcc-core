import dtcc_core
import numpy as np
from dtcc_core.model import City
from typing import Literal, Optional
from pydantic import Field

from .dataset import DatasetDescriptor, DatasetBaseArgs
from dtcc_core.common.progress import ProgressTracker


def _ground_level_from_raster(raster) -> float:
    valid_mask = np.isfinite(raster.data)
    if not np.isnan(raster.nodata):
        valid_mask &= raster.data != raster.nodata
    return float(raster.data[valid_mask].min())


class CitySurfaceMeshArgs(DatasetBaseArgs):
    max_mesh_size: float = Field(10.0, description="Maximum triangle size in meters")
    min_mesh_angle: float = Field(25.0, description="Minimum triangle angle in degrees")
    raster_cell_size: float = Field(
        2.0, description="Cell size for terrain raster in meters"
    )
    raster_radius: float = Field(
        3.0, description="Radius for terrain raster interpolation"
    )
    remove_outliers: bool = Field(
        True, description="Whether to remove global outliers from point cloud"
    )
    outlier_threshold: float = Field(
        3.0, description="Threshold for outlier removal (standard deviations)"
    )
    min_building_detail: float = Field(
        0.5, description="Minimum building feature size to resolve in meters"
    )
    min_building_area: float = Field(
        15.0, description="Smallest building footprint area to include (m²)"
    )
    merge_buildings: bool = Field(
        True, description="Whether to merge adjacent building footprints"
    )
    smoothing: int = Field(0, description="Number of terrain smoothing iterations")
    show_footprints: bool = Field(
        False,
        description="Whether to show a live matplotlib view of raw vs conditioned footprints before meshing",
    )
    footprint_cleaning_plot_block: bool = Field(
        True,
        description="Whether the optional footprint cleaning plot should block until the window is closed",
    )
    flat_ground: bool = Field(
        False,
        description="Whether to replace topography with a flat terrain raster",
    )
    ground_level: Optional[float] = Field(
        None,
        description="Ground level for flat terrain (defaults to minimum terrain elevation)",
    )
    format: Optional[Literal["obj", "stl", "vtu"]] = Field(
        None, description="Output file format"
    )


class CitySurfaceMeshDataset(DatasetDescriptor):
    name = "city_surface_mesh"
    description = (
        "Triangular surface mesh of a city with terrain and extruded buildings."
    )
    ArgsModel = CitySurfaceMeshArgs

    def build(self, args: CitySurfaceMeshArgs):
        """Build a city surface mesh from point cloud and building data.

        This method:
        1. Downloads point cloud and building footprints for the given bounds
        2. Removes outliers from the point cloud
        3. Builds a terrain raster
        4. Extracts roof points and computes building heights
        5. Creates a city model with terrain and buildings
        6. Generates a triangular surface mesh with extruded buildings

        Args:
            args: Validated arguments containing bounds and meshing parameters

        Returns:
            Mesh object or bytes if format is specified
        """
        progress_phases = {
            "download_pointcloud": 0.15,
            "download_footprints": 0.10,
            "remove_outliers": 0.05,
            "build_terrain": 0.10,
            "extract_roof_points": 0.05,
            "compute_building_heights": 0.05,
            "build_city": 0.03,
            "build_mesh": 0.42,
            "export": 0.05,
        }
        with ProgressTracker(total=1.0, phases=progress_phases) as progress:
            bounds = self.parse_bounds(args.bounds)

            with progress.phase(
                "download_pointcloud", "Downloading point cloud data..."
            ):
                pointcloud = dtcc_core.io.data.download_pointcloud(bounds=bounds)

            with progress.phase(
                "download_footprints", "Downloading building footprints..."
            ):
                buildings = dtcc_core.io.data.download_footprints(bounds=bounds)

            with progress.phase(
                "remove_outliers",
                (
                    "Removing outliers..."
                    if args.remove_outliers
                    else "Skipping outlier removal..."
                ),
            ):
                if args.remove_outliers:
                    pointcloud = pointcloud.remove_global_outliers(
                        args.outlier_threshold
                    )

            with progress.phase("build_terrain", "Building terrain raster..."):
                raster = dtcc_core.builder.build_terrain_raster(
                    pointcloud,
                    cell_size=args.raster_cell_size,
                    radius=args.raster_radius,
                    ground_only=True,
                )

            with progress.phase("extract_roof_points", "Extracting roof points..."):
                buildings = dtcc_core.builder.extract_roof_points(buildings, pointcloud)

            with progress.phase(
                "compute_building_heights", "Computing building heights..."
            ):
                buildings = dtcc_core.builder.compute_building_heights(
                    buildings, raster, overwrite=True
                )

            with progress.phase(
                "build_city",
                (
                    "Preparing flat-ground city model..."
                    if args.flat_ground
                    else "Assembling city model..."
                ),
            ):
                if args.flat_ground:
                    raster = dtcc_core.builder.flatten_terrain_raster(
                        raster, height=args.ground_level
                    )
                    buildings = dtcc_core.builder.set_building_heights_from_attribute(
                        buildings,
                        raster,
                        height_attribute="height",
                        default_ground_height=_ground_level_from_raster(raster),
                        always_use_default_ground=True,
                    )
                city = City()
                city.add_terrain(raster)
                city.add_buildings(buildings, remove_outside_terrain=True)

            with progress.phase("build_mesh", "Building city surface mesh..."):
                surface_mesh = dtcc_core.builder.build_city_surface_mesh(
                    city,
                    max_mesh_size=args.max_mesh_size,
                    min_mesh_angle=args.min_mesh_angle,
                    min_building_detail=args.min_building_detail,
                    min_building_area=args.min_building_area,
                    merge_buildings=args.merge_buildings,
                    smoothing=args.smoothing,
                    show_footprints=args.show_footprints,
                    footprint_cleaning_plot_block=args.footprint_cleaning_plot_block,
                )

            with progress.phase(
                "export",
                (
                    f"Exporting to {args.format}..."
                    if args.format
                    else "Preparing surface mesh..."
                ),
            ):
                if args.format is not None:
                    return self.export_to_bytes(surface_mesh, args.format)
                return surface_mesh
