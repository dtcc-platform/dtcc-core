import dtcc_core
from dtcc_core.model import City
from typing import Literal, Optional
from pydantic import Field

from .dataset import DatasetDescriptor, DatasetBaseArgs
from ._city_mesh_common import prepare_city_from_bounds
from dtcc_core.common.progress import ProgressTracker


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
    mesher: Optional[Literal["auto", "dtcc_mesher", "triangle"]] = Field(
        None,
        description="2D meshing backend to use for the ground triangulation",
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

    def _build_mesh_from_city(self, city: City, args: CitySurfaceMeshArgs):
        return dtcc_core.builder.build_city_surface_mesh(
            city,
            max_mesh_size=args.max_mesh_size,
            min_mesh_angle=args.min_mesh_angle,
            min_building_detail=args.min_building_detail,
            min_building_area=args.min_building_area,
            merge_buildings=args.merge_buildings,
            smoothing=args.smoothing,
            show_footprints=args.show_footprints,
            footprint_cleaning_plot_block=args.footprint_cleaning_plot_block,
            mesher=args.mesher,
        )

    def build_from_city(self, city: City, **kwargs):
        args = self.validate(kwargs)
        surface_mesh = self._build_mesh_from_city(city, args)
        if args.format is not None:
            return self.export_to_bytes(surface_mesh, args.format)
        return surface_mesh

    def build(self, args: CitySurfaceMeshArgs):
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
            city = prepare_city_from_bounds(
                bounds,
                raster_cell_size=args.raster_cell_size,
                raster_radius=args.raster_radius,
                remove_outliers=args.remove_outliers,
                outlier_threshold=args.outlier_threshold,
                flat_ground=args.flat_ground,
                ground_level=args.ground_level,
                progress=progress,
            )

            with progress.phase("build_mesh", "Building city surface mesh..."):
                surface_mesh = self._build_mesh_from_city(city, args)

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
