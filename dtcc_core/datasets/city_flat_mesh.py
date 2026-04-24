import dtcc_core
from dtcc_core.model import City, GeometryType
from typing import Literal, Optional
from pydantic import Field

from .dataset import DatasetDescriptor, DatasetBaseArgs
from ._city_mesh_common import prepare_city_from_bounds
from dtcc_core.common.progress import ProgressTracker


class CityFlatMeshArgs(DatasetBaseArgs):
    max_mesh_size: Optional[float] = Field(
        10.0, description="Maximum triangle size in meters"
    )
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
    show_footprints: bool = Field(
        False,
        description="Whether to show a live matplotlib view of raw vs conditioned footprints before meshing",
    )
    footprint_cleaning_plot_block: bool = Field(
        True,
        description="Whether the optional footprint cleaning plot should block until the window is closed",
    )
    mesher: Optional[Literal["auto", "dtcc_mesher", "triangle", "spade"]] = Field(
        None,
        description="2D meshing backend to use for the flat-mesh triangulation",
    )
    report_mesh_quality: bool = Field(
        True,
        description="Whether to log a mesh-quality summary after meshing",
    )
    stage_audit_enabled: bool = Field(
        False,
        description="Whether to attach per-stage audit data to the returned mesh",
    )
    pipeline_mode: Literal["strict"] = Field(
        "strict",
        description="Meshing pipeline mode",
    )
    format: Optional[Literal["obj", "stl", "vtu"]] = Field(
        None, description="Output file format"
    )


class CityFlatMeshDataset(DatasetDescriptor):
    name = "city_flat_mesh"
    description = (
        "Flat 2D triangular mesh at z=0 with building footprints marked as subdomains."
    )
    ArgsModel = CityFlatMeshArgs

    def _build_mesh_from_city(self, city: City, args: CityFlatMeshArgs):
        stage_audit = {} if args.stage_audit_enabled else None
        flat_mesh = dtcc_core.builder.build_city_flat_mesh(
            city,
            lod=GeometryType.LOD0,
            max_mesh_size=args.max_mesh_size,
            min_mesh_angle=args.min_mesh_angle,
            min_building_detail=args.min_building_detail,
            min_building_area=args.min_building_area,
            merge_buildings=args.merge_buildings,
            show_footprints=args.show_footprints,
            footprint_cleaning_plot_block=args.footprint_cleaning_plot_block,
            mesher=args.mesher,
            report_mesh_quality=args.report_mesh_quality,
            pipeline_mode=args.pipeline_mode,
            stage_audit=stage_audit,
        )
        return flat_mesh

    def build_from_city(self, city: City, **kwargs):
        args = self.validate(kwargs)
        flat_mesh = self._build_mesh_from_city(city, args)
        if args.format is not None:
            return self.export_to_bytes(flat_mesh, args.format)
        return flat_mesh

    def build(self, args: CityFlatMeshArgs):
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
                progress=progress,
            )

            with progress.phase("build_mesh", "Building city flat mesh..."):
                flat_mesh = self._build_mesh_from_city(city, args)

            with progress.phase(
                "export",
                (
                    f"Exporting to {args.format}..."
                    if args.format
                    else "Preparing flat mesh..."
                ),
            ):
                if args.format is not None:
                    return self.export_to_bytes(flat_mesh, args.format)
                return flat_mesh
