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


def _regular_tet_volume(edge_length: float) -> float:
    return (np.sqrt(2.0) / 12.0) * float(edge_length) ** 3


class CityVolumeMeshArgs(DatasetBaseArgs):
    max_mesh_size: float = Field(
        25.0,
        description="Maximum target edge size for the 2D ground and shell meshing stages in meters",
    )
    top_cap_max_mesh_size: Optional[float] = Field(
        None,
        description="Optional separate target edge size for the lifted top cap triangulation in meters",
    )
    domain_height: float = Field(
        80.0, description="Height of the computational domain (H parameter) in meters"
    )
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
    boundary_face_markers: bool = Field(
        True,
        description=(
            "Whether to add boundary face markers to the mesh "
            "(-1 ground, -2 top, -3 west/xmin, -4 east/xmax, "
            "-5 south/ymin, -6 north/ymax)"
        ),
    )
    min_building_detail: float = Field(
        0.5, description="Minimum building feature size to resolve in meters"
    )
    max_volume: Optional[float] = Field(
        None,
        description=(
            "Maximum tetrahedron volume used as TetGen's 3D size cap "
            "(defaults to the regular-tetrahedron volume implied by "
            "max_mesh_size if not set)"
        ),
    )
    mesher: Optional[Literal["auto", "dtcc_mesher", "triangle", "spade"]] = Field(
        None,
        description=(
            "2D meshing backend for the intermediate flat and shell meshes "
            "(defaults to dtcc_mesher when omitted)"
        ),
    )
    tetgen_extra: str = Field(
        "",
        description="Extra switches to pass to TetGen (e.g., 'VV' for verbose output)",
    )
    flat_ground: bool = Field(
        False,
        description="Whether to replace topography with a flat terrain raster",
    )
    show_footprints: bool = Field(
        False,
        description="Whether to show a live matplotlib view of raw vs conditioned footprints before meshing",
    )
    footprint_cleaning_plot_block: bool = Field(
        True,
        description="Whether the optional footprint cleaning plot should block until the window is closed",
    )
    ground_level: Optional[float] = Field(
        None,
        description="Ground level for flat terrain (defaults to minimum terrain elevation)",
    )
    format: Optional[Literal["xdmf", "vtu"]] = Field(
        None, description="Output file format"
    )


class CityVolumeMeshDataset(DatasetDescriptor):
    name = "city_volume_mesh"
    description = "Tetrahedral volume mesh from point cloud and building data, suitable for CFD/FEM simulations."
    ArgsModel = CityVolumeMeshArgs

    def build(self, args: CityVolumeMeshArgs):
        """Build a volume mesh from point cloud and building data.

        This method:
        1. Downloads point cloud and building footprints for the given bounds
        2. Removes outliers from the point cloud
        3. Builds a terrain raster
        4. Extracts roof points and computes building heights
        5. Creates a city model with terrain and buildings
        6. Generates a tetrahedral volume mesh suitable for simulations

        Args:
            args: Validated arguments containing bounds and meshing parameters

        Returns:
            VolumeMesh object or bytes if format is specified
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

            with progress.phase("build_mesh", "Building volume mesh..."):
                max_vol = (
                    args.max_volume
                    if args.max_volume is not None
                    else _regular_tet_volume(args.max_mesh_size)
                )
                volume_mesh = dtcc_core.builder.build_city_volume_mesh(
                    city,
                    max_mesh_size=args.max_mesh_size,
                    top_cap_max_mesh_size=args.top_cap_max_mesh_size,
                    domain_height=args.domain_height,
                    min_building_detail=args.min_building_detail,
                    boundary_face_markers=args.boundary_face_markers,
                    show_footprints=args.show_footprints,
                    footprint_cleaning_plot_block=args.footprint_cleaning_plot_block,
                    mesher=args.mesher,
                    max_volume=max_vol,
                    tetgen_switches={
                        "extra": args.tetgen_extra,
                    },
                )

            with progress.phase(
                "export",
                (
                    f"Exporting to {args.format}..."
                    if args.format
                    else "Preparing volume mesh..."
                ),
            ):
                if args.format is not None:
                    return self.export_to_bytes(volume_mesh, args.format)
                return volume_mesh
