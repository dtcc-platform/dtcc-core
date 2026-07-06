import dtcc_core
import numpy as np
from pathlib import Path
from dtcc_core.model import City, GeometryType
from typing import Any, Literal, Optional
from pydantic import Field

from .dataset import DatasetDescriptor, DatasetBaseArgs
from ._city_mesh_common import prepare_city_from_bounds
from .providers import provider_entry
from dtcc_core.common.progress import ProgressTracker


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
    lod: Optional[GeometryType] = Field(
        None,
        description="Optional building geometry level of detail for meshing",
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
    min_building_area: float = Field(
        15.0, description="Smallest building footprint area to include (m²)"
    )
    merge_buildings: bool = Field(
        True, description="Whether to merge adjacent building footprints"
    )
    smoothing: int = Field(0, description="Number of terrain smoothing iterations")
    max_volume: Optional[float] = Field(
        None,
        description=(
            "Maximum tetrahedron volume used as TetGen's 3D size cap "
            "(defaults to the regular-tetrahedron volume implied by "
            "max_mesh_size if not set)"
        ),
    )
    mesher: Optional[Literal["auto", "dtcc_mesher", "triangle"]] = Field(
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
    tetgen_switches: Optional[dict[str, Any]] = Field(
        None,
        description="Optional descriptive TetGen switch dictionary forwarded to the volume builder",
    )
    tetgen_debug_output_dir: Optional[str | Path] = Field(
        None,
        description="Optional directory for saved TetGen input artifacts",
    )
    tetgen_debug_output_stem: Optional[str] = Field(
        None,
        description="Optional filename stem for saved TetGen input artifacts",
    )
    tetgen_quality_failure_output_dir: Optional[str | Path] = Field(
        None,
        description="Optional directory for TetGen quality-failure reports",
    )
    tetgen_quality_failure_output_stem: Optional[str] = Field(
        None,
        description="Optional filename stem for TetGen quality-failure reports",
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
    format: Optional[Literal["xdmf", "vtu"]] = Field(
        None, description="Output file format"
    )


class CityVolumeMeshDataset(DatasetDescriptor):
    name = "city_volume_mesh"
    title = "City Volume Mesh"
    description = (
        "Tetrahedral volume mesh for FEM/CFD preprocessing, derived from point "
        "cloud terrain and default DTCC footprint sources."
    )
    ArgsModel = CityVolumeMeshArgs
    data_category = "derived"
    result_kind = "mesh"
    python_return_type = "dtcc_core.model.VolumeMesh"
    multi_file_formats = ("xdmf",)
    provider = [
        provider_entry("lantmateriet", role="source_provider"),
        provider_entry("dtcc-platform", role="processor"),
        {"name": "TetGen", "role": "tetrahedral_meshing_backend"},
    ]
    source = [
        {
            "name": "DTCC/Lantmäteriet point cloud backend/cache",
            "role": "upstream_dataset",
            "source_terms_status": "requires_review",
        },
        {
            "name": "DTCC/Lantmäteriet footprint backend/cache",
            "role": "upstream_dataset",
            "source_terms_status": "requires_review",
        },
    ]
    license = (
        "Requires review: derived from upstream point cloud and footprint data; "
        "verify source terms before redistribution."
    )
    collection_period = (
        "Requires review: inherits point cloud acquisition date and footprint "
        "provider/cache vintage, which are not currently surfaced in the result."
    )
    default_crs = "EPSG:3006"
    lod = "3D computational volume with terrain, building, and domain boundary surfaces"
    data_types = [
        "volume_mesh",
        "tetrahedral_mesh",
        "terrain",
        "buildings",
        "fem",
        "cfd",
        "boundary_markers",
    ]
    geographic_coverage = "Sweden, constrained by requested bounds and source coverage"
    update_frequency = "derived on demand from upstream source data"
    processing_steps = [
        "Download point cloud and DTCC footprint data for the requested bounds in EPSG:3006",
        "Optionally remove global point-cloud outliers using outlier_threshold",
        "Build a terrain raster with raster_cell_size and raster_radius",
        "Optionally replace topography with a flat terrain raster when flat_ground=True",
        "Extract roof points and compute building heights from the terrain raster",
        "Condition building footprints using lod, min_building_detail, min_building_area, merge_buildings, and pipeline_mode",
        "Build intermediate ground, shell, and top-cap surfaces using max_mesh_size, top_cap_max_mesh_size, min_mesh_angle, smoothing, and mesher",
        "Run TetGen with max_volume and tetgen_switches/tetgen_extra to produce tetrahedra",
        "Optionally attach boundary face markers (-1 ground, -2 top, -3 west, -4 east, -5 south, -6 north)",
        "Return the native VolumeMesh object or serialize VTU/XDMF output",
    ]
    derived_from = [
        {
            "name": "point_cloud",
            "relationship": "terrain and building-height estimation",
            "source_terms_status": "requires_review",
        },
        {
            "name": "building_footprints",
            "relationship": "building outlines and domain obstacles",
            "source_terms_status": "requires_review",
        },
    ]
    presentation_headline = "Tetrahedral City Volume Mesh"
    presentation_summary = (
        "A tetrahedral computational-domain mesh for downstream FEM or CFD "
        "experiments over the requested city bounds."
    )
    presentation_narrative = [
        {
            "heading": "What you are seeing",
            "body": (
                "The result is a 3D tetrahedral mesh bounded by terrain, the top "
                "cap, side walls, and building obstacle surfaces."
            ),
        },
        {
            "heading": "How it is made",
            "body": (
                "The dataset prepares a city model from point cloud and footprint "
                "sources, conditions building footprints, constructs intermediate "
                "surface meshes, and tetrahedralizes the domain with TetGen."
            ),
        },
        {
            "heading": "Limitations",
            "body": (
                "The mesh is a preprocessing artifact. Boundary markers are labels, "
                "not physical boundary conditions, and the dataset does not prove "
                "solver convergence or scientific validity."
            ),
        },
    ]
    key_points = [
        "domain_height sets the top of the computational domain in meters",
        "max_volume defaults to the regular-tetrahedron volume implied by max_mesh_size",
        "boundary_face_markers=True labels ground, top, and side boundary faces",
        "XDMF is a multi-file export format; VTU is single-file",
        "TetGen/dtcc-tetgen-wrapper availability is required for actual volume meshing",
    ]
    presentation_legend = {
        "title": "Volume mesh boundaries",
        "entries": [
            {"label": "-1 ground", "meaning": "terrain boundary face marker"},
            {"label": "-2 top", "meaning": "top cap boundary face marker"},
            {"label": "-3 west / -4 east", "meaning": "xmin/xmax side boundary markers"},
            {"label": "-5 south / -6 north", "meaning": "ymin/ymax side boundary markers"},
            {"label": "Tetrahedron", "meaning": "interior computational volume element"},
        ],
    }
    view_hints = {
        "preferred_geometry": "volume_mesh",
        "default_crs": "EPSG:3006",
        "mesh_role": "fem_cfd_preprocessing",
        "boundary_marker_convention": {
            "ground": -1,
            "top": -2,
            "west_xmin": -3,
            "east_xmax": -4,
            "south_ymin": -5,
            "north_ymax": -6,
        },
    }
    presentation_warnings = [
        "Source and license terms inherit upstream point-cloud and footprint review status.",
        "TetGen/dtcc-tetgen-wrapper must be available for actual volume meshing.",
        "Boundary markers do not define simulation boundary conditions by themselves.",
        "The generated mesh is not automatically validated for solver stability or scientific accuracy.",
    ]
    presentation_limitations = [
        "Point-cloud acquisition date, density, and footprint cache vintage are not currently surfaced in the result.",
        "Mesh quality depends on source data, conditioning choices, TetGen settings, and domain_height.",
        "No CFD/FEM solver assumptions, material parameters, or boundary conditions are included in this dataset.",
    ]

    @staticmethod
    def _tetgen_switch_payload(args: CityVolumeMeshArgs) -> dict[str, Any]:
        switches = dict(args.tetgen_switches or {})
        extra = str(switches.get("extra", ""))
        if args.tetgen_extra:
            extra = f"{extra}{args.tetgen_extra}"
        switches["extra"] = extra
        return switches

    def _build_mesh_from_city(self, city: City, args: CityVolumeMeshArgs):
        max_vol = (
            args.max_volume
            if args.max_volume is not None
            else _regular_tet_volume(args.max_mesh_size)
        )
        stage_audit = {} if args.stage_audit_enabled else None
        return dtcc_core.builder.build_city_volume_mesh(
            city,
            lod=args.lod,
            max_mesh_size=args.max_mesh_size,
            top_cap_max_mesh_size=args.top_cap_max_mesh_size,
            domain_height=args.domain_height,
            min_mesh_angle=args.min_mesh_angle,
            merge_buildings=args.merge_buildings,
            min_building_detail=args.min_building_detail,
            min_building_area=args.min_building_area,
            smoothing=args.smoothing,
            boundary_face_markers=args.boundary_face_markers,
            tetgen_switches=self._tetgen_switch_payload(args),
            report_mesh_quality=args.report_mesh_quality,
            show_footprints=args.show_footprints,
            footprint_cleaning_plot_block=args.footprint_cleaning_plot_block,
            mesher=args.mesher,
            tetgen_debug_output_dir=args.tetgen_debug_output_dir,
            tetgen_debug_output_stem=args.tetgen_debug_output_stem,
            tetgen_quality_failure_output_dir=args.tetgen_quality_failure_output_dir,
            tetgen_quality_failure_output_stem=args.tetgen_quality_failure_output_stem,
            stage_audit=stage_audit,
            pipeline_mode=args.pipeline_mode,
            max_volume=max_vol,
        )

    def build_from_city(self, city: City, **kwargs):
        args = self.validate(kwargs)
        volume_mesh = self._build_mesh_from_city(city, args)
        if args.format is not None:
            return self.export_to_bytes(volume_mesh, args.format)
        return volume_mesh

    def build(self, args: CityVolumeMeshArgs):
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

            with progress.phase("build_mesh", "Building volume mesh..."):
                volume_mesh = self._build_mesh_from_city(city, args)

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
