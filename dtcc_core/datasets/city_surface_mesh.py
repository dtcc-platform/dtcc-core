import dtcc_core
from dtcc_core.model import City
from typing import Literal, Optional
from pydantic import Field

from .dataset import DatasetDescriptor, DatasetBaseArgs
from ._city_mesh_common import prepare_city_from_bounds
from .providers import provider_entry
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
    merge_tolerance: float = Field(
        0.5, description="Distance tolerance for merging footprints in meters"
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


class CitySurfaceMeshDataset(DatasetDescriptor):
    name = "city_surface_mesh"
    title = "City Surface Mesh"
    description = (
        "Triangular surface mesh of terrain and extruded building surfaces, "
        "prepared from point cloud data and default DTCC footprint sources."
    )
    ArgsModel = CitySurfaceMeshArgs
    data_category = "derived"
    result_kind = "mesh"
    python_return_type = "dtcc_core.model.Mesh"
    provider = [
        provider_entry("lantmateriet", role="source_provider"),
        provider_entry("dtcc-platform", role="processor"),
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
    lod = "Terrain surface plus extruded building surface mesh"
    data_types = ["terrain", "buildings", "surface_mesh", "mesh", "lod1"]
    geographic_coverage = "Sweden, constrained by requested bounds and source coverage"
    update_frequency = "derived on demand from upstream source data"
    processing_steps = [
        "Download point cloud and DTCC footprint data for the requested bounds in EPSG:3006",
        "Optionally remove global point-cloud outliers using outlier_threshold",
        "Build a terrain raster with raster_cell_size and raster_radius",
        "Optionally replace topography with a flat terrain raster when flat_ground=True",
        "Extract roof points and compute building heights from the terrain raster",
        "Condition footprints using min_building_detail, min_building_area, merge_buildings, merge_tolerance, and pipeline_mode",
        "Generate a terrain plus extruded-building surface mesh using max_mesh_size, min_mesh_angle, smoothing, mesher, and mesh-quality reporting settings",
        "Return the native Mesh object or serialize the requested mesh format",
    ]
    derived_from = [
        {
            "name": "point_cloud",
            "relationship": "terrain and building-height estimation",
            "source_terms_status": "requires_review",
        },
        {
            "name": "building_footprints",
            "relationship": "building outlines for surface meshing",
            "source_terms_status": "requires_review",
        },
    ]
    presentation_headline = "Terrain and Building Surface Mesh"
    presentation_summary = (
        "A triangular shell mesh combining point-cloud terrain with generalized "
        "building surfaces."
    )
    presentation_narrative = [
        {
            "heading": "What you are seeing",
            "body": (
                "The mesh is a surface representation of the city envelope: "
                "terrain triangles plus extruded building walls and roofs."
            ),
        },
        {
            "heading": "How it is made",
            "body": (
                "The dataset prepares a City object from point cloud and footprint "
                "sources, conditions building footprints, and runs the DTCC surface "
                "meshing pipeline with the requested size, angle, smoothing, and backend settings."
            ),
        },
        {
            "heading": "Limitations",
            "body": (
                "The result is preprocessing geometry. It is not automatically "
                "validated as a watertight CFD/FEM boundary or as surveyed LoD geometry."
            ),
        },
    ]
    key_points = [
        "max_mesh_size and min_mesh_angle control terrain/building surface triangles",
        "min_building_area and min_building_detail condition small or complex footprints",
        "flat_ground=True replaces terrain topography before meshing",
        "report_mesh_quality controls logging; stage_audit_enabled requests diagnostic data from the builder",
    ]
    presentation_legend = {
        "title": "Surface mesh layers",
        "entries": [
            {"label": "Terrain triangle", "meaning": "point-cloud-derived ground surface facet"},
            {"label": "Building surface", "meaning": "extruded and conditioned footprint boundary"},
            {"label": "Conditioned footprint", "meaning": "merged or simplified source footprint used for meshing"},
        ],
    }
    view_hints = {
        "preferred_geometry": "surface_mesh",
        "default_crs": "EPSG:3006",
        "mesh_role": "visualization_or_surface_preprocessing",
        "building_lod": "generalized extruded surfaces",
    }
    presentation_warnings = [
        "Source and license terms inherit upstream point-cloud and footprint review status.",
        "Surface mesh quality depends on source data, footprint conditioning, and selected meshing parameters.",
        "The dataset does not prove watertightness, solver suitability, or scientific validity.",
    ]
    presentation_limitations = [
        "Point-cloud acquisition date, density, and footprint cache vintage are not currently surfaced in the result.",
        "Buildings are generalized from footprints and height estimates, not detailed roof/facade reconstructions.",
        "Mesh-quality reports are builder diagnostics and do not replace downstream solver validation.",
    ]

    def _build_mesh_from_city(self, city: City, args: CitySurfaceMeshArgs):
        stage_audit = {} if args.stage_audit_enabled else None
        return dtcc_core.builder.build_city_surface_mesh(
            city,
            max_mesh_size=args.max_mesh_size,
            min_mesh_angle=args.min_mesh_angle,
            min_building_detail=args.min_building_detail,
            min_building_area=args.min_building_area,
            merge_buildings=args.merge_buildings,
            merge_tolerance=args.merge_tolerance,
            smoothing=args.smoothing,
            show_footprints=args.show_footprints,
            footprint_cleaning_plot_block=args.footprint_cleaning_plot_block,
            mesher=args.mesher,
            report_mesh_quality=args.report_mesh_quality,
            pipeline_mode=args.pipeline_mode,
            stage_audit=stage_audit,
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
