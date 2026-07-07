import dtcc_core
from typing import Literal, Optional
from pydantic import Field

from .dataset import DatasetDescriptor, DatasetBaseArgs
from .providers import provider_entry
from dtcc_core.common.progress import ProgressTracker


class TerrainSurfaceMeshArgs(DatasetBaseArgs):
    raster_resolution: float = Field(
        2, description="Resolution of the terrain raster in meters"
    )
    max_mesh_size: float = Field(
        5,
        description="Maximum triangle size in meters (when not using adaptive meshing)",
    )

    adaptive_mesh: bool = Field(
        False,
        description="Whether to use adaptive meshing for the terrain surface mesh (dynmically adjusts mesh density based on terrain complexity)",
    )

    error_threshold: float = Field(
        0.5, description="Maximum allowed error (in meters) for adaptive meshing"
    )

    smoothing: int = Field(
        3, description="Number of smoothing iterations to apply to the terrain mesh"
    )
    mesher: Optional[Literal["auto", "dtcc_mesher", "triangle"]] = Field(
        None,
        description="2D meshing backend to use for the terrain triangulation",
    )

    remove_outliers: bool = Field(
        True, description="Whether to remove global outliers from the terrain raster"
    )
    remove_outlier_threshold: float = Field(
        3.0, description="Threshold for outlier removal"
    )
    format: Optional[Literal["tif", "obj", "stl"]] = Field(
        None, description="Output file format"
    )


class TerrainSurfaceMeshDataset(DatasetDescriptor):
    name = "terrain_surface_mesh"
    title = "Terrain Surface Mesh"
    description = (
        "Terrain raster or triangular surface mesh derived from point cloud "
        "ground data for the requested bounds."
    )
    ArgsModel = TerrainSurfaceMeshArgs
    data_category = "derived"
    result_kind = "mesh"
    python_return_type = "dtcc_core.model.Mesh | dtcc_core.model.Raster"
    provider = [
        provider_entry("lantmateriet", role="source_provider"),
        provider_entry("dtcc-platform", role="processor"),
    ]
    source = [
        {
            "name": "DTCC/Lantmäteriet point cloud backend/cache",
            "role": "upstream_dataset",
            "source_terms_status": "requires_review",
        }
    ]
    license = (
        "Requires review: derived from upstream point cloud data; verify "
        "source terms before redistribution."
    )
    collection_period = (
        "Requires review: inherits point cloud acquisition date/cache vintage, "
        "which is not currently surfaced in the terrain result."
    )
    default_crs = "EPSG:3006"
    lod = "Terrain raster or terrain surface mesh; no building LoD"
    data_types = ["terrain", "point_cloud_derived", "raster", "surface_mesh", "mesh"]
    geographic_coverage = "Sweden, constrained by requested bounds and source coverage"
    update_frequency = "derived on demand from upstream source data"
    processing_steps = [
        "Download point cloud data for the requested bounds in EPSG:3006",
        "Optionally remove global point-cloud outliers using remove_outlier_threshold",
        "For format='tif', build a terrain raster with raster_resolution cell size",
        "For adaptive_mesh=True, build an adaptive terrain mesh using error_threshold and raster_resolution",
        "Otherwise build a terrain surface mesh using max_mesh_size, smoothing, and the selected mesher backend",
        "Return the native mesh/raster object or serialize the requested terrain format",
    ]
    derived_from = [
        {
            "name": "point_cloud",
            "relationship": "terrain elevation source",
            "source_terms_status": "requires_review",
        }
    ]
    presentation_headline = "Point-Cloud Terrain Surface"
    presentation_summary = (
        "A terrain-only raster or triangular mesh built from point cloud data."
    )
    presentation_narrative = [
        {
            "heading": "What you are seeing",
            "body": (
                "The result represents terrain elevation over the requested "
                "bounds. It does not include buildings, vegetation semantics, "
                "or simulation boundary conditions."
            ),
        },
        {
            "heading": "How it is made",
            "body": (
                "The dataset downloads point cloud data, optionally removes "
                "global outliers, and then builds either a terrain raster or a "
                "triangular surface mesh with the requested size and mesher settings."
            ),
        },
        {
            "heading": "Limitations",
            "body": (
                "Mesh usefulness depends on point-cloud coverage, classification "
                "quality, raster resolution, smoothing, and later validation for "
                "the intended analysis."
            ),
        },
    ]
    key_points = [
        "format='tif' returns raster bytes instead of a Mesh object",
        "adaptive_mesh=True uses error_threshold rather than max_mesh_size",
        "remove_outliers=True applies a global point-cloud outlier filter before terrain construction",
        "mesher selects the 2D terrain triangulation backend when supported",
    ]
    presentation_legend = {
        "title": "Terrain outputs",
        "entries": [
            {"label": "Raster cell", "meaning": "terrain elevation sample in meters"},
            {"label": "Triangle", "meaning": "terrain surface facet"},
            {
                "label": "Mesh edge size",
                "meaning": "controlled by max_mesh_size or adaptive error settings",
            },
        ],
    }
    view_hints = {
        "preferred_geometry": "terrain_surface_mesh",
        "default_crs": "EPSG:3006",
        "mesh_role": "terrain_visualization_or_preprocessing",
        "z_units": "meters",
    }
    presentation_warnings = [
        "Source and license terms inherit upstream point-cloud review status.",
        "Terrain mesh quality is not automatically certified for FEM, CFD, or hydrology use.",
        "Outlier removal and smoothing can change terrain detail.",
    ]
    presentation_limitations = [
        "Point-cloud acquisition date, density, and classification quality are not currently surfaced in the result.",
        "The dataset does not add building, road, vegetation, or water semantics.",
        "No mesh-quality threshold is enforced beyond the selected builder parameters.",
    ]

    def build(self, args: TerrainSurfaceMeshArgs):
        progress_phases = {
            "download_pointcloud": 0.40,
            "remove_outliers": 0.10,
            "build_terrain": 0.40,
            "export": 0.10,
        }
        with ProgressTracker(total=1.0, phases=progress_phases) as progress:
            bounds = self.parse_bounds(args.bounds)

            with progress.phase(
                "download_pointcloud", "Downloading point cloud data..."
            ):
                pc = dtcc_core.io.data.download_pointcloud(bounds=bounds)

            with progress.phase(
                "remove_outliers",
                (
                    f"Removing outliers (threshold={args.remove_outlier_threshold})..."
                    if args.remove_outliers
                    else "Skipping outlier removal"
                ),
            ):
                if args.remove_outliers:
                    pc = pc.remove_global_outliers(args.remove_outlier_threshold)

            with progress.phase(
                "build_terrain",
                (
                    "Building terrain raster..."
                    if args.format == "tif"
                    else "Building terrain surface mesh..."
                ),
            ):
                if args.format == "tif":
                    result = dtcc_core.builder.build_terrain_raster(
                        pc, cell_size=args.raster_resolution
                    )
                else:
                    if args.adaptive_mesh:
                        result = dtcc_core.builder.adaptive_terrain_mesh(
                            pc, args.error_threshold, args.raster_resolution
                        )
                    else:
                        result = dtcc_core.builder.build_terrain_surface_mesh(
                            pc,
                            max_mesh_size=args.max_mesh_size,
                            smoothing=args.smoothing,
                            mesher=args.mesher,
                        )

            with progress.phase(
                "export",
                (
                    f"Exporting terrain to {args.format}..."
                    if args.format
                    else "Preparing terrain result..."
                ),
            ):
                if args.format == "tif":
                    return self.export_to_bytes(result, "tif")
                elif args.format is None:
                    return result
                elif args.format in ("obj", "stl"):
                    return self.export_to_bytes(result, args.format)
