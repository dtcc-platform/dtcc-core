import dtcc_core
from dtcc_core.model import BuildingCollection, City
from typing import Literal, Optional
from pydantic import Field

from .dataset import DatasetDescriptor, DatasetBaseArgs
from .footprints import _filter_small_buildings, _provider_for_source
from .providers import provider_entry
from dtcc_core.common.progress import ProgressTracker, report_progress


class BuildingArgs(DatasetBaseArgs):
    source: Literal["OSM", "LM"] = Field(
        "LM",
        description=(
            "Data source for building footprints. LM uses the DTCC/Lantmäteriet "
            "footprint backend/cache; OSM uses the OpenStreetMap/Overpass path."
        ),
    )
    smallest_building_size: float = Field(
        15.0, description="Smallest building size to include (in square meters)"
    )
    place_on_zero: bool = Field(
        False, description="Whether to place buildings on Z=0 plane"
    )
    format: Optional[Literal["obj", "stl"]] = Field(
        None, description="Output file format"
    )


class BuildingDataset(DatasetDescriptor):
    name = "buildings"
    title = "LoD1 Buildings"
    description = (
        "LoD1 building solids derived from source building footprints and point "
        "cloud height estimates for the requested bounds."
    )
    ArgsModel = BuildingArgs
    data_category = "derived"
    result_kind = "building_collection"
    python_return_type = "dtcc_core.model.BuildingCollection"
    provider = [
        provider_entry("lantmateriet", role="source_provider"),
        provider_entry("openstreetmap", role="source_provider"),
        provider_entry("dtcc-platform", role="processor"),
    ]
    source = [
        {
            "name": "DTCC/Lantmäteriet point cloud backend/cache",
            "role": "upstream_dataset",
            "source_terms_status": "requires_review",
        },
        {
            "name": "Selected building footprint provider",
            "role": "upstream_dataset",
            "selected_by": "source",
            "source_terms_status": "requires_review",
        },
    ]
    license = (
        "Requires review: derived from upstream point cloud and footprint data; "
        "verify all source terms before redistribution."
    )
    collection_period = (
        "Requires review: inherits point cloud acquisition date and footprint "
        "provider/cache vintage, which are not currently surfaced in the result."
    )
    default_crs = "EPSG:3006"
    lod = "LoD1 block buildings extruded from footprint geometry"
    data_types = ["city_model", "buildings", "lod1", "mesh"]
    geographic_coverage = "Sweden, constrained by requested bounds and source coverage"
    update_frequency = "derived on demand from upstream source data"
    processing_steps = [
        "Download point cloud data for the requested bounds in EPSG:3006",
        "Download selected building footprints for the requested bounds",
        "Filter footprints below smallest_building_size when requested",
        "Remove global point-cloud elevation outliers using a 3.0 standard-deviation threshold",
        "Build a terrain raster from ground-classified point-cloud returns",
        "Extract roof points inside each footprint",
        "Estimate building heights from roof points and terrain",
        "Extrude footprint geometry into LoD1 block buildings",
        "Return a BuildingCollection or serialize merged LoD1 meshes",
    ]
    derived_from = [
        {
            "name": "point_cloud",
            "relationship": "terrain and roof-height estimation",
            "source_terms_status": "requires_review",
        },
        {
            "name": "building_footprints",
            "relationship": "footprint outlines and source IDs",
            "source_terms_status": "requires_review",
        },
    ]
    presentation_headline = "LoD1 Building Blocks"
    presentation_summary = (
        "Blocky building solids that combine source footprints with estimated "
        "terrain-relative heights from point cloud data."
    )
    presentation_narrative = [
        {
            "heading": "What you are seeing",
            "body": (
                "Each object is a LoD1 building: a footprint extruded to an "
                "estimated height. Roof shape, facade detail, and semantic parts "
                "are not reconstructed."
            ),
        },
        {
            "heading": "How it is made",
            "body": (
                "The pipeline downloads footprints and point cloud data, filters "
                "small footprints, estimates terrain, extracts roof points, and "
                "extrudes each footprint into a block model."
            ),
        },
        {
            "heading": "Limitations",
            "body": (
                "Height quality depends on source point density, footprint quality, "
                "terrain estimation, and roof-point coverage. Treat outputs as "
                "generalized context geometry, not surveyed BIM."
            ),
        },
    ]
    key_points = [
        "source='LM' uses the DTCC/Lantmäteriet footprint backend/cache",
        "source='OSM' uses OpenStreetMap/Overpass footprints",
        "Heights are estimated from point-cloud roof points relative to terrain",
        "smallest_building_size filters footprints before height estimation",
        "place_on_zero affects exported merged meshes only",
    ]
    presentation_legend = {
        "title": "LoD1 buildings",
        "entries": [
            {"label": "Footprint", "meaning": "source building outline"},
            {"label": "Height", "meaning": "point-cloud-derived LoD1 extrusion height"},
            {"label": "Merged mesh", "meaning": "OBJ/STL export artifact when requested"},
        ],
    }
    view_hints = {
        "preferred_geometry": "solid_buildings",
        "default_crs": "EPSG:3006",
        "level_of_detail": "LoD1",
    }
    presentation_warnings = [
        "Source and license terms inherit upstream point cloud and footprint review status.",
        "LoD1 output does not contain roof geometry, facade detail, or BIM semantics.",
        "Height estimates can be wrong when roof points, terrain, or footprint geometry are incomplete.",
    ]
    presentation_limitations = [
        "Acquisition date and point density are not currently surfaced per building.",
        "Small or complex buildings can be filtered, merged, generalized, or missed.",
        "Exported OBJ/STL meshes are visualization geometry, not validated simulation meshes.",
    ]

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
                city.download_footprints(provider=_provider_for_source(args.source))
                if args.smallest_building_size > 0:
                    city.replace_buildings(
                        _filter_small_buildings(
                            city.buildings,
                            min_area=args.smallest_building_size,
                        )
                    )

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

    def prepare_result(self, result, validated_args: BuildingArgs):
        if validated_args.format is None and isinstance(result, list):
            return BuildingCollection(result)
        return result
