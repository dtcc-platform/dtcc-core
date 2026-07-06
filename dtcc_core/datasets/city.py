from dtcc_core.model import City
from typing import Literal, Optional
from pydantic import Field

from .dataset import DatasetDescriptor, DatasetBaseArgs
from .footprints import _filter_small_buildings, _provider_for_source
from .providers import provider_entry
from dtcc_core.common.progress import ProgressTracker


class CityArgs(DatasetBaseArgs):
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
    format: Optional[Literal["cityjson", "json"]] = Field(
        None, description="Output file format"
    )


class CityDataset(DatasetDescriptor):
    name = "city"
    title = "City"
    description = (
        "DTCC city model with terrain and LoD1 buildings derived from point cloud "
        "data and selected building footprints for the requested bounds."
    )
    ArgsModel = CityArgs
    data_category = "derived"
    result_kind = "city_model"
    python_return_type = "dtcc_core.model.City"
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
    lod = "Terrain mesh plus LoD1 block buildings"
    data_types = ["city_model", "terrain", "buildings", "lod1", "mesh"]
    geographic_coverage = "Sweden, constrained by requested bounds and source coverage"
    update_frequency = "derived on demand from upstream source data"
    processing_steps = [
        "Download point cloud data for the requested bounds in EPSG:3006",
        "Download selected building footprints for the requested bounds",
        "Filter footprints below smallest_building_size when requested",
        "Build terrain raster and terrain mesh from point cloud data",
        "Estimate building heights and extrude LoD1 building geometry",
        "Return a City object or serialize CityJSON-compatible JSON bytes",
    ]
    derived_from = [
        {
            "name": "point_cloud",
            "relationship": "terrain and building-height estimation",
            "source_terms_status": "requires_review",
        },
        {
            "name": "building_footprints",
            "relationship": "building outlines and source IDs",
            "source_terms_status": "requires_review",
        },
    ]
    presentation_headline = "Terrain and LoD1 City Model"
    presentation_summary = (
        "A DTCC city model that combines terrain from point cloud data with "
        "generalized LoD1 building solids."
    )
    presentation_narrative = [
        {
            "heading": "What you are seeing",
            "body": (
                "The result is a city container with terrain geometry and LoD1 "
                "buildings. Buildings are block extrusions, not detailed roof or "
                "facade reconstructions."
            ),
        },
        {
            "heading": "How it is made",
            "body": (
                "The dataset downloads point cloud and footprint sources, filters "
                "small footprints, builds terrain, estimates heights, and stores "
                "the outputs in a native DTCC City object."
            ),
        },
        {
            "heading": "Limitations",
            "body": (
                "The model is useful for context, visualization, and downstream "
                "meshing preparation, but it is not a surveyed digital twin or "
                "validated simulation mesh."
            ),
        },
    ]
    key_points = [
        "source='LM' uses the DTCC/Lantmäteriet footprint backend/cache",
        "source='OSM' uses OpenStreetMap/Overpass footprints",
        "Terrain and LoD1 buildings are derived on demand",
        "smallest_building_size filters footprints before terrain/building assembly",
        "cityjson and json currently use the same JSON byte export path",
    ]
    presentation_legend = {
        "title": "City model layers",
        "entries": [
            {"label": "Terrain", "meaning": "point-cloud-derived terrain raster/mesh"},
            {"label": "LoD1 building", "meaning": "extruded footprint with estimated height"},
            {"label": "City JSON", "meaning": "serialized city-model package when requested"},
        ],
    }
    view_hints = {
        "preferred_geometry": "city_model",
        "default_crs": "EPSG:3006",
        "terrain": "point_cloud_derived",
        "building_lod": "LoD1",
    }
    presentation_warnings = [
        "Source and license terms inherit upstream point cloud and footprint review status.",
        "City buildings are LoD1 block geometry and do not include detailed roofs or facades.",
        "Terrain and building quality depend on source coverage, classification quality, and cache vintage.",
    ]
    presentation_limitations = [
        "Acquisition date and point density are not currently surfaced in the city result.",
        "The CityJSON-compatible export path currently serializes both cityjson and json requests through JSON bytes.",
        "The generated city model is not automatically validated as an analysis or simulation mesh.",
    ]

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
                city.download_footprints(provider=_provider_for_source(args.source))
                if args.smallest_building_size > 0:
                    city.replace_buildings(
                        _filter_small_buildings(
                            city.buildings,
                            min_area=args.smallest_building_size,
                        )
                    )

            with progress.phase("build_terrain", "Building terrain mesh..."):
                city.build_terrain(build_mesh=True)

            with progress.phase("build_lod1", "Building LOD1 city model..."):
                city.build_lod1_buildings()

            with progress.phase("export", "Exporting city model..."):
                if args.format is None:
                    return city
                else:
                    return self.export_to_bytes(city, "json")
