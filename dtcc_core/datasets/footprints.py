import dtcc_core
from dtcc_core.model import Building, City, FootprintCollection
from typing import List, Literal, Optional
from pydantic import Field

from .dataset import DatasetDescriptor, DatasetBaseArgs
from .providers import provider_entry
from dtcc_core.common.progress import ProgressTracker


_SOURCE_TO_PROVIDER = {
    "LM": "dtcc",
    "OSM": "OSM",
}


class FootprintsArgs(DatasetBaseArgs):
    source: Literal["OSM", "LM"] = Field(
        "LM",
        description=(
            "Data source for building footprints. LM uses the DTCC/Lantmäteriet "
            "footprint backend/cache; OSM uses the OpenStreetMap/Overpass path."
        ),
    )
    smallest_building_size: float = Field(
        0.0, description="Smallest building size to include (in square meters)"
    )
    calculate_heights: bool = Field(
        False, description="Whether to calculate building heights from point cloud"
    )
    format: Optional[Literal["geojson", "gpkg", "shp.zip"]] = Field(
        None, description="Output file format"
    )
    crs: Optional[str] = Field(
        None,
        description=(
            "Output CRS for serialized formats (e.g. 'EPSG:3006'). If None, "
            "format conventions apply: GeoJSON is reprojected to EPSG:4326."
        ),
    )


class FootprintsDataset(DatasetDescriptor):
    name = "building_footprints"
    title = "Building Footprints"
    description = (
        "Building footprint polygons for the requested bounds, sourced from "
        "the selected live/cache provider and returned as table-alignment "
        "context geometry."
    )
    ArgsModel = FootprintsArgs
    data_category = "raw"
    result_kind = "building_footprints"
    python_return_type = "dtcc_core.model.FootprintCollection"
    provider = [
        provider_entry("lantmateriet", role="source_provider"),
        provider_entry("openstreetmap", role="source_provider"),
        provider_entry("dtcc-platform", role="processor"),
    ]
    source = [
        {
            "name": "DTCC/Lantmäteriet footprint backend/cache",
            "role": "source_provider",
            "selected_when": 'source="LM"',
            "source_terms_status": "requires_review",
        },
        {
            "name": "OpenStreetMap building footprints via Overpass",
            "role": "source_provider",
            "selected_when": 'source="OSM"',
            "source_terms_status": "requires_review",
        },
    ]
    license = (
        "Requires review: redistribution terms depend on the selected source "
        "(`LM`/DTCC-Lantmäteriet backend or `OSM`/OpenStreetMap)."
    )
    collection_period = (
        "Requires review: temporal coverage and update cadence depend on the "
        "selected provider/cache tile and are not currently surfaced in the result."
    )
    default_crs = "EPSG:3006"
    data_types = ["vector", "polygon", "building_footprints"]
    geographic_coverage = "Sweden, constrained by requested bounds and source coverage"
    update_frequency = "varies by selected upstream provider"
    processing_steps = [
        "Map source='LM' to the DTCC/Lantmäteriet footprint backend/cache or source='OSM' to the OpenStreetMap/Overpass path",
        "Download building footprints intersecting the requested bounds in EPSG:3006",
        "Filter footprints below smallest_building_size when requested",
        "Optionally download point cloud data and estimate building heights",
        "Return a FootprintCollection or serialize the selected vector format",
    ]
    derived_from = [
        {
            "name": "Selected upstream building footprint provider",
            "relationship": "source chosen by request parameter",
        },
        {
            "name": "Point cloud data",
            "relationship": "optional height enrichment when calculate_heights=True",
        },
    ]
    presentation_headline = "Building Footprint Alignment Layer"
    presentation_summary = (
        "Source building outlines for the requested area, useful as a real-world "
        "context layer and tangible-table alignment overlay."
    )
    presentation_narrative = [
        {
            "heading": "What you are seeing",
            "body": (
                "Each polygon is an upstream building footprint clipped or "
                "filtered to the requested bounds. For table use, export as "
                "GeoJSON with crs='EPSG:3006' so coordinates stay in meters."
            ),
        },
        {
            "heading": "How to interpret it",
            "body": (
                "The layer is best used as context geometry: footprints should "
                "align with the printed buildings and help spot projection, "
                "bounds, or CRS mismatches."
            ),
        },
        {
            "heading": "Limitations",
            "body": (
                "Coverage, update date, geometry quality, and redistribution "
                "terms depend on the selected provider. Height enrichment is "
                "optional and comes from point cloud processing, not the raw "
                "footprint source."
            ),
        },
    ]
    key_points = [
        "source='LM' uses the DTCC/Lantmäteriet footprint backend/cache",
        "source='OSM' uses the OpenStreetMap/Overpass footprint path",
        "GeoJSON defaults to WGS84 unless crs='EPSG:3006' is passed for table alignment",
        "calculate_heights=True adds point-cloud-derived height attributes",
    ]
    presentation_legend = {
        "title": "Building footprints",
        "entries": [
            {
                "label": "Polygon outline",
                "meaning": "building footprint from the selected provider",
            },
            {
                "label": "Height attribute",
                "meaning": "optional point-cloud-derived value when calculate_heights=True",
            },
        ],
    }
    view_hints = {
        "preferred_geometry": "polygons",
        "table_role": "alignment_context",
        "default_crs_for_table": "EPSG:3006",
        "default_style": "footprint_outline",
    }
    presentation_warnings = [
        "Source and license terms require review before redistributing generated footprint packages.",
        "GeoJSON table exports must pass crs='EPSG:3006'; the dataset-level GeoJSON default is WGS84.",
    ]
    presentation_limitations = [
        "Provider coverage and update dates are not currently surfaced per feature.",
        "Geometry can be missing, generalized, stale, or incomplete for some buildings.",
        "Height enrichment, when enabled, depends on point cloud coverage and processing assumptions.",
    ]

    def build(self, args: FootprintsArgs):
        progress_phases = {
            "download_footprints": 0.20,
            "download_pointcloud": 0.25,
            "compute_building_heights": 0.40,
            "export": 0.15,
        }

        with ProgressTracker(total=1.0, phases=progress_phases) as progress:
            bounds = self.parse_bounds(args.bounds)
            city = City()
            city.bounds = bounds
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

            if args.calculate_heights:
                with progress.phase(
                    "download_pointcloud", "Downloading point cloud data..."
                ):
                    city.download_pointcloud()
                with progress.phase(
                    "compute_building_heights", "Computing building heights..."
                ):
                    city.building_heights_from_pointcloud(keep_roof_points=False)

            else:
                progress.update(increment=65)

            with progress.phase("export", "Exporting footprints..."):
                if args.format is None:
                    return city.buildings
                else:
                    return self.export_to_bytes(
                        city,
                        args.format,
                        save_callable=dtcc_core.io.footprints.save,
                        output_crs=args.crs,
                    )

    def prepare_result(self, result, validated_args: FootprintsArgs):
        if validated_args.format is None and isinstance(result, list):
            return FootprintCollection.from_buildings(result)
        return result


def _provider_for_source(source: str) -> str:
    try:
        return _SOURCE_TO_PROVIDER[source]
    except KeyError as exc:
        raise ValueError(f"Unsupported building footprint source: {source!r}") from exc


def _filter_small_buildings(
    buildings: List[Building],
    *,
    min_area: float,
) -> List[Building]:
    filtered: List[Building] = []
    for building in buildings:
        footprint = building.footprint()
        if footprint is None:
            continue
        polygon = footprint.to_polygon(simplify=0.0)
        if not polygon.is_empty and polygon.area >= min_area:
            filtered.append(building)
    return filtered
