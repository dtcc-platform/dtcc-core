import dtcc_core
from dtcc_core.model import Building, City, FootprintCollection
from typing import List, Literal, Optional
from pydantic import Field

from .dataset import DatasetDescriptor, DatasetBaseArgs
from .providers import provider_entry
from dtcc_core.common.progress import ProgressTracker


class FootprintsArgs(DatasetBaseArgs):
    source: Literal["OSM", "LM"] = Field(
        "LM", description="Data source for building footprints"
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
    description = "Building footprints"
    ArgsModel = FootprintsArgs
    data_category = "raw"
    result_kind = "building_footprints"
    python_return_type = "dtcc_core.model.FootprintCollection"
    provider = [
        provider_entry("lantmateriet"),
        provider_entry("openstreetmap"),
    ]
    source = ["Building footprint source selected by the dataset request"]
    license = "Review selected upstream footprint source terms before redistribution."
    default_crs = "EPSG:3006"
    geographic_coverage = "Sweden, constrained by requested bounds and source coverage"
    update_frequency = "varies by selected upstream provider"
    processing_steps = [
        "Download building footprints for requested bounds",
        "Optionally estimate building heights from point cloud data",
    ]
    presentation_summary = (
        "Building footprints for the requested area, with optional height "
        "enrichment from point cloud data."
    )

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
                city.download_footprints()

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
