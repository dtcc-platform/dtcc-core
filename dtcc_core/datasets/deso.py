from pathlib import Path
import tempfile
from typing import Literal, Optional

import dtcc_core
from dtcc_core.model import DeSO
from pydantic import Field

from .dataset import DatasetBaseArgs, DatasetDescriptor
from .providers import provider_entry


class DeSOArgs(DatasetBaseArgs):
    source: Literal["SCB"] = Field("SCB", description="Data source")
    year: Literal[2018, 2025] = Field(2025, description="DeSO geometry vintage")
    statistics: Optional[
        list[Literal["population", "households", "cars", "employment"]]
    ] = Field(
        None,
        description=(
            "Optional SCB DeSO statistics to attach as area-aligned fields "
            "(population, households, cars, employment)."
        ),
    )
    statistics_year: Optional[int] = Field(
        None,
        description="Reference year for attached statistics; defaults to latest supported.",
    )
    format: Optional[Literal["pb", "geojson", "gpkg"]] = Field(
        None,
        description="Output format (pb for protobuf bytes, geojson/gpkg for vector bytes)",
    )


class DeSODataset(DatasetDescriptor):
    name = "deso"
    title = "DeSO"
    description = (
        "Swedish DeSO statistical area boundaries for the requested bounds, with "
        "optional SCB area-aligned statistics fields."
    )
    ArgsModel = DeSOArgs
    data_category = "raw"
    result_kind = "administrative_areas"
    python_return_type = "dtcc_core.model.DeSO"
    default_crs = "EPSG:3006"
    provider = [
        provider_entry("scb", role="source_provider"),
        provider_entry("dtcc-platform", role="processor"),
    ]
    source = [
        {
            "name": "SCB DeSO WFS boundaries",
            "role": "source_provider",
            "supported_years": [2018, 2025],
            "source_terms_status": "requires_review",
        },
        {
            "name": "SCB Statistikdatabasen DeSO statistics",
            "role": "source_provider",
            "selected_when": "statistics is not None",
            "topics": ["population", "households", "cars", "employment"],
            "source_terms_status": "requires_review",
        },
    ]
    license = "Requires review: verify SCB source terms before redistribution."
    collection_period = (
        "Geometry vintages supported by this dataset are 2018 and 2025. Optional "
        "statistics use the latest supported SCB API year per topic unless "
        "statistics_year is provided."
    )
    data_types = ["vector", "administrative_areas", "statistics"]
    geographic_coverage = "Sweden, constrained by requested bounds"
    update_frequency = "selected DeSO geometry/statistics vintage"
    processing_steps = [
        "Download DeSO geometry for requested bounds from the SCB WFS endpoint",
        "Convert SCB polygons to a dtcc-core DeSO object in EPSG:3006",
        "Normalize requested statistics topics when provided",
        "Query or reuse cached SCB Statistikdatabasen values per DeSO code",
        "Attach optional statistics as area-aligned fields",
        "Return a DeSO object or serialize the selected vector/protobuf format",
    ]
    derived_from = [
        {
            "name": "SCB DeSO boundary dataset",
            "relationship": "administrative area geometry",
            "source_terms_status": "requires_review",
        },
        {
            "name": "SCB Statistikdatabasen",
            "relationship": "optional area-aligned statistics fields",
            "source_terms_status": "requires_review",
        },
    ]
    presentation_headline = "Swedish DeSO Areas"
    presentation_summary = (
        "Administrative/statistical area polygons that can carry population, "
        "household, car, or employment attributes from SCB."
    )
    presentation_narrative = [
        {
            "heading": "What you are seeing",
            "body": (
                "Each polygon is a DeSO statistical area intersecting the requested "
                "bounds. Optional statistics are attached as fields aligned by DeSO code."
            ),
        },
        {
            "heading": "How to interpret it",
            "body": (
                "Use DeSO areas for aggregated socioeconomic context. Values describe "
                "the whole area, not individual buildings, parcels, or residents."
            ),
        },
        {
            "heading": "Limitations",
            "body": (
                "Geometry vintage and statistics year can differ. Statistics may be "
                "unavailable for some topic/year combinations and should not be "
                "downscaled without a separate method."
            ),
        },
    ]
    key_points = [
        "Supported geometry vintages are 2018 and 2025",
        "Optional statistics topics are population, households, cars, and employment",
        "statistics_year selects a supported SCB API year for requested topics",
        "Fields are area-aligned by DeSO code",
    ]
    presentation_legend = {
        "title": "DeSO areas",
        "entries": [
            {"label": "Polygon", "meaning": "SCB DeSO statistical area"},
            {"label": "population_total", "meaning": "persons by DeSO area when requested"},
            {"label": "households_total", "meaning": "households by DeSO area when requested"},
            {"label": "cars_total/cars_in_traffic", "meaning": "passenger-car counts when requested"},
            {"label": "employed_residents_total", "meaning": "employed residents when requested"},
        ],
    }
    view_hints = {
        "preferred_geometry": "polygons",
        "default_crs": "EPSG:3006",
        "statistics_topics": ["population", "households", "cars", "employment"],
        "geometry_years": [2018, 2025],
    }
    presentation_warnings = [
        "SCB source and redistribution terms require review before publishing packages.",
        "Statistics describe entire DeSO areas and are not building-level observations.",
        "Geometry vintage and statistics year may differ unless explicitly reviewed for the use case.",
    ]
    presentation_limitations = [
        "Statistics are attached only for requested topics and supported years.",
        "Missing or suppressed SCB values may appear as NaN in fields.",
        "The dataset does not currently expose per-field quality flags beyond source/year metadata.",
    ]

    def build(self, args: DeSOArgs):
        bounds = self.parse_bounds(args.bounds)
        if args.format in ("geojson", "gpkg"):
            if args.statistics:
                deso = dtcc_core.io.data.download_deso(
                    bounds=bounds,
                    year=args.year,
                    source=args.source,
                    statistics=args.statistics,
                    statistics_year=args.statistics_year,
                )
                gdf = deso.to_dataframe()
            else:
                gdf = dtcc_core.io.data.deso.download_deso_geodataframe(
                    bounds=bounds,
                    year=args.year,
                    source=args.source,
                )
            return self._export_gdf_to_bytes(gdf, args.format)

        deso: DeSO = dtcc_core.io.data.download_deso(
            bounds=bounds,
            year=args.year,
            source=args.source,
            statistics=args.statistics,
            statistics_year=args.statistics_year,
        )
        if args.format == "pb":
            return deso.to_proto().SerializeToString()
        return deso

    @staticmethod
    def _export_gdf_to_bytes(gdf, format: str) -> bytes:
        with tempfile.TemporaryDirectory() as tmpdir:
            if format == "geojson":
                path = Path(tmpdir) / "deso.geojson"
                gdf.to_file(path, driver="GeoJSON")
            elif format == "gpkg":
                path = Path(tmpdir) / "deso.gpkg"
                gdf.to_file(path, layer="deso", driver="GPKG")
            else:
                raise ValueError(f"Unsupported DeSO export format: {format}")
            return path.read_bytes()


__all__ = ["DeSOArgs", "DeSODataset"]
