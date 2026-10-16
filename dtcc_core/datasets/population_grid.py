from pathlib import Path
import tempfile
from typing import Literal, Optional

import dtcc_core
from dtcc_core.io.data.population_grid import POPULATION_GRID_DEFAULT_YEAR
from dtcc_core.model import Raster
from pydantic import Field

from .dataset import DatasetBaseArgs, DatasetDescriptor
from .providers import provider_entry


Year = Literal[2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025]
Variable = Literal[
    "total", "men", "women",
    "age_0_4", "age_5_9", "age_10_14", "age_15_19", "age_20_24", "age_25_29",
    "age_30_34", "age_35_39", "age_40_44", "age_45_49", "age_50_54", "age_55_59",
    "age_60_64", "age_65_69", "age_70_74", "age_75_79", "age_80_84", "age_85_89",
    "age_90_94", "age_95_99", "age_100_plus",
]


class PopulationGridArgs(DatasetBaseArgs):
    """Arguments for the SCB 1 km population grid dataset."""

    year: Year = Field(
        POPULATION_GRID_DEFAULT_YEAR,
        description="Reference year; population registered on 31 December.",
    )
    variable: Variable = Field(
        "total",
        description="Population variable for the raster result: total, men, women or a five-year age band.",
    )
    format: Optional[Literal["tif", "gpkg", "geojson"]] = Field(
        None,
        description=(
            "Output format (tif for the selected variable as GeoTIFF bytes; "
            "gpkg/geojson for all variables as 1 km cell polygons)"
        ),
    )


class PopulationGridDataset(DatasetDescriptor):
    """Population grid dataset, registered as ``population_grid``.

    SCB registered population per 1 km cell for the requested bounds, returned
    as a Raster in EPSG:3006 for one variable (total, sex or age band).

    Call it with keyword arguments defined by ``PopulationGridArgs``.
    """

    name = "population_grid"
    title = "Population Grid (1 km)"
    description = (
        "SCB registered population per 1 km cell for the requested bounds, "
        "returned as a Raster in EPSG:3006 for total population, sex or a "
        "five-year age band."
    )
    ArgsModel = PopulationGridArgs
    data_category = "raw"
    result_kind = "raster"
    python_return_type = "dtcc_core.model.Raster"
    default_crs = "EPSG:3006"
    provider = [
        provider_entry("scb", role="source_provider"),
        provider_entry("dtcc-platform", role="processor"),
    ]
    source = [
        {
            "name": "SCB WFS stat:befolkning_1km_{year}",
            "role": "source_provider",
            "service": "geodata.scb.se/geoserver/stat/wfs",
            "supported_years": list(range(2015, 2026)),
            "source_terms_status": "requires_review",
        }
    ]
    license = (
        "CC0 1.0 Universal (SCB open geodata). No attribution is required; "
        "citing \"Source: Statistics Sweden (SCB)\" is recommended."
    )
    collection_period = "Registered population on 31 December, 2015-2025"
    data_types = ["raster", "population_grid", "vector", "polygons"]
    geographic_coverage = "Sweden, constrained by requested bounds"
    update_frequency = "Yearly"
    processing_steps = [
        "Query the SCB WFS for 1 km cells overlapping the requested bounds, in pages",
        "Keep cells whose area overlaps the bounds, not those only touching an edge",
        "Rename SCB fields to total, men, women and age_<from>_<to> variables",
        "Place the selected variable on the 1 km SCB grid, extended outward to whole cells",
        "Set cells that SCB omits (no residents) to 0",
        "Return a Raster or serialize GeoTIFF, GeoPackage or GeoJSON bytes",
    ]
    presentation_headline = "Population per Square Kilometre"
    presentation_summary = (
        "Official registered population per 1 km cell from Statistics Sweden, "
        "with sex and five-year age bands, every year from 2015."
    )
    presentation_narrative = [
        {
            "heading": "What you are seeing",
            "body": (
                "Each 1 km cell holds the number of people registered as living "
                "there on 31 December of the chosen year, for the chosen variable."
            ),
        },
        {
            "heading": "How to interpret it",
            "body": (
                "Counts come from the Swedish population register, not a model. "
                "SCB adjusts small counts to protect privacy, so sex and age "
                "bands may not add up exactly to the total."
            ),
        },
        {
            "heading": "Limitations",
            "body": (
                "In cities a 1 km cell spans several DeSO areas; use the deso "
                "dataset for finer areas. Cells are whole, so the raster can "
                "extend slightly beyond the requested bounds."
            ),
        },
    ]
    key_points = [
        "Persons per 1 km cell from the Swedish population register",
        "Total, men, women and 21 five-year age bands",
        "Yearly 2015-2025; 2025 by default",
        "CC0 licence",
    ]
    presentation_legend = {
        "title": "Population grid",
        "entries": [
            {"label": "Raster value", "meaning": "Persons in the 1 km cell for the selected variable"},
            {"label": "0", "meaning": "No registered residents, or outside Sweden"},
        ],
    }
    view_hints = {
        "preferred_geometry": "raster",
        "alternate_geometry": "polygons",
        "default_crs": "EPSG:3006",
        "unit": "persons per cell",
    }
    presentation_warnings = [
        "Privacy adjustments mean sex and age bands may not sum exactly to the total.",
    ]
    presentation_limitations = [
        "1 km cells are coarser than DeSO areas in cities.",
        "Cells outside Sweden are 0, not missing.",
    ]

    def build(self, args: PopulationGridArgs):
        bounds = self.parse_bounds(args.bounds)
        if args.format in ("gpkg", "geojson"):
            gdf = dtcc_core.io.data.population_grid.download_population_grid_geodataframe(
                bounds=bounds, year=args.year
            )
            return self._export_gdf_to_bytes(gdf, args.format)

        raster: Raster = dtcc_core.io.data.download_population_grid(
            bounds=bounds, year=args.year, variable=args.variable
        )
        if args.format == "tif":
            return self.export_to_bytes(raster, "tif")
        return raster

    @staticmethod
    def _export_gdf_to_bytes(gdf, format: str) -> bytes:
        with tempfile.TemporaryDirectory() as tmpdir:
            if format == "geojson":
                path = Path(tmpdir) / "population_grid.geojson"
                gdf.to_file(path, driver="GeoJSON")
            else:
                path = Path(tmpdir) / "population_grid.gpkg"
                gdf.to_file(path, layer="population_grid", driver="GPKG")
            return path.read_bytes()


__all__ = ["PopulationGridArgs", "PopulationGridDataset"]
