from typing import Literal, Optional

import dtcc_core
from dtcc_core.io.data.ghsl import BUILT_UP_CITATION, BUILT_UP_DEFAULT_EPOCH
from dtcc_core.model import Raster
from pydantic import Field

from .dataset import DatasetBaseArgs, DatasetDescriptor
from .providers import provider_entry


Epoch = Literal[1975, 1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020, 2025, 2030]


class BuiltUpArgs(DatasetBaseArgs):
    """Arguments for the GHSL built-up surface dataset."""

    epoch: Epoch = Field(
        BUILT_UP_DEFAULT_EPOCH,
        description="Epoch from 1975 to 2030 in steps of 5; 2025 and 2030 are projections.",
    )
    component: Literal["total", "non_residential"] = Field(
        "total",
        description="All built-up surface, or only the non-residential part.",
    )
    format: Optional[Literal["tif"]] = Field(
        None, description="Output format (tif for GeoTIFF bytes)"
    )


class BuiltUpDataset(DatasetDescriptor):
    """Built-up surface dataset, registered as ``built_up``.

    GHSL built-up surface for the requested bounds and epoch, returned as a
    100 m Raster in EPSG:3006 with square metres of built-up surface per cell.

    Call it with keyword arguments defined by ``BuiltUpArgs``.
    """

    name = "built_up"
    title = "Built-up Surface"
    description = (
        "GHSL built-up surface for the requested bounds and epoch (1975-2030), "
        "returned as a 100 m Raster in EPSG:3006 with square metres of built-up "
        "surface per cell."
    )
    ArgsModel = BuiltUpArgs
    data_category = "raw"
    result_kind = "raster"
    python_return_type = "dtcc_core.model.Raster"
    default_crs = "EPSG:3006"
    provider = [
        provider_entry("jrc", role="source_provider"),
        provider_entry("dtcc-platform", role="processor"),
    ]
    source = [
        {
            "name": "GHS-BUILT-S R2023A, global 100 m Mollweide grids",
            "role": "source_provider",
            "service": "jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL",
            "doi": "10.2905/9F06F36F-4B11-47EC-ABB0-4F8B7B1D72EA",
            "source_terms_status": "requires_review",
        }
    ]
    license = (
        "CC BY 4.0 (Commission Decision 2011/833/EU). Credit \"© European Union\" "
        f"and cite: {BUILT_UP_CITATION}. Indicate changes such as clipping and "
        "reprojection. Requires review before redistribution."
    )
    collection_period = (
        "Epochs 1975-2030 every 5 years, interpolated from satellite observations "
        "around 1975, 1990, 2000, 2014 and 2018; 2025 and 2030 are projections"
    )
    data_types = ["raster", "built_up_surface"]
    geographic_coverage = "Global, constrained by requested bounds"
    update_frequency = "GHSL releases every few years (current: R2023A)"
    processing_steps = [
        "Read the requested window of the global 100 m Mollweide GeoTIFF directly from the remote zip",
        "Convert built-up m² per cell to a built-up share",
        "Average the share onto the 100 m EPSG:3006 grid, which preserves area totals",
        "Convert back to m² per EPSG:3006 cell",
        "Return a Raster or serialize it as GeoTIFF bytes",
    ]
    presentation_headline = "Built-up Surface Over Time"
    presentation_summary = (
        "How much of each 100 m cell is covered by buildings, for every five years "
        "from 1975 to 2030, from the EU Global Human Settlement Layer."
    )
    presentation_narrative = [
        {
            "heading": "What you are seeing",
            "body": (
                "Each 100 m cell holds the area covered by buildings in square "
                "metres, from 0 to 10,000. Comparing epochs shows where a city "
                "has grown or densified."
            ),
        },
        {
            "heading": "How to interpret it",
            "body": (
                "Built-up surface counts building roofs only, not roads or paved "
                "ground, which an imperviousness layer covers. "
                "The non-residential component isolates commercial, industrial "
                "and public buildings."
            ),
        },
        {
            "heading": "Limitations",
            "body": (
                "Values between observed years are interpolated, and 2025 and 2030 "
                "are projections. At 100 m the grid shows districts, not "
                "individual buildings; use building_footprints for those."
            ),
        },
    ]
    key_points = [
        "Square metres of built-up surface per 100 m cell",
        "Epochs 1975-2030 every 5 years; 2020 by default",
        "Total or non-residential built-up surface",
        "Reprojection to EPSG:3006 preserves area totals",
    ]
    presentation_legend = {
        "title": "Built-up surface",
        "entries": [
            {"label": "0", "meaning": "No built-up surface in the cell"},
            {"label": "10000", "meaning": "Cell fully built up (m² per 100 m cell)"},
            {"label": "NaN", "meaning": "No data"},
        ],
    }
    view_hints = {
        "preferred_geometry": "raster",
        "default_crs": "EPSG:3006",
        "value_range": [0, 10000],
        "unit": "m2 per cell",
    }
    presentation_warnings = [
        "GHSL credit and citation must accompany redistributed data.",
        "Epochs other than the observed years are interpolated or projected.",
    ]
    presentation_limitations = [
        "The 100 m grid cannot resolve individual buildings.",
        "Built-up surface excludes roads, squares and other paved ground.",
    ]

    def build(self, args: BuiltUpArgs):
        bounds = self.parse_bounds(args.bounds)
        raster: Raster = dtcc_core.io.data.download_built_up(
            bounds=bounds, epoch=args.epoch, component=args.component
        )
        if args.format == "tif":
            return self.export_to_bytes(raster, "tif")
        return raster


__all__ = ["BuiltUpArgs", "BuiltUpDataset"]
