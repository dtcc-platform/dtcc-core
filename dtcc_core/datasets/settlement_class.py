from typing import Literal, Optional

import dtcc_core
from dtcc_core.io.data.ghsl import (
    SETTLEMENT_CLASS_CITATION,
    SETTLEMENT_CLASS_DEFAULT_EPOCH,
    SETTLEMENT_CLASS_NODATA,
    SETTLEMENT_CLASSES,
)
from dtcc_core.model import Raster
from pydantic import Field

from .built_up import Epoch
from .dataset import DatasetBaseArgs, DatasetDescriptor
from .providers import provider_entry


class SettlementClassArgs(DatasetBaseArgs):
    """Arguments for the GHSL settlement class dataset."""

    epoch: Epoch = Field(
        SETTLEMENT_CLASS_DEFAULT_EPOCH,
        description="Epoch from 1975 to 2030 in steps of 5; 2025 and 2030 are projections.",
    )
    format: Optional[Literal["tif"]] = Field(
        None, description="Output format (tif for GeoTIFF bytes)"
    )


class SettlementClassDataset(DatasetDescriptor):
    """Settlement class dataset, registered as ``settlement_class``.

    GHSL settlement class (Degree of Urbanisation) for the requested bounds and
    epoch, returned as a 1 km Raster of class codes in EPSG:3006.

    Call it with keyword arguments defined by ``SettlementClassArgs``.
    """

    name = "settlement_class"
    title = "Settlement Class"
    description = (
        "GHSL settlement class (Degree of Urbanisation) for the requested bounds "
        "and epoch (1975-2030), returned as a 1 km Raster of class codes in "
        "EPSG:3006."
    )
    ArgsModel = SettlementClassArgs
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
            "name": "GHS-SMOD R2023A V2.0, global 1 km Mollweide grids",
            "role": "source_provider",
            "service": "jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL",
            "doi": "10.2905/A0DF7A6F-49DE-46EA-9BDE-563437A6E2BA",
            "source_terms_status": "requires_review",
        }
    ]
    license = (
        "CC BY 4.0 (Commission Decision 2011/833/EU). Credit \"© European Union\" "
        f"and cite: {SETTLEMENT_CLASS_CITATION}. Indicate changes such as clipping "
        "and reprojection. Requires review before redistribution."
    )
    collection_period = (
        "Epochs 1975-2030 every 5 years, derived from GHSL population and built-up "
        "grids; 2025 and 2030 are projections"
    )
    data_types = ["raster", "settlement_class", "categorical"]
    geographic_coverage = "Global, constrained by requested bounds"
    update_frequency = "GHSL releases every few years (current: R2023A)"
    processing_steps = [
        "Read the requested window of the global 1 km Mollweide GeoTIFF directly from the remote zip",
        "Assign each 1 km EPSG:3006 cell the class of the nearest source cell",
        f"Mark cells without data as {SETTLEMENT_CLASS_NODATA}",
        "Return a Raster or serialize it as GeoTIFF bytes",
    ]
    presentation_headline = "Settlement Class Over Time"
    presentation_summary = (
        "Whether each square kilometre is an urban centre, a suburb or rural land, "
        "using the EU Degree of Urbanisation, every five years from 1975 to 2030."
    )
    presentation_narrative = [
        {
            "heading": "What you are seeing",
            "body": (
                "Each 1 km cell holds a Degree of Urbanisation class code, from 30 "
                "(urban centre) down to 11 (very low density rural), with 10 for water."
            ),
        },
        {
            "heading": "How to interpret it",
            "body": (
                "Classes come from population and built-up density in clusters of "
                "cells, following the method recommended by the UN Statistical "
                "Commission. Comparing epochs shows where settlements grew."
            ),
        },
        {
            "heading": "Limitations",
            "body": (
                "At 1 km the classes describe districts, not streets. Epochs between "
                "observations are modelled, and 2025 and 2030 are projections."
            ),
        },
    ]
    key_points = [
        "Degree of Urbanisation class per 1 km cell",
        "Epochs 1975-2030 every 5 years; 2020 by default",
        "Same method worldwide, recommended by the UN Statistical Commission",
    ]
    presentation_legend = {
        "title": "Settlement class",
        "entries": [
            {"label": str(code), "meaning": name} for code, name in SETTLEMENT_CLASSES.items()
        ]
        + [{"label": str(SETTLEMENT_CLASS_NODATA), "meaning": "No data"}],
    }
    view_hints = {
        "preferred_geometry": "raster",
        "default_crs": "EPSG:3006",
        "value_type": "categorical",
        "classes": {str(code): name for code, name in SETTLEMENT_CLASSES.items()},
    }
    presentation_warnings = [
        "GHSL credit and citation must accompany redistributed data.",
        "Epochs other than the observed years are modelled or projected.",
    ]
    presentation_limitations = [
        "The 1 km grid describes districts, not individual streets or blocks.",
        "Cells near the requested edge take the class of the nearest source cell.",
    ]

    def build(self, args: SettlementClassArgs):
        bounds = self.parse_bounds(args.bounds)
        raster: Raster = dtcc_core.io.data.download_settlement_class(
            bounds=bounds, epoch=args.epoch
        )
        if args.format == "tif":
            return self.export_to_bytes(raster, "tif")
        return raster


__all__ = ["SettlementClassArgs", "SettlementClassDataset"]
