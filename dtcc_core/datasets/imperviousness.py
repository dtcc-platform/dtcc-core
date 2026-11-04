from typing import Literal

from pydantic import Field

import dtcc_core
from dtcc_core.io.data.imperviousness import (
    IMPERVIOUSNESS_DEFAULT_YEAR,
    IMPERVIOUSNESS_DOI,
    empty_imperviousness_raster,
)
from dtcc_core.model import Raster

from ..common import warning
from .context import attach_dataset_context
from .dataset import DatasetBaseArgs, DatasetDescriptor, DatasetUpstreamError
from .providers import provider_entry


class ImperviousnessArgs(DatasetBaseArgs):
    """Arguments for the Copernicus HRL Imperviousness Density dataset."""

    year: Literal[2018, 2021, 2024] = Field(
        IMPERVIOUSNESS_DEFAULT_YEAR,
        description="Reference year of the Imperviousness Density layer.",
    )
    format: Literal["tif"] | None = Field(
        None, description="Output format (tif for GeoTIFF bytes)"
    )


class ImperviousnessDataset(DatasetDescriptor):
    """Imperviousness dataset, registered as ``imperviousness``.

    Copernicus HRL Imperviousness Density for the requested bounds, returned
    as a 10 m Raster in EPSG:3006 with the sealed share of each cell in percent.

    Call it with keyword arguments defined by ``ImperviousnessArgs``.
    """

    name = "imperviousness"
    title = "Imperviousness Density"
    description = (
        "Copernicus HRL Imperviousness Density for the requested bounds, returned "
        "as a 10 m Raster in EPSG:3006 with the sealed share of each cell in percent."
    )
    ArgsModel = ImperviousnessArgs
    data_category = "raw"
    result_kind = "raster"
    python_return_type = "dtcc_core.model.Raster"
    default_crs = "EPSG:3006"
    provider = [
        provider_entry("copernicus-land", role="source_provider"),
        provider_entry("dtcc-platform", role="processor"),
    ]
    source = [
        {
            "name": "EEA ArcGIS ImageServer HRL_ImperviousnessDensity_2018",
            "role": "source_provider",
            "service": "image.discomap.eea.europa.eu",
            "supported_years": [2018],
            "doi": IMPERVIOUSNESS_DOI[2018].removeprefix("https://doi.org/"),
            "source_terms_status": "requires_review",
        },
        {
            "name": "CLMS GeoServer WMS nvlcc_2021 and nvlcc_2024",
            "role": "source_provider",
            "service": "geoserver.geoville.com",
            "supported_years": [2021, 2024],
            "doi": [
                IMPERVIOUSNESS_DOI[2021].removeprefix("https://doi.org/"),
                IMPERVIOUSNESS_DOI[2024].removeprefix("https://doi.org/"),
            ],
            "source_terms_status": "requires_review",
        },
    ]
    license = (
        "Copernicus full, free and open data policy (Regulation (EU) 2021/696). "
        "Attribution: \"Generated using European Union's Copernicus Land "
        'Monitoring Service information; <DOI>", with the DOI of the requested '
        "year. State any modification and do not "
        "imply EU endorsement. Requires review before redistribution."
    )
    collection_period = "Reference years 2018, 2021 and 2024 (three-yearly)"
    data_types = ["raster", "imperviousness_density"]
    geographic_coverage = (
        "EEA38 countries and the United Kingdom, constrained by requested bounds"
    )
    update_frequency = "Three-yearly reference years"
    processing_steps = [
        (
            "2018: export raw 10 m values in EPSG:3006 from the anonymous EEA "
            "ArcGIS ImageServer"
        ),
        (
            "2021 and 2024: request single-band palette GeoTIFFs in EPSG:3035 from "
            "the CLMS GeoServer WMS, where each pixel index is the density"
        ),
        "Split large areas into tiles within each service's request size limit",
        (
            "Check each image's size, CRS, grid and, for 2021 and 2024, the "
            "palette before reading it"
        ),
        "Set values above 100 (sea, and outside the product for 2018) to NaN",
        "Reproject EPSG:3035 data onto the 10 m EPSG:3006 grid with nearest neighbour",
        (
            "On service failure, raise a typed upstream error with strict_live, "
            "otherwise warn and return an all-NaN, uncached Raster marked partial"
        ),
        "Return a Raster or serialize it as GeoTIFF bytes",
    ]
    presentation_headline = "Imperviousness Density"
    presentation_summary = (
        "The share of sealed ground, such as roofs, asphalt and concrete, in each "
        "10 m cell, mapped from satellite imagery across Europe."
    )
    presentation_narrative = [
        {
            "heading": "What you are seeing",
            "body": (
                "Each 10 m cell holds the percentage of its area covered by "
                "artificial impervious surfaces, from 0 (fully permeable) to 100 "
                "(fully sealed)."
            ),
        },
        {
            "heading": "How to interpret it",
            "body": (
                "Values come from Sentinel-2 imagery classified by the Copernicus "
                "Land Monitoring Service. They describe sealing, not land use: a "
                "green roof or a gravel yard can read as partly permeable."
            ),
        },
        {
            "heading": "Limitations",
            "body": (
                "A 10 m cell cannot resolve narrow paths or single trees over "
                "pavement, and the newest layer may be a few years old."
            ),
        },
    ]
    key_points = [
        "10 m cells with sealed share in percent (0-100)",
        "Reference years 2018, 2021 and 2024; 2024 by default",
        "No account needed; sea is NaN",
        "Useful as a sealing input for land use, runoff and heat studies",
    ]
    presentation_legend = {
        "title": "Imperviousness Density",
        "entries": [
            {"label": "0", "meaning": "No sealed surface in the cell"},
            {"label": "100", "meaning": "Cell fully sealed"},
            {
                "label": "NaN",
                "meaning": "Sea, outside the product (2018 only), or no data",
            },
        ],
    }
    view_hints = {
        "preferred_geometry": "raster",
        "default_crs": "EPSG:3006",
        "value_range": [0, 100],
        "unit": "percent",
    }
    presentation_warnings = [
        "Copernicus attribution and the DOI of the requested year must accompany redistributed data.",
        "Neither service is documented for programmatic use or has an SLA.",
        (
            "2021 and 2024 values are read from WMS map images. The download "
            "checks the image palette and fails with an error if the map "
            "style changes."
        ),
    ]
    presentation_limitations = [
        "The 10 m grid blurs narrow features such as paths and single trees over pavement.",
        "2021 and 2024 are resampled from the EPSG:3035 product grid.",
        (
            "Outside the product footprint, 2021 and 2024 read as 0 % because "
            "the map images draw it the same as unsealed ground; only 2018 "
            "returns NaN there."
        ),
        "A request is limited to 50 km x 50 km (25 million 10 m cells).",
    ]

    def __call__(self, **kwargs):
        args = self.validate(kwargs)
        result, errors = self._realize(args)
        result = self.prepare_result(result, args)
        return attach_dataset_context(result, self._context(args, errors))

    def build(self, args: ImperviousnessArgs):
        return self._realize(args)[0]

    def _export_payload(self, args):
        result, errors = self._realize(args)
        context = self._context(args, errors)
        return result, {"dataset_context": context.model_dump(mode="json")}

    def _realize(self, args: ImperviousnessArgs):
        """Return the Raster or GeoTIFF bytes and any upstream errors."""
        bounds = self.parse_bounds(args.bounds)
        errors: list[DatasetUpstreamError] = []
        try:
            raster: Raster = dtcc_core.io.data.download_imperviousness(
                bounds=bounds, year=args.year
            )
        except DatasetUpstreamError as exc:
            if args.strict_live:
                raise
            warning(
                f"{exc} Returning an all-NaN, uncached imperviousness raster; "
                "pass strict_live=True to raise instead."
            )
            errors.append(exc)
            raster = empty_imperviousness_raster(bounds)
        if args.format == "tif":
            return self.export_to_bytes(raster, "tif"), errors
        return raster, errors

    def _context(self, args: ImperviousnessArgs, errors):
        context = self.create_context(args)
        # A failed tile discards the whole download, so any error means no data.
        health: dict = {"status": "failed" if errors else "complete"}
        self.apply_result_health_metadata(health, upstream_errors=errors)
        warnings = [str(exc) for exc in errors]
        presentation = context.presentation.model_copy(
            update={"warnings": warnings + list(context.presentation.warnings)}
        )
        return context.model_copy(
            update={
                "health": health,
                "warnings": warnings,
                "presentation": presentation,
            }
        )


__all__ = ["ImperviousnessArgs", "ImperviousnessDataset"]
