import re

import dtcc_core
from dtcc_core.model import Raster
from pydantic import Field, field_validator
from typing import Literal, Optional

from .dataset import DatasetDescriptor, DatasetBaseArgs
from .providers import provider_display_name, provider_entry


# Keep in step with dtcc_core.io.data.digitalearth.VALID_STYLES.
ImageryStyle = Literal[
    "rgb",
    "false_color",
    "false_color_urban",
    "swir",
    "ndvi",
    "ndmi",
    "ndwi",
    "nbi",
    "msi",
]

_DATE_RE = re.compile(r"^\d{4}-\d{2}-\d{2}$")


class SatelliteImageryArgs(DatasetBaseArgs):
    """Arguments for the Digital Earth Sweden satellite imagery dataset."""

    source: Literal["DES"] = Field(
        "DES",
        description=(
            "Imagery source. DES uses the anonymous Digital Earth Sweden "
            "Sentinel-2 services operated by RISE."
        ),
    )
    date: Optional[str] = Field(
        None,
        description=(
            "Acquisition date as YYYY-MM-DD. When omitted, the most recent "
            "date with at most max_cloud percent cloud cover is used."
        ),
    )
    style: ImageryStyle = Field(
        "rgb",
        description=(
            "Rendering: rgb true colour, false colour composites, or an index "
            "such as ndvi (vegetation) or ndwi (water)."
        ),
    )
    max_cloud: float = Field(
        10.0,
        ge=0.0,
        le=100.0,
        description="Highest cloud cover, in percent, accepted when picking a date.",
    )
    resolution: float = Field(
        10.0,
        gt=0.0,
        description="Ground sample distance in metres; Sentinel-2 is natively 10 m.",
    )
    search_days: int = Field(
        120,
        ge=1,
        description="How many days back to search for a clear acquisition.",
    )
    format: Optional[Literal["tif", "png"]] = Field(
        None, description="Output file format"
    )

    @field_validator("date")
    @classmethod
    def validate_date(cls, value):
        if value is None:
            return value
        if not _DATE_RE.match(value):
            raise ValueError("date must be formatted as YYYY-MM-DD.")
        return value


class SatelliteImageryDataset(DatasetDescriptor):
    """Satellite imagery dataset, registered as ``satellite_imagery``.

    Sentinel-2 imagery from Digital Earth Sweden for the requested bounds,
    defaulting to the most recent acquisition with little cloud cover.

    Call it with keyword arguments defined by ``SatelliteImageryArgs``.
    """

    name = "satellite_imagery"
    title = "Satellite Imagery"
    description = (
        "Sentinel-2 satellite imagery from Digital Earth Sweden (RISE) for the "
        "requested bounds, returned as an RGBA Raster in EPSG:3006 and "
        "defaulting to the most recent acquisition with little cloud cover."
    )
    ArgsModel = SatelliteImageryArgs
    data_category = "raw"
    result_kind = "raster"
    python_return_type = "dtcc_core.model.Raster"
    provider = [
        provider_entry("rise", role="source_provider"),
        provider_entry("dtcc-platform", role="processor"),
    ]
    source = [
        {
            "name": "Digital Earth Sweden Sentinel-2 WMS and STAC catalogue",
            "role": "source_provider",
            "service": "digital-earth-sweden",
            "selected_when": 'source="DES"',
            "source_terms_status": "requires_review",
        }
    ]
    license = (
        "Requires review: Copernicus Sentinel data is free and open, but confirm "
        f"the {provider_display_name('rise')} Digital Earth Sweden service terms "
        "before redistribution."
    )
    collection_period = (
        "Single Sentinel-2 acquisition date, either the requested date or the "
        "most recent clear date within search_days of the request."
    )
    default_crs = "EPSG:3006"
    data_types = ["raster", "imagery", "satellite", "sentinel-2"]
    geographic_coverage = "Sweden, constrained by requested bounds and source coverage"
    update_frequency = "new Sentinel-2 acquisitions every few days, subject to cloud"
    processing_steps = [
        "Validate requested bounds and imagery options",
        (
            "Without a date, query the STAC catalogue for acquisitions below "
            "max_cloud percent within search_days and probe the newest ones "
            "for real pixels"
        ),
        "Request the chosen date from the WMS server as 512-pixel tiles in EPSG:3006",
        "Stitch the tiles into one RGBA mosaic and cache it as a GeoTIFF",
        "Warn when the mosaic is empty or only partly covered",
        "Return a Raster or serialize it as tif or png",
    ]
    derived_from = [
        {
            "name": "Copernicus Sentinel-2 Level-2A surface reflectance",
            "relationship": "rendered by the Digital Earth Sweden imagery server",
            "source_terms_status": "requires_review",
        }
    ]
    presentation_headline = "Satellite Imagery"
    presentation_summary = (
        "A cloud-aware Sentinel-2 image of the requested area, as true colour "
        "or as a false colour or index rendering."
    )
    presentation_narrative = [
        {
            "heading": "What you are seeing",
            "body": (
                "Each pixel is a rendered Sentinel-2 value for one acquisition "
                "date. The fourth channel is transparency, marking pixels that "
                "the source had no data for."
            ),
        },
        {
            "heading": "How to use it",
            "body": (
                "Use rgb for a background map, and indices such as ndvi or ndwi "
                "to highlight vegetation or water. Pass date to repeat an "
                "earlier image exactly."
            ),
        },
        {
            "heading": "Limitations",
            "body": (
                "At 10 m per pixel single buildings are barely visible. Cloud "
                "cover is filtered per date, so thin or local cloud can still "
                "appear."
            ),
        },
    ]
    key_points = [
        "source='DES' uses the anonymous Digital Earth Sweden services",
        "Without date, the newest acquisition below max_cloud percent is used",
        "Imagery is returned in EPSG:3006 at 10 m unless resolution is changed",
        "Downloads are cached, so repeating a request does not download again",
    ]
    presentation_legend = {
        "title": "Rendering styles",
        "entries": [
            {"label": "rgb", "meaning": "true colour"},
            {"label": "false_color", "meaning": "near infrared, vegetation shows red"},
            {"label": "ndvi", "meaning": "vegetation index"},
            {"label": "ndwi", "meaning": "water index"},
        ],
    }
    view_hints = {
        "preferred_geometry": "raster",
        "default_crs": "EPSG:3006",
        "table_role": "basemap",
    }
    presentation_warnings = [
        "Digital Earth Sweden service terms require review before redistributing imagery.",
        "The chosen acquisition date is logged but not yet stored on the returned Raster.",
    ]
    presentation_limitations = [
        "Only Swedish coverage is available from Digital Earth Sweden.",
        "Cloud cover is judged per date for the whole area, not per pixel.",
        "Transparent pixels mark areas the chosen acquisition did not cover.",
    ]

    def build(self, args: SatelliteImageryArgs):
        bounds = self.parse_bounds(args.bounds)
        raster: Raster = dtcc_core.io.data.download_imagery(
            bounds=bounds,
            provider=args.source,
            epsg="3006",
            date=args.date,
            style=args.style,
            max_cloud=args.max_cloud,
            resolution=args.resolution,
            search_days=args.search_days,
        )
        if args.format is not None:
            return self.export_to_bytes(raster, args.format)
        return raster


__all__ = ["SatelliteImageryArgs", "SatelliteImageryDataset"]
