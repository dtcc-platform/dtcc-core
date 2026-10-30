import datetime
from typing import Literal

import requests
from pydantic import Field

import dtcc_core
from dtcc_core.common.progress import ProgressTracker
from dtcc_core.io.data.digitalearth import (
    DEFAULT_MAX_CLOUD,
    DEFAULT_RESOLUTION,
    DEFAULT_SEARCH_DAYS,
    STAC_URL,
    VALID_STYLES,
    WMS_URL,
    ServiceError,
)
from dtcc_core.model import Raster

from .dataset import DatasetBaseArgs, DatasetDescriptor, DatasetUpstreamError
from .providers import provider_display_name, provider_entry

ImageryStyle = Literal[tuple(VALID_STYLES)]


class SatelliteImageryArgs(DatasetBaseArgs):
    """Arguments for the Digital Earth Sweden satellite imagery dataset."""

    source: Literal["DES"] = Field(
        "DES",
        description=(
            "Imagery source. DES uses the anonymous Digital Earth Sweden "
            "Sentinel-2 services operated by RISE."
        ),
    )
    date: datetime.date | None = Field(
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
        DEFAULT_MAX_CLOUD,
        ge=0.0,
        le=100.0,
        description="Highest cloud cover, in percent, accepted when picking a date.",
    )
    resolution: float = Field(
        DEFAULT_RESOLUTION,
        gt=0.0,
        description="Ground sample distance in metres; Sentinel-2 is natively 10 m.",
    )
    search_days: int = Field(
        DEFAULT_SEARCH_DAYS,
        ge=1,
        description="How many days back to search for a clear acquisition.",
    )
    format: Literal["tif", "png"] | None = Field(None, description="Output file format")


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
    # The largest allowed request is about 400 tiles at about 2 s each.
    timeout_hint = 900
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
        "With an explicit date, return the cached GeoTIFF if there is one",
        "Probe the date in rgb and fail if the area has no imagery",
        "Request the chosen date from the WMS server as 512-pixel tiles in EPSG:3006",
        (
            "Retry failed requests; without strict_live, leave a tile that keeps "
            "failing transparent, otherwise raise a typed upstream error"
        ),
        "Stitch the tiles into one RGBA mosaic",
        "Fail if the mosaic is empty; warn and skip the cache if partly covered",
        "Cache complete mosaics as GeoTIFFs",
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
                "the source had no data for. Index styles also make pixels "
                "transparent where the index does not apply, such as land in "
                "ndwi."
            ),
        },
        {
            "heading": "How to use it",
            "body": (
                "Use rgb for a background map, and the index styles such as "
                "ndvi or ndwi to see where vegetation or water is. Index styles "
                "are colour ramps, not index values. Pass date to repeat an "
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
        (
            "Complete downloads are cached, so repeating a request with the same "
            "date does not download again; without date, the catalogue is still "
            "queried each time"
        ),
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
        (
            "Rendered 8-bit visualization, not surface reflectance or index "
            "values; do not use for quantitative analysis."
        ),
        "Digital Earth Sweden service terms require review before redistributing imagery.",
        "The chosen acquisition date is logged but not yet stored on the returned Raster.",
    ]
    presentation_limitations = [
        "Only Swedish coverage is available from Digital Earth Sweden.",
        "Cloud cover is judged per date for the whole area, not per pixel.",
        (
            "One request is limited to 100 million pixels, a 100 km square at "
            "10 m, and the area must be at least one pixel wide and tall."
        ),
        (
            "Transparent pixels mark areas the chosen acquisition did not cover, "
            "and in index styles also pixels the index masks out."
        ),
    ]

    def build(self, args: SatelliteImageryArgs):
        progress_phases = {"download": 0.95, "export": 0.05}
        with ProgressTracker(total=1.0, phases=progress_phases) as progress:
            bounds = self.parse_bounds(args.bounds)
            # Without strict_live a tile that keeps failing is left transparent
            # and the partial result is returned uncached with a warning.
            tile_errors = None if args.strict_live else []
            with progress.phase("download", "Downloading satellite imagery..."):
                try:
                    raster: Raster = dtcc_core.io.data.download_imagery(
                        bounds=bounds,
                        provider=args.source,
                        epsg="3006",
                        date=args.date,
                        style=args.style,
                        max_cloud=args.max_cloud,
                        resolution=args.resolution,
                        search_days=args.search_days,
                        tile_errors=tile_errors,
                    )
                except requests.RequestException as exc:
                    raise self._upstream_error(exc) from exc

            if args.format is None:
                return raster
            with progress.phase("export", f"Exporting imagery as {args.format}..."):
                return self.export_to_bytes(raster, args.format)

    def _upstream_error(self, exc: requests.RequestException) -> DatasetUpstreamError:
        """Turn a failed Digital Earth Sweden request into a typed error."""
        operation, target = _failed_operation(exc)
        if isinstance(exc, ServiceError):
            return DatasetUpstreamError(
                dataset=self.name,
                operation=operation,
                target=target,
                failure_class="invalid_payload",
                message=f"{self.name} {operation} failed for {target}: {exc}",
            )
        return self.build_upstream_error(self.name, operation, target, exc)


def _failed_operation(exc: requests.RequestException) -> tuple[str, str]:
    """Name the Digital Earth Sweden service a request exception came from."""
    request = getattr(exc, "request", None)
    url = getattr(request, "url", None) or ""
    if url.startswith(STAC_URL):
        return "stac_search", STAC_URL
    return "wms_getmap", WMS_URL


__all__ = ["SatelliteImageryArgs", "SatelliteImageryDataset"]
