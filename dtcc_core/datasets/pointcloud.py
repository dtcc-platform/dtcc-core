import dtcc_core
from dtcc_core.model import PointCloud
from typing import Literal, Optional, List, Union
from pydantic import Field, field_validator

from .dataset import DatasetDescriptor, DatasetBaseArgs
from .providers import provider_display_name, provider_entry
from dtcc_core.common.progress import ProgressTracker, report_progress


CLASSIFICATION_PRESETS = {
    "terrain": [2, 8],
    "buildings": [6, 9],
    "vegetation": [3, 4, 5, 7],
}


class PointCloudArgs(DatasetBaseArgs):
    classifications: Union[
        int, list[int], Literal["all", "terrain", "buildings", "vegetation"]
    ] = Field(
        "all",
        description="Which classifications to include (e.g., [2, 9] for ground and water, or 'vegetation' for all vegetation classes). 'all' includes all points.",
    )
    source: Literal["LM"] = Field(
        "LM",
        description=(
            "Data source. LM uses the DTCC/Lantmäteriet point cloud backend/cache."
        ),
    )
    format: Optional[Literal["copc", "las", "laz"]] = Field(
        None, description="Output file format"
    )
    remove_outliers: bool = Field(
        False, description="Whether to remove global outliers from the point cloud"
    )
    remove_outlier_threshold: float = Field(
        3.0, description="Threshold for outlier removal"
    )

    crs: Optional[str] = Field(
        None,
        description=(
            "Coordinate reference system for context/export. Point cloud "
            "coordinates are currently available only in EPSG:3006."
        ),
    )

    @field_validator("classifications")
    @classmethod
    def validate_classifications(cls, value):
        if isinstance(value, list) and not value:
            raise ValueError(
                "classifications list must contain at least one classification code."
            )
        return value

    @field_validator("crs")
    @classmethod
    def validate_crs(cls, value):
        if value is None:
            return value
        crs = str(value).strip().upper()
        if crs in {"3006", "EPSG:3006"}:
            return "EPSG:3006"
        raise ValueError(
            "point_cloud does not currently reproject coordinates; crs must be "
            "'EPSG:3006' or omitted."
        )


class PointCloudDataset(DatasetDescriptor):
    name = "point_cloud"
    title = "Point Cloud"
    description = (
        "Lantmäteriet point cloud tiles for the requested bounds, accessed through "
        "the DTCC backend/cache with optional classification filtering and global "
        "outlier removal."
    )
    ArgsModel = PointCloudArgs
    data_category = "raw"
    result_kind = "point_cloud"
    python_return_type = "dtcc_core.model.PointCloud"
    provider = [
        provider_entry("lantmateriet", role="source_provider"),
        provider_entry("dtcc-platform", role="processor"),
    ]
    source = [
        {
            "name": "DTCC/Lantmäteriet point cloud backend/cache",
            "role": "source_provider",
            "selected_when": 'source="LM"',
            "source_terms_status": "requires_review",
        }
    ]
    license = (
        f"Requires review: verify {provider_display_name('lantmateriet')} point "
        "cloud source terms before redistribution."
    )
    collection_period = (
        "Requires review: acquisition date, scan campaign, and cache tile vintage "
        "are not currently surfaced in the dataset result."
    )
    default_crs = "EPSG:3006"
    data_types = ["point_cloud", "lidar", "classified_points"]
    geographic_coverage = "Sweden, constrained by requested bounds and source coverage"
    update_frequency = (
        f"varies by {provider_display_name('lantmateriet')} source product"
    )
    processing_steps = [
        "Map source='LM' to the DTCC/Lantmäteriet point cloud backend/cache",
        "Download point cloud tiles intersecting the requested bounds in EPSG:3006",
        "Resolve classification presets: terrain=[2, 8], buildings=[6, 9], vegetation=[3, 4, 5, 7]",
        "Keep only requested classifications unless classifications='all'",
        "Optionally remove global Z outliers using remove_outlier_threshold",
        "Return a PointCloud or serialize the selected point-cloud format",
    ]
    derived_from = [
        {
            "name": f"{provider_display_name('lantmateriet')} point cloud source product",
            "relationship": "source data selected through the DTCC backend/cache",
            "source_terms_status": "requires_review",
        }
    ]
    presentation_headline = "Classified Point Cloud"
    presentation_summary = (
        "A bounded lidar point cloud that can be filtered to terrain, building, "
        "or vegetation classification groups before downstream city-model work."
    )
    presentation_narrative = [
        {
            "heading": "What you are seeing",
            "body": (
                "Each point is an upstream lidar return in the requested bounds. "
                "Classification presets are convenience groups over the source "
                "classification codes."
            ),
        },
        {
            "heading": "How to use it",
            "body": (
                "Use the full point cloud for inspection, terrain/building/"
                "vegetation presets for targeted extraction, and outlier removal "
                "only when extreme elevations are known to be noise."
            ),
        },
        {
            "heading": "Limitations",
            "body": (
                "Coverage, scan date, point density, classification quality, and "
                "redistribution terms depend on the selected upstream tile."
            ),
        },
    ]
    key_points = [
        "source='LM' uses the DTCC/Lantmäteriet point cloud backend/cache",
        "classifications='terrain' keeps classes [2, 8]",
        "classifications='buildings' keeps classes [6, 9]",
        "classifications='vegetation' uses the point-cloud vegetation extractor",
        "remove_outliers=True removes global Z outliers after classification filtering",
    ]
    presentation_legend = {
        "title": "Classification presets",
        "entries": [
            {"label": "terrain", "meaning": "classes 2 and 8"},
            {"label": "buildings", "meaning": "classes 6 and 9"},
            {"label": "vegetation", "meaning": "classes 3, 4, 5, and 7"},
            {"label": "all", "meaning": "no classification filter"},
        ],
    }
    view_hints = {
        "preferred_geometry": "points",
        "default_crs": "EPSG:3006",
        "recommended_filters": list(CLASSIFICATION_PRESETS),
    }
    presentation_warnings = [
        "Source and license terms require review before redistributing point cloud packages.",
        "Classification codes are treated as source-provided labels; the dataset does not validate classification accuracy.",
        "Global outlier removal can discard real high or low points when the requested bounds contain steep terrain or tall structures.",
    ]
    presentation_limitations = [
        "Acquisition date, density, and quality metadata are not currently surfaced per tile.",
        "Point cloud coverage can be incomplete or unavailable for some requested bounds.",
        "The crs argument is restricted to EPSG:3006 because point-cloud serializers do not currently reproject coordinates.",
    ]

    @staticmethod
    def _resolve_classifications(classifications) -> List[int]:
        """Resolve classification parameter to a list of class IDs."""
        if isinstance(classifications, int):
            return [classifications]
        elif isinstance(classifications, list):
            if not classifications:
                raise ValueError(
                    "classifications list must contain at least one classification code."
                )
            return classifications
        elif isinstance(classifications, str):
            try:
                return CLASSIFICATION_PRESETS[classifications]
            except KeyError as exc:
                raise ValueError(
                    f"Unsupported point cloud classification preset: {classifications!r}"
                ) from exc
        raise ValueError(
            "classifications must be an integer, non-empty integer list, or supported preset."
        )

    def build(self, args: PointCloudArgs):
        progress_phases = {
            "download_pointcloud": 0.50,
            "filter_classifications": 0.15,
            "remove_outliers": 0.15,
            "export": 0.20,
        }
        with ProgressTracker(total=1.0, phases=progress_phases) as progress:
            bounds = self.parse_bounds(args.bounds)

            with progress.phase(
                "download_pointcloud", "Downloading point cloud data..."
            ):
                pc: PointCloud = dtcc_core.io.data.download_pointcloud(
                    bounds=bounds,
                    provider="dtcc",
                )

            with progress.phase(
                "filter_classifications",
                (
                    "Filtering point cloud classifications..."
                    if args.classifications not in ("all", None)
                    else "Using all classifications"
                ),
            ):
                if args.classifications == "vegetation":
                    report_progress(
                        percent=30, message="Extracting vegetation points..."
                    )
                    pc = pc.get_vegetation()
                elif args.classifications not in ("all", None):
                    classifications = self._resolve_classifications(
                        args.classifications
                    )
                    if classifications:
                        report_progress(
                            percent=30,
                            message=f"Filtering to classes {classifications}...",
                        )
                        pc = pc.classification_filter(classifications, keep=True)

            with progress.phase(
                "remove_outliers",
                (
                    f"Removing outliers (threshold={args.remove_outlier_threshold})..."
                    if args.remove_outliers
                    else "Skipping outlier removal"
                ),
            ):
                if args.remove_outliers:
                    pc = pc.remove_global_outliers(args.remove_outlier_threshold)

            with progress.phase(
                "export",
                (
                    f"Exporting to {args.format}..."
                    if args.format
                    else "Preparing point cloud result..."
                ),
            ):
                if args.format is not None:
                    return self.export_to_bytes(pc, args.format)
                return pc
