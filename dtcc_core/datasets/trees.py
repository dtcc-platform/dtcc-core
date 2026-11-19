import warnings
from typing import Any, Literal

from pydantic import Field, model_validator

import dtcc_core
from dtcc_core.builder import tree_raster_from_pointcloud
from dtcc_core.datasets import DatasetBaseArgs, DatasetDescriptor
from dtcc_core.datasets.providers import provider_entry
from dtcc_core.io.trees import save_trees
from dtcc_core.model import City, PointCloud, TreeCollection


class TreeArgs(DatasetBaseArgs):
    """Arguments for the trees dataset derived from point cloud vegetation."""

    tree_type: Literal["urban", "mixed", "dense", "arid"] = Field(
        "urban", description="Type of trees to detect"
    )
    cell_size: float = Field(
        0.5, description="Cell size of the output tree raster (in meters)"
    )

    vector_geometry: Literal["point", "circle"] = Field(
        "point",
        description="Save trees as points or circles when exporting to vector format",
    )

    format: Literal["tif", "gpkg", "geojson"] | None = Field(
        None, description="Output file format"
    )

    product: Literal["trees", "raster"] = Field(
        "trees",
        description="Return detected trees or the tree-height raster",
    )

    @model_validator(mode="after")
    def validate_product(self):
        if self.product == "raster" and "vector_geometry" in self.model_fields_set:
            raise ValueError("vector_geometry applies only to product='trees'")
        if self.format == "tif" and self.product == "trees":
            if "product" in self.model_fields_set:
                raise ValueError("format='tif' requires product='raster'")
            warnings.warn(
                "format='tif' without product is deprecated; pass product='raster'",
                DeprecationWarning,
                stacklevel=2,
            )
            self.product = "raster"
        if self.product == "raster" and self.format not in (None, "tif"):
            raise ValueError("product='raster' supports only format='tif'")
        return self


class TreesDataset(DatasetDescriptor):
    """Trees dataset, registered as ``trees``.

    Tree locations or tree-height raster derived from vegetation returns in
    Lantmäteriet point cloud data for the requested bounds.

    Call it with keyword arguments defined by ``TreeArgs``.
    """

    name = "trees"
    title = "Trees"
    description = (
        "Tree locations or tree-height raster derived from vegetation returns in "
        "Lantmäteriet point cloud data for the requested bounds."
    )
    ArgsModel = TreeArgs
    data_category = "derived"
    result_kind = "tree_collection"
    python_return_type = "dtcc_core.model.TreeCollection"
    export_capabilities = [
        {
            "product": "trees",
            "canonical": True,
            "web_formats": ["geojson"],
            "constraints": {},
            "options": {},
            "exclusions": [],
        },
        {
            "product": "raster",
            "canonical": True,
            "web_formats": ["png", "tif"],
            "constraints": {},
            "options": {
                "colormap": {"formats": ["png"], "default": "viridis"},
                "value_range": {"formats": ["png"], "default": "valid min/max"},
                "unit": {"formats": ["png", "tif"], "default": "m"},
            },
            "exclusions": [],
        },
    ]
    provider = [
        provider_entry("lantmateriet", role="source_provider"),
        provider_entry("dtcc-platform", role="processor"),
    ]
    source = [
        {
            "name": "DTCC/Lantmäteriet point cloud backend/cache",
            "role": "upstream_dataset",
            "source_terms_status": "requires_review",
        },
        {
            "name": "dtcc-core tree extraction pipeline",
            "role": "processor",
        },
    ]
    license = (
        "Requires review: derived from upstream point cloud data, so "
        "Lantmäteriet/source redistribution terms must be verified."
    )
    collection_period = (
        "Requires review: inherits point cloud acquisition date/cache vintage, "
        "which is not currently surfaced in the result."
    )
    default_crs = "EPSG:3006"
    data_types = ["point_cloud_derived", "tree_points", "tree_height_raster"]
    geographic_coverage = "Sweden, constrained by requested bounds and source coverage"
    update_frequency = "derived on demand from upstream source data"
    processing_steps = [
        "Download point cloud data for requested bounds in EPSG:3006",
        "Extract vegetation-classified point returns from the point cloud",
        "Build or use a terrain raster to estimate above-ground canopy height",
        "Apply the selected tree_type threshold profile",
        "Return detected tree points or export a tree-height raster/vector format",
    ]
    derived_from = [
        {
            "name": "point_cloud",
            "relationship": "vegetation and terrain returns drive tree extraction",
            "source_terms_status": "requires_review",
        }
    ]
    presentation_headline = "Detected Tree Layer"
    presentation_summary = (
        "A derived vegetation layer that estimates tree positions or canopy-height "
        "raster cells from source point cloud classifications."
    )
    presentation_narrative = [
        {
            "heading": "What you are seeing",
            "body": (
                "The vector result is a collection of detected trees with position, "
                "height, and crown-radius estimates. The raster result is a "
                "canopy-height surface built from vegetation returns."
            ),
        },
        {
            "heading": "How it is made",
            "body": (
                "The pipeline separates vegetation returns, estimates terrain, "
                "builds a canopy-height model, and applies the selected tree_type "
                "threshold profile."
            ),
        },
        {
            "heading": "Limitations",
            "body": (
                "This is an automated extraction, not a field survey. Individual "
                "trees can be merged, missed, or shifted when source classifications "
                "or canopy segmentation are imperfect."
            ),
        },
    ]
    key_points = [
        "Derived from the same upstream point cloud source used by point_cloud",
        "tree_type selects detection thresholds for urban, mixed, dense, or arid settings",
        "product='raster' returns a tree-height raster; the default returns detected trees",
        "vector exports can use point or circle geometry",
    ]
    presentation_legend = {
        "title": "Tree extraction",
        "entries": [
            {"label": "Tree point", "meaning": "estimated tree top/base position"},
            {
                "label": "Circle",
                "meaning": "estimated crown radius when vector_geometry='circle'",
            },
            {
                "label": "Raster value",
                "meaning": "estimated canopy height above terrain",
            },
        ],
    }
    view_hints = {
        "preferred_geometry": "points",
        "alternate_geometry": "raster",
        "default_crs": "EPSG:3006",
        "tree_type_profiles": ["urban", "mixed", "dense", "arid"],
        "value_units": "m",
    }
    presentation_warnings = [
        "Tree extraction quality depends on upstream point cloud density, classification quality, and acquisition season.",
        "Derived tree objects should not be treated as surveyed individual-tree inventory.",
        "Source and license terms inherit the upstream point cloud review status.",
    ]
    presentation_limitations = [
        "The result does not currently expose source point density, scan date, or confidence per tree.",
        "Nearby crowns can merge and sparse canopies can be missed.",
        "Building and terrain interactions can affect canopy-height estimates.",
    ]

    def describe(self) -> dict[str, Any]:
        metadata = super().describe()
        metadata["products"] = [
            {
                "name": "trees",
                "python_return_type": "dtcc_core.model.TreeCollection",
                "formats": ["gpkg", "geojson"],
                "description": "Detected trees with position, height, and crown radius.",
            },
            {
                "name": "raster",
                "python_return_type": "dtcc_core.model.Raster",
                "formats": ["tif"],
                "description": "Canopy-height raster in meters built from vegetation returns.",
            },
        ]
        return metadata

    def build(self, args: TreeArgs):
        bounds = self.parse_bounds(args.bounds)
        pc: PointCloud = dtcc_core.io.data.download_pointcloud(bounds=bounds)
        if args.product == "raster":
            raster = tree_raster_from_pointcloud(
                pc,
                None,
                tree_type=args.tree_type,
                cell_size=args.cell_size,
            )
            raster.crs = self.default_crs
            if args.format == "tif":
                return self.export_to_bytes(raster, "tif")
            return raster
        else:
            city = City()
            city.add_point_cloud(pc)
            city.bounds = bounds
            trees = city.build_trees_from_pointcloud(tree_type=args.tree_type)
            if args.format is None:
                return trees
            elif args.format in ("gpkg", "geojson", "json"):
                if args.format in ("geojson", "json"):
                    args.format = "json"
                as_circles = args.vector_geometry == "circle"
                return self.export_to_bytes(
                    trees,
                    args.format,
                    save_callable=save_trees,
                    as_circles=as_circles,
                )

    def prepare_result(self, result, validated_args: TreeArgs):
        if validated_args.format is None and isinstance(result, list):
            return TreeCollection(result)
        return result
