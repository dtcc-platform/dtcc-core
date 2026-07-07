import dtcc_core
from dtcc_core.builder import tree_raster_from_pointcloud
from dtcc_core.model import PointCloud, City, TreeCollection
from dtcc_core.io.trees import save_trees
from pydantic import Field
from typing import Optional, Literal

from dtcc_core.datasets import DatasetDescriptor, DatasetBaseArgs
from dtcc_core.datasets.providers import provider_entry


class TreeArgs(DatasetBaseArgs):
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

    format: Optional[Literal["tif", "gpkg", "geojson"]] = Field(
        None, description="Output file format"
    )


class TreesDataset(DatasetDescriptor):
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
        "format='tif' exports a tree-height raster",
        "vector exports can use point or circle geometry",
    ]
    presentation_legend = {
        "title": "Tree extraction",
        "entries": [
            {"label": "Tree point", "meaning": "estimated tree top/base position"},
            {"label": "Circle", "meaning": "estimated crown radius when vector_geometry='circle'"},
            {"label": "Raster value", "meaning": "estimated canopy height above terrain"},
        ],
    }
    view_hints = {
        "preferred_geometry": "points",
        "alternate_geometry": "raster",
        "default_crs": "EPSG:3006",
        "tree_type_profiles": ["urban", "mixed", "dense", "arid"],
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

    def build(self, args: TreeArgs):
        bounds = self.parse_bounds(args.bounds)
        pc: PointCloud = dtcc_core.io.data.download_pointcloud(bounds=bounds)
        if args.format == "tif":
            raster = tree_raster_from_pointcloud(
                pc,
                None,
                tree_type=args.tree_type,
                cell_size=args.cell_size,
            )
            return self.export_to_bytes(raster, "tif")
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
