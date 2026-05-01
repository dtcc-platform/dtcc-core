import dtcc_core
from dtcc_core.model import RoadNetwork
from pydantic import Field
from typing import Literal, Optional

from .dataset import DatasetDescriptor, DatasetBaseArgs


class RoadsArgs(DatasetBaseArgs):
    source: Literal["OSM"] = Field("OSM", description="Data source")
    format: Optional[Literal["pb"]] = Field(
        None, description="Output format (pb for protobuf bytes)"
    )


class RoadsDataset(DatasetDescriptor):
    name = "roads"
    description = "Road network data from OSM/Overpass."
    ArgsModel = RoadsArgs

    def build(self, args: RoadsArgs):
        bounds = self.parse_bounds(args.bounds)
        roads: RoadNetwork = dtcc_core.io.data.download_roadnetwork(
            bounds=bounds,
            provider=args.source,
            epsg="3006",
        )
        if args.format == "pb":
            return roads.to_proto().SerializeToString()
        return roads


__all__ = ["RoadsArgs", "RoadsDataset"]
