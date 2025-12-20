import dtcc_core
from dtcc_core.model import City, Bounds
from typing import Literal, Optional, List, Tuple, Sequence, Union
from pydantic import BaseModel, Field
import tempfile

from .dataset import DatasetDescriptor, DatasetBaseArgs


class BuildingArgs(DatasetBaseArgs):
    source: Literal["OSM", "LM"] = Field(
        "LM", description="Data source for building footprints"
    )
    lod: Literal[1] = Field(1, description="Level of Detail for building")
    smallest_building_size: float = Field(
        15.0, description="Smallest building size to include (in square meters)"
    )
    place_on_zero: bool = Field(
        False, description="Whether to place buildings on Z=0 plane"
    )
    format: Optional[Literal["obj", "stl", "cityjson", "json"]] = Field(
        None, description="Output file format"
    )


class BuildingDataset(DatasetDescriptor):
    name = "buildings"
    description = "Generate 3D buildings from point cloud and building footprints."
    ArgsModel = BuildingArgs

    def build(self, args: BuildingArgs):
        bounds = self.parse_bounds(args.bounds)
        city = City()
        city.bounds = bounds
        city.download_pointcloud()
        city.download_footprints()
        city.add_point_cloud(city.pointcloud.remove_global_outliers(3.0))

        city.build_lod1_buildings()
        if args.format is None:
            return city.buildings
        elif args.format in ("cityjson", "json"):
            self.export_to_bytes(city, "json", as_text=True)
        elif args.format in ("obj", "stl"):
            building_meshes = [
                b.lod1.mesh(weld=True, snap=0.005) for b in city.buildings
            ]
            if args.place_on_zero:
                for b in building_meshes:
                    b.offset([0, 0, -b.bounds.zmin])
            merged_mesh = dtcc_core.builder.meshing.merge_meshes(building_meshes)
            return self.export_to_bytes(merged_mesh, args.format)
