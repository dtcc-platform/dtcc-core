import dtcc_core
from dtcc_core.model import City, Bounds
from typing import Literal, Optional, List, Tuple, Sequence, Union
from pydantic import BaseModel, Field
import tempfile

from .dataset import DatasetDescriptor


class BuildingArgs(BaseModel):
    bounds: Sequence[float] = Field(
        ...,
        description="Bounding box [minx, miny, [optional_minz], maxx, maxy, [optional_maxz]]",
    )

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
    format: Optional[Literal["obj", "stl", "cityjson"]] = Field(
        None, description="Output file format"
    )


class BuildingDataset(DatasetDescriptor):
    name = "buildings"
    ArgsModel = BuildingArgs

    def validate(self, kwargs) -> BuildingArgs:
        args = self.ArgsModel(**kwargs)
        if len(args.bounds) not in (4, 6):
            raise ValueError("Bounds must be a list of 4 or 6 floats.")

        return args

    def build(self, args: BuildingArgs):
        if len(args.bounds) == 4:
            bounds = Bounds(
                xmin=args.bounds[0],
                ymin=args.bounds[1],
                xmax=args.bounds[2],
                ymax=args.bounds[3],
            )
        else:
            bounds = Bounds(
                xmin=args.bounds[0],
                ymin=args.bounds[1],
                zmin=args.bounds[2],
                xmax=args.bounds[3],
                ymax=args.bounds[4],
                zmax=args.bounds[5],
            )

        city = City()
        city.bounds = bounds
        city.download_pointcloud()
        city.download_footprints()
        city.add_point_cloud(city.pointcloud.remove_global_outliers(3.0))

        city.build_lod1_buildings()
        if args.format is None:
            return city.buildings
        elif args.format == "cityjson":
            with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as tmpfile:
                city.save_cityjson(tmpfile.name)
                with open(tmpfile.name, "r") as f:
                    cityjson_data = f.read()
            return cityjson_data
        elif args.format in ("obj", "stl"):
            building_meshes = [
                b.lod1.mesh(weld=True, snap=0.005) for b in city.buildings
            ]
            if args.place_on_zero:
                for b in building_meshes:
                    b.offset([0, 0, -b.bounds.zmin])
            merged_mesh = dtcc_core.builder.meshing.merge_meshes(building_meshes)
            with tempfile.NamedTemporaryFile(
                suffix="." + args.format, delete=False
            ) as tmpfile:
                merged_mesh.save(tmpfile.name)
                with open(tmpfile.name, "rb") as f:
                    mesh_data = f.read()
            return mesh_data
