# Missing Docstrings

Found 92 functions/methods missing docstrings.

Run `python scripts/find_missing_docstrings.py` to regenerate this file.

## dtcc_core/builder/geometry/multisurface.py

- `dfs` (line 170)

## dtcc_core/builder/meshing/mesh_tiler.py

- `__init__` (line 166)

## dtcc_core/builder/meshing/tetgen.py

- `is_tetgen_available` (line 21)

## dtcc_core/common/dtcc_logging.py

- `error` (line 45)
- `critical` (line 49)

## dtcc_core/io/cityjson/converters.py

- `__post_init__` (line 24)
- `__init__` (line 114)
- `add_points` (line 132)
- `converter` (line 311)
- `converter` (line 320)

## dtcc_core/io/cityjson/utils.py

- `get_root_buildings` (line 223)
- `get_root_objects` (line 232)
- `set_buildings` (line 240)

## dtcc_core/io/container/pointcloud_container.py

- `bounds` (line 56)
- `__str__` (line 85)
- `__len__` (line 88)
- `to_proto` (line 91)
- `from_proto` (line 94)

## dtcc_core/io/data/lidar.py

- `add_bbox_coords` (line 66)

## dtcc_core/io/data/overpass.py

- `load_cache_metadata` (line 53)
- `save_cache_metadata` (line 59)

## dtcc_core/io/footprints.py

- `list_io` (line 332)
- `print_io` (line 336)

## dtcc_core/io/generic.py

- `save` (line 8)
- `load` (line 32)
- `list_io` (line 56)
- `print_io` (line 63)

## dtcc_core/io/landuse.py

- `load` (line 71)

## dtcc_core/io/meshes.py

- `has_assimp` (line 31)
- `load_mesh` (line 327)
- `load_volume_mesh` (line 331)
- `load_mesh_as_city` (line 335)

## dtcc_core/io/pointcloud.py

- `bounds_filter_poinst` (line 63)
- `load_list` (line 130)
- `save` (line 234)
- `list_io` (line 268)
- `print_io` (line 272)

## dtcc_core/io/roadnetwork.py

- `load` (line 89)
- `to_dataframe` (line 107)

## dtcc_core/io/utils.py

- `get_epsg` (line 9)
- `protobuf_to_json` (line 26)
- `save_to_pb` (line 31)

## dtcc_core/model/geometry/geometry.py

- `calculate_bounds` (line 44)
- `bounds` (line 48)
- `bounds` (line 54)

## dtcc_core/model/geometry/grid.py

- `__str__` (line 28)
- `calculate_bounds` (line 33)
- `__str__` (line 159)
- `calculate_bounds` (line 162)

## dtcc_core/model/geometry/linestring.py

- `from_proto` (line 89)
- `from_proto` (line 163)

## dtcc_core/model/geometry/pointcloud.py

- `__repr__` (line 54)

## dtcc_core/model/geometry/polygon.py

- `to_proto` (line 12)
- `from_proto` (line 15)
- `shapely` (line 19)
- `vertices` (line 23)
- `holes` (line 27)
- `area` (line 31)

## dtcc_core/model/geometry/surface.py

- `xmin` (line 43)
- `ymin` (line 47)
- `zmin` (line 51)
- `xmax` (line 55)
- `ymax` (line 59)
- `zmax` (line 63)
- `centroid` (line 67)
- `__str__` (line 215)
- `__len__` (line 241)
- `zmax` (line 266)
- `copy` (line 290)
- `__str__` (line 338)

## dtcc_core/model/mixins/city/io_mixins.py

- `load_terrain_raster` (line 85)

## dtcc_core/model/model.py

- `to_proto` (line 17)
- `from_proto` (line 21)

## dtcc_core/model/object/city.py

- `has_terrain` (line 51)
- `get_building_attributes` (line 72)
- `remove_terrain` (line 101)
- `remove_buildings` (line 105)

## dtcc_core/model/object/landuse.py

- `surfaces` (line 62)
- `to_proto` (line 68)
- `from_proto` (line 71)

## dtcc_core/model/object/object.py

- `from_class_name` (line 92)
- `from_class` (line 96)

## dtcc_core/model/object/roadnetwork.py

- `linestrings` (line 54)
- `multilinestrings` (line 61)
- `bounds` (line 66)
- `to_shapely` (line 76)
- `to_proto` (line 82)
- `from_proto` (line 94)

## dtcc_core/model/object/terrain.py

- `__str__` (line 52)

## dtcc_core/model/object/tree.py

- `to_proto` (line 31)
- `from_proto` (line 34)

## dtcc_core/model/values/raster.py

- `copy` (line 248)

