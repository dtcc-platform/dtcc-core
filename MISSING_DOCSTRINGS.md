Found 153 functions/methods missing docstrings:

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/builder/building/modify.py
  - get_footprint (line 22)
  - merge_building_footprints (line 51)
  - merge_building_attributes (line 111)
  - simplify_building_footprints (line 120)
  - fix_building_footprint_clearance (line 143)
  - split_footprint_walls (line 166)
  - clean_building_geometry (line 196)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/builder/city/modify.py
  - clean_building_surfaces (line 175)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/builder/geometry/multisurface.py
  - find_edge_connections (line 79)
  - group_coplanar_surfaces (line 99)
  - dfs (line 107)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/builder/geometry_builders/buildings.py
  - compute_building_heights (line 58)
  - extract_roof_points (line 141)
  - building_heights_from_pointcloud (line 185)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/builder/geometry_builders/terrain.py
  - build_terrain_mesh (line 24)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/builder/meshing/convert.py
  - mesh_to_raster (line 7)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/builder/meshing/meshing.py
  - mesh_multisurfaces (line 80)
  - merge_meshes (line 107)
  - merge (line 114)
  - snap_vertices (line 122)
  - disjoint_meshes (line 129)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/builder/pointcloud/filter.py
  - find_statistical_outliers (line 46)
  - get_vegetation (line 162)
  - pts_in_bounds (line 218)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/builder/polygons/polygons.py
  - merge_polygons_convexhull (line 29)
  - merge_polygon_hulls (line 34)
  - merge_polygons_buffering (line 45)
  - merge_polygons_snapping (line 57)
  - merge_polygons (line 83)
  - simplify_polygon (line 99)
  - remove_slivers (line 108)
  - remove_holes (line 121)
  - merge_multipolygon (line 125)
  - merge_list_of_polygons (line 213)
  - buffer_intersect_bounds (line 266)
  - widen_gaps (line 275)
  - lengthen_edges (line 303)
  - flatten_sharp_angles (line 355)
  - remove_short_edges (line 386)
  - fix_clearance (line 408)
  - split_polygon_sides (line 486)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/builder/raster/analyse.py
  - slope_aspect (line 9)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/builder/raster/filter.py
  - erode_small_lines (line 55)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/builder/register.py
  - register_model_method (line 5)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/builder/roadnetwork/convert.py
  - to_matrix (line 12)
  - to_surfaces (line 30)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/builder/surface/modify.py
  - clean_surface (line 6)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/common/dtcc_logging.py
  - error (line 45)
  - critical (line 49)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/io/city.py
  - load (line 68)
  - save (line 72)
  - buildings_to_df (line 76)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/io/cityjson/cityjson.py
  - setup_city (line 12)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/io/cityjson/utils.py
  - tin_geom_to_mesh (line 11)
  - get_terrain_mesh (line 31)
  - get_building_geometry (line 100)
  - build_dtcc_building (line 118)
  - get_root_buildings (line 143)
  - get_root_objects (line 152)
  - set_buildings (line 160)
  - set_cityobject (line 176)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/io/cityjson/write_cityjson.py
  - to_cityjson_multisurface (line 9)
  - to_cityjson (line 27)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/io/container/pointcloud_container.py
  - bounds (line 56)
  - __str__ (line 85)
  - __len__ (line 88)
  - to_proto (line 91)
  - from_proto (line 94)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/io/data/cache.py
  - empty_cache (line 8)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/io/data/lidar.py
  - add_bbox_coords (line 66)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/io/data/overpass.py
  - load_cache_metadata (line 53)
  - save_cache_metadata (line 59)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/io/data/wrapper.py
  - download_pointcloud (line 148)
  - download_footprints (line 154)
  - download_roadnetwork (line 162)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/io/footprints.py
  - list_io (line 333)
  - print_io (line 337)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/io/generic.py
  - save (line 8)
  - load (line 20)
  - list_io (line 44)
  - print_io (line 51)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/io/landuse.py
  - load (line 71)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/io/meshes.py
  - has_assimp (line 25)
  - load_mesh (line 272)
  - load_volume_mesh (line 276)
  - load_mesh_as_city (line 279)
  - save (line 283)
  - list_io (line 287)
  - print_io (line 291)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/io/pointcloud.py
  - bounds_filter_poinst (line 63)
  - load_list (line 130)
  - save (line 234)
  - list_io (line 268)
  - print_io (line 272)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/io/roadnetwork.py
  - load (line 89)
  - to_dataframe (line 107)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/io/utils.py
  - get_epsg (line 9)
  - protobuf_to_json (line 26)
  - save_to_pb (line 31)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/model/geometry/geometry.py
  - calculate_bounds (line 44)
  - bounds (line 48)
  - bounds (line 54)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/model/geometry/grid.py
  - __str__ (line 28)
  - calculate_bounds (line 33)
  - __str__ (line 159)
  - calculate_bounds (line 162)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/model/geometry/linestring.py
  - from_proto (line 78)
  - from_proto (line 144)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/model/geometry/pointcloud.py
  - __repr__ (line 54)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/model/geometry/polygon.py
  - to_proto (line 12)
  - from_proto (line 15)
  - shapely (line 19)
  - vertices (line 23)
  - holes (line 27)
  - area (line 31)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/model/geometry/surface.py
  - xmin (line 43)
  - ymin (line 47)
  - zmin (line 51)
  - xmax (line 55)
  - ymax (line 59)
  - zmax (line 63)
  - centroid (line 67)
  - __str__ (line 211)
  - __len__ (line 237)
  - zmax (line 262)
  - copy (line 286)
  - __str__ (line 334)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/model/mixins/city/io_mixins.py
  - load_terrain_raster (line 85)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/model/model.py
  - to_proto (line 17)
  - from_proto (line 21)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/model/object/building.py
  - height (line 23)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/model/object/city.py
  - has_terrain (line 38)
  - get_building_attributes (line 59)
  - remove_terrain (line 88)
  - remove_buildings (line 92)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/model/object/landuse.py
  - surfaces (line 33)
  - to_proto (line 39)
  - from_proto (line 42)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/model/object/object.py
  - from_str (line 58)
  - bounds (line 176)
  - get_children (line 243)
  - set_child_attributues (line 246)
  - get_child_attributes (line 256)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/model/object/roadnetwork.py
  - linestrings (line 38)
  - multilinestrings (line 45)
  - bounds (line 50)
  - to_shapely (line 60)
  - to_proto (line 66)
  - from_proto (line 78)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/model/object/terrain.py
  - __str__ (line 52)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/model/object/tree.py
  - to_proto (line 18)
  - from_proto (line 21)

File: /Users/dwastberg/repos/dtcc2/dtcc-core/dtcc_core/model/values/raster.py
  - copy (line 248)

