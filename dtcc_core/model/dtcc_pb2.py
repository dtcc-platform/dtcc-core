# -*- coding: utf-8 -*-
# Generated by the protocol buffer compiler.  DO NOT EDIT!
# NO CHECKED-IN PROTOBUF GENCODE
# source: dtcc.proto
# Protobuf Python Version: 5.27.3
"""Generated protocol buffer code."""
from google.protobuf import descriptor as _descriptor
from google.protobuf import descriptor_pool as _descriptor_pool
from google.protobuf import runtime_version as _runtime_version
from google.protobuf import symbol_database as _symbol_database
from google.protobuf.internal import builder as _builder
_runtime_version.ValidateProtobufRuntimeVersion(
    _runtime_version.Domain.PUBLIC,
    5,
    27,
    3,
    '',
    'dtcc.proto'
)
# @@protoc_insertion_point(imports)

_sym_db = _symbol_database.Default()




DESCRIPTOR = _descriptor_pool.Default().AddSerializedFile(b'\n\ndtcc.proto\x12\x04\x44TCC\"\xe2\x03\n\x06Object\x12\n\n\x02id\x18\x01 \x01(\t\x12\x12\n\nattributes\x18\x02 \x01(\t\x12\x1e\n\x08\x63hildren\x18\x03 \x03(\x0b\x32\x0c.DTCC.Object\x12,\n\x08geometry\x18\x04 \x03(\x0b\x32\x1a.DTCC.Object.GeometryEntry\x12\x1c\n\x06\x62ounds\x18\x05 \x01(\x0b\x32\x0c.DTCC.Bounds\x12\x1a\n\x04\x63ity\x18\x06 \x01(\x0b\x32\n.DTCC.CityH\x00\x12\"\n\x08\x62uilding\x18\x07 \x01(\x0b\x32\x0e.DTCC.BuildingH\x00\x12 \n\x07terrain\x18\x08 \x01(\x0b\x32\r.DTCC.TerrainH\x00\x12\'\n\x0b\x63ity_object\x18\t \x01(\x0b\x32\x10.DTCC.CityObjectH\x00\x12+\n\rbuilding_part\x18\n \x01(\x0b\x32\x12.DTCC.BuildingPartH\x00\x12)\n\x0croad_network\x18\x0b \x01(\x0b\x32\x11.DTCC.RoadNetworkH\x00\x12 \n\x07landuse\x18\x0c \x01(\x0b\x32\r.DTCC.LanduseH\x00\x1a?\n\rGeometryEntry\x12\x0b\n\x03key\x18\x01 \x01(\t\x12\x1d\n\x05value\x18\x02 \x01(\x0b\x32\x0e.DTCC.Geometry:\x02\x38\x01\x42\x06\n\x04type\"\x06\n\x04\x43ity\"\n\n\x08\x42uilding\"\t\n\x07Terrain\"\x0c\n\nCityObject\"\x0e\n\x0c\x42uildingPart\"L\n\x0bRoadNetwork\x12\x10\n\x08vertices\x18\x01 \x03(\x02\x12\x0b\n\x03\x64im\x18\x02 \x01(\r\x12\r\n\x05\x65\x64ges\x18\x03 \x03(\r\x12\x0f\n\x07lengths\x18\x04 \x03(\x02\"\x1b\n\x07Landuse\x12\x10\n\x08landuses\x18\x01 \x01(\r\"\\\n\x06\x42ounds\x12\x0c\n\x04xmin\x18\x01 \x01(\x02\x12\x0c\n\x04ymin\x18\x02 \x01(\x02\x12\x0c\n\x04zmin\x18\x03 \x01(\x02\x12\x0c\n\x04xmax\x18\x04 \x01(\x02\x12\x0c\n\x04ymax\x18\x05 \x01(\x02\x12\x0c\n\x04zmax\x18\x06 \x01(\x02\"(\n\tTransform\x12\x0b\n\x03srs\x18\x01 \x01(\t\x12\x0e\n\x06\x61\x66\x66ine\x18\x02 \x03(\x02\"+\n\nLineString\x12\x10\n\x08vertices\x18\x01 \x03(\x02\x12\x0b\n\x03\x64im\x18\x02 \x01(\x05\"9\n\x0fMultiLineString\x12&\n\x0cline_strings\x18\x01 \x03(\x0b\x32\x10.DTCC.LineString\"\xd0\x03\n\x08Geometry\x12\x1c\n\x06\x62ounds\x18\x01 \x01(\x0b\x32\x0c.DTCC.Bounds\x12\"\n\ttransform\x18\x02 \x01(\x0b\x32\x0f.DTCC.Transform\x12\x1b\n\x06\x66ields\x18\x03 \x03(\x0b\x32\x0b.DTCC.Field\x12 \n\x07surface\x18\x04 \x01(\x0b\x32\r.DTCC.SurfaceH\x00\x12+\n\rmulti_surface\x18\x05 \x01(\x0b\x32\x12.DTCC.MultiSurfaceH\x00\x12\'\n\x0bpoint_cloud\x18\x06 \x01(\x0b\x32\x10.DTCC.PointCloudH\x00\x12\x1a\n\x04mesh\x18\x07 \x01(\x0b\x32\n.DTCC.MeshH\x00\x12\'\n\x0bvolume_mesh\x18\x08 \x01(\x0b\x32\x10.DTCC.VolumeMeshH\x00\x12\x1a\n\x04grid\x18\t \x01(\x0b\x32\n.DTCC.GridH\x00\x12\'\n\x0bvolume_grid\x18\n \x01(\x0b\x32\x10.DTCC.VolumeGridH\x00\x12\'\n\x0bline_string\x18\x0b \x01(\x0b\x32\x10.DTCC.LineStringH\x00\x12\x32\n\x11multi_line_string\x18\x0c \x01(\x0b\x32\x15.DTCC.MultiLineStringH\x00\x42\x06\n\x04type\"L\n\x07Surface\x12\x10\n\x08vertices\x18\x01 \x03(\x02\x12\x0e\n\x06normal\x18\x02 \x03(\x02\x12\x1f\n\x05holes\x18\x03 \x03(\x0b\x32\x10.DTCC.LineString\"/\n\x0cMultiSurface\x12\x1f\n\x08surfaces\x18\x01 \x03(\x0b\x32\r.DTCC.Surface\"s\n\nPointCloud\x12\x0e\n\x06points\x18\x01 \x03(\x02\x12\x16\n\x0e\x63lassification\x18\x02 \x03(\r\x12\x11\n\tintensity\x18\x03 \x03(\r\x12\x15\n\rreturn_number\x18\x04 \x03(\r\x12\x13\n\x0bnum_returns\x18\x05 \x03(\r\"\'\n\x04Mesh\x12\x10\n\x08vertices\x18\x01 \x03(\x02\x12\r\n\x05\x66\x61\x63\x65s\x18\x02 \x03(\r\"-\n\nVolumeMesh\x12\x10\n\x08vertices\x18\x01 \x03(\x02\x12\r\n\x05\x63\x65lls\x18\x02 \x03(\r\"C\n\x04Grid\x12\r\n\x05width\x18\x01 \x01(\x05\x12\x0e\n\x06height\x18\x02 \x01(\x05\x12\r\n\x05xstep\x18\x03 \x01(\x02\x12\r\n\x05ystep\x18\x04 \x01(\x02\"g\n\nVolumeGrid\x12\r\n\x05width\x18\x01 \x01(\x05\x12\x0e\n\x06height\x18\x02 \x01(\x05\x12\r\n\x05\x64\x65pth\x18\x03 \x01(\x05\x12\r\n\x05xstep\x18\x04 \x01(\x02\x12\r\n\x05ystep\x18\x05 \x01(\x02\x12\r\n\x05zstep\x18\x06 \x01(\x02\"U\n\x05\x46ield\x12\x0c\n\x04name\x18\x01 \x01(\t\x12\x0c\n\x04unit\x18\x02 \x01(\t\x12\x13\n\x0b\x64\x65scription\x18\x03 \x01(\t\x12\x0e\n\x06values\x18\x04 \x03(\x02\x12\x0b\n\x03\x64im\x18\x05 \x01(\x05\"c\n\x06Raster\x12\x0c\n\x04name\x18\x01 \x01(\t\x12\x0c\n\x04unit\x18\x02 \x01(\t\x12\x13\n\x0b\x64\x65scription\x18\x03 \x01(\t\x12\x0e\n\x06values\x18\x04 \x03(\x02\x12\x18\n\x04grid\x18\x05 \x01(\x0b\x32\n.DTCC.Gridb\x06proto3')

_globals = globals()
_builder.BuildMessageAndEnumDescriptors(DESCRIPTOR, _globals)
_builder.BuildTopDescriptorsAndMessages(DESCRIPTOR, 'dtcc_pb2', _globals)
if not _descriptor._USE_C_DESCRIPTORS:
  DESCRIPTOR._loaded_options = None
  _globals['_OBJECT_GEOMETRYENTRY']._loaded_options = None
  _globals['_OBJECT_GEOMETRYENTRY']._serialized_options = b'8\001'
  _globals['_OBJECT']._serialized_start=21
  _globals['_OBJECT']._serialized_end=503
  _globals['_OBJECT_GEOMETRYENTRY']._serialized_start=432
  _globals['_OBJECT_GEOMETRYENTRY']._serialized_end=495
  _globals['_CITY']._serialized_start=505
  _globals['_CITY']._serialized_end=511
  _globals['_BUILDING']._serialized_start=513
  _globals['_BUILDING']._serialized_end=523
  _globals['_TERRAIN']._serialized_start=525
  _globals['_TERRAIN']._serialized_end=534
  _globals['_CITYOBJECT']._serialized_start=536
  _globals['_CITYOBJECT']._serialized_end=548
  _globals['_BUILDINGPART']._serialized_start=550
  _globals['_BUILDINGPART']._serialized_end=564
  _globals['_ROADNETWORK']._serialized_start=566
  _globals['_ROADNETWORK']._serialized_end=642
  _globals['_LANDUSE']._serialized_start=644
  _globals['_LANDUSE']._serialized_end=671
  _globals['_BOUNDS']._serialized_start=673
  _globals['_BOUNDS']._serialized_end=765
  _globals['_TRANSFORM']._serialized_start=767
  _globals['_TRANSFORM']._serialized_end=807
  _globals['_LINESTRING']._serialized_start=809
  _globals['_LINESTRING']._serialized_end=852
  _globals['_MULTILINESTRING']._serialized_start=854
  _globals['_MULTILINESTRING']._serialized_end=911
  _globals['_GEOMETRY']._serialized_start=914
  _globals['_GEOMETRY']._serialized_end=1378
  _globals['_SURFACE']._serialized_start=1380
  _globals['_SURFACE']._serialized_end=1456
  _globals['_MULTISURFACE']._serialized_start=1458
  _globals['_MULTISURFACE']._serialized_end=1505
  _globals['_POINTCLOUD']._serialized_start=1507
  _globals['_POINTCLOUD']._serialized_end=1622
  _globals['_MESH']._serialized_start=1624
  _globals['_MESH']._serialized_end=1663
  _globals['_VOLUMEMESH']._serialized_start=1665
  _globals['_VOLUMEMESH']._serialized_end=1710
  _globals['_GRID']._serialized_start=1712
  _globals['_GRID']._serialized_end=1779
  _globals['_VOLUMEGRID']._serialized_start=1781
  _globals['_VOLUMEGRID']._serialized_end=1884
  _globals['_FIELD']._serialized_start=1886
  _globals['_FIELD']._serialized_end=1971
  _globals['_RASTER']._serialized_start=1973
  _globals['_RASTER']._serialized_end=2072
# @@protoc_insertion_point(module_scope)
