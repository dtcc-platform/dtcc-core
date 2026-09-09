# DTCC Model: contract and design proposal

Status: first hardening milestone implemented for [dtcc-core #85](https://github.com/dtcc-platform/dtcc-core/issues/85).
This document separates the existing, partial exchange contract from proposed work.
It does not declare a new file format or introduce an ontology implementation.

DTCC Model should remain a small, generic representation of spatial objects,
geometry and values. Its usefulness as a city model also depends on agreed meanings
for objects and relationships. We should strengthen both: reliable storage and
exchange in the core, with a small optional city semantic profile built from the
domain classes that already exist.

## Public types and current exchange coverage

The inventory follows the exports in `dtcc_core.model` and also includes the base
`Model` and module-visible `dtcc_core.model.geometry.Polygon`. The generated `proto`
module is an implementation interface, not another model class.

In the tables, **partial PB** means that `to_proto`/`from_proto` encode the existing
schema; it does **not** promise lossless Python state, dtype or coordinate precision.
**JSON out** means inherited protobuf JSON output through `to_json`, subject to the
same omissions; there is no matching public `from_json` contract. It is different
from CityJSON and GeoJSON. **Unsupported** means explicit `NotImplementedError`.
Native I/O entries list available entry points, not whole-model round-trip claims.

| Object/value type | PB | Protobuf JSON | Native I/O beyond low-level PB methods | CityJSON |
| --- | --- | --- | --- | --- |
| `Object` | Partial PB; untyped root | JSON out | No dedicated file loader/saver | No generic object adapter |
| `City` | Partial PB | JSON out | City PB files; mesh import; footprint vector import/export | Restricted import/export |
| `CityObject` | Partial PB; type tag | JSON out | Through containing City PB | Not currently mapped |
| `Building` | Partial PB; type tag | JSON out | Through City; footprint vector adapters | Through City |
| `BuildingPart` | Partial PB; type tag | JSON out | Through City PB | Through parent Building |
| `Terrain` | Partial PB; type tag | JSON out | Through City; raster/mesh adapters | Mesh as TINRelief |
| `Tree` | Partial PB; position, height, crown radius, type tag | JSON out | `io.trees.save_trees` exports a list to vector files | Not currently mapped |
| `Landuse` | Unsupported | Unsupported | Fiona vector import | Not currently mapped |
| `RoadNetwork` | Partial PB; type tag | JSON out | Fiona vector import; dataframe conversion | Not currently mapped |
| `SensorCollection` | Partial PB as known root class; nested export rejected | JSON out as known root | Dataset PB output; arrays/plotting | Not currently mapped |
| `VehicleCollection` | Partial PB as known root class; nested export rejected | JSON out as known root | Dataset PB output; arrays/plotting | Not currently mapped |
| `DeSO` | Partial PB as known root class; nested export rejected | JSON out as known root | Dataset PB output; arrays/dataframe/plotting | Not currently mapped |
| `BuildingCollection` | Unsupported | Unsupported | List/footprint conversion | No direct adapter; use a City |
| `FootprintCollection` | Unsupported | Unsupported | GeoJSON output; arrays/Shapely conversion | No direct adapter |
| `TreeCollection` | Unsupported | Unsupported | List conversion; tree-list vector exporter | No direct adapter |
| `CalibrationGrid` | Unsupported | Unsupported | GeoJSON import/output; Python mapping conversion | No direct adapter |
| `DatasetCollection` | Unsupported | Unsupported | Transitional sequence/list conversion | No adapter |
| `DatasetValue` | Unsupported | Unsupported | Transitional Python-value conversion | No adapter |
| `Field` | Partial PB; float values and `dim` | JSON out | Through geometry/native mesh I/O where supported | No field adapter |
| `Raster` | Partial PB as separate Raster message | JSON out | PB, GeoTIFF/images, CSV; additional read formats | No adapter |

| Geometry/support type | PB | Protobuf JSON | Native I/O beyond low-level PB methods | CityJSON |
| --- | --- | --- | --- | --- |
| `Model` | Abstract interface | Delegates to concrete PB | No generic load contract | No adapter |
| `Geometry` | Abstract base; shared message fields only | Through concrete types | Through concrete types | No generic adapter |
| `Bounds` | Partial PB; float coordinates | JSON out | Through owning model | City extent metadata only |
| `Transform` | Partial PB; SRS and float affine values | JSON out | Through Geometry; Object transform omitted | Coordinate conversion, not preservation of DTCC transforms |
| `Point` | PB coordinates are double; shared metadata limitations | JSON out | Through Object PB | Not currently mapped |
| `LineString` | Partial PB | JSON out | Through RoadNetwork/native vector adapters | No direct adapter |
| `MultiLineString` | Partial PB | JSON out | Through RoadNetwork/native vector adapters | No direct adapter |
| `Surface` | Partial PB | JSON out | Footprint/Shapely conversions | Through supported city objects |
| `MultiSurface` | Partial PB | JSON out | Through objects; surface conversions | Through supported city objects |
| `Polygon` | Abstract (no bounds implementation); serializer methods explicitly unsupported | Unsupported | Shapely-backed methods; not a concrete public model | No adapter |
| `Mesh` | Partial PB; markers/normals omitted | JSON out | PB and mesh formats including OBJ/PLY/STL/VTK/VTU/XDMF; glTF output | Through supported city objects |
| `VolumeMesh` | Partial PB; markers omitted | JSON out | PB and mesh formats, including XDMF | No adapter |
| `PointCloud` | Partial PB | JSON out | PB, LAS/LAZ, CSV | No adapter |
| `Grid` | Partial PB; empty dimensions supported | JSON out | Through Geometry PB | No adapter |
| `VolumeGrid` | Partial PB; empty dimensions supported | JSON out | Through Geometry PB | No adapter |
| `FieldSlice` | Unsupported | Unsupported | Explicit GeoJSON/PNG artifacts | No adapter |
| `StreamlineCollection` | Unsupported | Unsupported | Explicit GeoJSON/PNG artifacts | No adapter |

`GeometryType`, `RoadType` and `LanduseClasses` are public enums rather than
standalone serializable models. Their interpretation belongs to the containing
model's contract. Geometry keys currently mix representation roles (`LOD0`) and
geometry kinds (`MESH`); existing enum keys and intentional custom names must
remain usable. A nonempty custom name is different from a malformed key.

The authoritative file-format dispatch tables live in
[`dtcc_core/io`](../../dtcc_core/io). Availability of a format is not evidence that
every attribute survives it. Footprint `.pb`/`.pb2` handlers now fail explicitly;
use the City protobuf I/O route for City files.
Raster is a separate protobuf message, so storing it in an Object's Geometry map
is not supported by the current wire schema even though `add_raster` exists.
Object export rejects Raster and Bounds geometry with `NotImplementedError`.

CityJSON import/export currently covers buildings, building parts and terrain.
The importer warns and skips unsupported root object types and catches some
processing errors; the writer does not emit other object classes. Thus a CityJSON
load/save cycle is not a general DTCC round trip. Define supported/lossy cases
explicitly before extending this adapter; preserve existing supported use cases.

## Phase 1: harden the existing contract

These changes preserve the protobuf schema and established public interfaces.
The test matrix exercises minimal/populated instances, nested instances
where supported, byte and message decoding, and decoding into an already populated
instance. Unsupported classes must fail explicitly, including when nested.

- [x] Implement Tree's existing protobuf representation; make Landuse and Polygon
  serialization explicitly unsupported. Landuse stores a Python list, while its
  current protobuf field is one `uint32`; silently choosing one value loses data.
- [x] Make repeated Object/Geometry decoding replace prior state, and reject invalid
  geometry keys and non-JSON attributes with actionable errors.
- [x] Fix independent Transform instances, existing Grid/VolumeGrid behavior and
  Field/Raster validation without changing their wire representations.
- [x] Complete the public-type round-trip/unsupported matrix and targeted I/O tests;
  document each represented subset and remaining omission.
- [x] Check downstream contracts, including viewer grid access and simulation
  scalar fields shaped both `(n,)` and `(n, 1)`.

The executable inventory is
[`tests/model/test_serialization_matrix.py`](../../tests/model/test_serialization_matrix.py).
It checks every exported concrete model, protobuf JSON, nested objects and geometries,
and reuse of decoded instances. Passing tests establish the documented subset;
they do not establish lossless serialization of fields absent from the schema.

Intentional boundary changes and compatibility rules:

- Geometry keys are enum members or nonempty string roles. Known enum aliases
  normalize consistently; `location`, `grid`, `volume_grid` and `centerlines`
  remain strings. Malformed reserved names, non-string keys and alias collisions
  fail. Nested objects/geometries whose wire tag would erase their concrete type
  now fail instead of exporting a misleading payload.
- Attributes must be JSON objects containing string keys and JSON values.
  Tuples, non-string mapping keys, unsupported NumPy values, cycles and nonfinite
  numbers now fail with a path. Convert tuples/arrays to lists and NumPy integers
  to Python integers explicitly; use `None` for a missing attribute.
- Field values accept `(N,)` for scalar fields and `(N, dim)` for scalar/vector
  fields; the legacy decoder returns `(N, dim)` float arrays. NaN measurements are
  permitted in Field/Raster values. Association and source dtype are not encoded.
- Empty Grid dimensions remain valid; asking for a step along a zero-cell axis
  raises `ValueError`. Transform arrays are finite affine 4-by-4 matrices, owned
  independently by each new Transform.
- Raster decoding validates dimensions, pixel counts, dtype and georeferencing.
  Legacy grid dimensions, absent channel count/dtype/georeferencing, and the old
  zero-dimensional empty raster's stray scalar have explicit compatibility paths.
  Empty rasters now emit no pixel values. Single-channel rasters decode as `(H, W)`.
  Integer pixels that float32 cannot preserve are rejected before serialization;
  use GeoTIFF for such data. Floating pixel values still use legacy float32.
- VolumeMesh offsets invalidate its bounds, and mesh bounds cover every stored
  vertex. General Object/Geometry serialization still trusts cached bounds after
  direct array or child changes; explicit `calculate_bounds()` is still needed.
  General shape/connectivity validation, point-cloud metadata ranges, and wrong
  concrete geometry decoder/oneof combinations remain follow-up work.

## Remaining exchange and model milestones

These are required design work for the broader issue, not fixes already delivered
by phase 1. Keep #85 open until its agreed acceptance criteria are met.

- [ ] **Precision and compatibility:** preserve spatial coordinates as float64.
  Most current coordinates, bounds and affine matrices use protobuf float32;
  this can lose decimetres in projected coordinate systems. Retain legacy fixture
  decoding. Never change existing float field numbers to double: the binary
  representations differ, including packed arrays. Add new fields/messages or
  introduce an explicitly versioned format. [Protobuf encoding](https://protobuf.dev/programming-guides/encoding/),
  [schema evolution](https://protobuf.dev/programming-guides/proto3/#updating).
- [ ] **Self-describing exchange:** specify a versioned root with concrete root and
  nested model identity. Preserve Object transforms, mesh markers/normals and
  per-node metadata/context. Root metadata alone cannot preserve different data
  sources attached to children. Use a small explicit set of supported types;
  avoid a plugin registry or opaque fallback payloads without a demonstrated need.
  MultiSurface and MultiLineString currently store bare child shapes, omitting
  each child's transform and fields; that needs an explicit encoding too.
- [ ] **Fields and arrays:** specify dtype, shape and association
  (vertex/edge/face/cell/sample). Preserve accepted scalar shapes or document a
  migration. Validate association against geometry at attachment/export boundaries;
  matching array length alone is ambiguous. Preserve Raster values as well as its
  dtype label. Decide handling of missing/NaN measurements separately from invalid
  geometry coordinates.
- [ ] **Collections and unsupported types:** give semantic collections distinct
  encodings and settle Landuse's per-surface classification representation. A raw
  Object envelope cannot reconstruct an untagged collection class automatically.
- [ ] **Mutable bounds and transforms:** define local/global coordinate semantics,
  containment traversal and when cached bounds become invalid after public-array
  or child mutation. Fixing individual serializers does not solve arbitrary
  in-place NumPy mutations; avoid an observer framework without evidence it is needed.
- [ ] **Artifact identity:** record the actual exported model type in artifact
  metadata, including derived products, and validate it against the payload.
  Dataset-level `python_return_type` is insufficient when exports change type.

A small explicit exchange envelope is a candidate, not an approved API. Keep
low-level `to_proto`/`from_proto` compatibility distinct from any new high-level
file API. Choose migration/version behavior before changing producer output.

## A small optional city profile

DTCC already has useful semantic classes: Building, BuildingPart, Terrain, Tree,
Landuse and RoadNetwork. What is missing is an explicit shared contract for their
meaning, permitted relationships and important attributes. Keep generic Object,
Geometry and Field usable for data outside that vocabulary.

Proposed first step: document a DTCC city profile and implement only rules needed
by existing workflows using ordinary Python validation functions. For example,
within this profile a `BuildingPart` must be contained by a `Building` or another
`BuildingPart`; placing it directly under `City` would produce an error naming the
part and containment path. The generic Object container can still represent other
hierarchies. Agree the rule before enforcing it; no profile is implemented here.

Treat `children` as containment/parts. Do not overload it with unrelated meanings
such as a sensor observing a building. Add a named cross-reference only when an
actual workflow needs it, with validation of referenced IDs. Likewise, roof/wall
labels describe geometric primitives; they are separate from object type and from
numeric simulation fields. Start with vocabulary, a few rules and examples, not
an exhaustive ontology hierarchy, inference engine or configurable schema framework.

This direction follows the separation in CityGML 3.0 between a semantic conceptual
model and its encodings. [OGC CityGML overview](https://www.ogc.org/standards/citygml/).
CityJSON provides prescribed object types, containment and geometry restrictions,
while allowing flexible ordinary attributes and extensions. It separately models
the semantics of geometric primitives. These are useful reference contracts for
adapters; DTCC need not adopt every class to benefit from them.
[CityJSON specification, sections 2, 3.3 and 8](https://www.cityjson.org/specs/2.0.2/).

## Cross-package rollout

The following references are to sibling repositories in the DTCC workspace.

| Consumer | Observed dependency | Migration/check |
| --- | --- | --- |
| `dtcc-viewer` | `src/dtcc_viewer/opengl/wrp_city.py` reads literal `grid`/`volume_grid`; `scripts/main.py` writes them. Other wrappers use GeometryType enums and iterate `.fields`. | Update reads/writes together if core canonicalizes grid names. Later teach field consumers explicit association; preserve scalar display and LOD behavior. |
| `dtcc-sim` | `dtcc_sim/urban_wind.py` writes vector `(n,3)` and scalar `(n,)` fields; `smooth_reconstruction.py` writes `(n,1)` scalars. `datasets.py` emits raw traffic-result protobuf bytes. | Run `tests/test_dtcc_core_contract.py`, field/traffic tests; migrate producers to explicit vertex association and any agreed exchange envelope together. |
| `dtcc-atlas` | `server/jobs/worker.py` passes dataset bytes through; upload/catalog paths consume artifact manifests. | No model parser migration for phase 1. Verify new manifest/type metadata and existing download behavior when exchange changes. |
| `dtcc-tangible-twin` | `scripts/generate_table_catalog.py` consumes DatasetArtifact manifests. | Check catalog compatibility with added artifact type/version information. |
| `dtcc-mesher` | Targeted Python-package audit found no direct model protobuf, geometry-key or field dependency. | No immediate change identified; retain core meshing integration checks. |

Roll out in order: phase 1 regression fixes; agreed schema/version and legacy
fixtures; core readers/writers and artifact metadata; affected simulation/viewer
producers and consumers; remote/catalog integration tests. Discuss and implement
the optional city profile separately so ontology decisions do not delay basic
serialization correctness.

## Verification of this milestone

Commands run from the named repository, using the shared workspace virtualenv:

| Repository | Command | Result |
| --- | --- | --- |
| dtcc-core | `../venv/bin/python -m pytest tests/model tests/io tests/reproject tests/datasets -q` | 1121 passed, 1 skipped, 50 live tests deselected |
| dtcc-core | `../venv/bin/python -m pytest tests/builder/test_builder_datamodel.py tests/builder/test_gridfield.py -q` | 3 passed |
| dtcc-sim | `../venv/bin/python -m pytest tests/test_dtcc_core_contract.py tests/test_traffic.py -q` | 12 passed |
| dtcc-viewer | `../venv/bin/python -m pytest tests/test_roadnetwork_wrapper.py tests/test_deso_wrapper.py -q` | 7 passed |
| dtcc-core | `../venv/bin/python -m compileall -q dtcc_core/model dtcc_core/io/footprints.py` | Passed |
| dtcc-core | `git diff --check` | Passed |
| dtcc-core | `../venv/bin/python -m black --version` | Unavailable: no module named black; formatting inspected manually |

No protobuf schema, generated bindings or production dependencies changed.
No sibling-package source changes were needed for this milestone. Independent
reviews found and prompted fixes for subclass type erasure, old empty-raster
compatibility and integer raster conversion; no blocking findings remain for
this milestone. The broader exchange limitations above remain unresolved.
