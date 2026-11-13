# Canonical model exchange v1

This milestone document records historical implementation evidence. For current
wire support and schema coverage, see [the native inventory](model-inventory.md)
and [standard I/O contract](standard-schema-io.md). Old wire readers and legacy
serialization paths described below have since been removed.

Historical initial subset. New writers now emit model wire version 3, adding
ordered representations and Solid after v2 polygon surfaces/regions. Readers
retain v1/v2 support. See [representations](geometry-representations.md).
See [the buildings milestone](buildings-profile.md#versioning-and-validation) for
current version selection and the expanded admission boundary. Rules below remain
applicable to v1 data; the v1-only writer/reader statements record the initial slice.

This is the first implemented subset of the [model contract](model-contract.md),
not completion of issue #85. The authority is [Core Design](../../DESIGN.md).
Its wire definition is the former `model_exchange.proto` (now superseded by [dtcc.proto](../../dtcc_core/schemas/dtcc.proto));
admission and mapping live in [exchange.py](../../dtcc_core/model/exchange.py).

## Public workflow

```python
from dtcc_core import io
from dtcc_core.datasets import load_model_package

io.save_model(city, "city.dtcc")  # Also city.save("city.dtcc").
restored = io.load_model("city.dtcc")  # Root type comes from the artifact.

# Requires city.dataset_context, ordinarily attached by a dataset call.
package = city.export("city.dtccpkg", canonical=True)
restored = load_model_package(package.path)
```

The same file functions handle every supported root. `io.load_city` and
`io.load_mesh` also accept `.dtcc` and reject an incompatible root type. A package
path without `.dtccpkg` produces a directory. In canonical mode, `format="vtu"`
on a Mesh requests a supplemental artifact alongside the canonical model.

## Supported facts and admission

Exact supported classes are Object, City, Building, BuildingPart, Mesh, Point and
Field. Other classes and subclasses fail explicitly. Object children retain their
concrete types. Geometry representations retain their normalized keys. Empty child
groups are not semantic state. Bounds are derived caches and are omitted; explicit
`Object.calculate_bounds()` refreshes geometry and descendant bounds after array
mutation. This does not introduce automatic mutation tracking.

Objects carry model-scoped, nonempty, unique IDs. Containment is a tree; cycles and
shared child instances are rejected. Optional keyword-only `semantic_type` is an
absolute URI independent of the Python class. Optional `profile_id` (absolute URI)
and `profile_version` are supplied together. Named `relations` map nonempty names
to ordered lists of IDs in this model. Duplicate and dangling targets are rejected.
`parent` and `children` are reserved: containment remains `Object.children`.
Semantic type ranges, relation cardinalities and domain-specific cycles belong to
an optional profile. Setting a profile identity does not load or enforce it.

Attributes support string-keyed dictionaries, lists, null, booleans, finite Python
floats, integers and strings. The wire uses an explicit value union: integers are
canonical decimal strings (including values outside int64), floats are doubles,
and empty containers and null remain distinct. It is not Protobuf Struct's
double-only numeric representation. Arbitrary Python objects remain unsupported.
Canonical attribute maps do not preserve key insertion order across exchange;
consumers must not rely on the order of restored keys.

Mesh arrays preserve numerical dtype, shape and C-order values. The allowed wire
dtypes are `|b1`, `|i1`, `|u1`, `<i2`, `<i4`, `<i8`, `<u2`, `<u4`, `<u8`, `<f4`,
and `<f8`. Byte order normalizes to little endian; memory strides and aliasing are
not transported. Arrays have at most two dimensions and exact byte-length checks.
Vertices are finite float coordinates of shape `(N,3)`; faces are integer triangles
with valid indices. Default empty `(0,)` vertices/faces remain supported. Nonempty
markers contain one signed or unsigned integer per face. Nonempty normals are
finite float face normals of shape `(F,3)`; vertex normals are outside this subset.
Point coordinates and Object/Geometry affine coefficients use doubles. Every
transform preserves its own SRS string and validated 4×4 affine matrix. Native
coordinates or mutated affine coefficients that cannot be represented exactly as
doubles are rejected rather than rounded.

Fields preserve name, description, units, component count, dtype, original shape
and values, including NaN measurements. Association must be explicit: `vertex` or
`face` on a Mesh, `sample` on a Point. Counts must match the owner. Standalone fields
can declare any of these three associations without an owner count check. General
temporal/scenario axes and other associations await a concrete supported workflow.

Dataset Context belongs in the package manifest. Standalone model files preserve
intrinsic state only. Context attached to a nested object, geometry or field is
rejected because per-entity provenance associations are not specified in v1.

## Wire and package boundary

`ModelFile` declares `format="dtcc-model"`, version `1` and a typed root union.
Readers reject unsupported versions, unknown fields, unspecified types, invalid
transforms, malformed arrays and invalid native invariants. Raw legacy bytes are
never used to guess a root class. Containment/attribute nesting is limited to
32 levels. The initial implementation also imposed a 256 MiB payload cap; that
arbitrary cap has been removed. Native payloads must fit in one serialized
Protobuf message, strictly smaller than 2 GiB (2,147,483,648 bytes), as required
for [cross-implementation support](https://protobuf.dev/programming-guides/proto-limits/#total-size-of-the-message).
This is an in-memory implementation; peak memory can exceed the file size.

Canonical packages use `dtcc-dataset-manifest-v3`. They include exactly one artifact
with role `canonical_model`, format `dtcc`, data kind `model`, media type
`application/vnd.dtcc.model+protobuf`, concrete `model_type` and
`model_schema_version`. The writer uses `artifacts/model.dtcc`. Supplemental files
have role `derived` and `derived_from` pointing to that artifact. Every artifact
has a size and SHA-256 digest. Root Dataset Context, including health/warnings,
is restored from the manifest.

The canonical artifact's CRS describes the canonical root's own transform SRS,
when present. Its bounds remain unset: request bounds and stored local coordinate bounds
do not establish global bounds after composing nested transforms. A derived
artifact may carry `crs` and `bounds` when its serializer computes them from its
own output: canonical `format="geojson"`, admitted only for model types whose
writer follows the RFC 7946 web profile (class attribute `web_geojson = True`),
records `crs` `EPSG:4326` and WGS84 `bounds` `[minlon, minlat, maxlon, maxlat]`,
or six values `[minlon, minlat, minz, maxlon, maxlat, maxz]` when positions have
z. For other types canonical `format="geojson"` is refused before anything is
written. Map-placed `format="png"` and `format="tif"` artifacts, admitted for a
`Raster` or a model with a `to_raster(field=None)` hook returning
`(Raster, rendering)`, leave `crs` and `bounds` unset and carry `placement` and
`rendering` instead. Canonical PNG and GeoTIFF require a Raster with a CRS;
`format="tiff"` is written as `tif`, and other raster formats are refused before
anything is written. A colour-mapped PNG of integer data requires valid values
within `[-2**53, 2**53]`, the integers float64 holds exactly, and is refused
before anything is written otherwise.
`placement` records the artifact grid: `crs`, `transform`
(affine coefficients `a, b, c, d, e, f`), positive `width` and `height`, `bounds`
in that CRS, `bounds_wgs84` and the source `nodata` as a finite number or the
string `"NaN"`. The PNG is warped to the north-up EPSG:3857 grid from
`rasterio.warp.calculate_default_transform` (nearest for integer data, bilinear
for floating data, validity mask nearest as alpha); the GeoTIFF keeps the source
grid with band unit and description. `rendering` records the written `bands`
(uint8 RGB(A) or GeoTIFF bands) or the hook's `field`, the colour-mapped `band`,
`unit`, `value_range` and `colormap`, `transparent_nodata`, and the hook's
`components` when given; re-export of a loaded package needs the export options
again. Other derived artifacts leave all four unset. The manifest version stays
`dtcc-dataset-manifest-v3`: both records are optional and omitted when absent,
and packages without them load unchanged in the current reader. Readers older
than this change reject a package carrying them (unknown manifest fields fail
validation), so packages with map-placed raster artifacts are not backward
compatible. The reader accepts a derived `crs` only as a nonempty string,
`bounds` only as four or six finite ordered numbers, and `placement` only with
finite numbers, ordered bounds and positive sizes. Dataset discovery/request
metadata remains Context, not an intrinsic artifact claim. This slice preserves
transforms; it does not settle global transform composition or computed spatial
summaries.

The reader checks manifest version and structure, safe paths, artifact integrity,
canonical type/version/CRS agreement and derived relationships. ZIPs are read
without extraction; duplicate or undeclared ZIP members fail. Limits are 4 MiB
for the manifest and 9,999 artifacts. Packages have no fixed aggregate byte cap;
the Protobuf ceiling applies only to the native model artifact. Derived artifacts
are integrity-checked in chunks without loading them into memory in full.
Applications accepting untrusted packages can pass `max_bytes` to
`load_model_package` to bound total uncompressed artifact bytes before reading
them. CityJSON file readers similarly accept `max_bytes` to bound input bytes,
including decompressed ZIP contents, independently of Protobuf. Both budgets
default to `None`; upload services retain their own configured upload limits.
File and package writers stage output before replacing a destination;
failed admission or supplemental serialization preserves existing output. This
provides local atomic replacement, not a cross-process transaction or a promise
of crash durability for the whole package.

### Tile packages

A `RasterTileCollection` (root-only kind 19) carries tile metadata and, per tile,
the package-relative reference `artifacts/tiles/<collection>/<index>-<id>.tif`
(collection and id reduced to `A-Za-z0-9_-`), never a local cache path.
`export(path, canonical=True)` writes the canonical model and copies each tile file
to its reference as an artifact with role `tile`, format `tif`, `size` and `sha256`,
and no other metadata; `format=` and the legacy package are refused. Publishing
keeps refusing: `publish()`, the exported package's `publish()` and the upload
admission (`DatasetUploadClient.upload_package`) reject a tile package before any
upload configuration is read or anything is sent. In a `.dtccpkg` the tile members
are stored uncompressed (`ZIP_STORED`); the model and manifest keep their
compression. The package is staged and moved into place, so a failed tile copy
leaves the destination untouched.

`load_model_package` binds each reference only to a manifest artifact with role
`tile` that passed the size and SHA-256 checks; a reference without such an artifact
(including one pointing at a `derived` artifact), a `tile` artifact no reference
names, and a `tile` artifact with any other metadata all fail. In a directory
package the restored tile's `path` is the artifact file. In a `.dtccpkg`, `path` is
the member name and `archive` the local archive path; nothing is extracted, and an
archive whose resolved path holds a brace is refused. `load_tile` checks the binding
again because a `RasterTile` can be built directly: a tile restored from a
standalone model file has no binding and refuses to load or export, and is never
resolved against the working directory. An archive tile needs an existing local
archive (no GDAL virtual path, no braces) and a member that passes the portable
package-path check; the reader opens `/vsizip/{archive}/{member}`, which it builds
itself, and admits the header and the memory budget before any pixel is read. A
local tile `path` is an ordinary filesystem path, opened in its absolute form; only
GDAL virtual paths (`/vsi...`) are refused. `build_mosaic` reads local cache tiles
only and refuses unbound and archive-bound tiles before opening any file. Re-export
of a collection loaded from an archive reads its members through `zipfile`. The
tiles load and re-export after the original cache is gone.

## Compatibility and migration

The original `dtcc.proto`, generated legacy bindings, `.pb`/`.pb2` routes and
default v2 package selection remain in place. They are explicitly legacy and do
not claim canonical conformance. Legacy Object/Field writers reject the new facts
they cannot preserve. Reading legacy bytes leaves these facts absent; migration
must supply required associations explicitly, never infer them from array length.

Consumers opt into `.dtcc` and `canonical=True` while this subset is evaluated.
Canonical package publication fails clearly until the catalog/upload boundary
supports v3. Do not switch defaults or retire legacy readers until required types
and Python/browser/C++ consumers have migration evidence. Unknown new fields are
rejected by v1, so an extension needs an explicit compatibility/version decision.

Derived-artifact `crs` and `bounds` are a recorded compatibility decision of the
same kind: the manifest version stays `dtcc-dataset-manifest-v3` because the
manifest fields already exist, are optional and are null when absent,
and packages without them load unchanged in the current reader. Readers older
than this change reject a package whose derived artifact carries `crs` or
`bounds` ("Supplemental artifact must identify its canonical source"), so
packages enriched with derived spatial metadata are not backward compatible.

The generated Python binding requires existing dependency `protobuf>=5.29.0,<6`.
Regenerate only this schema with the development compiler `grpcio-tools==1.71.0`:

```bash
python -m grpc_tools.protoc -I dtcc_core/schemas --python_out=dtcc_core/model dtcc_core/schemas/dtcc.proto
```

To verify the generated file by hand, hash the committed `dtcc_pb2.py`, regenerate it
with the pinned compiler and compare the two hashes; they must be equal.

`Object.Kind` carries seven root-only result-wrapper kinds: `TREE_COLLECTION` (13),
`BUILDING_COLLECTION` (14), `FOOTPRINT_COLLECTION` (15), `CALIBRATION_GRID` (16),
`FIELD_SLICE` (17), `STREAMLINE_COLLECTION` (18) and `RASTER_TILE_COLLECTION` (19).
This kind extension is a recorded compatibility decision: wire version stays **6** and
schema 0.9.0 is unchanged. The kinds are rejected anywhere below the root. A kind
retired later is marked `reserved` in `dtcc.proto` rather than reused.

`PLACE_COLLECTION` (20) is a generic object kind under the same decision: wire
version 6 and schema 0.9.0 are unchanged. A reader without kind 20 rejects a file
that contains one ("Unsupported canonical object kind 20"); files without it
load as before.

LinkML remains isolated in [the profile experiment](../../sandbox/model_profiles/README.md).
It is neither a runtime dependency nor an alternative wire authority. Python
round-trip and rejection checks are in
[test_canonical_exchange.py](../../tests/model/test_canonical_exchange.py).
No browser/C++ interoperability or representative large-city benchmark is claimed.
