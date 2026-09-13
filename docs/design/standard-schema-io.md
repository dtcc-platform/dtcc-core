# Standard schema at save/load boundaries

The current standard DTCC semantic schema is
[dtcc.yaml](../../dtcc_core/schemas/dtcc.yaml), alongside the Protobuf wire
definition [dtcc.proto](../../dtcc_core/schemas/dtcc.proto). Both specifications
have stable paths; semantic and wire versions remain independent.
It is a self-contained LinkML file shipped inside the installed Core package.
It declares the semantic rules applied by default to canonical I/O and
[strict CityJSON I/O](cityjson-schema-io.md). Native
numerical structure remains governed by [model-contract.md](model-contract.md)
and Core's existing admission code.

Schema 0.9.0 also declares `RasterElevation` records on Terrain. Its
`geometry_id` uses `dtcc_representation_reference: Raster` to name an existing
Raster representation on the owning Object. The semantic projection checks the
local ID and carrier type using the schema's existing native bindings; LinkML
checks record shape and values. This is a narrow inline-record annotation, not a
second metadata or geometry registry. Bypass skips this semantic link check while
retaining all native numerical and persistence integrity checks.

The current [exterior-city profile](exterior-city-profile.md) defines the finite
feature and geometry coverage. Schema 0.9.0 uses the logical semantic namespace
`https://github.com/dtcc-platform/dtcc-core/schemas/model#`; schema identity and
version remain separate. These URIs identify local definitions and are not fetched.
Earlier experimental namespaces are not automatically translated. Unknown URIs
continue to receive generic native validation rather than inferred meanings.

## Public workflow

```python
from dtcc_core import io
from dtcc_core.datasets import load_model_package

city.save("city.dtcc")
city = io.load_city("city.dtcc")

city.save("unfinished.dtcc", validate_schema=False)
city = io.load_city("unfinished.dtcc", validate_schema=False)

city.export("city.dtccpkg", canonical=True)
city = load_model_package("city.dtccpkg")
city.export("unfinished.dtccpkg", canonical=True, validate_schema=False)
city = load_model_package("unfinished.dtccpkg", validate_schema=False)
```

The same keyword works in `io.save_model`, `io.load_model`, the City/Mesh/
VolumeMesh/PointCloud/Raster `.dtcc` wrappers, `Model.to_proto`/`from_proto`,
and the canonical `exchange.dumps`/`exchange.loads` APIs.
Values must be literal booleans; a string such as `"false"` fails clearly.
There is no global bypass setting. Attribute/geometry access, construction and
edits perform no schema evaluation. An incomplete intermediate model can be edited
normally; its next default save checks its current data.

Schema violations raise `ValueError` with model/property paths before an existing
file/archive is replaced or a loaded model is returned. For example:

```text
objects['building'].attributes['measured_height']: -1 is less than the minimum of 0
```

`validate_schema=False` bypasses **only semantic evaluation**. It does not bypass
array/connectivity checks, native types, containment and reference integrity,
canonical wire-version checks, size limits, valid schema declarations, package
hashes, manifest consistency or safe package paths. Bypassed canonical writes still declare
a schema version. That declaration identifies the intended contract; it is not
a claim that the data was successfully validated.

Coverage includes **canonical .dtcc files and canonical packages**, including
directory packages, and **strict CityJSON JSON/ZIP I/O**. The sole wire definition is `dtcc_core/schemas/dtcc.proto`; a `.dtcc` file is one
binary ModelFile message. `to_proto` and `from_proto` use the same boundary.
Legacy `.pb`/`.pb2` file handlers and old model layouts have been removed. Other
format adapters and default artifact-only packages are not universally validated. Continue
passing `canonical=True` for canonical packages. Passing the bypass option to a
legacy package export is rejected, rather than suggesting that it changed checks.
Strict CityJSON now follows the same semantic contract with mode-dependent flag
defaults and explicit schema-selection limits; see [its contract](cityjson-schema-io.md).
Mesh interchange, raster and point-cloud adapters remain subsequent work.

## Coverage and generic data

The user selected allowing unfamiliar classifications as generic data. The schema
records this policy in `dtcc_open_semantics: true`, with explicit class annotations
mapping native type names to schema classes.

| Native data | Standard semantic projection |
| --- | --- |
| Object | Generic object, containment and declared attributes/references. |
| City, Building, BuildingPart | Native bindings supply the corresponding validation class when no declared semantic URI matches. [Optional building attributes](building-attributes.md), including separate measured/estimated heights; sibling Building/BuildingPart schema classes and declared containment restrictions. |
| Tree | Native height and crown_radius project as height and crown_radius; other declared plant attributes retain their meaning. |
| Landuse | Native enum-name list projects as native_landuses. Source class/function/usage codes are not translated into native codes. |
| Mesh, Point, Surface, MultiSurface, Solid | Supported standalone roots; semantic regions and field metadata are evaluated. Numerical arrays, topology, transforms and shell structure stay in Core. |
| SemanticRegion | All eleven declared building surface classifications receive owner rules; Window/Door retain their optional host rules. Unfamiliar classifications use the generic region contract. |
| Field | Native name, unit, description and association metadata are exposed to schema rules, including fields on nested polygon surfaces. Values and element counts stay in Core. |

Known semantic URIs select their declared schema class. Otherwise the validation
class comes from the native binding: an unfamiliar URI on a generic Object receives
the generic Object contract; an unclassified native Building still receives the
Building rules. **Stored semantic URIs are never rewritten or inferred.** This
does not establish the meaning of an external vocabulary or infer a missing
specialized class. Explicit validation against a stricter domain profile remains
available through `SemanticProfile`.

Undeclared attributes are allowed. Unknown named Object relationships still undergo
Core's ID/reference checks, but no undeclared target-type restriction is invented.
Only declared semantic slots enter the flattened projection; generic metadata named
`id` or `parent` cannot replace native identity/containment. Native measurements
cannot be shadowed by attributes. Objects retain their full original attributes
and relations in memory and on the wire.

The projection uses namespaced internal identifiers for object/geometry/region/
field records. Its identifier slot is not a lexical constraint on the original
Object ID. Original ID validity and uniqueness remain Core checks. No coordinate
or field-value buffers are converted into LinkML objects. Plain nested polygon
surfaces without semantic metadata do not become extra validation records.

BuildingPart containment remains restricted to Building/BuildingPart when an owner
is present. A standalone BuildingPart root is valid. Building/city-feature owners
may be generic Object containers; the standard does not insist that every useful
native hierarchy be a complete CityJSON document. Known boundary regions belong
to Building/BuildingPart or a standalone Geometry. Window/Door hosts, when present,
must be wall or roof regions. Missing parents/hosts are not inferred.

Schema 0.7.0 also binds CityObject, Terrain, RoadNetwork, SensorCollection,
VehicleCollection, DeSO, PointCloud, VolumeMesh, LineString, MultiLineString,
Grid, VolumeGrid, Raster, Bounds and Transform. Their native state is encoded
by the shared codec; generic bindings do not claim new CityGML equivalences.
The complete supported/unsupported classification is in [model-inventory.md](model-inventory.md). This schema
does not claim full CityGML conformance, geometric watertightness, unrestricted
LinkML language support, or that every semantic URI must correspond to one concrete
Python class.

## Version selection and evolution

Schema 0.6.0 introduced [qualified height measurements and classification codes](qualified-values.md)
through generic identifier-free nested records. Their declared members receive
the same default evaluation at canonical and strict CityJSON boundaries; ordinary
scalar height access remains unchanged. No wire revision accompanies this change.
Schema 0.7.0 adds signed [qualified elevations and explicit 3DBAG mapping](3dbag-attribute-mapping.md).
Schema 0.8.0 declares TINRelief as a specialization of Terrain and accompanies
[strict triangular CityJSON terrain exchange](strict-tin-relief.md).

The standard schema ID is
`https://github.com/dtcc-platform/dtcc-core/schemas/model`, currently version
`0.9.0`. Earlier standard schemas 0.1.0 through 0.8.0 are preserved in Git
history, outside the checkout and runtime bundle; default
validation rejects those unavailable versions. The explicit bypass preserves
declarations and data without migration. See [naming](model-naming.md).
The schema ID is a logical identifier, not a URL fetched during I/O. The existing
experimental semantic class URIs are retained to avoid silently reclassifying data.

Canonical wire **v6** adds required `schema_id` and `schema_version` fields to the
root ModelFile envelope. They work for Object, Geometry, Field, Raster, Bounds and Transform roots and are
independent of the wire version. They are also independent of existing
`Object.profile_id/profile_version` domain-profile labels: those labels survive
exchange but do not select or replace the standard contract.

- New roots with neither schema property set use the bundled default.
- Loaded roots expose `model.schema_id` and `model.schema_version`. A subsequent
  save preserves that selection. Saving does not alter the source model's metadata.
- Only wire version 6 is supported. Older layouts are rejected, including with
  validation bypassed. No legacy baseline or migration reader remains.
- Default reads/writes reject an unavailable schema ID/version. An explicit
  bypass permits handling such data while retaining its declaration for re-save.
  Missing/malformed v6 declarations and unsupported wire versions still fail.
- The package bundles one active schema at `dtcc_core/schemas/dtcc.yaml`.
  Contract changes update its declared version and the runtime default together;
  Git retains earlier definitions. Older files are not silently validated against
  newer rules: their requested version must match the bundled declaration.

These root properties are I/O selection metadata. Only the serialized root's
properties select an artifact's schema; nesting an independently loaded mesh in a
city does not select a second schema for that representation. Numerical structure
and native facts remain independent of the schema-loading machinery.

The schema loader admits only the bundled schema identity and version,
verifies the loaded YAML's declaration and rejects imports. It performs no remote
lookup. At most eight standard profile instances are retained by the selector;
their compiled validators are reused. Model validity is never cached. Canonical
admission runs once per codec call; schema evaluation follows within that same
operation, without re-entering the public standalone validator.

## Installation and evidence

Core now declares LinkML **1.11.1** and linkml-runtime **1.11.1** as dependencies.
They are loaded lazily when evaluation is needed. An environment installed before
this dependency change needs its Core dependencies updated, normally with
`python -m pip install -e .` from the checkout. A normal wheel installation includes
the standard YAML and declares the same requirements. The development requirements
also pin jsonschema 4.26.0 for repeatable experiment comparisons.

The built macOS wheel was installed and exercised outside the source checkout.
Valid default save/load, semantic failure and explicit bypass passed using that
installed package and its bundled schema. The offline test environment reused
the existing native/tooling dependency directories without upgrading either.
This verifies the installed artifact and public entry point; it is not a fresh
network dependency-resolution test or a cross-platform certification.

The project's existing virtual environment was then updated with LinkML/runtime
1.11.1 using constraints that fixed every previously installed package version.
Installation added 54 distributions, changed no existing versions, and the 110
focused regressions passed directly in that environment. This dependency footprint
is a cost of using the full LinkML toolchain for runtime evaluation. `pip check`
still reports the pre-existing, unrelated dtcc-atlas/FastAPI version mismatch.
The isolated wheel environment additionally retains an existing h5py 3.14.0
installation below Core's already-declared 3.16.0 requirement; wheel smoke coverage
does not establish that the complete environment satisfies every Core dependency.

The original standard-schema boundary tests covered semantic failures, generic data,
standalone roots, bypass limits, version preservation, then-supported v1–v5 fixtures, package
integrity and schema-only Field metadata changes. The broader affected suite passed
**1,349 tests**, with **1 skipped and 53 deselected**. The final plan records later
focused checks and artifact evidence:
[standard-schema I/O plan](../../.agent/plans/2026-09-11-standard-schema-io.md).

For new integrated measurements, supply a current-format model/package pair and
run with an installed Core environment:

```bash
python sandbox/model_profiles/benchmark_schema_io.py /path/to/recorded-artifacts /tmp/dtcc-schema-io-benchmark
```

The script verifies the recorded complete 3DBAG model and its package payload,
prepares v6 inputs, and measures public file/package load/save with validation
enabled and bypassed. It runs three fresh sequential processes for each case,
each with one cold sample and three warm samples. No native checks are disabled.
Preparation, Python/Core imports, garbage collection and RSS snapshots are outside
timings; the first enabled operation includes lazy LinkML loading/compilation.

Before format unification, on the complete 3DBAG tile (2,222 objects, 4,443
representations and 68,399 surfaces), the installed-wheel run with schema 0.1.0 measured these medians over nine warm calls
per operation/mode:

| Public operation | Semantic bypass | Default validation | Additional time |
| --- | ---: | ---: | ---: |
| File save | 1.153 s | 1.496 s | 0.344 s / 29.8% |
| File load | 2.156 s | 2.575 s | 0.419 s / 19.4% |
| Package save | 1.462 s | 1.742 s | 0.280 s / 19.1% |
| Package load | 2.213 s | 2.626 s | 0.413 s / 18.7% |

The first enabled operation added **1.14–1.25 s** relative to its bypassed cold
counterpart, including dependency/schema setup. Reusing the loaded schema removes
that startup work from subsequent operations. These integrated costs supersede
the earlier approximately 0.25-second incremental estimate; they are complete
public I/O measurements, not timings of only the semantic backend.

Tested Python 3.12.12, macOS arm64, LinkML/runtime 1.11.1, jsonschema 4.26.0 and
the installed native dependencies reused from the existing environments. Inputs
were unchanged; benchmark reads used prepared v6 files/packages. The recorded v3
source payload is 24,504,507 bytes, SHA-256
`e8d9b7a88754bab5b1c87ba1f43bcb211498e539e4ff4d2dc871294af9d3e909`.
The source package SHA-256 is
`955d07dcf92510aa1e2c8c0d29aec79ecf6f43b71789243272cdce58009a6079`.
Evidence is in `/private/tmp/dtcc-standard-schema-evidence/benchmark/report.json`,
with raw samples and process logs alongside it. This is one workload/platform;
it does not establish cost for larger cities, heavily invalid data or other formats.

The earlier measurements above are historical; they are not measurements of the
extended schema 0.2.0. The old v3 source artifact is no longer accepted by the reader.
The benchmark now fingerprints its supplied current-format inputs and verifies
that the package contains the same model. Unified-format completion and checks
are recorded in [the plan](../../.agent/plans/2026-09-11-unified-protobuf.md).
