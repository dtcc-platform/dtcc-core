# DTCC Model: exchange contract and implementation status

Status: unified Protobuf wire version **6** and standard schema **0.9.0** implemented
for the bounded model contract in [dtcc-core #85](https://github.com/dtcc-platform/dtcc-core/issues/85).
Current contract reconciled on 12 September 2026; historical phase-1 evidence is
labeled separately below.

[`DESIGN.md`](../../DESIGN.md) is the authority for the target behavior. This
document makes its model-exchange requirements concrete, records implementation
gaps, and separates open implementation choices from requirements. The sole wire definition is now `dtcc.proto`, with the shared codec in
`model/exchange.py`. Legacy model layouts, old experiment readers and per-class
serializers have been removed at the user's explicit direction. A `.dtcc` file
is one ordinary binary Protobuf ModelFile message. See the complete current
[native/schema inventory](model-inventory.md) and [I/O contract](standard-schema-io.md).

DTCC Model remains the shared, typed representation of spatial objects, geometry
and values. Valid domain data accepted or produced by the platform must be
representable faithfully. Missing representations require extending the model;
application-specific DTOs or lossy presentation formats are not substitutes.
Implementation may proceed in small supported subsets, but those subsets must
not silently discard facts or redefine the platform requirement.

## Requirements and implementation status

Requirements derive from the Design sections on
[one semantic model](../../DESIGN.md#one-semantic-model),
[model exchange](../../DESIGN.md#model-exchange-contract),
[object-first workflows](../../DESIGN.md#object-first-workflows),
[Dataset packages](../../DESIGN.md#dataset-contract),
[boundary validation](../../DESIGN.md#trust-and-quality), and
[evolution](../../DESIGN.md#evolution).

**Implemented** means the stated, bounded requirement has evidence below.
**Partial** means some behavior exists but the requirement is not met in full.
**Missing** means the required exchange mechanism or policy is not implemented;
related Python properties or manifest fields do not establish it. Status describes
the current implementation, not whether a requirement is optional.

| Requirement | Status | Evidence or remaining gap |
| --- | --- | --- |
| Public types and their existing serialization support are inventoried. | Implemented | The complete native/schema inventory classifies all 36 exported Model subclasses: 27 supported concrete roots, eight unsupported wrappers/results and abstract Geometry. The test matrix checks this public boundary. |
| Public workflows use native typed DTCC models, including semantic collections, with explicit save/export boundaries. | Implemented | Existing model classes, Dataset results, and `Model.export()` provide this workflow; exchange completeness is assessed separately below. |
| Every type admitted to canonical Protobuf exchange preserves all facts needed for interpretation and computation. | Partial | The unified codec covers the 27 concrete native types listed in the inventory, including geometry, rasters and semantic collections. Eight wrappers/results are explicitly not admitted; they are not silently downcast. Intrinsic grid domains, including zero-area domains, survive refresh and exchange. Global-frame semantics, undeclared domain concepts and downstream integration remain limited by the documented scope. |
| Canonical artifacts identify their schema version and concrete root type; readers validate compatibility without dataset-name-to-decoder mappings. | Implemented | ModelFile declares format, wire version 6, concrete root and semantic schema. Readers accept only the current wire version and default to schema 0.9.0 for new roots. |
| Concrete nested types, containment and other required relationships survive exchange. | Partial | All admitted Object kinds, including SensorCollection, VehicleCollection and DeSO, preserve type when nested. IDs, containment and named references are checked. Domain meanings beyond declared schema rules remain open. |
| Geometry, topology, coordinate precision, CRS, transforms, markers and classifications survive exchange. | Partial | Arrays, Raster georeferencing, Bounds and transforms have explicit native storage. Derived caches are omitted and ordinary Grid domains are preserved. Zero-area intrinsic-grid domains, raw-envelope aggregation edge cases and global transform composition remain qualified in the spatial contract. |
| Numerical values preserve dtype, shape, units, components, missingness and association, and relevant temporal/scenario/other axes. | Partial | Field/Raster arrays preserve dtype, shape and numerical values. Explicit vertex/edge/face/cell/sample/geometry associations are checked. General temporal/scenario axes remain missing. |
| Attributes preserve the types and meanings required by platform data across consumers. | Partial | Wire 6 uses a typed value union for null, booleans, arbitrary-precision decimal integers, finite doubles, strings, lists and maps. The schema 0.9.0 browser run exercises these meanings; production consumer and C++ certification remain separate. |
| Intrinsic model facts remain authoritative in the model; Dataset Context is attached in memory and snapshotted in the package manifest. | Partial | Wire 6 stores supported intrinsic facts in the payload; v3 packages snapshot and restore root Context including health/warnings. Nested Context is rejected pending a per-entity provenance association. |
| Every semantic Dataset Package includes a canonical Protobuf model artifact; supplemental artifacts have explicit roles and relationships. | Partial | `Model.export(canonical=True)` guarantees the canonical artifact plus explicitly related supplements. Default artifact-only v2 packages remain pending canonical package adoption and consumer changes. |
| Package creation, reading and validation preserve the realization and context, verify integrity, and keep manifest summaries consistent with the canonical model. | Partial | V3 local readers reconstruct the model/Context and verify paths, hashes, sizes, model type/version and root CRS. Bounds summaries remain unset pending global transform semantics; remote consumers are not migrated. |
| Admission and exchange reject invalid required data and unsupported representations clearly; documented partial imports expose health/warnings. | Partial | The shared codec checks connectivity, field association, IDs and type/version; packages verify their independent boundary. Caches are omitted and descendant bounds can be explicitly refreshed. Strict CityJSON admits the 12 documented exterior feature types. Permissive CityJSON and other external adapters do not inherit a complete strict semantic contract. |
| Evolution has explicit reader compatibility, writer migration and retirement rules. | Implemented for this break | The user explicitly retired old formats without backward compatibility. Only the unified version 6 is read/written. Unknown fields are rejected before reconstruction can lose them. Downstream consumers must adopt dtcc.proto; none are silently migrated. |

Current implementation evidence is in the public-type inventory and linked workflow
contracts; the phase-1 checks below are historical. Relevant implementation paths include
[`Model`](../../dtcc_core/model/model.py),
[`Dataset schemas`](../../dtcc_core/datasets/schema.py), and
[`package export`](../../dtcc_core/datasets/package.py). This is a model-exchange
assessment, not certification of every platform workflow or package consumer.

The completed [exterior-city profile](exterior-city-profile.md) and
[native DEM workflow](terrain-dem.md) define the current bounded standards alignment.
Schema 0.9.0 adds exterior transport, water, vegetation cover and explicit elevation
raster interpretation while retaining generic Object and numeric carriers.

## Current implementation and earlier milestones

The current [inventory](model-inventory.md) maps every exported Model class to
its wire support, schema binding and standards coverage. It supersedes the former
legacy serialization tables. `to_proto`/`from_proto`, `.dtcc` file handlers and
canonical package payloads share admission, encoding and schema evaluation.
The standard schema is self-contained and versioned; numerical validity stays in
Core. Direct construction, edits and data access perform no schema evaluation.

The [real model workflow](real-model-workflow.md) exercises the current contract
on a complete cached 3DBAG tile, including edits, boundary meshing, attached fields,
native exchange and canonical package context preservation.

Earlier bounded milestones established [buildings](buildings-profile.md),
[representations/Solid](geometry-representations.md), [openings](building-openings.md),
[mixed-city semantics](model-semantic-coverage.md) and [explicit profiles](profile-validation.md).
Their recorded results remain historical evidence. Earlier wire layouts and legacy
compatibility promises are superseded by the unified-format decision, not supported
runtime paths. The [standard I/O contract](standard-schema-io.md) describes current
flags, schema identity and limits. The [browser example](../../sandbox/protobuf_interop/README.md)
records a passing Python → Chromium 152 JavaScript → Python run using wire 6 and
schema 0.9.0, including exact numeric values, browser edits and rejected receiver
replacement for invalid semantics, arrays and unknown fields. Production TypeScript
integration, C++ consumers and cross-browser certification remain outside that
demonstration. The current issue-closeout check is recorded in
[the closeout assessment](model-issue-85-closeout.md).

## Canonical exchange requirements

These requirements apply to the target contract; they do not describe guarantees
already established for every platform type and consumer. The
[issue-closeout assessment](model-issue-85-closeout.md) evaluates #85 against its
agreed bounded acceptance criteria and records separately scoped follow-ups.
Unimplemented whole-platform requirements are not silently claimed complete.

- **Spatial state:** preserve spatial coordinates as float64, topology, CRS,
  transforms, normals, signed markers (including terrain marker `-2`) and
  classifications. Preserve these facts on nested objects and geometry children.
  The current wire preserves array dtype and shape, including float64 coordinates,
  while Point coordinates, Bounds and affine matrices use Protobuf doubles. It
  does not reduce projected coordinates to float32. Local-coordinate interpretation,
  explicit bounds refresh after public-array/child mutation and remaining global
  composition limits are described in the [spatial contract](model-spatial-contract.md).
  Validate connectivity and relevant metadata ranges at admission/exchange
  boundaries. A stale cache must not become authoritative exported spatial state.
- **Identity and relationships:** canonical artifacts carry schema version and
  concrete root identity. Preserve nested types, IDs, containment and required
  named relationships. Validate references according to their declared meaning,
  including target existence/type, ID uniqueness within the declared scope, and
  prohibited cycles. Not every relationship is containment or must be acyclic.
  Preserve semantic collections and Landuse classifications without generic-object
  substitution or dropping values. A temporary unsupported-type error is a clear
  limitation, not fulfillment of the requirement to represent valid platform data.
- **Fields and arrays:** preserve dtype, shape, components, units, missingness and
  association with vertices, edges, faces, cells, objects or samples as applicable.
  Matching array length alone cannot establish association. Preserve Raster pixel
  values as well as their dtype label. Distinguish missing/NaN measurements from
  invalid geometry coordinates. Temporal, scenario, ensemble or other axes needed
  to interpret actual platform results must have explicit meanings and survive
  exchange; this does not require inventing every possible axis in advance.
- **Typed attributes:** define the accepted value types and their cross-language
  meaning, including numeric ranges, integer/float distinctions and null/missing
  behavior where meaningful. Reject values outside that contract clearly. The
  current native contract is the finite JSON-like typed union defined in
  `dtcc.proto`, not an unrestricted Python-object serializer. Attribute integers
  retain exact decimal identity beyond JavaScript Number and int64. Large numerical
  arrays belong in numerical payloads rather than attribute dictionaries.
- **Normalization and errors:** document any canonical normalization and preserve
  meaning and computational usefulness. Scalar Field values may have shape `(N,)`
  or `(N, 1)`; native exchange preserves the supplied shape without converting
  between them. Unsupported versions, types and required values fail clearly. Validation
  belongs at independent admission, persistence and exchange boundaries; this does
  not require continuous enforcement on every mutable Python container.

### Model facts, Dataset Context and packages

The model artifact is authoritative for facts intrinsic to interpretation or
computation: geometry, CRS/transforms, object identity/relationships, field units
and association, classifications, and relevant semantic axes. These facts must
remain usable when the canonical model artifact is read independently of a
Dataset Package.

Dataset identity, request, acquisition/processing lineage, source terms,
presentation and health remain Dataset Context, attached to the realization in
memory and snapshotted in the package manifest. A fact that is intrinsic to the
model belongs in the model even if a provider supplied it. Different child sources
require preserving the association between provenance and the relevant model
entities in the package contract; they do not imply copying full Dataset Context
into every node. The representation of that association remains open.

Every Dataset Package containing semantic DTCC data must include a canonical
Protobuf model artifact. Derived display/export formats and source/provenance
artifacts supplement it. The manifest must identify their roles, relationships,
semantic dimensions, capabilities, formats, coordinate context and integrity
information. It must describe the actual exported artifacts; dataset-level
`python_return_type` is insufficient when an export changes type. Any manifest
summary of an intrinsic model fact must agree with the canonical payload, which
remains authoritative.

Core owns package creation, validation and reading. The ordinary package workflow
must reconstruct the model and its Dataset Context, validate compatibility and
integrity, and reject contradictory required metadata. Catalog storage,
publication and audience policy remain Twin responsibilities. The current
canonical package path meets the declared local v3 boundary. Default artifact-only
v2 packages and downstream consumers have not all adopted this broader package
requirement; there is no retained legacy native model codec.

### Format evolution and retirement

The user explicitly chose a breaking replacement on 11 September 2026. The sole
`dtcc.proto` defines ModelFile version 6; old dtcc.proto payloads and experimental
wire versions 1–5 are rejected, including with semantic validation bypassed.
Model file suffixes are `.dtcc`; the old `.pb`/`.pb2` handlers have been removed.
External dataset selectors called `format="pb"` still mean Protobuf output, now
using ModelFile. Third-party Protobuf formats such as GTFS are unrelated.

Current schema 0.9.0 uses snake_case DTCC properties and semantic identifiers under
`https://github.com/dtcc-platform/dtcc-core/schemas/model#`; see the
[naming contract](model-naming.md). Building/BuildingPart.height reads the stored
measured_height (metres), returning None when missing. Builders write
estimated_height without replacing measurements; see the [height contract](building-height.md).
Geometry extent requires an explicit geometry selection and retains its coordinate units.
The schema declares [basic building attributes and eleven surface types](building-attributes.md),
with Building and BuildingPart sharing an abstract schema parent, plus the bounded
[exterior-city profile](exterior-city-profile.md) and [native DEM interpretation](terrain-dem.md).

Earlier milestones added the [qualified-value contract](qualified-values.md) in
schema 0.6.0, signed [qualified elevations and explicit 3DBAG import](3dbag-attribute-mapping.md)
in 0.7.0, and [strict triangular terrain exchange](strict-tin-relief.md) in 0.8.0.
These capabilities remain in 0.9.0. Qualified records use native attribute
dictionaries and one generic LinkML backend; they do not override scalar measurements.
The [terrain/transportation crosswalk](terrain-transportation-crosswalk.md) records
both the earlier comparison and the subsequent bounded implementation without
replacing native graph/raster carriers.
Earlier standard schemas 0.1.0 through 0.8.0 are preserved in Git history, outside
the checkout and runtime bundle. No old-name migration or alias
is applied. The explicit semantic bypass can still preserve their declarations.
Semantic versions are independent of the wire version. A valid declaration selects
a local bundled schema; an explicit bypass can preserve an unavailable semantic
version, but never permits an unknown wire layout. Future field-number changes
still require deliberate wire decisions. Unknown fields are currently rejected
because native reconstruction does not preserve them. Review future evolution
against browser and downstream consumer evidence; do not add a fallback codec.

## Acceptance evidence for canonical exchange

For each implemented subset, declare supported types and semantic properties and
prove the following through the ordinary public entry points. Expand coverage as
additional required platform types are admitted; do not label the entire model
complete based on one subset.

1. Save and load representative populated models, including nested cases, using
   the artifact's type/version. Compare the original and restored semantic state,
   not only protobuf equality. Exercise projected float64 coordinates and all
   required metadata, numerical values and relationships for the admitted types.
2. Prove focused boundary failures for invalid required input, unsupported
   type/version, and any new reference/association rules. Verify valid public
   mutations cannot export stale authoritative bounds. Unsupported required state
   must fail rather than disappear from an otherwise successful export.
3. Export and read a semantic Dataset Package through the public workflow. Verify
   the canonical artifact is present even when supplemental formats are requested,
   the model and context are restored, and mismatched metadata or failed integrity
   checks are rejected. Include a derived artifact whose type differs from its
   dataset result when that export path is supported.
4. Verify direct decoding with the published dtcc.proto and rejection of retired
   layouts. Keep representative fixtures for the supported contract. Check Python,
   browser and C++ consumers when their exchange boundary is introduced or changed;
   a Python-only result must not be reported as cross-language verification.

Use focused invariant and workflow checks rather than a Cartesian product of every
type, format and value. Historical phase-1 results prove only their then-documented
subset. The [current inventory](model-inventory.md), exterior examples and
[issue-closeout assessment](model-issue-85-closeout.md) identify the current acceptance
boundary; they do not establish complete platform or CityGML coverage.

## Open implementation choices

The requirements above are settled by Design. The initial
[canonical v1 design](canonical-model-v1.md) is historical background; `dtcc.proto`,
the standard I/O contract and the current inventory describe the implemented wire.
Remaining decisions include explicitly unsupported wrappers/results, general axes,
global transform composition and per-entity provenance associations. Choose each against a current workflow and its
acceptance evidence; avoid speculative registries or parallel models.

## Standard schema and optional domain profiles

DTCC now has an explicit, editable LinkML contract for its declared semantic types,
attributes and relationships. Schema 0.9.0 is selected by default at canonical native
and package I/O and strict CityJSON boundaries. Direct Python construction, edits
and quick access do not trigger schema evaluation. `validate_schema=False` bypasses
semantic evaluation only; native and external-format integrity checks remain.
Generic Object, Geometry and Field remain usable outside the declared vocabulary,
and unfamiliar semantic URIs retain generic data meaning.

The [isolated profile experiment](../../sandbox/model_profiles/README.md), completed
on 9 September 2026, is historical evidence: all 18 expected outcomes matched the
baseline, a schema-added semantic type needed no new native class, and two profiles
coexisted. LinkML's generated validation alone did not establish DTCC graph integrity.
The current shared backend combines schema rules with the required reference checks;
Core retains numerical and structural admission. There is no OWL inference engine
or alternative model runtime.

[Standard-schema I/O](standard-schema-io.md) and
[strict CityJSON I/O](cityjson-schema-io.md) apply the selected standard schema.
The [native validation experiment](profile-validation.md) provides an optional
explicit domain-profile API. Its disposable projection is a validation tool, not a
second wire authority. Schema selection is local and explicit; stored domain-profile
labels do not trigger fetching or replace the standard schema. Native wire state
preserves those paired profile labels; strict CityJSON rejects profile identity
without an extension mapping rather than dropping it.

The current standard allows a BuildingPart as a standalone native root; when it has
a containment parent, that parent must be Building or BuildingPart. Strict CityJSON
requires exactly one such parent and consistent source children/parents lists.
These are implemented rules, not outstanding profile proposals.

Existing typed public classes and geometry operations remain intact. New exterior
classes use generic Object and the existing geometry/region carriers, with rules in
YAML rather than a CityGML-shaped Python hierarchy. `children` means containment/parts;
named references represent other relationships where required by a workflow. Surface
labels and source grouping remain distinct from object type and numeric fields.

CityGML 3.0 separates a semantic conceptual model from its encodings, while CityJSON
prescribes a bounded external feature and geometry vocabulary. DTCC applies those
concepts without claiming complete coverage. The
[exterior-city profile](exterior-city-profile.md) is the current finite boundary;
additional themes, relationships or inference are separate scope decisions.

## Cross-package rollout

The following sibling-repository observations come from the September 2026 phase-1
audit. They are a rollout checklist, not verification of those consumers against
wire 6 and schema 0.9.0. No current consumer guarantee follows from an old passing test.

| Consumer | Observed dependency | Migration/check |
| --- | --- | --- |
| `dtcc-viewer` | `src/dtcc_viewer/opengl/wrp_city.py` reads literal `grid`/`volume_grid`; `scripts/main.py` writes them. Other wrappers use GeometryType enums and iterate `.fields`. | Update reads/writes together if core canonicalizes grid names. Later teach field consumers explicit association; preserve scalar display and LOD behavior. |
| `dtcc-sim` | `dtcc_sim/urban_wind.py` writes vector `(n,3)` and scalar `(n,)` fields; `smooth_reconstruction.py` writes `(n,1)` scalars. `datasets.py` emits raw traffic-result protobuf bytes. | Run `tests/test_dtcc_core_contract.py`, field/traffic tests; migrate producers to explicit vertex association and any agreed exchange envelope together. |
| `dtcc-atlas` | `server/jobs/worker.py` passes dataset bytes through; upload/catalog paths consume artifact manifests. | No model parser migration for phase 1. Verify new manifest/type metadata and existing download behavior when exchange changes. |
| `dtcc-tangible-twin` | `scripts/generate_table_catalog.py` consumes DatasetArtifact manifests. | Check catalog compatibility with added artifact type/version information. |
| `dtcc-mesher` | Targeted Python-package audit found no direct model protobuf, geometry-key or field dependency. | No immediate change identified; retain core meshing integration checks. |
| `dtcc-twin` | Design assigns Twin application orchestration and catalog responsibilities; the phase 1 audit did not certify its model reader. | Verify canonical type/version decoding and package consumption through the Core-owned contract, without dataset-name-to-decoder mappings. |

The browser demonstration has its own versioned evidence; production consumers and
C++ still require checks at their actual boundaries. Unsupported wrappers/results,
canonical package adoption and other adapters are separately scoped follow-ups, not
an instruction to implement every wrapper or external format to close #85. The
[inventory](model-inventory.md) records those limits. No backward-compatible legacy
serializer is retained for downstream rollout.

## Historical phase-1 verification record

The following phase 1 results were recorded on 7 September 2026 in the source
workspace, before integration into current `develop`. They are historical evidence,
not results rerun by this document revision. Commands are relative to each named
repository and use that workspace's shared virtualenv:

| Repository | Command | Result |
| --- | --- | --- |
| dtcc-core | `../venv/bin/python -m pytest tests/model tests/io tests/reproject tests/datasets -q` | 1121 passed, 1 skipped, 50 live tests deselected |
| dtcc-core | `../venv/bin/python -m pytest tests/builder/test_builder_datamodel.py tests/builder/test_gridfield.py -q` | 3 passed |
| dtcc-sim | `../venv/bin/python -m pytest tests/test_dtcc_core_contract.py tests/test_traffic.py -q` | 12 passed |
| dtcc-viewer | `../venv/bin/python -m pytest tests/test_roadnetwork_wrapper.py tests/test_deso_wrapper.py -q` | 7 passed |
| dtcc-core | `../venv/bin/python -m compileall -q dtcc_core/model dtcc_core/io/footprints.py` | Passed |
| dtcc-core | `git diff --check` | Passed |
| dtcc-core | `../venv/bin/python -m black --version` | Unavailable: no module named black; formatting inspected manually |

In that phase-1 milestone, no Protobuf schema, generated bindings or production
dependencies changed. Later wire/schema milestones explicitly superseded those
serialization paths. No sibling-package source changes were needed in phase 1. Independent
reviews found and prompted fixes for subclass type erasure, old empty-raster
compatibility and integer raster conversion; no blocking findings remain for
that phase-1 milestone. Current remaining limits are described above, not inferred
from this historical result.

After merge `a8e6059` on 9 September 2026, a focused local run of the four new
contract suites, `test_dataset_semantic_collections.py`,
`test_sensor_collection_pb_roundtrip.py`, and `tests/io/test_write_footprints.py`
passed: **252 tests**. `git diff --check` also passed. The broader suites and
cross-package checks above were not rerun in that integration check.
