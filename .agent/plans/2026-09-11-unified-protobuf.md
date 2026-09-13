# One DTCC Protobuf contract

Status: complete, 11 September 2026.
Authority: the user's explicit decision to retire legacy serialization without
backward compatibility, DESIGN.md and docs/design/model-contract.md.

## Acceptance boundary

Make dtcc_core/schemas/dtcc.proto the sole model wire definition. A .dtcc file is
one ordinary ModelFile Protobuf message. Remove the parallel definition, legacy
class serializers and old wire readers. Keep one shared native admission/codec
and default semantic validation, including to_proto/from_proto. Preserve existing
serializable native capabilities by porting their actual state; do not introduce
a full CityGML hierarchy. Unsupported convenience/result wrappers fail explicitly.
Retain external data-service transport contracts unless separately authorized;
never publish incompatible packages to existing consumers.

## Decisions

- Keep the current envelope version 6; accept only that version. Earlier local
  experiment versions and the old dtcc.proto layout are retired. Semantic schema
  versioning remains independent. There is no migration or fallback reader.
- Route public Protobuf methods and .dtcc file I/O through the same codec.
  Retire .pb/.pb2 model file handlers. Native arrays retain dtype/shape/raw bytes.
- Port prior concrete native types: CityObject, Terrain, RoadNetwork,
  SensorCollection, VehicleCollection, DeSO; PointCloud, VolumeMesh, LineString,
  MultiLineString, Grid, VolumeGrid; Raster, Bounds and Transform. Generic schema
  bindings do not claim CityGML semantic equivalence.
- Standard schema 0.2.0 adds native bindings; the prior semantic file 0.1.0 is
  unchanged. Geometry-wide Field association explicitly preserves DeSO area
  statistics; dataset point measurements and smoke vertex fields now declare
  their existing sample locations at their producers.
- Record a complete public-model inventory against wire support, schema binding,
  and CityGML relevance. Identify unimplemented semantic modules explicitly.
- Do not expand this step into browser bindings or all external format validation;
  make those subsequent checkpoints of the wider overhaul.

## Checkpoints

- [x] Consolidate schema/generated bindings and remove obsolete codecs/readers.
- [x] Preserve existing native data through the unified codec and public I/O.
- [x] Update focused invariants and run ordinary save/load and failure workflows.
- [x] Record inventory, semantic gaps, installed artifact evidence and next steps.

## Verification

Exercise public to_proto/from_proto and .dtcc save/load, exact numerical state,
malformed arrays, atomic failed reads/writes, schema bypass limits and packages.
Update or delete tests for intentionally retired behavior. Run affected model,
I/O, dataset and reproject tests; build and smoke the installed wheel because
schema and generated-module ownership change. Preserve all unrelated edits.

## Completion

Implemented the sole `dtcc.proto` / `dtcc_pb2.py` definition and removed
`model_exchange.proto`, its generated module, legacy per-class serializers,
unused polygon adapters, old .pb/.pb2 model handlers and v1–v5 readers/fixtures.
Model.to_proto/from_proto now use the same admitted ModelFile codec as file and
canonical package I/O; failed decoding preserves the receiver. The protobuf
package is `DTCC`, and the root message is `DTCC.ModelFile`.

Ported the listed concrete native classes, retaining array dtype/shape, exact
raster integers and channels, mesh markers, graph arrays, domain bounds and
concrete collection types. Added .dtcc Raster/PointCloud/VolumeMesh wrappers.
No new production dependency was added by this step.

The complete `docs/design/model-inventory.md` accounts for 36 exported Model
subclasses: 27 supported roots, eight unsupported result/container types and
abstract Geometry. All admitted types have explicit schema 0.2.0 bindings.
Bindings for newly admitted computational types are generic; the CityGML
property-level crosswalk is still partial. Updated the controlling contract and
I/O documentation, marked historical milestone documents accordingly, and fixed
the compiler script's missing output-directory creation and error handling.

Observed verification:

- Affected model, I/O, datasets, reproject and meshing suite: **1,314 passed,
  1 skipped, 53 deselected** in 38.32 s. The smaller count than the prior milestone
  reflects removal/replacement of obsolete legacy-format assertions.
- After public-spec formatting and test-inventory cleanup, focused public class,
  schema and unified-format checks: **74 passed** in 4.36 s.
- Added a focused malformed Field/Raster metadata check to retain protection at
  the replacement boundary: unified-format file **5 passed** in 2.98 s.
- Regenerated Python bindings from dtcc.proto. Built the macOS arm64 wheel and
  verified it includes only dtcc.proto, dtcc_pb2.py and versioned semantic schemas;
  no model_exchange module or schema remains.
- Installed the wheel outside the checkout and verified the imported package path.
  Direct generated-Protobuf decoding, schema 0.2.0 selection, exact uint64 raster
  pixels, default semantic failure, atomic file/receiver preservation, explicit
  bypass and retired-version rejection passed. The installed smoke reused the
  existing native dependency directories; it is not cross-platform or clean
  full-dependency-resolution certification.
- compileall, bash syntax check for scripts/build_proto, current documentation
  link checks and final git diff --check passed. No browser/C++ reader, external
  publication, commit or sibling-repository update was performed.

Evidence: `/private/tmp/dtcc-unified-final-tests.log`,
`/private/tmp/dtcc-unified-build.log`, and
`/private/tmp/dtcc-unified-dist/dtcc_core-0.9.8.dev0-cp312-cp312-macosx_26_0_arm64.whl`.

Remaining wider-plan work is recorded in the inventory: browser/TypeScript and
C++ interchange proof; property/relationship/units crosswalk and permanent semantic
URIs; the eight result/container decisions; canonical package/downstream adoption;
other adapter validation and representative performance checks. Existing service
format="pb" selectors now produce ModelFile bytes; they are not legacy codecs.
The older artifact-only manifest-v2 package contract remains separately supported,
and the existing gate against publishing canonical packages to unmigrated catalog
consumers remains intact. Prior performance measurements are historical, not
new measurements of schema 0.2.0.

## Independent-agent handoff

Implement `.agent/plans/2026-09-11-unified-protobuf.md` through its checkpoints,
keep it updated as material decisions or status change, preserve unrelated work,
and run the specified verification. Retire old bytes and duplicate serializers;
keep Protobuf language-neutral and distinguish schema coverage from CityGML
conformance. Do not add backward compatibility or publish external artifacts.
