# Optional semantic profiles: experiment and decision

The current finite acceptance examples are:

```bash
../venv/bin/python -m sandbox.model_profiles.exterior_city_example /tmp/dtcc-exterior
../venv/bin/python -m sandbox.model_profiles.dem_example /tmp/dtcc-dem
../venv/bin/python -m sandbox.model_profiles.native_profile_example /tmp/dtcc-schema-edit
```

See [exterior coverage](../../docs/design/exterior-city-profile.md) and
[DEM interpretation](../../docs/design/terrain-dem.md). Historical profile examples
below retain their explicitly selected experimental schemas; they are not the
standard save/load contract.

Current wire definition: `dtcc_core/schemas/dtcc.proto`; standard schema: model 0.9.0.

Run `python -m sandbox.model_profiles.tin_relief_example OUTPUT` for the synthetic
mixed building/terrain strict CityJSON and native/package workflow. See
[the strict TIN contract](../../docs/design/strict-tin-relief.md).

For the audited real 3DBAG attribute convention, run
`python -m sandbox.model_profiles.three_d_bag_example SOURCE OUTPUT` from the repo.
See [the mapping contract](../../docs/design/3dbag-attribute-mapping.md) for the
pinned source, qualified elevations, source preservation and limitations.
Current optional profiles: city 0.2.0 and buildings 0.3.0; current example scripts
use their snake_case properties. Earlier profile files and recorded evidence below
are historical. See [the naming contract](../../docs/design/model-naming.md).
The old codec/readers are removed. Historical benchmark artifacts must be recreated
with the current writer before reuse; the recorded old measurements remain historical.
See [the current inventory](../../docs/design/model-inventory.md).

The [qualified-values example](../../docs/design/qualified-values.md) authors
height records and namespace-qualified codes, then verifies native/package and
strict CityJSON exchange under the current default schema. Run
`../venv/bin/python sandbox/model_profiles/qualified_values_example.py /tmp/dtcc-qualified-values`.

Current canonical behavior is documented in
[standard-schema I/O](../../docs/design/standard-schema-io.md): a bundled standard
schema runs by default, with `validate_schema=False` as an explicit semantic bypass.
LinkML/runtime are now Core dependencies. The comparison and dependency-isolation
notes below describe the earlier optional experiments, which remain reproducible.

Status: experiment completed on 9 September 2026 against Core merge `a8e6059`;
then rerun through the canonical v1 implementation in this working tree.
This is executable design evidence for [issue #85](https://github.com/dtcc-platform/dtcc-core/issues/85),
governed by [Core Design](../../DESIGN.md) and the
[model contract](../../docs/design/model-contract.md). The original comparison
below is historical experiment evidence, not a canonical exchange format or a
CityGML/CityJSON conformance profile. The subsequent
[native validation experiment](../../docs/design/profile-validation.md) adds an
explicit experimental API and real-tile measurements; LinkML remains optional.

The next [buildings milestone](../../docs/design/buildings-profile.md) has a
[self-contained versioned schema](../../schemas/profiles/buildings/0.2.0/schema.yaml)
and an end-to-end CityJSON/native/canonical example. The shared optional validator
is now `dtcc_core/model/_profile_backend.py`, loaded by the `linkml_profile.py`
development shim. The original city/low-rise files below remain comparison
fixtures; they are not the new buildings profile.

The [real-building evaluation](../../docs/design/real-buildings.md) audits unchanged
public CityJSON samples and exercises region-preserving native triangulation.
`audit_real_buildings.py` takes local paths and performs no network operations.

## Decision

Keep LinkML as the preferred candidate for optional semantic-profile authoring and
validation. The experiment demonstrates schema-only changes to constraints and
semantic types while preserving native DTCC classes. It also demonstrates that
the default validator alone is insufficient for the graph rules exercised here.
Retain a small graph-validation boundary, with domain target rules read from the
profile rather than hardcoded again in production Python.

Do not replace DTCC computational classes with generated classes or derive a new
wire schema from every profile. The [first canonical subset](../../docs/design/canonical-model-v1.md) now makes
the required intrinsic facts explicit and preserves them in a canonical exchange
workflow. This experiment makes no production dependency decision; it installed
LinkML only in a separate temporary development virtualenv.

## Representative model and adapter

`example.py` constructs real `City`, `Building`, `BuildingPart`, `Mesh`, `Field`
and `Point` instances. A sensor is an existing generic `Object`, consistent with
the sensor representation used by `SensorCollection`. The city contains a building
and sensor; the building contains its part. The part has a projected EPSG:3006 mesh,
one signed marker, a face normal, and a scalar temperature field in kelvin. Marker
`-2` exercises signed transport; this synthetic value does not classify the part
as terrain.

The adapter walks actual `Object.children` and projects IDs, absolute semantic
URIs, attributes and named ID-list relationships from native properties. `observes`
links the sensor to the building; it does not change containment. The projected
scalar `parent` comes from containment. Field association is stored on the native
Field. No semantic type or association is inferred from Python classes or lengths.
The example saves and loads `.dtcc` and a canonical `.dtccpkg`, checks exact arrays
and semantic facts, then projects the restored model for validation. Root profile
identity/version and package Dataset Context also survive.

The output JSON is a disposable experiment input, not a parallel DTCC model. Mesh
vertices and numerical field arrays stay in native DTCC objects. The adapter checks
missing semantic types, repeated containment objects and duplicate IDs before
projection, and the example verifies arrays are unchanged. The evaluator separately
checks duplicate IDs and cycles in projected data crossing the process boundary.

`city.yaml` declares the four semantic types and their attribute/reference rules.
`low-rise.yaml` imports it, limits building height to 10 m and adds `WeatherStation`
as a Sensor subtype with required nonnegative accuracy. These are illustrative
profile rules, not universal building or sensor requirements. The extension example
still uses a generic Python `Object` for that new semantic type.

## Reproduce

Run from the `dtcc-core` repository root with the existing DTCC environment. The
example uses it for native model imports; the evaluator uses a separate environment
so LinkML cannot upgrade Core dependencies. Python 3.12.12 was used for this run.

```bash
../venv/bin/python sandbox/model_profiles/example.py /tmp/dtcc-profile-example.json
../venv/bin/python -m venv /tmp/dtcc-model-profile-85-env
/tmp/dtcc-model-profile-85-env/bin/python -m pip install -r sandbox/model_profiles/requirements.txt
PYSTOW_HOME=/tmp/dtcc-profile-pystow /tmp/dtcc-model-profile-85-env/bin/python sandbox/model_profiles/evaluate.py /tmp/dtcc-profile-example.json --output /tmp/dtcc-profile-results.json
```

The install requires network access; subsequent evaluation uses local schemas and
the installed `linkml:types` import. `PYSTOW_HOME` confines a transitive dependency's
import-time cache setup to a temporary directory. LinkML is absent from Core's production dependencies and shared virtualenv.
The canonical binding raises the existing Protobuf minimum to 5.29.0. The requirements pin the evaluated validator/runtime
versions, not every transitive dependency or platform.

The full run exits 0 only when all expected outcomes match. For a single case it
acts as a validator and exits 1 for invalid data, with an object/property diagnostic:

```bash
PYSTOW_HOME=/tmp/dtcc-profile-pystow /tmp/dtcc-model-profile-85-env/bin/python sandbox/model_profiles/evaluate.py /tmp/dtcc-profile-example.json --case dangling_id
```

Observed output: `objects['sensor-1'].observes: No object 'missing'`, exit 1.
Use `--case valid` for the corresponding successful invocation.

## Observed results

Tested LinkML **1.11.1**, linkml-runtime **1.11.1**, and jsonschema **4.26.0**.
All **18 expected outcomes** matched. Both profile instances were reused in the
same process: the original example validated under city, failed under low-rise,
then validated under city again. Validation did not mutate the input projection.

| Cases | LinkML schema stage | Schema plus graph pass | Python baseline |
| --- | --- | --- | --- |
| Valid city; valid low-rise building | Accept | Accept | Accept |
| Missing height; negative height; string height; invalid usage; unknown attribute | Reject | Reject | Reject |
| Missing observes; multiple targets in single-target observes | Reject | Reject | Reject |
| Wrong observes target type; dangling ID; wrong parent type; containment cycle; duplicate ID | Accept | Reject | Reject |
| Unknown semantic type; tall building under low-rise | Reject | Reject | Reject |
| Schema-added WeatherStation with accuracy | Accept | Accept | Reject: baseline has no such type |
| Schema-added WeatherStation missing accuracy | Reject | Reject | Reject: baseline has no such type |

The schema stage maps each native semantic URI to the corresponding class URI
declared by the selected schema through SchemaView, then uses LinkML's JSON Schema plugin with unknown properties forbidden. Thus an
unknown type is rejected by the adapter before calling the plugin. The experiment
does not demonstrate automatic type inference or a universal schema loader.

The graph pass checks model-scoped ID uniqueness, reference existence, allowed
target types (including inheritance), and containment cycles. Target types come
from induced LinkML slot ranges. The `containment` annotation is an experiment
convention interpreted by our code, not a built-in LinkML acyclicity guarantee.
The `observes` ID list is constrained to exactly one target; containment has one
scalar parent. Arbitrary multivalued graphs, multiple parents and inverse
consistency have not been evaluated.

The Python baseline deliberately duplicates the original four-type profile rules
for comparison. Both backends use the same graph algorithm, with independently
supplied target declarations; this is not two independent proofs of graph-algorithm
correctness. The named graph cases exercise its required behavior directly. Only
the LinkML path reads new classes/rules from YAML. The baseline is frozen experiment
code and must not become a second production authority.

## What this establishes about model and exchange design

1. **Semantic type is separate from Python class.** A generic Object can represent
   the schema-added WeatherStation and retain its geometry operations. Native
   Objects now carry stable semantic URIs and optional profile identity/version;
   both survive the canonical wire contract.
2. **Containment and references are distinct.** `children` supplies containment;
   `observes` supplies a cross-reference. The native contract now stores named
   lists of model-scoped IDs; canonical boundaries reject duplicate/dangling IDs.
3. **Numerical validity stays in Core.** Profile success says nothing about mesh
   connectivity, CRS/transform correctness, marker meaning or field association.
   Association is now supplied on Field and validated at canonical exchange
   boundaries, including the count on its owning Mesh or Point.
4. **Profiles can share one transport.** Preserve the intrinsic types, attributes
   and references that profiles constrain. Profile changes should not each require
   new handwritten DTCC classes or a separately generated Protobuf schema.
5. **Boundary validation is sufficient for this experiment.** Operations select
   their profile explicitly, and two profiles coexist without global replacement.
   Continuous enforcement on mutable NumPy arrays/dictionaries is not implemented.

The original experiment measured the legacy City Protobuf round trip before the
native properties were added (historical evidence, no longer rerun by the example):

| Fact | Observed legacy result |
| --- | --- |
| Projected float64 coordinates | Maximum absolute error about **0.235 m** |
| Signed mesh markers and normals | Both omitted |
| Scalar field shape | Changed from `(3,)` to `(3, 1)` |
| Temperature values | Maximum absolute error about **6.10e-6 K** |
| Explicit semantic types, observes reference and field association | Experiment side inputs have no native wire representation |

The canonical rerun observed **zero coordinate and field-value error**, preserved
markers, normals, dtype and shape, and preserved native semantic facts and package
Context. Legacy losses are historical evidence, not assertions requiring loss to
persist. Canonical conformance remains limited to the explicitly supported subset.

## Next milestone and limits

The native facts and local canonical file/package milestone are implemented; see
[the exact scope and migration policy](../../docs/design/canonical-model-v1.md).
The next useful step is to exercise a representative city workflow and extend the
required type coverage before migrating package defaults or remote consumers.

The original comparison supports proceeding with LinkML as a candidate. The later
[native validation experiment](../../docs/design/profile-validation.md) now measures
memory/latency on a complete real city tile and provides an opt-in native API.
The production dependency decision remains open. RDF/OWL inference, SHACL
comparison, CityGML conformance and browser/C++ wire validation were not performed.
An alternative should be evaluated against
the same cases if a concrete requirement exceeds this approach's useful scope.

Primary tool references:
[LinkML validation](https://linkml.io/linkml/data/validating-data.html),
[JSON Schema reference encoding](https://linkml.io/linkml/generators/json-schema.html#inlining),
[SchemaView](https://linkml.io/linkml/developers/schemaview.html).


## Representations and the complete 3DBAG tile

The [implemented representation contract](../../docs/design/geometry-representations.md)
adds exact LoDs, ordered attachment records and Solid shells using canonical wire
v3. `representations_example.py SOURCE OUTPUT_DIRECTORY` checks the unchanged,
checksum-pinned tile, exercises public canonical file/package and strict CityJSON
round trips, and writes timings, size and full-workflow peak memory to `report.json`.
It uses explicit extent recomputation with original discrepant summaries retained
in package Dataset Context. The command needs only the ordinary DTCC environment.
No source geometry is repaired and no geometric validity certification is claimed.


## Real-model performance profiling

`profile_model.py SOURCE ARTIFACT_DIRECTORY OUTPUT_DIRECTORY` measures import,
codec and package operations, and CityJSON export, each in fresh sequential
processes. It records phase timings, RSS snapshots and process peak memory. Add
`--profile --repeats 1` for a separate cProfile run. See the
[performance report](../../docs/design/model-performance.md) for measured
before/after results and memory-accounting limits. The existing environment's
psutil is a development measurement aid, not a new production dependency.


`profile_model_memory.py ARTIFACT OUTPUT_JSON` measures live memory after canonical
loading and after release. Use `--copies 3` for a controlled collection of independent
copies, and run `--trace` separately for tracked live allocations and a native-state
inventory. See the [live-memory report](../../docs/design/model-live-memory.md) for
the accounting boundaries and the measured reader-local string-reuse improvement.


## Opening surfaces and hosts

[Buildings profile 0.2.0](../../schemas/profiles/buildings/0.2.0/schema.yaml) adds
Window/Door surface regions and optional wall/roof host rules. Native parent
indices, canonical v4, CityJSON normalization, version compatibility and all
reproduction commands are documented in [the openings contract](../../docs/design/building-openings.md).
`openings_example.py OUTPUT_DIRECTORY` exercises public access, file/package/
CityJSON/mesh workflows; `evaluate_openings.py PROFILE_JSON BUILDINGS_EXAMPLE_JSON`
checks the portable schema in the isolated LinkML environment. The original
0.1.0 file and examples remain unchanged in meaning and continue to validate.


## Mixed-city coverage

The [coverage map and contract](../../docs/design/model-semantic-coverage.md)
distinguish semantic-only extensions from native numerical gaps. The single
[city/0.1.0 schema](../../schemas/profiles/city/0.1.0/schema.yaml) adds land use,
individual vegetation, explicit native trees and furniture to the buildings
vocabulary. Existing building profiles remain available.

`mixed_city_example.py OUTPUT_DIRECTORY` exercises source and explicitly authored
native cities through canonical v5 files/packages, plus strict CityJSON for the
source city. `evaluate_mixed_city.py MIXED_PROFILE_JSON BUILDINGS_EXAMPLE_JSON`
checks the portable schema, older buildings, missing/invalid plant measurements
and a schema-only ChargingBench subtype. Native Tree state is not silently
converted to an unqualified source plant; see the explicit mapping boundaries
in the coverage document.
