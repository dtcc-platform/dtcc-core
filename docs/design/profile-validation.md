# An explicit, editable semantic contract

Status: experimental native validation API, 11 September 2026.
Authority: DESIGN.md and [model-contract.md](model-contract.md).

This report records the earlier opt-in experiment. The subsequent
[standard-schema I/O milestone](standard-schema-io.md) now declares LinkML as a Core
dependency and enables standard schema validation by default on canonical I/O.
The timing results and dependency-isolation discussion below describe the earlier
experiment; the new document controls current default behavior.
The later [qualified-values milestone](qualified-values.md) also extends this
backend with explicit identifier-free inlined records and string/record unions.
The rejected inline encodings described below refer to the original experiment;
entity-ID graphs inside records remain unsupported.

The LinkML schema describes what projected DTCC semantic data means and which
constraints it must satisfy. It can now be applied directly to native objects:

```python
from dtcc_core.model.profiles import SemanticProfile

profile = SemanticProfile("sandbox/model_profiles/profiles/city/0.1.0/schema.yaml")
report = profile.validate(city)
print(report.valid)
for issue in report.issues:
    print(f"{issue.path}: {issue.message}")
```

For example:

```text
objects['building-1'].attributes['measuredHeight']: -1 is less than the minimum of 0
objects['bench-1'].attributes['chargingPower']: Missing required value 'chargingPower'
```

The entry point is in [profiles.py](../../dtcc_core/model/profiles.py), in the Core
source tree. It imports LinkML only when a profile is explicitly loaded. The
shared backend is [_profile_backend.py](../../dtcc_core/model/_profile_backend.py).
The earlier isolated evaluators load that same backend without importing Core;
their native building/opening/mixed-city examples use the same projection.
There is no second live implementation of the semantic rules.

LinkML remains optional experiment tooling; this step changes no production
dependencies, model classes, wire schema, or automatic construction/load/save
behavior. Installing Core alone does not install LinkML. The versioned schemas
remain repository files selected by explicit path; automatic package discovery
and a supported dependency extra are not part of this experiment.

## What validation means

1. Core's existing canonical admission checks native structure, IDs/references,
   numerical arrays and supported geometry. Malformed state raises the existing
   `ValueError`; unsupported native types raise `NotImplementedError`. A root must
   be an Object in the current canonical subset. This reuses Core's authority
   rather than introducing another numerical or structural validator.
2. A disposable projection carries explicit semantic URIs, IDs, attributes,
   containment, named relationships, region ownership/hosts and existing typed
   Tree/Landuse facts. Numerical geometry and field arrays stay native.
3. LinkML's JSON Schema plugin checks attribute types, required values, ranges,
   enums and cardinalities. The shared graph pass reads allowed reference targets
   and inheritance from the schema. Missing/unknown/abstract semantic types and
   undeclared relationships fail. Unspecified ordinary attributes remain allowed.

The result is an immutable `ValidationReport` with the **selected** profile ID,
version and a tuple of `ValidationIssue(path, rule, message)`. `valid` means no
semantic issues at that instant. It does not certify the whole CityGML standard,
geometric solid closure, every LinkML language construct, or every possible native
class/semantic-URI combination. No semantic class is inferred from a Python class.
The selected schema can intentionally differ from stored profile labels to compare
versions; validation neither verifies those labels nor rewrites them. A saved
label is not proof that validation has occurred.

Diagnostics use `objects['id']` as a model-scoped locator, not a new Python lookup
API. Native Tree errors point to `.height`/`.crown_radius`; Landuse code errors to
`.landuses`. Containment errors point to the owner's `.children[Type][index]`.
Opening host errors point to
`.geometry['representation-id'].geometry.regions[index].parent`.

## Projection and schema boundaries

- The projection namespaces Object IDs and generated region locators to keep them
  distinct in its graph. Its `id` slot is that internal projection identifier;
  lexical constraints on original source IDs are not currently exposed separately.
  Native Object IDs themselves remain unchanged and appear in diagnostic paths.
- Native `Tree.height`, `Tree.crown_radius` and `Landuse.landuses` project as
  `height`, `crownRadius` and `nativeLanduses`. Enum names and measurements are not
  translated or inferred. These are projection bindings for existing native
  facts, not mutable duplicates in `attributes`.
- Named Object references belong in `Object.relations`, including schema-added
  relationships. Projected `parent` comes from `children`. Region `parent` means
  the owning Object; region `host` comes from its existing local parent index.
  Other region references currently fail explicitly. Colliding attributes fail
  with their native path instead of overriding an intrinsic value/reference.
- The public loader accepts a local self-contained YAML file and rejects imports.
  It does not fetch URIs or resolve schemas from model metadata. Profile ID,
  version and classes are required. A loaded instance retains its definition;
  edit the file and construct a new instance to apply changes.
- The evaluated reference encoding is class-valued ID references, optionally
  multivalued, including `any_of` class ranges. Inline references and mixed
  class/scalar alternatives are rejected. The experimental `containment: true`
  annotation is supported on the scalar `parent` slot only. General ontology
  inference, inverse consistency and unrestricted LinkML constructs are outside
  the supported contract; adding language features needs focused admission and
  validation evidence. LinkML itself notes that its default JSON Schema strategy
  does not express every LinkML construct. [LinkML validation](https://linkml.io/linkml/data/validating-data.html).

The profile and projection must evolve together when a *new native fact* needs to
be exposed. New semantic subclasses and constraints on already exposed facts can
be changed in YAML alone. That is the useful separation demonstrated here.

## Editable-contract experiment

[native_profile_example.py](../../sandbox/model_profiles/native_profile_example.py)
loads the mixed-city fixture, validates it and saves an extended native model. It
writes self-contained city schemas 0.1.1 and 0.1.2 to the chosen output directory:

- 0.1.1 adds `ChargingBench` below `CityFurniture`, with required nonnegative
  `chargingPower`. Missing power fails at the bench's attribute path; adding
  100 succeeds. The bench remains a generic Object through a canonical round trip.
- 0.1.2 adds a maximum of 50. The same native data fails the new contract and
  remains valid against the already loaded 0.1.1 contract. The original 0.1.0
  instance still validates the original mixed city and rejects the new subtype.

The existing published-in-repository profile files are unchanged. These example
versions and `example.org` identities are experimental, not released standards.

## Measured cost on the complete 3DBAG tile

Three fresh sequential processes, followed by one final verification process,
loaded the unchanged canonical artifact used in
the preceding performance milestones: **24,504,507 bytes**, SHA-256
`e8d9b7a88754bab5b1c87ba1f43bcb211498e539e4ff4d2dc871294af9d3e909`.
It contains **2,222 native objects**, **4,443 representations**, **68,399 surfaces**
and **15,554 projected records** including regions. All validations passed the
city 0.1.0 schema, SHA-256
`80cc5bd04281b5fd5444b691327b96324b0fdd41df14d7709859b44030f5e975`.

| Operation | Observed seconds |
| --- | ---: |
| Native file load, separate from validation | 2.06–2.23 |
| Profile construction, including first LinkML import | 0.65–0.69 |
| First explicit validation, including native admission and lazy schema compilation | 0.92–0.98 |
| Repeated explicit validation, twelve calls across the four processes | 0.77–0.81 |

Separate diagnostic breakdowns after the timed calls measured native admission
at **0.52–0.54 s**, projection at **0.063–0.067 s**, record flattening at **0.003–0.004 s**,
and schema plus graph checks at **0.18–0.19 s**. This is consistent with the complete public-call
measurements; it is not an API for bypassing native admission. No new optimization
machinery was justified. Reusing the profile already reuses LinkML's compiled
per-class validators.

Process RSS after repeated validation was **65–68 MiB above the loaded-model
baseline**. The process peak rose by **84–88 MiB** above that baseline. These are
whole-process deltas including tooling, retained caches and allocator effects,
not exact LinkML heap sizes. GC occurred outside timed calls. The benchmark did
not run tracing or simultaneous workload generation.

This supports explicit validation at ingestion/review/export checkpoints for this
tile size. The dependency and memory overhead are meaningful, and the experiment
does not justify validating every edit or imposing LinkML on every Core import.
The projection and full error list are materialized, not streamed; very large or
heavily invalid datasets and long-lived loading of many schema versions have not
been profiled. Reuse a loaded profile and revalidate after relevant data changes.

## Reproduction and verification

Tested: Python 3.12.12, macOS arm64, LinkML 1.11.1, linkml-runtime 1.11.1,
jsonschema 4.26.0; native NumPy 2.5.0, Protobuf 5.29.6 and Pydantic 2.13.4.
The pinned optional tooling is in
[requirements.txt](../../sandbox/model_profiles/requirements.txt).
Both Core and that tooling must be available in the interpreter running the new
API. Normal Core imports require no optional tooling.

For this local experiment the existing native environment and the earlier
isolated tooling environment were composed explicitly, without installing or
upgrading anything in either. This development-only runner reproduces that setup:

```bash
cat > /tmp/dtcc-profile-python.py <<'PY'
import pathlib, runpy, site, sys
site.addsitedir('/private/tmp/dtcc-model-profile-85-env/lib/python3.12/site-packages')
if sys.argv[1] == '-m':
    module = sys.argv[2]
    sys.argv = sys.argv[2:]
    runpy.run_module(module, run_name='__main__', alter_sys=True)
else:
    sys.argv = sys.argv[1:]
    sys.path.insert(0, str(pathlib.Path(sys.argv[0]).resolve().parent))
    runpy.run_path(sys.argv[0], run_name='__main__')
PY

PYSTOW_HOME=/tmp/dtcc-profile-pystow ../venv/bin/python /tmp/dtcc-profile-python.py sandbox/model_profiles/native_profile_example.py /tmp/dtcc-native-profile
PYSTOW_HOME=/tmp/dtcc-profile-pystow ../venv/bin/python /tmp/dtcc-profile-python.py sandbox/model_profiles/validate_native.py /tmp/dtcc-native-profile/charging-bench.dtcc /tmp/dtcc-native-profile/city-0.1.1.yaml
PYSTOW_HOME=/tmp/dtcc-profile-pystow ../venv/bin/python /tmp/dtcc-profile-python.py sandbox/model_profiles/validate_native.py /tmp/dtcc-native-profile/missing-power.dtcc /tmp/dtcc-native-profile/city-0.1.1.yaml
PYSTOW_HOME=/tmp/dtcc-profile-pystow ../venv/bin/python /tmp/dtcc-profile-python.py -m pytest tests/model/test_profiles.py -q
```

The CLI exits 0 for valid data, 1 for semantic violations and 2 when validation
cannot run because of configuration or native admission failures. The schema is
explicitly selected even when the data stores another profile version.
`PYSTOW_HOME` contains a transitive dependency's import-time cache setup.
The composed development environment is not a reproducible distribution lock or
a decision about a supported production dependency combination.

For the real tile, run `sandbox/model_profiles/profile_validation.py MODEL SCHEMA
OUTPUT_JSON` with the same interpreter in three fresh sequential processes.
It measures setup, first/repeated calls and process memory; `--breakdown` adds a
separate phase diagnostic. Evidence from this run is under
`/private/tmp/dtcc-profile-validation-evidence/`, with final measurements in
`tile-2.json`, `tile-3.json`, `tile-4.json`, and `tile-final.json`. The
[implementation plan](../../.agent/plans/2026-09-11-profile-validation.md) records
the observed verification and completion status.

The final focused native/profile/IO suite passed **101 tests** in the composed
environment. The native environment alone passed **95**, skipping the **6**
optional integration tests. The original **18** profile cases, buildings/openings/
mixed-city evaluators and their public examples also passed. CLI success,
missing-required/range violations and missing/malformed schema failures returned
the expected exit codes. Compilation and `git diff --check` passed.
