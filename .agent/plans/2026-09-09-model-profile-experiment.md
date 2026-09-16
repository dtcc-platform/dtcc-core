# DTCC Model semantic-profile experiment

Status: completed (bounded experiment; production exchange remains open)
Created: 2026-09-09

## Objective and authority

Implement the bounded experiment approved in the issue #85 discussion. Follow
`DESIGN.md` and `docs/design/model-contract.md`. Determine whether an optional
LinkML profile can constrain existing DTCC objects through a small adapter, and
which model facts must become explicit before canonical exchange is implemented.
No applicable `PLAN_TEMPLATE.md` was found in the repository or workspace ancestors.

## Scope

Create a City, Building, BuildingPart and sensor represented by an existing Object,
with a projected mesh, signed markers, normals and a numerical field. Derive
containment from the actual model, supply an explicit `observes` cross-reference,
and keep proposed semantic types and field associations visibly experimental.

Compare two independently scoped LinkML profiles with a small handwritten Python
baseline. Test required/invalid attributes, unknown types/attributes, cardinality,
wrong/dangling targets, duplicate IDs, containment cycles and profile coexistence.
Distinguish schema validation from graph checks; do not claim inference support.
Record legacy exchange losses on the representative model without changing its wire
format. Keep production code/dependencies and unrelated work unchanged.

## Checkpoints

- [x] Build the native example and a small semantic projection without numerical arrays.
- [x] Run LinkML in a separate development virtualenv; record the tested versions.
- [x] Compare declarative, graph and Python-baseline outcomes on focused cases.
- [x] Record usability, coverage, limitations and the next implementation decision.

## Verification

Run the documented example-generation and experiment commands from the repository
root. Require a nonzero exit for a deliberately invalid case and for an unexpected
comparison result. Verify profile selection A/B/A has no global leakage, the same
generic adapter handles a schema-only rule/type addition, and numerical arrays
are neither copied into the projection nor modified by semantic validation.
Observe legacy precision and metadata losses separately from profile validity.
Run `git diff --check` and verify production dependencies/schema are unchanged.

## Completion evidence and decisions

Reproduction commands, observed comparison table, decision and limits are in
`sandbox/model_profiles/README.md`. Native generation succeeded. All 18 expected
comparison outcomes matched using LinkML/linkml-runtime 1.11.1 and jsonschema
4.26.0. Profile A/B/A produced valid/invalid/valid. A schema-added WeatherStation
worked through the unchanged adapter and existing generic Object class. Array and
projection mutation checks passed. The dangling-ID CLI case returned exit 1 and
identified `objects['sensor-1'].observes`; native adapter smoke checks rejected
missing semantic types, duplicate IDs and repeated containment objects.
The valid CLI case returned exit 0. Python syntax, 15 local documentation links
and anchors, `pip check`, and `git diff --check` passed. Core's shared environment
has no LinkML installation; production code, schema, dependencies and `DESIGN.md`
are unchanged. Broad Core/sibling suites were not rerun for this isolated experiment.

The legacy round-trip probe observed about 0.235 m coordinate error, omitted mesh
markers/normals, scalar shape conversion and about 6.10e-6 K field-value error.
No wire/API changes were made. LinkML is a useful optional candidate with a small
graph pass; default JSON Schema validation alone missed five graph cases. The
handwritten baseline is comparison evidence only, not a second production authority.

Next: specify native semantic identity, named references and field association for
this example, then implement its canonical exchange/package workflow with explicit
versioning and compatibility. This experiment does not complete canonical exchange
or issue #85. Production dependency adoption remains a separate decision.

The subsequent native/canonical milestone is recorded separately in
`.agent/plans/2026-09-09-canonical-model-example.md`. It replaces the example's side
inputs with native properties and reruns the experiment after canonical exchange.
The measurements above describe the original isolated experiment.

## Independent-agent handoff

Implement the plan `.agent/plans/2026-09-09-model-profile-experiment.md` through its
checkpoints. Keep this plan updated as material decisions or status change. Read
`DESIGN.md` and `docs/design/model-contract.md`, preserve unrelated work (including
the existing contract edit and untracked PNG), and run the specified verification.
Keep all prototype code and dependencies isolated from production. Complete the
decision record; do not start a production profile framework or wire migration.
