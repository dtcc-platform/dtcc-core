# Canonical exchange for the representative model

Status: completed (bounded local subset, 9 September 2026)
Authority: `DESIGN.md`, `docs/design/model-contract.md`, and the completed
`.agent/plans/2026-09-09-model-profile-experiment.md`.

## Scope and decisions

Implement the example's native semantic facts and a faithful, self-describing
Protobuf exchange subset: Object, City, Building, BuildingPart, Mesh, Point and
Field. Other concrete types fail explicitly. No LinkML production dependency.
Use keyword-only optional Object semantic type/profile ID/profile version and
named ID-list relations; retain children as containment. Field gets an optional
explicit association. Validate native/decoded data at canonical boundaries.

Use a separately named/versioned Protobuf message set; keep existing field numbers
and legacy decoder behavior. Preserve typed attributes, numerical dtype/shape,
float64 geometry, transforms and signed markers/normals. Bounds are derived caches,
not authority in the new payload. Reject unrepresented required state explicitly.

Provide `io.save_model` / `io.load_model` for canonical `.dtcc` files. Keep legacy
`.pb` calls explicit; legacy writers reject newly added facts they would lose.
Canonical package export is an explicit `canonical=True` transition on the existing
`Model.export` workflow. Existing package writers stay legacy during this bounded
rollout; they do not claim canonical conformance. New packages always contain the
canonical model and may include the requested supplemental format. Add Core package
reading with compatibility, integrity and intrinsic-metadata agreement checks.
No publication or remote consumer migration is authorized by this local milestone.

## Checkpoints

- [x] Native facts and precise format/migration specification.
- [x] Canonical encoding/decoding and ordinary file I/O, with focused failure tests.
- [x] Canonical package export/read through existing package authority.
- [x] Move the example off side inputs, verify profile projection after round trip,
      and update contract status and completion evidence.

## Verification

Run native-file/package workflows on the representative City. Assert exact values,
dtype, shape, metadata and semantic relationships; validate missing/dangling IDs,
bad connectivity/association, incompatible versions, unsupported types, malformed
arrays, cache mutation and package hash/metadata mismatches. Preserve a small legacy
fixture and existing model/I/O/package regression suites. Exercise the actual public
APIs early. Keep all data local; no live datasets or publication.

## Evidence and limitations

Implemented native semantic/profile URIs, model-scoped ID-list references and Field
association, the separate v1 Protobuf schema/codec and public `.dtcc` I/O, and opt-in
v3 package export/read using the existing exporter. The exact specification and
migration boundary are in `docs/design/canonical-model-v1.md`. Legacy schema and
bindings remain unchanged; the existing Protobuf runtime floor is now 5.29.0 to
match the generated binding. LinkML remains outside production dependencies.

The representative City saved/read through both public file and package APIs with
zero coordinate/value error, exact array dtype/shape, preserved markers/normals,
semantic facts and root Dataset Context. The LinkML experiment now reads those
native facts after round trip. All 18 expected outcomes and A/B/A profile
coexistence passed; the valid CLI exited 0 and dangling-reference CLI exited 1 with
an object/property diagnostic. Temporary evidence is `/private/tmp/dtcc-profile-results.json`;
reproduction commands are in the experiment README.

Verification: `../venv/bin/python -m pytest tests/model tests/io tests/datasets
tests/reproject -q` completed with **1,227 passed, 1 skipped, 53 deselected**,
including 26 focused canonical checks and fixed legacy/v1 Point fixtures.
A first reproducibility test caught a concurrent docstring edit between its two
subprocesses; with stable source the targeted checks and complete rerun passed.
Final precision review also reproduced and fixed rounding of Point coordinates
and publicly mutated integer affine coefficients outside exact double precision;
both now fail before replacement. The complete suite above was rerun afterwards.
`git diff --check`, 25 local documentation file links and byte-identical regeneration
of the new Protobuf binding passed. No live datasets, uploads or publication were
performed.

This is not completion of issue #85. Remaining limits include unsupported model
classes, general field axes, nested provenance, global transform/bounds semantics,
large-city performance and browser/C++ consumer verification. Default v2 producers
and publication remain unmigrated. No global profile registry or runtime LinkML
integration was added.

## Independent-agent handoff

Implement `.agent/plans/2026-09-09-canonical-model-example.md` through its checkpoints.
Keep the plan updated as material decisions/status change, preserve unrelated work
and the completed profile experiment, and run the specified verification. Follow
`DESIGN.md` and `docs/design/model-contract.md`. Do not silently widen the supported
subset, switch legacy producers, add production dependencies or publish artifacts.
