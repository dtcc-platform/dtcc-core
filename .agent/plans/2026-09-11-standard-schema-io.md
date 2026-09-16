# Standard schema at canonical save/load boundaries

Status: complete, 11 September 2026
Authority: DESIGN.md, model-contract.md and the user's approved standard-schema milestone.

## Acceptance boundary

Ship the first standard versioned DTCC schema and apply it by default to canonical
`.dtcc` file and `.dtccpkg`/directory package reads and writes, including the
existing City/Mesh wrappers. One keyword, `validate_schema=False`, bypasses only
semantic evaluation. Structural, numerical, wire-version and package integrity
checks remain mandatory. No validation on getters, edits or construction. Other
file adapters and legacy package defaults are a later integration step.

## Decisions

- Reuse the existing LinkML backend and canonical admission, with one structural
  admission per codec operation. Bundle the YAML under the installed Core package;
  use the measured LinkML tooling for evaluation, without runtime schema fetching.
- Standard schema coverage follows the current canonical subset: Object, City,
  Building, BuildingPart, Tree, Landuse; Mesh, Point, Surface, MultiSurface, Solid;
  Field and semantic regions. Numerical arrays/geometry validity remain Core's
  authority. Additional native types are still explicitly unsupported by canonical
  exchange rather than silently claimed as schema-covered.
- Preserve ordinary standalone roots and unclassified native data. The user
  selected allowing unfamiliar types as generic data. Keep domain rules for known
  classes in YAML; distinguish native generic admission from validation of an
  external vocabulary's meaning. Native bindings and this policy are explicit
  schema annotations, without rewriting stored semantic URIs.
- Add required standard schema ID/version to canonical wire v6's root envelope,
  independently of existing Object domain-profile labels. Freeze v5 evidence;
  continue reading v1–v5 using a fixed legacy schema baseline. Unknown semantic
  schema versions can only bypass evaluation explicitly; malformed declarations
  and unknown wire versions remain errors with bypass enabled.
- Retain root schema selection on loaded models for faithful re-save, including
  bypassed unknown versions. New roots select the bundled default. Never stamp a
  successful-validation claim or cache model validity. Cache only loaded schemas.
- Validate before file replacement/package publication and before returning a
  decoded model. Package metadata continues to identify its payload wire version.

## Checkpoints

- [x] Record native/schema coverage and generic-data policy; package the schema
      and declared validation dependency with a working ordinary installation.
- [x] Implement versioned codec validation and forward the single bypass keyword
      through canonical file, package and existing object-oriented wrappers.
- [x] Verify schema failures, bypass limits, old versions, unknown schema versions,
      package integrity and preservation of existing output on save failures.
- [x] Benchmark default/bypassed public file/package operations on the real tile,
      check distribution contents, run affected regressions and document limits.

## Verification

Start with public canonical save/load on a small valid/invalid city, standalone
geometry/field and an existing package fixture. Add focused regression checks for
the new semantic/default/bypass/version boundaries. Reuse earlier codec and package
integrity coverage. Exercise schema-only constraints through the existing explicit
profile API and rerun the affected optional profile examples. Validate a built
distribution contains the schema and declares the runtime needed for default I/O.
Use the unchanged local 3DBAG tile for sequential default/bypass measurements.
No remote publication, commit or unrelated edits are authorized.

## Completion evidence

Implemented in Core, with the self-contained schema at
`dtcc_core/schemas/model/0.1.0/schema.yaml`, shared projection/backend, a local
version selector, v6 envelope declarations and default validation in canonical
file/package entry points. Added `tests/io/test_standard_schema.py` and a frozen
v5 fixture; adapted older wire-version tests to keep their original boundaries.
Public workflow, schema coverage and limitations are documented in
`docs/design/standard-schema-io.md`.

Verification observed:

- Affected model, I/O, reproject, datasets and meshing suite: **1,349 passed,
  1 skipped, 53 deselected** (41.92 s).
- Final focused schema/profile/codec/representations/mixed-city/openings/CityJSON
  tests: **110 passed** (8.09 s) with the existing tooling environment. After
  installing dependencies in the project environment, the same tests passed
  directly with `../venv/bin/python -m pytest` (**110 passed**, 9.75 s).
- All five public examples (original, buildings, openings, mixed city, native
  profile) and the earlier profile comparison evaluators passed. Historical
  strict profiles keep their declared-type policy; the standard schema implements
  the user's explicit generic-data choice.
- Built the macOS arm64 wheel, verified bundled YAML and dependency metadata,
  installed it outside the checkout and exercised default save/load, semantic
  rejection before replacement, and explicit bypass. Repeated the smoke after
  the final wheel rebuild. Wheel:
  `/private/tmp/dtcc-standard-schema-dist/dtcc_core-0.9.8.dev0-cp312-cp312-macosx_26_0_arm64.whl`.
- Installed LinkML/runtime 1.11.1 in the project virtual environment while fixing
  every existing package version with a constraints file. Verified no existing
  versions changed; 54 distributions were added. The pre-existing dtcc-atlas /
  FastAPI mismatch remains the only project-environment `pip check` finding.
  The isolated wheel smoke reused native/tooling dependency directories and also
  retains h5py 3.14.0 below Core's pre-existing 3.16.0 requirement. These are not
  claims of a clean full dependency resolution or cross-platform verification.
- Integrated installed-wheel benchmark: 24 fresh sequential workers, three cold
  and nine warm samples per operation/mode, on the unchanged complete 3DBAG tile.
  Default validation added warm medians of **0.344 s file save, 0.419 s file load,
  0.280 s package save and 0.413 s package load** (18.7–29.8%). Cold additions
  were 1.14–1.25 s including lazy dependency/schema setup. These supersede the
  earlier approximate 0.25 s estimate. Raw samples and inputs are identified in
  `/private/tmp/dtcc-standard-schema-evidence/benchmark/report.json` and the
  design document. The subsequent schema edits changed descriptions only.

This completes the canonical boundary milestone. Other interchange adapters and
legacy package defaults remain subsequent work. Schema evaluation has a measurable
cost and a substantial dependency footprint; numerical checks remain native and
mandatory, and direct model access remains free of schema evaluation. No commit,
publication or unrelated cleanup was performed.

## Independent-agent handoff

Implement `.agent/plans/2026-09-11-standard-schema-io.md` through its checkpoints.
Keep the plan updated as material decisions or status change, preserve unrelated
work and all previous uncommitted milestones, and run the specified verification.
Keep the validation boundary shared, the schema local and versioned, numerical
checks mandatory, and model construction/data access free of schema evaluation.
