# Consistent model display and detailed inspection

Status: completed, 16 September 2026
Authority: user-approved display contract; issue #33; DESIGN.md.
No PLAN_TEMPLATE.md exists in this repository or the workspace. This follows
the established concise plan structure in `.agent/plans`.

## Acceptance boundary

Compact `ClassName(name=value, ...)` representations throughout the public model,
with `<...>` for summaries that omit reconstructible state. `str(obj)` and
`print(obj)` use the same representation. Small complete values may use executable
constructor expressions; dataset context or omitted intrinsic state requires the
summary form. Repr must not scan geometry, compute statistics or expand children.

Detailed `.info()` uses the existing Rich table styling where appropriate. One
model implementation assembles subclass content and optional dataset context,
prints by default, or returns plain text with `print=False`. Preserve
`presentation=False`, `print_info(file=...)`, domain statistics and useful status
messages. Dataset descriptors follow the same compact/detail split. Explicit
dataset catalogue inspection must not depend on log level. Keep `.tree()`, file
metadata dictionaries and operational logging distinct.

No new dependencies, serialization changes, network datasets, publication, or
unrelated native-binding fixes. Preserve the pre-existing untracked history PNG.
Deliver all task changes in one commit and show actual output examples.

## Checkpoints

- [x] Implement shared display/report helpers and compact model representations;
      account explicitly for dataclass-generated repr in subclasses.
- [x] Unify detailed model reports, sensor/vehicle/road/DeSO reports, dataset
      descriptor help and catalogue inspection; update affected documentation.
- [x] Verify ordinary display, return/print/file behavior, empty and populated
      models, reconstruction markers, bounded cheap repr and literal table text.
      Run focused tests and affected regression suites; self-review.
- [x] Record evidence and runnable examples for the single-commit delivery.

## Implementation decisions

Use small private summary/report hooks and the existing `make_table()` helper;
do not introduce a renderer registry, public report schema or new display API.
Limit displayed collection rows with explicit omission counts. Keep concrete
class names in inherited summaries. Reconstructible forms are a narrow opt-in
for complete finite small values, not an inference from arbitrary dataclasses.

Model tables show at most 20 rows with omission counts; dataset parameter help,
catalogue listings and attached context remain complete. Remote parameter help
uses the existing `show_options()` authority without contacting the service.
Native geometry-representation and semantic-region records use compact repr
without being promoted to Model subclasses. Native extension internals remain
outside this Python model/display change.

## Verification

Start with `tests/model/test_model_info.py` and focused dataset display tests.
Exercise all exported concrete model types, plus representative populated city,
mesh, fields, sensors and dataset context. Check a large city and point cloud
without recursive expansion or bounds calculation; test malformed field input
at an existing validating boundary. Run model, datasets (excluding live), common,
I/O and import-time tests as warranted by affected boundaries. No live provider
calls. Inspect actual reports and `git diff --check` before committing.

## Completion evidence

- Broad affected regression run: **1,322 passed, 2 skipped**, in 48.20 s:
  `.venv/bin/python -m pytest -q tests/model tests/datasets tests/common tests/io tests/test_import_time.py --ignore=tests/datasets/live --maxfail=5`.
- After final report-limit and reconstruction-marker corrections, **143 focused
  tests passed** across model display, roads, DeSO, collections, dataset context,
  registration, remote descriptors and vehicle guidance. Subsequent small
  collection-omission and native schema/landuse checks also passed.
- Runtime coverage exercises all **35 exported concrete Model types**, checks
  identical repr/str, literal markup, print/return/file equality, finite small
  value reconstruction and omitted-state markers. A 10,000-child city and
  100,000-point cloud retain compact repr without invoking child repr or bounds
  calculation. A malformed Field passed to `to_proto()` still fails clearly at
  the existing shape-validation boundary.
- `.venv/bin/python sandbox/model_display.py` ran successfully with real Point,
  Bounds, Mesh, Field, City, dataset descriptor and dataset-context examples.
  It runs offline. The contract and migration notes are in `docs/model-display.md`.
- `git diff --check` passed. No new dependencies or serialization behavior.
  Existing environment/dependency warnings remain; live provider tests were not
  run. The unrelated history PNG is excluded from the commit.

## Independent-agent handoff

Implement `.agent/plans/2026-09-16-model-display.md` through its checkpoints.
Keep the plan updated as material decisions or status change, preserve unrelated
work, and run the specified verification. Follow the acceptance boundary above.
Deliver the task in one commit and show runnable examples of the resulting
`repr`, `str`, and `.info()` output.
