# Building height authority

Status: complete, 11 September 2026.
Authority: user-approved height milestone, DESIGN.md and model contract.

## Acceptance boundary

Building/BuildingPart.height is a convenience for attributes["measured_height"],
with an equally explicit measured_height property. Missing measurements return
None, zero remains zero, access does not validate or derive geometry. Builder
outputs use estimated_height and never overwrite measured_height. Both stored
scalars are optional nonnegative metres in standard schema 0.4.0. Geometry extent
is explicit through selected geometry bounds.depth, in local coordinate units.
Existing metre-coordinate builder preconditions remain; no unit conversion,
qualified CityGML Height records, new geometry API or Protobuf layout change.
No old height-key interpretation or migration. Tree height is unchanged.

## Checkpoints

- [x] Add direct Python access and the versioned estimated-height rule.
- [x] Align point-cloud, attribute, conditioning and meshing consumers; reuse the
  existing attribute builder for the public City LoD1 path. Ground elevation must
  not become a missing building height. Preserve measurements when clamping;
  multi-building aggregates must not inherit an arbitrary source measurement.
- [x] Verify measured versus computed data through the public building workflow,
  CityJSON and default .dtcc I/O, including a rejected invalid estimate.
- [x] Update current docs, audit remaining references and record focused checks.

## Verification

Focused model/access and real City LoD1/point-cloud/conditioning tests, then the
model/I/O and affected builder/dataset suites. Check wheel schema resources if
packaging changes. Preserve unrelated work. No remote dataset/network runs.

## Completion and evidence

Implemented direct Building/BuildingPart measurement and estimate properties,
standard schema 0.4.0, aligned builder/dataset consumers, and shared public City
attribute-height processing. Geometry extent uses the existing selected geometry
bounds.depth API; no new geometry API, runtime dependency or Protobuf layout.
The previous standard 0.3.0 definition is archived byte-for-byte. Updated current
contract/inventory/naming/I/O docs and added docs/design/building-height.md with
the audit and remaining boundaries. Current browser fixture expects 0.4.0;
optional domain-profile versions remain independent and unchanged.

Observed verification:

- Public City attribute LoD1 smoke: 1 m measurement preserved while the model uses
  2.5 m above 100 m ground; both values survive default canonical I/O.
- Focused tests cover Building and BuildingPart missing/zero/no-fallback access,
  invalid-estimate save/load rejection and failed-write preservation, point-cloud
  estimates and missing-point defaults, single-source versus aggregate measurements,
  and the CityJSON/native naming round trip with separate measurement and estimate.
- Combined model/I/O, affected building/LoD1/LoD2/conditioning/meshing and dataset
  tests: **1,026 passed, 3 skipped in 37.93 s**. No live/network dataset run.
- Built wheel contains only standard model/0.4.0/schema.yaml. Archived 0.3.0 is
  byte-identical to the preceding naming wheel.
- Installed outside the checkout: verified import location, resource selection,
  public City LoD1, default .dtcc round trip and invalid-estimate failed-save
  preservation. Existing dependency directories were reused; this is not a fresh
  dependency-resolution certification.
- JavaScript syntax check, updated local documentation links and git diff --check
  passed. The browser was not rerun for this unchanged wire layout; earlier real
  browser evidence remains labelled with its original schema version.

Evidence: /private/tmp/dtcc-building-height-tests.log,
/private/tmp/dtcc-building-height-build.log, /private/tmp/dtcc-building-height-dist,
/private/tmp/dtcc-building-height-installed-smoke.py and
/private/tmp/dtcc-building-height-installed.log.

Initial test fixtures were corrected to request footprint merging with a positive
merge distance and zero minimum area. Strict CityJSON rejects native builder
representation IDs under its existing explicit mapping boundary; the public native
workflow and bounded CityJSON workflow are therefore verified separately, with the
restriction documented. No adapter expansion was added to this height milestone.
The old C++ city helper's height source is aligned, but its unrelated legacy
footprint/UUID/ground-level interface is not certified. Qualified CityGML Height
records, vertical-reference/unit conversion and default schema evaluation in external
adapters remain separate milestones. Next: basic exterior attributes and the
remaining boundary-surface declarations.

## Independent-agent handoff

Implement `.agent/plans/2026-09-11-building-height.md` through its checkpoints,
keep the plan updated as material decisions or status change, preserve unrelated
work, and run the specified verification. Keep measurement authority in the
attribute map and semantic constraints in the schema; do not infer measurements
from modelling estimates or geometry, and do not expand to qualified records.
