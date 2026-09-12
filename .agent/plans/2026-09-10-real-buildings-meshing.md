# Real building evidence and region-preserving triangulation

Status: completed, 10 September 2026
Authority: DESIGN.md, model-contract.md, and the user's request to continue the
generic, efficient, intuitive, independently versioned DTCC Model work.
The earlier buildings-profile milestone is completed and remains historical.

## Acceptance boundary

Audit public CityJSON samples unchanged, recording provenance and unsupported
content. Complete the ordinary `building.lod2.mesh()` workflow for Building_1
within the unchanged Montréal tutorial sample, and report Building_2's invalid
polygon rings explicitly. Preserve roof/wall/ground region identity, attributes
and membership when polygons become triangles; save/read the resulting canonical
Mesh. No new model classes, wire version or production dependency.

## Decisions

- Audit the published two-building Montréal tutorial, full Montréal VM05 tile,
  and 3DBAG 9-284-556 tile. Download only public source data to temporary files;
  do not vendor datasets or strip unsupported source content to claim success.
- Use the existing per-surface meshers and merge path, recording triangle offsets
  before attaching copied regions. Keep LinkML and semantic vocabulary out of
  meshing. Preserve the enclosing geometry transform and do not mutate input.
- Support ordinary triangulation for regions. Reject cleaning, welding/snapping,
  nested transform interpretation and field interpolation when they would require
  mappings this slice does not implement. Protect the batch meshing entry point
  and both backend routes from silently dropping regions. Only dtcc_mesher is
  installed locally; no Triangle runtime verification is claimed.
- Record Solid/shell topology, fractional/repeated LoDs and appearance as actual
  data requirements. Do not disguise a Solid as a MultiSurface or round LoDs.
- Source audit found two zero-area triangles in tutorial Building_2. Canonical
  preservation remains valid, but polygon operations reject that geometry. Do not
  repair it or report the entire tutorial as successfully meshed/exported.

## Checkpoints

- [x] Source audit, checksums, public load results and next representation decisions.
- [x] Region-preserving meshing through existing public single/batch entry points.
- [x] Focused regression for grouped polygons, holes, unclassified faces, copied
      metadata, input immutability and explicit unsupported-operation failure.
- [x] Real sample canonical workflow, timings and affected regressions.

## Verification

Use an offline sandbox command taking explicit local CityJSON paths. Record
unchanged source hashes and support counts; do not make tests network-dependent.
Exercise `MultiSurface.mesh()` and batch meshing with each installed backend.
Check region areas/membership against source polygons and exact canonical reload.
Run focused builder and CityJSON/model regressions; preserve earlier tests/work.

## Completion evidence

`docs/design/real-buildings.md` records source URLs, checksums, observed coverage,
limitations, measurements, reproduction and next representation decisions.

- `audit_real_buildings.py` completed on all three unchanged sources. Full tutorial
  canonical round trip was exact; Building_1 produced 110 triangles with region
  counts 9/61/40 and matching semantic areas. Its Mesh round trip was exact and
  input unchanged. Building_2 meshing and whole-tutorial strict CityJSON export
  were explicitly rejected for collapsed rings, without modifying the source.
- All three sources passed the official bundled CityJSON 2.0.2 JSON Schema;
  this establishes structural schema validation, not geometric certification.
- `../venv/bin/python -m pytest tests/model tests/io tests/datasets tests/reproject
  tests/builder/test_meshing.py tests/builder/test_semantic_meshing.py -q` passed
  **1,283 tests**, with 1 skipped and 53 deselected. Five focused semantic meshing
  tests cover the named invariants/failures. Only dtcc_mesher was installed.
- `git diff --check` and compileall on the edited Python modules passed.
- Tutorial model: 113,545 canonical bytes; three-run warmed medians 10.49 ms encode
  and 16.65 ms decode. Building_1 meshing took 18.39 ms in one run. No city-scale
  claim; source vertex sharing versus native per-surface storage is explicitly
  documented. No schema/base-class/wire-version or dependency changes were needed.

## Independent-agent handoff

Implement `.agent/plans/2026-09-10-real-buildings-meshing.md` through its checkpoints.
Keep the plan updated as material decisions or status change, preserve unrelated
work and prior uncommitted milestones, and run its specified verification. Follow
DESIGN.md and the user's simplicity, performance, Python usability and portable
schema constraints. Do not publish data or add production dependencies.
