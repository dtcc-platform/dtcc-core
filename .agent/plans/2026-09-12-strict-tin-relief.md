# Strict triangular TINRelief exchange

Status: complete, 12 September 2026. Authority: user-approved next slice in
docs/design/terrain-transportation-crosswalk.md and the model contract.
No applicable repository plan template exists.

## Acceptance boundary

Explicit TINRelief semantics compatible with existing Terrain, and strict CityJSON
triangular CompositeSurface ↔ Mesh. Preserve all feature IDs, attributes, ordered
representations, exact LoD strings, triangle order, coordinates and CRS. Native
round trips remain exact; CityJSON comparison respects explicit quantization.
Generic Terrain is not automatically classified as TINRelief. No new geometry
classes, wire revision, dependency, raster conversion or transportation expansion.

## Checkpoints

- [x] Add standard schema 0.8.0 and bounded import/export using existing adapters.
- [x] Exercise a mixed building/two-terrain/multiple-representation example via
      public load, direct NumPy access, native and canonical package exchange.
- [x] Check nontriangles, holes, unsupported semantics/mesh state and collapsing
      quantization fail without overwriting existing output; run affected checks.
- [x] Document scope and evidence; verify the installed wheel outside checkout.

## Verification

Use a small synthetic source with shared indexed vertices, two terrain features,
two representations on one terrain and an existing building fixture. Check every
triangle and metadata field. Import compacts the document vertex table per mesh;
CityJSON global vertex numbering is not native identity. Reject unsupported mesh
fields/markers/normals and unused vertices on strict export to prevent silent loss.
No geometric validity or terrain accuracy certification is implied by a triangle
carrier. Preserve unrelated dirty-tree work and archived schema bytes.

## Completion

Implemented the strict adapter and schema vocabulary without new native classes,
wire changes or dependencies. Generic Terrain requires explicit TINRelief semantic
classification for strict export. Import preserves first-use source index sharing
in a compact Mesh. Strict export rejects unmapped mesh state and emits triangle
boundaries without generated surface labels. The existing shared vertex indexer
and quantization policy are reused.

The synthetic public example passed: building and part, two terrain features,
three terrain representations, five triangles, direct NumPy access, exact native
and canonical package round trips, preserved package Dataset Context and mixed
strict CityJSON comparison. Marker rejection leaves an existing output unchanged.
Artifacts and report: `/private/tmp/dtcc-tin-relief` (native payload 5,576 bytes).

Focused checks: 29 passed in 3.63 s. The first test invocation found a missing
DatasetContext in the new package test fixture; the fixture now supplies the
required context, with no production package change. Full model/I/O checks:
690 passed, 1 skipped in 30.28 s. No source geometry is repaired or triangulated.
Malformed later terrain input rejects the whole source; marker/normal/field/region,
unused-vertex and quantization failures protect output files, including semantic
bypass. Generic native Terrain continues to serialize without TIN inference.

A synthetic 20,000-triangle/10,201-vertex grid measured warm dictionary import/export
medians 0.4162/0.0697 s over three runs, with default semantic validation. Excludes
initial schema loading, JSON parsing/stringification and disk I/O. No non-obvious
optimization was introduced. Raw samples: `/private/tmp/dtcc-tin-timing.json`.

Built the wheel, verified final adapter/schema bytes and that only standard 0.8.0
is bundled, and installed it outside checkout. Archived 0.7.0 matches the preceding
milestone wheel exactly. Installed public strict import/access/native/package/
CityJSON and rejected-write preservation all passed. Existing dependencies were
reused without resolution. No browser or other platform certification is claimed.
Documentation links and `git diff --check` pass; unrelated dirty-tree work remains.

The resulting contract is `docs/design/strict-tin-relief.md`. The crosswalk and
current reference pages link to it; historical evidence remains labeled and pinned.
Limits: triangular CompositeSurface only, no source semantics, no raster conversion,
full GML hierarchy, transport mapping or geometric/terrain accuracy certificate.
CityJSON global vertex numbering and duplicate-vertex identity are not preserved;
ordered triangle coordinates are compared within quantization, while native arrays
remain exact. Standalone ModelFile still excludes root Dataset Context.

Logs: `/private/tmp/dtcc-tin-{focused,regression,example,timing,build,install,installed}.log`.
Wheel: `/private/tmp/dtcc-tin-dist`.
Installed smoke: `/private/tmp/dtcc-tin-installed-smoke.py`.

## Independent-agent handoff

Implement `.agent/plans/2026-09-12-strict-tin-relief.md` through its checkpoints,
keep the plan updated as material decisions or status change, preserve unrelated
work, and run the specified verification. Use existing Terrain/Mesh and the single
schema validator; finish public strict CityJSON and native/package workflows,
including meaningful failures, without expanding to raster or transport semantics.
