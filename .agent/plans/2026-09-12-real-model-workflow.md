# Real-data model workflow with the standard schema

Status: complete, 12 September 2026.
Authority: user-approved real-data integration milestone; DESIGN.md and the model contract.

## Acceptance boundary

Use the unchanged, checksum-pinned cached 3DBAG tile. Strictly import under the
current standard schema, edit building metadata, round-trip native data, mesh
one explicitly selected building part, and export/reload a canonical package.
Preserve source identity, representations, shell topology, attributes and regions;
retain derived mesh fields and package Dataset Context. Measure public operations
with default validation enabled. No source repair, implicit height inference,
new dependency, wire/schema revision, or claim of whole-tile geometric validity.

## Decisions

- Reuse the existing representations example and meshing implementation.
- Close the observed missing Solid.mesh entry point with boundary triangulation
  through the existing surface/region mapping. Include every shell's surfaces;
  keep the Solid as the authority for shell partition. The derived Mesh is a
  surface mesh, with no volume cells or shell partition. Reject unsupported
  remapping and field transfer, including for solids without semantic regions.
- Check strict CityJSON exchange before adding the derived Mesh, because that
  external adapter has no supported mapping for the added representation/field.
- Standalone native files retain model data; canonical packages additionally
  preserve Dataset Context. Report this existing distinction explicitly.
- The selected 271-polygon part exposed a normal-calculation bug: the first
  almost-collinear vertex triple tilted a roof's projection, losing 841 m².
  Fix Surface.calculate_normal at its existing authority using the full ring's
  translated area vector. Verify winding, large offsets and degenerate failure;
  rerun relevant surface/builder tests and real per-region area checks.

## Checkpoints

- [x] Exercise the actual strict import and identify workflow gaps.
- [x] Add the bounded Solid meshing entry point and focused invariant coverage.
- [x] Extend and run the existing complete-tile example with edit, mesh, field,
      default-validation failures, round trips and time/memory evidence.
- [x] Document results, review changes and run affected checks plus installed smoke.

## Verification

Test all-shell triangulation, region mapping, transform and input preservation,
and rejection of invalid rings or unmapped fields. Run existing semantic meshing
and geometry representation tests, then relevant builder/I/O tests if changes
warrant them. Run the checksum-pinned example sequentially, with public-operation
timings and clearly labelled process RSS. Build and smoke the wheel outside the
checkout because the new builder method depends on module registration.

## Completion and evidence

The full checksum-pinned workflow passed with schema 0.5.0 and wire 6. It retains
all 4,443 source representations and adds one Mesh (904 triangles, four semantic
regions and a face-associated area field). Canonical comparisons are exact;
package Dataset Context matches, including the 1,111 extent discrepancies.
The formerly failing first source part also meshes after the normal fix (3,186
triangles). Four collapsed source rings remain in the unchanged model; no full-tile
meshing or geometric validity claim is made. A boolean storey-count edit fails
default save validation and leaves the existing native file unchanged.

23 focused checks passed. The full model/I/O/builder/reprojection run produced
1,226 passed, 5 skipped and one stale test failure in 38.76 s. The failure still
looked up an enum directly in the superseded geometry dictionary layout. Updated
that assertion to the existing public get_geometry selector; all three flattening
checks then passed in 0.53 s. No runtime change was needed for that failure, and
the already-passing broad checks were not repeated for the assertion-only edit.

Built the wheel without dependency resolution, installed it into the existing
isolated environment, and ran the smoke outside the checkout. Verified actual
installed import location, schema 0.5.0/wire 6, full real-data native/package load,
public Solid.mesh registration, identical remeshed geometry, source preservation,
the 904-value float64 face field and package context. Only the installed
dtcc_mesher backend/macOS environment was exercised; this is not cross-platform,
Triangle-backend, dependency-resolution or volume-meshing certification.

Source self-review and git diff --check passed. No new dependency, schema/wire
revision, source repair or external publication. Existing unrelated work remains.

Measurements and source evidence: docs/design/real-model-workflow.md,
/private/tmp/dtcc-real-model-workflow/report.json and
/private/tmp/dtcc-real-model-memory/{import,decode}-0.json. Public native save/load:
1.456/2.386 s; selected triangulation: 44 ms; fresh import/decode peaks:
456.6/468.5 MiB including dependencies. Schema/wire/dependencies are unchanged.

Verification logs: /private/tmp/dtcc-real-model-{focused,regression,flatten,build,installed}.log.
Installed smoke: /private/tmp/dtcc-real-model-installed-smoke.py.
Wheel: /private/tmp/dtcc-real-model-dist.

## Independent-agent handoff

Implement `.agent/plans/2026-09-12-real-model-workflow.md` through its checkpoints,
keep the plan updated as material decisions or status change, preserve unrelated
work, and run the specified verification. Keep this to the real-data workflow
and the observed Solid boundary meshing gap; do not repair source geometry or
expand schema scope to make the example succeed.
