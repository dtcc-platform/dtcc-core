# Building openings and boundary relationships

Status: completed, 10 September 2026
Authority: DESIGN.md, model-contract.md and the user's approved next step after
the completed live-memory milestone.

## Acceptance boundary

Carry window/door surface regions and their wall/roof relationships through the
ordinary strict CityJSON, canonical file/package and native meshing workflows.
Keep domain rules in a new self-contained buildings profile, retain 0.1.0, and
keep native geometry generic. No physical window/door object hierarchy, geometric
repair, full CityGML conformance, new runtime dependency or production validator.

## Decisions

- One optional region-local parent index; derive inverse children. Indices match
  the existing ordered region store and need no generated IDs or object graph.
  Reordering regions requires remapping parent indices; merge offsets them.
- Canonical v4 adds this fact; v1–v3 remain readable. Freeze a v3 fixture before
  updating the writer. Reject new facts under old version declarations.
- CityJSON parent-only or children-only relations normalize to one parent.
  Reject malformed, dangling, cyclic or conflicting references. Export both
  directions deterministically; redundant key presence/child order is not a fact.
- Profile 0.2.0 retains region ownership and adds an optional host reference for
  Window/Door surface regions. The vocabulary means CityGML WindowSurface and
  DoorSurface, not the separate physical filling elements. Missing host remains
  missing. Wall/roof target rules are schema facts, not native class checks.

## Checkpoints

- [x] Verify CityJSON 2.0.2 semantics and CityGML 3.0 filling-surface distinction.
- [x] Implement generic region relationships, wire evolution and strict mapping.
- [x] Add profile 0.2.0 and one facade fixture with window and door; exercise
      ordinary load/access/save/package/mesh and portable profile validation.
- [x] Verify malformed relationships, older wires/profiles, merge remapping and
      affected regressions; document scope, evidence and remaining limitations.

## Verification

Focused openings tests first; existing model, IO, datasets, reproject and semantic
meshing tests after changes. Validate the small source/export against the locally
recorded official CityJSON 2.0.2 schema in the isolated tooling environment. Run
existing profile examples/evaluators and a representative complete-tile codec
comparison to check the ordinary no-opening path and avoid a performance regression.
No network/download requirement for reproducing checked-in synthetic examples.

## Completion evidence

`docs/design/building-openings.md` records the implemented contract, primary
sources, exact scope and reproduction commands. The native change is one optional
integer; domain types and host target ranges live in profile 0.2.0. Region merging
remaps parent indices, and existing semantic meshing preserves them automatically.
The optional profile validator now rejects undeclared relationships even when
attributes are open, so a supplied host cannot silently escape graph checking.

Public file/package/CityJSON/mesh example passes, including independently editable
metadata and front-wall/window/door areas 50/4/6 m². Source and exported fixtures
pass the recorded official CityJSON 2.0.2 schema. The new portable profile check
and previous 7-case/18-case evaluators pass; old 0.1.0 remains unchanged.

Affected regression run: 1,315 passed, 1 skipped, 53 deselected (39.48 s). The
complete 24,504,507-byte 3DBAG tile rewrites identically except for the v4 envelope
version. A same-process timing control against reconstructed pre-parent codec
paths initially found 2.7% overhead; reading the envelope version once reduced
the final median difference to 0.4% (2.098 -> 2.106 s). Both variants share current
native dataclasses; this does not measure old-runtime memory. Raw evidence is in
`/private/tmp/dtcc-openings-evidence/`.

Final focused model/version/package/openings checks and `git diff --check` passed;
changed Python modules compiled. Self-review found no blocking issue. Index-based
relationships require explicit remapping after manual region-list reordering;
this editing obligation is documented. Physical filling objects, broader geometry
operations, production LinkML integration and downstream consumers remain outside
this completed slice. Prior uncommitted work and unrelated files are preserved.

## Independent-agent handoff

Implement `.agent/plans/2026-09-10-building-openings.md` through its checkpoints.
Keep the plan updated as material decisions/status change, preserve unrelated work
and all previous uncommitted milestones, and run its specified verification.
Follow DESIGN.md and the native model contract; keep one geometry/relationship
authority and keep LinkML outside production dependencies.
