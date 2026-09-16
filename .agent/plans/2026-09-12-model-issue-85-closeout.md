# Close out model hardening issue 85

Status: complete locally (12 September 2026). Authority: user agreed to the final issue-closeout milestone after
reviewing the original issue's acceptance criteria. Existing DESIGN.md and the
model contract control behavior; this plan does not expand the exterior profile.

## Acceptance boundary

Reconcile the original issue's five acceptance criteria with implemented behavior,
refresh and run Python/browser/Python interchange under schema0.9.0, and make
bounds/cache/local-coordinate limits explicit with a focused workflow check.
Known remaining concerns must have concrete follow-up descriptions and acceptance
boundaries; supporting every wrapper or all CityGML is not required for issue85.
Review the accumulated implementation for integration, preserve unrelated work,
and prepare a local integration commit and copy-ready issue closure/follow-up text.
Do not publish, push, or post GitHub messages as part of this local milestone.

## Checkpoints

- [x] Reconcile current model contract and inventory, distinguishing historical evidence.
- [x] Refresh browser fixture expectations without duplicating schema authority;
      run a real browser and preserve numerical/semantic/unknown-field failures.
- [x] Record bounds and transform behavior, exercise mutation/recalculation and
      exchange, and define separately scoped follow-ups for remaining contracts.
- [x] Review accumulated implementation, run affected checks and normal examples,
      resolve blockers, and integrate task changes locally while preserving unrelated work.
- [x] Record exact verification and completion evidence; prepare issue85 closure
      summary and actionable follow-up drafts with no implied full CityGML claim.

## Verification

Browser round trip checks actual coordinates, numeric precision, metadata, nested
identity, browser edits and rejection before receiver replacement. Bounds proof
checks public mutation, explicit refresh, and persisted-state reconstruction;
no observer framework or speculative global-transform algorithm. Reuse the model,
I/O and affected builder/package regressions. Review and independently inspect
changes before committing. Do not claim unrun C++/production-consumer certification.

## Decisions and evidence

The authoritative closeout summary is docs/design/model-issue-85-closeout.md.
Three subagents reviewed documentation, spatial integrity and integration. The
browser visibly passed under schema0.9.0/Chromium152 after moving evolving fixture
expectations into Python-generated fixture.json. Review found and resolved
Grid/VolumeGrid explicit-domain overwrite (includingzeroareas); public reprojection
source mutation/metadata loss now has a bounded copying/fail-closed contract.
These are corrections to current integrity requirements, not new city features.
No general transform hierarchy, observer framework or metadata transformation
engine was introduced. Four follow-up drafts bound the remaining broader work.
Existing milestone plans remain completed historical evidence.

## Completion evidence

- Implementation integrated on `develop` as `3bbafbe3ba90d38ea89dbf69d13db92df2156dc5`; this final
  plan/closeout record is documentation-only. All reviewed task files were staged
  explicitly; unrelated `git-loc-history-2026-08-01_21-06-37.png` was excluded.
- Current acceptance checklist, final contracts and copy-ready closure/follow-up
  text: `docs/design/model-issue-85-closeout.md`. Four bounded follow-up drafts
  cover global-frame semantics, empty/intrinsic aggregate extents, production
  consumers/package adoption and simulation/wrapper meaning. They are local
  drafts, not published issue numbers.
- Final affected suite: **1,811 passed, 5 skipped, 53 live tests deselected**,
  52.84s. Command and log are recorded in the closeout document. Earlier two
  obsolete normal-preservation expectations were replaced with meaningful
  fail-closed checks; no final test failures remain.
- Actual Chromium152 Python/browser/Python under schema0.9.0 passed unchanged and
  edited values plus semantic, array and unknown-wire failures. JavaScript syntax,
  Python compilation and staged diff checks passed. Browser fixture expectations
  now come from the selected Python schema, avoiding duplicated version/URI constants.
- Grid/VolumeGrid retain explicit and degenerate domains through recursive refresh,
  resolution edits and persistence. Public reprojection preserves supported bare
  source state, fails ambiguous/unmapped cases, and cannot silently change integer
  Z values. Focused spatial/reprojection proofs and independent review passed.
- Built `/private/tmp/dtcc-issue85-dist/` wheel and installed without new dependencies.
  Exterior/DEM examples and spatial/reprojection checks passed from outside the
  checkout; changed installed production files and schema match source bytes.
- No push, release, external message or GitHub issue mutation. The local test
  server was stopped after browser verification. Publication/issue closure is
  a separate authorized action, not an unreported part of this milestone.

## Independent-agent handoff

Implement `.agent/plans/2026-09-12-model-issue-85-closeout.md` through its checkpoints,
keep the plan updated as material decisions or status change, preserve unrelated
work, and run the specified verification. Finish the finite acceptance and local
integration milestone; do not expand CityGML scope or publish external messages.
