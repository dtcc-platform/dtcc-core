# Explicit native semantic-profile validation

Status: completed, 11 September 2026
Authority: DESIGN.md, model-contract.md and the user's approved validation experiment.

## Acceptance boundary

Make the versioned schema an explicit, editable contract for ordinary native DTCC
data: load a local profile once, validate a native city, and receive useful data
paths. Demonstrate a schema-only change and measure a complete real city tile.
Keep the API experimental and opt-in. No production dependency change, automatic
validation during model construction/exchange, remote schema discovery, generated
native classes, or full LinkML/CityGML conformance claim.

## Decisions

- Reuse canonical native admission before semantic projection; malformed native
  state raises the existing structural error. This limits the experiment to the
  current canonical object/geometry subset and includes that cost in measurements.
- Consolidate the existing optional LinkML/graph backend and buildings projection;
  the older process-separated examples must use the same implementation.
- Accept an explicitly selected self-contained local schema. Validate against
  that selection, independently of stored profile labels, without modifying data.
  No schema changes take effect on an already loaded profile instance.
- Preserve open attributes and reject undeclared relationships. Read semantic
  class/attribute/target constraints from YAML, with native containment and region
  hosts projected from their existing authorities. Numerical arrays stay native.
- Keep semantic diagnostics distinct from structural/configuration errors. Paths
  identify object IDs, attributes, typed native fields and geometry-region hosts.

## Checkpoints

- [x] Expose one reusable optional native validation entry point and consolidate
      the existing backend/projection without adding mandatory imports.
- [x] Exercise valid and invalid native examples, an edited local schema, existing
      profile coexistence and absence of LinkML in normal Core usage.
- [x] Measure setup, first and repeated validation, plus process memory on the
      unchanged complete 3DBAG canonical tile; investigate material overhead.
- [x] Run focused regressions, document scope, reproduction and measured adoption
      tradeoffs; self-review the final implementation.

## Verification

Use the existing DTCC environment for optional-dependency absence and native
regressions. For integration only, explicitly append the isolated LinkML tooling
site-packages to that interpreter; do not change the shared environment or embed
environment discovery in library code. Exercise the same public API in examples,
tests and the development benchmark. Re-run the original profile cases and current
buildings/openings/mixed-city evaluators. Performance runs use the recorded local
24,504,507-byte 3DBAG artifact and fresh processes, with no network or data repair.

## Completion evidence

`dtcc_core.model.profiles.SemanticProfile` now loads a local self-contained schema
and validates native Object roots, returning immutable reports with data paths.
It imports optional tooling on profile construction only. The existing backend
moved into `dtcc_core/model/_profile_backend.py`; the isolated evaluators use a
small import shim, and native building/opening/mixed-city examples share the new
projection. No production dependencies or normal model/exchange paths changed.

`native_profile_example.py` demonstrated a schema-only ChargingBench, required and
maximum-value errors, canonical file preservation and independent loaded versions.
The `validate_native.py` CLI passed valid (exit 0), missing required/range violations
(exit 1), and missing/malformed local schema (exit 2) checks. Missing LinkML is
reported with an actionable optional-tooling message; normal imports do not load it.

The final focused command used the explicitly composed development interpreter:

```bash
PYSTOW_HOME=/private/tmp/dtcc-profile-pystow ../venv/bin/python /private/tmp/dtcc-profile-python.py -m pytest tests/model/test_profiles.py tests/model/test_canonical_exchange.py tests/model/test_geometry_representations.py tests/io/test_mixed_city.py tests/io/test_building_openings.py tests/io/test_cityjson_buildings_profile.py -q
```

Result: **101 passed in 4.47 s**. The same selection in the normal native
environment passed **95**, skipping the **6** optional integration checks.
All original **18** comparison cases, the **7** buildings checks and the openings
and mixed-city evaluators passed with the shared backend. Their public native,
file/package, CityJSON and opening-mesh examples passed. Compilation and
`git diff --check` passed. No unrelated broader suite was needed for this optional
entry point; preceding milestones' broader evidence remains in their own plans.

Three sequential fresh measurement processes and a final verification process
validated all **15,554 semantic records** from **2,222 native objects** in the
unchanged **24,504,507-byte** 3DBAG artifact. Profile setup took **0.65–0.69 s**;
first validation **0.92–0.98 s**; twelve repeated calls **0.77–0.81 s**, including
native admission. Diagnostic native admission took **0.52–0.54 s** versus
**0.18–0.19 s** for schema plus graph checks. Whole-process RSS rose approximately
**65–68 MiB** above the loaded-model baseline, with peak increases of **84–88 MiB**.
Reusing existing compiled validators was sufficient; no new optimization layer
was warranted.

`docs/design/profile-validation.md` records the supported projection, interpretation
of valid, measured tradeoffs and reproducible development setup. Artifacts are in
`/private/tmp/dtcc-profile-validation-evidence/`; measured final reports are
`tile-2.json`, `tile-3.json`, `tile-4.json`, and `tile-final.json`. The schema files
and source artifact were not changed. Temporary extended schemas are local examples.

Self-review tightened diagnostic paths for actual containment/region-host/native
field authorities, rejected projection collisions and unsupported region reference
mappings, and normalized malformed YAML into a clear configuration failure.
There is no blocking issue within this acceptance boundary. Remaining limits are
explicit: canonical native type coverage; evaluated LinkML subset; projected IDs
rather than original-ID lexical constraints; no streaming/bounded-error API or
many-version soak measurement; and no supported production dependency distribution
decision. These costs support opt-in validation at selected checkpoints, not
automatic validation on every edit. All prior and unrelated work is preserved.

## Independent-agent handoff

Implement `.agent/plans/2026-09-11-profile-validation.md` through its checkpoints.
Keep the plan updated as material decisions or status change, preserve unrelated
work and all previous uncommitted milestones, and run the specified verification.
Follow DESIGN.md; keep LinkML optional and domain rules authoritative in the local
schema, with one shared projection and validation backend.
