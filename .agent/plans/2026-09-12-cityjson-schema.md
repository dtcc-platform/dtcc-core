# Standard schema at strict CityJSON boundaries

Status: complete, 12 September 2026.
Authority: user-approved next milestone; model contract and completed building
attributes/crosswalk plans.

## Acceptance boundary

Strict CityJSON load/save validates native structure and the selected standard
schema by default. validate_schema=False skips only semantics: malformed geometry,
references, unsupported mappings, duplicate keys and invalid schema identities
still fail. Reuse native admission and the existing semantic evaluator once each.

CityJSON has no DTCC schema declaration in the supported mapping. Strict import
selects the bundled default and records its ID/version on the returned City even
when evaluation is bypassed. Strict export selects the root's schema declaration
or default, without mutating that declaration. No CityJSON extension or schema
version change; CityJSON does not persist the DTCC declaration. No implicit schema
migration on export. Permissive behavior stays optional and is not certified;
explicit validate_schema=True requires strict=True. An omitted flag follows mode.

## Checkpoints

- [x] Connect public/direct JSON and ZIP I/O and dictionary conversion to default
  semantic admission, sharing canonical schema selection.
- [x] Verify valid exchange, invalid attributes/hosts, explicit bypass, structural
  failures despite bypass, unavailable/partial declarations and file preservation.
- [x] Update docs and any intentionally invalid semantic test inputs to request
  explicit bypass where their test is about a later boundary.
- [x] Run focused and full model/I/O tests and an installed-package smoke.

## Verification

Exercise io.load_city and City.save for .json and .json.zip, direct dictionary
conversion, native round trips and meaningful failure paths. Verify the boundary
invokes native and semantic admission once without serialization as a validator.
Run model/I/O tests after focused checks; no browser or live dataset run. Build and
smoke the installed wheel because public dispatch and shared schema selection change.
Preserve unrelated work. No new dependency, backend validator or wire format.

## Completion and evidence

Implemented default semantic evaluation for strict CityJSON file/dictionary
imports and dictionary/file exports, with JSON/ZIP public dispatch and the
City.save_cityjson convenience wrapper forwarding the option. The former duplicate
public JSON/ZIP reader delegates to the existing CityJSON loader. Canonical and
strict CityJSON writers share root schema selection via exchange._schema_for_model;
semantic evaluation remains the existing admitted-model path. No domain/schema
version change, Protobuf change, backend validator or dependency was introduced.

Strict imports record the bundled schema selection (currently 0.5.0), even with
bypass. Exports honor root selections without mutating them. CityJSON does not
encode that DTCC declaration; this limitation and mode-dependent option defaults
are explicit in docs/design/cityjson-schema-io.md. Updated current standard I/O,
model contract, inventory, building attributes and historical crosswalk notes.

Observed verification:

- Existing affected building/opening/mixed-city I/O checks: 56 passed.
- New boundary checks: 7 passed. They cover public JSON/ZIP, direct dictionaries,
  invalid attributes and semantic hosts, bypass, malformed indices/hosts,
  duplicate keys, unsupported mappings, declaration selection, invalid option
  values and preservation of existing files on failure.
- Full model/I/O suite: **679 passed, 1 skipped in 37.09 s**.
- A final public-entry audit found City.save_cityjson lacked option forwarding.
  Added it and reran the focused checks: **7 passed in 3.08 s**. The final wheel
  was rebuilt and checked to include that wrapper change.
- Instrumented real strict import and dictionary export: one exchange.validate
  call and one validate_admitted call each; bypassed export: one native call and
  zero semantic calls. No serialization or duplicate geometry admission added.
- Wheel build and installed-package smoke passed outside the checkout. Confirmed
  actual installed import location, JSON load failure/bypass, ZIP convenience
  save/load, schema selection, boolean-count rejection and failed-write preservation,
  including malformed geometry despite bypass. Existing dependencies were reused;
  this is not a fresh dependency-resolution or cross-platform certification.
- Updated documentation links, source self-review and git diff --check passed.
  No existing tests needed semantic bypass changes; no browser or live dataset run.

Evidence: /private/tmp/dtcc-cityjson-schema-initial.log,
/private/tmp/dtcc-cityjson-schema-focused.log,
/private/tmp/dtcc-cityjson-schema-tests.log,
/private/tmp/dtcc-cityjson-schema-build.log,
/private/tmp/dtcc-cityjson-schema-dist,
/private/tmp/dtcc-cityjson-schema-installed-smoke.py and
/private/tmp/dtcc-cityjson-schema-installed.log.

Remaining limits: strict mode is still explicitly selected; permissive CityJSON
has no semantic validation claim. CityJSON has no DTCC schema declaration mapping,
and its physical writer has not been changed into a crash-atomic writer. Unsupported
format mappings remain unsupported. Other external adapters and artifact-only
packages need their own boundary adoption; qualified records and additional urban
themes remain separate schema work.

## Independent-agent handoff

Implement `.agent/plans/2026-09-12-cityjson-schema.md` through its checkpoints,
keep the plan updated as material decisions or status change, preserve unrelated
work, and run the specified verification. Keep format admission and semantic
admission distinct and once per operation; preserve the explicit bypass without
weakening native integrity or unsupported-content failures.
