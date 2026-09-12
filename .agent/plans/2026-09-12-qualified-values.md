# Qualified height measurements and classification codes

Status: complete, 12 September 2026.
Authority: user-approved next step; DESIGN.md, model contract and building crosswalk.

## Acceptance boundary

Introduce standard schema 0.6.0 with optional building/part height_measurements
and qualified classification codes. Use native nested dictionaries/lists and the
existing typed Protobuf values. No Python record classes, wire changes, new
dependency, network-dependent validator or implicit conversion/synchronization.

## Decisions

- Preserve scalar measured_height/estimated_height and building.height unchanged.
  height_measurements is a list of independent qualified vertical distances.
  Each requires value, unit, high_reference, low_reference and status. Status is
  measured or estimated, following CityGML 3.0 Height. Unit uses a small explicit
  length-unit vocabulary (m, cm, mm, km); no conversion or unit guessing.
- A QualifiedCode requires value and code_space, with optional label. Existing
  class, function, usage and roof_type slots accept either their plain strings
  or these records, keeping the value and its namespace together. Height references
  use the same forms. No universal code list or remote membership check is invented.
- An optional source string on a height record records author-supplied provenance;
  no claim of a complete provenance ontology, acquisition date or nil-reason model.
- Extend the existing LinkML backend for explicit identifier-free inlined records.
  Keep entity-ID graph checks separate; reject inline entity graphs or references
  inside value records until supported. Record classes cannot classify native Objects.
- Declared records are closed using LinkML's existing nested-record validation;
  misspelled keys fail. Ordinary undeclared Object attributes remain open and can
  carry extra supplier metadata. Required code_space prevents an unqualified
  dictionary or misspelled required namespace from passing as a QualifiedCode.
- String/record unions use LinkML's documented Any range constrained by any_of;
  the self-contained schema declares the linkml prefix and Any helper explicitly.
- CityJSON can carry these DTCC records as JSON attributes. No new official
  CityJSON/CityGML mapping is claimed and nested source keys are not renamed.
- Archive 0.5.0 unchanged; update the default and current documentation/examples.

## Checkpoints

- [x] Review CityGML Height, GML CodeType/MeasureType and LinkML inlining.
- [x] Implement the generic backend extension and schema, then exercise public I/O.
- [x] Add focused record/cardinality/diagnostic and unsupported-graph regressions.
- [x] Check native/package/strict CityJSON preservation, bypass and failed writes;
      run affected model/I/O checks and installed-wheel smoke.
- [x] Document scope, sources, Python use, measurements and remaining limitations.

## Verification

Use a small existing building CityJSON fixture augmented with explicitly authored
qualified values, including two reference pairs and units. No fabricated survey
facts attached to real data. Check ordinary load/save/package, value/namespace
preservation, unchanged scalar height, invalid required members/status/units/types,
unknown metadata and immutable old schema. A custom local schema must demonstrate
the record support is generic. Measure representative metadata validation with the
existing real tile to detect material regression; no large new benchmark harness.

## Completion and evidence

Implemented and example/public boundaries pass. Full model/I/O checks:
683 passed, 1 skipped in 36.11 s. A final generic-backend review also checks
identifier-bearing descendants of inline record ranges; the affected record and
profile checks then passed (11 tests, 5.23 s). No second validator or native record
classes were added.

The standalone authored example passed with exact native/package comparisons and
exact strict CityJSON comparison after aligning the aggregate ID. It retains the
centimetre-valued estimate, qualified namespaces, package context and a missing
scalar height. Artifacts and report: /private/tmp/dtcc-qualified-values.

Warm real-tile schema evaluation medians: 0.2888 s for archived 0.5.0 versus
0.2910 s for 0.6.0 on the same admitted source (five alternating evaluations).
Default native encode/decode medians: 1.4425/2.3029 s (three runs). A separate
synthetic 1,000-building/2,000-height-record model evaluates in 0.2849 s. These
measurements exclude first-use schema initialization and do not establish broad
scaling guarantees; they show no material regression on the exercised workload.
Evidence: /private/tmp/dtcc-qualified-performance.json and its temporary script.

Built the wheel without resolving dependencies, verified it contains the final
backend and 0.6.0 schema and excludes archived 0.5.0, then installed and exercised
it outside the checkout. Native/package/strict CityJSON round trips, version
selection, missing-unit diagnostic, file preservation and explicit bypass all
passed. Existing dependencies were reused; no browser, C++ consumer or other
platform certification is claimed. No new dependencies, wire changes, source
conversions or publication. Source review, documentation links and git diff --check
passed; unrelated work is preserved.

Verification logs: /private/tmp/dtcc-qualified-{boundaries,regression,final-focused,example,build,installed}.log.
Wheel: /private/tmp/dtcc-qualified-dist.
Installed smoke: /private/tmp/dtcc-qualified-installed-smoke.py.

Remaining scope: m/cm/mm/km only for qualified heights; code-space identifiers
are preserved but neither fetched nor membership-validated. No automatic scalar
height selection/conversion, GML adapter, qualified elevation/datum, uncertainty,
acquisition-time, nil-reason or full provenance model. The generic record support
is the foundation for later schema-only additions when those meanings are chosen.

Sources: OGC 20-010 Height tables 619/623; published construction/3.0
construction.xsd and GML basicTypes.xsd; LinkML inlining documentation. Cached OGC
copies from the previous crosswalk were inspected; web rendering cannot handle the
large conceptual HTML or XML schemas. See the resulting design note for links.

## Independent-agent handoff

Implement `.agent/plans/2026-09-12-qualified-values.md` through its checkpoints,
keep the plan updated as material decisions or status change, preserve unrelated
work, and run the specified verification. Keep scalar height access simple, use
one generic LinkML validation path and existing native metadata carriers, and do
not infer source measurements or introduce implicit unit/code conversion.
