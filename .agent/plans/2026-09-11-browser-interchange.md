# Prove browser interchange with the public Protobuf contract

Status: complete, 11 September 2026.
Authority: DESIGN.md, docs/design/model-contract.md, and the user's request to
continue the overhaul checkpoints after the unified-format milestone.

## Acceptance boundary

One useful slice: Python saves a mixed model using ordinary .dtcc I/O, a real
browser reads it with protobuf.js and the repository's dtcc.proto, edits it and
encodes it, and ordinary Python loading verifies the returned model. Demonstrate
building footprint access, identities, float64 coordinates, uint64 arrays,
arbitrary integer attributes, fields, opening relationships, unfamiliar semantic
types and rejection of invalid browser output. Preserve all unrelated changes.

This is a development interoperability example, not a second model runtime or a
production JavaScript SDK. Keep its npm dependency local and development-only.
Do not duplicate the LinkML validator in JavaScript, add another wire schema,
change sibling consumers, publish artifacts or restore legacy compatibility.
Full CityGML crosswalk, wrapper decisions, consumer adoption and performance
remain subsequent vertical slices of the wider roadmap.

## Checkpoints

- [x] Build a small reproducible local browser example and Python fixture/checker.
- [x] Run Python → real browser → Python, including one semantic and one numerical
  failure; preserve unknown wire fields through the relay so that Python rejects
  them rather than accepting silently stripped data. Inspect the actual browser
  outcome and record the evidence.
- [x] Document the language-boundary lessons and update the inventory checkpoint.

## Verification

Use public save/load entry points with default schema evaluation. Compare the
decoded Protobuf messages for the untouched round trip and an independently
constructed expected native edit, rather than requiring byte-identical encoding.
Check native footprint access after the browser edit. Invalid messages must fail
through the ordinary Python reader and leave an existing receiver unchanged.
Exercise browser array access and bounded malformed-array failure. Run syntax
checks and the existing focused unified-format tests; broader native tests only
if production code changes. No large-data throughput claim from this small fixture.

## Completion evidence and remaining scope

Implemented `sandbox/protobuf_interop/` with a pinned development-only protobuf.js
dependency, local Python fixture/server and browser example. The real Chromium 152
browser visibly reported PASS for the untouched round trip, independently checked
native edit and three rejected inputs (negative measuredHeight, malformed array,
unknown wire field). All rejected inputs preserved the existing Python receiver.
The browser also checked exact float64/uint64/arbitrary integer values, footprint
access, fields, UTF-8 identities, named references, typed defaults and opening parent
presence. The generic unfamiliar bench URI passed the default semantic boundary.

Observed the browser library's default unknown-field discard behavior in its source
and documentation. Configured reader.discardUnknown=false and verified relay to
Python rejection rather than silent stripping. No production codec/schema/API
change or JavaScript semantic validator was required. README records reproducible
commands, results and integration implications; native inventory now links it.

Verification: real browser workflow passed; existing unified-format tests **5
passed in 3.18 s**; JavaScript/Python syntax checks, offline npm ci from the lockfile
and git diff --check passed. Evidence files are in
`/private/tmp/dtcc-browser-interchange/` (Python source/returned binary files and
report.json). Browser assertions and Chromium version were observed directly in
the page. The synthetic return is 4,955 bytes, edited return 4,956 bytes.

This completes the browser interchange milestone, not the entire model overhaul.
Next is the urban property/relationship/units crosswalk. Production TypeScript
consumer and C++ integration, the eight unsupported wrapper decisions, broader
package/adapter adoption and representative performance work remain in the
inventory. No external publication, sibling edits or new production dependency.

## Independent-agent handoff

Implement `.agent/plans/2026-09-11-browser-interchange.md` through its checkpoints,
keep the plan updated as material decisions or status change, preserve unrelated
work, and run its specified verification. Use the single public dtcc.proto and
real browser execution; distinguish Protobuf decoding from semantic validation.
