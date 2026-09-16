# Adopt the revised model in downstream DTCC packages

Status: complete locally (13 September 2026). Authority: the user requested principled downstream fixes
following the local audit of seven sibling repositories. No PLAN_TEMPLATE.md was
found in this workspace. Core DESIGN.md and the model/package contracts control
behavior; this milestone does not add a viewer, solver, or CityGML capabilities.

## Acceptance boundary

Sim produces and delivers `.dtcc` results with explicit field associations.
Core's publication client and Upload transport canonical manifest-v3 packages
without losing model bytes or context, retaining path/hash/size and authorization
checks. Atlas and Tangible Twin recognize the new format and consume supported
display artifacts in canonical packages, keeping native models distinct from
display derivatives. Existing artifact-oriented package workflows remain valid.
Twin has no implementation to migrate; meshing libraries have no model-format
dependency. No legacy Protobuf reader is reintroduced. No production dependency
is added, no remote publication is performed, and no production deployment is
claimed. Numerical solver verification is limited to available runtimes.

## Checkpoints

- [x] Update Sim formats, field construction, result delivery, demos, and contract checks.
- [x] Enable canonical package publication and Upload storage with a real Core
      package, exact download/reload proof, and corrupt-input rejection.
- [x] Update Atlas format metadata/package consumption and Tangible Twin catalogue,
      importer, and browser package handling using real canonical fixtures.
- [x] Coordinate dependency revisions, document consumer responsibilities, review
      all changes, and run focused cross-repository regression checks.

## Verification

Use current Core source explicitly. Exercise synthetic traffic and volume fields
through public dataset/result boundaries without network providers. Export a real
canonical package, submit it through Upload's local test client, download it and
reload through Core. Verify semantic/numerical facts and context, plus rejection
of altered artifacts. Run existing Upload safety tests, Atlas package/download
tests, Tangible Twin generator/importer/browser tests, and affected Core package
tests. Browser display consumes supported supplemental media, not the native
model itself. Record exact passed/failed/skipped checks and remaining limits here.

## Decisions and completion evidence

Preserve unrelated untracked history images. The initial audit found
that current Sim/Upload contract checks pass while missing native result delivery
and canonical package adoption; those authoritative checks are now extended.

## Independent-agent handoff

Implement `/Users/logg/scratch/dtcc/dtcc-core/.agent/plans/2026-09-12-model-downstream-adoption.md`
through its checkpoints across the sibling repositories, keep this plan updated
as material decisions or status change, preserve unrelated work, and run the
specified verification. Do not expand model scope or publish remote changes.


### Implementation evidence (13 September)

- Sim: 45 passed, 5 skipped (optional FEniCSx environment absent). Native result
  delivery preserves mesh arrays, field association/unit and values; ambiguous
  association fails before producing a result file. Traffic uses the new format.
- Core: 583 affected dataset/canonical checks passed, 53 live checks deselected;
  after the final publication tests and idempotency change, 91 focused checks
  passed. The model publication workflow and rejection of extra ZIP members are
  covered in Core itself as well as the cross-repository Upload tests.
- Upload: the final full 165-test suite passed, including the catalogue summary
  selecting the canonical model rather than its alphabetically first sidecar. Real Core packages, both v2 and v3, traverse the
  authenticated API with exact byte preservation. Real Core publication client
  paths cover directories and ZIP packages, download/reload and corrupted input.
- Atlas: 226 tests passed after integrating listing, format-choice schema, native
  and preview downloads and jobs. The real-Core subprocess verifies native
  corruption rejection; normal Atlas unit tests otherwise mock Core.
- Tangible: 21 Python generator tests and the actual synthetic smoke catalogue
  command passed. It produced one .dtcc package without network simulation data
  or publication. 276 browser/CLI tests passed; the additional online canonical
  preview case and affected control-panel tests passed (74 focused checks).
  The final complete browser/CLI run passed 277 tests. TypeScript checking and
  production Vite build passed; the existing large-bundle advisory remains.
- Real canonical fixture: /private/tmp/dtcc-downstream-contract-20260912.
  Smoke workflow: /private/tmp/dtcc-downstream-smoke-catalog.
- Bounded work: archive member/size limits and bounded decompression protect the
  newly admitted native packages. Listings avoid model-array decoding; native
  artifact consumption uses Core. Online projector viewing fetches only the
  manifest and selected preview, checking its digest; it does not claim to have
  decoded the native model. Native-only packages have no projector preview.
- No new production dependencies or production browser model SDK. Existing v2
  artifact workflows remain supported. No service deployment or remote writes.


### Local integration

Core runtime/publication revision: `5ca2ca410f24763591dc62c7b61f870cef13f717`.
Sim (`e24a1f2`) and Atlas (`b201ff8`) both pin that exact revision; Atlas's UV
source declaration agrees with its project dependency. Upload is `26a53e3` and
Tangible Twin is `4f49558`. Core's final completion record is documentation-only.
Upload and Atlas had no configured Git identity, so their commits used the
identity already configured in Core for those commands only; no global or local
Git configuration was changed.

All changed-repository diff checks passed. Generated untracked Python caches were
removed from Tangible Twin; unrelated history images were retained. These are
local commits only. Publish the Core revision before installing the downstream
pins from GitHub, then roll out the consumers together. No remote push, service
publication, deployment, data deletion, or legacy-file migration occurred.
