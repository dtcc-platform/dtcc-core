# Model hardening: issue 85 closeout

Status: complete locally, 12 September 2026. This assessment applies
to [issue 85](https://github.com/dtcc-platform/dtcc-core/issues/85), its review
comments and the user's subsequent explicit format and schema decisions.
It is the issue acceptance checklist, not a promise to implement all CityGML.

## Original acceptance criteria

| Criterion | Implementation and proof |
|---|---|
| A model serialization/invariant matrix exists | `tests/model/test_serialization_matrix.py` inventories all 36 exported Model subclasses: 27 supported concrete roots, eight deliberately unsupported wrappers/results, and abstract Geometry. `model-inventory.md` states the supported native facts and adapter boundaries. |
| Supported public types, including nested city objects, round-trip | Minimal Protobuf bytes/JSON and representative nested geometry/object tests cover the admitted types. The unified codec preserves typed arrays, transforms/CRS, fields and associations, markers/normals, region membership, IDs and concrete collection identity. Real-data and package workflows supplement the matrix. |
| Unsupported paths fail explicitly | Unsupported wrappers, unknown concrete Python subclasses, unrepresentable native state and unsupported strict CityJSON mappings fail rather than downcast or disappear. A schema-unknown semantic URI on a supported generic Object is deliberately accepted as generic data; this does not invent semantic validation or lose its original URI. |
| Known quick-review gaps are fixed or assigned concrete follow-up scope | Tree, Landuse, nested collections, markers/normals, child geometry metadata, key normalization, attribute and array admission, precision and self-description are implemented. Cache behavior is explicit in `model-spatial-contract.md`; remaining global-frame and aggregate-coverage work is scoped below. |
| Existing CityJSON and Protobuf I/O tests continue to pass | Final affected regression and real-browser results are recorded below. Strict CityJSON is a documented subset; legacy permissive import is not a lossless exchange claim. |

The issue originally suggested compatibility where possible. The user later
explicitly retired legacy formats. `dtcc.proto` now defines the sole version 6
`DTCC.ModelFile` message, normally named `.dtcc`; there is no old `.pb/.pb2`
model decoder or automatic migration. Protobuf field/value types and schema
identity/version are explicit. A filename is not the format definition.

## Scope achieved beyond the original issue

The standard [schema 0.9.0](../../dtcc_core/schemas/dtcc.yaml) is
self-contained and evaluated by default at native/canonical-package and strict
CityJSON boundaries. `validate_schema=False` skips semantic evaluation, not
numerical, graph, file or package integrity checks. Direct Python access remains
ordinary arrays, dictionaries and native geometry accessors.

The [exterior-city profile](exterior-city-profile.md) covers 12 external feature
types and their deliberately restricted geometry/metadata mappings. Native DEMs
retain explicit elevation interpretation. Qualified heights/elevations/codes,
multiple geometry representations, solids, semantic regions and schema-only
subtypes use generic carriers. These achievements do not establish full CityGML
conformance, a CityGML XML adapter or a production JavaScript/C++ model SDK.

## Closeout fixes

The browser example no longer embeds a stale schema version or URI suffix.
Its local fixture expectations come from the Python-selected schema; the browser
still independently checks the actual decoded values, representation identity,
float64/uint64/arbitrary integer precision and edits. The public dtcc.proto is the
only wire definition. No JavaScript LinkML validator is introduced.

Spatial review found that refreshing Object bounds could overwrite the physical
domain of Grid/VolumeGrid with their cell counts. The domain is intrinsic state,
so refresh now retains it, including explicitly degenerate domains; lazy defaults
remain available. A focused workflow verifies mutation, refresh, source transform
preservation and native reconstruction. General world-frame composition remains
outside this fix.

Reprojection review also identified source-CRS mutation and silent metadata loss.
The closeout hardens the existing bounded reprojection path and rejects state
whose coordinate-dependent meaning it cannot preserve. It does not claim a
complete scene hierarchy, field-component or vertical-datum transformation engine.

Horizontal reprojection supports PointCloud, Surface (including holes), Mesh and
VolumeMesh (including connectivity/markers), MultiSurface and one non-nested object
whose spatial state is entirely in its supported representations. A CRS change
returns an independent copy with the target CRS and no stale derived cache; Z is
retained. A same-CRS no-op returns the original object. An actual two-axis source
CRS must be supplied or declared; conflicts fail unless the PointCloud override
explicitly selects the supplied source. Following #112, fields retain their values,
associations, units and component conventions on the copied geometry. Vector
fields warn on a CRS change because their components are not rotated or rescaled.
See [field-preserving reprojection](../reprojection.md) for usage and limitations.
Regions, stored normals, local affines, DatasetContext, nested objects, intrinsic
Tree/graph coordinates and unmapped geometry fail before a misleading result is
returned. Full vector-field,
vertical-datum and provenance transformation requires a separately scoped task.

## Separately scoped follow-up drafts

These are local, copy-ready issue drafts, not published GitHub issue numbers.
They prevent the original hardening issue from becoming an unlimited model backlog.

### F1 — Define world-frame bounds and hierarchy transform semantics

Problem: Object bounds currently aggregate raw local coordinate envelopes and do
not compose parent/child affine transforms. A bounds result cannot automatically
be treated as a world-space footprint or canonical package extent.

Acceptance: define frame inheritance and CRS conflict rules; specify matrix
composition order and transform-versus-reprojection semantics; implement one
nested transformed-object workflow with a correctly framed envelope; reject
ambiguous frames; populate package bounds only when their frame is explicit.
Keep public NumPy mutation and explicit refresh inexpensive. No implicit observer
framework is required.

### F2 — Define aggregate coverage for empty and intrinsic-only geometry

Problem: empty coordinate geometries can contribute an origin envelope, and
intrinsic spatial state such as Tree.position or RoadNetwork graph arrays is not
uniformly included in Object.calculate_bounds.

Acceptance: define when an object has no extent, which intrinsic representations
contribute, and how mixed children combine; exercise an empty child, native tree,
and graph-only road network without fabricated origin bounds. Preserve declared
Grid/VolumeGrid domains and the explicit frame boundary from F1.

### F3 — Adopt the model contract at production consumer boundaries

Problem: a successful development browser relay is not production Twin/TypeScript
or C++ consumer adoption. Canonical manifest-v3 packages are opt-in; default
artifact-only packages and remote consumers have separate contracts.

Acceptance: choose one actual consumer, decode the public dtcc.proto, compare
coordinates and all required semantic/numerical facts, reject unsupported
versions/state without silent field stripping, and exercise canonical package
context/integrity through that consumer. Add C++ fixtures at its actual boundary
when needed. Address deployed producer revision identification explicitly;
`dtcc_core.__version__` is not currently a provided API. Record installed package
version and VCS provenance using supported package metadata, without silently
recording an unknown placeholder as verified provenance.

### F4 — Define durable simulation values and wrapper scope from real workflows

Problem: Field supports numeric arrays and explicit associations but not general
named component/time/scenario/ensemble axes. Eight result/container types remain
intentionally unsupported; serializing only their geometry would discard context.

Acceptance: select an actual simulation/result workflow; identify durable facts
versus convenience wrappers; preserve required axis/component meaning and source
identity or formally leave the wrapper nonserializable. Do not add a universal
Python-object serializer or implement every wrapper solely to increase coverage.

Further CityGML themes, appearance/instancing and supplier mappings require their
own concrete use cases. They are outside issue 85's original acceptance boundary.

## Verification and integration evidence

Final source regression: **1,811 passed, 5 skipped, 53 live tests deselected** in
52.84s. Command: `python -m pytest tests/model tests/io tests/builder tests/datasets
tests/reproject tests/test_top_level_api.py -q`. Log:
`/private/tmp/dtcc-issue85-final-regression.log`.

Real Chromium 152 browser: all round-trip, edit, semantic/array/unknown-field
rejection cases passed under schema 0.9.0; receiver state survived failures.
Browser-source syntax and Python compilation checks passed. Artifacts:
`/private/tmp/dtcc-issue85-browser/`.

The final wheel built successfully and was installed with existing dependencies
in `/private/tmp/dtcc-standard-schema-env/`. Exterior-city and DEM examples passed
outside the source checkout; the spatial-domain and reprojection workflows also
passed against that installed package. Its changed Grid/Object/reprojection files
and schema were byte-compared with source. Wheel:
`/private/tmp/dtcc-issue85-dist/`; artifacts:
`/private/tmp/dtcc-issue85-wheel-artifacts/`.

Three subagents reviewed current documentation, spatial invariants and selected
exchange/profile/package/reprojection paths. Identified blockers were fixed;
`git diff --check` passed. This is focused integration evidence, not an exhaustive
proof of every algorithm or downstream consumer. Earlier timing and real 3DBAG
evidence remain in the completed exterior-city milestone.

The reviewed implementation is integrated locally on `develop` in commit
`3bbafbe3ba90d38ea89dbf69d13db92df2156dc5`. The unrelated untracked history image was
preserved. The completion record is a subsequent documentation-only commit.
No release, remote push, GitHub comment or issue creation/closure was performed;
the issue is ready for publication of the changes and the closure summary below.

## Copy-ready GitHub closure summary

The model hardening acceptance criteria are implemented and verified. A public-type
matrix classifies all 36 exported Model subclasses, with 27 supported concrete roots,
eight explicit unsupported wrappers/results and abstract Geometry. Supported
Protobuf/native round trips preserve nested identity, coordinates/typed arrays,
transforms/CRS, fields/associations and declared metadata; unsupported paths fail
clearly. The user-approved replacement is one self-describing `DTCC.ModelFile`
message in `dtcc.proto`, wire 6, normally named `.dtcc`.

Schema 0.9.0 supplies a self-contained standard semantic contract, enabled by
default at native/canonical-package and strict CityJSON boundaries. Python access
remains generic and direct. Bounds/cache limits and unsupported transformations
are explicit. Closeout review fixed intrinsic grid-domain overwrite and public
reprojection source mutation/metadata loss. Final tests: 1,811 passed, 5 skipped,
53 live tests deselected; actual Python/browser/Python and installed-wheel examples
also passed.

Global-frame/aggregate extent semantics, production consumer/package adoption and
richer simulation/wrapper contracts are separately scoped in follow-up drafts
F1–F4 in `docs/design/model-issue-85-closeout.md`. Full CityGML support is outside
this issue's acceptance boundary. Link selected follow-up issues when publishing
this summary; do not imply that these local drafts are already GitHub issues.
