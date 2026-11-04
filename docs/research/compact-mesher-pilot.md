> Historical experiment. Its options and recovery behavior are superseded by the
> [current LOD2 pipeline](../lod2_pipeline.md). Saved results are historical evidence.

# Compact catalogue surface meshing

LOD2 calls using the resolved `dtcc_mesher` backend construct minimal planar
roof/wall domains and apply bounded improvement by default. Both public entry
points use the same path:

```python
from dtcc_core import builder
from dtcc_core.builder import CompactSurfaceOptions

mesh = builder.build_city_surface_mesh(city, lod="lod2",
                                       mesher="dtcc_mesher")
mesh = city.build_surface_mesh(lod="lod2", mesher="dtcc_mesher")

# Disable all optional point insertion; legal zero-point flips remain enabled.
zero = CompactSurfaceOptions(budgeted=False)
mesh = city.build_surface_mesh(lod="lod2", mesher="dtcc_mesher",
                               compact_surface=zero)

# Initial bounded improvement only (Stage A).
baseline = CompactSurfaceOptions(extension_source_vertices=0,
                                  extension_group_vertices=0)
# Fixed standard control (also the default).
standard = CompactSurfaceOptions(adaptive_extension=False)
# Preserve the standard prefix, then assess guarded extra edge points.
adaptive = CompactSurfaceOptions(adaptive_extension=True)
# Explicit doubled extension control; maxima, never insertion targets.
doubled = CompactSurfaceOptions(extension_source_vertices=8,
                                 extension_group_vertices=16,
                                 adaptive_extension=False)
```

`compact_surface=None` selects the normal compact path for catalogue/dtcc_mesher
calls. RANSAC is selectable with `lod2_method="ransac"` and uses the same
compact mesher. Default `lod="lod1"` and ordinary LOD1 behavior are unchanged.
Both LOD2 methods support isolated and joint surface layouts. Isolated compact
skins use the same planar faces, shared-edge registry and bounded improvement,
with an explicit base height instead of a reference terrain attachment. Their
only permitted open edges lie at that base. Terrain is meshed separately and no
vertices are welded between buildings or terrain. The `volume_boundary=0` guard
remains for isolated and separated outputs. Joint merged LOD2 shells can pass
through PLC validation into the TetGen volume pipeline.
The `City` method also accepts keyword-only `report_mesh_quality` for parity with
the builder entry point.

## Options, backends and dependencies

An explicit `CompactSurfaceOptions` requires `lod="lod2"` and excludes an explicit
`surface_refinement`. An explicit `surface_refinement` with `compact_surface=None`
selects that existing experiment instead of the automatic compact path. These
are mutually exclusive geometric representations; they are not combined.

When the resolved coverage backend is Triangle, `compact_surface=None` retains
that backend's existing catalogue representation. It does not request compact
improvement. Explicit budgeted shared-edge optimization requires both resolved
`mesher="dtcc_mesher"` and `CompactSurfaceOptions.backend="dtcc_mesher"`; unsupported
combinations raise before fitting or terrain meshing. No backend is switched.
Triangle/GEOS face adapters remain explicit zero-point or interior-only comparison
choices. Their polygon-only retessellation cannot retain an existing interior
point during shared-edge insertion, so budgeted shared-edge combinations fail at
option construction. `shared_edge_points=False` retains interior-only comparison.

Every compact mode needs the native `improve` operation, including zero-point
flips. Before reconstruction, a small valid-mesh probe checks the actual loaded
operation and bounded `segment_trials` result. Missing/incompatible capability
raises an actionable error; requested optimization is never silently disabled.
The portable dependency pin has not been invented or replaced with a local path.
**Release gate:** the required Mesher changes remain uncommitted. A clean install
at the pinned revision does not provide this capability. See
[installation](installation.md#compact-catalogue-meshing) for the exact local build.
Flat/LOD1 and ordinary Mesher calls do not require this additional capability.

Public mesh-size/angle parameters govern the reference coverage/terrain. They do
not impose a building-skin size or a global angle-refinement target. Changing
those inputs can change terrain attachment and thus T0. Comparisons must hold
reference settings and input terrain fixed. No corner shells, building-size
subdivision or automatic refinement cascades run in compact improvement.

## Two stages and connected groups

Core constructs minimal source-owned roof/wall polygons directly from cleaned
cells, fitted planes and the established native shell base. A canonical 3D edge
graph removes jointly redundant collinear samples while retaining corners,
creases, ownership changes, height crossings and necessary junctions. Equal XY
at different heights is not a shared 3D point.

Groups follow shared building faces/edges and connected cleaned roof projections,
including unequal-height neighbours. Terrain never joins independent buildings.
Each group has independent source/point/triangle/search budgets; a preceding or
distant group cannot spend them. Deterministic face/segment ordering and cached
incident versions govern selection. Source keys identify a group deterministically.

| Limit | Stage A | Stage B standard | Stage B doubled |
|---|---:|---:|---:|
| Optional points per source | 8 | 4 additional | 8 additional |
| Unique points per group | 32 | 8 additional | 16 additional |
| Interior points per face | 2 | 0 additional | 0 additional |

`max_face_vertices`, `max_source_vertices` and `max_group_vertices` may reduce
Stage A caps; they cannot exceed 2/8/32. `extension_source_vertices` and
`extension_group_vertices` control Stage B up to 8/16; defaults are 4/8.
`adaptive_extension=False` is the fixed policy. With `adaptive_extension=True`
and the standard 4/8 extension caps, Core completes that exact fixed-standard
pass first, then assesses capped source/edge intervals in the same search state.
Each effective source/group allowance increases by exactly the amount required
for an accepted transaction, up to **total Stage B** hard limits of 8/source and
16/group. All incident sources are charged. Failed geometry/quality trials or
remaining poor triangles alone never raise an allowance. Explicit nonstandard
numeric extension caps retain fixed semantics, including 0/0 (Stage A only) and
8/16 (fixed doubled), regardless of this boolean. There is no implicit cap reset.

The fixed standard policy remains the default: adaptive accepted 17 additional
points in cases 4, 6, 10, 11 and 12, but its final meshes matched fixed doubled
bitwise in all 29 cases, with no change in the 18 fully retained LOD2 cases.
This sample shows no final-mesh advantage over the simpler fixed doubled control. Adaptive is an explicit opt-in; its standard-result
prefix can help policy analysis even when its final mesh matches fixed doubled. A shared
point counts once per group and against **every** incident source. A single-source
standard group can therefore add at most four edge points. Interior and shared
points are reported separately. No case-specific budgets appear in production.

T0 is the fresh zero-optional-point building triangle count. Stage A is limited to
`floor(1.5*T0)` and produces baseline B with TB building triangles. Stage B starts
from that same in-memory state and uses edge points only. Its ceiling is
`min(floor(1.5*T0), TB + floor(.20*TB))`. These limits apply per connected group;
their sum also bounds total building growth. Terrain counts never dilute them.
Complete incident roof/wall costs are charged before commit. These are engineering
limits, not geometric laws. Allowances are ceilings, never targets.

Each complete stage/group has at most 512 native one-point requests, 256 joint trials and 60 seconds
of optional search, finishing an in-flight native operation. The finite search
stops on guarded finite-set exhaustion or binding budgets. Zero-point quality
measurements/retessellation flips do not count as proposal requests. Adaptive
continuation reuses valid native/proposal caches and incident versions; its
counters and original TB-based ceiling span both phases. No saved mesh is loaded
and Stage A is never restarted. A private standard-state recovery copy also
preserves the prefix if final whole-shell validation rejects the extension. No scheduler or
minimum-angle-to-completion loop is introduced.

## Native improvement and atomicity

Mesher owns generic planar CDT, legal flips, bounded candidate generation,
insertion, native topology and quality measurements. Its stateless
`improve(mesh, max_added_points=0|1, splittable_segments=(), allow_interior=True)`
returns an owned mesh, stable input point identities, split-segment provenance
and up to four validated segment trial parameters. The input remains borrowed
and unchanged; C callers free result-owned allocations, Python receives NumPy-owned
copies. Quality optimization can be non-Delaunay; ordinary CDT retains its contract.
Read [Mesher's operation documentation](../../dtcc-mesher/docs/improve.md).

The single energy is `E = sum((1-q)**2)`, where
`q = 4*sqrt(3)*area/sum(edge_length**2)`. The public `penalty` field reports E.
The severe-only historical score is retained in reports for comparison. Insertion
requires gain greater than .015 per **total** added triangle. Every incident face
must retain its minimum q and absolute below-5°/below-10° counts and areas within
the unchanged numerical policy. A remote unavoidable worst triangle need not
improve. Incidence is diagnostic, not the objective.

Core maps one returned parameter onto one original 3D segment, projects the same
point into every incident face and retessellates each with no extra points. It
tries lower native ranks only after joint rejection, then selects among requests
by combined benefit/cost. All incident caches, edge sequences and source charges
commit atomically. A failed optional trial leaves the validated state intact.
Unexpected native errors, geometric rejection and budget termination have distinct
diagnostics. Stage-boundary validation can restore its previous compact state;
optional errors cannot flatten a block or alter reconstruction decisions. Invalid
initial compact geometry raises explicitly instead of being presented as success.

Roof domains/planes, cleaned footprints, exposed wall polygons, source heights,
markers and fitting/fallback choices are fixed. No terrain points are proposed.
Reference terrain outside the contact star is bitwise preserved; the compact
terrain tessellation is invariant across all five optional-point controls. Plane,
coverage and shared-identity tolerances stay 1e-7 m, max(1e-7, area*1e-9) m² and
1e-8 m. Shell orientation, closure, vertex links, buried walls and ownership are
validated. Merged and separate outputs carry the same source/face semantics.

Normal diagnostics contain group/stage counts, source charges, effective
allowances, individual escalation reasons/costs/G, the standard-result snapshot,
quality gains and rejection summaries. Every stage is recorded, including disabled
Stage B. Disabled stages do not copy or rebuild the shell or assess edge candidates.
`stopping` is the authoritative additive record:

- `primary_reason`: `disabled`, `no_eligible_shared_edges`,
  `evaluated_candidate_exhaustion`, `source_point_allowance`,
  `group_point_allowance`, `hard_source_point_limit`, `hard_group_point_limit`,
  `building_triangle_ceiling`, `native_request_limit`, `joint_trial_limit`,
  `runtime_limit` or `unexpected_numerical_failure`.
- `binding_limits` retains simultaneous limitations. Search limits take precedence
  in the primary reason, followed by point, triangle and numerical limitations;
  the complete list and rejection counts explain the result.
- `evaluated_finite_set` is false for budget/failure/incomplete assessment stops.
  A finite-set plateau never proves geometric impossibility or global optimality.
- `assessment` records passed geometry/incident-quality trials, fully evaluated
  triangle-blocked candidates, unassessed allowance-capped faces and incomplete
  search. Triangle rejections have `passed_geometry_quality=True`; quality benefit
  is not mislabeled as absent. `best_fully_evaluated_allowance_blocked` gives
  measured G and cost for a still-uncommitted, version-valid candidate blocked by
  the final point allowance (the hard cap during adaptive continuation). It is
  null when no such fully evaluated candidate is known. Accepted historical gains
  belong to `escalations` and `optional_vertices`, not this stopping field.
  Controls do not spend extra requests merely to populate this field.

The legacy `termination` remains, mapping search/point limits to existing budget
strings. Triangle rejection now reports `triangle_budget`; finite evaluated
exhaustion remains `no_beneficial_transaction`. Consumers should use `stopping`
for precise interpretation. Rejection counts summarize attempts, not unique
geometric positions. Accepted transactions record `phase` and `gain_per_triangle`.
E can decrease while below-20° counts/areas or vertex incidence worsen; those
metrics remain diagnostics, separate from the unchanged 5°/10° incident guards. `CompactSurfaceOptions(detailed_audit=True)` also stores
native outcomes, reduced candidate inputs, full edge lineage and face rings.
Auditing does not alter meshing decisions.

## Reproduction and evaluation

The fixed manifest is `.agent/work/compact-integration/case-manifest.json`. It
contains all original 13 cases plus 16 training-only accepted shape-coverage cases,
selected before meshing-policy evaluation. Inputs and source/native binaries are
hashed. Holdout observations are never used for fitting or mesh policy selection.
Fitting fallbacks remain visible and are not counted as successful LOD2 cases.
The matrix uses public entry points from city inputs, with no optimizer checkpoint,
monkey-patching, intermediate mesh reload or sandbox import in production.

```sh
uv pip install --python .venv/bin/python --no-deps ../dtcc-mesher
cmake -S ../dtcc-mesher -B .agent/work/compact-adaptive/native-build \
  -DDTCC_MESHER_BUILD_TESTS=ON -DDTCC_MESHER_BUILD_DEMOS=OFF \
  -DDTCC_MESHER_BUILD_PYTHON=OFF
cmake --build .agent/work/compact-adaptive/native-build -j 4
ctest --test-dir .agent/work/compact-adaptive/native-build --output-on-failure
export MPLCONFIGDIR=/tmp/dtcc-pilot-mpl
export DTCC_MESHER_TEST_CLI="$PWD/.agent/work/compact-adaptive/native-build/dtcc_mesher"
export DTCC_MESHER_TEST_C_API_DUMP="$PWD/.agent/work/compact-adaptive/native-build/dtcc_mesher_c_api_dump"
export DTCC_MESHER_TEST_CASE_DIR="$PWD/../dtcc-mesher/tests/cases"
.venv/bin/python -m pytest ../dtcc-mesher/tests/python
.venv/bin/python -m pytest tests/builder
```

The gallery/report is `.agent/work/compact-integration/report.html`. Earlier
minimal, interface and edge-extension reports/meshes remain read-only historical
references; revised fresh insertion order is not expected to reproduce them.
This validation covers representative cases, not a full-city run.


### Adaptive comparison reproduction

The adaptive evaluation reuses the same manifest byte-for-byte and preserves the
integration directory as a frozen baseline. The five controls are `zero`, `A`,
`standard`, `doubled` and `adaptive`. Both public entry points reconstruct directly
from city inputs. Detailed auditing repeats adaptive on all 29 cases; representative
repeats check determinism. Fixed controls are compared bitwise with frozen meshes.
Use `--jobs 1` for isolated runtime measurements; parallel jobs are evaluation
workers only and are unrelated to the production search.

```sh
export MPLCONFIGDIR=/tmp/dtcc-pilot-mpl
.venv/bin/python sandbox/compact_adaptive_evaluation.py --all --jobs 1 \
  --output .agent/work/compact-adaptive/reproduced
.venv/bin/python sandbox/compact_adaptive_evaluation.py --all --variant adaptive \
  --detailed --jobs 1 --output .agent/work/compact-adaptive/reproduced-detailed
# Repeat representative cases with the same inputs/binary (example).
.venv/bin/python sandbox/compact_adaptive_evaluation.py --case 11 --variant adaptive \
  --output .agent/work/compact-adaptive/reproduced-repeat
.venv/bin/python sandbox/compact_adaptive_validate.py \
  --input-dir .agent/work/compact-adaptive/reproduced \
  --detailed-dir .agent/work/compact-adaptive/reproduced-detailed \
  --repeat-dir .agent/work/compact-adaptive/reproduced-repeat
.venv/bin/python sandbox/compact_adaptive_report.py \
  --input-dir .agent/work/compact-adaptive/reproduced \
  --detailed-dir .agent/work/compact-adaptive/reproduced-detailed
```

The new gallery is `.agent/work/compact-adaptive/report.html`; its policy summary
separates 29 valid compact cases from fully retained LOD2 fits, includes every
attempted case, accepted-step quality/triangle tables and trade-offs. The API
rationale is `.agent/work/compact-adaptive/api-rationale.md`. Native Mesher APIs,
borrowed inputs/owned result semantics, objective and thresholds are unchanged.
No generic candidate or optimization implementation is added to Core.
