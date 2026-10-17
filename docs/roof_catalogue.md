# LOD2 city reconstruction

`builder.build_city_surface_mesh(..., lod="lod2")` reconstructs
source roofs from observations, transfers them to cleaned source cells, and
builds roofs, exposed walls and terrain as one conforming surface. The default
`lod="lod1"` keeps the existing LOD1 pipeline. `lod` requests the
output representation; input footprints and observations need no LOD selection.
`lod2_method="catalogue"` is the default LOD2 algorithm. Select
`lod2_method="ransac"` to compare the existing RANSAC reconstructor through the
same city surface meshing path. The standalone `build_lod2_buildings` API is
unchanged.

```python
from dtcc_core import builder
from dtcc_core.builder.roofs import RoofFitOptions, supported_shapes

print(supported_shapes())  # ('flat', 'shed', 'gable', 'hip')
audit = {}
mesh = builder.build_city_surface_mesh(
    city,
    lod="lod2",
    roof_options=RoofFitOptions(shapes=("flat", "shed", "gable", "hip")),
    stage_audit=audit,
)
```

`city.build_surface_mesh` accepts the same reconstruction options. `lod2_method`
and the fitting/meshing options are keyword-only. `roof_options=None`
uses the defaults below. Unknown methods or shape keys fail immediately.
Passing catalogue fit options with LOD1 or RANSAC is an error.
Both methods support `isolated=True` and `isolated=False`. Joint cleaning builds
one conforming terrain/building surface. Isolation cleans buildings independently,
meshes each fitted skin down to its lowest supporting terrain level, and combines
it with a separate terrain mesh. Buildings may overlap; their vertices are never
welded to neighbouring buildings or terrain. Both layouts use the same roof
fitting, transfer, fallback and bounded meshing routines.

With `merge_meshes=False`, the joint layout returns terrain and patch shells;
validation uses their welded union. Isolation returns terrain followed by one
skin per retained cleaned building, with the `building` face field identifying
its output index. Compact skins are constructed directly from fitted domains and
a base height,
without triangulating a temporary terrain ring. Validation permits open edges
only at that base and still checks shared edges, orientation and manifoldness.
Diagnostics under `roof_catalogue.building_assemblies`
record those individual assemblies; there is no claim of whole-city closure.
Isolated meshes remain unsuitable as volume boundaries.

With the resolved `dtcc_mesher` backend, catalogue calls use minimal compact
roof/wall polygons and bounded two-stage improvement by default. Use
`CompactSurfaceOptions(budgeted=False)` to disable optional points. Mesh sizes
govern the reference terrain, not building-skin refinement. Read the
[compact meshing contract](compact_mesher_interface.md) for group budgets,
triangle ceilings, backend restrictions, detailed auditing and the current
compatible-Mesher release gate. These controls do not change roof fitting.

## Selecting output LOD and comparing methods

```python
lod1 = city.build_surface_mesh(lod="lod1")
catalogue = city.build_surface_mesh(lod="lod2")
ransac = city.build_surface_mesh(lod="lod2", lod2_method="ransac")
```

Both LOD2 methods use the source-specific observations, cleaned cells, shared
interfaces, compact mesher and fallback validation. RANSAC uses the existing
building reconstructor unchanged. Its planar roof domains retain their original
planes and interior ownership. On cleaned footprint fringes, internal domain
boundaries extend as straight lines and exterior fragments use the nearest
original roof domain. This is a geometric extrapolation, not additional evidence.
Stepped roofs retain distinct domain heights. Unsupported fits use explicit LOD1
fallbacks; invalid transferred partitions use the existing block fallback.
The `roof_catalogue` diagnostics contain fitting and meshing results for both
methods. The audit fields `requested_lod` and `lod2_method` identify the request.

`GeometryType.LOD1` and `GeometryType.LOD2` have the same output meaning as the
strings. The previous `lod` input selector is no longer exposed on surface/volume-mesh
calls: migrate LOD0 and per-building selectors to the desired output `lod` and let
source geometry resolve internally. Flat meshing still selects stored geometry.

The [comparison demo](../demos/compare_lod2.py) makes two dataset requests
for the same 200 m Gothenburg block (13 source buildings in the validation dataset)
and saves catalogue and RANSAC surface meshes for ParaView:

```bash
python demos/compare_lod2.py
```

Outputs are `output/lod2_catalogue.vtu` and `output/lod2_ransac.vtu`.
It downloads footprints and LiDAR through the dataset API. A compatible local
Mesher build is required; see [installation](installation.md).

## Supported catalogue

The ordered registry in `dtcc_core.builder.roofs.catalogue.REGISTRY` is explicit.
Each family has its own module, parameter definitions, initializations, bounds
and plane conversion.

| Shape key | Initial family | Free parameters |
| --- | --- | --- |
| `flat` | Horizontal roof | Elevation, 1 |
| `shed` | One inclined plane | Reference elevation and two slope components, 3 |
| `gable` | Symmetric sides, centred ridge in an oriented bounding rectangle | Orientation, eave elevation, pitch, 3 |
| `hip` | Common rectangle eave elevation, symmetric opposing faces | Orientation, eave elevation, side pitch, end pitch, 4 |

Elevations are in metres in the input vertical reference system. Slopes/pitches
are rise/run. Orientations are radians counterclockwise from the projected x
axis, along the ridge/rectangle axis. The reference point is the source
footprint centroid. Hip ridge length follows from rectangle dimensions and
pitches. A point ridge is a geometric limiting case, not another roof family.

Gable and hip models use the source footprint's oriented bounding rectangle;
their plane domains are then clipped to the actual cleaned cell, including
holes and disconnected pieces. On an irregular outline, this can give nonlevel
eaves. Rectangularity below 0.85 produces an applicability warning for these
families. This warning is distinct from residual or support rejection.

Asymmetric gables, mansards, gambrels, dormers, overhangs, curved roofs and
compound-roof decomposition are outside this catalogue. A successful optimizer
or complete footprint coverage alone does not establish a supported fit.

## Evidence and fitting

Cleaning runs once at the requested detail scale, or reuses a matching cleaning
record. Selection by `min_building_area` happens on cleaned buildings before any
roof subdivision. Labels, source footprints, cell precedence and discarded
partitions are respected. A merged parent's pooled point cloud is never fitted
to each part: observations must belong to the corresponding original building
or BuildingPart. Missing observations, source correspondence or ground height
produce explicit LOD1 fallback diagnostics. Original footprints preserved on
parts are preferred for orientation and fitting. A previously cleaned source
without a preserved original footprint uses LOD1 fallback and records the
missing orientation evidence; the cleaner does not retain original outlines
for every standalone building.

The fitter is independent of the city mesher. It uses local XY coordinates for
projected-coordinate stability, equal weights per occupied 1 m cell, bounded
least squares with a soft-L1 loss of scale 0.2 m, and deterministic orientation
candidates from footprint edges, covariance and the bounding rectangle. Gable
and hip orientations optimize within 7.5 degrees of each candidate; pitches
range from 0.035 to tan(60 degrees). Shed slope components are bounded by
tan(60 degrees)/sqrt(2). Reference/eave heights lie between ground + 0.5 m and
the larger of ground + 1.5 m and the observed 98th percentile + 2 m.

When enough observations exist, deterministic 2 m spatial blocks supply an
internal validation set. Otherwise selection uses penalized training error and
records `inner_spatial_split=False`. Model selection adds 0.025 m per free
parameter to weighted pseudo-Huber error. External benchmark holdouts never
enter the production fitter.

`RoofFitOptions` defaults are:

| Option | Default | Meaning |
| --- | --- | --- |
| `shapes` | all four registry keys | Ordered candidate families |
| `min_points` | 12 | Minimum number of source observations |
| `median_error` | 0.35 m | Maximum training median absolute error |
| `p90_error` | 1.5 m | Maximum training 90th percentile absolute error |
| `face_inlier_distance` | 0.6 m | Residual limit for face support |
| `minimum_clearance` | 0.25 m | Minimum roof-to-source-ground clearance |
| `max_nfev` | 100 | Maximum function evaluations per optimization |

Every face occupying at least 10% of the source footprint requires at least
max(3, floor(0.01 * observation count)) inliers. Convergence, finite parameters,
residuals, support and clearance are checked separately. Thresholds are starting
engineering defaults, not calibrated confidence probabilities.

## Geometry and validation

`RoofModel` retains source identity, parameters, local origin, planes and
JSON-ready diagnostics. `model.patches(target)` returns `RoofPatch` objects:
polygons with holes, affine coefficients and source/patch/block identifiers.
The lower-envelope representation of these four families stays inside the roof
package. Neither the native shell builder nor interface meshing knows catalogue
shape names. No fitted callables or pickled models are stored on Buildings.

Patch boundaries form a coordinated XY segment graph. Ridge and hip edges are
constraints. Shared boundaries receive the union of breakpoints from both
sides, including every interior root of the difference of adjacent affine
height profiles. Coincident intervals produce no step wall. At a crossing,
the higher owner switches and a wall panel can collapse to a triangle. Outer
walls descend to the existing flattened block base on terrain; internal walls
cover only the exposed step. Continuous interfaces share 3D vertices, while
height discontinuities retain separate Z vertices joined by walls.

Both existing `dtcc_mesher` and optional `triangle` backends are supported. The
complete partition is triangulated together; backend Steiner vertices evaluate
the marked patch plane. Triangle containment and area coverage verify that
creases survive triangulation. Native wall generation shares the legacy
vertical-line reconciliation and quality-aware zipper triangulation, generalized
to bottom/top heights at both endpoints. This supplies the validated reference
shell. The normal dtcc_mesher path then constructs compact geometry; Triangle
retains the connected planar-face wall method below.

Numerical XY noding uses a 1e-7 m grid, independently of footprint detail and
fit tolerances. Native shared-vertex heights within 1e-6 m are reconciled for
roundoff. Both backend coverages use a local coordinate frame; roof and terrain
boolean operations share one precision model. Numerical faces created by
rounding use the original patch ownership, while real courtyards remain ground.
Area coverage tolerances include the explicit perimeter times coordinate-grid
roundoff bound; overlap checks retain their independent tolerance. Constraints
closer than the native 1e-6 m coincidence scale reject the affected block before
triangulation, rather than being silently merged. Physical
`min_building_detail` never snaps fitted elevations. A
collapsed patch, ambiguous very short interval, lost constraint or unsupported
thin geometry must fall back or raise; it is never silently removed. Generic
shell construction bounds vertical-line subdivision to 1024 levels, including
when physical detail is zero, to bound wall-triangulation memory. Exceeding this
budget invokes the same fallback/error policy.

Shell validation checks interior cracks, paired-edge orientation, nonmanifold
edges, vertex links, zero-area triangles, upward roof orientation and buried
walls. The buried-wall audit evaluates the roof triangle plane at the probe.
Triangle quality minima and percentiles are reported separately from validity;
validity does not guarantee a useful observation fit or an engineering-quality
mesh in every application. The existing surface quality contract is applied:
zero-area triangles are errors, quality below 0.02 and extreme aspect ratios
produce explicit warnings. These thresholds are unchanged. Audits retain the
worst face, its vertices, edge lengths and aspect ratio. Remaining constrained
wall faces retain explicit quality warnings; outputs are not certified for
engineering use or volume meshing.

## Connected planar wall faces

This section describes the legacy/Triangle wall representation. Normal compact
catalogue meshing uses the separate bounded contract linked above.

`meshing.wall_faces` consumes the valid exposed panels from the shared native
shell builder, with source ownership. It has no catalogue-family dispatch.
Edge-connected, equally oriented panels merge only within the same source and
one plane at the existing 1e-7 m numerical tolerance. A fixed reference plane
prevents tolerance drift. Actual corners, different owners, disconnected faces,
courtyards, steps and required surface junctions remain constraints. Concave
outlines and holes use their true polygon union.

Only faces containing triangles below quality 0.02, and their incident faces
when a junction is subdivided, need replacement. DTCC Mesher supplies an initial
constrained triangulation. Interior off-centres refine it using existing SciPy
Delaunay tessellation, accepting a result only if every boundary segment and
the complete polygon coverage survive. Optional Triangle support is not needed
by this wall stage, whichever backend built the roof coverage.

Eave and terrain boundary discretizations remain fixed. A genuinely narrow
face may require shorter vertical corner edges; subdivisions are made once on
the global wall edge and used by every incident face, including across owner
boundaries. They are allowed only on edges shared exclusively by two wall
triangles. Proposed interior junction vertices are graded and reconciled;
existing boundary vertices are never removed. Roof and terrain triangles stay
unchanged. Each coupled group is a transaction: an unsupported boundary,
quality failure or exceeded budget retains its original valid panels.

The internal budgets are 32 refinement rounds, 4096 points per wall face and
100000 added points per city shell. A replacement must meet quality 0.02,
preserve coverage and winding, and pass the unchanged strict shell validator.
A wall-quality rejection cannot reject a fitted roof or invoke reconstruction
fallback. It records its reason under
`stage_audit["roof_catalogue"]["wall_triangulation"]`, independently of source
and block fallbacks. Face records include original and candidate quality,
coverage error, refinement effort, contributing markers and input-angle bounds.
A bound below the target proves a geometric limitation; other rejections are
limitations of this bounded refinement method, not proofs of impossibility.

A triangle spanning former patch seams uses a representative wall marker from
its source. The face field `catalogue_wall_face` and semantic regions with URI
`urn:dtcc:wall-face` map accepted replacements to their source and all contributing
patch markers. Region attributes persist in canonical `.dtcc` files; export the
audit alongside formats that do not preserve semantic regions. Value -1 in the
field denotes an unchanged triangle, whose original marker remains exact.
Separate-mesh output assigns a spanning face to its representative patch shell
and preserves its region metadata. Its welded union is the same conforming
surface. Roof marker semantics and default LOD1 markers are unchanged.

## Rejection and bounded fallback

Rejected individual fits use the ordinary grouped LOD1 levels for their source
cells. This is different from an accepted `flat` candidate estimated from
observations. Block fallback calls the same `_building_level_regions` policy
as ordinary LOD1, including its level grouping and cleaner decisions, while
retaining source ownership and step walls. A single roof is used only where
the ordinary cleaner/LOD1 policy supplies one; an 8 m/14 m block retains both
heights. A previously cleaned city follows the same policy.

Each block receives a numerical coverage admission check. Local failures record
its sources and reason. A city assembly failure may rebuild attributed blocks
once with their height-preserving LOD1 coverage. There are at most two city
assembly attempts, and each block may fall back only once. An unattributed
backend error raises with the audit rather than flattening the city. Failure of
a flat replacement raises `RoofAssemblyError`; invalid shells are never returned.
Triangle catalogue calls retain all input constraints and cap refinement at
50,000 Steiner points. Refinement budget and resulting quality are audited;
the cap does not certify that the requested angle was achieved.

`stage_audit["roof_catalogue"]` contains source outcomes and candidate scores,
block/source mappings, partition diagnostics, assembly attempts, topology and
quality metrics, whole-block fallback reasons, and final patch-marker mappings.
Roof marker `N + m` refers to patch `m`. Wall marker `m` refers to that patch
for unchanged panels and to a representative patch of the same source for
connected-face replacements, with the explicit contributor mapping above;
`N` is the number of meshing patch regions in this path. The mapping records
source ownership separately: a roof patch is not a new building. Separate
meshes are terrain followed by active patch shells, as in the existing region
shell interface. A failed API call retains its audit but returns no mesh.

## Volume meshing

```python
import dtcc_core as dtcc

volume = city.build_volume_mesh(lod="lod2")
# The dataset entry point has the same output-LOD and method selection:
volume = dtcc.datasets.city_volume_mesh(bounds=bounds, lod="lod2")
```

All city surface and volume entry points default to `lod="lod1"`.
LOD2 volume meshing calls the joint surface builder directly, with catalogue
reconstruction by default and `lod2_method="ransac"` available. Fitting,
footprint cleaning, source fallbacks, shared-edge discretization and compact
roof/wall triangulation are identical. `roof_options` and `compact_surface`
are available on the builder and City volume methods.

The existing TetGen stages add the domain cap and side walls, validate the
closed PLC, then tetrahedralize the air domain outside the buildings. There
is no second footprint cleaning or flat-roof extrusion. Cap size seeding
retains its interior points but re-triangulates against the prescribed boundary
if the mesher subdivides it. No additional Mesher API is required.

`stage_audit["roof_catalogue"]` retains roof outcomes and marker-to-source
mappings; the volume attempt records the surface pipeline, PLC and volume
metrics. Roof/wall markers remain patch-based, with terrain `-1`, top `-2`,
and domain sides `-3` through `-6`.

Joint merged LOD2 surfaces are admitted to PLC validation. Isolated meshes
and separate patch outputs retain `volume_boundary=0` and are refused. Use
native `.dtcc` serialization to preserve this field; VTU may omit geometry
fields and is not a volume-validation certificate.

The real hip-building and merged-group checks verify positive tetrahedron
volumes, conforming closed boundaries, preserved roof planes and boundary
areas, and agreement between tetrahedron volume and the closed PLC volume.
These are geometry checks, not a solver-quality guarantee. Small constraints
in the surface can cause sliver tetrahedra and substantial TetGen refinement.
In the hip fixture, a 60 m³ volume cap increased the LOD2 mesh from 1,009 to
56,450 tetrahedra while its worst quality remained poor; explicit quality
refinement did not eliminate that limitation. The integration preserves the
existing TetGen defaults rather than silently changing geometry or applying
an aggressive refinement policy. Inspect quality reports before simulation.

## Adding a family

Implement a module under `roofs/catalogue/` defining `PARAMETERS`,
`candidates(context)` yielding initial/lower/upper vectors, and
`planes(parameters, local_vertices)`. Add its entry to the ordered registry,
document its parameters and add noisy-fit, support and patch-coverage tests.
The shared fitter owns optimization and diagnostics. No meshing changes should
be necessary for another family representable by the roof package's current
plane-domain converter. A family with different domains should implement its
patch conversion within the roof package while preserving `RoofPatch`.

## Tests and saved-data comparison

```bash
uv run pytest tests/builder/test_roof_catalogue.py tests/builder/test_building_shells.py
```

Backend-specific tests explicitly skip an unavailable optional backend. Native
changes require rebuilding (`uv sync --reinstall-package dtcc-core`), not just
an editable Python import. See the native/installation documentation for the
optional Triangle build configuration.

The separate benchmark runner consumes the existing local 341-building cache
and fixed thirteen-case list, including both five-building groups:

```bash
.venv/bin/python sandbox/evaluate_roof_catalogue_integration.py \
  --input .agent/work/roofer-evaluation \
  --cases .agent/work/catalogue-evaluation/cases.json \
  --output .agent/work/catalogue-integration --mesher dtcc_mesher
```

Use `--mesher triangle` and a separate output directory to repeat both complete
areas and the thirteen cases with the optional backend.

It reconstructs both complete areas and compares assembled default LOD1 and
catalogue shells in every case, with held-out LiDAR overlays. The local terrain
is a 4 m raster sampled from classified ground in the manifest's LAZ files.
JSON audits, VTU surfaces, per-source held-out errors, failed assemblies, fallback
reasons and triangle-quality metrics accompany `report.html`. A failed assembled
candidate may be retained by this diagnostic runner and is explicitly labelled;
production never returns it. Caches, figures and large mesh files stay outside
installation and CI. Missing local inputs raise a clear error; no substitute
measurements are invented. The runner depends on the production API, not the
sandbox catalogue prototype or Roofer.


For a wall-only benchmark against preserved catalogue outputs, use the same
entry point with a baseline root containing both backend report directories and
their complete-city `.dtcc` meshes:

```sh
python sandbox/evaluate_roof_catalogue_integration.py \
  --input .agent/work/roofer-evaluation \
  --cases .agent/work/catalogue-evaluation/cases.json \
  --wall-baseline .agent/work/catalogue-wall-baseline \
  --output .agent/work/catalogue-walls-dtcc_mesher --mesher dtcc_mesher
```

Repeat with `--mesher triangle` and a separate output directory. This compares
old catalogue, new catalogue and ordinary LOD1 on both complete cities and all
thirteen cases, retaining roof/wall/terrain metrics, reconstruction outcomes,
geometry checks and paired views. Quality distributions include absolute poor
triangle counts, used and stored vertices, growth and runtime; process memory
high-water marks are labelled separately from isolated allocation measurements.
It preserves the saved training/holdout split and does not rerun Roofer.
