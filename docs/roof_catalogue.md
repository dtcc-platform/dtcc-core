# LOD2 city reconstruction

`builder.build_city_surface_mesh(..., lod="lod2")` reconstructs
source roofs from observations, transfers them to cleaned source cells, and
builds roofs, exposed walls and terrain as one conforming surface. The default
`lod="lod1"` keeps the existing LOD1 pipeline. `lod` requests the
output representation; input footprints and observations need no LOD selection.
`lod2_method="catalogue"` is the default LOD2 algorithm. Select
`lod2_method="ransac"` to compare the existing RANSAC reconstructor through the
same city surface meshing path.

For a city model with individual closed building geometries, use
`datasets.city(bounds=bounds, lod="lod2")`. It uses the same catalogue fitter;
`lod2_method="ransac"` selects RANSAC. The dataset still defaults to LOD1.
`city.build_lod2_buildings()` and `builder.build_lod2_buildings(buildings)` also
default to catalogue reconstruction. Individual city-model shells retain the
original footprints, without joint cleaning or terrain attachment. Rejected
fits remain LOD1, including footprints with holes unsupported by this shell
builder. Catalogue outcomes are stored in `building.attributes["roof_reconstruction"]`.
CityJSON export retains both the reconstructed LOD2 and fallback LOD1 geometries.

To mesh supplied CityJSON geometry, use
`city.flatten_geometry(GeometryType.LOD2).mesh()`. This triangulates the existing
surfaces directly; it does not reconstruct roofs or create a joint terrain
boundary. See `demos/build_cityjson_surface_mesh.py`.

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
interfaces, compact mesher and stage contracts. RANSAC uses the shared
building reconstructor. Its planar roof domains retain their fitted
planes and interior ownership. On cleaned footprint fringes, internal domain
boundaries extend as straight lines and exterior fragments use the nearest
original roof domain. This is a geometric extrapolation, not additional evidence.
Stepped roofs retain distinct domain heights. Unsupported fits use explicit LOD1
models; invalid transfer is rejected during roof admission. Assembly errors raise
without flattening admitted roofs.
The `roof_catalogue` diagnostics contain fitting and meshing results for both
methods. The audit fields `requested_lod` and `lod2_method` identify the request.

RANSAC least-squares plane refinement uses centered coordinates, so dense
observations in large projected map coordinates retain both slope components.
When more than two planes are detected and footprint decomposition does not
succeed, `geometry_builders/roof_plane_selection.py` searches subsets of the
maximum six detected planes. Lower envelopes support continuous gable/hip roofs;
upper envelopes support valleys. Selection requires at least 80% of observations
within 0.4 m vertically, at least 20 supporting observations and 2 m² per face,
and at least half of each selected plane's inliers inside its assigned domain.
Candidates minimize mean absolute residual clipped at 1 m plus 0.01 m per plane.
The search preserves fitted planes, the footprint, and shared junctions; it
does not introduce catalogue templates. Near-horizontal subsets remain with the
stepped-flat path. Accepted subsets are recorded in each source's
`diagnostics.plane_selection`, including detected/selected counts and support.
These envelopes do not cover arbitrary combinations of intersecting or stepped
roofs: unsupported cases still use decomposition, stepped-flat reconstruction,
or explicit LOD1 fallback. Face-specific triangulation failures are attributed
to their building so they cannot silently flatten the entire city.

`GeometryType.LOD1` and `GeometryType.LOD2` have the same output meaning as the
strings. The previous `lod` input selector is no longer exposed on surface/volume-mesh
calls: migrate LOD0 and per-building selectors to the desired output `lod` and let
source geometry resolve internally. Flat meshing still selects stored geometry.

The [Höganäs comparison script](../sandbox/compare_lod2.py) reconstructs the same
LM footprints with catalogue, RANSAC, and standalone Roofer roofs, each using
either LM LiDAR or the supplied dense point cloud. It exports six building
meshes and the supplied reference models to VTU. Roofer must be on `PATH` or
configured through `ROOFER_EXECUTABLE`; see the [setup instructions](../sandbox/compare_lod2.md#roofer-setup).
Configure `DATA_DIR` and `BOUNDS` at the top:

```bash
python sandbox/compare_lod2.py
```

See [comparison setup, provenance and metrics](../sandbox/compare_lod2.md).
A compatible local Mesher build is required; see [installation](installation.md).

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

## Geometry, admission and validation

`RoofModel` retains source identity, parameters, local origin, planes and
JSON-ready fitting diagnostics. `model.patches(target)` returns polygons with
holes, affine coefficients and source/patch/block identifiers. No fitted
callables or pickled models are stored on Buildings.

Both methods use the [shared LOD2 pipeline](lod2_pipeline.md). Source-labelled
cleaned regions, terrain attachment, candidate roof admission, common 3D face
boundaries, generic triangulation and bounded refinement have distinct contracts.
Cleaned outlines are authoritative. Roof ridges, steps and affine height crossings
are mandatory constraints; physical height steps retain separate vertices.
Exposed walls are constructed directly from the incident affine profiles and
block foundation. A complete intermediate building shell is not required.

One bounded quality policy coordinates optional roof-edge samples across every
incident roof/wall face. Its caps are 2 interior points per face, 8 points per
source, 32 per connected group, and `floor(1.5*T0)` building triangles per group.
Mesh-size and angle parameters govern terrain. See the [options and audit
contract](compact_mesher_interface.md) for quality guards and stopping reasons.

Numerical calculations use local metre coordinates, 100 nm XY overlay precision
and 10 nm coincidence. These scales do not define physical cleaning detail.
Plane, coverage, winding, edge pairing, vertex-link and degeneracy checks apply
before and after meshing. The completed world-coordinate surface is checked again.
Quality and fit accuracy are separate from topological validity.

## Roof admission and errors

Rejected individual fits use the source's declared grouped LOD1 level. This is
different from an accepted `flat` candidate estimated from observations. Candidate
transfer is validated against its cleaned domain. Admission checks all affine
extrema against the actual terrain foundation. Inadmissible candidates select
the declared source model before shared-face construction, recording the reason;
planes are never lifted or clamped.

Accepted roofs cannot change because triangulation or assembly fails. There is
one assembly attempt, no whole-block flat retry and no alternate recovery chain.
Failures carry their responsible stage and available block/source/face/edge IDs.
Invalid meshes are never returned as successful refinement stops.

`stage_audit["roof_catalogue"]` contains source fitting and admission outcomes,
subdivision changes, partition diagnostics, terrain bases, assembly timings,
strict topology checks, quality by face kind and final marker mappings. Roof
marker `N+m` refers to patch `m`; wall marker `m` refers to its owner patch or a
representative contributor of the same source. Shared-face semantic regions
retain contributor markers. A patch is not a new building. Separate outputs
retain these ownership mappings and are excluded from volume meshing.

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
