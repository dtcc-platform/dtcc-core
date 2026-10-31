# Tagged volume-mesh facets

`VolumeMesh.boundary_faces` is an integer (N, 3) array of vertex indices;
`boundary_markers` is an aligned integer (N,) array. A mapping from three vertex
indices to an integer marker is also accepted. Native loading normalizes either
input to aligned arrays. Cell markers remain separate.

Canonical `.dtcc` and `.dtccpkg` artifacts preserve these arrays, including
negative marker values. XDMF/HDF5 stores the same facet-label pairs. Every
facet must contain three distinct, in-range vertices; duplicate facets,
noninteger labels and mismatched lengths are rejected before a native file is
replaced. This validates the representation, not whether a tagged face is an
exterior boundary or whether a particular PDE boundary condition is appropriate.

The wire schema adds optional fields 4 and 5 to VolumeMesh. Untagged artifacts
omit both fields and remain readable. Tagged artifacts require this reader;
older strict Core readers reject their unknown fields. Regenerate Python
bindings from `dtcc_core/schemas/dtcc.proto` with a Protobuf 5.29-compatible
compiler. Do not edit generated bindings manually.

```python
mesh.boundary_markers = {(0, 1, 2): -2, (0, 1, 3): 17}
mesh.save("mesh.dtcc")
restored = dtcc_core.io.load_volume_mesh("mesh.dtcc")
labels = dict(zip(map(tuple, restored.boundary_faces), restored.boundary_markers))
```

A dataset package additionally retains its DatasetContext. Native mesh files
remain standalone artifacts, rather than numerical restart checkpoints.

## Singleton boundary corners

`dtcc_core.builder.find_boundary_corners(mesh, *, boundary_markers=None)` returns
`BoundaryCornerReport`: aligned one-dimensional integer arrays `vertex_indices`
and `cell_indices`, sorted by vertex index, and a `num_corners` property. A corner
is a vertex incident to exactly one tetrahedron, with all three faces meeting
that vertex selected. On a conforming mesh these faces are exterior. Several
vertices may belong to the same cell. Reports are transient snapshots and do not
track subsequent mesh edits.

With `boundary_markers=None`, detection uses geometric connectivity and needs no
tags. An explicit collection selects a vertex only when all three incident face
labels belong to that collection. An empty collection selects nothing. Nonempty
filters require an unambiguous label on every candidate face; missing labels
raise an error even when other candidate faces would exclude the vertex. Partial
tags away from candidate faces are allowed. Facet matching ignores stored order
and orientation. Malformed shapes, noninteger or out-of-range connectivity,
repeated cell vertices, and malformed or duplicate facet labels are rejected.

The filter's physical meaning belongs to the caller. For strongly imposed P2
velocity on three incident faces, all ten original tetrahedron velocity nodes
are constrained; the singleton P1 pressure basis has no free velocity test
support. A flow solver should select its actual strong full-vector Dirichlet
faces, including prescribed inlets where appropriate. Slip, weak enforcement,
and open boundaries do not meet that definition. Marker sign alone has no
physical meaning. The detector does not validate manifoldness, geometric quality,
Taylor–Hood stability, pressure accuracy, or nonlinear convergence. An interior
edge can join boundary vertices; endpoint membership is not a boundary-edge test.

```python
from dtcc_core.builder import find_boundary_corners, repair_boundary_corners

strong_velocity_markers = {-1, 17}
corners = find_boundary_corners(mesh, boundary_markers=strong_velocity_markers)
if corners.num_corners:
    mesh = repair_boundary_corners(mesh, boundary_markers=strong_velocity_markers)

# Equivalent model operations:
corners = mesh.find_boundary_corners()
prepared = mesh.repair_boundary_corners()
```

## Nonmutating local repair

`repair_boundary_corners(mesh, *, boundary_markers=None)` detects its own current
selection and returns an independent `VolumeMesh`. It accepts no inspection
report. Each distinct selected parent gets one interior barycenter and four
children, retaining one child at the parent's original index and appending three
in deterministic order. Original vertices and unaffected cell indices remain
intact. All original triangular faces remain intact, including shared interfaces;
no neighboring refinement or hanging interfaces are introduced. Growth is one
vertex and three cells per repaired parent, regardless of its corner count.
Detection counts vertex incidence in a linear scan without constructing all mesh
faces or finite-element spaces.

Input arrays remain unchanged on success and failure. A no-op returns an
independent equivalent mesh with a zero-count repair record; attached context,
health and processing history remain unchanged. Repair preserves original
coordinates, boundary face order/orientation/labels, bounds, CRS, transform,
ordinary model metadata, and cell markers (copied from each parent to every
child). Empty cell markers stay empty. Signed child volumes must have the
parent's orientation and sum to its volume; negative parent orientation is
preserved. Checks use local coordinate
scales and reject degenerate, nonfinite, unresolved, collapsed, or unrepresentable
volumes. Large coordinate offsets can make a barycenter unrepresentable and
cause an explicit failure.

Topology-changing repair rejects meshes carrying fields or semantic regions:
repair before attaching these data. Field transfer is not inferred. No-op repair
preserves attached data. Attached Dataset Context is copied with its original
request and dataset identity; repair appends actual processing counts, refreshes
the generation time, and clears prior health findings. The repair record includes
the marker selection and final counts; an attached stage audit retains its
generation history and gains a processing event. Package export computes
fresh artifact sizes and hashes for the returned geometry. Canonical and XDMF
exchange already represent the repaired topology and tags.

## Generation option

The keyword-only boolean `repair_boundary_corners=False` is available on
`builder.meshing.tetgen.build_volume_mesh`, `builder.build_city_volume_mesh`,
`City.build_volume_mesh`, and the `city_volume_mesh` dataset request. It defaults
to false. When true, generation removes all geometric singleton corners,
independent of marker labels and even when boundary-face output is disabled:

```python
from dtcc_core import builder

mesh = builder.build_city_volume_mesh(
    city,
    repair_boundary_corners=True,
    boundary_face_markers=False,
)
```

The city builder forwards the flag to the low-level builder, which applies
repair once after tetrahedralization and before final quality summaries, stage
audit completion, Dataset Context, and export. When repair is requested, the mesh's
`boundary_corner_repair` record contains the request flag, detected corner/parent
counts, repaired counts, and final vertex/cell counts. City stage audits and
dataset processing history include these execution results. With the default
false flag, generation skips corner inspection and adds no repair execution
record. The flag remains in the city audit configuration and validated dataset
request, so true and false requests remain distinct. City volume generation requires
`isolated=False`; independently cleaned building shells are unsupported.
Structured `BoxMesh` has no generation flag and can use the repair operation.

Generation deliberately includes open corners because it cannot know future
solver boundary choices. Selective preprocessing uses the marker-filtered model
operation. The option guarantees the specified corner configuration is absent;
it does not certify a valid or stable solver mesh.

## Measured behavior and limits

On a frozen city mesh with 29,086 cells and 6,398 vertices, seven singleton
vertices in seven cells produce 29,107 cells and 6,405 vertices. Selected-parent
minimum element quality decreases from 0.356 to 0.105 in the children. The global
minimum stays 0.0007468, and mean quality changes from 0.733175 to 0.732812.
Barycenter subdivision can worsen shape quality and is not general remeshing.

Operation timings exclude loading and quality calculation. Values are medians
of three timed operations after one warmup. RSS is sampled at 1 ms intervals;
process peaks include the loaded input, while increments describe operation
allocation. The million-cell fixture comprises 35 disjoint translated copies
of the city, exercising the same repair density without adjacency propagation.

| Fixture | Detection | Repair | Detection peak/increment RSS | Repair peak/increment RSS |
| --- | ---: | ---: | ---: | ---: |
| 29,086 cells | 4.29 ms | 8.72 ms | 177 / 2.93 MB | 178 / 4.73 MB |
| 1,018,010 cells | 156 ms | 314 ms | 371 / 119 MB | 465 / 214 MB |

In FEniCSx 0.11.0 with real PETSc/MUMPS, the repaired all-Dirichlet single and
two-tetrahedron cases each have exactly one pressure null mode (the constant
gauge), and recover matching-force linear pressure with L2 error below 3e-16.
The original meshes trigger the existing uncoupled-pressure rejection.

City comparisons use identical physical settings, P2/P1 spaces, one warmup and
three measured solves per method, on an Apple M5 Pro with 48 GiB memory. The
numerical-library thread limits are one; MPI progress threads remain active.
The reference uses accepted local cell-star refinement to 30,251 cells.

| City solve | Barycenter repair | Local refinement reference |
| --- | ---: | ---: |
| Warm total median [range], s | 53.17 [52.77–53.46] | 57.14 [57.00–57.17] |
| Jacobian assembly total, s | 5.03 | 5.24 |
| Numerical LU total, s | 43.43 | 46.95 |
| Triangular solves total, s | 0.317 | 0.333 |
| Process peak RSS, GB | 2.94 | 3.12 |
| Full-convection Newton iterations | 5 | 5 |

Solver totals exclude mesh preparation and field export. One measured preparation
took 0.071 s for barycenter repair and 1.420 s for local refinement; these include
loading, and the reference also assembles the algebraic pressure check.

Both methods pass the unchanged residual, mass and boundary criteria. Repaired
city residual is 5.09e-6 against a 1.14e-4 limit, dimensionless divergence is
0.00110 against 0.1, flux imbalance is 1.51e-16 against 1e-4, and wall-normal
velocity RMS is zero against 1e-8 m/s. Inlet flux agrees to roundoff at
-100070.660 m³/s. On 3,823 common probes, sampled relative velocity/pressure
differences are 0.450% / 0.820%, with speed RMS difference 0.0186 m/s. Pressure
is traction-determined and no mean shift is removed. These are differences from
a numerical reference, not exact city errors or physical validation.

These checks support this subdivision on the tested meshes. They do not establish
general inf-sup stability, accuracy on other cities, million-cell solver
performance, turbulence validity, or scalable iterative/MPI solver behavior.
Sim integration should repair the native mesh before finite-element conversion
and partitioning, using actual strong velocity markers, and retain algebraic
pressure/gauge and full-field checks. Distributed meshes need an explicit native
preprocessing path. Core introduces no solver dependency or implicit gather.
