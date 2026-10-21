"""Selectable shared-edge quality-refinement pilot for catalogue surfaces."""

from dataclasses import dataclass, asdict
import time

import numpy as np
import shapely

from ...model import Field, Mesh, SemanticRegion
from .shared_edges import SharedEdges, seed_corner_shells
from .surface_faces import extract_surface_faces
from .surface_refinement_backend import mesh_face
from .wall_faces import _angle_constraint
from .face_validation import validate_face_geometry


@dataclass(frozen=True)
class SurfaceRefinementOptions:
    """Explicit pilot controls; None leaves catalogue meshing unchanged.

    min_angle is the face-refinement target, independent of min_mesh_angle
    used by the city coverage triangulation. Nominal spacing defaults to the
    building mesh size (terrain uses max_mesh_size). Quality exceptions must
    be reported against actual constraints or backend protection, not hidden.
    """

    backend: str = "dtcc_mesher"
    min_angle: float = 20.0
    nominal_spacing: float | None = None
    max_rounds: int = 12
    max_vertices: int = 300000
    max_face_vertices: int = 50000
    max_refine_steps: int = 100000
    max_seconds: float = 180.0

    def __post_init__(self):
        if self.backend not in ("dtcc_mesher", "triangle"):
            raise ValueError("Face backend must be dtcc_mesher or triangle")
        if not 0 < self.min_angle <= 28:
            raise ValueError("Face minimum angle must be in (0, 28] degrees")
        if self.nominal_spacing is not None and (
            not np.isfinite(self.nominal_spacing) or self.nominal_spacing <= 0
        ):
            raise ValueError("Nominal face spacing must be positive")
        for key in (
            "max_rounds",
            "max_vertices",
            "max_face_vertices",
            "max_refine_steps",
        ):
            value = getattr(self, key)
            if (
                not isinstance(value, (int, np.integer))
                or isinstance(value, bool)
                or value <= 0
            ):
                raise ValueError(f"{key} must be a positive integer")
        if not np.isfinite(self.max_seconds) or self.max_seconds <= 0:
            raise ValueError("Refinement time budget must be positive and finite")


def _boundary_requests(face, points, triangles, registry):
    """Map backend boundary points to canonical 3D edges and their parameters."""
    edges = np.sort(
        np.concatenate(
            [triangles[:, [0, 1]], triangles[:, [1, 2]], triangles[:, [2, 0]]]
        ),
        axis=1,
    )
    unique, counts = np.unique(edges, axis=0, return_counts=True)
    boundary = unique[counts == 1]
    if np.any(counts > 2):
        raise ValueError("Face mesher returned a nonmanifold edge")
    ids = np.unique(boundary)
    refs = sorted(set(ei for ei, _ in registry.face_edges[face.id]))
    endpoints = np.array(
        [[registry.vertices[v] for v in registry.edges[ei].endpoints] for ei in refs]
    )
    xy = (endpoints - face.origin) @ face.basis.T
    direction = xy[:, 1] - xy[:, 0]
    length2 = np.sum(direction * direction, axis=1)
    mapped = {}
    for pi in ids:
        p = points[pi]
        parameters = np.sum((p - xy[:, 0]) * direction, axis=1) / length2
        residual = np.linalg.norm(
            p - (xy[:, 0] + parameters[:, None] * direction), axis=1
        )
        possible = np.flatnonzero(
            (parameters >= -1e-10)
            & (parameters <= 1 + 1e-10)
            & (residual <= registry.tolerance)
        )
        if not len(possible):
            raise ValueError(
                f"Boundary point is off every original edge (distance {residual.min():g})"
            )
        k = min(possible, key=lambda j: (residual[j], refs[j]))
        mapped[int(pi)] = (refs[k], float(np.clip(parameters[k], 0.0, 1.0)))
    return mapped, boundary


def _cache_matches(face, result, registry):
    mapped = result["boundary_map"]
    try:
        ids = {pi: registry.locate(ei, t) for pi, (ei, t) in mapped.items()}
    except ValueError:
        return False
    actual = {
        tuple(sorted((ids[int(a)], ids[int(b)])))
        for a, b in result["boundary_segments"]
    }
    expected = {tuple(sorted(e)) for e in registry.boundary(face.id)}
    return actual == expected


def refine_shared_surface(
    mesh, owners, options, *, building_spacing=4.0, terrain_spacing=8.0, adapter=None
):
    """Transactionally refine a shell using one shared discretization per edge.

    No geometry, fitting or reconstruction fallback is performed. Failure
    returns the entire original valid shell and an explicit failed pilot audit.
    """
    started = time.perf_counter()
    audit = dict(
        method="shared_edge_quality_pilot",
        options=asdict(options),
        accepted=False,
        rounds=[],
        face_refinement_backend=options.backend,
        geometry_tolerance_m=1e-7,
        shared_edge_tolerance_m=1e-8,
        faces=[],
    )
    refine = mesh_face if adapter is None else adapter
    current_face = None
    try:
        if len(mesh.vertices) > options.max_vertices:
            raise ValueError(f"Initial vertex budget exceeded ({options.max_vertices})")
        faces = extract_surface_faces(mesh, owners)
        registry = SharedEdges(mesh.vertices, {f.id: f.boundary for f in faces})
        audit["initial_corner_discretization"] = seed_corner_shells(
            registry,
            faces,
            nominal_spacing=options.nominal_spacing or building_spacing,
            max_vertices=options.max_vertices,
        )
        cache = {}
        dirty = set(range(len(faces)))
        for round_number in range(1, options.max_rounds + 1):
            row = dict(
                round=round_number, dirty_faces=len(dirty), new_boundary_vertices=0
            )
            audit["rounds"].append(row)
            requests = []
            for fi in sorted(dirty):
                current_face = faces[fi]
                if time.perf_counter() - started > options.max_seconds:
                    raise ValueError(
                        f"Pilot runtime budget exceeded ({options.max_seconds:g} s)"
                    )
                boundary = registry.boundary(fi)
                ids = np.unique(boundary)
                lookup = {int(v): i for i, v in enumerate(ids)}
                segments = np.array([[lookup[a], lookup[b]] for a, b in boundary])
                points = (
                    np.array(registry.vertices)[ids] - current_face.origin
                ) @ current_face.basis.T
                # Evaluate optional edge points in local coordinates directly;
                # subtracting large world coordinates would add fresh roundoff
                # at every synchronization round.
                for ei, _ in registry.face_edges[fi]:
                    edge = registry.edges[ei]
                    a, b = (
                        np.array([registry.vertices[v] for v in edge.endpoints])
                        - current_face.origin
                    ) @ current_face.basis.T
                    for t, vertex in zip(edge.parameters[1:-1], edge.vertices[1:-1]):
                        points[lookup[vertex]] = a + t * (b - a)
                spacing = options.nominal_spacing or (
                    terrain_spacing
                    if current_face.kind == "terrain"
                    else building_spacing
                )
                if len(points) > options.max_face_vertices:
                    raise ValueError(
                        f"Face input point budget exceeded ({options.max_face_vertices})"
                    )
                try:
                    vertices, triangles, backend = refine(
                        points,
                        segments,
                        current_face.polygon,
                        options,
                        spacing=spacing,
                        reuse_acute_corners=True,
                    )
                except (ValueError, RuntimeError, shapely.errors.GEOSException) as exc:
                    audit["reproducer"] = dict(
                        face_id=fi,
                        kind=current_face.kind,
                        points=points.tolist(),
                        segments=segments.tolist(),
                        polygon=shapely.to_wkt(
                            current_face.polygon, rounding_precision=15
                        ),
                        spacing=spacing,
                    )
                    raise ValueError(
                        f"Face {fi} ({current_face.kind}) mesher failed: {exc}"
                    ) from exc
                if len(vertices) > options.max_face_vertices:
                    raise ValueError(
                        f"Face output point budget exceeded ({options.max_face_vertices})"
                    )
                mapped, output_boundary = _boundary_requests(
                    current_face, vertices, triangles, registry
                )
                # Every input constraint must be present as a union of refined
                # segments; full polygon equality is checked before adoption.
                polygons = shapely.polygons(vertices[triangles])
                union = shapely.union_all(polygons)
                tolerance = max(1e-7, current_face.polygon.area * 1e-9)
                if (
                    current_face.polygon.symmetric_difference(union).area > tolerance
                    or abs(float(shapely.area(polygons).sum()) - union.area) > tolerance
                ):
                    raise ValueError(f"Face {fi} backend lost polygon coverage")
                cache[fi] = dict(
                    points=vertices,
                    triangles=triangles,
                    boundary_map=mapped,
                    backend=backend,
                    spacing=spacing,
                    round=round_number,
                    boundary_segments=output_boundary,
                )
                requests.extend(mapped.values())
            affected, added = registry.reconcile(
                requests, max_vertices=options.max_vertices
            )
            row.update(
                new_boundary_vertices=added, registry_vertices=len(registry.vertices)
            )
            dirty = {
                fi for fi in cache if not _cache_matches(faces[fi], cache[fi], registry)
            }
            row["unpropagated_faces"] = len(dirty)
            row["dirty_face_ids"] = sorted(dirty)
            if not dirty:
                break
        else:
            raise ValueError(
                f"Shared-edge synchronization budget exceeded ({options.max_rounds} rounds)"
            )
        vertices = list(registry.vertices)
        triangle_arrays, markers, face_ids, exempt_flags = [], [], [], []
        for face in faces:
            current_face = face
            result = cache[face.id]
            local, triangles, mapped = (
                result["points"],
                result["triangles"].copy(),
                result["boundary_map"],
            )
            global_ids = np.full(len(local), -1, dtype=np.int64)
            for pi, (ei, t) in mapped.items():
                global_ids[pi] = registry.locate(ei, t)
            # Boundary points use registry geometry, never an independent lift
            # or a final proximity weld. Interior points stay on this plane.
            inside = np.flatnonzero(global_ids < 0)
            if len(vertices) + len(inside) > options.max_vertices:
                raise ValueError(
                    f"Total surface vertex budget exceeded ({options.max_vertices})"
                )
            global_ids[inside] = np.arange(len(vertices), len(vertices) + len(inside))
            vertices.extend(face.origin + local[inside] @ face.basis)
            xy = local[triangles]
            reverse = (xy[:, 1, 0] - xy[:, 0, 0]) * (xy[:, 2, 1] - xy[:, 0, 1]) - (
                xy[:, 1, 1] - xy[:, 0, 1]
            ) * (xy[:, 2, 0] - xy[:, 0, 0]) < 0
            triangles[reverse] = triangles[reverse][:, [0, 2, 1]]
            tri = global_ids[triangles]
            geometry = validate_face_geometry(face, np.asarray(vertices), tri)
            triangle_arrays.append(tri)
            markers.append(np.full(len(tri), face.markers[0], dtype=np.int64))
            face_ids.append(np.full(len(tri), face.id, dtype=np.int64))
            exempt = np.zeros(len(tri), dtype=np.int64)
            exempt[result["backend"].get("quality_exempt_triangles", [])] = 1
            exempt_flags.append(exempt)
            audit["faces"].append(
                dict(
                    face_id=face.id,
                    kind=face.kind,
                    source_key=face.source_key,
                    markers=face.markers,
                    original_triangles=len(face.original_faces),
                    triangles=len(tri),
                    added_interior_vertices=len(inside),
                    spacing=result["spacing"],
                    last_refinement_round=result["round"],
                    backend=result["backend"],
                    angle_constraint=_angle_constraint(face.polygon),
                    geometry=geometry,
                    boundary_edges=[
                        dict(
                            edge_id=ei,
                            forward=forward,
                            vertices=registry.chain(ei, forward),
                        )
                        for ei, forward in registry.face_edges[face.id]
                    ],
                )
            )
        candidate = Mesh(
            vertices=np.asarray(vertices),
            faces=np.concatenate(triangle_arrays),
            markers=np.concatenate(markers),
        )
        ids = np.concatenate(face_ids)
        candidate.fields.append(
            Field(name="catalogue_surface_face", values=ids, association="face")
        )
        candidate.fields.append(
            Field(
                name="refinement_quality_exempt",
                values=np.concatenate(exempt_flags),
                association="face",
            )
        )
        for face in faces:
            candidate.regions.append(
                SemanticRegion(
                    "urn:dtcc:surface-face",
                    id=str(face.id),
                    indices=np.flatnonzero(ids == face.id),
                    attributes=dict(
                        kind=face.kind,
                        source_key=face.source_key,
                        contributing_patch_markers=face.markers,
                    ),
                )
            )
        mismatches = sum(
            not _cache_matches(face, cache[face.id], registry) for face in faces
        )
        if mismatches:
            raise ValueError(
                f"{mismatches} faces have inconsistent shared edge sequences"
            )
        from .roof_validation import validate_shell
        from ..geometry_builders.meshes import _shell_closure_audit

        audit["topology"] = validate_shell(
            candidate, num_patches=len(owners), closure_audit=_shell_closure_audit
        )
        if not audit["topology"]["valid"]:
            raise ValueError(
                "Pilot replacement failed strict shell topology validation"
            )
        audit.update(
            accepted=True,
            shared_edges=len(registry.edges),
            added_boundary_vertices=len(registry.vertices) - len(mesh.vertices),
            added_interior_vertices=len(vertices) - len(registry.vertices),
            shared_boundary_mismatches=mismatches,
            original_vertices=len(mesh.vertices),
            vertices=len(vertices),
            edge_registry=[
                dict(
                    edge_id=i,
                    endpoints=e.endpoints,
                    parameters=e.parameters,
                    vertices=e.vertices,
                    incident_faces=sorted(e.incidents),
                )
                for i, e in enumerate(registry.edges)
            ],
        )
        audit["seconds"] = time.perf_counter() - started
        return candidate, audit
    except (ValueError, RuntimeError, shapely.errors.GEOSException) as exc:
        audit.update(
            error=str(exc),
            failed_face=current_face.id
            if current_face and "reproducer" in audit
            else None,
            last_processed_face=current_face.id if current_face else None,
            seconds=time.perf_counter() - started,
            transaction="original complete shell retained",
        )
        return mesh, audit
