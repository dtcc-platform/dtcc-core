"""Minimal catalogue skins and bounded improvement on connected building groups."""

from dataclasses import dataclass, asdict, replace
import time

import numpy as np

from ...model import Field, Mesh, SemanticRegion
from ...common import info
from ...common.progress import report_progress, progress_range
from .compact_geometry import (
    building_domains,
    geometric_graph,
    face_coordinates,
    terrain_domains,
)
from .compact_triangulation import triangulate
from .face_validation import validate_face_geometry
from .wall_faces import _angle_constraint
from .roof_interfaces import RoofAssemblyError


def quality_by_kind(mesh, kinds):
    """Report shape, severe angles and cost separately for each geometric face kind."""
    xyz = mesh.vertices[mesh.faces]
    cross = np.cross(xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0])
    twice_area = np.linalg.norm(cross, axis=1)
    denominator = sum(
        np.sum((xyz[:, i] - xyz[:, (i + 1) % 3]) ** 2, axis=1) for i in range(3)
    )
    q = 2 * np.sqrt(3) * twice_area / denominator
    angles = np.min(
        [
            np.degrees(
                np.arctan2(
                    twice_area,
                    np.sum(
                        (xyz[:, (i + 1) % 3] - xyz[:, i])
                        * (xyz[:, (i + 2) % 3] - xyz[:, i]),
                        axis=1,
                    ),
                )
            )
            for i in range(3)
        ],
        axis=0,
    )
    ids = next(
        f.values for f in mesh.fields if f.name == "catalogue_surface_face"
    ).astype(int)
    labels = np.asarray(kinds)[ids]
    result = {}
    for kind in ("roof", "wall", "terrain"):
        selected = labels == kind
        values, a, area = q[selected], angles[selected], twice_area[selected] / 2
        result[kind] = dict(
            triangles=int(selected.sum()),
            minimum_q=float(values.min()) if len(values) else None,
            q_p01=float(np.quantile(values, 0.01)) if len(values) else None,
            q_p05=float(np.quantile(values, 0.05)) if len(values) else None,
            minimum_angle_deg=float(a.min()) if len(a) else None,
            defect_energy=float(np.sum((1 - values) ** 2)),
            **{
                f"below_q_{threshold:g}": int(np.sum(values < threshold))
                for threshold in (0.02, 0.1, 0.25)
            },
            **{
                f"below_{threshold}_degrees": int(np.sum(a < threshold))
                for threshold in (5, 10)
            },
            **{
                f"area_below_{threshold}_degrees_m2": float(area[a < threshold].sum())
                for threshold in (5, 10)
            },
        )
    return result


@dataclass(frozen=True)
class CompactSurfaceOptions:
    """Compact skin policy: size parameters affect sampled terrain only.

    Refinement, corner protection and building-size subdivision are disabled.
    Budgeted mode permits coordinated roof-edge points and at most two interior
    points per face, eight unique points per source and 32 per group. Each
    connected group has one fixed triangle ceiling of 1.5*T0. Invalid geometry
    or topology raises; it never becomes a successful stopping condition.
    detailed_audit stores candidate inputs/lineage without changing decisions.
    Coordinated boundary insertion requires the dtcc_mesher backend; other
    backends support zero-point and interior-only comparison modes.
    """

    budgeted: bool = True
    shared_edge_points: bool = True
    backend: str = "dtcc_mesher"
    max_face_vertices: int = 2
    max_source_vertices: int = 8
    max_group_vertices: int = 32
    detailed_audit: bool = False

    def __post_init__(self):
        if not isinstance(self.shared_edge_points, bool):
            raise TypeError("shared_edge_points must be boolean")
        if not isinstance(self.budgeted, bool):
            raise TypeError("budgeted must be boolean")
        if not isinstance(self.detailed_audit, bool):
            raise TypeError("detailed_audit must be boolean")
        if self.backend not in ("dtcc_mesher", "geos", "triangle"):
            raise ValueError("Compact backend must be dtcc_mesher, geos or triangle")
        if self.budgeted and self.shared_edge_points and self.backend != "dtcc_mesher":
            raise ValueError(
                'Budgeted shared-edge optimization requires backend="dtcc_mesher"; '
                "use budgeted=False or shared_edge_points=False for comparison backends"
            )
        for name, limit in [
            ("max_face_vertices", 2),
            ("max_source_vertices", 8),
            ("max_group_vertices", 32),
        ]:
            value = getattr(self, name)
            if (
                isinstance(value, bool)
                or not isinstance(value, int)
                or not 0 <= value <= limit
            ):
                raise ValueError(f"{name} must be an integer in [0, {limit}]")


def require_compact_mesher():
    """Check the actual loaded capability before fitting or terrain meshing."""
    from .dtcc_mesher_backend import _load_dtcc_mesher

    dm = _load_dtcc_mesher()
    message = (
        "Compact catalogue meshing requires native dtcc-mesher improve with "
        "bounded segment_trials. Rebuild the compatible local checkout; "
        "see docs/installation.md#compact-catalogue-meshing"
    )
    if not callable(getattr(dm, "improve", None)):
        raise RuntimeError(message)
    try:
        points = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
        triangles = np.array([[0, 1, 2]], dtype=np.uint32)
        segments = np.array([[0, 1], [1, 2], [2, 0]], dtype=np.uint32)
        probe = dm.improve(
            dm.Mesh(points, triangles, segments, summary=None), max_added_points=0
        )
        # A Python result default can mask an older loaded extension. Check
        # the native provenance payload as well as the public adapter.
        raw = dm._core._improve_raw(
            points, triangles, segments, 0, False, np.empty(0, dtype=np.uintp)
        )
    except (TypeError, ValueError, RuntimeError, AttributeError) as exc:
        raise RuntimeError(message) from exc
    if not hasattr(probe, "segment_trials") or "segment_trials" not in raw:
        raise RuntimeError(message)


def compact_surface(
    reference,
    patches,
    options,
    *,
    allow_source_steps=False,
    base_z=None,
    block_bases=None,
):
    """Mesh shared planar faces using a city attachment or an isolated base.

    The sampled terrain fixes attachment. With reference=None, base_z fixes the open bottom;
    analytic face coverage and the open-base topology contract govern the skin.
    """
    started = time.perf_counter()
    stage_seconds = {}
    if (reference is None) == (base_z is None):
        raise ValueError(
            "Provide either a reference city surface or an isolated base height"
        )
    if reference is not None:
        info("LOD2 mesh | constructing compact roof, wall and terrain faces")
    report_progress(percent=0, message="Constructing compact LOD2 faces...")
    owners = [p.source_key for p in patches]
    building = building_domains(
        patches,
        reference,
        allow_source_steps=allow_source_steps,
        base_z=base_z,
        block_bases=block_bases,
    )
    # Intrinsic complexity is measured before any terrain tessellation enters.
    intrinsic = [replace(f, rings=[], boundary=[]) for f in building]
    for f, source in zip(intrinsic, building):
        f._geometric_rings = [r.copy() for r in face_coordinates(source)]
    with progress_range(5, 10, "Building intrinsic LOD2 edge graph..."):
        _, intrinsic_audit = geometric_graph(intrinsic)
    report_progress(percent=10, message="Extracting terrain attachment faces...")
    if reference is None:
        terrain, terrain_locality = (
            [],
            dict(
                applicable=False,
                attachment_star_triangles=0,
                outside_star_triangles=0,
                outside_star_face_ids=[],
            ),
        )
    else:
        terrain, terrain_locality = terrain_domains(reference)
    original_terrain_rings = {
        id(f): [r.copy() for r in face_coordinates(f)] for f in terrain
    }
    faces = building + terrain
    stage_seconds["domains_and_intrinsic_graph"] = time.perf_counter() - started
    clock = time.perf_counter()
    with progress_range(20, 40, "Reconciling shared LOD2 edges..."):
        registry, graph = geometric_graph(faces)
    stage_seconds["shared_boundary_graph"] = time.perf_counter() - clock
    clock = time.perf_counter()
    vertices = np.asarray(registry.vertices)
    cache, rows = {}, []
    for face_index, f in enumerate(faces, 1):
        ids = sorted(set(x for e in f.boundary for x in e))
        xy = (np.asarray(vertices)[ids] - f.origin) @ f.basis.T
        from shapely.geometry import Polygon

        triangulation_polygon = Polygon(
            (np.asarray(vertices)[f.rings[0]] - f.origin) @ f.basis.T,
            [(np.asarray(vertices)[r] - f.origin) @ f.basis.T for r in f.rings[1:]],
        )
        inverse = {v: i for i, v in enumerate(ids)}
        segments = np.array([[inverse[a], inverse[b]] for a, b in f.boundary])
        unchanged = f.kind == "terrain" and all(
            len(a) == len(b)
            and np.max(np.min(np.linalg.norm(a[:, None] - b[None, :], axis=2), axis=1))
            <= 1e-8
            for a, b in zip(
                original_terrain_rings[id(f)],
                [np.asarray(vertices)[r] for r in f.rings],
            )
        )
        if unchanged:
            tri = reference.faces[f.original_faces]
            # Reuse the exact original vertices and triangles. Original IDs are
            # remapped once during assembly; no raster sampling is performed.
            points, triangles, flips = None, None, 0
        else:
            points = xy
            try:
                triangles, flips = triangulate(
                    points, segments, triangulation_polygon, backend=options.backend
                )
            except (ValueError, RuntimeError) as exc:
                blocks = {
                    patches[m % len(patches)].block_id
                    for m in f.markers
                    if 0 <= m < 2 * len(patches)
                }
                raise RoofAssemblyError(
                    f"Initial compact face {f.id} ({f.kind}, {f.source_key}): {exc}",
                    block_ids=blocks,
                ) from exc
        cache[f.id] = dict(
            ids=ids,
            points=points,
            triangles=triangles,
            segments=segments,
            unchanged=unchanged,
            optional=0,
            polygon=triangulation_polygon,
        )
        rows.append(
            dict(
                face_id=f.id,
                kind=f.kind,
                source_key=f.source_key,
                markers=f.markers,
                area_m2=float(f.polygon.area),
                boundary_vertices=len(ids),
                holes=len(f.rings) - 1,
                raw_boundary_vertices=graph["raw_face_boundary_vertices"][str(f.id)],
                removed_redundant_vertices=graph["raw_face_boundary_vertices"][
                    str(f.id)
                ]
                - len(ids),
                zero_steiner_triangles=(
                    len(f.original_faces) if unchanged else len(triangles)
                ),
                boundary_only_expected=len(ids) + 2 * (len(f.rings) - 1) - 2,
                fixed_terrain_triangles=bool(unchanged),
                legal_flips=flips,
                angle_constraint=_angle_constraint(f.polygon),
                optional_vertices=[],
            )
        )
        if face_index % 100 == 0 or face_index == len(faces):
            report_progress(
                percent=40 + 15 * face_index / len(faces),
                message=f"Triangulated LOD2 face {face_index}/{len(faces)}",
            )
    stage_seconds["initial_triangulation"] = time.perf_counter() - clock
    from copy import deepcopy

    initial_rows, initial_graph = rows, graph
    # Terrain attachment and geometric domains are immutable during optional
    # refinement. Reuse their verification while rechecking changed faces and
    # the complete shell topology for every accepted stage.
    face_validation = {}
    reference_map = {}
    from scipy.spatial import cKDTree

    tree = cKDTree(vertices)
    for f in terrain:
        if cache[f.id]["unchanged"]:
            for x in np.unique(reference.faces[f.original_faces]):
                if int(x) not in reference_map:
                    point = reference.vertices[x]
                    distance, index = tree.query(point)
                    if distance > 1e-8:
                        index = len(registry.vertices)
                        registry.vertices.append(point.copy())
                    reference_map[int(x)] = int(index)

    def assemble(optimization, *, final=False):
        vertices = np.asarray(registry.vertices)
        if final:
            rows, graph = deepcopy(initial_rows), deepcopy(initial_graph)
            for record in optimization["optional_vertices"]:
                for fi in record["affected_faces"]:
                    rows[fi]["optional_vertices"].append(record)
            additions = optimization["optional_vertices"]
            for record in additions:
                if record["kind"] == "segment":
                    edge = registry.edges[record["edge_id"]]
                    graph["retained_vertices"].append(
                        dict(
                            vertex=record["vertex"],
                            xyz=record["xyz"],
                            reason="optional quality improvement",
                            incident_faces=record["affected_faces"],
                            corner_faces=[],
                            fixed_terrain_faces=[],
                            neighbours=list(edge.endpoints),
                            collinear_retained_because=record["reason"],
                        )
                    )
            source_used = optimization["source_optional_vertices"]
        arrays, markers, face_ids = [], [], []
        shared_boundary_mismatches = 0
        for f in faces:
            c = cache[f.id]
            if c["unchanged"]:
                tri = np.array(
                    [
                        [reference_map[int(x)] for x in t]
                        for t in reference.faces[f.original_faces]
                    ]
                )
            else:
                tri = np.array(c["ids"])[c["triangles"]]
            previous = face_validation.get(f.id)
            if previous is not None and previous[0] is c:
                geometry = previous[1]
            else:
                from collections import Counter

                occurrences = Counter(
                    tuple(sorted((int(a), int(b))))
                    for t in tri
                    for a, b in zip(t, np.roll(t, -1))
                )
                actual = {edge for edge, count in occurrences.items() if count == 1}
                expected = {
                    tuple(sorted((int(a), int(b)))) for a, b in registry.boundary(f.id)
                }
                if actual != expected:
                    raise ValueError(
                        f"Face {f.id} does not consume its canonical shared boundary"
                    )
                geometry = validate_face_geometry(f, vertices, tri)
                face_validation[f.id] = (c, geometry)
            if final:
                rows[f.id].update(
                    triangles=len(tri),
                    geometry=geometry,
                    boundary_vertices=len(set(x for r in f.rings for x in r)),
                    rings_xyz=[np.asarray(vertices)[r].tolist() for r in f.rings],
                )
            arrays.append(tri)
            markers.extend([f.markers[0]] * len(tri))
            face_ids.extend([f.id] * len(tri))
        triangles = np.concatenate(arrays)
        used, inverse = np.unique(triangles, return_inverse=True)
        result = Mesh(
            vertices=np.asarray(vertices)[used],
            faces=inverse.reshape(-1, 3),
            markers=np.array(markers),
        )
        result.fields.append(
            Field(
                name="catalogue_surface_face",
                values=np.array(face_ids),
                association="face",
            )
        )
        if final:
            face_ids_array = np.asarray(face_ids)
            for f in faces:
                result.regions.append(
                    SemanticRegion(
                        "urn:dtcc:surface-face",
                        id=str(f.id),
                        indices=np.flatnonzero(face_ids_array == f.id),
                        attributes=dict(
                            kind=f.kind,
                            source_key=f.source_key,
                            contributing_patch_markers=f.markers,
                        ),
                    )
                )
        from .roof_validation import validate_shell
        from ..geometry_builders.meshes import _shell_closure_audit

        topology = validate_shell(
            result,
            num_patches=len(patches),
            closure_audit=_shell_closure_audit,
            base_z=base_z,
        )
        if not topology["valid"]:
            raise ValueError(f"Compact shell topology failed: {topology}")
        wall_errors, ownership_areas = [], []
        if reference is not None:

            def triangle_keys(points, triangles):
                return {
                    tuple(sorted(tuple(map(float, p)) for p in points[t]))
                    for t in triangles
                }

            outside = reference.faces[terrain_locality["outside_star_face_ids"]]
            fixed_keys = triangle_keys(reference.vertices, outside)
            final_keys = triangle_keys(
                result.vertices, result.faces[result.markers < 0]
            )
            if not fixed_keys <= final_keys:
                raise ValueError(
                    "Compact attachment changed terrain outside the contact star"
                )
            terrain_locality["outside_star_bitwise_preserved"] = True
        if not final:
            return result, dict(topology=topology)
        graph["face_edge_sequences"] = {
            str(f.id): [
                dict(edge_id=int(i), forward=bool(forward))
                for i, forward in registry.face_edges[f.id]
            ]
            for f in faces
        }
        graph["edges"] = [
            dict(
                edge_id=i,
                endpoints=list(e.endpoints),
                incident_faces=sorted(e.incidents),
                parameters=e.parameters,
                vertices=e.vertices,
                xyz=vertices[e.vertices].tolist(),
            )
            for i, e in enumerate(registry.edges)
        ]
        intrinsic_triangles = sum(
            len(set(x for r in f.rings for x in r)) + 2 * (len(f.rings) - 1) - 2
            for f in intrinsic
        )
        attached = sum(
            r["zero_steiner_triangles"] for r in rows if r["kind"] != "terrain"
        )
        return result, dict(
            method="minimal_catalogue_geometry",
            accepted=True,
            options=asdict(options),
            skin_size_policy="boundary only; public mesh sizes apply to sampled terrain, not compact skin",
            geometry_tolerance_m=1e-7,
            shared_edge_tolerance_m=1e-8,
            intrinsic=dict(
                faces=len(intrinsic),
                vertices=len(intrinsic_audit["retained_vertices"]),
                triangles=intrinsic_triangles,
                graph=intrinsic_audit,
            ),
            graph=graph,
            faces=rows,
            analytic_wall_areas=[
                dict(
                    source_key=source,
                    area_m2=sum(
                        f.polygon.area
                        for f in building
                        if f.kind == "wall" and f.source_key == source
                    ),
                )
                for source in sorted(set(owners))
            ],
            optimization=optimization,
            optional_vertices=additions,
            source_optional_vertices=source_used,
            terrain_attachment_triangle_overhead=attached - intrinsic_triangles,
            zero_steiner_building_triangles=attached,
            terrain_locality=terrain_locality,
            terrain_triangles=sum(
                r["triangles"] for r in rows if r["kind"] == "terrain"
            ),
            unchanged_terrain_triangles=sum(
                len(f.original_faces) for f in terrain if cache[f.id]["unchanged"]
            ),
            topology=topology,
            shared_boundary_mismatches=shared_boundary_mismatches,
            seconds=time.perf_counter() - started,
        )

    from .compact_improvement import building_groups, improve_shared

    groups = building_groups(registry, faces)
    optimization = dict(
        optional_vertices=[],
        source_optional_vertices=dict.fromkeys(sorted(set(owners)), 0),
        groups=[],
        stages=[],
    )
    report_progress(percent=55, message="Validating compact LOD2 shell...")
    clock = time.perf_counter()
    result, audit = assemble(optimization)
    stage_seconds["minimal_validation"] = time.perf_counter() - clock
    quality_before = quality_by_kind(result, [f.kind for f in faces])
    # Invalid minimal geometry raises above; optional stages cannot convert it
    # into a successful reconstruction or trigger a roof-fitting fallback.
    if reference is not None:
        info(
            f"LOD2 mesh | improving {len(groups)} building groups ({len(faces)} faces)"
        )
    clock = last_report = time.perf_counter()
    for group_index, group in enumerate(groups, 1):
        report_progress(
            percent=60 + 35 * (group_index - 1) / len(groups),
            message=f"Improving LOD2 building group {group_index}/{len(groups)}...",
        )
        ids = group["face_ids"]
        t0 = sum(len(cache[i]["triangles"]) for i in ids)
        ceiling = 3 * t0 // 2
        registry, measured = improve_shared(
            registry,
            cache,
            faces,
            options,
            max_building_triangles=ceiling,
            face_ids=set(ids),
        )
        for point in measured["optional_vertices"]:
            point["group_id"] = group["sources"][0]
        optimization["optional_vertices"].extend(measured["optional_vertices"])
        for source, count in measured["source_optional_vertices"].items():
            optimization["source_optional_vertices"][source] += count
        optimization["groups"].append(
            dict(
                group_id=group["sources"][0],
                **group,
                T0=t0,
                overall_triangle_ceiling=ceiling,
                end_triangles=measured["building_triangles"],
                interior_points=sum(
                    p["kind"] == "interior" for p in measured["optional_vertices"]
                ),
                shared_points=sum(
                    p["kind"] == "segment" for p in measured["optional_vertices"]
                ),
                source_charges=measured["source_optional_vertices"],
                termination=measured["termination"],
                stopping=measured["stopping"],
            )
        )
        optimization["stages"].append(dict(group_id=group["sources"][0], **measured))
        report_progress(
            percent=60 + 35 * group_index / len(groups),
            message=f"Improved LOD2 building group {group_index}/{len(groups)}",
        )
        if reference is not None and (
            time.perf_counter() - last_report >= 5 or group_index == len(groups)
        ):
            info(
                f"LOD2 mesh | improved {group_index}/{len(groups)} building groups "
                f"({time.perf_counter() - started:.1f} s)"
            )
            last_report = time.perf_counter()
    report_progress(percent=95, message="Finalizing compact LOD2 mesh...")
    stage_seconds["bounded_refinement"] = time.perf_counter() - clock
    clock = time.perf_counter()
    result, audit = assemble(optimization, final=True)
    stage_seconds["final_validation"] = time.perf_counter() - clock
    audit.update(
        quality_before=quality_before,
        quality_after=quality_by_kind(result, [f.kind for f in faces]),
        stage_seconds=stage_seconds,
    )
    if reference is not None:
        info(
            f"LOD2 mesh | compact surface ready: {len(result.faces)} triangles "
            f"({time.perf_counter() - started:.1f} s)"
        )
    audit["T0"] = sum(g["T0"] for g in optimization["groups"])
    audit["building_triangles"] = sum(
        r["triangles"] for r in audit["faces"] if r["kind"] != "terrain"
    )
    audit["optional_counts"] = {
        kind: sum(p["kind"] == kind for p in optimization["optional_vertices"])
        for kind in ("interior", "segment")
    }
    if not options.detailed_audit:
        audit["intrinsic"].pop("graph", None)
        audit["graph"] = {
            key: value
            for key, value in audit["graph"].items()
            if key
            in (
                "raw_vertices",
                "removed_redundant_vertices",
                "raw_faces",
                "raw_face_boundary_vertices",
            )
        }
        for row in audit["faces"]:
            row.pop("rings_xyz", None)
    report_progress(percent=100, message="Compact LOD2 mesh complete")
    return result, audit
