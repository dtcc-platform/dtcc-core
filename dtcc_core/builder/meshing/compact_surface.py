"""Minimal catalogue skins and bounded improvement on connected building groups."""
from dataclasses import dataclass, asdict, replace
import time

import numpy as np

from ...model import Field, Mesh, SemanticRegion
from .compact_geometry import building_domains, geometric_graph, extract_surface_faces, face_coordinates, terrain_domains
from .compact_triangulation import triangulate
from .surface_refinement import _validate_geometry
from .wall_faces import _angle_constraint


@dataclass(frozen=True)
class CompactSurfaceOptions:
    """Compact skin policy: size parameters affect reference terrain only.

    Refinement, corner protection and building-size subdivision are disabled.
    Budgeted mode permits coordinated roof-edge points and at most two interior
    points per face, eight unique points per source and 32 per group in Stage A.
    Stage B defaults to four extra shared points/source and eight/group;
    extension caps may be doubled explicitly. adaptive_extension opts into
    guarded continuation of standard 4/8 caps to hard Stage B totals of 8/16;
    other explicit caps retain fixed semantics. Both stages obey group triangle
    ceilings of 1.5*T0 overall and 0.20*TB extra in Stage B.
    detailed_audit stores candidate inputs/lineage without changing decisions.
    Coordinated boundary insertion requires the dtcc_mesher backend; other
    backends support zero-point and interior-only comparison modes.
    """
    budgeted: bool = True
    shared_edge_points: bool = True
    backend: str = 'dtcc_mesher'
    max_face_vertices: int = 2
    max_source_vertices: int = 8
    max_group_vertices: int = 32
    extension_source_vertices: int = 4
    extension_group_vertices: int = 8
    detailed_audit: bool = False
    adaptive_extension: bool = False

    def __post_init__(self):
        if not isinstance(self.shared_edge_points, bool):
            raise TypeError("shared_edge_points must be boolean")
        if not isinstance(self.budgeted, bool):
            raise TypeError('budgeted must be boolean')
        if not isinstance(self.adaptive_extension, bool):
            raise TypeError('adaptive_extension must be boolean')
        if not isinstance(self.detailed_audit, bool):
            raise TypeError('detailed_audit must be boolean')
        if self.backend not in ('dtcc_mesher', 'geos', 'triangle'):
            raise ValueError('Compact backend must be dtcc_mesher, geos or triangle')
        if self.budgeted and self.shared_edge_points and self.backend != 'dtcc_mesher':
            raise ValueError('Budgeted shared-edge optimization requires backend="dtcc_mesher"; '
                             'use budgeted=False or shared_edge_points=False for comparison backends')
        for name, limit in [('max_face_vertices',2), ('max_source_vertices',8), ('max_group_vertices',32),
                            ('extension_source_vertices',8),('extension_group_vertices',16)]:
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or not 0 <= value <= limit:
                raise ValueError(f'{name} must be an integer in [0, {limit}]')


def require_compact_mesher():
    """Check the actual loaded capability before fitting or terrain meshing."""
    from .dtcc_mesher_backend import _load_dtcc_mesher
    dm = _load_dtcc_mesher()
    message = ('Compact catalogue meshing requires native dtcc-mesher improve with '
               'bounded segment_trials. Rebuild the compatible local checkout; '
               'see docs/installation.md#compact-catalogue-meshing')
    if not callable(getattr(dm,'improve',None)):
        raise RuntimeError(message)
    try:
        points=np.array([[0.,0.],[1.,0.],[0.,1.]])
        triangles=np.array([[0,1,2]],dtype=np.uint32)
        segments=np.array([[0,1],[1,2],[2,0]],dtype=np.uint32)
        probe=dm.improve(dm.Mesh(points,triangles,segments,summary=None),max_added_points=0)
        # A Python result default can mask an older loaded extension. Check
        # the native provenance payload as well as the public adapter.
        raw=dm._core._improve_raw(points,triangles,segments,0,False,np.empty(0,dtype=np.uintp))
    except (TypeError,ValueError,RuntimeError,AttributeError) as exc:
        raise RuntimeError(message) from exc
    if not hasattr(probe,'segment_trials') or 'segment_trials' not in raw:
        raise RuntimeError(message)


def compact_surface(reference, patches, options, *, allow_source_steps=False, base_z=None):
    """Mesh shared planar faces using a city attachment or an isolated base.

    A reference surface fixes terrain attachment and supplies an independent
    exposed-wall check. With reference=None, base_z fixes the open bottom;
    analytic face coverage and the open-base topology contract govern the skin.
    """
    started = time.perf_counter()
    if (reference is None) == (base_z is None):
        raise ValueError('Provide either a reference city surface or an isolated base height')
    owners = [p.source_key for p in patches]
    building = building_domains(patches, reference, allow_source_steps=allow_source_steps, base_z=base_z)
    # Intrinsic complexity is measured before any terrain tessellation enters.
    intrinsic = [replace(f, rings=[], boundary=[]) for f in building]
    for f, source in zip(intrinsic, building):
        f._geometric_rings = [r.copy() for r in face_coordinates(source)]
    _, intrinsic_audit = geometric_graph(intrinsic)
    if reference is None:
        terrain, terrain_locality = [], dict(applicable=False, attachment_star_triangles=0,
            outside_star_triangles=0, outside_star_face_ids=[])
    else:
        old_faces = extract_surface_faces(reference, owners)
        terrain, terrain_locality = terrain_domains(reference, old_faces)
    original_terrain_rings = {id(f):[r.copy() for r in face_coordinates(f)] for f in terrain}
    faces = building+terrain
    registry, graph = geometric_graph(faces)
    vertices = registry.vertices
    cache, rows = {}, []
    for f in faces:
        ids = sorted(set(x for e in f.boundary for x in e))
        xy = (np.asarray(vertices)[ids]-f.origin) @ f.basis.T
        from shapely.geometry import Polygon
        triangulation_polygon = Polygon(
            (np.asarray(vertices)[f.rings[0]]-f.origin) @ f.basis.T,
            [(np.asarray(vertices)[r]-f.origin) @ f.basis.T for r in f.rings[1:]])
        inverse = {v:i for i,v in enumerate(ids)}
        segments = np.array([[inverse[a],inverse[b]] for a,b in f.boundary])
        unchanged = f.kind == 'terrain' and all(
            len(a)==len(b) and np.max(np.min(np.linalg.norm(a[:,None]-b[None,:], axis=2), axis=1)) <= 1e-8
            for a,b in zip(original_terrain_rings[id(f)],
                           [np.asarray(vertices)[r] for r in f.rings]))
        if unchanged:
            tri = reference.faces[f.original_faces]
            # Reuse the exact original vertices and triangles. Original IDs are
            # remapped once during assembly; no raster sampling is performed.
            points, triangles, flips = None, None, 0
        else:
            points = xy
            try:
                triangles, flips = triangulate(points, segments, triangulation_polygon, backend=options.backend)
            except (ValueError,RuntimeError) as exc:
                raise ValueError(f'Initial compact face {f.id} ({f.kind}, {f.source_key}): {exc}') from exc
        cache[f.id] = dict(ids=ids, points=points, triangles=triangles, segments=segments,
                           unchanged=unchanged, optional=0, polygon=triangulation_polygon)
        rows.append(dict(face_id=f.id,kind=f.kind,source_key=f.source_key,markers=f.markers,
                         area_m2=float(f.polygon.area),boundary_vertices=len(ids),holes=len(f.rings)-1,
                         raw_boundary_vertices=graph['raw_face_boundary_vertices'][str(f.id)],
                         removed_redundant_vertices=graph['raw_face_boundary_vertices'][str(f.id)]-len(ids),
                         zero_steiner_triangles=(len(f.original_faces) if unchanged else len(triangles)),
                         boundary_only_expected=len(ids)+2*(len(f.rings)-1)-2,
                         fixed_terrain_triangles=bool(unchanged),legal_flips=flips,
                         angle_constraint=_angle_constraint(f.polygon),optional_vertices=[]))
    from copy import deepcopy
    initial_rows, initial_graph = rows, graph

    def assemble(optimization):
        vertices = registry.vertices
        rows, graph = deepcopy(initial_rows), deepcopy(initial_graph)
        for record in optimization['optional_vertices']:
            for fi in record['affected_faces']:
                rows[fi]['optional_vertices'].append(record)
        additions = optimization['optional_vertices']
        for record in additions:
            if record['kind']=='segment':
                edge=registry.edges[record['edge_id']]
                graph['retained_vertices'].append(dict(vertex=record['vertex'],xyz=record['xyz'],
                    reason='optional quality improvement',incident_faces=record['affected_faces'],
                    corner_faces=[],fixed_terrain_faces=[],neighbours=list(edge.endpoints),
                    collinear_retained_because=record['reason']))
        source_used = optimization['source_optional_vertices']
        arrays, markers, face_ids = [], [], []
        shared_boundary_mismatches = 0
        reference_map = {}
        # Canonicalize unchanged terrain boundary IDs through the geometric graph.
        from scipy.spatial import cKDTree
        tree = cKDTree(np.asarray(vertices))
        for f in faces:
            c = cache[f.id]
            if c['unchanged']:
                tri = []
                for t in reference.faces[f.original_faces]:
                    output = []
                    for x in t:
                        if int(x) not in reference_map:
                            p = reference.vertices[x]
                            distance, index = tree.query(p)
                            if distance <= 1e-8:
                                reference_map[int(x)] = int(index)
                            else:
                                reference_map[int(x)] = len(vertices)
                                vertices.append(p.copy())
                        output.append(reference_map[int(x)])
                    tri.append(output)
                tri = np.array(tri, dtype=int)
            else:
                tri = np.array(c['ids'])[c['triangles']]
            from collections import Counter
            occurrences = Counter(tuple(sorted((int(a),int(b)))) for t in tri
                                  for a,b in zip(t,np.roll(t,-1)))
            actual = {edge for edge,count in occurrences.items() if count==1}
            expected = {tuple(sorted((int(a),int(b)))) for a,b in registry.boundary(f.id)}
            if actual != expected:
                shared_boundary_mismatches += 1
                raise ValueError(f'Face {f.id} does not consume its canonical shared boundary')
            geometry = _validate_geometry(f, np.asarray(vertices), tri)
            rows[f.id].update(triangles=len(tri),geometry=geometry,
                              boundary_vertices=len(set(x for r in f.rings for x in r)),
                              rings_xyz=[np.asarray(vertices)[r].tolist() for r in f.rings])
            arrays.append(tri)
            markers.extend([f.markers[0]]*len(tri))
            face_ids.extend([f.id]*len(tri))
        triangles = np.concatenate(arrays)
        used, inverse = np.unique(triangles, return_inverse=True)
        result = Mesh(vertices=np.asarray(vertices)[used], faces=inverse.reshape(-1,3), markers=np.array(markers))
        result.fields.append(Field(name='catalogue_surface_face',values=np.array(face_ids),association='face'))
        for f in faces:
            result.regions.append(SemanticRegion('urn:dtcc:surface-face',id=str(f.id),
                indices=np.flatnonzero(np.array(face_ids)==f.id), attributes=dict(kind=f.kind,
                    source_key=f.source_key,contributing_patch_markers=f.markers)))
        from .roof_validation import validate_shell
        from ..geometry_builders.meshes import _shell_closure_audit
        topology = validate_shell(result,num_patches=len(patches),closure_audit=_shell_closure_audit,base_z=base_z)
        if not topology['valid']:
            raise ValueError(f'Compact shell topology failed: {topology}')
        wall_errors, ownership_areas = [], []
        if reference is not None:
            def triangle_keys(points, triangles):
                return {tuple(sorted(tuple(map(float,p)) for p in points[t])) for t in triangles}
            outside = reference.faces[terrain_locality['outside_star_face_ids']]
            fixed_keys = triangle_keys(reference.vertices, outside)
            final_keys = triangle_keys(result.vertices, result.faces[result.markers<0])
            if not fixed_keys <= final_keys:
                raise ValueError('Compact attachment changed terrain outside the contact star')
            terrain_locality['outside_star_bitwise_preserved'] = True
            # Prove exposed wall coverage against the established native interpretation.
            wall_errors = []
            import shapely
            for f in building:
                if f.kind != 'wall':
                    continue
                wall_ids = np.flatnonzero((reference.markers >= 0) & (reference.markers < len(owners)))
                wall_ids = np.array([i for i in wall_ids if owners[int(reference.markers[i])]==f.source_key], dtype=int)
                xyz = reference.vertices[reference.faces[wall_ids]]
                residuals = np.max(abs((xyz-f.origin) @ f.normal), axis=1)
                normals = np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0])
                normals /= np.linalg.norm(normals,axis=1)[:,None]
                selected = (residuals <= 1e-7) & (normals @ f.normal > 1-1e-8)
                projected = (reference.vertices-f.origin) @ f.basis.T
                polygons = shapely.polygons(projected[reference.faces[wall_ids[selected]]])
                union = shapely.union_all(polygons)
                error = f.polygon.difference(union).area
                if error > max(1e-7,f.polygon.area*1e-9):
                    raise ValueError(f'Compact wall differs from exposed native wall by {error} m2')
                wall_errors.append(float(error))
            ownership_areas = []
            xyz = reference.vertices[reference.faces]
            native_area = np.linalg.norm(np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]),axis=1)/2
            for owner in sorted(set(owners)):
                native_wall = sum(native_area[i] for i,m in enumerate(reference.markers)
                                  if 0 <= m < len(owners) and owners[int(m)]==owner)
                direct_wall = sum(f.polygon.area for f in building if f.kind=='wall' and f.source_key==owner)
                error = abs(native_wall-direct_wall)
                if error > max(1e-7,native_wall*1e-9):
                    raise ValueError(f'Exposed wall area/ownership mismatch for {owner}: {error}')
                ownership_areas.append(dict(source_key=owner,native_wall_area_m2=float(native_wall),
                                           direct_wall_area_m2=float(direct_wall),error_m2=float(error)))
        graph['face_edge_sequences'] = {str(f.id):[dict(edge_id=int(i),forward=bool(forward))
            for i,forward in registry.face_edges[f.id]] for f in faces}
        graph['edges'] = [dict(edge_id=i,endpoints=list(e.endpoints),incident_faces=sorted(e.incidents),
                               parameters=e.parameters,vertices=e.vertices,
                               xyz=np.asarray(registry.vertices)[e.vertices].tolist()) for i,e in enumerate(registry.edges)]
        intrinsic_triangles = sum(len(set(x for r in f.rings for x in r))+2*(len(f.rings)-1)-2 for f in intrinsic)
        attached = sum(r['zero_steiner_triangles'] for r in rows if r['kind']!='terrain')
        return result, dict(method='minimal_catalogue_geometry',accepted=True,options=asdict(options),
            skin_size_policy='boundary only; public mesh sizes apply to reference terrain, not compact skin',
            geometry_tolerance_m=1e-7, shared_edge_tolerance_m=1e-8,
            intrinsic=dict(faces=len(intrinsic),vertices=len(intrinsic_audit['retained_vertices']),
                           triangles=intrinsic_triangles,graph=intrinsic_audit),
            graph=graph,faces=rows,ownership_areas=ownership_areas,optimization=optimization,optional_vertices=additions,source_optional_vertices=source_used,
            terrain_attachment_triangle_overhead=attached-intrinsic_triangles,
            zero_steiner_building_triangles=attached,
            terrain_locality=terrain_locality,
            terrain_triangles=sum(r['triangles'] for r in rows if r['kind']=='terrain'),
            unchanged_terrain_triangles=sum(len(f.original_faces) for f in terrain if cache[f.id]['unchanged']),
            wall_reference_coverage_max_m2=max(wall_errors,default=0.) if reference is not None else None,
            topology=topology,shared_boundary_mismatches=shared_boundary_mismatches,seconds=time.perf_counter()-started)

    from .compact_improvement import building_groups, improve_shared
    groups = building_groups(registry, faces)
    optimization = dict(optional_vertices=[], source_optional_vertices=dict.fromkeys(sorted(set(owners)),0),
                        groups=[], stages=[])
    result, audit = assemble(optimization)
    # Invalid minimal geometry raises above; optional stages cannot convert it
    # into a successful reconstruction or trigger a roof-fitting fallback.
    for group in groups:
        ids=group['face_ids']
        t0=sum(len(cache[i]['triangles']) for i in ids)
        overall=3*t0//2
        group_audit=dict(group_id=group['sources'][0],**group,T0=t0,overall_triangle_ceiling=overall)
        optimization['groups'].append(group_audit)
        for stage in ('A','B'):
            tb=sum(len(cache[i]['triangles']) for i in ids)
            policy=(options if stage=='A' else replace(options,
                budgeted=options.budgeted and options.shared_edge_points,max_face_vertices=0,
                max_source_vertices=options.extension_source_vertices,max_group_vertices=options.extension_group_vertices))
            ceiling=overall if stage=='A' else min(overall,tb+tb//5)
            enabled=policy.budgeted and policy.max_group_vertices>0 and policy.max_source_vertices>0
            before=deepcopy((registry,cache,faces)) if enabled else None
            standard_checkpoint={}
            measured=None

            def validated_trial(measured):
                for point in measured['optional_vertices']:
                    point.update(stage=stage,group_id=group_audit['group_id'])
                trial=deepcopy(optimization)
                trial['optional_vertices'].extend(measured['optional_vertices'])
                for source,count in measured['source_optional_vertices'].items():
                    trial['source_optional_vertices'][source]+=count
                candidate,candidate_audit=assemble(trial)
                return trial,candidate,candidate_audit

            try:
                registry, measured=improve_shared(registry,cache,faces,policy,
                    max_building_triangles=ceiling,face_ids=set(ids),
                    adaptive=(stage=='B' and policy.budgeted and options.adaptive_extension and
                              options.extension_source_vertices==4 and options.extension_group_vertices==8),
                    _standard_checkpoint=standard_checkpoint)
                if enabled:
                    trial,candidate,candidate_audit=validated_trial(measured)
            except (ValueError,RuntimeError) as exc:
                if not enabled:
                    raise
                recovered=False
                if standard_checkpoint:
                    registry=standard_checkpoint['registry'];cache=standard_checkpoint['caches']
                    faces[:]=standard_checkpoint['faces']
                    measured=dict(measured or {},optional_vertices=standard_checkpoint['optional_vertices'],
                        source_optional_vertices=standard_checkpoint['source_optional_vertices'],
                        standard_result=standard_checkpoint['standard_result'],escalations=[],
                        effective_allowances=dict(source=dict.fromkeys(group['sources'],policy.max_source_vertices),
                                                  group=policy.max_group_vertices),
                        building_triangles=standard_checkpoint['standard_result']['building_triangles'],
                        termination='unexpected_numerical_failure',error=str(exc),adaptive_rollback=True,
                        stopping=dict(primary_reason='unexpected_numerical_failure',
                            binding_limits=['unexpected_numerical_failure'],evaluated_finite_set=False,
                            rejection_counts={'whole shell validation':1}))
                    try:
                        trial,candidate,candidate_audit=validated_trial(measured)
                    except (ValueError,RuntimeError):
                        pass
                    else:
                        optimization=trial;result,audit=candidate,candidate_audit;recovered=True
                if not recovered:
                    registry,cache,restored=before
                    faces[:]=restored
                    measured=dict(optional_vertices=[],source_optional_vertices=dict.fromkeys(group['sources'],0),
                        termination='unexpected_numerical_failure',error=str(exc),
                        stopping=dict(primary_reason='unexpected_numerical_failure',
                            binding_limits=['unexpected_numerical_failure'],evaluated_finite_set=False,
                            rejection_counts={'whole shell validation':1}))
            else:
                if enabled:
                    optimization=trial;result,audit=candidate,candidate_audit
            group_audit[stage]=dict(start_triangles=tb,triangle_ceiling=ceiling,
                end_triangles=sum(len(cache[i]['triangles']) for i in ids),
                interior_points=sum(p['kind']=='interior' for p in measured['optional_vertices']),
                shared_points=sum(p['kind']=='segment' for p in measured['optional_vertices']),
                source_charges=measured['source_optional_vertices'],termination=measured['termination'],
                stopping=measured.get('stopping'),standard_result=measured.get('standard_result'),
                effective_allowances=measured.get('effective_allowances'),escalations=measured.get('escalations',[]))
            # Audit containers are independent of transactional state copies.
            optimization['groups']=[group_audit if g['group_id']==group_audit['group_id'] else g for g in optimization['groups']]
            optimization['stages'].append(dict(group_id=group_audit['group_id'],stage=stage,**measured))
    result,audit=assemble(optimization)
    audit['T0']=sum(g['T0'] for g in optimization['groups'])
    audit['TB']=sum(g.get('A',{}).get('end_triangles',g['T0']) for g in optimization['groups'])
    audit['building_triangles']=sum(r['triangles'] for r in audit['faces'] if r['kind']!='terrain')
    audit['optional_counts']={stage:{kind:sum(p['stage']==stage and p['kind']==kind for p in optimization['optional_vertices'])
                                   for kind in ('interior','segment')} for stage in ('A','B')}
    if not options.detailed_audit:
        audit['intrinsic'].pop('graph',None)
        audit['graph']={key:value for key,value in audit['graph'].items() if key in
            ('raw_vertices','removed_redundant_vertices','raw_faces','raw_face_boundary_vertices')}
        for row in audit['faces']:
            row.pop('rings_xyz',None)
    return result,audit
