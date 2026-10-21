"""Bounded source-budgeted transactions on the compact 3D face graph."""
from collections import Counter, defaultdict
from copy import deepcopy
from types import SimpleNamespace
import time

import numpy as np
from shapely.geometry import Polygon

from .compact_triangulation import improve, triangulate
from .surface_faces import boundary_rings
from .face_validation import validate_face_geometry

MAX_NATIVE_REQUESTS = 512
MAX_JOINT_TRIALS = 256
MAX_GROUP_SECONDS = 60.


def building_groups(registry, faces):
    """Connected source domains and shared building faces, excluding terrain."""
    sources=sorted({f.source_key for f in faces if f.kind!='terrain'})
    parent={s:s for s in sources}
    def root(s):
        while parent[s]!=s:
            parent[s]=parent[parent[s]];s=parent[s]
        return s
    def join(a,b):
        a,b=root(a),root(b)
        if a!=b:parent[max(a,b)]=min(a,b)
    for edge in registry.edges:
        owned=sorted({faces[i].source_key for i in edge.incidents if faces[i].kind!='terrain'})
        for source in owned[1:]:join(owned[0],source)
    # Unequal-height neighbours can have no common roof-height edge. Their
    # cleaned roof projections still share a source-cell boundary.
    import shapely
    domains=[]
    for source in sources:
        polygons=[Polygon(np.asarray(registry.vertices)[f.rings[0],:2],
                          [np.asarray(registry.vertices)[r,:2] for r in f.rings[1:]])
                  for f in faces if f.kind=='roof' and f.source_key==source]
        domains.append(shapely.union_all(polygons))
    tree=shapely.STRtree(domains)
    for i,domain in enumerate(domains):
        for j in sorted(tree.query(domain,predicate='intersects')):
            if j<i and domain.intersection(domains[j]).length>1e-8:
                join(sources[i],sources[j])
    groups=defaultdict(list)
    for f in faces:
        if f.kind!='terrain':groups[root(f.source_key)].append(f.id)
    return [dict(sources=sorted({faces[i].source_key for i in ids}),face_ids=sorted(ids))
            for _,ids in sorted(groups.items())]


def quality(points, triangles, segments):
    """Use native returned metrics for application-level joint acceptance."""
    measured = improve(points, triangles, segments, max_added_points=0)
    return dict(measured.before, min_angle=measured.quality_before.min_angle_deg_min,
                poor_20_count=measured.quality_before.count_min_angle_lt_20)


def aggregate(rows):
    """Combine face metrics using worst quality and summed defect measures."""
    return {key:(min(r[key] for r in rows) if key in ('minimum_q','min_angle')
                 else sum(r[key] for r in rows)) for key in rows[0]}


def joint_acceptance(before, after, triangle_cost):
    """Require measured combined gain and guard every incident face's defects."""
    a,b = aggregate(before),aggregate(after)
    if triangle_cost <= 0:
        return False,'nonpositive triangle cost',a,b
    for old,new in zip(before,after):
        if new['minimum_q'] < old['minimum_q']-1e-10:
            return False,'incident face worst q regressed',a,b
        for angle in (5,10):
            if (new[f'poor_{angle}_count']>old[f'poor_{angle}_count'] or
                new[f'poor_{angle}_area']>old[f'poor_{angle}_area']+1e-10*max(1.,old[f'poor_{angle}_area'])):
                return False,f'incident face {angle}-degree defects regressed',a,b
    if a['penalty']-b['penalty'] <= .015*triangle_cost:
        return False,'insufficient combined gain per total triangle',a,b
    return True,'combined defect reduction per total triangle cost',a,b


def permitted_segments(face, cache, registry, faces):
    """Only non-terrain roof/roof or roof/wall intervals may split."""
    eligible = {}
    for ei,_ in registry.face_edges[face.id]:
        e = registry.edges[ei]
        incident = [faces[i] for i in sorted(e.incidents)]
        if len(incident)!=2 or not any(f.kind=='roof' for f in incident) or any(f.kind=='terrain' for f in incident):
            continue
        for a,b in zip(e.vertices,e.vertices[1:]):
            eligible[tuple(sorted((a,b)))]=ei
    return {i:eligible[key] for i,(a,b) in enumerate(cache['segments'])
            if (key:=tuple(sorted((cache['ids'][int(a)],cache['ids'][int(b)])))) in eligible}


def _retessellate(face, cache, registry, backend):
    boundary = registry.boundary(face.id)
    ids = cache['ids'].copy()
    required = sorted(set(v for edge in boundary for v in edge)-set(ids))
    ids.extend(required)
    inverse = {v:i for i,v in enumerate(ids)}
    coordinates = np.asarray([registry.vertices[i] for i in ids])
    points = (coordinates-face.origin) @ face.basis.T
    segments = np.array([[inverse[a],inverse[b]] for a,b in boundary],dtype=int)
    rings = boundary_rings(boundary)
    rings.sort(key=lambda r:-Polygon(points[[inverse[v] for v in r]]).area)
    local = [points[[inverse[v] for v in r]] for r in rings]
    polygon = Polygon(local[0],local[1:])
    triangles,flips = triangulate(points,segments,polygon,backend=backend)
    validate_face_geometry(face,coordinates,triangles)
    return dict(cache,ids=ids,points=points,triangles=triangles,segments=segments,
                polygon=polygon,rings=rings,flips=flips)


def edge_transaction(registry, caches, faces, proposer, result, segment_map, *, backend):
    """Stage every incident face on a private graph; no mutation before acceptance."""
    c = caches[proposer.id]
    parent = int(result.parent_segment)
    ei = segment_map[parent]
    e = registry.edges[ei]
    ia,ib = map(int,c['segments'][parent])
    ga,gb = c['ids'][ia],c['ids'][ib]
    ta,tb = e.parameters[e.vertices.index(ga)],e.parameters[e.vertices.index(gb)]
    t = ta+float(result.segment_parameter)*(tb-ta)
    staged = registry.copy_for_split(ei)
    affected,added = staged.reconcile([(ei,t)],max_vertices=len(staged.vertices)+1)
    if added != 1:
        return None,dict(accepted=False,reason='duplicate canonical split',edge_id=ei,parameter=t)
    updates = {}
    for fi in sorted(affected):
        if faces[fi].kind=='terrain' or caches[fi]['unchanged']:
            return None,dict(accepted=False,reason='fixed terrain edge',edge_id=ei,parameter=t)
        updates[fi] = _retessellate(faces[fi],caches[fi],staged,backend)
    before = [quality(caches[i]['points'],caches[i]['triangles'],caches[i]['segments']) for i in sorted(affected)]
    after = [quality(updates[i]['points'],updates[i]['triangles'],updates[i]['segments']) for i in sorted(affected)]
    cost = sum(len(updates[i]['triangles'])-len(caches[i]['triangles']) for i in affected)
    accepted,reason,a,b = joint_acceptance(before,after,cost)
    record = dict(kind='segment',edge_id=ei,parameter=t,local_parent_segment=parent,
                  local_parameter=float(result.segment_parameter),proposer=proposer.id,
                  affected_faces=sorted(affected),unique_point_cost=1,triangle_cost=cost,
                  before=a,after=b,per_face_before=dict(zip(map(str,sorted(affected)),before)),
                  per_face_after=dict(zip(map(str,sorted(affected)),after)),
                  accepted=accepted,reason=reason,xyz=staged.vertices[-1].tolist())
    return (staged,updates) if accepted else None,record


def improve_shared(registry, caches, faces, options, *, max_building_triangles=None, face_ids=None):
    """Commit conforming quality gains under one fixed set of point/work budgets."""
    started = time.monotonic()
    active = [f for f in faces if f.kind!='terrain' and (face_ids is None or f.id in face_ids)]
    building_triangles = sum(len(caches[f.id]['triangles']) for f in active)
    if max_building_triangles is not None and (isinstance(max_building_triangles, bool)
            or not isinstance(max_building_triangles, int) or max_building_triangles < building_triangles):
        raise ValueError('Building triangle ceiling must be an integer at least the current count')
    used = dict.fromkeys(sorted({f.source_key for f in active}), 0)
    additions, rejections, outcomes = [], [], []
    versions = defaultdict(int)
    memo, native_memo = {}, {}
    rejection_counts = Counter()
    requests = trials = 0
    triangle_blocked = False
    detailed = options.detailed_audit
    enabled = options.budgeted

    def limits():
        return [name for name, hit in (
            ('native_request_limit', requests>=MAX_NATIVE_REQUESTS),
            ('joint_trial_limit', trials>=MAX_JOINT_TRIALS),
            ('runtime_limit', time.monotonic()-started>=MAX_GROUP_SECONDS)) if hit]

    def reject(record):
        rejection_counts[record['reason']] += 1
        if detailed or rejection_counts[record['reason']]==1:
            rejections.append(record if detailed else {k:v for k,v in record.items()
                if k in ('kind','proposer','reason','error','accepted')})

    while enabled and len(additions)<options.max_group_vertices and not limits():
        proposals = []
        triangle_blocked = False
        for f in active:
            if limits(): break
            if used[f.source_key]>=options.max_source_vertices: continue
            c = caches[f.id]
            allowed = permitted_segments(f,c,registry,faces) if options.shared_edge_points else {}
            modes = [(None,True)] if c['optional']<options.max_face_vertices else []
            modes += [(segment,False) for segment in sorted(allowed)]
            for segment, interior in modes:
                if limits(): break
                charged = ({f.source_key} if segment is None else
                    {faces[i].source_key for i in registry.edges[allowed[segment]].incidents})
                if any(used[source]>=options.max_source_vertices for source in charged): continue
                neighbours = (() if segment is None else tuple((i,versions[i]) for i in
                    sorted(registry.edges[allowed[segment]].incidents)))
                native_key = (f.id,versions[f.id],segment,interior)
                key = (*native_key,neighbours)
                if key not in memo:
                    memo[key] = None
                    if native_key not in native_memo:
                        requests += 1
                        native_memo[native_key] = improve(c['points'],c['triangles'],c['segments'],
                            max_added_points=1, splittable_segments=[] if segment is None else [segment],
                            allow_interior=interior)
                        result = native_memo[native_key]
                        if detailed:
                            outcomes.append(dict(face_id=f.id, version=versions[f.id], segment=segment,
                                insertion_kind=result.insertion_kind, before=result.before,
                                after=result.after, termination=result.termination))
                    result = native_memo[native_key]
                    if result.inserted_point is None: continue
                    if result.insertion_kind=='segment':
                        candidates = result.segment_trials or (dict(parent_segment=result.parent_segment,
                                                                    segment_parameter=result.segment_parameter),)
                        for rank, candidate in enumerate(candidates,1):
                            if limits(): break
                            trials += 1
                            staged, record = edge_transaction(registry,caches,faces,f,
                                SimpleNamespace(**candidate),allowed,backend=options.backend)
                            record['native_trial_rank'] = rank
                            if staged is None:
                                reject(record)
                                continue
                            record['charged_sources'] = sorted(charged)
                            memo[key] = (staged[0],staged[1],record)
                            break
                    else:
                        trials += 1
                        updated = dict(c,points=np.asarray(result.mesh.points),
                            triangles=np.asarray(result.mesh.triangles,dtype=int),optional=c['optional']+1,
                            ids=c['ids']+[len(registry.vertices)])
                        xyz = f.origin+updated['points'][-1]@f.basis
                        validate_face_geometry(f,np.vstack([registry.vertices,xyz]),
                                           np.asarray(updated['ids'])[updated['triangles']])
                        before = quality(c['points'],c['triangles'],c['segments'])
                        after = quality(updated['points'],updated['triangles'],updated['segments'])
                        cost = len(updated['triangles'])-len(c['triangles'])
                        accepted,reason,a,b = joint_acceptance([before],[after],cost)
                        record = dict(kind='interior',proposer=f.id,affected_faces=[f.id],
                            charged_sources=[f.source_key],unique_point_cost=1,triangle_cost=cost,
                            before=a,after=b,accepted=accepted,reason=reason,xyz=xyz.tolist())
                        if not accepted:
                            reject(record)
                            continue
                        memo[key] = (None,{f.id:updated},record)
                proposal = memo[key]
                if proposal is None: continue
                record = proposal[2]
                if (max_building_triangles is not None and
                        building_triangles+record['triangle_cost']>max_building_triangles):
                    triangle_blocked = True
                    reject(dict(record,accepted=False,reason='building triangle ceiling',
                                passed_geometry_quality=True))
                    continue
                gain = (record['before']['penalty']-record['after']['penalty'])/record['triangle_cost']
                proposals.append((gain,-f.id,-(segment if segment is not None else -1),proposal))
        if not proposals: break
        gain,_,_,(staged,updates,record) = max(proposals,key=lambda row:row[:3])
        updates, record = deepcopy(updates), deepcopy(record)
        if staged is None:
            fi = record['affected_faces'][0]
            updates[fi]['ids'] = caches[fi]['ids']+[len(registry.vertices)]
            registry.vertices.append(np.asarray(record['xyz']))
        else:
            previous_id = len(staged.vertices)-1
            current = registry.copy_for_split(record['edge_id'])
            _, added = current.reconcile([(record['edge_id'],record['parameter'])],
                                        max_vertices=len(current.vertices)+1)
            if added!=1: raise ValueError('Cached edge transaction lost its unique canonical split')
            new_id = len(current.vertices)-1
            for c in updates.values():
                c['ids'] = [new_id if i==previous_id else i for i in c['ids']]
                c['rings'] = [[new_id if i==previous_id else i for i in ring] for ring in c['rings']]
            registry = current
        for fi,c in updates.items():
            caches[fi] = c
            versions[fi] += 1
            if 'rings' in c:
                faces[fi].rings = c['rings']
                faces[fi].boundary = registry.boundary(fi)
        for source in record['charged_sources']: used[source] += 1
        record.update(vertex=len(registry.vertices)-1,gain_per_triangle=gain,
                      reason='optional quality improvement: '+record['reason'])
        additions.append(record)
        building_triangles += record['triangle_cost']
        memo = {key:value for key,value in memo.items()
                if not ({key[0],*[i for i,_ in key[4]]} & set(updates))}
    bindings = limits() if enabled else []
    if enabled and len(additions)>=options.max_group_vertices: bindings.append('group_point_allowance')
    if enabled and any(used[s]>=options.max_source_vertices for s in used): bindings.append('source_point_allowance')
    if enabled and triangle_blocked: bindings.append('building_triangle_ceiling')
    # Describe the final search state without evaluating excluded insertions.
    # A cap excludes a mode, not a proven quality gain or geometric possibility.
    mode_counts = Counter()
    excluded_modes = []
    unevaluated_modes = []
    if enabled:
        for f in active:
            c = caches[f.id]
            allowed = permitted_segments(f,c,registry,faces) if options.shared_edge_points else {}
            modes = [(None,True), *[(segment,False) for segment in sorted(allowed)]]
            for segment, interior in modes:
                charged = ({f.source_key} if interior else
                    {faces[i].source_key for i in registry.edges[allowed[segment]].incidents})
                restrictions = []
                if len(additions)>=options.max_group_vertices:
                    restrictions.append('group_point_allowance')
                if any(used[s]>=options.max_source_vertices for s in charged):
                    restrictions.append('source_point_allowance')
                if interior and c['optional']>=options.max_face_vertices:
                    restrictions.append('face_point_allowance')
                mode = dict(face_id=f.id,kind='interior' if interior else 'segment',segment=segment)
                if restrictions:
                    excluded_modes.append(dict(mode,binding_limits=restrictions,
                                               additional_gain='unknown'))
                    mode_counts['excluded_by_allowance'] += 1
                elif (f.id,versions[f.id],segment,interior) in native_memo:
                    mode_counts['evaluated'] += 1
                else:
                    unevaluated_modes.append(mode)
                    mode_counts['not_evaluated'] += 1
        if any('face_point_allowance' in m['binding_limits'] for m in excluded_modes):
            bindings.append('face_point_allowance')
    primary = bindings[0] if bindings else ('constraints_prevent_improvement' if enabled and not requests
              else 'marginal_gain_exhausted' if enabled else 'disabled')
    stopping = dict(primary_reason=primary,binding_limits=bindings,
                    category='budget_reached' if bindings else primary,
                    rejection_counts=dict(rejection_counts),
                    evaluated_finite_set=enabled and not bindings and not unevaluated_modes,
                    finite_set_scope='all final-state interior and permitted segment modes, before point allowances',
                    mode_counts=dict(mode_counts),excluded_modes=excluded_modes,
                    unevaluated_modes=unevaluated_modes,
                    mixed_mode_restrictions=bool(excluded_modes and mode_counts['evaluated']),
                    restriction_scope='listed excluded modes only; other modes retain separate evaluation status',
                    gain_exhaustion_scope='evaluated eligible modes only')
    return registry, dict(optional_vertices=additions,source_optional_vertices=used,
        rejected_transactions=rejections,native_outcomes=outcomes,rejection_counts=dict(rejection_counts),
        stopping=stopping,native_requests=requests,joint_trials=trials,
        limits=dict(native_requests=MAX_NATIVE_REQUESTS,joint_trials=MAX_JOINT_TRIALS,
                    building_triangles=max_building_triangles,seconds=MAX_GROUP_SECONDS),
        building_triangles=building_triangles,seconds=time.monotonic()-started,termination=primary)
