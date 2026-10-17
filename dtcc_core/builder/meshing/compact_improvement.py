"""Bounded source-budgeted transactions on the compact 3D face graph."""
from collections import Counter, defaultdict
from copy import deepcopy
from types import SimpleNamespace
import time

import numpy as np
from shapely.geometry import Polygon

from .compact_triangulation import improve, triangulate
from .surface_faces import boundary_rings
from .surface_refinement import _validate_geometry

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
    points = (np.asarray(registry.vertices)[ids]-face.origin) @ face.basis.T
    segments = np.array([[inverse[a],inverse[b]] for a,b in boundary],dtype=int)
    rings = boundary_rings(boundary)
    rings.sort(key=lambda r:-Polygon((np.asarray(registry.vertices)[r]-face.origin)@face.basis.T).area)
    local = [(np.asarray(registry.vertices)[r]-face.origin)@face.basis.T for r in rings]
    polygon = Polygon(local[0],local[1:])
    triangles,flips = triangulate(points,segments,polygon,backend=backend)
    _validate_geometry(face,np.asarray(registry.vertices),np.array(ids)[triangles])
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
    staged = deepcopy(registry)
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


def improve_shared(registry, caches, faces, options, *, max_building_triangles=None, face_ids=None,
                   adaptive=False, _standard_checkpoint=None):
    """Select guarded transactions; optionally continue the fixed allowance in place.

    Adaptive continuation is edge-only, with total hard caps of 8/source and
    16/group. Search counters, versions, caches and the triangle ceiling span
    both phases. Effective allowances increase only with an accepted point.
    """
    if adaptive and (options.max_face_vertices or not options.shared_edge_points):
        raise ValueError("Adaptive continuation requires shared-edge-only improvement")
    started=time.monotonic()
    active=[f for f in faces if f.kind!='terrain' and (face_ids is None or f.id in face_ids)]
    building_triangles = sum(len(caches[f.id]['triangles']) for f in active)
    if max_building_triangles is not None:
        if (isinstance(max_building_triangles,bool) or not isinstance(max_building_triangles,int)
                or max_building_triangles < building_triangles):
            raise ValueError('Building triangle ceiling must be an integer at least the current count')
    used = dict.fromkeys(sorted({f.source_key for f in active}),0)
    additions,rejections = [],[]
    rejection_counts=Counter()
    detailed=getattr(options,'detailed_audit',False)
    def reject(record):
        record.setdefault('passed_geometry_quality',False)
        rejection_counts[record['reason']]+=1
        if detailed or rejection_counts[record['reason']]==1:
            rejections.append(record if detailed else {k:record[k] for k in
                ('kind','proposer','reason','error','accepted','passed_geometry_quality') if k in record})
    versions = defaultdict(int)
    memo = {}
    native_memo = {}
    outcomes = []
    reproducers = []
    requests = 0
    joint_trials = 0
    phase = 'standard'
    effective_source = dict.fromkeys(used, options.max_source_vertices)
    effective_group = options.max_group_vertices
    hard_source = 8 if adaptive else options.max_source_vertices
    hard_group = 16 if adaptive else options.max_group_vertices
    escalations = []
    standard_result = None
    assessments = Counter()
    # Records refer to current incident versions; stale constraints never decide
    # final stopping. Historical rejections remain separately counted.
    triangle_blocked = {}
    allowance_unassessed = set()
    search_incomplete = False
    enabled = options.budgeted and options.max_group_vertices > 0 and options.max_source_vertices > 0
    eligible = any(permitted_segments(f,caches[f.id],registry,faces) for f in active) if enabled and options.shared_edge_points else False

    def search_limits():
        return ([name for name, hit in (
            ('native_request_limit', requests >= MAX_NATIVE_REQUESTS),
            ('joint_trial_limit', joint_trials >= MAX_JOINT_TRIALS),
            ('runtime_limit', time.monotonic()-started >= MAX_GROUP_SECONDS)) if hit])

    def snapshot():
        return dict(optional_points=len(additions), source_charges=dict(used),
                    building_triangles=building_triangles, native_requests=requests,
                    joint_trials=joint_trials)

    def begin_continuation():
        # A private recovery copy lets the caller retain the complete standard
        # prefix if final whole-shell validation rejects the optional extension.
        if _standard_checkpoint is not None:
            _standard_checkpoint.update(deepcopy(dict(registry=registry,caches=caches,
                faces=faces,optional_vertices=additions,source_optional_vertices=used,
                standard_result=standard_result)))

    while enabled:
        cap_source = hard_source if phase == 'adaptive' else options.max_source_vertices
        cap_group = hard_group if phase == 'adaptive' else options.max_group_vertices
        if len(additions)>=cap_group or search_limits():
            if search_limits(): search_incomplete = True
            if phase=='standard':
                standard_result=snapshot()
                if adaptive and len(additions)<hard_group and not search_limits():
                    begin_continuation()
                    phase='adaptive'
                    continue
            break
        proposals = []
        for f in active:
            if used[f.source_key]>=cap_source:
                allowance_unassessed.add(f.id)
                continue
            if phase=='adaptive':allowance_unassessed.discard(f.id)
            c = caches[f.id]
            allowed = permitted_segments(f,c,registry,faces) if options.shared_edge_points else {}
            # Separate bounded native calls make edge proposals visible even
            # when a cheaper local interior winner has a different joint cost.
            modes = [(None,True)] if c['optional']<options.max_face_vertices else []
            modes += [(si,False) for si in sorted(allowed)]
            for segment,interior in modes:
                charged = ({f.source_key} if segment is None else
                    {faces[i].source_key for i in registry.edges[allowed[segment]].incidents})
                if phase=='adaptive' and not (len(additions)>=effective_group or
                        any(used[s]>=effective_source[s] for s in charged)):
                    continue
                neighbours = (() if segment is None else tuple((i,versions[i]) for i in
                    sorted(registry.edges[allowed[segment]].incidents)))
                native_key = (f.id,versions[f.id],segment,interior)
                key = (*native_key,neighbours)
                if key not in memo:
                    if (requests>=MAX_NATIVE_REQUESTS or joint_trials>=MAX_JOINT_TRIALS
                            or time.monotonic()-started>=MAX_GROUP_SECONDS):
                        search_incomplete=True
                        break
                    if native_key not in native_memo:
                        requests += 1
                        try:
                            native_memo[native_key] = improve(c['points'],c['triangles'],c['segments'],max_added_points=1,
                                splittable_segments=[] if segment is None else [segment],allow_interior=interior)
                        except (ValueError,RuntimeError) as exc:
                            native_memo[native_key]=None
                            reject(dict(kind='native',proposer=f.id,reason='unexpected native numerical failure',
                                        error=str(exc),accepted=False))
                        native = native_memo[native_key]
                        if native is None:
                            memo[key]=None
                            continue
                        if detailed:outcomes.append(dict(face_id=f.id,version=versions[f.id],segment=segment,
                            edge_id=allowed.get(segment),allow_interior=interior,kind=native.insertion_kind,
                            before=native.before,after=native.after,termination=native.termination,
                            candidate_evaluations=native.candidate_evaluations,
                            segment_trials=list(getattr(native,'segment_trials',()))))
                        if detailed and segment is not None and versions[f.id]==0:
                            reproducers.append(dict(face_id=f.id,edge_id=allowed[segment],segment=int(segment),
                                points=c['points'].tolist(),triangles=c['triangles'].tolist(),
                                segments=c['segments'].tolist(),incident_faces=sorted(registry.edges[allowed[segment]].incidents)))
                    result = native_memo[native_key]
                    if result is None or result.inserted_point is None:
                        memo[key]=None
                        continue
                    if result.insertion_kind=='segment':
                        affected_sources = {faces[i].source_key for i in registry.edges[allowed[segment]].incidents
                                            if faces[i].source_key is not None}
                        if any(used[s]>=cap_source for s in affected_sources):
                            allowance_unassessed.add(f.id)
                            # Keep the native result for a later allowance assessment.
                            reject(dict(kind='segment',proposer=f.id,edge_id=allowed[segment],reason='affected source budget',accepted=False))
                            continue
                        trials=getattr(result,'segment_trials',()) or (
                            dict(parent_segment=result.parent_segment,segment_parameter=result.segment_parameter),)
                        memo[key]=None
                        for rank,trial in enumerate(trials,1):
                            if joint_trials>=MAX_JOINT_TRIALS or time.monotonic()-started>=MAX_GROUP_SECONDS:
                                search_incomplete=True
                                break
                            joint_trials+=1
                            try:
                                staged,record=edge_transaction(registry,caches,faces,f,
                                    SimpleNamespace(**trial),allowed,backend=options.backend)
                            except (ValueError,RuntimeError) as exc:
                                staged,record=None,dict(kind='segment',proposer=f.id,
                                    edge_id=allowed[segment],local_parameter=trial['segment_parameter'],
                                    reason='invalid incident retessellation: '+str(exc),accepted=False)
                            record['native_trial_rank']=rank
                            if staged is None:
                                reject(record)
                                continue
                            record['charged_sources']=sorted(affected_sources)
                            assessments['passed_geometry_quality']+=1
                            memo[key]=(staged[0],staged[1],record)
                            # Alternatives rescue an incompatible local winner.
                            # Stop at the first guarded trial to preserve the
                            # joint-search budget for other permitted intervals.
                            break
                    else:
                        joint_trials += 1
                        updated = dict(c,points=np.asarray(result.mesh.points),triangles=np.asarray(result.mesh.triangles,dtype=int),optional=c['optional']+1)
                        xyz=f.origin+updated['points'][-1]@f.basis
                        updated['ids']=c['ids']+[len(registry.vertices)]
                        try:
                            _validate_geometry(f,np.vstack([registry.vertices,xyz]),np.array(updated['ids'])[updated['triangles']])
                            before = quality(c['points'],c['triangles'],c['segments'])
                            after = quality(updated['points'],updated['triangles'],updated['segments'])
                        except (ValueError,RuntimeError) as exc:
                            reject(dict(kind='interior',proposer=f.id,reason='invalid interior transaction',error=str(exc),accepted=False))
                            memo[key]=None
                            continue
                        cost = len(updated['triangles'])-len(c['triangles'])
                        accepted,reason,a,b = joint_acceptance([before],[after],cost)
                        record = dict(kind='interior',proposer=f.id,affected_faces=[f.id],charged_sources=[f.source_key],
                                      unique_point_cost=1,triangle_cost=cost,before=a,after=b,accepted=accepted,reason=reason,
                                      xyz=xyz.tolist())
                        if not accepted:
                            reject(record);memo[key]=None;continue
                        assessments['passed_geometry_quality']+=1
                        updated['ids']=c['ids']+[len(registry.vertices)]
                        memo[key]=(None,{f.id:updated},record)
                proposal = memo.get(key)
                if proposal is not None:
                    record=proposal[2]
                    if (max_building_triangles is not None and
                            building_triangles+record['triangle_cost']>max_building_triangles):
                        reject(dict(record,accepted=False,reason='building triangle ceiling',passed_geometry_quality=True))
                        triangle_blocked[key]=record
                        continue
                    if any(used[s]>=cap_source for s in record['charged_sources']):
                        allowance_unassessed.add(f.id)
                        continue
                    if phase=='adaptive':
                        assessments['passed_blocked_only_by_allowance']+=1
                    proposals.append(((record['before']['penalty']-record['after']['penalty'])/record['triangle_cost'],
                                      -f.id,-(segment if segment is not None else -1),key,proposal))
        if not proposals:
            if phase=='standard':
                standard_result=snapshot()
                capped=(len(additions)>=effective_group or any(used[s]>=effective_source[s] for s in used))
                if adaptive and capped and len(additions)<hard_group and not search_limits():
                    begin_continuation()
                    phase='adaptive'
                    continue
            break
        _,_,_,key,(staged,updates,record) = max(proposals,key=lambda x:x[:3])
        updates = deepcopy(updates)
        record = deepcopy(record)
        record['phase']=phase
        record['gain_per_triangle']=(record['before']['penalty']-record['after']['penalty'])/record['triangle_cost']
        if phase=='adaptive':
            source_changes={s:dict(before=effective_source[s],after=used[s]+1)
                for s in record['charged_sources'] if used[s]>=effective_source[s]}
            group_after=max(effective_group,len(additions)+1)
            escalations.append(dict(point_index=len(additions),sources=source_changes,
                group_before=effective_group,group_after=group_after,
                affected_faces=record['affected_faces'],edge_id=record['edge_id'],
                gain_per_triangle=record['gain_per_triangle'],triangle_cost=record['triangle_cost'],
                reason='guarded complete transaction fits unchanged triangle ceiling'))
            for source,change in source_changes.items():effective_source[source]=change['after']
            effective_group=group_after
        if staged is None:
            fi=record['affected_faces'][0]
            updates[fi]['ids']=caches[fi]['ids']+[len(registry.vertices)]
            registry.vertices.append(np.array(record['xyz']))
        else:
            # Cached geometry/quality is reusable while incident versions agree.
            # Rebase only its appended identity onto the current complete graph.
            previous_id=len(staged.vertices)-1
            current=deepcopy(registry)
            _,added=current.reconcile([(record['edge_id'],record['parameter'])],
                                      max_vertices=len(current.vertices)+1)
            if added!=1:
                raise ValueError('Cached edge transaction no longer has a unique split')
            new_id=len(current.vertices)-1
            for c in updates.values():
                c['ids']=[new_id if i==previous_id else i for i in c['ids']]
                c['rings']=[[new_id if i==previous_id else i for i in r] for r in c['rings']]
            registry=current
        for fi,c in updates.items():
            caches[fi]=c
            versions[fi]+=1
            if 'rings' in c:
                faces[fi].rings=c['rings']
                faces[fi].boundary=registry.boundary(fi)
        for source in record['charged_sources']:
            used[source]+=1
        record['vertex']=len(registry.vertices)-1
        record['reason']='optional quality improvement: '+record['reason']
        additions.append(record)
        building_triangles += record['triangle_cost']
        # Incident versions define geometric validity; unrelated appends only
        # require identity rebasing, not another native search/quality trial.
        triangle_blocked={k:v for k,v in triangle_blocked.items()
              if not ({k[0],*[i for i,_ in k[4]]} & set(updates))}
        allowance_unassessed.difference_update(updates)
        memo={k:v for k,v in memo.items()
              if not ({k[0],*[i for i,_ in k[4]]} & set(updates))}
    if standard_result is None:standard_result=snapshot()
    # Only current, uncommitted proposals explain the final point limit.
    # Accepted adaptive gains remain in escalations and optional_vertices.
    best_allowance_blocked = None
    for proposal in memo.values():
        if proposal is None:continue
        record=proposal[2]
        blocked=(len(additions)>=cap_group or
                 any(used[s]>=cap_source for s in record['charged_sources']))
        fits=max_building_triangles is None or building_triangles+record['triangle_cost']<=max_building_triangles
        if blocked and fits:
            gain=(record['before']['penalty']-record['after']['penalty'])/record['triangle_cost']
            if best_allowance_blocked is None or gain>best_allowance_blocked['gain_per_triangle']:
                best_allowance_blocked=dict({k:v for k,v in record.items() if k not in ('per_face_before','per_face_after')},
                    gain_per_triangle=gain,passed_geometry_quality=True,
                    accepted=False,reason='point allowance at termination',
                    blocked_only_by_allowance=True,at_point_count=len(additions))
    bindings=search_limits() if enabled else []
    source_reached=[s for s in used if used[s]>=cap_source] if enabled else []
    if not enabled:
        primary='disabled'
    else:
        if len(additions)>=cap_group:
            bindings.append('hard_group_point_limit' if phase=='adaptive' or cap_group>=16 else 'group_point_allowance')
        if source_reached:
            bindings.append('hard_source_point_limit' if phase=='adaptive' or cap_source>=8 else 'source_point_allowance')
        if triangle_blocked:bindings.append('building_triangle_ceiling')
        if rejection_counts.get('unexpected native numerical failure'):
            bindings.append('unexpected_numerical_failure')
        primary=next(iter(bindings), 'evaluated_candidate_exhaustion' if eligible or options.max_face_vertices else 'no_eligible_shared_edges')
    stopping=dict(primary_reason=primary,binding_limits=bindings,
        rejection_counts=dict(rejection_counts),
        evaluated_finite_set=bool(enabled and not (bindings or search_incomplete or allowance_unassessed)),
        assessment=dict(assessments,triangle_blocked=len(triangle_blocked),
            allowance_unassessed_faces=sorted(allowance_unassessed),
            search_incomplete=bool(enabled and (search_incomplete or search_limits()))),
        best_fully_evaluated_allowance_blocked=best_allowance_blocked)
    compatibility={'native_request_limit':'request_budget','joint_trial_limit':'joint_trial_budget',
        'runtime_limit':'runtime_budget','group_point_allowance':'point_budget',
        'hard_group_point_limit':'point_budget','source_point_allowance':'source_budget',
        'hard_source_point_limit':'source_budget','building_triangle_ceiling':'triangle_budget',
        'evaluated_candidate_exhaustion':'no_beneficial_transaction'}
    return registry,dict(optional_vertices=additions,source_optional_vertices=used,
        rejected_transactions=rejections,native_outcomes=outcomes,initial_edge_reproducers=reproducers,
        rejection_counts=dict(rejection_counts),stopping=stopping,
        standard_result=standard_result,escalations=escalations,
        effective_allowances=dict(source=effective_source,group=effective_group),
        hard_limits=dict(source=hard_source,group=hard_group),
        native_requests=requests,joint_trials=joint_trials,
        limits=dict(native_requests=MAX_NATIVE_REQUESTS,joint_trials=MAX_JOINT_TRIALS,
                    building_triangles=max_building_triangles,seconds=MAX_GROUP_SECONDS),building_triangles=building_triangles,
        seconds=time.monotonic()-started,termination=compatibility.get(primary,primary))
