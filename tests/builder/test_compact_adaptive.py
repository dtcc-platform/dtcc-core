"""Cumulative adaptive allowances and evidence-based stopping on native fixtures."""
from collections import Counter
from copy import deepcopy
import importlib

import numpy as np
import pytest

from .test_compact_mesher_interface import fixture
from .test_compact_integration import two_groups
from dtcc_core.builder import CompactSurfaceOptions

ci=importlib.import_module('dtcc_core.builder.meshing.compact_improvement')


def policy(source=1,group=1,**kwargs):
    return CompactSurfaceOptions(max_face_vertices=0,max_source_vertices=source,
                                 max_group_vertices=group,**kwargs)


def paired():
    r,c,f=two_groups()
    # Separate geometric intervals with shared source identities expose reusable
    # disjoint cached work after both sources reach their first allowance.
    f[2].source_key='a';f[3].source_key='b'
    return r,c,f


def run(options=None,*,adaptive=False,ceiling=None,paired_input=False):
    r,c,f=paired() if paired_input else fixture()
    r,a=ci.improve_shared(r,c,f,options or policy(),adaptive=adaptive,
                          max_building_triangles=ceiling)
    return r,c,f,a


def test_triangle_budget_has_positive_quality_evidence_and_atomic_rejection():
    r,c,f=fixture();initial=deepcopy((r,c));t=sum(len(v['triangles']) for v in c.values())
    r,a=ci.improve_shared(r,c,f,policy(4,8),max_building_triangles=t+1)
    stop=a['stopping']
    assert stop['primary_reason']=='building_triangle_ceiling'
    assert stop['assessment']['triangle_blocked']>0
    assert a['termination']=='triangle_budget'
    assert all(x['passed_geometry_quality'] for x in a['rejected_transactions'])
    assert not stop['evaluated_finite_set']
    assert np.array_equal(r.vertices,initial[0].vertices)
    assert all(np.array_equal(v['triangles'],initial[1][i]['triangles']) for i,v in c.items())


@pytest.mark.parametrize('source,group,binding',[(1,8,'source_point_allowance'),(4,1,'group_point_allowance')])
def test_separate_soft_allowances(source,group,binding):
    *_,a=run(policy(source,group))
    assert binding in a['stopping']['binding_limits']
    assert not a['stopping']['evaluated_finite_set']


def test_adaptive_preserves_standard_prefix_and_reuses_disjoint_cached_work():
    *_,fixed=run(paired_input=True)
    r,c,f,a=run(adaptive=True,paired_input=True)
    n=a['standard_result']['optional_points']
    assert n==1 and len(a['optional_vertices'])==2
    assert a['optional_vertices'][:n]==fixed['optional_vertices']
    assert a['standard_result']['native_requests']==fixed['native_requests']
    assert a['standard_result']['joint_trials']==fixed['joint_trials']
    assert a['native_requests']>=fixed['native_requests']
    assert a['source_optional_vertices']=={'a':2,'b':2}
    assert a['effective_allowances']=={'source':{'a':2,'b':2},'group':2}
    assert len(a['escalations'])==1
    event=a['escalations'][0]
    assert event['sources']=={'a':{'before':1,'after':2},'b':{'before':1,'after':2}}
    assert event['gain_per_triangle']>.015 and event['triangle_cost']==2
    assert a['stopping']['best_fully_evaluated_allowance_blocked'] is None
    assert a['building_triangles']==a['standard_result']['building_triangles']+2
    for point in a['optional_vertices']:
        assert point['affected_faces']==sorted(r.edges[point['edge_id']].incidents)
        assert all(point['vertex'] in c[i]['ids'] for i in point['affected_faces'])


def test_stopping_candidate_is_current_and_uncommitted():
    r,_,_,fixed=run(paired_input=True)
    best=fixed['stopping']['best_fully_evaluated_allowance_blocked']
    assert best is not None and best['passed_geometry_quality']
    assert best['blocked_only_by_allowance'] and not best['accepted']
    assert best['at_point_count']==len(fixed['optional_vertices'])
    assert best['parameter'] not in r.edges[best['edge_id']].parameters
    assert best['gain_per_triangle']>.015

    # Continuation consumes the cached opportunity; its gain is historical now.
    *_,adaptive=run(adaptive=True,paired_input=True)
    assert adaptive['stopping']['primary_reason']=='evaluated_candidate_exhaustion'
    assert adaptive['stopping']['best_fully_evaluated_allowance_blocked'] is None
    point=adaptive['optional_vertices'][fixed['standard_result']['optional_points']]
    assert (point['edge_id'],point['parameter'])==(best['edge_id'],best['parameter'])
    assert adaptive['escalations'][0]['gain_per_triangle']==best['gain_per_triangle']


def test_no_escalation_when_unavoidable_poor_triangles_remain():
    *_,a=run(adaptive=True)
    assert len(a['optional_vertices'])==1
    assert not a['escalations']
    assert a['stopping']['primary_reason']=='evaluated_candidate_exhaustion'
    assert a['stopping']['evaluated_finite_set']
    assert a['rejection_counts']  # The same incident guards reject further trials.


def test_adaptive_keeps_original_triangle_ceiling():
    r,c,f=paired();t=sum(len(c[i]['triangles']) for i in range(4))
    r,a=ci.improve_shared(r,c,f,policy(),adaptive=True,max_building_triangles=t+3)
    assert len(a['optional_vertices'])==1 and not a['escalations']
    assert a['limits']['building_triangles']==t+3
    assert a['building_triangles']==t+2
    assert a['stopping']['primary_reason']=='building_triangle_ceiling'


@pytest.mark.parametrize('budget,reason',[('MAX_NATIVE_REQUESTS','native_request_limit'),('MAX_JOINT_TRIALS','joint_trial_limit')])
def test_continuation_does_not_reset_search_budget(monkeypatch,budget,reason):
    *_,fixed=run(paired_input=True)
    limit=fixed['native_requests' if budget=='MAX_NATIVE_REQUESTS' else 'joint_trials']
    monkeypatch.setattr(ci,budget,limit)
    *_,a=run(adaptive=True,paired_input=True)
    assert len(a['optional_vertices'])==1
    assert reason in a['stopping']['binding_limits']
    assert a['stopping']['assessment']['search_incomplete']
    assert not a['escalations']
    assert a['native_requests']<=ci.MAX_NATIVE_REQUESTS and a['joint_trials']<=ci.MAX_JOINT_TRIALS
    assert 'group_point_allowance' in a['stopping']['binding_limits']


def test_simultaneous_triangle_point_search_limits(monkeypatch):
    r,c,f=paired();t=sum(len(c[i]['triangles']) for i in range(4))
    # Fully assess a cheap and an expensive proposal before the joint limit.
    original=ci.edge_transaction
    def costed(*args,**kwargs):
        staged,record=original(*args,**kwargs)
        if staged is not None and record['proposer']>=2:record['triangle_cost']=4
        return staged,record
    monkeypatch.setattr(ci,'edge_transaction',costed)
    monkeypatch.setattr(ci,'MAX_JOINT_TRIALS',4)
    r,a=ci.improve_shared(r,c,f,policy(),max_building_triangles=t+2)
    assert {'joint_trial_limit','group_point_allowance','building_triangle_ceiling'}<=set(a['stopping']['binding_limits'])


def test_controlled_clock_exhaustion(monkeypatch):
    now=[0.];original=ci.improve
    def elapsed(*args,**kwargs):
        result=original(*args,**kwargs)
        if kwargs.get('max_added_points')==1:now[0]=ci.MAX_GROUP_SECONDS+1
        return result
    monkeypatch.setattr(ci.time,'monotonic',lambda:now[0])
    monkeypatch.setattr(ci,'improve',elapsed)
    *_,a=run(policy(4,8),adaptive=True)
    assert a['stopping']['primary_reason']=='runtime_limit'
    assert not a['stopping']['evaluated_finite_set'] and not a['escalations']


def test_disabled_and_no_eligible_edges():
    *_,a=run(policy(budgeted=False))
    assert a['stopping']['primary_reason']=='disabled' and a['native_requests']==0
    assert not a['stopping']['evaluated_finite_set']
    *_,a=run(policy(shared_edge_points=False))
    assert a['stopping']['primary_reason']=='no_eligible_shared_edges'
    assert a['native_requests']==0


@pytest.mark.parametrize('group_count',[1,2])
@pytest.mark.parametrize('options,enabled_stages',[
    (CompactSurfaceOptions(budgeted=False),0),
    (CompactSurfaceOptions(max_source_vertices=0,max_group_vertices=0,
        extension_source_vertices=0,extension_group_vertices=0),0),
    (CompactSurfaceOptions(extension_source_vertices=0,extension_group_vertices=0),1),
    (CompactSurfaceOptions(shared_edge_points=False),1),
])
def test_disabled_stages_do_not_revalidate_whole_shell(monkeypatch,group_count,options,enabled_stages):
    from shapely.geometry import box
    from .test_roof_catalogue import assemble
    from dtcc_core.builder.roofs import RoofPatch
    from dtcc_core.builder.meshing.compact_surface import compact_surface
    validation=importlib.import_module('dtcc_core.builder.meshing.roof_validation')
    patches=[RoofPatch(box(30*i-18,-5,30*i-3,5),np.zeros(2),np.array([0.,0.,8.]),
                      str(i),block_id=i) for i in range(group_count)]
    reference,patches=assemble(patches,slope=.4)
    calls=[];original=validation.validate_shell
    def counted(*args,**kwargs):
        calls.append(True)
        return original(*args,**kwargs)
    monkeypatch.setattr(validation,'validate_shell',counted)
    if not enabled_stages:
        def unexpected_assessment(*args,**kwargs):
            pytest.fail('Disabled refinement assessed shared-edge candidates')
        monkeypatch.setattr(ci,'permitted_segments',unexpected_assessment)
    _,audit=compact_surface(reference,patches,options)
    stages=audit['optimization']['stages']
    assert len(audit['optimization']['groups'])==group_count
    assert len(stages)==2*group_count
    assert sum(s['stopping']['primary_reason']=='disabled' for s in stages)==(2-enabled_stages)*group_count
    # Initial/final validation is constant; only enabled stages add shell work.
    assert len(calls)==2+enabled_stages*group_count
    assert audit['topology']['valid']
    for stage in stages:
        if stage['stopping']['primary_reason']=='disabled':
            assert stage['native_requests']==stage['joint_trials']==0
            assert not stage['optional_vertices']


def test_unexpected_native_failure_is_not_search_exhaustion(monkeypatch):
    original=ci.improve
    def failed(*args,**kwargs):
        if kwargs.get('max_added_points')==1:raise RuntimeError('controlled native failure')
        return original(*args,**kwargs)
    monkeypatch.setattr(ci,'improve',failed)
    *_,a=run(policy(4,8),adaptive=True)
    assert a['stopping']['primary_reason']=='unexpected_numerical_failure'
    assert not a['escalations'] and not a['optional_vertices']


def test_adaptive_atomic_failure_keeps_standard_prefix(monkeypatch):
    *_,fixed=run()
    original=ci.edge_transaction
    def fail_second(registry,*args,**kwargs):
        if any(len(e.vertices)>2 for e in registry.edges):raise ValueError('controlled geometry rejection')
        return original(registry,*args,**kwargs)
    monkeypatch.setattr(ci,'edge_transaction',fail_second)
    r,c,f,a=run(adaptive=True)
    assert a['optional_vertices']==fixed['optional_vertices']
    assert a['source_optional_vertices']==fixed['source_optional_vertices']
    assert not a['escalations']


def test_detailed_audit_has_identical_adaptive_geometry():
    states=[]
    for detailed in (False,True):
        r,c,f,a=run(policy(detailed_audit=detailed),adaptive=True,paired_input=True)
        states.append((np.array(r.vertices),[v['triangles'] for v in c.values()],a))
    assert np.array_equal(states[0][0],states[1][0])
    assert all(np.array_equal(a,b) for a,b in zip(states[0][1],states[1][1]))
    assert states[0][2]['optional_vertices']==states[1][2]['optional_vertices']
    assert not states[0][2]['native_outcomes'] and states[1][2]['native_outcomes']


def test_private_standard_checkpoint_retains_exact_prefix():
    r,c,f=paired();checkpoint={}
    _,a=ci.improve_shared(r,c,f,policy(),adaptive=True,_standard_checkpoint=checkpoint)
    r0,c0,f0=paired();r0,fixed=ci.improve_shared(r0,c0,f0,policy())
    assert checkpoint['optional_vertices']==fixed['optional_vertices']
    assert checkpoint['source_optional_vertices']==fixed['source_optional_vertices']
    assert np.array_equal(checkpoint['registry'].vertices,r0.vertices)
    assert all(np.array_equal(c0[i]['triangles'],checkpoint['caches'][i]['triangles']) for i in c0)
    assert len(a['optional_vertices'])>len(checkpoint['optional_vertices'])


@pytest.mark.parametrize('shared_sources,count,reason',[(True,9,'hard_source_point_limit'),(False,18,'hard_group_point_limit')])
def test_hard_total_point_limits(shared_sources,count,reason):
    from dtcc_core.builder.meshing.shared_edges import SharedEdges
    vertices=[];faces=[];caches={};boundaries={}
    for j in range(count):
        r,c,f=fixture();offset=len(vertices);face_offset=len(faces)
        vertices.extend(np.array(r.vertices)+[j*100,0,0])
        for face in f:
            old=face.id;face.id+=face_offset;face.origin=face.origin+[j*100,0,0]
            if not shared_sources:face.source_key+=str(j)
            face.rings=[[v+offset for v in ring] for ring in face.rings]
            face.boundary=[(a+offset,b+offset) for a,b in face.boundary]
            boundaries[face.id]=face.boundary
            caches[face.id]=dict(c[old],ids=[v+offset for v in c[old]['ids']])
        faces.extend(f)
    registry=SharedEdges(vertices,boundaries)
    _,a=ci.improve_shared(registry,caches,faces,policy(),adaptive=True)
    assert reason in a['stopping']['binding_limits']
    assert len(a['optional_vertices'])==(8 if shared_sources else 16)
    assert all(n<=8 for n in a['source_optional_vertices'].values())
    assert len(a['optional_vertices'])<=16
    assert dict(Counter(s for p in a['optional_vertices'] for s in p['charged_sources']))=={
        s:n for s,n in a['source_optional_vertices'].items() if n}


def test_whole_shell_rejection_recovers_standard_prefix(monkeypatch):
    from dataclasses import replace
    from dtcc_core import builder
    from .test_compact_integration import sample,signature
    validation=importlib.import_module('dtcc_core.builder.meshing.roof_validation')
    original_search=ci.improve_shared;original_validation=validation.validate_shell
    def small_allowance(r,c,f,options,**kwargs):
        if options.max_face_vertices==0:
            options=replace(options,max_source_vertices=1,max_group_vertices=1)
        return original_search(r,c,f,options,**kwargs)
    monkeypatch.setattr(ci,'improve_shared',small_allowance)
    def build(adaptive):
        audit={}
        mesh=builder.build_city_surface_mesh(sample(),lod='lod2',mesher='dtcc_mesher',
            compact_surface=CompactSurfaceOptions(max_source_vertices=0,max_group_vertices=0,
                adaptive_extension=adaptive),report_mesh_quality=False,stage_audit=audit)
        return mesh,audit['roof_catalogue']
    fixed,reference=build(False)
    t0=reference['compact_surface']['T0'];failed=[]
    def reject_extension_once(mesh,**kwargs):
        result=original_validation(mesh,**kwargs)
        if sum(mesh.markers>=0)>t0 and not failed:
            failed.append(True)
            return dict(result,valid=False,controlled_failure=True)
        return result
    monkeypatch.setattr(validation,'validate_shell',reject_extension_once)
    mesh,cat=build(True);audit=cat['compact_surface'];stage=audit['optimization']['stages'][-1]
    assert failed and stage['adaptive_rollback']
    assert stage['termination']=='unexpected_numerical_failure'
    assert len(stage['optional_vertices'])==stage['standard_result']['optional_points']==1
    assert signature(mesh)==signature(fixed)
    assert audit['topology']['valid'] and not cat['block_fallbacks']


@pytest.mark.parametrize('budgeted,shared',[(False,True),(True,False)])
def test_public_disabled_stage_b_with_adaptive_option(budgeted,shared):
    from dtcc_core import builder
    from .test_compact_integration import sample
    audit={}
    builder.build_city_surface_mesh(sample(),lod='lod2',mesher='dtcc_mesher',
        compact_surface=CompactSurfaceOptions(budgeted=budgeted,shared_edge_points=shared,
            adaptive_extension=True),report_mesh_quality=False,stage_audit=audit)
    stage=audit['roof_catalogue']['compact_surface']['optimization']['stages'][-1]
    assert stage['stage']=='B' and stage['stopping']['primary_reason']=='disabled'
    assert not stage['optional_vertices'] and not stage['stopping']['evaluated_finite_set']


def test_adaptive_boolean_is_validated():
    with pytest.raises(TypeError,match='adaptive_extension'):
        CompactSurfaceOptions(adaptive_extension='adaptive')


def test_adaptive_groups_keep_independent_budgets_and_order():
    results=[]
    for reverse in (False,True):
        r,c,f=two_groups();groups=ci.building_groups(r,f)
        if reverse:groups.reverse()
        charges={}
        for group in groups:
            start=sum(len(c[i]['triangles']) for i in group['face_ids'])
            r,a=ci.improve_shared(r,c,f,policy(),face_ids=set(group['face_ids']),
                max_building_triangles=start+2,adaptive=True)
            assert a['standard_result']['optional_points']==1
            assert len(a['optional_vertices'])==1
            charges.update(a['source_optional_vertices'])
        assert charges=={'a':1,'b':1,'c':1,'d':1}
        xyz=np.asarray(r.vertices)
        results.append(sorted(tuple(sorted(map(tuple,xyz[np.array(cache['ids'])[tri]])))
            for cache in c.values() for tri in cache['triangles']))
    assert results[0]==results[1]
