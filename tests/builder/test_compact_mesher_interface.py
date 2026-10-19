"""Canonical roof-edge transactions and incident-face quality/ownership guards."""
from copy import deepcopy
import numpy as np
import pytest
from shapely.geometry import Polygon

from dtcc_core.builder.meshing.compact_geometry import planar_face,geometric_graph
from dtcc_core.builder.meshing.compact_triangulation import improve,triangulate
from dtcc_core.builder.meshing.compact_improvement import (
    permitted_segments,edge_transaction,improve_shared,joint_acceptance,quality)
from dtcc_core.builder.meshing.compact_surface import CompactSurfaceOptions


@pytest.mark.parametrize('backend',['geos','triangle'])
def test_shared_edge_backend_contract_is_validated_at_option_construction(backend):
    with pytest.raises(ValueError,match='requires backend="dtcc_mesher"'):
        CompactSurfaceOptions(backend=backend)
    assert CompactSurfaceOptions(backend=backend,budgeted=False).backend==backend
    assert CompactSurfaceOptions(backend=backend,shared_edge_points=False).backend==backend


def fixture(*,ridge=False,neighbour_width=10.,height_offset=0.):
    bottom=np.array([[x,0,height_offset] for x in np.linspace(0,20,21)])
    wall=np.vstack([bottom,[[20,0,8+height_offset],[0,0,8+height_offset]]])
    if ridge:
        wall[:,1]=8-wall[:,2]+height_offset
        wall[:,2]=8+height_offset-.5*wall[:,1]
        normal=[0,.5,1]
    else:normal=[0,-1,0]
    donor=planar_face([wall],kind='roof' if ridge else 'wall',source='a',markers=[0],normal=normal)
    if ridge:
        donor._geometric_rings=[donor._geometric_rings[0][::-1]]
    donor._fixed_boundary=True  # Explicit required opposite-edge junctions.
    y=-neighbour_width if ridge else neighbour_width
    z=8+height_offset-(.5*neighbour_width if ridge else 0.)
    other=planar_face([np.array([[0,0,8+height_offset],[20,0,8+height_offset],[20,y,z],[0,y,z]])],
                     kind='roof',source='b',markers=[1],normal=[0,-.5,1] if ridge else [0,0,1])
    if ridge:
        other._geometric_rings=[other._geometric_rings[0][::-1]]
    faces=[donor,other]
    registry,_=geometric_graph(faces)
    caches={}
    for f in faces:
        ids=sorted(set(x for e in f.boundary for x in e))
        xy=(np.asarray(registry.vertices)[ids]-f.origin)@f.basis.T
        inverse={v:i for i,v in enumerate(ids)}
        segments=np.array([[inverse[a],inverse[b]] for a,b in f.boundary])
        local=[(np.asarray(registry.vertices)[r]-f.origin)@f.basis.T for r in f.rings]
        polygon=Polygon(local[0],local[1:]);t,_=triangulate(xy,segments,polygon)
        caches[f.id]=dict(ids=ids,points=xy,triangles=t,segments=segments,polygon=polygon,optional=0,unchanged=False)
    return registry,caches,faces


def propose(registry,caches,faces):
    f=faces[0];c=caches[0];mapping=permitted_segments(f,c,registry,faces)
    assert len(mapping)==1
    r=improve(c['points'],c['triangles'],c['segments'],max_added_points=1,
              splittable_segments=list(mapping),allow_interior=False)
    assert r.insertion_kind=='segment'
    return r,mapping


@pytest.mark.parametrize('ridge',[False,True])
def test_shared_eave_or_ridge_split_has_exact_point_and_total_triangle_cost(ridge):
    registry,caches,faces=fixture(ridge=ridge)
    result,mapping=propose(registry,caches,faces)
    staged,record=edge_transaction(registry,caches,faces,faces[0],result,mapping,backend='dtcc_mesher')
    assert record['accepted'],record
    new,updates=staged
    assert record['affected_faces']==[0,1] and record['triangle_cost']==2
    assert record['before']['penalty']>record['after']['penalty']
    assert len(new.vertices)==len(registry.vertices)+1
    p=np.array(record['xyz']);ei=record['edge_id'];a,b=np.array(registry.vertices)[list(registry.edges[ei].endpoints)]
    assert np.array_equal(p,a+record['parameter']*(b-a))
    for fi,c in updates.items():
        assert c['ids'][-1]==len(registry.vertices)
        assert np.linalg.norm(c['points'][-1]-(p-faces[fi].origin)@faces[fi].basis.T)<1e-12
    # Opposite local edge orientations reconcile to the same canonical parameter.
    parent=result.parent_segment
    s=caches[0]['segments'][parent]
    a,b=[caches[0]['ids'][int(i)] for i in s]
    assert (a,b) in registry.boundary(0) and (b,a) in registry.boundary(1)
    assert registry.edges[ei].parameters==[0.,1.]


def test_bad_neighbour_rejects_entire_transaction_without_mutation():
    registry,caches,faces=fixture(neighbour_width=.5)
    before=deepcopy((registry,caches,[f.rings for f in faces]))
    r,mapping=propose(registry,caches,faces)
    staged,record=edge_transaction(registry,caches,faces,faces[0],r,mapping,backend='dtcc_mesher')
    assert staged is None and not record['accepted']
    assert 'incident face' in record['reason']
    assert np.array_equal(registry.vertices,before[0].vertices)
    assert [e.parameters for e in registry.edges]==[e.parameters for e in before[0].edges]
    assert [f.rings for f in faces]==before[2]
    assert all(np.array_equal(caches[i]['triangles'],before[1][i]['triangles']) for i in caches)
    assert record['triangle_cost']==2 and len(record['per_face_after'])==2


def test_shared_point_charges_each_source_once_and_respects_interior_cap():
    registry,caches,faces=fixture()
    options=CompactSurfaceOptions(max_face_vertices=0,max_source_vertices=1,max_group_vertices=1)
    r,a=improve_shared(registry,caches,faces,options)
    assert len(a['optional_vertices'])==1
    assert a['optional_vertices'][0]['kind']=='segment'
    assert a['source_optional_vertices']=={'a':1,'b':1}
    assert all(c['optional']==0 for c in caches.values())
    assert a['optional_vertices'][0]['unique_point_cost']==1
    assert a['optional_vertices'][0]['triangle_cost']==2
    r,c,f=fixture();old=np.array(r.vertices)
    r,a=improve_shared(r,c,f,CompactSurfaceOptions(max_source_vertices=0))
    assert not a['optional_vertices'] and np.array_equal(r.vertices,old)


def test_distinct_height_edges_remain_distinct():
    r0,c0,f0=fixture();r1,c1,f1=fixture(height_offset=10.)
    r,_=geometric_graph(f0+f1)
    candidates=[e for e in r.edges if len(e.incidents)==2]
    assert len(candidates)==2
    points=[np.array(r.vertices)[list(e.endpoints)] for e in candidates]
    assert np.allclose(points[0][:,:2],points[1][:,:2])
    assert np.allclose(abs(points[1][:,2]-points[0][:,2]),10.)


def test_remote_global_minimum_is_not_an_acceptance_veto():
    old=dict(penalty=1.,minimum_q=.00001,min_angle=.001,poor_5_count=5,poor_5_area=1.,poor_10_count=6,poor_10_area=2.,poor_20_count=8,poor_20_area=4.)
    new=dict(old,penalty=.5,poor_5_count=3,poor_5_area=.5)
    yes,_,_,_=joint_acceptance([old],[new],2)
    assert yes
    bad=dict(new,poor_5_count=6)
    assert not joint_acceptance([old],[bad],2)[0]
    assert not joint_acceptance([old],[new],100)[0]


def test_thin_neighbour_trial_never_mutates_the_original_registry():
    registry,caches,faces=fixture(neighbour_width=.01)
    before=np.array(registry.vertices)
    result,mapping=propose(registry,caches,faces)
    staged, record = edge_transaction(registry,caches,faces,faces[0],result,mapping,backend='dtcc_mesher')
    assert staged is None and not record['accepted']
    assert 'defects regressed' in record['reason']
    assert np.array_equal(registry.vertices,before)
    assert all(e.parameters==[0.,1.] for e in registry.edges)


def test_cached_disjoint_transactions_rebase_ids_without_losing_prior_splits():
    from dtcc_core.builder.meshing.shared_edges import SharedEdges
    r0,c0,f0=fixture();r1,c1,f1=fixture(height_offset=10.)
    offset=len(r0.vertices)
    for f in f1:
        f.id+=2;f.source_key='c' if f.id==2 else 'd'
        f.rings=[[i+offset for i in ring] for ring in f.rings]
    r=SharedEdges(np.vstack([r0.vertices,r1.vertices]),
        {0:r0.boundary(0),1:r0.boundary(1),2:[(a+offset,b+offset) for a,b in r1.boundary(0)],
         3:[(a+offset,b+offset) for a,b in r1.boundary(1)]})
    caches=dict(c0)
    for i,c in c1.items():caches[i+2]=dict(c,ids=[j+offset for j in c['ids']])
    new,a=improve_shared(r,caches,f0+f1,CompactSurfaceOptions(max_face_vertices=0,max_source_vertices=1,max_group_vertices=2))
    assert len(a['optional_vertices'])==2
    assert a['source_optional_vertices']=={'a':1,'b':1,'c':1,'d':1}
    ids=[p['vertex'] for p in a['optional_vertices']]
    assert len(set(ids))==2 and ids[1]==ids[0]+1
    for p in a['optional_vertices']:
        assert p['vertex'] in new.edges[p['edge_id']].vertices
        assert np.array_equal(new.vertices[p['vertex']],p['xyz'])
        for fi in p['affected_faces']:assert p['vertex'] in caches[fi]['ids']


def test_additional_edge_budget_triangle_ceiling_and_determinism():
    options=CompactSurfaceOptions(max_face_vertices=0,max_source_vertices=4,max_group_vertices=4)
    outcomes=[]
    for _ in range(2):
        r,c,f=fixture()
        start=sum(len(x['triangles']) for x in c.values())
        fixed=np.array(r.vertices)
        r,a=improve_shared(r,c,f,options,max_building_triangles=start+2)
        assert len(a['optional_vertices'])==1 and a['building_triangles']==start+2
        assert all(x['kind']=='segment' for x in a['optional_vertices'])
        assert a['source_optional_vertices']=={'a':1,'b':1}
        assert np.array_equal(r.vertices[:len(fixed)],fixed)
        assert all(x['optional']==0 for x in c.values())
        outcomes.append((np.array(r.vertices),[x['triangles'] for x in c.values()]))
    assert np.array_equal(outcomes[0][0],outcomes[1][0])
    assert all(np.array_equal(a,b) for a,b in zip(outcomes[0][1],outcomes[1][1]))


def test_triangle_ceiling_rejects_before_atomic_commit():
    r,c,f=fixture();before=deepcopy((r,c,[x.rings for x in f]))
    start=sum(len(x['triangles']) for x in c.values())
    result,a=improve_shared(r,c,f,CompactSurfaceOptions(max_face_vertices=0),max_building_triangles=start+1)
    assert not a['optional_vertices']
    assert 'building_triangle_ceiling' in a['stopping']['binding_limits']
    assert np.array_equal(result.vertices,before[0].vertices)
    assert [x.rings for x in f]==before[2]
    assert all(np.array_equal(c[i]['triangles'],before[1][i]['triangles']) for i in c)
    with pytest.raises(ValueError,match='ceiling'):
        improve_shared(r,c,f,CompactSurfaceOptions(),max_building_triangles=start-1)


def test_joint_hip_trial_can_use_native_alternative_without_loosening_guards():
    import json
    from pathlib import Path
    from dtcc_core.builder.meshing.surface_faces import SurfaceFace
    from dtcc_core.builder.meshing.shared_edges import SharedEdges
    data=json.loads((Path(__file__).parent/'data/compact_hip_joint_trial.json').read_text())
    vertices=np.array(data['vertices']);faces=[];caches={}
    for i,row in enumerate(data['faces']):
        origin,basis,normal=(np.array(row[k]) for k in ('origin','basis','normal'))
        local=[(vertices[r]-origin)@basis.T for r in row['rings']]
        polygon=Polygon(local[0],local[1:])
        boundary=[(a,b) for r in row['rings'] for a,b in zip(r,np.roll(r,-1))]
        faces.append(SurfaceFace(i,row['kind'],'a',[i],np.array([],int),origin,basis,normal,
                                 polygon,boundary,row['rings'],0.))
        caches[i]=dict(ids=row['ids'],points=np.array(row['points']),triangles=np.array(row['triangles']),
            segments=np.array(row['segments']),optional=row['optional'],unchanged=False,polygon=polygon)
    registry=SharedEdges(vertices,{f.id:f.boundary for f in faces})
    c=caches[1];allowed=permitted_segments(faces[1],c,registry,faces)
    native=improve(c['points'],c['triangles'],c['segments'],max_added_points=1,
                   splittable_segments=list(allowed),allow_interior=False)
    assert len(native.segment_trials)==4
    staged,record=edge_transaction(registry,caches,faces,faces[1],native,allowed,backend='dtcc_mesher')
    assert staged is None and 'incident face' in record['reason']
    before=np.array(registry.vertices)
    result,a=improve_shared(registry,caches,faces,CompactSurfaceOptions(max_face_vertices=0,
        max_source_vertices=1,max_group_vertices=1),max_building_triangles=sum(len(c['triangles']) for c in caches.values())+2)
    assert len(a['optional_vertices'])==1
    record=a['optional_vertices'][0]
    assert record['native_trial_rank']==4 and record['triangle_cost']==2
    assert record['after']['poor_20_count']<record['before']['poor_20_count']
    assert np.array_equal(result.vertices[:-1],before)
    for i,old in record['per_face_before'].items():
        new=record['per_face_after'][i]
        assert new['minimum_q']>=old['minimum_q']-1e-10
        for angle in (5,10):
            assert new[f'poor_{angle}_count']<=old[f'poor_{angle}_count']
            assert new[f'poor_{angle}_area']<=old[f'poor_{angle}_area']+1e-10*max(1,old[f'poor_{angle}_area'])
