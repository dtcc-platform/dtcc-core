"""Public compact catalogue routing, group isolation and two-stage accounting."""
from copy import deepcopy
from types import SimpleNamespace
import importlib
import numpy as np
import pytest
from shapely.geometry import box

from dtcc_core import builder
from dtcc_core.builder import CompactSurfaceOptions
from dtcc_core.builder.geometry_builders.meshes import _welded_mesh
from dtcc_core.builder.meshing.compact_improvement import building_groups,improve_shared
from dtcc_core.builder.meshing.shared_edges import SharedEdges
from .test_compact_mesher_interface import fixture
from .test_roof_catalogue import building,city_of,model


def sample():
    polygon=box(-10,-5,10,5)
    return city_of([building(polygon,model('gable',[0.,8.,.6],polygon))])


def signature(mesh):
    rows=[]
    for triangle,marker in zip(mesh.vertices[mesh.faces],mesh.markers):
        rows.append((int(marker),tuple(sorted(map(tuple,triangle)))))
    return sorted(rows)


@pytest.mark.parametrize('entry',['builder','city'])
def test_fresh_public_default_and_zero_optional_points(entry):
    outputs=[]
    for options in (CompactSurfaceOptions(budgeted=False),None):
        city=sample();audit={}
        kwargs=dict(lod='lod2',mesher='dtcc_mesher',compact_surface=options,
                    report_mesh_quality=False,stage_audit=audit)
        mesh=builder.build_city_surface_mesh(city,**kwargs) if entry=='builder' else city.build_surface_mesh(**kwargs)
        compact=audit['roof_catalogue']['compact_surface'];outputs.append(compact)
        assert compact['topology']['valid'] and compact['shared_boundary_mismatches']==0
        assert audit['roof_catalogue']['retained_nonflat_fits']==1
        assert compact['building_triangles']<=3*compact['T0']//2
        assert not audit['roof_catalogue']['block_fallbacks']
        for group in compact['optimization']['groups']:
            for stage in ('A','B'):
                if stage in group:
                    x=group[stage]
                    assert x['end_triangles']<=x['triangle_ceiling']
                    assert max(x['source_charges'].values())<=(8 if stage=='A' else 4)
                    if stage=='B':
                        assert x['interior_points']==0
                        assert x['end_triangles']-x['start_triangles']<=x['start_triangles']//5
        assert all(r.attributes.get('source_key') is not None for r in mesh.regions if r.attributes['kind']!='terrain')
    assert not outputs[0]['optional_vertices']
    assert outputs[0]['T0']==outputs[1]['T0']


def test_auditing_and_merged_separate_outputs_preserve_geometry_and_markers():
    results=[];audits=[]
    for detailed,separate,external_audit in ((False,False,True),(True,False,True),
                                           (False,True,True),(False,False,False)):
        a={};mesh=builder.build_city_surface_mesh(sample(),lod='lod2',mesher='dtcc_mesher',
            compact_surface=CompactSurfaceOptions(detailed_audit=detailed),merge_meshes=not separate,
            stage_audit=a if external_audit else None,report_mesh_quality=False)
        if not external_audit:a=mesh.stage_audit
        if separate:mesh=_welded_mesh(mesh)
        results.append(signature(mesh));audits.append(a['roof_catalogue']['compact_surface'])
    assert results[0]==results[1]==results[2]==results[3]
    assert 'edges' not in audits[0]['graph'] and 'edges' in audits[1]['graph']
    assert not any(s['initial_edge_reproducers'] for s in audits[0]['optimization']['stages'])
    assert audits[0]['optional_counts']==audits[1]['optional_counts']


def test_backend_and_capability_errors_precede_city_work(monkeypatch):
    import dtcc_mesher as dm
    meshers=importlib.import_module('dtcc_core.builder.meshing.backends')
    monkeypatch.setattr(meshers,'available_2d_meshers',lambda:('triangle','dtcc_mesher'))
    with pytest.raises(ValueError,match='resolved mesher'):
        builder.build_city_surface_mesh(None,lod='lod2',mesher='triangle',compact_surface=CompactSurfaceOptions())
    backend=importlib.import_module('dtcc_core.builder.meshing.dtcc_mesher_backend')
    monkeypatch.setattr(backend,'_load_dtcc_mesher',lambda:SimpleNamespace())
    with pytest.raises(RuntimeError,match='bounded segment_trials'):
        builder.build_city_surface_mesh(None,lod='lod2',mesher='dtcc_mesher')
    # A new Python dataclass default must not hide an older native payload.
    monkeypatch.setattr(backend,'_load_dtcc_mesher',lambda:SimpleNamespace(Mesh=dm.Mesh,
        improve=lambda *a,**kw:SimpleNamespace(segment_trials=()),
        _core=SimpleNamespace(_improve_raw=lambda *a:{})))
    with pytest.raises(RuntimeError,match='bounded segment_trials'):
        builder.build_city_surface_mesh(None,lod='lod2',mesher='dtcc_mesher')
    monkeypatch.setattr(backend,'_load_dtcc_mesher',lambda:SimpleNamespace(Mesh=dm.Mesh,improve=lambda *a,**kw:SimpleNamespace()))
    with pytest.raises(RuntimeError,match='bounded segment_trials'):
        builder.build_city_surface_mesh(None,lod='lod2',mesher='dtcc_mesher')


def test_explicit_refinement_conflicts_and_extension_bounds():
    from dtcc_core.builder.meshing.surface_refinement import SurfaceRefinementOptions
    with pytest.raises(ValueError,match='excludes surface_refinement'):
        builder.build_city_surface_mesh(None,lod='lod2',compact_surface=CompactSurfaceOptions(),
                                        surface_refinement=SurfaceRefinementOptions())
    for key,value in [('extension_source_vertices',9),('extension_group_vertices',17),('extension_source_vertices',True)]:
        with pytest.raises(ValueError):CompactSurfaceOptions(**{key:value})


def two_groups():
    a,ca,fa=fixture();b,cb,fb=fixture()
    for i in range(len(b.vertices)):b.vertices[i]=b.vertices[i]+[100,0,0]
    offset=len(a.vertices)
    for f in fb:
        f.id+=2;f.source_key='c' if f.id==2 else 'd';f.origin=f.origin+[100,0,0]
        f.rings=[[i+offset for i in ring] for ring in f.rings]
        f.boundary=[(i+offset,j+offset) for i,j in f.boundary]
    vertices=np.vstack([a.vertices,b.vertices]);faces=fa+fb;caches=dict(ca)
    for i,c in cb.items():caches[i+2]=dict(c,ids=[j+offset for j in c['ids']])
    boundaries={f.id:f.boundary for f in faces}
    # A common terrain interface is deliberately incident to both components.
    boundaries[4]=[(u,v) for f in (faces[0],faces[2]) for u,v in f.boundary
                   if vertices[u,2]==vertices[v,2]==0]
    faces.append(SimpleNamespace(id=4,kind='terrain',source_key=None))
    return SharedEdges(vertices,boundaries),caches,faces


def test_independent_group_budgets_and_order_do_not_starve_distant_group():
    results=[]
    for reverse,alone in ((False,False),(True,False),(False,True)):
        registry,caches,faces=two_groups();groups=building_groups(registry,faces)
        assert [g['sources'] for g in groups]==[['a','b'],['c','d']]
        if alone:groups=groups[:1]
        if reverse:groups.reverse()
        charges={}
        for g in groups:
            start=sum(len(caches[i]['triangles']) for i in g['face_ids'])
            registry,a=improve_shared(registry,caches,faces,CompactSurfaceOptions(max_face_vertices=0,
                max_source_vertices=1,max_group_vertices=1),face_ids=set(g['face_ids']),max_building_triangles=start+2)
            assert len(a['optional_vertices'])==1
            charges.update(a['source_optional_vertices'])
        assert charges['a']==charges['b']==1
        if not alone:assert charges['c']==charges['d']==1
        result=[]
        for i in (0,1):
            c=caches[i]
            result+=sorted(tuple(sorted(map(tuple,np.array(registry.vertices)[np.array(c['ids'])[t]]))) for t in c['triangles'])
        results.append(sorted(result))
    assert results[0]==results[1]==results[2]


def test_stage_a_baseline_identical_under_standard_and_doubled_allowances():
    histories=[]
    for source,group in ((0,0),(4,8),(8,16)):
        a={}
        builder.build_city_surface_mesh(sample(),lod='lod2',mesher='dtcc_mesher',
            compact_surface=CompactSurfaceOptions(extension_source_vertices=source,extension_group_vertices=group),
            report_mesh_quality=False,stage_audit=a)
        compact=a['roof_catalogue']['compact_surface']
        histories.append([(p['kind'],p['xyz'],p['triangle_cost']) for p in compact['optional_vertices'] if p['stage']=='A'])
    assert histories[0]==histories[1]==histories[2]


def test_group_order_and_distant_group_on_fixed_reference_shell(monkeypatch):
    from .test_roof_catalogue import assemble
    from dtcc_core.builder.roofs import RoofPatch
    from dtcc_core.builder.meshing.compact_surface import compact_surface
    ci=importlib.import_module('dtcc_core.builder.meshing.compact_improvement')
    patches=[RoofPatch(box(-18,-5,-3,5),np.zeros(2),np.array([0.,0.,8.]),'a'),
             RoofPatch(box(8,-5,28,5),np.zeros(2),np.array([0.,0.,8.]),'b',block_id=1)]
    reference,patches=assemble(patches,slope=.4)
    original=ci.building_groups
    results=[]
    policy=CompactSurfaceOptions(max_source_vertices=1,max_group_vertices=1,
                                 extension_source_vertices=0,extension_group_vertices=0)
    for mode in ('normal','reverse','first'):
        def selected(registry,faces):
            groups=original(registry,faces)
            assert len(groups)==2
            return groups[::-1] if mode=='reverse' else groups[:1] if mode=='first' else groups
        monkeypatch.setattr(ci,'building_groups',selected)
        mesh,audit=compact_surface(reference,patches,policy)
        assert audit['topology']['valid'] and audit['terrain_locality']['outside_star_bitwise_preserved']
        assert all(g.get('A',{}).get('source_charges',{}).get('a',0)<=1 for g in audit['optimization']['groups'])
        selected_mesh=deepcopy(mesh)
        selected_mesh.faces=mesh.faces[np.isin(mesh.markers,[0,len(patches)])]
        selected_mesh.markers=mesh.markers[np.isin(mesh.markers,[0,len(patches)])]
        results.append((signature(selected_mesh),signature(SimpleNamespace(vertices=mesh.vertices,
            faces=mesh.faces[mesh.markers<0],markers=mesh.markers[mesh.markers<0]))))
    assert results[0]==results[1]==results[2]
