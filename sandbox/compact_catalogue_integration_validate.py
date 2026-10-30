"""Validate the public case matrix and its geometry/dependency provenance."""
import ast
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
from dtcc_core import io
import dtcc_mesher._core as native
from dtcc_core.builder import _dtcc_builder
from compact_catalogue_integration import ROOT,OUT,VARIANTS,digest


def signature(mesh,mask):
    return sorted(tuple(sorted(map(tuple,t))) for t in mesh.vertices[mesh.faces[mask]])


def point_history(a,stage):
    return [(p['kind'],p['xyz'],p['affected_faces'],p['charged_sources'],p['triangle_cost'])
            for p in a['optional_vertices'] if p['stage']==stage]


def main(input_dir=OUT/'public',detailed_dir=OUT/'public-detailed'):
    baseline=json.loads((OUT/'baseline-hashes.json').read_text())
    mismatches=[p for p,h in baseline.items() if digest(p)!=h]
    assert not mismatches,mismatches
    manifest=json.loads((OUT/'case-manifest.json').read_text())
    loaded={m.__file__:digest(m.__file__) for m in (native,_dtcc_builder)}
    outcomes={};success=0;lod2=0;failures=[];double_changes=[]
    for case in manifest['cases']:
        i=case['case'];records={};audits={};models={};detail={};checks=[]
        for p,h in case['inputs'].items():assert digest(p)==h,p
        for v in VARIANTS:
            p=input_dir/f'case-{i:02d}-{v}'
            records[v]=r=json.loads(p.with_suffix('.json').read_text())
            assert r['manifest_sha256']==digest(OUT/'case-manifest.json')
            assert r['native_binaries']==loaded
            assert r['peak_rss_bytes']<=2*1024**3
            if not r['valid_compact']:
                assert r['error'];failures.append(dict(case=i,variant=v,error=r['error']))
                continue
            models[v]=m=io.load_mesh(p.with_suffix('.dtcc'))
            audits[v]=a=r['audit']['roof_catalogue']['compact_surface']
            dpath=detailed_dir/f'case-{i:02d}-{v}'
            d=json.loads(dpath.with_suffix('.json').read_text())
            assert d['valid_compact'] and d['native_binaries']==loaded
            dm=io.load_mesh(dpath.with_suffix('.dtcc'))
            assert np.array_equal(m.vertices,dm.vertices),(i,v,'auditing changed vertices')
            assert np.array_equal(m.faces,dm.faces),(i,v,'auditing changed triangles')
            assert np.array_equal(m.markers,dm.markers)
            detail[v]=d['audit']['roof_catalogue']['compact_surface']
            assert a['topology']['valid'] and a['shared_boundary_mismatches']==0
            assert a['terrain_locality']['outside_star_bitwise_preserved']
            assert a['building_triangles']==np.count_nonzero(m.markers>=0)
            assert a['building_triangles']<=3*a['T0']//2
            assert a['building_triangles']-a['TB']<=a['TB']//5
            for f in a['faces']:
                g=f['geometry'];tol=max(1e-7,f['area_m2']*1e-9)
                assert g['plane_error_m']<=1e-7 and g['coverage_error_m2']<=tol and g['overlap_m2']<=tol
            for g in a['optimization']['groups']:
                assert all(source in g['sources'] for p in a['optional_vertices'] if p['group_id']==g['group_id'] for source in p['charged_sources'])
                assert g['sources']
                assert g.get('B',g.get('A',{})).get('end_triangles',g['T0'])<=3*g['T0']//2
                for stage in ('A','B'):
                    points=[p for p in a['optional_vertices'] if p['stage']==stage and p['group_id']==g['group_id']]
                    if stage not in g:
                        assert not points;continue
                    entry=g[stage];charges=Counter(s for p in points for s in p['charged_sources'])
                    source_cap=(8 if stage=='A' or v=='doubled' else 4)
                    group_cap=(32 if stage=='A' else 16 if v=='doubled' else 8)
                    assert len(points)<=group_cap and all(n<=source_cap for n in charges.values())
                    assert dict(charges)=={s:n for s,n in entry['source_charges'].items() if n}
                    assert entry['end_triangles']<=entry['triangle_ceiling']
                    if stage=='B':
                        assert all(p['kind']=='segment' for p in points)
                        assert entry['end_triangles']-entry['start_triangles']<=entry['start_triangles']//5
                    assert entry['end_triangles']-entry['start_triangles']==sum(p['triangle_cost'] for p in points)
            interior=Counter(i for p in a['optional_vertices'] if p['kind']=='interior' for i in p['affected_faces'])
            assert all(n<=2 for n in interior.values())
            assert len({p['vertex'] for p in a['optional_vertices']})==len(a['optional_vertices'])
            graph=detail[v]['graph']
            for p in a['optional_vertices']:
                assert p['unique_point_cost']==1 and p['triangle_cost']>0
                assert p['before']['penalty']-p['after']['penalty']>.015*p['triangle_cost']
                if p['kind']=='segment':
                    edge=graph['edges'][p['edge_id']]
                    x,y=np.array(edge['xyz'])[[0,-1]]
                    assert np.array_equal(p['xyz'],x+p['parameter']*(y-x))
                    assert p['affected_faces']==edge['incident_faces']
                pairs=([(old,p['per_face_after'][fi]) for fi,old in p['per_face_before'].items()]
                       if p['kind']=='segment' else [(p['before'],p['after'])])
                for old,new in pairs:
                    assert new['minimum_q']>=old['minimum_q']-1e-10
                    for angle in (5,10):
                        assert new[f'poor_{angle}_count']<=old[f'poor_{angle}_count']
                        assert new[f'poor_{angle}_area']<=old[f'poor_{angle}_area']+1e-10*max(1,old[f'poor_{angle}_area'])
            for stage in a['optimization']['stages']:
                assert stage.get('native_requests',0)<=512 and stage.get('joint_trials',0)<=256
            checks.append(v)
        if len(models)==4:
            success+=1
            lod2+=int(all(r['lod2_validated'] for r in records.values()))
            zero=models['zero'];basecat=records['zero']['audit']['roof_catalogue']
            if i<13:
                fields=('source_key','shape','outcome','final_outcome','fallback_area')
                prior=[tuple(s.get(k) for k in fields) for s in case['historical_selected_sources']]
                current=[tuple(s.get(k) for k in fields) for s in basecat['sources']]
                assert current==prior,(i,'historical fit/fallback decisions changed')
            assert set(map(tuple,zero.vertices))<=set(map(tuple,models['A'].vertices))
            for v,m in models.items():
                assert signature(zero,zero.markers<0)==signature(m,m.markers<0),(i,v,'terrain changed')
                cat=records[v]['audit']['roof_catalogue']
                for key in ('patch_markers','blocks','block_fallbacks'):assert cat[key]==basecat[key],(i,v,key)
                fields=lambda mesh:{r.id:r.attributes for r in mesh.regions}
                assert fields(m)==fields(zero)
            assert not audits['zero']['optional_vertices']
            for v in ('standard','doubled'):
                assert audits[v]['TB']==audits['A']['building_triangles']
                assert point_history(audits[v],'A')==point_history(audits['A'],'A')
                assert set(map(tuple,models['A'].vertices))<=set(map(tuple,models[v].vertices))
            if signature(models['standard'],models['standard'].markers>=0)!=signature(models['doubled'],models['doubled'].markers>=0):
                double_changes.append(i)
        outcomes[str(i)]=dict(valid_controls=checks,lod2_validated=all(r['lod2_validated'] for r in records.values()),
                              auditing_bitwise_equal=True if checks else None)
    assert success>=10 and lod2>=10,(success,lod2)
    untouched=[]
    for root,name,paths in [(ROOT,'core',['pyproject.toml','uv.lock',
        'dtcc_core/builder/meshing/compact_geometry.py','dtcc_core/builder/meshing/shared_edges.py',
        'dtcc_core/builder/meshing/roof_regions.py','dtcc_core/builder/meshing/roof_interfaces.py',
        'dtcc_core/cpp/include/MeshBuilder.h','dtcc_core/cpp/include/Triangulate.h']),
        (ROOT.parent/'dtcc-mesher','mesher',['include/dtcc_mesher/dtcc_mesher.h','src/core/mesh.c','src/core/mesh.h',
         'tests/test_mesh.c','tests/python/test_corner_reuse.py'])]:
        if name=='core':paths += [str(p.relative_to(ROOT)) for p in (ROOT/'dtcc_core/builder/roofs').rglob('*.py')]
        for path in paths:
            assert digest(root/path)==digest(OUT/'starting-source'/name/path),path
            untouched.append(str(root/path))
    for path in (ROOT/'dtcc_core/builder/meshing').glob('compact*.py'):
        tree=ast.parse(path.read_text())
        for node in ast.walk(tree):
            if isinstance(node,ast.ImportFrom):assert not (node.module or '').startswith(('sandbox','unittest','pickle'))
            elif isinstance(node,ast.Import):assert all(not item.name.startswith(('sandbox','unittest','pickle')) for item in node.names)
    logs={}
    for path,expected in [('builder-tests-complete.log','761 passed, 2 skipped'),('mesher-python-final.log','62 passed'),
                          ('c-tests-final.log','100% tests passed, 0 tests failed out of 4'),
                          ('public-api-final.log','8 passed')]:
        assert expected in (OUT/path).read_text(),path
        logs[path]=dict(result=expected,sha256=digest(OUT/path))
    source_hashes={}
    for root in (ROOT,ROOT.parent/'dtcc-mesher'):
        names=subprocess.check_output(['git','ls-files','--cached','--others','--exclude-standard','-z'],cwd=root).decode().split('\0')
        for name in names:
            p=root/name
            if name and p.is_file() and p.suffix in ('.py','.c','.h','.cpp','.md','.toml','.json'):
                source_hashes[str(p)]=digest(p)
    previous_attempts={}
    for folder in ('matrix','final','detailed'):
        attempts=[]
        for file in sorted((OUT/folder).glob('case-*-*.json')):
            if file.name.endswith('-worker.json'):continue
            x=json.loads(file.read_text());attempts.append(dict(case=x['case'],variant=x['variant'],error=x['error'],valid_compact=x['valid_compact']))
        previous_attempts[folder]=dict(records=attempts,status='development control; not counted in settled validation')
    result=dict(valid=True,case_count=len(manifest['cases']),valid_compact_cases=success,lod2_cases=lod2,
        failures=failures,cases=outcomes,doubled_changes_mesh_cases=double_changes,
        normal_vs_detailed_audit='All valid control meshes bitwise equal in both full matrices',
        prior_baseline_files_unchanged=len(baseline),source_hashes=source_hashes,unrelated_work_preserved=untouched,
        historical_thirteen_fit_and_fallback_decisions_preserved=True,
        binaries=loaded,tests=logs,previous_attempts=previous_attempts,
        dependency_gate=dict(status='required Mesher changes are uncommitted; compatible release and portable pin update remain',
            base_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT.parent/'dtcc-mesher',text=True).strip(),
            fabricated_revision=False,portable_pin_and_lockfile_unchanged=True),
        recommendation='Standard additional 4/source and 8/group, doubled 8/16 explicit; identical triangle and severe-defect guards',
        scope='Representative fixed cases only; no full city, commit or publishing')
    (OUT/'validation.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({k:result[k] for k in ('valid','case_count','valid_compact_cases','lod2_cases','prior_baseline_files_unchanged','doubled_changes_mesh_cases')}))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir',type=Path,default=OUT/'public')
    parser.add_argument('--detailed-dir',type=Path,default=OUT/'public-detailed')
    args=parser.parse_args();main(args.input_dir,args.detailed_dir)
