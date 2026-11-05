"""Check the fixed five-way public matrix, frozen meshes and cumulative policies."""
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
from compact_adaptive_evaluation import ROOT,OUT,VARIANTS,digest
from compact_catalogue_integration_validate import signature,point_history
from compact_adaptive_report import scores,arrays


def same(a,b):
    return all(np.array_equal(getattr(a,k),getattr(b,k)) for k in ('vertices','faces','markers'))


def main(input_dir,detailed_dir,repeat_dir):
    baseline=json.loads((OUT/'baseline-hashes.json').read_text())
    assert all(digest(p)==h for p,h in baseline.items()),'Historical baseline modified'
    manifest=json.loads((OUT/'case-manifest.json').read_text())
    assert digest(OUT/'case-manifest.json')==digest(OUT.parent/'compact-integration/case-manifest.json')
    loaded={m.__file__:digest(m.__file__) for m in (native,_dtcc_builder)}
    prior=json.loads((OUT.parent/'compact-integration/validation.json').read_text())
    assert loaded==prior['binaries']
    source_hashes={}
    # No fitting, geometry, terrain, native algorithm or dependency changes.
    untouched=[]
    for repo,name in ((ROOT,'core'),(ROOT.parent/'dtcc-mesher','mesher')):
        names=subprocess.check_output(['git','ls-files','--cached','--others','--exclude-standard','-z'],cwd=repo).decode().split('\0')
        for name0 in names:
            p=repo/name0
            if not name0 or not p.is_file() or p.suffix not in ('.py','.c','.h','.cpp','.md','.toml','.lock'):continue
            source_hashes[str(p)]=digest(p)
            before=OUT/'starting-source'/name/name0
            allowed=name=='core' and (name0 in ('dtcc_core/builder/meshing/compact_surface.py',
                'dtcc_core/builder/meshing/compact_improvement.py','docs/compact_mesher_interface.md') or
                name0.startswith(('tests/','sandbox/compact_adaptive')))
            if before.exists() and not allowed:
                assert digest(p)==digest(before),str(p)
                untouched.append(str(p))
    results={};failures=[];changes=[];double_equal=[];escalation_count=0;tradeoffs=[];fixed_equal=0;repeated=[]
    totals={v:Counter() for v in VARIANTS};lod2totals={v:Counter() for v in VARIANTS};valid=lod2=0
    for case in manifest['cases']:
        i=case['case'];models={};audits={};records={};qualities={};checks={}
        for p,h in case['inputs'].items():assert digest(p)==h
        for v in VARIANTS:
            prefix=input_dir/f'case-{i:02d}-{v}';r=records[v]=json.loads(prefix.with_suffix('.json').read_text())
            assert r['manifest_sha256']==digest(OUT/'case-manifest.json')
            assert r['native_binaries']==loaded
            if not r['valid_compact']:
                failures.append(dict(case=i,variant=v,error=r['error']));continue
            models[v]=m=io.load_mesh(prefix.with_suffix('.dtcc'));a=audits[v]=r['audit']['roof_catalogue']['compact_surface']
            assert a['topology']['valid'] and a['shared_boundary_mismatches']==0
            assert a['terrain_locality']['outside_star_bitwise_preserved']
            assert a['building_triangles']==sum(m.markers>=0)
            assert a['building_triangles']<=3*a['T0']//2
            assert a['building_triangles']-a['TB']<=a['TB']//5
            for f in a['faces']:
                g=f['geometry'];tol=max(1e-7,f['area_m2']*1e-9)
                assert g['plane_error_m']<=1e-7 and g['coverage_error_m2']<=tol and g['overlap_m2']<=tol
            n=len(r['audit']['roof_catalogue']['patch_markers']);q=qualities[v]=scores(m,n)['building']
            q['energy']=float(np.square(1-arrays(m,m.markers>=0)[3]).sum())
            for population in (totals, lod2totals if r['lod2_validated'] else None):
                if population is None:continue
                population[v].update(triangles=q['triangles'],energy=q['energy'],seconds=r['seconds'])
                for angle in (5,10,20):
                    population[v][f'poor_{angle}_count']+=q['poor'][str(angle)]['count']
                    population[v][f'poor_{angle}_area']+=q['poor'][str(angle)]['area_m2']
            for g in a['optimization']['groups']:
                for stage in ('A','B'):
                    entry=g[stage];points=[p for p in a['optional_vertices'] if p['stage']==stage and p['group_id']==g['group_id']]
                    charges=Counter(s for p in points for s in p['charged_sources'])
                    assert dict(charges)=={s:n for s,n in entry['source_charges'].items() if n}
                    assert entry['end_triangles']-entry['start_triangles']==sum(p['triangle_cost'] for p in points)
                    assert entry['end_triangles']<=entry['triangle_ceiling']
                    assert len(points)<=(32 if stage=='A' else 16 if v in ('adaptive','doubled') else 8)
                    assert all(n<=(8 if stage=='A' or v in ('adaptive','doubled') else 4) for n in charges.values())
                    if stage=='B':
                        tb=g['A']['end_triangles']
                        assert entry['start_triangles']==tb
                        assert entry['triangle_ceiling']==min(3*g['T0']//2,tb+tb//5)
                        assert all(p['kind']=='segment' for p in points)
            for p in a['optional_vertices']:
                assert p['unique_point_cost']==1 and p['triangle_cost']>0
                assert p['before']['penalty']-p['after']['penalty']>.015*p['triangle_cost']
                pairs=([(old,p['per_face_after'][fi]) for fi,old in p['per_face_before'].items()]
                    if p['kind']=='segment' else [(p['before'],p['after'])])
                for old,new in pairs:
                    assert new['minimum_q']>=old['minimum_q']-1e-10
                    for angle in (5,10):
                        assert new[f'poor_{angle}_count']<=old[f'poor_{angle}_count']
                        assert new[f'poor_{angle}_area']<=old[f'poor_{angle}_area']+1e-10*max(1,old[f'poor_{angle}_area'])
                if p['after']['poor_20_count']>p['before']['poor_20_count']:
                    tradeoffs.append(dict(case=i,variant=v,stage=p['stage'],phase=p['phase'],vertex=p['vertex'],
                        metric='incident below-20 count',before=p['before']['poor_20_count'],after=p['after']['poor_20_count']))
            for s in a['optimization']['stages']:
                assert s.get('native_requests',0)<=512 and s.get('joint_trials',0)<=256
                stop=s['stopping'];assert stop['primary_reason']
                if stop['primary_reason']=='disabled':assert not stop['evaluated_finite_set']
                if stop['primary_reason'] in ('request_budget','native_request_limit','joint_trial_limit','runtime_limit'):
                    assert not stop['evaluated_finite_set']
                if stop.get('assessment',{}).get('triangle_blocked'):
                    assert 'building_triangle_ceiling' in stop['binding_limits']
            if v!='adaptive':
                frozen=io.load_mesh(OUT.parent/f'compact-integration/public/case-{i:02d}-{v}.dtcc')
                assert same(m,frozen),(i,v,'frozen fixed mesh changed')
                fixed_equal+=1
            else:
                d=json.loads((detailed_dir/f'case-{i:02d}-adaptive.json').read_text())
                assert d['valid_compact'] and d['native_binaries']==loaded
                dm=io.load_mesh(detailed_dir/f'case-{i:02d}-adaptive.dtcc')
                assert same(m,dm),(i,'detailed audit changed mesh')
                graph=d['audit']['roof_catalogue']['compact_surface']['graph']
                for p in a['optional_vertices']:
                    if p['kind']!='segment':continue
                    edge=graph['edges'][p['edge_id']];x,y=np.array(edge['xyz'])[[0,-1]]
                    assert np.array_equal(p['xyz'],x+p['parameter']*(y-x))
                    assert p['affected_faces']==edge['incident_faces']
                if (repeat_dir/f'case-{i:02d}-adaptive.dtcc').exists():
                    repeated.append(i);assert same(m,io.load_mesh(repeat_dir/f'case-{i:02d}-adaptive.dtcc'))
            checks[v]=True
        if len(models)!=5:
            results[str(i)]=dict(valid_controls=checks,errors=[x for x in failures if x['case']==i]);continue
        valid+=1;lod2+=int(all(r['lod2_validated'] for r in records.values()))
        z=models['zero'];base=records['zero']['audit']['roof_catalogue']
        for v,m in models.items():
            assert signature(z,z.markers<0)==signature(m,m.markers<0),(i,v,'terrain changed')
            cat=records[v]['audit']['roof_catalogue']
            for k in ('sources','blocks','block_fallbacks','patch_markers'):assert cat[k]==base[k],(i,v,k)
            assert {r.id:r.attributes for r in m.regions}=={r.id:r.attributes for r in z.regions}
        assert set(map(tuple,models['zero'].vertices))<=set(map(tuple,models['A'].vertices))
        for v in ('standard','doubled','adaptive'):
            assert set(map(tuple,models['A'].vertices))<=set(map(tuple,models[v].vertices))
            assert point_history(audits[v],'A')==point_history(audits['A'],'A')
            assert audits[v]['TB']==audits['A']['building_triangles']
        adaptive=audits['adaptive'];std=audits['standard']
        for s in adaptive['optimization']['stages']:
            if s['stage']!='B':continue
            fixed=next(x for x in std['optimization']['stages'] if x['stage']=='B' and x['group_id']==s['group_id'])
            stdpts=fixed['optional_vertices'];stdresult=s['standard_result']
            assert s['optional_vertices'][:len(stdpts)]==stdpts,(i,'standard prefix changed')
            assert stdresult['optional_points']==len(stdpts)
            assert stdresult['building_triangles']==fixed['building_triangles']
            assert stdresult['source_charges']==fixed['source_optional_vertices']
            assert stdresult['native_requests']==fixed['native_requests'] and stdresult['joint_trials']==fixed['joint_trials']
            expected_source=dict.fromkeys(s['source_optional_vertices'],4);expected_group=8
            for e in s['escalations']:
                p=s['optional_vertices'][e['point_index']]
                assert p['phase']=='adaptive' and p['gain_per_triangle']>.015
                for source,change in e['sources'].items():
                    assert change['before']==expected_source[source] and change['after']==change['before']+1
                    assert source in p['charged_sources'];expected_source[source]=change['after']
                assert e['group_before']==expected_group and expected_group<=e['group_after']<=expected_group+1
                expected_group=e['group_after']
            assert s['effective_allowances']==dict(source=expected_source,group=expected_group)
            escalation_count+=len(s['escalations'])
        if not same(models['standard'],models['adaptive']):changes.append(i)
        if same(models['doubled'],models['adaptive']):double_equal.append(i)
        for part in ('building',):
            for metric in ('maximum_vertex_incidence',):
                if qualities['adaptive'][metric]>qualities['standard'][metric]:
                    tradeoffs.append(dict(case=i,metric=metric,before=qualities['standard'][metric],after=qualities['adaptive'][metric]))
            for angle in (5,10,20):
                for key in ('count','area_m2'):
                    x=qualities['standard']['poor'][str(angle)][key];y=qualities['adaptive']['poor'][str(angle)][key]
                    if y>x+1e-8:tradeoffs.append(dict(case=i,metric=f'below-{angle} {key}',before=x,after=y))
        results[str(i)]=dict(valid_controls=checks,lod2_validated=records['adaptive']['lod2_validated'],
            adaptive_changed=i in changes,adaptive_equals_doubled=i in double_equal,
            fixed_standard_prefix_equal=True,detailed_bitwise_equal=True,
            quality=qualities,adaptive_policy=adaptive['optimization'])
    assert not failures,failures
    assert valid==29 and lod2==18
    recommendation=('Adaptive continuation retains the standard prefix and accepts '+
        f'{escalation_count} guarded extra points in cases {changes}, charging only the sources/groups that need them. '+
        'Keep fixed standard for the smallest meshes and fixed doubled as an explicit comparison. '+
        'The fully retained LOD2 cases are unchanged; additional gains here concern existing fallback/mixed-fit controls. '+
        'No general optimality is claimed.')
    if len(double_equal)==29:
        recommendation=('Keep fixed standard 4/8 as the normal default and fixed doubled 8/16 as the simpler explicit policy for measured extra quality. '+
            f'Adaptive accepts {escalation_count} guarded extra points in cases {changes}, but its final meshes equal fixed doubled on all 29 cases. '+
            'The fully retained LOD2 cases remain unchanged. Adaptive remains opt-in for its standard-prefix and per-source allowance accounting; this sample does not demonstrate a final-mesh advantage over fixed doubled.')
    policy=dict(recommendation=recommendation,adaptive_changed_cases=changes,
        previously_doubled_helped_cases=[4,6,10,11,12],adaptive_equals_doubled_cases=double_equal,
        escalations=escalation_count,totals={k:dict(v) for k,v in totals.items()},
        retained_lod2_totals={k:dict(v) for k,v in lod2totals.items()},tradeoffs=tradeoffs,cases=results)
    (OUT/'policy-summary.json').write_text(json.dumps(policy,indent=2))
    attempts={}
    for folder in ('public','final','settled','detailed','repeat'):
        rows=[]
        for p in sorted((OUT/folder).glob('case-*-*.json')):
            if p.stem.endswith('worker'):continue
            r=json.loads(p.read_text());rows.append(dict(case=r['case'],variant=r['variant'],valid_compact=r['valid_compact'],error=r['error']))
        attempts[folder]=rows
    required_tests={'builder-verified-tests.log':'783 passed, 2 skipped',
        'adaptive-tests.log':'22 passed','mesher-python.log':'62 passed',
        'c-tests.log':'100% tests passed, 0 tests failed out of 4'}
    for file,expected in required_tests.items():assert expected in (OUT/file).read_text(),file
    assert subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT.parent/'dtcc-mesher',text=True).strip()==prior['dependency_gate']['base_commit']
    result=dict(valid=True,case_count=29,valid_compact_cases=valid,lod2_cases=lod2,failures=failures,
        preserved_baseline_files=len(baseline),frozen_fixed_meshes_bitwise_equal=fixed_equal,
        adaptive_detailed_meshes_bitwise_equal=29,repeated_deterministic_cases=repeated,
        source_hashes=source_hashes,unrelated_sources_unchanged=untouched,binaries=loaded,
        native_source_provenance='Loaded binaries match frozen verified build; all Mesher/native source unchanged; C/CLI rebuilt and tested from current sources',
        tests={p.name:dict(sha256=digest(p),summary=next((line for line in p.read_text().splitlines() if ' passed' in line or '100% tests passed' in line),None)) for p in OUT.glob('*.log') if 'passed' in p.read_text()},
        attempted_cases=attempts,
        required_tests=required_tests,
        artifacts={str(p.relative_to(OUT)):digest(p) for p in [OUT/'report.html',OUT/'api-rationale.md',OUT/'case-manifest.json',*OUT.glob('case-*.png')]},
        dependency_gate=prior['dependency_gate'],recommendation=recommendation)
    (OUT/'validation.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({k:result[k] for k in ('valid','case_count','lod2_cases','frozen_fixed_meshes_bitwise_equal','repeated_deterministic_cases')}))
    print(recommendation)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir',type=Path,default=OUT/'settled')
    parser.add_argument('--detailed-dir',type=Path,default=OUT/'detailed')
    parser.add_argument('--repeat-dir',type=Path,default=OUT/'repeat')
    a=parser.parse_args();main(a.input_dir,a.detailed_dir,a.repeat_dir)
