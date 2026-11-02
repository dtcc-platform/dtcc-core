"""Verify fixed geometry, bounded transactions, binaries and preserved baselines."""
from collections import Counter
import hashlib
import json
from pathlib import Path
import pickle
import re
import subprocess
import numpy as np
from dtcc_core import io
import dtcc_mesher._core as native
from compact_edge_extension import ROOT,OUT
from roof_wall_evaluation import _triangle_signature


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    baseline=json.loads((OUT/'baseline-hashes.json').read_text())
    assert all(digest(p)==sha for p,sha in baseline.items())
    inherited=[]
    for manifest in (OUT.parent/'compact-mesher-interface/baseline-hashes.json',
                     OUT.parent/'shared-edge-pilot/baseline-manifest.json'):
        data=json.loads(manifest.read_text())
        if isinstance(data,dict):
            for path,sha in data.items():
                candidate=Path(path) if Path(path).is_absolute() else ROOT/path
                if isinstance(sha,str) and re.fullmatch('[0-9a-f]{64}',sha):
                    assert digest(candidate)==sha,(manifest,candidate)
                    inherited.append(str(candidate))
    binaries={str(Path(native.__file__)):digest(native.__file__)}
    symbols=subprocess.check_output(['nm','-g',native.__file__],text=True)
    for symbol in ('_dtcc_mesher_improve','_dtcc_mesher_improvement_init','_dtcc_mesher_improvement_free'):
        assert symbol in symbols
    cases={}
    for case in (1,11):
        state=pickle.load(open(OUT/f'case-{case:02d}-checkpoint.pkl','rb'))
        base=io.load_mesh(OUT.parent/f'compact-mesher-interface/case-{case:02d}-shared-dtcc_mesher.dtcc')
        rows=[]
        for budget in ([0,2,4] if case==1 else [0,2,4,8]):
            p=OUT/'joint-fallback'/f'case-{case:02d}-extra-{budget}'
            mesh=io.load_mesh(p.with_suffix('.dtcc'));record=json.loads(p.with_suffix('.json').read_text())
            a,e=record['audit'],record['extension'];points=e['optional_vertices']
            assert record['binary_sha256']==binaries[str(Path(native.__file__))]
            assert e['old_points_bitwise_fixed'] and e['old_points_used'] and e['no_new_interiors']
            assert a['topology']['valid'] and not a['shared_boundary_mismatches']
            assert a['terrain_locality']['outside_star_bitwise_preserved']
            assert record['terrain_bitwise_fixed']
            assert np.array_equal(_triangle_signature(mesh,mesh.markers<0),_triangle_signature(base,base.markers<0))
            assert set(map(tuple,base.vertices))<=set(map(tuple,mesh.vertices))
            assert len(points)<=budget and max(e['source_optional_vertices'].values())<=4
            assert np.count_nonzero(mesh.markers>=0)<= (85 if case==1 else 307)
            charges=Counter()
            for point in points:
                assert point['kind']=='segment' and point['triangle_cost']==2 and point['unique_point_cost']==1
                edge=state['registry'].edges[point['edge_id']]
                x,y=np.array(state['registry'].vertices)[list(edge.endpoints)]
                assert np.array_equal(point['xyz'],x+point['parameter']*(y-x))
                assert len(point['affected_faces'])==2 and all(state['faces'][i].kind!='terrain' for i in point['affected_faces'])
                charges.update(point['charged_sources'])
                for i,before in point['per_face_before'].items():
                    after=point['per_face_after'][i]
                    assert after['minimum_q']>=before['minimum_q']-1e-10
                    for angle in (5,10):
                        assert after[f'poor_{angle}_count']<=before[f'poor_{angle}_count']
                        assert after[f'poor_{angle}_area']<=before[f'poor_{angle}_area']+1e-10*max(1,before[f'poor_{angle}_area'])
                assert point['before']['penalty']-point['after']['penalty']>.015*point['triangle_cost']
            assert dict(charges)=={s:n for s,n in e['source_optional_vertices'].items() if n}
            assert e['native_requests']<=512 and e['joint_trials']<=256
            assert sum(p['triangle_cost'] for p in points)==np.count_nonzero(mesh.markers>=0)-np.count_nonzero(base.markers>=0)
            for face in a['faces']:
                g=face['geometry'];area=face['area_m2']
                assert g['plane_error_m']<=1e-7 and g['coverage_error_m2']<=max(1e-7,area*1e-9)
                assert g['overlap_m2']<=max(1e-7,area*1e-9)
            assert record['peak_rss_bytes']<=2*1024**3 and record['seconds']<600
            rows.append(dict(budget=budget,used=len(points),triangles=int(np.count_nonzero(mesh.markers>=0)),
                             geometry=True,topology=True,ownership=True,terrain_bitwise=True,all_old_points_fixed=True,
                             termination=e['termination'],charges=dict(charges)))
        cases[str(case)]=rows
    logs={}
    for file,count in [('c-tests.log',4),('mesher-python.log',61),('core-targeted-final.log',131),('builder-tests-final.log',753)]:
        text=(OUT/file).read_text()
        assert (f'{count} passed' in text if file!='c-tests.log' else '100% tests passed, 0 tests failed out of 4' in text)
        logs[file]=dict(sha256=digest(OUT/file),passed=count)
    hashes={}
    for repo in (ROOT,ROOT.parent/'dtcc-mesher'):
        names=subprocess.check_output(['git','ls-files','--cached','--others','--exclude-standard','-z'],cwd=repo).decode().split('\0')
        for name in names:
            p=repo/name
            if name and p.is_file() and p.suffix in ('.py','.c','.cpp','.h','.md','.json'):
                hashes[str(p)]=digest(p)
    untouched=[]
    for repo,paths in [(ROOT,['dtcc_core/builder/geometry_builders/meshes.py',
        'dtcc_core/builder/meshing/catalogue_surface.py','dtcc_core/builder/meshing/compact_geometry.py',
        'dtcc_core/builder/meshing/shared_edges.py','dtcc_core/builder/meshing/compact_triangulation.py',
        'dtcc_core/cpp/include/MeshBuilder.h','dtcc_core/cpp/include/Triangulate.h','pyproject.toml','uv.lock']),
        (ROOT.parent/'dtcc-mesher',['include/dtcc_mesher/dtcc_mesher.h','src/core/mesh.c','src/core/mesh.h',
         'src/core/cdt.c','src/core/cdt.h','tests/test_mesh.c','tests/python/test_corner_reuse.py'])]:
        for name in paths:
            snapshot=OUT/'starting-source'/('core' if repo==ROOT else 'mesher')/name
            assert digest(repo/name)==digest(snapshot)
            untouched.append(str(repo/name))
    assert 'Bitwise deterministic 1' in (OUT/'determinism.log').read_text()
    assert 'Bitwise deterministic 11' in (OUT/'determinism.log').read_text()
    result=dict(valid=True,cases=cases,native_binaries=binaries,native_symbols_verified=True,
        preserved_current_baseline_files=len(baseline),preserved_inherited_manifest_entries=len(inherited),
        frozen_starting_binary=json.loads((OUT/'starting-binaries.json').read_text()),source_hashes=hashes,tests=logs,
        unrelated_work_preserved=untouched,pilot_repeat_meshes_bitwise_deterministic=True,
        sanitizer='Not rerun: the previous host smoke tests stalled before main, including instrumented hello world; no sanitizer success claimed.',
        public_api='Same improve arguments; penalty now squared q-defect energy; result adds at most four segment provenance records, C consumers rebuild.',
        scope='Only case-01 and case-11; no commit, publishing, default enablement or geometric approximation')
    (OUT/'validation.json').write_text(json.dumps(result,indent=2))
    print('Validated both pilots;',len(baseline),'current baseline files unchanged;',len(inherited),'inherited hashes verified')


if __name__=='__main__':main()
