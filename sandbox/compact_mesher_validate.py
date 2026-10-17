"""Verify bounded native/shared pilot outputs and preservation manifests."""
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'.agent/work/compact-mesher-interface';MESHER=ROOT.parent/'dtcc-mesher'


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    baseline=json.loads((OUT/'baseline-hashes.json').read_text())
    assert all(digest(Path(p))==h for p,h in baseline.items())
    previous=json.loads((ROOT/'.agent/work/shared-edge-pilot/baseline-manifest.json').read_text())
    assert all(digest(ROOT/p)==h for p,h in previous.items())
    starting=json.loads((OUT/'starting-source-hashes.json').read_text())
    untouched=[MESHER/p for p in ['include/dtcc_mesher/dtcc_mesher.h','src/core/mesh.c','src/core/mesh.h','tests/test_mesh.c','tests/python/test_corner_reuse.py']]
    untouched += [ROOT/p for p in ['dtcc_core/builder/meshing/catalogue_surface.py','dtcc_core/builder/geometry_builders/meshes.py',
                                  'dtcc_core/builder/meshing/roof_regions.py','dtcc_core/builder/meshing/roof_interfaces.py']]
    assert all(digest(p)==starting[str(p)] for p in untouched)
    output=dict(frozen_compact_files=len(baseline),previous_baseline_files=len(previous),baseline_mismatches=[],
                preserved_unrelated_experiment_files=[str(p) for p in untouched],cases={})
    for case in (1,11):
        variants={}
        for variant in ('zero','interior','shared'):
            path=OUT/f'case-{case:02d}-{variant}-dtcc_mesher.json';r=json.loads(path.read_text())
            assert r['error'] is None and all(r['preservation'].values())
            for p,h in r['inputs'].items():assert digest(ROOT/p)==h,p
            for p,h in r['source_hashes'].items():assert digest(ROOT/p)==h,p
            for p,h in r['mesher_source_hashes'].items():assert digest(Path(p))==h,p
            for p,h in r['native_binaries'].items():assert digest(Path(p))==h,p
            a=r['audit']['roof_catalogue']['compact_surface'];assert a['topology']['valid'] and a['shared_boundary_mismatches']==0
            assert a['terrain_locality']['outside_star_bitwise_preserved']
            assert len(a['optional_vertices'])<=32
            assert len({p['vertex'] for p in a['optional_vertices']})==len(a['optional_vertices'])
            used=dict.fromkeys(a['source_optional_vertices'],0)
            for point in a['optional_vertices']:
                assert point['unique_point_cost']==1 and point['triangle_cost']==2
                assert point['before']['penalty']-point['after']['penalty']>.015*point['triangle_cost']
                for s in point['charged_sources']:used[s]+=1
            assert used==a['source_optional_vertices'] and all(n<=8 for n in used.values())
            for f in a['faces']:
                assert f['geometry']['plane_error_m']<=1e-7
                assert f['geometry']['coverage_error_m2']<=max(1e-7,f['area_m2']*1e-9)
                interior=sum(p['kind']=='interior' for p in f['optional_vertices']);boundary=sum(p['kind']=='segment' for p in f['optional_vertices'])
                assert interior<=2
                assert f['triangles']==f['zero_steiner_triangles']+2*interior+boundary
            if variant=='zero':assert not a['optional_vertices']
            if variant=='interior':assert all(p['kind']=='interior' for p in a['optional_vertices'])
            resource=json.loads(path.with_name(path.stem+'-resource.json').read_text())
            assert resource['status']=='completed' and resource['within_memory_budget']
            variants[variant]=dict(optional_points=len(a['optional_vertices']),edge_points=sum(p['kind']=='segment' for p in a['optional_vertices']),
                source_charges=used,topology=a['topology'],terrain_locality=a['terrain_locality'],resource=resource)
        comparison=json.loads((OUT/f'case-{case:02d}-comparison.json').read_text())
        assert comparison['terrain_triangles_bitwise_equal_across_new_variants']
        for variant in variants:assert comparison['records'][variant]==json.loads((OUT/f'case-{case:02d}-{variant}-dtcc_mesher.json').read_text())
        output['cases'][str(case)]=variants
    evidence=json.loads((OUT/'case-01-edge-evidence.json').read_text())
    assert len(evidence['trials'])==91 and not any(p['joint']['accepted'] for p in evidence['trials'])
    expected={'core-targeted-final.log':'126 passed','builder-tests-final.log':'748 passed, 2 skipped',
              'mesher-python-tests.log':'58 passed','c-tests.log':'100% tests passed, 0 tests failed out of 4'}
    for p,text in expected.items():assert text in (OUT/p).read_text(),p
    output['tests']=expected
    output['optional_sanitizer']=json.loads((OUT/'sanitizer-environment.json').read_text())
    output['optional_sanitizer']['status']='unavailable: instrumented hello world stalls before main'
    code=[*sorted((ROOT/'dtcc_core/builder/meshing').glob('compact*.py')),
          ROOT/'tests/builder/test_compact_surface.py',ROOT/'tests/builder/test_compact_mesher_interface.py',
          *sorted((ROOT/'sandbox').glob('compact_*.py')),ROOT/'sandbox/compact_edge_probe.c',
          MESHER/'include/dtcc_mesher/dtcc_mesher_improve.h',MESHER/'src/api/improve.c',
          MESHER/'python/bindings.cpp',MESHER/'python/dtcc_mesher/__init__.py',
          MESHER/'tests/c/test_improve.c',MESHER/'tests/python/test_improve.py',
          ROOT/'docs/compact_mesher_interface.md',MESHER/'docs/improve.md']
    output['source_manifest']={str(p):digest(p) for p in code}
    import dtcc_mesher,dtcc_mesher._core as native
    output['loaded_mesher']=dict(version=dtcc_mesher.__version__,python=dtcc_mesher.__file__,binary=native.__file__,sha256=digest(Path(native.__file__)),
                                 starting_binary_sha256=digest(OUT/'starting-mesher.so'))
    assert output['loaded_mesher']['sha256']!=output['loaded_mesher']['starting_binary_sha256']
    output['dependency_pin_and_lockfile_changes']=False
    (OUT/'validation.json').write_text(json.dumps(output,indent=2))
    print(json.dumps(dict(valid=True,frozen_compact_files=len(baseline),previous_baseline_files=len(previous),tests=expected)))

if __name__=='__main__':main()
