"""Verify the two saved compact cases, test logs and preserved baseline hashes."""
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
WORK=ROOT/'.agent/work'
OUT=WORK/'minimal-catalogue'


def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    previous=json.loads((WORK/'shared-edge-pilot/baseline-manifest.json').read_text())
    mismatches=[p for p,h in previous.items() if not (ROOT/p).exists() or digest(ROOT/p)!=h]
    assert not mismatches,mismatches
    broad=json.loads((OUT/'baseline-hashes.json').read_text())
    broad_mismatches=[p for p,h in broad.items() if not (ROOT/p).exists() or digest(ROOT/p)!=h]
    assert not broad_mismatches,broad_mismatches
    result=dict(previous_manifest_files=len(previous),previous_manifest_mismatches=mismatches,
                baseline_files=len(broad),baseline_mismatches=broad_mismatches,cases={})
    for case in [1,11]:
        entries={}
        for variant,backend in [('zero','dtcc_mesher'),('budgeted','dtcc_mesher'),('zero','geos'),('zero','triangle')]:
            path=OUT/f'case-{case:02d}-{variant}-{backend}.json'
            r=json.loads(path.read_text());assert r['error'] is None,r['error']
            assert all(r['preservation'].values()),r['preservation']
            a=r['audit']['roof_catalogue']['compact_surface']
            assert a['accepted'] and a['topology']['valid']
            assert a['shared_boundary_mismatches']==0
            assert a['terrain_locality']['outside_star_bitwise_preserved']
            assert len(a['optional_vertices'])<=32
            assert all(n<=8 for n in a['source_optional_vertices'].values())
            assert all(len(f['optional_vertices'])<=2 for f in a['faces'])
            assert not a['optional_vertices'] if variant=='zero' else True
            for f in a['faces']:
                assert f['geometry']['plane_error_m']<=1e-7
                assert f['geometry']['coverage_error_m2']<=max(1e-7,f['area_m2']*1e-9)
                if not f['fixed_terrain_triangles']:
                    assert f['zero_steiner_triangles']==f['boundary_only_expected']
                    assert f['triangles']==f['boundary_only_expected']+2*len(f['optional_vertices'])
            for p,h in r['source_hashes'].items():assert digest(ROOT/p)==h,p
            for p,h in r['inputs'].items():assert digest(ROOT/p)==h,p
            for p,h in r['native_binaries'].items():assert digest(Path(p))==h,p
            paired=json.loads((WORK/f'shared-edge-pilot/case-{case:02d}-dtcc_mesher.json').read_text())
            assert r['native_binaries']==paired['native_binaries']
            early=json.loads((WORK/'shared-edge-pilot/native-baseline.json').read_text())
            result['earlier_pre_rebuild_native_snapshot']=early
            result['paired_shared_pilot_binaries']=paired['native_binaries']
            resource=json.loads(path.with_name(path.stem+'-resource.json').read_text())
            assert resource['status']=='completed' and resource['within_memory_budget']
            entries[f'{variant}-{backend}']=dict(preservation=r['preservation'],topology=a['topology'],
                terrain_locality=a['terrain_locality'],optional_vertices=len(a['optional_vertices']),
                optional_per_source=a['source_optional_vertices'],seconds=r['seconds'],peak_rss_bytes=r['peak_rss_bytes'])
        result['cases'][str(case)]=entries
    tests={}
    for name,expected in [('targeted-tests-final.log','118 passed'),('builder-tests-delivery.log','740 passed, 2 skipped')]:
        text=(OUT/name).read_text();assert expected in text,name
        tests[name]=dict(result=expected,sha256=digest(OUT/name))
    result.update(tests=tests,native_code_changed_by_compact_pilot=False,
                  native_binary_matches_final_shared_pilot=True,no_commit=True,
                  source_manifest={str(p.relative_to(ROOT)):digest(p) for p in [
                      *sorted((ROOT/'dtcc_core/builder/meshing').glob('compact*.py')),
                      ROOT/'dtcc_core/builder/meshing/catalogue_surface.py',
                      ROOT/'dtcc_core/builder/geometry_builders/meshes.py',
                      ROOT/'dtcc_core/model/mixins/city/builder_mixin.py',
                      ROOT/'dtcc_core/builder/__init__.py',
                      *sorted((ROOT/'sandbox').glob('minimal_catalogue_*.py')),
                      ROOT/'tests/builder/test_compact_surface.py',ROOT/'docs/minimal_catalogue_geometry.md']})
    (OUT/'validation.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(dict(valid=True,baseline_files=len(broad),baseline_mismatches=0,tests=tests)))

if __name__=='__main__':main()
