"""Two fixed shared-edge pilot cases using frozen training data and terrain.

Run one worker per case/backend, then generate the paired report. Baselines
are read-only. Face backend and original coverage backend are separate.
"""
import argparse
import hashlib
import json
from pathlib import Path
import pickle
import resource
import subprocess
import sys
import time
from unittest.mock import patch

import numpy as np
import shapely
from shapely.geometry import box

from dtcc_core import builder,io
from dtcc_core.model import GeometryType,PointCloud
from dtcc_core.builder.geometry_builders import meshes
from dtcc_core.builder.meshing.surface_refinement import SurfaceRefinementOptions
from evaluate_roof_catalogue_integration import terrain_raster,crop_raster,city_of
from roof_wall_evaluation import load_vtu,_triangle_signature

ROOT=Path(__file__).resolve().parents[1]
WORK=ROOT/'.agent/work'


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def evaluate(case,backend,output,variant="pilot"):
    output.mkdir(parents=True,exist_ok=True)
    prefix=f'case-{case:02d}-{backend}'
    data=WORK/'roofer-evaluation/gothenburg'
    originals=pickle.load(open(data/'buildings.pkl','rb'))
    definition=json.loads((WORK/'catalogue-evaluation/cases.json').read_text())[case]
    indices=definition.get('indices',[definition.get('index')])
    assert indices==([14] if case==1 else [0,3,23,32,46])
    buildings=[originals[i] for i in indices]
    for i,b in zip(indices,buildings):
        b.add_geometry(PointCloud(points=np.load(data/f'b{i:03d}/points.npz')['train']),GeometryType.POINT_CLOUD)
    raster=terrain_raster(json.loads((data/'manifest.json').read_text()),WORK/'catalogue-integration')
    extent=box(*shapely.union_all([b.lod0.to_polygon() for b in buildings]).bounds).buffer(6).intersection(box(*meshes._raster_bounds_tuple(raster))).bounds
    city=city_of(buildings,crop_raster(raster,extent))
    audit={};before=[]
    from dtcc_core.builder.meshing import catalogue_surface
    import dtcc_core.builder.meshing.surface_refinement as sr
    original=sr.refine_shared_surface
    from dtcc_core.builder.meshing import wall_faces as wf
    original_wall=wf.refine_wall_faces
    def capture(mesh,*args,**kwargs):
        before.append(mesh)
        return original(mesh,*args,**kwargs)
    def baseline_capture(mesh,*args,**kwargs):
        before.append(mesh)
        return (mesh,dict(method='native_catalogue_baseline')) if variant=='original' else original_wall(mesh,*args,**kwargs)
    started=time.perf_counter();error=None;mesh=None
    options=SurfaceRefinementOptions(backend=backend,max_seconds=180,max_rounds=12)
    try:
        with (patch.object(sr,'refine_shared_surface',capture) if variant=='pilot' else patch.object(wf,'refine_wall_faces',baseline_capture)):
            mesh=builder.build_city_surface_mesh(city,lod='lod2',mesher='dtcc_mesher',
                surface_refinement=options if variant=='pilot' else None,stage_audit=audit,report_mesh_quality=False,
                max_mesh_size=8.,building_mesh_triangle_size=4.,min_building_area=15.)
    except Exception as exc:error=str(exc)
    cat=audit.get('roof_catalogue',{})
    old_cat=json.loads((WORK/f'catalogue-integration/case-{case:02d}.json').read_text())['catalogue']['audit']['roof_catalogue']
    old=load_vtu(WORK/f'catalogue-integration/case-{case:02d}-catalogue.vtu')
    comparison={k:cat.get(k)==old_cat.get(k) for k in ['patch_markers','blocks','block_fallbacks']}
    comparison['source_outcomes']=([(s.get('source_key'),s.get('shape'),s.get('outcome'),s.get('final_outcome'),s.get('fallback_area')) for s in cat.get('sources',[])]==[(s.get('source_key'),s.get('shape'),s.get('outcome'),s.get('final_outcome'),s.get('fallback_area')) for s in old_cat['sources']])
    comparison['native_shell_bitwise_equal']=bool(before and np.array_equal(_triangle_signature(old,np.ones(len(old.faces),bool)),_triangle_signature(before[0],np.ones(len(before[0].faces),bool))))
    import dtcc_mesher as dm
    from dtcc_core.builder import _dtcc_builder
    import dtcc_mesher._core as native_dm
    record=dict(case=case,indices=indices,face_backend=backend,coverage_backend='dtcc_mesher',error=error,
        variant=variant,seconds=time.perf_counter()-started,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024),
        limits=dict(seconds=600,rss_bytes=2*1024**3),audit=audit,preservation=comparison,
        native_binaries={str(m.__file__):digest(m.__file__) for m in [native_dm,_dtcc_builder]},
        mesher_version=dm.__version__,holdout_used=False,
        inputs={str(p.relative_to(ROOT)):digest(p) for p in [data/'buildings.pkl',WORK/'catalogue-evaluation/cases.json',*[data/f'b{i:03d}/points.npz' for i in indices]]})
    if mesh is not None:
        if variant!='pilot':
            reference=old if variant=='original' else io.load_mesh(WORK/f'catalogue-walls-dtcc_mesher/case-{case:02d}-catalogue.dtcc')
            comparison['final_shell_bitwise_equal_to_saved_baseline']=bool(np.array_equal(_triangle_signature(reference,np.ones(len(reference.faces),bool)),_triangle_signature(mesh,np.ones(len(mesh.faces),bool))))
        mesh.save(output/f'{prefix}.dtcc')
    (output/f'{prefix}.json').write_text(json.dumps(record,indent=2))
    print(json.dumps({k:record[k] for k in ['case','face_backend','error','seconds','peak_rss_bytes','preservation']}),flush=True)
    print('accepted',cat.get('surface_refinement',{}).get('accepted'),cat.get('surface_refinement',{}).get('error'),flush=True)


def supervised(args):
    """Enforce a process deadline; audit measured peak RSS against the budget."""
    output=args.output.resolve();output.mkdir(parents=True,exist_ok=True)
    prefix=f'case-{args.case:02d}-{args.backend}'
    command=[sys.executable,str(Path(__file__).resolve()),'--case',str(args.case),
             '--backend',args.backend,'--variant',args.variant,'--output',str(output)]
    started=time.monotonic();status='completed'
    with (output/f'{prefix}.log').open('w') as log:
        process=subprocess.Popen(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
        try:code=process.wait(timeout=600)
        except subprocess.TimeoutExpired:
            process.kill();code=process.wait();status='incomplete_timeout'
    record=dict(status=status if code==0 else ('failed' if status=='completed' else status),
        exit_code=code,seconds=time.monotonic()-started,wall_limit_seconds=600,
        rss_limit_bytes=2*1024**3,memory_limit_enforcement='post-run peak audit; native face point/step budgets')
    result=output/f'{prefix}.json'
    if code==0 and result.exists():
        measured=json.loads(result.read_text())['peak_rss_bytes']
        record.update(peak_rss_bytes=measured,within_memory_budget=measured<=2*1024**3)
    (output/f'{prefix}-resource.json').write_text(json.dumps(record,indent=2))
    print(json.dumps(record),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--case',type=int,choices=[1,11],required=True);parser.add_argument('--backend',choices=['dtcc_mesher','triangle'],required=True);parser.add_argument('--output',type=Path,default=WORK/'shared-edge-pilot');parser.add_argument('--variant',choices=['pilot','original','fan'],default='pilot');parser.add_argument('--supervised',action='store_true',help='Enforce 600 s external process timeout and record peak-memory budget audit');args=parser.parse_args()
    if args.supervised:supervised(args)
    else:evaluate(args.case,args.backend,args.output,args.variant)
