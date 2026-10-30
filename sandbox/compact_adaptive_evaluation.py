"""Fresh public-pipeline comparison on a fixed training-only case manifest."""
import argparse
import hashlib
import json
from pathlib import Path
import pickle
import resource
import subprocess
import sys
import time
import traceback
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import shapely
from shapely.geometry import box
from dtcc_core import builder
from dtcc_core.model import Bounds,City,GeometryType,PointCloud,Raster,Terrain
from dtcc_core.builder.geometry_builders import meshes
from dtcc_core.builder import CompactSurfaceOptions

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'.agent/work/compact-adaptive'
VARIANTS=('zero','A','standard','doubled','adaptive')


def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def city_input(case):
    data=ROOT/'.agent/work/roofer-evaluation'/case['area']
    originals=pickle.load(open(data/'buildings.pkl','rb'))
    buildings=[originals[i] for i in case['indices']]
    for i,b in zip(case['indices'],buildings):
        b.add_geometry(PointCloud(points=np.load(data/f'b{i:03d}/points.npz')['train']),GeometryType.POINT_CLOUD)
    manifest=json.loads((data/'manifest.json').read_text())
    raster=Raster(data=np.load(ROOT/f'.agent/work/catalogue-integration/{case["area"]}-terrain.npz')['data'])
    raster.set_bounds(Bounds(*manifest['bounds']))
    bounds=meshes._raster_bounds_tuple(raster)
    extent=box(*shapely.union_all([b.lod0.to_polygon() for b in buildings]).bounds).buffer(6).intersection(box(*bounds)).bounds
    xmin,ymin,xmax,ymax=extent
    values=np.array([[raster.get_value(np.clip(x,bounds[0]+.001,bounds[2]-.001),
                                      np.clip(y,bounds[1]+.001,bounds[3]-.001))
                      for x in np.linspace(xmin,xmax,32)] for y in np.linspace(ymax,ymin,32)])
    terrain=Terrain();cropped=Raster(data=values);cropped.set_bounds(Bounds(*extent))
    terrain.add_geometry(cropped,GeometryType.RASTER)
    city=City();city.add_terrain(terrain);city.add_buildings(buildings)
    return city


def evaluate(case,variant,output,detailed=False):
    output.mkdir(parents=True,exist_ok=True)
    options=dict(zero=CompactSurfaceOptions(budgeted=False,detailed_audit=detailed),
        A=CompactSurfaceOptions(extension_source_vertices=0,extension_group_vertices=0,detailed_audit=detailed),
        standard=CompactSurfaceOptions(detailed_audit=detailed,adaptive_extension=False),
        adaptive=CompactSurfaceOptions(detailed_audit=detailed,adaptive_extension=True),
        doubled=CompactSurfaceOptions(extension_source_vertices=8,extension_group_vertices=16,detailed_audit=detailed,adaptive_extension=False))[variant]
    city=city_input(case);audit={};started=time.perf_counter();mesh=None;error=None;trace=None
    entry='builder' if case['case']%2 else 'city'
    try:
        kwargs=dict(lod='lod2',mesher='dtcc_mesher',compact_surface=options,
            max_mesh_size=8.,building_mesh_triangle_size=4.,min_building_area=15.,
            report_mesh_quality=False,stage_audit=audit)
        mesh=(builder.build_city_surface_mesh(city,**kwargs) if entry=='builder' else city.build_surface_mesh(**kwargs))
    except Exception as exc:
        error=f'{type(exc).__name__}: {exc}';trace=traceback.format_exc()
    import dtcc_mesher._core as native
    from dtcc_core.builder import _dtcc_builder
    cat=audit.get('roof_catalogue',{});compact=cat.get('compact_surface',{})
    retained=bool(cat.get('sources')) and all(s.get('final_outcome')=='retained_fit' for s in cat['sources']) and not cat.get('block_fallbacks')
    record=dict(case=case['case'],definition=case,variant=variant,entry_point=entry,detailed_audit=detailed,
        error=error,traceback=trace,seconds=time.perf_counter()-started,
        peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024),
        audit=audit,valid_compact=bool(not error and compact.get('topology',{}).get('valid')),
        lod2_validated=bool(not error and compact.get('topology',{}).get('valid') and retained),
        native_binaries={m.__file__:digest(m.__file__) for m in (native,_dtcc_builder)},
        manifest_sha256=digest(OUT/'case-manifest.json'),holdout_used=False)
    prefix=f'case-{case["case"]:02d}-{variant}'
    if mesh is not None:mesh.save(output/f'{prefix}.dtcc')
    (output/f'{prefix}.json').write_text(json.dumps(record,indent=2))
    print(case['case'],variant,'T',compact.get('building_triangles'),'LOD2',record['lod2_validated'],
          'error',error,'seconds',round(record['seconds'],3),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',type=int)
    parser.add_argument('--all',action='store_true')
    parser.add_argument('--variant',choices=VARIANTS)
    parser.add_argument('--output',type=Path,default=OUT/'matrix')
    parser.add_argument('--detailed',action='store_true')
    parser.add_argument('--jobs',type=int,default=1)
    args=parser.parse_args();manifest=json.loads((OUT/'case-manifest.json').read_text())
    if args.all:
        def worker(case):
            output=args.output.resolve();output.mkdir(parents=True,exist_ok=True)
            command=[sys.executable,str(Path(__file__).resolve()),'--case',str(case['case']),'--output',str(output)]
            if args.detailed:command+=['--detailed']
            if args.variant:command+=['--variant',args.variant]
            tag=f'case-{case["case"]:02d}'+(f'-{args.variant}' if args.variant else '')
            with (output/f'{tag}.log').open('w') as log:
                try:
                    result=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,timeout=600)
                    status=dict(exit_code=result.returncode,termination='completed')
                except subprocess.TimeoutExpired:
                    status=dict(exit_code=None,termination='case_runtime_budget')
            (output/f'{tag}-worker.json').write_text(json.dumps(status,indent=2))
            print(case['case'],status,flush=True)
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            list(pool.map(worker,manifest['cases']))
    elif args.case is not None:
        case=next(c for c in manifest['cases'] if c['case']==args.case)
        for p,h in case['inputs'].items():assert digest(p)==h,p
        for variant in ([args.variant] if args.variant else VARIANTS):evaluate(case,variant,args.output,args.detailed)
    else:parser.error('Choose --case or --all')
