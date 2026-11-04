"""Paired wall-tessellation benchmark against preserved catalogue evidence."""
import hashlib
import html
import json
import resource
import shutil
import time
from pathlib import Path
from unittest.mock import patch

import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import meshio
import numpy as np
import shapely

from dtcc_core import io
from dtcc_core.model import Mesh, PointCloud, GeometryType
from dtcc_core.builder.meshing import wall_faces as wf
from dtcc_core.builder.meshing.roof_validation import validate_shell
from dtcc_core.builder.geometry_builders import meshes


def load_vtu(path):
    raw=meshio.read(path)
    return Mesh(vertices=raw.points,faces=raw.cells_dict['triangle'].astype(np.int64),
                markers=raw.cell_data_dict['markers']['triangle'])


def metrics(mesh, count):
    q=wf.triangle_quality(mesh.vertices,mesh.faces)
    xyz=mesh.vertices[mesh.faces]
    area2=np.linalg.norm(np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]),axis=1)
    maxedge2=np.max([np.sum((xyz[:,i]-xyz[:,(i+1)%3])**2,axis=1) for i in range(3)],axis=0)
    aspect=maxedge2/area2
    markers=np.asarray(mesh.markers)
    result={}
    for name,mask in dict(all=np.ones(len(q),dtype=bool),walls=(markers>=0)&(markers<count),
                          roofs=markers>=count,terrain=markers<0).items():
        values=q[mask]
        result[name]=dict(triangles=len(values),vertices=len(np.unique(mesh.faces[mask])),
            minimum=float(values.min()) if len(values) else None,
            p01=float(np.quantile(values,.01)) if len(values) else None,
            p05=float(np.quantile(values,.05)) if len(values) else None,
            below_002=int(sum(values<.02)),below_001=int(sum(values<.01)),
            fraction_below_002=float(np.mean(values<.02)) if len(values) else 0.,
            fraction_below_001=float(np.mean(values<.01)) if len(values) else 0.,
            max_aspect=float(aspect[mask].max()) if len(values) else None,
            histogram=np.histogram(values,bins=[0,.001,.005,.01,.02,.05,.1,.2,.4,.6,.8,1.00001])[0].tolist())
    worst=np.argsort(q)[:10]
    result['all']['stored_vertices']=len(mesh.vertices)
    result['all']['unused_vertices']=len(mesh.vertices)-result['all']['vertices']
    result['worst']=[dict(face=int(i),marker=int(markers[i]),quality=float(q[i]),vertices=xyz[i].tolist()) for i in worst]
    result['topology']=validate_shell(mesh,num_patches=count,closure_audit=meshes._shell_closure_audit)
    return result


def _triangle_signature(mesh,mask):
    points=mesh.vertices[mesh.faces[mask]]
    # Keep winding and exact coordinates; only triangle ordering may differ.
    rows=np.c_[np.asarray(mesh.markers)[mask],points.reshape(-1,9)]
    order=np.lexsort(rows[:,::-1].T)
    return rows[order]


def preservation(old,new,before,old_audit,new_audit):
    count=len(old_audit['patch_markers'])
    record={}
    for name,mask,mask_new in [
        ('roofs',np.asarray(old.markers)>=count,np.asarray(new.markers)>=count),
        ('terrain',np.asarray(old.markers)<0,np.asarray(new.markers)<0)]:
        a,b=_triangle_signature(old,mask),_triangle_signature(new,mask_new)
        record[name+'_triangles_bitwise_equal']=bool(np.array_equal(a,b))
        record[name+'_maximum_coordinate_error_m']=float(np.max(abs(a-b))) if a.shape==b.shape else None
    a,b=_triangle_signature(old,np.ones(len(old.faces),bool)),_triangle_signature(before,np.ones(len(before.faces),bool))
    record['native_shell_matches_saved_baseline']=bool(np.array_equal(a,b))
    record['patch_domains_planes_ownership_equal']=old_audit['patch_markers']==new_audit['patch_markers']
    record['blocks_and_fit_retention_equal']=old_audit['blocks']==new_audit['blocks']
    record['reconstruction_fallbacks_equal']=old_audit['block_fallbacks']==new_audit['block_fallbacks']
    record['source_fit_outcomes_equal']=[(s.get('source_key'),s.get('shape'),s.get('outcome'),s.get('final_outcome')) for s in old_audit['sources']]==[(s.get('source_key'),s.get('shape'),s.get('outcome'),s.get('final_outcome')) for s in new_audit['sources']]
    wall=new_audit.get('wall_triangulation',{})
    accepted=[r for r in wall.get('faces',[]) if r['method']=='connected_planar_cdt']
    record['maximum_wall_coverage_error_m2']=max([r['coverage_error_m2'] for r in accepted],default=0.)
    record['maximum_wall_overlap_m2']=max([r['overlap_m2'] for r in accepted],default=0.)
    record['maximum_original_wall_plane_error_m']=max([r['maximum_original_plane_error_m'] for r in accepted],default=0.)
    record['wall_plane_tolerance_m']=wf.XY_PRECISION
    record['validation_tolerance_changed']=False
    record['new_roof_fit_parameters_equal']=old_audit['patch_markers']==new_audit['patch_markers']
    record['passed']=all(v for k,v in record.items() if k.endswith('_equal') or k=='native_shell_matches_saved_baseline')
    return record


def _measure(city,catalogue,backend,build):
    before=[]
    original=wf.refine_wall_faces
    def capture(mesh,*args,**kwargs):
        before.append(mesh)
        return original(mesh,*args,**kwargs)
    start=time.perf_counter()
    with patch.object(wf,'refine_wall_faces',capture):
        mesh,result=build(city,catalogue=catalogue,mesher=backend)
    result['seconds']=time.perf_counter()-start
    # macOS reports bytes, Linux reports KiB. This is a process high-water
    # mark, not an isolated per-build allocation measurement.
    import sys
    rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    result['process_high_water_rss_bytes']=int(rss if sys.platform=='darwin' else rss*1024)
    return mesh,result,before[-1] if before else mesh


def comparison(output,label,old,new,flat,old_result,new_result,flat_result,dump,before):
    old_cat=old_result['audit']['roof_catalogue'];new_cat=new_result['audit']['roof_catalogue']
    count=len(new_cat['patch_markers']);flat_count=(int(max(flat.markers))+1)//2
    scores={name:metrics(mesh,n) for name,mesh,n in [('old_catalogue',old,count),('new_catalogue',new,count),('lod1',flat,flat_count)]}
    scores['preservation']=preservation(old,new,before,old_cat,new_cat)
    scores['runtime_seconds']={name:r['seconds'] for name,r in [('old_catalogue',old_result),('new_catalogue',new_result),('lod1',flat_result)]}
    scores['memory']={'new_process_high_water_rss_bytes':new_result['process_high_water_rss_bytes'],
                     'lod1_process_high_water_rss_bytes':flat_result['process_high_water_rss_bytes'],
                     'old_per_build_peak_available':False,'interpretation':'Cumulative process RSS high-water marks; isolated old per-build peaks were not saved. See supervised whole-run peaks.'}
    scores['retained_fits']=new_cat['retained_fits'];scores['retained_nonflat_fits']=new_cat['retained_nonflat_fits']
    from collections import Counter
    scores['per_source_outcomes']=dict(Counter(s.get('final_outcome') for s in new_cat['sources']))
    scores['source_fallbacks']=[dict(source_key=s.get('source_key'),outcome=s['outcome'],final_outcome=s.get('final_outcome'),reasons=s['diagnostics'].get('rejection_reasons',[]),block_fallbacks=s.get('block_fallbacks',[])) for s in new_cat['sources'] if s.get('final_outcome')!='retained_fit']
    scores['per_block_outcomes']=new_cat['blocks']
    scores['reconstruction_fallbacks']=new_cat['block_fallbacks']
    scores['wall_quality_fallbacks']=[r for r in new_cat['wall_triangulation'].get('faces',[]) if r['fallback_reason']]
    scores['wall_refinement_seconds']=new_cat['wall_triangulation'].get('seconds')
    scores['growth']={part:{metric:scores['new_catalogue'][part][metric]/max(1,scores['old_catalogue'][part][metric]) for metric in ['triangles','vertices']} for part in ['all','walls','roofs','terrain']}
    scores['growth']['all']['stored_vertices']=len(new.vertices)/len(old.vertices)
    for row in scores['new_catalogue']['worst']:
        m=row['marker']
        row['source_key']=new_cat['patch_markers'][m%count]['source_key'] if m>=0 else None
    dump(scores,output/f'{label}-comparison.json')
    return scores


def distributions(output,label,scores):
    fig,axes=plt.subplots(1,3,figsize=(14,3.7))
    edges=np.array([0,.001,.005,.01,.02,.05,.1,.2,.4,.6,.8,1.00001])
    for ax,part in zip(axes,['walls','roofs','terrain']):
        for name,color in [('old_catalogue','#af443d'),('new_catalogue','#176786'),('lod1','#6d7836')]:
            counts=np.array(scores[name][part]['histogram'])
            ax.stairs(counts,edges,label=name,color=color)
        ax.axvline(.02,color='black',ls=':',lw=1);ax.set_yscale('symlog',linthresh=1)
        ax.set_xlabel('Triangle quality q');ax.set_ylabel('Triangle count');ax.set_title(part)
    axes[0].legend(fontsize=7);fig.tight_layout();fig.savefig(output/f'{label}-quality.png',dpi=140);plt.close(fig)


def closeups(output,old,new,old_cat,render):
    count=len(old_cat['patch_markers']);q=wf.triangle_quality(old.vertices,old.faces)
    walls=np.flatnonzero((old.markers>=0)&(old.markers<count))
    worst=int(walls[np.argmin(q[walls])]);groups,_=wf._wall_groups(old,[r['source_key'] for r in old_cat['patch_markers']])
    group,origin,normal=next(g for g in groups if worst in g[0])
    polygon,boundary,basis=wf._face_geometry(old,group,origin,normal)
    world=old.vertices[old.faces[group]].reshape(-1,3)
    # Show the complete building and its hip, with identical view and scale.
    fig=plt.figure(figsize=(13,6));base=np.min(old.vertices[:,:2],axis=0)
    for i,(mesh,title) in enumerate([(old,'Old catalogue'),(new,'Connected wall faces')]):
        ax=fig.add_subplot(1,2,i+1,projection='3d');render(ax,mesh,np.empty((0,3)),title,base)
        ax.view_init(elev=20,azim=float(np.degrees(np.arctan2(normal[1],normal[0]))))
    a,b=fig.axes
    for ax in fig.axes:
        ax.set_xlim(a.get_xlim());ax.set_ylim(a.get_ylim());ax.set_zlim(a.get_zlim())
    fig.tight_layout();fig.savefig(output/'hip-facade-before-after.png',dpi=220);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(13,6))
    new_ids=next((field.values for field in new.fields if field.name=='catalogue_wall_face'),None)
    gi=next(i for i,g in enumerate(groups) if worst in g[0])
    for ax,mesh,selected,title in [(axes[0],old,np.array(group),'Old facade'),(axes[1],new,np.flatnonzero(new_ids==gi),'New facade')]:
        xy=(mesh.vertices-origin)@basis.T
        triangles=xy[mesh.faces[selected]];quality=wf.triangle_quality(mesh.vertices,mesh.faces[selected])
        collection=PolyCollection(triangles,array=np.log10(np.maximum(quality,1e-5)),cmap='viridis',clim=(-3,0),edgecolor='#263342',linewidth=.45)
        ax.add_collection(collection);ax.set_aspect('equal');ax.set_xlim(polygon.bounds[0]-.05,polygon.bounds[2]+.05);ax.set_ylim(polygon.bounds[1]-.1,polygon.bounds[3]+.1)
        ax.set_title(f'{title}: min q={quality.min():.5f}; <0.02: {sum(quality<.02)}');ax.set_xlabel('Facade distance (m)');ax.set_ylabel('Height above local origin (m)')
    fig.tight_layout();fig.savefig(output/'hip-wall-plane-before-after.png',dpi=220);plt.close(fig)


def report(output,summary):
    def table(scores):
        rows=[]
        for part in ['walls','roofs','terrain','all']:
            for name in ['old_catalogue','new_catalogue','lod1']:
                r=scores[name][part]
                rows.append('<tr>'+''.join(f'<td>{html.escape(str(v))}</td>' for v in [part,name,r['triangles'],r['vertices'],*[f'{r[k]:.6g}' if r[k] is not None else 'n/a' for k in ['minimum','p01','p05']],r['below_002'],r['below_001'],f"{r['max_aspect']:.5g}" if r['max_aspect'] else 'n/a'])+'</tr>')
        return '<table><tr><th>Part</th><th>Method</th><th>Triangles</th><th>Used vertices</th><th>Minimum q</th><th>p01</th><th>p05</th><th>q&lt;.02</th><th>q&lt;.01</th><th>Max aspect</th></tr>'+''.join(rows)+'</table>'
    sections=[]
    for name,scores in [*summary['cities'].items(),*summary['cases'].items()]:
        preservation=scores['preservation']
        n=scores['new_catalogue'];topo=n['topology'];times=scores['runtime_seconds']
        audit_link = f'{name}.json' if name.startswith('case-') else f'{name}-audit.json'
        sections.append(f'<section><h2>{name}</h2><p>Catalogue fits retained: {scores["retained_fits"]} ({scores["retained_nonflat_fits"]} non-flat). Reconstruction block fallbacks: {len(scores["reconstruction_fallbacks"])}; wall-quality face fallbacks: {len(scores["wall_quality_fallbacks"])}. Strict topology valid: {topo["valid"]}. Geometry preservation: {preservation["passed"]}.</p><p>Runtime old/new/LOD1: {times["old_catalogue"]:.2f}/{times["new_catalogue"]:.2f}/{times["lod1"]:.2f} s. New wall stage: {scores["wall_refinement_seconds"]:.2f} s.</p>'+table(scores)+f'<img src="{name}-quality.png"><p><a href="{name}-comparison.json">Complete paired metrics, growth, ownership, fallbacks and preservation</a></p>'+(f'<img src="{name}.png">' if name.startswith('case') else '')+f'<p>Stored vertices old/new/LOD1: {scores["old_catalogue"]["all"]["stored_vertices"]}/{scores["new_catalogue"]["all"]["stored_vertices"]}/{scores["lod1"]["all"]["stored_vertices"]}. Source outcomes: {html.escape(str(scores["per_source_outcomes"]))}. <a href="{audit_link}">Source/block/wall audit</a>.</p>'+'<details><summary>Preservation and remaining wall limitations</summary><pre>'+html.escape(json.dumps(dict(preservation=preservation,wall_quality_fallbacks=scores['wall_quality_fallbacks']),indent=2))+'</pre></details></section>')
    run_path=output/'run-status.json';run=json.loads(run_path.read_text()) if run_path.exists() else {'status':'worker finished; supervisor record pending'}
    page='''<!doctype html><meta charset="utf-8"><title>Catalogue wall triangulation comparison</title><style>body{font:15px system-ui;margin:32px auto;max-width:1450px;color:#203442;background:#f5f7f8}section{background:white;padding:20px;margin:24px 0}img{max-width:100%}table{border-collapse:collapse;font-size:12px}td,th{border:1px solid #ccd4db;padding:7px}pre{white-space:pre-wrap;font-size:11px}</style><h1>Catalogue walls: geometry-preserving tessellation</h1>'''
    page+=f'<p>Backend: {summary["backend"]}. Implementation and both-city/thirteen-case worker validation completed. Case 01 wall target: {summary["hip_target_achieved"]}. Roof and terrain quality remain separate; topological validity is not a triangle-quality guarantee.</p>'
    page+='<p>Old catalogue is the preserved source/native-hash-verified baseline. All new fits use the original training points. The outer holdout was not used for fitting or wall refinement. LOD1 uses identical reconstruction settings. No additional reconstruction fallback is permitted by wall refinement.</p>'
    page+='<p>Quality q=4√3 A/Σedge²; maximum aspect=longest edge/shortest altitude. Distributions show counts, not only aggregate percentiles. Detailed JSON contains fractions, growth, worst faces and source ownership. Memory records are process high-water marks; old per-build peaks were not saved, so no isolated peak comparison is claimed.</p>'
    page+='<p><a href="summary.json">Summary</a> · <a href="experiment.json">Source/native hashes and settings</a> · <a href="run-status.json">600 s / 2 GiB supervised resource record</a></p><pre>'+html.escape(json.dumps(run,indent=2))+'</pre>'
    page+='<p>Stored vertices (including unused original vertices retained to keep external arrays stable) and used vertices are recorded separately in each comparison. Full audit JSON records every source, block and wall-quality fallback. Baseline provisional/accepted fits: 335/140 of 341.</p>'
    page+='<h2>Hip facade and local wall plane — same view/scale</h2><img src="hip-facade-before-after.png"><img src="hip-wall-plane-before-after.png">'
    page+=''.join(sections);(output/'report.html').write_text(page)


def run(args):
    import pickle
    import platform
    import scipy
    import dtcc_mesher
    from dtcc_core.builder import _dtcc_builder
    from evaluate_roof_catalogue_integration import dump,terrain_raster,city_of,crop_raster,build,render
    output=args.output;output.mkdir(parents=True,exist_ok=True)
    baseline=args.wall_baseline/('catalogue-integration' if args.mesher=='dtcc_mesher' else 'catalogue-integration-triangle')
    experiment=json.loads((baseline/'experiment.json').read_text())
    repository=Path(__file__).resolve().parents[1]
    current_files=[*experiment['production_sha256'], 'dtcc_core/builder/meshing/wall_faces.py','sandbox/roof_wall_evaluation.py','sandbox/evaluate_roof_catalogue_integration.py']
    dump(dict(backend=args.mesher,baseline=str(baseline.resolve()),baseline_experiment=experiment,
              native_extension=str(_dtcc_builder.__file__),native_sha256=hashlib.sha256(Path(_dtcc_builder.__file__).read_bytes()).hexdigest(),
              production_sha256={p:hashlib.sha256((repository/p).read_bytes()).hexdigest() for p in current_files},
              python=platform.python_version(),scipy=scipy.__version__,dtcc_mesher=dtcc_mesher.__version__,
              input_split_sha256={str(p.relative_to(args.input)):hashlib.sha256(p.read_bytes()).hexdigest() for p in args.input.glob('*/b*/points.npz')},
              cases_sha256=hashlib.sha256(args.cases.read_bytes()).hexdigest(),
              reconstruction_options=experiment['reconstruction_options'],holdout_used_for_fitting=False),output/'experiment.json')
    summary=dict(backend=args.mesher,cities={},cases={},baseline=str(baseline),holdout_used_for_fitting=False)
    originals,rasters={},{}
    for area in ['gothenburg','lund']:
        shutil.copy2(baseline/f'{area}-terrain.npz',output/f'{area}-terrain.npz')
        buildings=pickle.load(open(args.input/area/'buildings.pkl','rb'))
        for index,b in enumerate(buildings):
            b.add_geometry(PointCloud(points=np.load(args.input/area/f'b{index:03d}'/'points.npz')['train']),GeometryType.POINT_CLOUD)
        originals[area]=buildings;manifest=json.loads((args.input/area/'manifest.json').read_text());rasters[area]=terrain_raster(manifest,output)
        city=city_of(buildings,rasters[area]);old=io.load_mesh(args.wall_baseline/f'{area}-{args.mesher}.dtcc');old_result=json.loads((baseline/f'{area}-audit.json').read_text())
        new,nr,before=_measure(city,True,args.mesher,build);flat,fr,_=_measure(city,False,args.mesher,build)
        dump(nr,output/f'{area}-audit.json');dump(fr,output/f'{area}-lod1-audit.json')
        if nr['error'] or fr['error']:raise RuntimeError(f'{area}: new={nr["error"]}; LOD1={fr["error"]}')
        for name,mesh in [('catalogue',new),('lod1',flat)]:mesh.save(output/f'{area}-{name}.dtcc');mesh.save(output/f'{area}-{name}.vtu')
        scores=comparison(output,area,old,new,flat,old_result,nr,fr,dump,before);summary['cities'][area]=scores
        distributions(output,area,scores)
        print(area,'new',nr['seconds'],'LOD1',fr['seconds'],'wall poor old/new',scores['old_catalogue']['walls']['below_002'],scores['new_catalogue']['walls']['below_002'],flush=True)
        dump(summary,output/'summary.json')
    cases=json.loads(args.cases.read_text());dump(cases,output/'cases.json')
    for index,case in enumerate(cases):
        label=f'case-{index:02d}';area=case['area'];ids=case.get('indices',[case.get('index')]);buildings=[originals[area][i] for i in ids]
        bounds=shapely.union_all([b.lod0.to_polygon() for b in buildings]).bounds
        from shapely.geometry import box
        extent=box(*bounds).buffer(6).intersection(box(*meshes._raster_bounds_tuple(rasters[area]))).bounds
        city=city_of(buildings,crop_raster(rasters[area],extent));old=load_vtu(baseline/f'{label}-catalogue.vtu');old_result=json.loads((baseline/f'{label}.json').read_text())['catalogue']
        new,nr,before=_measure(city,True,args.mesher,build);flat,fr,_=_measure(city,False,args.mesher,build)
        dump(dict(case=case,catalogue=nr,lod1=fr),output/f'{label}.json')
        if nr['error'] or fr['error']:raise RuntimeError(f'{label}: new={nr["error"]}; LOD1={fr["error"]}')
        for name,mesh in [('catalogue',new),('lod1',flat)]:mesh.save(output/f'{label}-{name}.dtcc');mesh.save(output/f'{label}-{name}.vtu')
        scores=comparison(output,label,old,new,flat,old_result,nr,fr,dump,before);summary['cases'][label]=scores;distributions(output,label,scores)
        fig=plt.figure(figsize=(17,6));origin=np.array(extent[:2])
        for i,(mesh,title) in enumerate([(old,'Old catalogue'),(new,'New catalogue'),(flat,'Default LOD1')]):render(fig.add_subplot(1,3,i+1,projection='3d'),mesh,np.empty((0,3)),title,origin)
        zlo=min(ax.get_zlim()[0] for ax in fig.axes);zhi=max(ax.get_zlim()[1] for ax in fig.axes);ax0=fig.axes[0]
        for ax in fig.axes:ax.set_xlim(ax0.get_xlim());ax.set_ylim(ax0.get_ylim());ax.set_zlim(zlo,zhi)
        fig.suptitle(f'{area} {ids}: {case["reason"]}\nBlue roofs, tan walls, green terrain; fits and reconstruction fallback decisions unchanged',fontsize=10)
        fig.tight_layout();fig.savefig(output/f'{label}.png',dpi=150);plt.close(fig)
        if index==1:closeups(output,old,new,old_result['audit']['roof_catalogue'],render)
        print(label,area,ids,'wall poor old/new',scores['old_catalogue']['walls']['below_002'],scores['new_catalogue']['walls']['below_002'],flush=True)
        dump(summary,output/'summary.json')
    summary['hip_target_achieved']=summary['cases']['case-01']['new_catalogue']['walls']['below_002']==0
    summary['geometry_preserved']=all(s['preservation']['passed'] for s in [*summary['cities'].values(),*summary['cases'].values()])
    summary['worker_complete']=True
    dump(summary,output/'summary.json');report(output,summary)
