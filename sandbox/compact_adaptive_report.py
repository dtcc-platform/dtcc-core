"""Gallery and quality/growth accounting for the frozen public-pipeline matrix."""
from collections import Counter
import argparse
import html
import json
import os
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import numpy as np
from dtcc_core import io
from compact_adaptive_evaluation import OUT,VARIANTS
from compact_mesher_report import table
from shared_edge_pilot_report import arrays,metrics,render3d

LABELS={'zero':'Zero points','A':'Stage A baseline B','standard':'Standard extension','doubled':'Fixed doubled','adaptive':'Adaptive continuation'}


def scores(mesh,n):
    return {name:metrics(mesh,n,mask) for name,mask in dict(building=mesh.markers>=0,
        roofs=mesh.markers>=n,walls=(mesh.markers>=0)&(mesh.markers<n),terrain=mesh.markers<0).items()}


def points(ax,records,origin,face=None,basis=None):
    for i,p in enumerate(records):
        if face is not None and face not in p['affected_faces']:continue
        xyz=np.array(p['xyz'])-origin
        color=('#169c69' if p.get('phase')=='adaptive' else '#ee8200') if p['stage']=='B' else '#9f3980'
        if basis is None:
            ax.scatter(*xyz,s=18,c=color,depthshade=False)
        else:
            xy=xyz@basis.T;ax.scatter(*xy,s=25,c=color,zorder=4)
            ax.annotate(f'{p["stage"]}{i+1}',xy,xytext=(2,-10),textcoords='offset points',fontsize=7,va='top')


def gallery(case,models,records,quality):
    valid=next(iter(models.values()),None)
    if valid is None:
        fig,ax=plt.subplots(figsize=(12,3));ax.axis('off');ax.text(.02,.8,f'Case {case}: no valid initial compact model',fontsize=16)
        ax.text(.02,.6,next(iter(records.values()))['error'],wrap=True,fontsize=8)
    else:
        origin=valid.vertices.min(axis=0);xyz=valid.vertices-origin;limits=np.array([xyz.min(axis=0),xyz.max(axis=0)]).T
        fig=plt.figure(figsize=(25,5))
        for i,v in enumerate(VARIANTS):
            ax=fig.add_subplot(1,5,i+1,projection='3d')
            if v not in models:
                ax.text2D(.02,.5,records[v]['error'],transform=ax.transAxes,fontsize=7,wrap=True);ax.set_title(LABELS[v]);continue
            q=quality[v]['building'];cat=records[v]['audit']['roof_catalogue'];n=len(cat['patch_markers'])
            render3d(ax,models[v],n,origin,limits,f'{LABELS[v]}: {q["triangles"]} T; <20° {q["poor"]["20"]["count"]}')
            points(ax,cat['compact_surface']['optional_vertices'],origin)
        fig.suptitle(f'Case {case:02d}: identical orthographic camera (24°, −55°), true XYZ aspect; magenta A, orange standard B, green adaptive B')
    fig.tight_layout();fig.savefig(OUT/f'case-{case:02d}-whole.png',dpi=140);plt.close(fig)


def details(case,models,records,quality):
    if len(models)!=5:return
    base=models['A'];cat=records['A']['audit']['roof_catalogue'];a=cat['compact_surface'];n=len(cat['patch_markers'])
    faceids=next(f.values for f in base.fields if f.name=='catalogue_surface_face')
    # Keep the same wall as the historical two-case controls when its owner,
    # area and domain identity still agree; other cases use their largest wall.
    wall=max((f for f in a['faces'] if f['kind']=='wall'),key=lambda f:f['area_m2'])
    if case in (1,11):
        old=json.loads((OUT.parent/f'compact-mesher-interface/case-{case:02d}-shared-dtcc_mesher.json').read_text())['audit']['roof_catalogue']['compact_surface']
        target=next(f for f in old['faces'] if f['face_id']==(9 if case==1 else 33))
        matches=[f for f in a['faces'] if f['kind']=='wall' and f['source_key']==target['source_key'] and abs(f['area_m2']-target['area_m2'])<1e-7]
        if matches:wall=matches[0]
    selected=faceids==wall['face_id'];p=base.vertices[np.unique(base.faces[selected])]
    start=p[0];d=p[np.argmax(np.linalg.norm(p[:,:2]-start[:2],axis=1))]-start
    basis=np.array([np.r_[d[:2]/np.linalg.norm(d[:2]),0.],[0,0,1.]])
    local=(p-start)@basis.T;bounds=np.array([local.min(axis=0),local.max(axis=0)]).T
    fig,axes=plt.subplots(5,1,figsize=(17,15))
    for ax,v in zip(axes,VARIANTS):
        mesh=models[v];ids=next(f.values for f in mesh.fields if f.name=='catalogue_surface_face')==wall['face_id']
        xy=((mesh.vertices-start)@basis.T)[mesh.faces[ids]];angles=arrays(mesh,ids)[0].min(axis=1)
        ax.add_collection(PolyCollection(xy,facecolors=np.where((angles<20)[:,None],[[.98,.72,.66,1]],[[.82,.93,.98,1]]),edgecolors='#243d4b',linewidths=.4))
        ax.set_xlim(bounds[0]+[-.1,.1]);ax.set_ylim(bounds[1]+[-.1,.1]);ax.set_aspect('equal',adjustable='box')
        ax.set_title(f'{LABELS[v]}, wall {wall["face_id"]}: {sum(ids)} T; <5° {sum(angles<5)}, <10° {sum(angles<10)}, <20° {sum(angles<20)}, min {min(angles):.5f}°')
        points(ax,records[v]['audit']['roof_catalogue']['compact_surface']['optional_vertices'],start,face=wall['face_id'],basis=basis)
        ax.set_xlabel('Along wall (m)');ax.set_ylabel('Height offset (m)')
    fig.suptitle('True-aspect fixed wall domain and owner. Red: minimum angle <20°. Magenta: A; orange: standard B; green: adaptive B.')
    fig.tight_layout();fig.savefig(OUT/f'case-{case:02d}-wall.png',dpi=160);plt.close(fig)
    origin=base.vertices.min(axis=0);xyz=base.vertices-origin;limits=np.array([xyz.min(axis=0),xyz.max(axis=0)]).T
    fig,axes=plt.subplots(1,5,figsize=(25,5))
    for ax,v in zip(axes,VARIANTS):
        mesh=models[v];ids=mesh.markers>=n;xy=mesh.vertices[mesh.faces[ids],:2]-origin[:2];angles=arrays(mesh,ids)[0].min(axis=1)
        ax.add_collection(PolyCollection(xy,facecolors=np.where((angles<20)[:,None],[[.98,.72,.66,1]],[[.72,.86,.94,1]]),edgecolors='#243d4b',linewidths=.35))
        for i,p in enumerate(records[v]['audit']['roof_catalogue']['compact_surface']['optional_vertices']):
            xy=np.array(p['xyz'])[:2]-origin[:2]
            ax.scatter(*xy,s=18,c='#ee8200' if p['stage']=='B' else '#9f3980',zorder=5)
            if p['stage']=='B':ax.annotate(f'B{i+1}',xy,xytext=(2,4),textcoords='offset points',fontsize=7)
        ax.set_xlim(limits[0]);ax.set_ylim(limits[1]);ax.set_aspect('equal');ax.set_title(f'{LABELS[v]}, roof <20° {sum(angles<20)}');ax.set_xlabel('x (m)');ax.set_ylabel('y (m)')
    fig.suptitle('Roof-only comparison exposes transferred defects; identical XY limits and aspect')
    fig.tight_layout();fig.savefig(OUT/f'case-{case:02d}-roof.png',dpi=160);plt.close(fig)
    additions=records['standard']['audit']['roof_catalogue']['compact_surface']['optional_vertices']
    focus=next((np.array(p['xyz']) for p in additions if p['stage']=='B'),start)
    center=focus-origin;window=np.array([center-3,center+3]).T
    fig=plt.figure(figsize=(25,5))
    for i,v in enumerate(VARIANTS):
        ax=fig.add_subplot(1,5,i+1,projection='3d');render3d(ax,models[v],n,origin,window,LABELS[v])
        points(ax,records[v]['audit']['roof_catalogue']['compact_surface']['optional_vertices'],origin)
    fig.suptitle('Identical 6 m corner window and physical aspect')
    fig.tight_layout();fig.savefig(OUT/f'case-{case:02d}-corner.png',dpi=160);plt.close(fig)


def adaptive_focus(case,models,records):
    a=records['adaptive']['audit']['roof_catalogue']['compact_surface']
    extra=[p for p in a['optional_vertices'] if p.get('phase')=='adaptive']
    if not extra:return
    point=max(extra,key=lambda p:p['before']['penalty']-p['after']['penalty'])
    kinds={f['face_id']:f['kind'] for f in a['faces']}
    candidates=[i for i in point['affected_faces'] if kinds[i]=='wall'] or point['affected_faces']
    face=max(candidates,key=lambda i:point['per_face_before'][str(i)]['penalty']-point['per_face_after'][str(i)]['penalty'])
    base=models['adaptive'];ids=next(f.values for f in base.fields if f.name=='catalogue_surface_face')==face
    xyz=base.vertices[np.unique(base.faces[ids])];origin=xyz.mean(axis=0)
    _,_,basis=np.linalg.svd(xyz-origin,full_matrices=False);basis=basis[:2]
    local=(xyz-origin)@basis.T;bounds=np.array([local.min(axis=0),local.max(axis=0)]).T
    margin=np.maximum(np.diff(bounds,axis=1).ravel()*.01,.1)
    fig,axes=plt.subplots(5,1,figsize=(17,15))
    for ax,v in zip(axes,VARIANTS):
        mesh=models[v];mask=next(f.values for f in mesh.fields if f.name=='catalogue_surface_face')==face
        xy=((mesh.vertices-origin)@basis.T)[mesh.faces[mask]];angles=arrays(mesh,mask)[0].min(axis=1)
        ax.add_collection(PolyCollection(xy,facecolors=np.where((angles<20)[:,None],[[.98,.72,.66,1]],[[.82,.93,.98,1]]),edgecolors='#243d4b',linewidths=.45))
        points(ax,records[v]['audit']['roof_catalogue']['compact_surface']['optional_vertices'],origin,face=face,basis=basis)
        ax.set_xlim(bounds[0]+[-margin[0],margin[0]]);ax.set_ylim(bounds[1]+[-margin[1],margin[1]])
        ax.set_aspect('equal',adjustable='box');ax.set_xlabel('Local planar coordinate (m)');ax.set_ylabel('Local planar coordinate (m)')
        e=float(np.square(1-arrays(mesh,mask)[3]).sum())
        ax.set_title(f'{LABELS[v]}: {kinds[face]} {face}; {sum(mask)} T, E={e:.6f}, <20° {sum(angles<20)}')
    fig.suptitle(f'Face affected by adaptive continuation: full transaction G={point["gain_per_triangle"]:.6f}, +{point["triangle_cost"]} total triangles, incident faces {point["affected_faces"]}\nGreen: adaptive points; orange: standard B; magenta: A. Same domain and physical aspect.')
    fig.tight_layout();fig.savefig(OUT/f'case-{case:02d}-adaptive-focus.png',dpi=160);plt.close(fig)


def main(input_dir=OUT/'public',detailed_dir=OUT/'public-detailed'):
    manifest=json.loads((OUT/'case-manifest.json').read_text());summary={};rows=[];sections=[];totals=Counter();shapes=Counter()
    detailed_cases={1,4,6,8,10,11,12,15,18,22,26}
    for case in manifest['cases']:
        index=case['case'];records={v:json.loads((input_dir/f'case-{index:02d}-{v}.json').read_text()) for v in VARIANTS}
        models={v:io.load_mesh(input_dir/f'case-{index:02d}-{v}.dtcc') for v,r in records.items() if r['valid_compact']}
        quality={v:scores(m,len(records[v]['audit']['roof_catalogue']['patch_markers'])) for v,m in models.items()}
        for v,mesh in models.items():
            for part,mask in dict(building=mesh.markers>=0,roofs=mesh.markers>=len(records[v]['audit']['roof_catalogue']['patch_markers']),
                    walls=(mesh.markers>=0)&(mesh.markers<len(records[v]['audit']['roof_catalogue']['patch_markers'])),terrain=mesh.markers<0).items():
                quality[v][part]['energy']=float(np.square(1-arrays(mesh,mask)[3]).sum())
        gallery(index,models,records,quality)
        if index in detailed_cases:details(index,models,records,quality)
        if len(models)==5:adaptive_focus(index,models,records)
        selected=[dict(source=s.get('source_key'),shape=s.get('shape'),outcome=s.get('final_outcome'),
                       reasons=s.get('diagnostics',{}).get('rejection_reasons')) for s in records['standard']['audit'].get('roof_catalogue',{}).get('sources',[])]
        if all(r['lod2_validated'] for r in records.values()):
            totals['LOD2_cases']+=1;shapes.update(s['shape'] for s in selected)
        if all(r['valid_compact'] for r in records.values()):totals['compact_cases']+=1
        countrows=[];poor={};per_face={};per_source={};pointrows=[];terminations={};budgetrows=[]
        for v,mesh in models.items():
            r=records[v];cat=r['audit']['roof_catalogue'];a=cat['compact_surface'];n=len(cat['patch_markers'])
            for part,q in quality[v].items():
                countrows.append([LABELS[v],part,q['vertices'],q['triangles'],f'{q["minimum_angle"]["min"]:.6f}',f'{q["q"]["min"]:.6g}',f'{q["q"]["p50"]:.4f}',
                    *[q['poor'][str(k)]['count'] for k in (5,10,20)],*[f'{q["poor"][str(k)]["area_m2"]:.3f}' for k in (5,10,20)],q['maximum_vertex_incidence']])
            fids=next(f.values for f in mesh.fields if f.name=='catalogue_surface_face')
            per_face[v]={str(f['face_id']):metrics(mesh,n,fids==f['face_id']) for f in a['faces']}
            per_source[v]={s:{kind:metrics(mesh,n,np.isin(fids,[f['face_id'] for f in a['faces'] if f['source_key']==s and f['kind']==kind])) for kind in ('roof','wall')} for s in a['source_optional_vertices']}
            angle=arrays(mesh)[0].min(axis=1)
            poor[v]=[dict(face=int(fids[i]),triangle=int(i),angle=float(angle[i]),xyz=mesh.vertices[mesh.faces[i]].tolist()) for i in np.flatnonzero(angle<20)]
            terminations[v]=[dict(group=g['group_id'],stage=s['stage'],termination=s['termination'],
                stopping=s.get('stopping'),effective_allowances=s.get('effective_allowances'),
                standard_result=s.get('standard_result'),escalations=s.get('escalations',[]),
                rejection_counts=s.get('rejection_counts',{}),requests=s.get('native_requests'),trials=s.get('joint_trials'))
                for g in a['optimization']['groups'] for s in a['optimization']['stages'] if s['group_id']==g['group_id']]
            for group in a['optimization']['groups']:
                for stage in ('A','B'):
                    if stage not in group:continue
                    entry=group[stage]
                    rejected=sum(s.get('rejection_counts',{}).get('building triangle ceiling',0)
                        for s in a['optimization']['stages'] if s['group_id']==group['group_id'] and s['stage']==stage)
                    budgetrows.append([LABELS[v],group['group_id'],stage,group['T0'],entry['start_triangles'],
                        entry['end_triangles'],entry['triangle_ceiling'],entry['interior_points'],entry['shared_points'],
                        entry['source_charges'],rejected,entry['termination']])
            for i,p in enumerate(a['optional_vertices']):
                pointrows.append([LABELS[v],f'{p["stage"]}{i+1} '+p.get('phase','standard'),p['kind'],p.get('edge_id','—'),p.get('parameter','—'),p['xyz'],p['affected_faces'],p['charged_sources'],p['triangle_cost'],
                    f'{p["before"]["penalty"]:.5f} → {p["after"]["penalty"]:.5f}',p['reason']])
        row=[index,case['area'],', '.join(f'{s["shape"]}: {s["outcome"]}' for s in selected),
             'LOD2' if all(r['lod2_validated'] for r in records.values()) else 'fit fallback' if len(models)==5 else 'FAILED']
        row += [quality[v]['building']['triangles'] if v in quality else 'failed' for v in VARIANTS]
        for v in ('standard','doubled','adaptive'):
            if v in quality:
                q=quality[v]['building'];a=records[v]['audit']['roof_catalogue']['compact_surface']
                row += [a['optional_counts']['B']['segment'],q['triangles']-a['TB'],f'{q["triangles"]/a["T0"]:.3f}',q['poor']['20']['count'],f'{q["poor"]["20"]["area_m2"]:.2f}']
            else:row += ['—']*5
        rows.append(row)
        result=dict(records=records,scores=quality,selected_sources=selected,per_face=per_face,per_source=per_source,
                    poor_triangles=poor,terminations=terminations)
        (OUT/f'case-{index:02d}-comparison.json').write_text(json.dumps(result,indent=2))
        summary[str(index)]=dict(scores=quality,selected_sources=selected,lod2_validated=all(r['lod2_validated'] for r in records.values()),
                               valid_compact=len(models)==5,stages={v:r['audit'].get('roof_catalogue',{}).get('compact_surface',{}).get('optimization',{}) for v,r in records.items()})
        resource_rows=[[LABELS[v],f'{r["seconds"]:.3f}',f'{r["peak_rss_bytes"]/1024**2:.1f}',r['error'] or 'valid compact; '+('all fits retained' if r['lod2_validated'] else 'reconstruction fallback visible')] for v,r in records.items()]
        constraints=records['standard']['audit'].get('roof_catalogue',{}).get('compact_surface',{}).get('faces',[])
        corner_rows=[[f['face_id'],f['kind'],f['source_key'],f['angle_constraint']['minimum_input_angle_degrees'],
                      f['angle_constraint']['corner_quality_upper_bound']]
                     for f in constraints if f['kind']!='terrain' and f['angle_constraint']['minimum_input_angle_degrees']<20]
        corner_evidence=(table(['Face','Kind','Source','Fixed corner °','Corner q upper bound'],corner_rows)
                         if corner_rows else '<p>No roof/wall polygon corner below 20° in this case. Remaining poor triangles are not proved geometrically unavoidable.</p>')
        steps=[]
        if all(v in models for v in ('A','standard','doubled','adaptive')):
            e0=quality['A']['building']['energy']
            for v in ('standard','adaptive','doubled'):
                a=records[v]['audit']['roof_catalogue']['compact_surface'];t=a['TB'];e=e0
                steps.append([LABELS[v],'Stage A baseline',t,f'{e:.8f}','—','—','—'])
                for j,p in enumerate(p for p in a['optional_vertices'] if p['stage']=='B'):
                    delta=p['before']['penalty']-p['after']['penalty'];e-=delta;t+=p['triangle_cost']
                    steps.append([LABELS[v],str(j+1)+' '+p.get('phase','standard'),t,f'{e:.8f}',f'{delta/p["triangle_cost"]:.6f}',p['triangle_cost'],p['charged_sources']])
        step_html='<details open><summary>Quality versus triangle count: exact marginal transaction steps</summary>'+table(
            ['Policy','Step','Building T','E','G','Total +T','Sources'],steps)+'</details>' if steps else ''
        sections.append(f'<section id="case-{index}"><h2>Case {index:02d}: {case["area"]} {case["indices"]}</h2><p>{html.escape(case["reason"])}</p><p>'+html.escape(str(selected))+'</p>'+
            f'<a href="case-{index:02d}-whole.png"><img src="case-{index:02d}-whole.png"></a>'+
            table(['Variant','Surface','V','T','Min °','Min q','Median q','<5°','<10°','<20°','Area <5','Area <10','Area <20','Max incidence'],countrows)+
            table(['Variant','Build seconds','Process peak MiB','Outcome'],resource_rows)+
            ''.join(f'<a href="case-{index:02d}-{view}.png"><img src="case-{index:02d}-{view}.png"></a>' for view in ('wall','roof','corner','adaptive-focus') if (OUT/f'case-{index:02d}-{view}.png').exists())+
            step_html+'<details><summary>Inserted points: exact positions, reasons, incident faces, source charges and complete triangle costs</summary>'+table(['Variant','Point','Kind','Edge','Canonical t','XYZ','Incident faces','Charged sources','Total +T','E before → after','Reason'],pointrows)+'</details>'+
            '<details><summary>Search termination and rejection reasons</summary><pre>'+html.escape(json.dumps(terminations,indent=2))+'</pre></details>'+
            '<details><summary>Per-group triangle ceilings, actual source charges and binding limits</summary>'+table(
                ['Variant','Group','Stage','T0','Start T','End T','Ceiling T','Interior points','Shared points','Source charges','Ceiling rejections','Termination'],budgetrows)+'</details>'+
            '<details><summary>Fixed corner limits versus finite-search limitations</summary>'+corner_evidence+
            '<p>A fixed acute polygon corner imposes an angle bound on triangles incident to it. Other residual fans can reflect finite candidate/rank selection, neighbour guards, numerical retessellation or engineering budgets. Rejection does not prove geometric impossibility; no fitting or geometry approximation is used to improve the statistics.</p></details>'+
            f'<p><a href="case-{index:02d}-comparison.json">Full distributions, per-face/source metrics, poor triangle locations and normal diagnostics</a> · '+
            ' · '.join(f'<a href="{os.path.relpath(detailed_dir,OUT)}/case-{index:02d}-{v}.json">{v}: detailed edge graph and candidate evidence</a>' for v in ('adaptive',))+'</p></section>')
    intro=f'''<!doctype html><meta charset="utf-8"><title>Adaptive compact Stage B</title><style>body{{font:14px system-ui;margin:24px;background:#eef3f6;color:#203c4a}}section{{padding:22px;margin:24px 0;background:white}}table{{border-collapse:collapse;font-size:11px}}td,th{{border:1px solid #c8d2da;padding:5px;max-width:280px;overflow-wrap:anywhere}}img{{width:100%}}p{{max-width:1300px;line-height:1.6}}pre{{white-space:pre-wrap}}</style><h1>Adaptive compact meshing: fixed 29-case public pipeline</h1>
    <p>Five fresh controls use the same native binary, city inputs, fitting and terrain settings. Both public APIs are exercised. No optimizer checkpoints, monkey-patching or sandbox code is used in production. All attempted cases appear below: {totals['compact_cases']} valid compact cases and {totals['LOD2_cases']} fully retained LOD2 cases. Existing fitting fallbacks remain visible and are excluded from LOD2 credit.</p>
    <p>Adaptive first completes the fixed standard 4/source,8/group pass, then assesses capped intervals in the same state. Each effective allowance increases only as needed for one accepted guarded split, to hard totals 8/source,16/group. Stage B's original min(floor(1.5T0),TB+floor(.20TB)) ceiling and cumulative 512 native requests,256 joint trials,60 seconds never reset. No extra interior points or refinement cascade runs.</p>
    <p>E=Σ(1−q)² and G=ΔE/complete incident triangle cost &gt;.015 remain unchanged. Every incident face retains minimum q and absolute below-5°/10° counts and areas. A lower E does not imply improvement in every other metric. Full distributions and absolute counts/areas are reported for roofs,walls,building and terrain, along with incidence, runtime and process peak RSS. Magenta points: Stage A; orange points: standard Stage B; green points: adaptive continuation. Detailed point records distinguish phases. Build runtimes were measured with parallel case workers and do not establish runtime superiority. Process peak RSS includes earlier controls in the same case worker.</p>
    <p><a href="public/">Initial diagnostic-development matrix (preserved)</a> · <a href="settled/">Settled five-way matrix and worker logs</a> · <a href="final/">Prior diagnostic development control</a> · <a href="detailed/">29 adaptive detailed runs</a> · <a href="repeat/">Representative repeats</a>. All attempted cases are listed in validation.json; only the settled matrix supplies policy totals.</p>
    <p>Roof planes/domains, footprints, exposed walls, source ownership, original coordinates, fitting decisions and terrain attachment locality are validated. Stopping records distinguish measured quality rejection from allowances, hard limits, independent triangle budgets, runtime/search exhaustion and numerical failure. None found refers only to the evaluated finite candidate set; no global optimum or geometric impossibility is claimed.</p>'''
    headers=['Case','City','Selected shapes/outcomes','Validation','T0','TB','Standard T','Doubled T','Adaptive T']
    for label in ('Std','Double','Adaptive'):headers += [label+' B points',label+' +T from B',label+' T/T0',label+' <20',label+' area <20']
    conclusion='''<section><h2>Scope and reproduction</h2><p>Adaptive retains the standard-result prefix. Fixed doubled starts with the larger allowance and can choose a different insertion order; see transaction step tables. Generic native algorithms and the public Mesher API are unchanged. Core adds one options boolean and a compact structured stopping record.</p><p><a href="validation.json">Validation, source/native provenance, tests and frozen baseline checks</a> · <a href="policy-summary.json">Policy comparison and trade-offs</a> · <a href="api-rationale.md">API rationale</a> · <a href="case-manifest.json">Unchanged manifest</a> · <a href="../../../docs/compact_mesher_interface.md">API and reproduction</a> · <a href="../compact-integration/report.html">Frozen integration report</a></p></section>'''
    aggregate_rows=[]
    doubled_changed=[]
    for group in ('All compact cases','Fully retained LOD2 cases'):
        subset=[c for c in summary.values() if c['valid_compact'] and
                (group=='All compact cases' or c['lod2_validated'])]
        baseline_t=sum(c['scores']['A']['building']['triangles'] for c in subset)
        for v in ('A','standard','doubled','adaptive'):
            metrics_rows=[c['scores'][v]['building'] for c in subset]
            t=sum(q['triangles'] for q in metrics_rows)
            aggregate_rows.append([group,len(subset),LABELS[v],t,t-baseline_t,
                sum(p['stage']=='B' for c in subset for p in c['stages'][v]['optional_vertices']),
                *[sum(q['poor'][str(k)]['count'] for q in metrics_rows) for k in (5,10,20)],
                *[f'{sum(q["poor"][str(k)]["area_m2"] for q in metrics_rows):.3f}' for k in (5,10,20)]])
    for i,c in summary.items():
        if c['valid_compact'] and (c['scores']['standard']['building']['triangles']!=c['scores']['doubled']['building']['triangles']):
            doubled_changed.append(int(i))
    aggregate_html=('<h2>Policy totals</h2>'+table(
        ['Population','Cases','Variant','Building T','+T versus B','B edge points','<5°','<10°','<20°','Area <5','Area <10','Area <20'],aggregate_rows))
    policy=json.loads((OUT/'policy-summary.json').read_text())
    aggregate_html+='<p>'+html.escape(policy['recommendation'])+'</p><pre>'+html.escape(json.dumps({k:v for k,v in policy.items() if k not in ('cases','recommendation')},indent=2))+'</pre>'
    limitations='''<section><h2>Limits</h2><p>The compatible Mesher changes remain uncommitted. A real release/revision and portable dependency update are required before a clean install can provide the capability. The existing dependency pin/lock are unchanged. Eleven cases include fitting fallbacks; this fixed sample does not establish full-city behavior. The native finite candidate/rank search, incident guards and engineering budgets leave poor triangles; their plateau is not proof of optimality. Quality optimization need not preserve Delaunay. No commit, publishing, geometric approximation or full-city rollout is included.</p></section>'''
    (OUT/'report.html').write_text(intro+aggregate_html+table(headers,rows)+''.join(sections)+conclusion+limitations)
    (OUT/'summary.json').write_text(json.dumps(dict(totals=dict(totals),shapes=dict(shapes),cases=summary),indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir',type=Path,default=OUT/'public')
    parser.add_argument('--detailed-dir',type=Path,default=OUT/'public-detailed')
    args=parser.parse_args();main(args.input_dir,args.detailed_dir)
