"""Fixed-camera, read-only comparison of the two frozen edge-extension pilots."""
from collections import Counter
import csv
import html
import json
from pathlib import Path
import pickle

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import numpy as np
from dtcc_core import io
from shared_edge_pilot_report import arrays,metrics,render3d,BINS
from compact_mesher_report import table
from compact_edge_extension import OUT

BASE=OUT.parent/'compact-mesher-interface'
EXPERIMENTS=('budget-only','local-candidates','energy','joint-alternatives','joint-fallback')


def measured(mesh,n,mask):
    result=metrics(mesh,n,mask)
    q=arrays(mesh,mask)[3]
    result.update(legacy_penalty=float(np.maximum(0,.25-q).sum()),energy=float(((1-q)**2).sum()))
    return result


def measures(mesh,n):
    m=mesh.markers
    return {name:measured(mesh,n,mask) for name,mask in dict(
        building=m>=0,roofs=m>=n,walls=(m>=0)&(m<n),terrain=m<0,all=np.ones(len(m),bool)).items()}


def labels(ax,points,origin,basis=None):
    for i,p in enumerate(points,1):
        xyz=np.array(p['xyz'])-origin
        if basis is not None:
            xy=xyz@basis.T;ax.scatter(*xy,c='#ef8500',edgecolors='black',s=40,zorder=5)
            ax.annotate(f'E{i}',xy,xytext=(3,-12),textcoords='offset points',fontsize=8,va='top')
        else:
            ax.scatter(*xyz,c='#ef8500',edgecolors='black',s=26,depthshade=False)
            ax.text(*xyz,f'E{i}',fontsize=8,color='#7d2f00')


def images(case,models,records,n,face):
    names=list(models);baseline=models[names[0]]
    origin=baseline.vertices.min(axis=0);xyz=baseline.vertices-origin
    limits=np.array([xyz.min(axis=0),xyz.max(axis=0)]).T
    additions=next((r['extension']['optional_vertices'] for r in records.values() if r['extension']['optional_vertices']),[])
    center=(np.array(additions[0]['xyz'])-origin if additions else np.array(face['rings_xyz'][0]).mean(axis=0)-origin)
    for view in ('whole','corner'):
        window=limits if view=='whole' else np.array([center-4,center+4]).T
        fig=plt.figure(figsize=(5.3*len(names),5.5))
        for i,(name,m) in enumerate(models.items()):
            ax=fig.add_subplot(1,len(names),i+1,projection='3d')
            q=measured(m,n,m.markers>=0)
            render3d(ax,m,n,origin,window,f'{name}: {q["triangles"]} T, <20° {q["poor"]["20"]["count"]}')
            labels(ax,records[name]['extension']['optional_vertices'],origin)
        fig.suptitle(f'Case {case:02d}: same orthographic camera (24°, −55°), true XYZ aspect; orange = additional edge points')
        fig.tight_layout();fig.savefig(OUT/f'case-{case:02d}-{view}.png',dpi=155);plt.close(fig)
    ring=np.array(face['rings_xyz'][0]);start=ring[0]
    delta=ring[np.argmax(np.linalg.norm(ring[:,:2]-start[:2],axis=1))]-start
    basis=np.array([np.r_[delta[:2]/np.linalg.norm(delta[:2]),0.],[0,0,1.]])
    local=(ring-start)@basis.T;bounds=np.array([local.min(axis=0),local.max(axis=0)]).T
    fig,axes=plt.subplots(len(names),1,figsize=(17,3*len(names)),squeeze=False)
    for ax,(name,m) in zip(axes[:,0],models.items()):
        ids=next(f.values for f in m.fields if f.name=='catalogue_surface_face')==face['face_id']
        xy=((m.vertices-start)@basis.T)[m.faces[ids]];a=arrays(m,ids)[0].min(axis=1)
        colors=np.where((a<20)[:,None],[[.98,.72,.66,1]],[[.84,.94,.98,1]])
        ax.add_collection(PolyCollection(xy,facecolors=colors,edgecolors='#223e4e',linewidths=.5))
        # Preserve global point labels even if this wall only uses a subset.
        for i,p in enumerate(records[name]['extension']['optional_vertices'],1):
            if face['face_id'] in p['affected_faces']:
                point=(np.array(p['xyz'])-start)@basis.T
                ax.scatter(*point,c='#ef8500',edgecolors='black',s=40,zorder=5)
                ax.annotate(f'E{i}',point,xytext=(3,-12),textcoords='offset points',fontsize=8,va='top')
        ax.set_xlim(bounds[0]+[-.1,.1]);ax.set_ylim(bounds[1]+[-.1,.1]);ax.set_aspect('equal',adjustable='box')
        ax.set_title(f'{name}, fixed wall {face["face_id"]}: {sum(ids)} T; <5° {sum(a<5)}, <10° {sum(a<10)}, <20° {sum(a<20)}, min {min(a):.4f}°')
        ax.set_xlabel('Along facade (m)');ax.set_ylabel('Height offset (m)')
    fig.suptitle('True 1:1 wall-plane aspect. Red = below 20°. Same domain, support plane and owner as wall-native-control.')
    fig.tight_layout();fig.savefig(OUT/f'case-{case:02d}-wall.png',dpi=170);plt.close(fig)
    fig,axes=plt.subplots(1,len(names),figsize=(5.3*len(names),6),squeeze=False)
    for ax,(name,m) in zip(axes[0],models.items()):
        mask=m.markers>=n;xy=m.vertices[m.faces[mask],:2]-origin[:2]
        a=arrays(m,mask)[0].min(axis=1)
        ax.add_collection(PolyCollection(xy,facecolors=np.where((a<20)[:,None],[[.98,.72,.66,1]],[[.75,.87,.94,1]]),edgecolors='#223e4e',linewidths=.5))
        for i,p in enumerate(records[name]['extension']['optional_vertices'],1):
            v=np.array(p['xyz'])[:2]-origin[:2];ax.scatter(*v,c='#ef8500',edgecolors='black',s=35,zorder=4)
            ax.annotate(f'E{i}',v,xytext=(3,4),textcoords='offset points',fontsize=8)
        ax.set_xlim(limits[0]);ax.set_ylim(limits[1]);ax.set_aspect('equal');ax.set_title(f'{name}, roof <20°: {sum(a<20)}')
        ax.set_xlabel('x (m)');ax.set_ylabel('y (m)')
    fig.suptitle('Roof-only view exposes transferred defects. Same XY limits and aspect; red = below 20°.')
    fig.tight_layout();fig.savefig(OUT/f'case-{case:02d}-roofs.png',dpi=170);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(16,4))
    for ax,part in zip(axes,('walls','roofs','terrain')):
        for name,m in models.items():
            mask=dict(walls=(m.markers>=0)&(m.markers<n),roofs=m.markers>=n,terrain=m.markers<0)[part]
            ax.stairs(measured(m,n,mask)['angle_histogram'],BINS,label=name)
        ax.set_title(part);ax.set_xlabel('Minimum angle (°)');ax.set_ylabel('Absolute triangle count');ax.set_yscale('symlog',linthresh=1)
    axes[0].legend(fontsize=7);fig.tight_layout();fig.savefig(OUT/f'case-{case:02d}-angles.png',dpi=170);plt.close(fig)


def main():
    summary={};sections=[]
    for case in (1,11):
        old=json.loads((BASE/f'case-{case:02d}-shared-dtcc_mesher.json').read_text())
        cat=old['audit']['roof_catalogue'];audit=cat['compact_surface'];n=len(cat['patch_markers'])
        # Match precisely the selection used by the earlier wall-native-control.
        old_edges=[p for p in audit['optional_vertices'] if p['kind']=='segment']
        wall=next((f for p in old_edges for f in audit['faces'] if f['face_id'] in p['affected_faces'] and f['kind']=='wall'),
                  max((f for f in audit['faces'] if f['kind']=='wall'),key=lambda f:f['area_m2']))
        state=pickle.load(open(OUT/f'case-{case:02d}-checkpoint.pkl','rb'));fixed=state['faces'][wall['face_id']]
        assert fixed.kind=='wall' and fixed.source_key==wall['source_key']
        assert abs(fixed.polygon.area-wall['area_m2'])<1e-7
        assert max(abs((np.array(wall['rings_xyz'][0])-fixed.origin)@fixed.normal))<1e-7
        budgets=[0,2,4] if case==1 else [0,2,4,8]
        all_scores={};all_records={};rows=[]
        for experiment in EXPERIMENTS:
            for b in budgets:
                file=OUT/experiment/f'case-{case:02d}-extra-{b}.json'
                r=json.loads(file.read_text());m=io.load_mesh(file.with_suffix('.dtcc'));s=measures(m,n)
                key=f'{experiment} +{b}';all_scores[key]=s;all_records[key]=r
                q=s['building'];e=r['extension']
                rows.append([experiment,b,len(e['optional_vertices']),q['vertices'],q['triangles'],
                    *[q['poor'][str(k)]['count'] for k in (5,10,20)],
                    *[f'{q["poor"][str(k)]["area_m2"]:.3f}' for k in (5,10,20)],
                    f'{q["legacy_penalty"]:.5f}',f'{q["energy"]:.5f}',q['maximum_vertex_incidence'],
                    f'{r["seconds"]:.3f}',f'{r["peak_rss_bytes"]/1024**2:.1f}',e['native_requests'],e['joint_trials'],e['termination']])
        models={};records={}
        for b in budgets:
            label='Current shared' if b==0 else f'Allowance +{b}'
            p=OUT/'joint-fallback'/f'case-{case:02d}-extra-{b}'
            models[label]=io.load_mesh(p.with_suffix('.dtcc'));records[label]=all_records[f'joint-fallback +{b}']
        images(case,models,records,n,wall)
        countrows=[];per_face={};per_source={};poor={};newpoints=[]
        for label,mesh in models.items():
            scores=measures(mesh,n)
            for part,q in scores.items():
                countrows.append([label,part,q['vertices'],q['triangles'],f'{q["minimum_angle"]["min"]:.6f}',
                    f'{q["minimum_angle"]["p50"]:.3f}',f'{q["q"]["min"]:.7f}',f'{q["q"]["p50"]:.4f}',
                    *[q['poor'][str(k)]['count'] for k in (5,10,20)],
                    *[f'{q["poor"][str(k)]["area_m2"]:.3f}' for k in (5,10,20)],q['maximum_vertex_incidence']])
            face_ids=next(f.values for f in mesh.fields if f.name=='catalogue_surface_face')
            per_face[label]={str(f['face_id']):measured(mesh,n,face_ids==f['face_id']) for f in audit['faces']}
            per_source[label]={source:{part:measured(mesh,n,np.isin(face_ids,[f['face_id'] for f in audit['faces']
                if f['source_key']==source and f['kind']==part])) for part in ('roof','wall')}
                for source in audit['source_optional_vertices']}
            angle=arrays(mesh)[0].min(axis=1)
            poor[label]=[dict(triangle=int(i),face=int(face_ids[i]),angle=float(angle[i]),xyz=mesh.vertices[mesh.faces[i]].tolist()) for i in np.flatnonzero(angle<20)]
            for i,p in enumerate(records[label]['extension']['optional_vertices'],1):
                newpoints.append([label,f'E{i}',p['edge_id'],f'{p["parameter"]:.12f}',p['xyz'],p['affected_faces'],p['charged_sources'],
                    p['native_trial_rank'],p['triangle_cost'],f'{p["before"]["penalty"]:.5f} → {p["after"]["penalty"]:.5f}',
                    f'{p["before"]["poor_20_count"]} → {p["after"]["poor_20_count"]}',p['reason']])
        # Keep the zero-point and native interior controls visible, unchanged.
        baseline_rows=[]
        for variant in ('zero','interior','shared'):
            m=io.load_mesh(BASE/f'case-{case:02d}-{variant}-dtcc_mesher.dtcc');q=measures(m,n)['building']
            baseline_rows.append([variant,q['vertices'],q['triangles'],*[q['poor'][str(k)]['count'] for k in (5,10,20)]])
        reasons={key:dict(Counter(r['reason'] for r in value['extension']['rejected_transactions'])) for key,value in all_records.items()}
        result=dict(scores=all_scores,records=all_records,per_face=per_face,per_source=per_source,poor_triangles=poor,
            fixed_wall=dict(face_id=wall['face_id'],source=wall['source_key'],area_m2=wall['area_m2'],
                            origin=fixed.origin.tolist(),basis=fixed.basis.tolist(),normal=fixed.normal.tolist(),verified=True),
            rejected_reasons=reasons)
        (OUT/f'case-{case:02d}-comparison.json').write_text(json.dumps(result,indent=2))
        with (OUT/f'case-{case:02d}-additional-points.csv').open('w') as stream:
            writer=csv.writer(stream);writer.writerow(['variant','point','edge','t','xyz','faces','sources','native_rank','total_triangle_cost','energy_before_after','poor20_before_after','reason']);writer.writerows(newpoints)
        evidence=json.loads((OUT/f'case-{case:02d}-extension-edge-evidence.json').read_text())
        diagnostics=Counter(t['joint']['reason'] for t in evidence['trials'])
        sections.append(f'<section><h2>Case {case:02d}: '+('hip building 14' if case==1 else 'merged [0, 3, 23, 32, 46]')+'</h2>'+
            f'<p>Fixed comparison wall: <b>{wall["face_id"]}</b>, source {html.escape(wall["source_key"])}; area {wall["area_m2"]:.6f} m². Identity checked against the frozen checkpoint using kind, owner, plane and polygon area. This is the wall in the previous <a href="../compact-mesher-interface/case-{case:02d}-wall-native-control.png">wall-native-control</a>.</p>'+
            f'<p>Fixed geometric complexity: {audit["intrinsic"]["triangles"]} intrinsic building triangles plus {audit["terrain_attachment_triangle_overhead"]} imposed by exact terrain contact gives the {audit["zero_steiner_building_triangles"]}-triangle minimal skin. Terrain itself remains {audit["terrain_triangles"]} triangles. Optional triangles are accounted for separately below; no geometric feature is removed to improve quality.</p>'+
            table(['Frozen control','Building V','Building T','<5°','<10°','<20°'],baseline_rows)+
            '<h3>Controlled budget/search comparison</h3><p>Allowances are relative to the current shared mesh. Every new point charges each incident source; per-source cap 4, no new interiors. Search caps stay 512 native requests / 256 joint trials. The exhaustive-alternative control hits the joint cap; the final fallback checks lower ranks only after rejection. No cap was increased. Time measures compact reconstruction/extension/assembly from saved fitted patches; RSS includes imports and cached terrain. Historical full-build times are in the linked baseline report and are not directly equivalent.</p>'+
            table(['Experiment','Budget','Used','V','T','<5°','<10°','<20°','Area <5','Area <10','Area <20','Old S','Energy E','Max incidence','Seconds','MiB','Native','Joint','Stop'],rows)+
            '<h3>Final bounded fallback: all surfaces and absolute defects</h3>'+table(['Variant','Surface','V','T','Min °','Median °','Min q','Median q','<5°','<10°','<20°','Area <5','Area <10','Area <20','Max incidence'],countrows)+
            ''.join(f'<a href="case-{case:02d}-{v}.png"><img src="case-{case:02d}-{v}.png"></a>' for v in ('whole','corner','wall','roofs','angles'))+
            '<h3>Every additional point</h3>'+table(['Variant','Point','Edge','Canonical t','XYZ','Incident faces','Charged sources','Native rank','Total +T','E before → after','<20 before → after','Reason'],newpoints)+
            '<p>Original shared/interior points remain fixed, with their original reasons in the frozen report. Every new boundary vertex is a single exact parent-segment interpolation, accepted for the measured incident-face energy gain above. There are no new terrain points or interior points.</p>'+
            f'<p>Uncommitted finite native diagnostic trials: {len(evidence["trials"])}; joint-compatible: {sum(t["joint"]["accepted"] for t in evidence["trials"])}. Rejections: {html.escape(str(dict(diagnostics)))}. These raw diagnostic trials also include candidates filtered by native guards; they are evidence, not additional committed search. <a href="case-{case:02d}-extension-edge-evidence.json">Reduced paired inputs and every native parameter/joint result</a>.</p>'+
            '<details><summary>Rejected transaction reasons by control</summary><pre>'+html.escape(json.dumps(reasons,indent=2))+'</pre></details>'+
            f'<p><a href="case-{case:02d}-comparison.json">Full distributions, per-source/per-face metrics, poor-triangle locations, preserved geometry checks, runtime, source charges and search diagnostics</a> · <a href="case-{case:02d}-additional-points.csv">Point accounting CSV</a></p></section>')
        summary[str(case)]={key:dict(building=s['building'],roofs=s['roofs'],walls=s['walls'],terrain=s['terrain']) for key,s in all_scores.items()}
    introduction='''<!doctype html><meta charset="utf-8"><title>Additional compact shared-edge pilot</title><style>body{font:14px system-ui;color:#193441;background:#eef3f6;margin:24px}section{padding:22px;margin:22px 0;background:white}table{border-collapse:collapse;font-size:11px}th,td{border:1px solid #c9d3db;padding:5px;max-width:290px;overflow-wrap:anywhere}img{width:100%;margin:10px 0}p{max-width:1250px;line-height:1.6}pre{white-space:pre-wrap}</style><h1>Additional shared roof-edge points: controlled two-case pilot</h1><p>The starting meshes reproduce the saved shared outputs bit-for-bit: hip 77 building triangles, group 291. All existing points are fixed. Only additional roof-edge points are permitted: ≤4/source, ≤4 unique hip / ≤8 group, with independent building ceilings 85 / 307. Geometry, ownership, fitting/fallback choices and all terrain triangles remain unchanged. Every accepted incident face passes its original coverage/plane checks, minimum-q and absolute &lt;5°/&lt;10° count/area guards. Shell closure, orientation, vertex links and shared sequences pass. No refinement cascades, terrain resampling, approximation or changed defaults.</p><p><b>Attribution:</b> budget alone helps the group but leaves the hip unchanged. Six bounded endpoint-star candidates fix a generic remote-triangle search weakness, with no extra gain in these real budget-only controls. E=Σ(1−q)² rewards moderate fans that old S=Σmax(0,.25−q) misses. Returning up to four validated segment parameters fixes the demonstrated hip local-winner/incident-roof conflict; Core tries lower ranks only after rejection. The stateless zero/one-point operation and arguments remain unchanged. The public penalty meaning changes to E, and the unreleased C result gains a small inline provenance array/count (rebuild required). Generic search and quality stay native; Core only coordinates 3D identities and complete transactions.</p><p>A generic 12×4 rectangle with four bottom intervals has old S=0 before and after. Under equal one-point budgets the old criterion does nothing; E reduces &lt;20° triangles 2→0 for +1 triangle with no &lt;10° defects. A generic asymmetric fan plus 12 remote acute components gains E≈0.307 with only the midpoint, versus ≈1.051 with segment-local proposals. <a href="synthetic-old.json">Old generic controls</a> · <a href="synthetic-local.json">Local-candidate controls</a> · <a href="synthetic-energy.json">Energy controls</a> · <a href="local-fixture-energy.json">Asymmetric/remote fixtures</a>.</p><p>The previous initial-geometry hip evidence rejected all 91 joint trials; it cannot establish impossibility for the frozen shared state with its fixed interior points. New diagnostics isolate the latter. Unsuccessful controls remain visible below.</p>'''
    conclusion='''<section><h2>Recommendation</h2><p>Choose the smallest useful allowance: +2 shared points for each case. The hip stops at two points: 77→81 building triangles, &lt;20° count 51→48 and area 226.586→152.964 m². The group at +2 uses 295 triangles (&lt;20° count 173, area 5937.569 m²); +4 uses 299 (168, 5502.356 m²); allowance +8 uses only six points and 303 triangles (167, 5398.416 m²). Each successive pair costs four triangles while reducing &lt;20° area by about 1019, 435 and 104 m²: diminishing returns. The group +4 option also removes all &lt;10° triangles on the fixed long wall, versus seven initially. Do not fill a ceiling merely because it exists.</p><p>The energy path does not dominate every control: with four points, the old-score group control has 29 triangles below 5°, versus 31 for energy; energy has lower poor area and fewer moderate defects. Both preserve the incident-face severe guards relative to their own starting state. Severe geometric slivers remain, including exact terrain-driven junctions. Neither the finite search nor these measurements prove a globally best triangulation or a minimum-angle guarantee. Keep this explicitly selected two-case pilot; do not enable production defaults or roll out cities.</p><p><a href="api-rationale.md">API rationale and reduced hip evidence</a> · <a href="validation.json">Validation, binary/source hashes and untouched baseline manifests</a> · <a href="../../../docs/compact_mesher_interface.md">API/geometry contract and reproduction commands</a> · <a href="../compact-mesher-interface/report.html">Frozen previous report and original point reasons</a></p></section>'''
    (OUT/'report.html').write_text(introduction+''.join(sections)+conclusion)
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2))


if __name__=='__main__':main()
