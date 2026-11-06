"""Two-case measured shared-edge pilot, with frozen compact baselines."""
import csv
import hashlib
import html
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import numpy as np

from dtcc_core import io
from shared_edge_pilot_report import arrays,metrics,score,render3d,BINS
from roof_wall_evaluation import _triangle_signature

ROOT=Path(__file__).resolve().parents[1];WORK=ROOT/'.agent/work';OUT=WORK/'compact-mesher-interface';OLD=WORK/'minimal-catalogue'


def table(headers,rows):
    return '<table><tr>'+''.join('<th>'+html.escape(str(x))+'</th>' for x in headers)+'</tr>'+''.join('<tr>'+''.join('<td>'+html.escape(str(x))+'</td>' for x in row)+'</tr>' for row in rows)+'</table>'


def images(case,versions,record,n,scores):
    a=record['audit']['roof_catalogue']['compact_surface'];points=a['optional_vertices']
    main=['Zero (native)','Frozen interior','Shared budget']
    origin=versions[main[0]].vertices.min(axis=0)
    xyz=versions[main[0]].vertices-origin
    limits=np.array([xyz.min(axis=0),xyz.max(axis=0)]).T
    for view in ['whole','corner']:
        window=limits.copy()
        if view=='corner':
            edge=next((p for p in points if p['kind']=='segment'),None)
            center=(np.array(edge['xyz'])-origin if edge else
                    (versions[main[0]].vertices[versions[main[0]].faces[np.flatnonzero((versions[main[0]].markers>=0)&(versions[main[0]].markers<n))[0]]]-origin).mean(axis=0))
            radius=5 if case==11 else 2
            window=np.array([center-radius,center+radius]).T
        fig=plt.figure(figsize=(18,6))
        for i,label in enumerate(main):
            mesh=versions[label];s=scores[label]['building']
            ax=fig.add_subplot(1,3,i+1,projection='3d')
            render3d(ax,mesh,n,origin,window,f'{label}: skin {s["vertices"]} V / {s["triangles"]} T; terrain {scores[label]["terrain"]["triangles"]} T')
            if label=='Shared budget':
                for j,p in enumerate(points):
                    q=np.array(p['xyz'])-origin
                    if view=='corner' and np.any((q<window[:,0])|(q>window[:,1])):continue
                    ax.scatter(*q,s=32 if p['kind']=='segment' else 12,c='#ed8d14' if p['kind']=='segment' else '#a52882',depthshade=False)
                    if p['kind']=='segment':ax.text(*q,str(j+1),fontsize=9,color='#7b3500')
        fig.suptitle(f'Case {case:02d}: identical orthographic cameras and XYZ aspect. Orange: shared edge points; magenta: interiors.')
        fig.tight_layout();fig.savefig(OUT/f'case-{case:02d}-{view}.png',dpi=175);plt.close(fig)
    edge_records=[p for p in points if p['kind']=='segment']
    face=next((f for p in edge_records for f in a['faces'] if f['face_id'] in p['affected_faces'] and f['kind']=='wall'),
              max((f for f in a['faces'] if f['kind']=='wall'),key=lambda f:f['area_m2']))
    rings=[np.array(r) for r in face['rings_xyz']];ring=rings[0];start=ring[0]
    delta=ring[np.argmax(np.linalg.norm(ring[:,:2]-start[:2],axis=1))]-start
    basis=np.array([np.r_[delta[:2]/np.linalg.norm(delta[:2]),0.],[0,0,1.]])
    local=(ring-start)@basis.T;bounds=np.array([local.min(axis=0),local.max(axis=0)]).T
    for view,labels in [('wall-plane',main),('wall-native-control',['Native interior','Shared budget'])]:
        fig,axes=plt.subplots(len(labels),1,figsize=(17,3*len(labels)),squeeze=False)
        for ax,label in zip(axes[:,0],labels):
            mesh=versions[label];ids=next(f.values for f in mesh.fields if f.name=='catalogue_surface_face')==face['face_id']
            xy=((mesh.vertices-start)@basis.T)[mesh.faces[ids]];angles=arrays(mesh,ids)[0].min(axis=1)
            ax.add_collection(PolyCollection(xy,facecolors=np.where((angles<10)[:,None],[[.96,.65,.60,1]],[[.86,.94,.98,1]]),edgecolors='#223e4e',linewidths=.5))
            if label=='Shared budget':
                for j,p in enumerate(points):
                    if face['face_id'] not in p['affected_faces']:continue
                    q=(np.array(p['xyz'])-start)@basis.T
                    ax.scatter(*q,s=40,color='#ed8d14' if p['kind']=='segment' else '#a52882',zorder=5)
                    ax.annotate(str(j+1),q,xytext=(3,6),textcoords='offset points',fontsize=9)
            ax.set_xlim(bounds[0]+[-.1,.1]);ax.set_ylim(bounds[1]+[-.1,.1]);ax.set_aspect('equal',adjustable='box')
            ax.set_title(f'{label}, wall {face["face_id"]}: {sum(ids)} T; <5° {sum(angles<5)}, <10° {sum(angles<10)}, min {angles.min():.5f}°',fontsize=10)
            ax.set_xlabel('Along facade (m)');ax.set_ylabel('Height offset (m)')
        fig.suptitle('True 1:1 wall-plane aspect. Red: minimum angle below 10°. Labels match insertion accounting.')
        fig.tight_layout();fig.savefig(OUT/f'case-{case:02d}-{view}.png',dpi=175);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(18,6))
    for ax,label in zip(axes,main):
        mesh=versions[label];ids=mesh.markers>=n
        ax.add_collection(PolyCollection(mesh.vertices[mesh.faces[ids],:2]-origin[:2],facecolors='#b8d7e7',edgecolors='#294652',linewidths=.4))
        for f in a['faces']:
            if f['kind']!='roof':continue
            for r in f['rings_xyz']:
                p=np.array(r)[:,:2]-origin[:2];p=np.vstack([p,p[:1]])
                ax.plot(p[:,0],p[:,1],color='#aa4143',lw=.8)
        if label=='Shared budget':
            for j,p in enumerate(points):
                if p['kind']!='segment':continue
                q=np.array(p['xyz'])[:2]-origin[:2];ax.scatter(*q,s=45,color='#ed8d14',zorder=5);ax.annotate(str(j+1),q,xytext=(3,3),textcoords='offset points')
        ax.autoscale();ax.set_aspect('equal');ax.set_title(label);ax.set_xlabel('x (m)');ax.set_ylabel('y (m)')
    fig.suptitle('Roofs, unchanged source/patch interfaces (red), and inserted shared edge points (orange)')
    fig.tight_layout();fig.savefig(OUT/f'case-{case:02d}-interfaces.png',dpi=170);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(18,4))
    for ax,part in zip(axes,['roofs','walls','terrain']):
        for label in versions:ax.stairs(scores[label][part]['angle_histogram'],BINS,label=label)
        ax.set_title(part);ax.set_yscale('symlog',linthresh=1);ax.set_xlabel('Minimum angle (°)');ax.set_ylabel('Absolute triangle count')
    axes[0].legend(fontsize=7);fig.tight_layout();fig.savefig(OUT/f'case-{case:02d}-angles.png',dpi=170);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(18,6))
    for ax,label in zip(axes,main):
        mesh=versions[label];angles=arrays(mesh)[0].min(axis=1);t=mesh.vertices[mesh.faces]
        terrain=mesh.markers<0
        ax.add_collection(PolyCollection(t[terrain,:,:2]-origin[:2],facecolors='#f0f3eb',edgecolors='#b1bfae',linewidths=.15))
        for lo,hi,color in [(10,20,'#266695'),(5,10,'#cb931f'),(0,5,'#bd3339')]:
            selected=(angles>=lo)&(angles<hi);c=t[selected].mean(axis=1)[:,:2]-origin[:2]
            ax.scatter(c[:,0],c[:,1],s=9,color=color,label=f'{lo}–{hi}°: {sum(selected)}')
        ax.set_xlim(limits[0]);ax.set_ylim(limits[1]);ax.set_aspect('equal');ax.set_title(label);ax.legend(fontsize=7)
    fig.suptitle('All remaining poor-triangle locations, including terrain; wall centroids project onto their outlines')
    fig.tight_layout();fig.savefig(OUT/f'case-{case:02d}-poor-locations.png',dpi=170);plt.close(fig)


def main():
    sections=[];summary={}
    for case in (1,11):
        records={v:json.loads((OUT/f'case-{case:02d}-{v}-dtcc_mesher.json').read_text()) for v in ['zero','interior','shared']}
        assert all(r['error'] is None and all(r['preservation'].values()) for r in records.values())
        cat=records['shared']['audit']['roof_catalogue'];n=len(cat['patch_markers']);owners=[p['source_key'] for p in cat['patch_markers']]
        a=cat['compact_surface'];opt=a['optimization']
        versions={'Frozen zero':io.load_mesh(OLD/f'case-{case:02d}-zero-dtcc_mesher.dtcc'),
                  'Zero (native)':io.load_mesh(OUT/f'case-{case:02d}-zero-dtcc_mesher.dtcc'),
                  'Frozen interior':io.load_mesh(OLD/f'case-{case:02d}-budgeted-dtcc_mesher.dtcc'),
                  'Native interior':io.load_mesh(OUT/f'case-{case:02d}-interior-dtcc_mesher.dtcc'),
                  'Shared budget':io.load_mesh(OUT/f'case-{case:02d}-shared-dtcc_mesher.dtcc')}
        scores={name:score(m,n) for name,m in versions.items()}
        for name,m in versions.items():scores[name]['building']=metrics(m,n,m.markers>=0)
        per_face={};per_source={};poor={}
        for label,m in versions.items():
            ids=next(f.values for f in m.fields if f.name=='catalogue_surface_face')
            per_face[label]=[dict(face_id=f['face_id'],kind=f['kind'],source_key=f['source_key'],quality=metrics(m,n,ids==f['face_id'])) for f in a['faces']]
            per_source[label]={}
            for source in sorted(set(owners)):
                mask=np.array([x>=0 and owners[int(x)%n]==source for x in m.markers])
                per_source[label][source]={part:metrics(m,n,mask&selection) for part,selection in dict(roofs=m.markers>=n,walls=(m.markers>=0)&(m.markers<n)).items()}
            angle=arrays(m)[0].min(axis=1)
            poor[label]=[dict(triangle=int(i),face_id=int(ids[i]),angle=float(angle[i]),xyz=m.vertices[m.faces[i]].tolist()) for i in np.flatnonzero(angle<20)]
        terrain_equal=all(np.array_equal(_triangle_signature(versions['Zero (native)'],versions['Zero (native)'].markers<0),
                                        _triangle_signature(versions[label],versions[label].markers<0)) for label in ['Native interior','Shared budget'])
        assert terrain_equal
        result=dict(records=records,scores=scores,per_face=per_face,per_source=per_source,poor_triangles=poor,
                    terrain_triangles_bitwise_equal_across_new_variants=terrain_equal)
        (OUT/f'case-{case:02d}-comparison.json').write_text(json.dumps(result,indent=2))
        images(case,versions,records['shared'],n,scores)
        rows=[]
        for label,s in scores.items():
            for part in ['roofs','walls','terrain','all']:
                x=s[part];rows.append([label,part,x['vertices'],x['triangles'],f'{x["minimum_angle"]["min"]:.6f}',
                    *[x['poor'][str(k)]['count'] for k in [5,10,20]],f'{x["poor"]["5"]["area_m2"]:.6f}',
                    f'{x["poor"]["10"]["area_m2"]:.6f}',f'{x["poor"]["20"]["area_m2"]:.6f}',f'{x["q"]["min"]:.6g}',x['maximum_vertex_incidence']])
        insertion_rows=[]
        with (OUT/f'case-{case:02d}-insertions.csv').open('w') as stream:
            writer=csv.writer(stream);writer.writerow(['point','kind','x','y','z','edge','parameter','affected_faces','charged_sources','unique_points','triangle_cost','penalty_before','penalty_after','poor5_before','poor5_after','poor5_area_before','poor5_area_after','reason'])
            for j,p in enumerate(a['optional_vertices']):
                writer.writerow([j+1,p['kind'],*p['xyz'],p.get('edge_id'),p.get('parameter'),p['affected_faces'],p['charged_sources'],1,p['triangle_cost'],p['before']['penalty'],p['after']['penalty'],p['before']['poor_5_count'],p['after']['poor_5_count'],p['before']['poor_5_area'],p['after']['poor_5_area'],p['reason']])
                insertion_rows.append([j+1,p['kind'],p.get('edge_id','—'),p['affected_faces'],p['charged_sources'],p['triangle_cost'],
                    f'{p["before"]["penalty"]:.4f} → {p["after"]["penalty"]:.4f}',
                    f'{p["before"]["poor_5_count"]} → {p["after"]["poor_5_count"]}',p['reason']])
        resource_rows=[]
        for label,variant in [('Zero (native)','zero'),('Native interior','interior'),('Shared budget','shared')]:
            r=records[variant];resource_rows.append([label,f'{r["seconds"]:.3f}',f'{r["peak_rss_bytes"]/1024**2:.1f}'])
        old=json.loads((OLD/f'case-{case:02d}-budgeted-dtcc_mesher.json').read_text())
        resource_rows.append(['Frozen interior',f'{old["seconds"]:.3f}',f'{old["peak_rss_bytes"]/1024**2:.1f}'])
        edges=[p for p in a['optional_vertices'] if p['kind']=='segment']
        note=(f'<p>Native bounded calls: {opt["native_requests"]}; joint trials: {opt["joint_trials"]}; termination: <b>{opt["termination"]}</b>. Unique optional points: {len(a["optional_vertices"])} ({len(edges)} shared edge, {len(a["optional_vertices"])-len(edges)} interior); total added triangles: {sum(p["triangle_cost"] for p in a["optional_vertices"])}. Charges per source: {html.escape(str(a["source_optional_vertices"]))}. Original limits remain 8/source, 32/group and 2 interiors/face. Separate request/trial caps are 512/256 and are not raised. Cached transactions reuse measured geometry while incident versions agree and rebase only the appended identity after unrelated commits.</p>')
        if not edges:
            evidence=json.loads((OUT/'case-01-edge-evidence.json').read_text());reasons={}
            for r in evidence['trials']:reasons[r['joint']['reason']]=reasons.get(r['joint']['reason'],0)+1
            note+=f'<p>No hip edge split survived the guards. A diagnostic native-only generator supplied {len(evidence["trials"])} finite candidates to the ordinary joint transaction check; {sum(t["joint"]["accepted"] for t in evidence["trials"])} passed. Reasons: {html.escape(str(reasons))}. <a href="case-01-edge-evidence.json">Reduced paired inputs, native candidate parameters and complete joint cost/quality evidence</a>. These diagnostic transactions are never committed. This is a bounded-search result, not proof that every possible edge point is useless.</p>'
        sections.append(f'<section><h2>Case {case:02d}: '+('hip building 14' if case==1 else 'merged [0, 3, 23, 32, 46]')+'</h2>'+note+
            '<p>Fits, cleaned coverage, patch markers, ownership, block fallback decisions and the generated reference shell match the frozen reconstruction inputs. Strict shell topology and canonical edge sequences pass. Roof/wall support planes and per-face coverage use the unchanged 1e-7 m / max(1e-7, area×1e-9) m² policy. All new-variant terrain triangles are bitwise identical; outside-contact terrain is also checked against the original reference.</p>'+
            table(['Variant','Surface','V','T','Min °','<5°','<10°','<20°','Area <5° m²','Area <10° m²','Area <20° m²','Min q','Max incidence'],rows)+
            table(['Run','Build seconds','Process peak MiB'],resource_rows)+'<p>Time excludes imports/cache loading; RSS includes them. Historical frozen runs are independent saved worker measurements, not allocation deltas.</p>'+
            table(['Point','Kind','Edge','Affected faces','Charged sources','Total +T','S before → after','<5° before → after','Reason'],insertion_rows)+
            ''.join(f'<a href="case-{case:02d}-{view}.png"><img src="case-{case:02d}-{view}.png"></a>' for view in ['whole','corner','wall-plane','wall-native-control','interfaces','angles','poor-locations'])+
            f'<p><a href="case-{case:02d}-comparison.json">Full q/angle distributions, per-source/per-face metrics, all poor-triangle coordinates, source/binary/input provenance and rejected transactions</a> · <a href="case-{case:02d}-insertions.csv">Insertion accounting CSV</a></p></section>')
        summary[str(case)]=dict(scores=scores,unique_optional_points=len(a['optional_vertices']),edge_points=len(edges),termination=opt['termination'],resources=resource_rows)
    intro='''<!doctype html><meta charset="utf-8"><title>Bounded shared roof-edge pilot</title><style>body{font:15px system-ui;max-width:1750px;margin:30px auto;background:#f3f6f8;color:#20343f}section{background:white;padding:24px;margin:24px 0}table{border-collapse:collapse;font-size:12px}th,td{border:1px solid #ccd5db;padding:6px}img{width:100%}p{line-height:1.55}</style><h1>Bounded native improvement and shared roof-edge pilot</h1><p>Generic flips, candidate generation, point insertion and triangulation improvement now run in dtcc-mesher's explicit stateless improve operation. Core owns the minimal 3D geometry, canonical shared edges, source budgets, combined roof/wall acceptance and atomic commit/rollback. No public mutable topology, callback/objective framework, scheduler or 3D mesher layer is added.</p><p>The frozen interior-only baseline remains unchanged. The native interior-only control separates the effect of the new objective/implementation from shared boundary points. Native zero-point optimization may produce valid non-Delaunay triangles; normal CDT generation/defaults are unchanged. S = Σmax(0, .25−q), q = 4√3 area/Σedge². Point trials must lower S by >.015 per total added triangle, with incident-face q and absolute &lt;5°/&lt;10° counts/areas guarded. An unchanged remote global minimum is allowed. No global angle refinement, corner shells, nominal-size subdivision or terrain insertion runs.</p><p>Shared points lie on one canonical 3D segment, use returned parent/t provenance, update every incident face via native zero-extra-point meshing/optimization, and charge all affected sources. Complete candidates commit atomically or are rejected. Terrain geometry/locality, fitting and reconstruction choices remain unchanged. Only the two requested cases are evaluated.</p>'''
    end='''<section><h2>Assessment and recommendation</h2><p>The architecture now separates generic numerical improvement from city geometry. The merged group demonstrates useful shared eave points with full roof/wall cost accounting; the hip case does not justify edge insertion under the unchanged conservative guards and finite candidate set. Interior points remain useful there. The measured alternatives trade quality against triangle count; no global optimality or universal minimum angle is claimed. Thin geometry, exact terrain contact and finite-search restrictions remain visible in per-face evidence.</p><p>Use the explicitly selected bounded pilot when its measured quality gains justify the extra elements. Preserve zero and frozen interior baselines for comparison. Do not enable this by default for production catalogue meshing or roll out cities. No geometry approximation, changed fit/fallback, terrain resampling, publishing or commits are included.</p><p>Validation: native C/CLI and Python regressions, Core targeted/shared-interface and builder suites pass. Optional ASan/UBSan execution is unavailable on this host: instrumented binaries time out before main, including an independent hello-world smoke test. Repeated native allocation/error/free tests pass normally; do not treat sanitizer timeouts as successful checks.</p></section>'''
    (OUT/'report.html').write_text(intro+''.join(sections)+end+'<p><a href="api-rationale.md">API rationale written before implementation</a> · <a href="validation.json">Validation and preserved baseline hashes</a> · <a href="../../../docs/compact_mesher_interface.md">Reproduction and API accounting</a></p>')
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2))

if __name__=='__main__':main()
