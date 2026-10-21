"""Frozen two-case compact geometry, triangulation and constraint report."""
import csv
import hashlib
import html
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
import shapely

from dtcc_core import io
from roof_wall_evaluation import load_vtu, _triangle_signature
from shared_edge_pilot_report import metrics, score, arrays, BINS

ROOT=Path(__file__).resolve().parents[1]
WORK=ROOT/'.agent/work'
OUT=WORK/'minimal-catalogue'


def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def table(headers, rows):
    return '<table><tr>'+''.join('<th>'+html.escape(str(x))+'</th>' for x in headers)+'</tr>'+''.join(
        '<tr>'+''.join('<td>'+html.escape(str(x))+'</td>' for x in row)+'</tr>' for row in rows)+'</table>'


def signature(mesh,mask):return _triangle_signature(mesh,mask)


def labels_for(mesh,n):
    m=np.asarray(mesh.markers)
    return dict(roofs=m>=n,walls=(m>=0)&(m<n),terrain=m<0,building=m>=0)


def render(ax,mesh,n,origin,limits,title,polygons=None):
    m=np.asarray(mesh.markers)
    if polygons is None:
        xyz=mesh.vertices[mesh.faces]-origin
        colors=np.tile([.35,.65,.82,1],(len(m),1))
        colors[(m>=0)&(m<n)]=[.85,.72,.56,1];colors[m<0]=[.78,.83,.72,1]
        ax.add_collection3d(Poly3DCollection(xyz,facecolors=colors,edgecolors=(.1,.16,.2,.65),linewidths=.32))
    else:
        rings=[];colors=[]
        for f in polygons:
            if f['kind']=='terrain':continue
            for ring in f['rings_xyz']:
                p=np.array(ring)-origin
                rings.append(p);colors.append('#70abc8' if f['kind']=='roof' else '#dcc3a0')
                ax.scatter(p[:,0],p[:,1],p[:,2],s=4,color='#233d48')
        ax.add_collection3d(Poly3DCollection(rings,facecolors=colors,edgecolors='#233d48',linewidths=.8))
    for k,func in enumerate([ax.set_xlim,ax.set_ylim,ax.set_zlim]):func(*limits[k])
    ax.set_box_aspect(np.diff(limits,axis=1).ravel());ax.view_init(24,-55);ax.set_proj_type('ortho')
    ax.set_title(title,fontsize=10);ax.tick_params(labelsize=6)
    ax.set_xlabel('x (m)',fontsize=8);ax.set_ylabel('y (m)',fontsize=8);ax.set_zlabel('z (m)',fontsize=8)


def figures(case,zero,budget,a,n,versions,scores):
    origin=zero.vertices.min(axis=0)
    skin=zero.vertices[np.unique(zero.faces[zero.markers>=0])]-origin
    limits=np.array([skin.min(axis=0)-.5,skin.max(axis=0)+.5]).T
    items=[('Minimal polygons',zero,a['faces']),('Zero Steiner',zero,None),('Budgeted',budget,None)]
    for view in ['whole','corner']:
        window=limits.copy()
        if view=='corner':
            angle=arrays(zero)[0].min(axis=1);wall=(zero.markers>=0)&(zero.markers<n)
            worst=np.flatnonzero(wall)[np.argmin(angle[wall])]
            center=(zero.vertices[zero.faces[worst]]-origin).mean(axis=0)
            window=np.array([center-2,center+2]).T
        fig=plt.figure(figsize=(17,6))
        for i,(label,mesh,polys) in enumerate(items):
            count=metrics(mesh,n,mesh.markers>=0)
            title=(f'{label}: {len([f for f in a["faces"] if f["kind"]!="terrain"])} faces; '
                   f'{a["zero_steiner_building_triangles"]} boundary-only T' if polys else
                   f'{label}: skin {count["vertices"]} V / {count["triangles"]} T; terrain {int(sum(mesh.markers<0))} T')
            render(fig.add_subplot(1,3,i+1,projection='3d'),mesh,n,origin,window,title,polys)
        fig.suptitle(f'Case {case:02d}: '+('same orthographic camera and XYZ scale' if view=='whole' else 'same 4 m corner window, true XYZ aspect'))
        fig.tight_layout();fig.savefig(OUT/f'case-{case:02d}-{view}.png',dpi=165);plt.close(fig)
    walls=[f for f in a['faces'] if f['kind']=='wall']
    def horizontal(f):
        v=np.concatenate([np.array(r) for r in f['rings_xyz']]);return np.linalg.norm(v[:,None,:2]-v[None,:,:2],axis=2).max()
    selected=[('wall-plane',max(walls,key=lambda f:f['area_m2'])),('thin-wall',min(walls,key=horizontal)),('thin-detail',min(walls,key=horizontal))]
    face_ids=[next(x.values for x in mesh.fields if x.name=='catalogue_surface_face') for mesh in [zero,budget]]
    for view,f in selected:
        ring=np.array(f['rings_xyz'][0]);start=ring[0]
        delta=ring[np.argmax(np.linalg.norm(ring[:,:2]-start[:2],axis=1))]-start
        tangent=np.r_[delta[:2]/np.linalg.norm(delta[:2]),0.]
        basis=np.array([tangent,[0.,0,1.]])
        boundary=(ring-start)@basis.T
        bounds=np.array([boundary.min(axis=0),boundary.max(axis=0)]).T
        fig,axes=plt.subplots(1,3,figsize=(17,7 if view=='thin-wall' else 5))
        for i,(ax,(label,mesh,polys)) in enumerate(zip(axes,items)):
            if polys:
                ax.add_patch(plt.Polygon(boundary,facecolor='#e1cbaa',edgecolor='#243d48',linewidth=1))
                ax.scatter(boundary[:,0],boundary[:,1],s=12,color='#243d48')
                title=f'{label}: {len(boundary)} boundary V'
            else:
                ids=face_ids[i-1]==f['face_id']
                local=((mesh.vertices-start)@basis.T)[mesh.faces[ids]]
                angle=arrays(mesh,ids)[0].min(axis=1)
                ax.add_collection(PolyCollection(local,facecolors=np.where((angle<20)[:,None],[[.96,.65,.6,1]],[[.86,.93,.97,1]]),edgecolors='#233d48',linewidths=.5))
                title=f'{label}: {sum(ids)} T; min {angle.min():.5f}°'
            ax.set_xlim(bounds[0] + [-.01,.01]);ax.set_ylim(bounds[1]+[-.05,.05])
            if view=='thin-detail':
                radius=max(.01,min(.5,horizontal(f)*1.5))
                center=np.array([boundary[:,0].mean(),boundary[:,1].min()])
                ax.set_xlim(center[0]-radius,center[0]+radius)
                ax.set_ylim(center[1]-radius/5,center[1]+radius*1.8)
            ax.set_aspect('equal',adjustable='box');ax.set_title(title,fontsize=10)
            ax.set_xlabel('Along facade (m)');ax.set_ylabel('Height (m)')
        fig.suptitle(f'Case {case:02d}, wall {f["face_id"]}: {horizontal(f):.6f} m wide, {np.ptp(ring[:,2]):.6f} m high; true 1:1 aspect. Red: <20°.')
        fig.tight_layout();fig.savefig(OUT/f'case-{case:02d}-{view}.png',dpi=190);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(17,6))
    sources=list(dict.fromkeys(f['source_key'] for f in a['faces'] if f['kind']=='roof'))
    for i,(ax,(label,mesh,polys)) in enumerate(zip(axes,items)):
        for f in a['faces']:
            if f['kind']!='roof':continue
            color=plt.get_cmap('tab10')(sources.index(f['source_key']))
            for r in f['rings_xyz']:
                xy=np.array(r)[:,:2]-origin[:2]
                ax.add_patch(plt.Polygon(xy,facecolor=color,edgecolor='#b12d31',alpha=.5,lw=1))
            if not polys:
                ids=next(x.values for x in mesh.fields if x.name=='catalogue_surface_face')==f['face_id']
                ax.add_collection(PolyCollection(mesh.vertices[mesh.faces[ids],:2]-origin[:2],facecolors='none',edgecolors='#263e48',linewidths=.35))
        ax.autoscale();ax.set_aspect('equal');ax.set_title(label);ax.set_xlabel('x (m)');ax.set_ylabel('y (m)')
    fig.suptitle('Stable source colours and shared interfaces (red); all roof diagonals shown')
    fig.tight_layout();fig.savefig(OUT/f'case-{case:02d}-interfaces.png',dpi=170);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(17,5))
    extent=np.array([zero.vertices.min(axis=0)[:2]-origin[:2],zero.vertices.max(axis=0)[:2]-origin[:2]]).T
    for ax,(label,mesh,polys) in zip(axes,items):
        terrain=mesh.vertices[mesh.faces[mesh.markers<0],:2]-origin[:2]
        ax.add_collection(PolyCollection(terrain,facecolors='#f0f3ec',edgecolors='#b3c0b3',linewidths=.15))
        angles=arrays(mesh)[0].min(axis=1)
        for lo,hi,color in [(10,20,'#316d9f'),(5,10,'#cf8e19'),(0,5,'#bd3438')]:
            selection=(angles>=lo)&(angles<hi)
            center=mesh.vertices[mesh.faces[selection]].mean(axis=1)[:,:2]-origin[:2]
            ax.scatter(center[:,0],center[:,1],s=9,color=color,label=f'{lo}–{hi}°: {sum(selection)}',alpha=.8)
        ax.set_xlim(extent[0]);ax.set_ylim(extent[1]);ax.set_aspect('equal');ax.set_title(label)
        ax.legend(fontsize=7);ax.set_xlabel('x (m)');ax.set_ylabel('y (m)')
    fig.suptitle('Locations of all poor triangles: plan view, exact terrain extent (wall locations project onto their outlines)')
    fig.tight_layout();fig.savefig(OUT/f'case-{case:02d}-poor-locations.png',dpi=170);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(17,4))
    for ax,part in zip(axes,['roofs','walls','terrain']):
        for label in versions:
            ax.stairs(scores[label][part]['angle_histogram'],BINS,label=label)
        ax.set_yscale('symlog',linthresh=1);ax.set_title(part);ax.set_xlabel('Minimum angle (°)');ax.set_ylabel('Absolute triangle count')
    axes[0].legend(fontsize=6);fig.tight_layout();fig.savefig(OUT/f'case-{case:02d}-angles.png',dpi=160);plt.close(fig)


def main():
    sections=[];summary={}
    baseline_hashes={}
    for name in ['roofer-evaluation','catalogue-evaluation','catalogue-integration','catalogue-integration-triangle','catalogue-walls-dtcc_mesher','catalogue-walls-triangle','shared-edge-pilot']:
        baseline_hashes.update({str(p.relative_to(ROOT)):digest(p) for p in sorted((WORK/name).rglob('*')) if p.is_file()})
    (OUT/'baseline-hashes.json').write_text(json.dumps(baseline_hashes,indent=2))
    for case in [1,11]:
        records={variant:json.loads((OUT/f'case-{case:02d}-{variant}-dtcc_mesher.json').read_text()) for variant in ['zero','budgeted']}
        assert all(r['error'] is None and all(r['preservation'].values()) for r in records.values())
        cat=records['zero']['audit']['roof_catalogue'];a=cat['compact_surface'];n=len(cat['patch_markers']);owners=[p['source_key'] for p in cat['patch_markers']]
        zero=io.load_mesh(OUT/f'case-{case:02d}-zero-dtcc_mesher.dtcc');budget=io.load_mesh(OUT/f'case-{case:02d}-budgeted-dtcc_mesher.dtcc')
        versions={
            'Original DTCC':load_vtu(WORK/f'catalogue-integration/case-{case:02d}-catalogue.vtu'),
            'Fan DTCC':io.load_mesh(WORK/f'catalogue-walls-dtcc_mesher/case-{case:02d}-catalogue.dtcc'),
            'Shared refinement DTCC':io.load_mesh(WORK/f'shared-edge-pilot/case-{case:02d}-dtcc_mesher.dtcc'),
            'Shared Triangle (rejected)':io.load_mesh(WORK/f'shared-edge-pilot/case-{case:02d}-triangle.dtcc'),
            'Original Triangle':load_vtu(WORK/f'catalogue-integration-triangle/case-{case:02d}-catalogue.vtu'),
            'Fan Triangle':io.load_mesh(WORK/f'catalogue-walls-triangle/case-{case:02d}-catalogue.dtcc'),
            'Zero Steiner':zero,'Budgeted':budget}
        scores={label:score(mesh,n) for label,mesh in versions.items()}
        for label,mesh in versions.items():scores[label]['building']=metrics(mesh,n,mesh.markers>=0)
        perface={};persource={};poor={}
        for variant,mesh in [('zero',zero),('budgeted',budget)]:
            audit=records[variant]['audit']['roof_catalogue']['compact_surface']
            ids=next(x.values for x in mesh.fields if x.name=='catalogue_surface_face')
            perface[variant]=[]
            for f in audit['faces']:
                row=dict(f,quality=metrics(mesh,n,ids==f['face_id']))
                xyz=np.concatenate([np.array(r) for r in f['rings_xyz']])
                row['horizontal_span_m']=float(np.linalg.norm(xyz[:,None,:2]-xyz[None,:,:2],axis=2).max())
                row['vertical_span_m']=float(np.ptp(xyz[:,2]))
                row['retained_boundary_reasons']=[r for r in audit['graph']['retained_vertices'] if f['face_id'] in r['incident_faces']]
                perface[variant].append(row)
            persource[variant]={}
            for source in dict.fromkeys(owners):
                selection=np.array([m>=0 and owners[int(m)%n]==source for m in mesh.markers])
                persource[variant][source]={part:metrics(mesh,n,selection&mask) for part,mask in labels_for(mesh,n).items() if part in ['roofs','walls','building']}
            angles=arrays(mesh)[0].min(axis=1)
            poor[variant]=[dict(triangle=int(i),kind=('terrain' if mesh.markers[i]<0 else 'roof' if mesh.markers[i]>=n else 'wall'),
                face_id=int(ids[i]),source_key=(None if mesh.markers[i]<0 else owners[int(mesh.markers[i])%n]),
                min_angle_degrees=float(angles[i]),xyz=mesh.vertices[mesh.faces[i]].tolist()) for i in np.flatnonzero(angles<20)]
        for label,mesh in versions.items():
            if label in ('Zero Steiner','Budgeted'):
                continue
            persource[label]={}
            for source in dict.fromkeys(owners):
                selection=np.array([m>=0 and owners[int(m)%n]==source for m in mesh.markers])
                persource[label][source]={part:metrics(mesh,n,selection&mask)
                    for part,mask in labels_for(mesh,n).items() if part in ['roofs','walls','building']}
        # Fixed terrain outside the contact star must survive as exact triangles.
        old=versions['Original DTCC'];native_terrain=signature(old,old.markers<0)
        terrain_matches={}
        for variant,mesh in [('zero',zero),('budgeted',budget)]:
            signature_new=signature(mesh,mesh.markers<0)
            # The compact audit also tracks the original contact star directly.
            terrain_matches[variant]=dict(reference_triangles=len(native_terrain),final_triangles=len(signature_new),
                exact_retained=records[variant]['audit']['roof_catalogue']['compact_surface']['unchanged_terrain_triangles'])
        with (OUT/f'case-{case:02d}-boundary-vertices.csv').open('w') as stream:
            writer=csv.writer(stream);writer.writerow(['vertex','x','y','z','reason','incident_faces','corner_faces','fixed_terrain_faces','retention_detail'])
            for r in a['graph']['retained_vertices']:
                writer.writerow([r['vertex'],*r['xyz'],r['reason'],r['incident_faces'],r['corner_faces'],r.get('fixed_terrain_faces',[]),r['collinear_retained_because']])
        figures(case,zero,budget,a,n,versions,scores)
        resources={}
        for key in ['original','fan']:
            r=json.loads((WORK/f'shared-edge-pilot/isolated-{key}/case-{case:02d}-dtcc_mesher.json').read_text())
            resources[key]=dict(seconds=r['seconds'],peak_rss_bytes=r['peak_rss_bytes'],measurement='saved isolated build')
        for backend in ['dtcc_mesher','triangle']:
            r=json.loads((WORK/f'shared-edge-pilot/case-{case:02d}-{backend}.json').read_text())
            resources['shared '+backend]=dict(seconds=r['seconds'],peak_rss_bytes=r['peak_rss_bytes'],accepted=r['audit']['roof_catalogue']['surface_refinement']['accepted'])
        for label,directory in [('original Triangle','catalogue-integration-triangle'),('fan Triangle','catalogue-walls-triangle')]:
            historic=json.loads((WORK/directory/f'case-{case:02d}.json').read_text())['catalogue']
            resources[label]=dict(seconds=historic['seconds'],peak_rss_bytes=historic.get('process_high_water_rss_bytes'),
                                  measurement='historical build time; RSS cumulative across evaluation, or unavailable')
        for variant,r in records.items():resources[variant]=dict(seconds=r['seconds'],peak_rss_bytes=r['peak_rss_bytes'],measurement='isolated worker; imports excluded from build time')
        backend_validation={}
        for backend in ['geos','triangle']:
            r=json.loads((OUT/f'case-{case:02d}-zero-{backend}.json').read_text())
            backend_validation[backend]=dict(error=r['error'],preservation=r['preservation'],seconds=r['seconds'])
            if r['error'] is None:
                mesh=io.load_mesh(OUT/f'case-{case:02d}-zero-{backend}.dtcc')
                backend_validation[backend]['scores']=score(mesh,n)
        result=dict(scores=scores,per_face=perface,per_source=persource,poor_triangles=poor,records=records,
                    resources=resources,terrain_preservation=terrain_matches,independent_backends=backend_validation)
        (OUT/f'case-{case:02d}-comparison.json').write_text(json.dumps(result,indent=2))
        summary[str(case)]=dict(scores=scores,resources=resources,intrinsic={k:v for k,v in a['intrinsic'].items() if k!='graph'},
            attachment_overhead=a['terrain_attachment_triangle_overhead'],optional=len(records['budgeted']['audit']['roof_catalogue']['compact_surface']['optional_vertices']))
        count_rows=[]
        for label,s in scores.items():
            for part in ['building','roofs','walls','terrain']:
                x=s[part];count_rows.append([label,part,x['vertices'],x['triangles'],f'{x["minimum_angle"]["min"]:.5f}',
                    *[x['poor'][str(k)]['count'] for k in [5,10,20]],f'{x["poor"]["20"]["area_m2"]:.5f}',
                    f'{x["q"]["min"]:.6g}',x['maximum_vertex_incidence']])
        budget_a=records['budgeted']['audit']['roof_catalogue']['compact_surface']
        skin_rows=[]
        for part in ['roofs','walls']:
            fs=[f for f in perface['zero'] if f['kind']==part[:-1]]
            skin_rows.append([part,len(fs),sum(f['raw_boundary_vertices'] for f in fs),sum(f['removed_redundant_vertices'] for f in fs),
                              sum(f['boundary_vertices'] for f in fs),sum(f['triangles'] for f in fs)])
        resource_rows=[[k,f'{v["seconds"]:.3f}',('unavailable' if v['peak_rss_bytes'] is None else f'{v["peak_rss_bytes"]/1024**2:.1f}'),v.get('measurement',f'accepted={v.get("accepted")}')] for k,v in resources.items()]
        thin=min([f for f in perface['zero'] if f['kind']=='wall'],key=lambda f:f['horizontal_span_m'])
        width,height=thin['horizontal_span_m'],thin['vertical_span_m']
        # With an unsplit vertical edge H and all opposite points within width W,
        # at least one endpoint angle is <= atan(2 W/H), regardless of interiors.
        bound=float(np.degrees(np.arctan2(2*width,height))) if height else 90.
        limitations=f'<p>Thinnest wall {thin["face_id"]}: {width:.6f} m × {height:.6f} m; minimum zero-Steiner angle {thin["quality"]["minimum_angle"]["min"]:.6f}°. With an unsplit full-height side and all opposite points inside that width, one endpoint angle is at most atan(2W/H) ≈ {bound:.6f}°. Interior points alone cannot eliminate this fixed-boundary restriction. This bound is conditional on the full-height side; general polygon optimality is not claimed. Per-face input corners, dimensions and remaining triangle coordinates distinguish geometry from heuristic limitations.</p>'
        sections.append(f'<section><h2>Case {case:02d}: '+('building 14, hip' if case==1 else 'merged sources [0, 3, 23, 32, 46]')+'</h2>'+f'<p>Intrinsic: {a["intrinsic"]["faces"]} faces, {a["intrinsic"]["vertices"]} shared vertices, {a["intrinsic"]["triangles"]} triangles. Exact terrain attachment adds {a["terrain_attachment_triangle_overhead"]} building triangles. Final zero-Steiner skin: {a["zero_steiner_building_triangles"]} triangles; budgeted skin: {scores["Budgeted"]["building"]["triangles"]}. Terrain: {a["terrain_triangles"]} triangles, {a["unchanged_terrain_triangles"]} retained unchanged; {a["terrain_locality"]["outside_star_triangles"]} outside the contact star. Optional interior points: {len(budget_a["optional_vertices"])}; optional boundary points: 0.</p>'+table(['Surface','Faces','Raw loop V','Removed loop V','Retained loop V','Zero T'],skin_rows)+'<p>Loop vertex counts count face incidences; intrinsic vertices count unique shared identities. The CSV accounts for every boundary vertex, including terrain junctions and fixed outside-contact triangles.</p>'+table(['Variant','Surface','V','T','Min °','<5°','<10°','<20°','Area <20° m²','Min q','Max incidence'],count_rows)+table(['Method','Build seconds','Process peak MiB','Measurement'],resource_rows)+f'<p>Independent no-refinement backends: {html.escape(str(backend_validation.keys()))}; see comparison JSON for acceptance, geometry, quality and timing. Historical Triangle coverage produces a different reference terrain tessellation and is an auxiliary comparison. The failed Triangle shared-refinement attempt retains its original shell; it is not a successful improved result. Original/fan DTCC timings are the saved isolated reruns. Historical full-evaluation Triangle peaks are cumulative and are not comparable per-build allocations.</p>'+limitations+''.join(f'<a href="case-{case:02d}-{view}.png"><img src="case-{case:02d}-{view}.png"></a>' for view in ['whole','corner','wall-plane','thin-wall','thin-detail','interfaces','poor-locations','angles'])+f'<p><a href="case-{case:02d}-comparison.json">Complete per-source/per-face metrics, poor triangle coordinates, source/binary/input hashes, topology and ownership checks</a> · <a href="case-{case:02d}-boundary-vertices.csv">Every retained boundary vertex and reason</a></p></section>')
    intro='''<!doctype html><meta charset="utf-8"><title>Minimal catalogue geometry — two-case pilot</title><style>body{font:15px system-ui;max-width:1650px;margin:30px auto;background:#f3f6f8;color:#203540}section{background:white;padding:24px;margin:24px 0}img{width:100%}table{border-collapse:collapse;font-size:12px}th,td{border:1px solid #ccd5da;padding:6px}p{line-height:1.5}</style><h1>Minimal catalogue geometry — two-case pilot</h1><p>Roofs and exposed walls are constructed from the jointly noded cleaned cells, fitted affine planes and established native block bases. Equivalent source-owned planes are merged. All incident faces jointly remove redundant collinear vertices under the existing numerical policy (1e-8 m shared identity, 1e-7 m plane residual, max(1e-7, area×1e-9) m² coverage). Canonical shared 3D endpoints are retained during local projection. No fitting thresholds, footprint features or reconstruction fallbacks are changed.</p><p>Zero Steiner: constrained polygon triangulation with refinement and acute protection disabled, followed by bounded legal diagonal flips maximizing the worse q of each pair. For each remeshed connected face, T = n + 2h − 2. Unchanged terrain faces deliberately retain their exact original triangles and interior vertices; they are exempt from the boundary-only count. No global 20° target or building-size subdivision runs.</p><p>Budgeted alternative: at most two interior vertices per face, eight per source and 32 per group. Finite centroid/off-centre trials are accepted only when the face worst q improves, aggregate serious-defect penalty Σmax(0, 0.25−q) falls by at least max(0.015, 10%), and the measured improvement justifies the two-triangle cost. No boundary proposals, corner seeding or propagation. Rejected candidates consume no budget. This is a bounded heuristic, not an optimal triangulation claim.</p><p>The original native shell is generated solely to preserve baseline bases and exact terrain, and as a geometry reference; its roof and wall tessellation does not constrain the new skin. This retains some reference-generation runtime/memory cost. Only the original contact-vertex star may be retriangulated. All triangles outside it remain unchanged. Per-face coverage, support planes, ownership area, edge pairing, orientation, buried walls and vertex links are checked. Source outcomes and native reference shells match saved DTCC baselines bitwise. No holdout observations are used.</p>'''
    conclusion='''<section><h2>Recommendation</h2><p>Keep the compact mode explicitly selectable. Prefer the zero-Steiner result for catalogue storage and visualisation; choose the small-budget result only for its measured local gains. The large reduction comes from removing redundant roof/wall tessellation before meshing. Remaining growth is mostly exact terrain contact; remaining bad angles include genuine short facades and immutable boundary constraints. The finite candidate search can miss useful interior or shared-edge improvements, so unchanged defects are not all proved unavoidable. Boundary-point optimisation is intentionally deferred, rather than spending the budget on unmeasured propagation.</p><p>Compactness trades away angle quality: in these cases, budgeted insertion reduces the severe &lt;5° defects, but it does not uniformly improve every angle band or poor-area measure. Original meshes can have better worst angles at substantially greater triangle counts. This pilot does not establish universally high-quality coarse meshes.</p><p>No geometric approximation, full-city rollout or volume support is included. Geometry fidelity takes priority over a universal angle target. The original native reference stage remains a performance limitation and could later be extracted into a terrain/base-only stage without changing its rules.</p></section>'''
    (OUT/'report.html').write_text(intro+''.join(sections)+conclusion+'<p><a href="baseline-hashes.json">Preserved baseline hashes</a> · <a href="validation.json">Validation</a> · <a href="../../../docs/minimal_catalogue_geometry.md">Reproduction and policy</a></p>')
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2))

if __name__=='__main__':main()
