"""Read-only paired geometry and spatial-quality report for the two pilot cases."""
import hashlib,html,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
import shapely
from dtcc_core import io
from dtcc_core.builder.meshing.surface_faces import extract_surface_faces
from dtcc_core.builder.meshing.wall_faces import triangle_quality
from dtcc_core.builder.meshing.roof_validation import validate_shell
from dtcc_core.builder.geometry_builders import meshes
from roof_wall_evaluation import load_vtu

ROOT=Path(__file__).resolve().parents[1];WORK=ROOT/'.agent/work';OUT=WORK/'shared-edge-pilot'
BINS=[0,1,5,10,15,20,25,30,40,60.0001]


def arrays(mesh,mask=None):
    faces=mesh.faces if mask is None else mesh.faces[mask]
    t=mesh.vertices[faces];u=np.roll(t,-1,axis=1)-t;v=np.roll(t,1,axis=1)-t
    a=np.degrees(np.arctan2(np.linalg.norm(np.cross(u,v),axis=2),np.sum(u*v,axis=2)))
    lengths=np.linalg.norm(u,axis=2);area=np.linalg.norm(np.cross(u[:,0],-v[:,0]),axis=1)/2
    return a,lengths,area,triangle_quality(mesh.vertices,faces)


def distribution(values):
    return dict(zip(['min','p01','p05','p50','p95','p99','max'],map(float,np.quantile(values,[0,.01,.05,.5,.95,.99,1])))) if len(values) else None


def metrics(mesh,n,mask=None,spacing=4):
    if mask is None:mask=np.ones(len(mesh.faces),dtype=bool)
    a,l,area,q=arrays(mesh,mask);f=mesh.faces[mask]
    if not len(f):return dict(triangles=0)
    groups={}
    for fi,t in enumerate(f):
        for e in [t[[0,1]],t[[1,2]],t[[2,0]]]:groups.setdefault(tuple(sorted(e)),[]).append(fi)
    pairs=np.array([v for v in groups.values() if len(v)==2],dtype=int).reshape(-1,2)
    exempt=next((field.values.astype(bool)[mask] for field in mesh.fields if field.name=='refinement_quality_exempt'),np.zeros(len(f),bool))
    incident_exempt=exempt[pairs].any(axis=1) if len(pairs) else np.array([],bool)
    size=np.sqrt(area);ratios=np.max(size[pairs],axis=1)/np.min(size[pairs],axis=1) if len(pairs) else np.array([])
    return dict(triangles=len(f),vertices=len(np.unique(f)),area_m2=float(area.sum()),minimum_angle=distribution(a.min(axis=1)),maximum_angle=distribution(a.max(axis=1)),q=distribution(q),
        angle_histogram=np.histogram(a.min(axis=1),BINS)[0].tolist(),
        poor={str(k):dict(count=int(sum(a.min(axis=1)<k-1e-6)),fraction=float(np.mean(a.min(axis=1)<k-1e-6)),area_m2=float(area[a.min(axis=1)<k-1e-6].sum())) for k in [5,10,20]},
        q_below_002=dict(count=int(sum(q<.02)),area_m2=float(area[q<.02].sum())),edge_length_m=distribution(l.ravel()),
        adjacent_size_ratio=distribution(ratios),adjacent_ratio_above_2=int(sum(ratios>2)),adjacent_ratio_above_4=int(sum(ratios>4)),
        adjacent_ratio_above_4_incident_backend_exempt=int(sum((ratios>4)&incident_exempt)),
        adjacent_ratio_above_4_nonexempt=int(sum((ratios>4)&~incident_exempt)),
        vertex_incidence=distribution(np.bincount(f.ravel())[np.unique(f)]),
        triangle_area_m2=distribution(area),area_above_nominal_cap=int(sum(area>.5*spacing**2*(1+1e-6))),
        maximum_vertex_incidence=int(np.bincount(f.ravel()).max()),longest_edge_above_nominal_sqrt2=int(sum(l.max(axis=1)>spacing*np.sqrt(2)*(1+1e-6))))


def score(mesh,n):
    m=np.asarray(mesh.markers);exempt=next((f.values.astype(bool) for f in mesh.fields if f.name=='refinement_quality_exempt'),np.zeros(len(m),bool))
    result={name:metrics(mesh,n,mask,8 if name=='terrain' else 4) for name,mask in dict(all=np.ones(len(m),bool),walls=(m>=0)&(m<n),roofs=m>=n,terrain=m<0).items()}
    for name,mask in dict(walls=(m>=0)&(m<n),roofs=m>=n,terrain=m<0).items():
        result[name]['backend_exempt']=metrics(mesh,n,mask&exempt)
        result[name]['not_backend_exempt']=metrics(mesh,n,mask&~exempt)
    result['all'].update(stored_vertices=len(mesh.vertices),unused_vertices=len(mesh.vertices)-result['all']['vertices'])
    result['topology']=validate_shell(mesh,num_patches=n,closure_audit=meshes._shell_closure_audit)
    return result


def render3d(ax,mesh,n,origin,limits,title,*,zoom=None):
    xyz=mesh.vertices[mesh.faces]-origin;m=np.asarray(mesh.markers)
    colors=np.tile([.38,.67,.83,1.],(len(m),1));colors[(m>=0)&(m<n)]=[.85,.73,.57,1];colors[m<0]=[.75,.81,.70,1]
    if zoom is not None:
        center,radius=zoom;mask=np.linalg.norm(xyz.mean(axis=1)[:,:2]-center[:2],axis=1)<radius*2
        xyz,colors=xyz[mask],colors[mask]
    ax.add_collection3d(Poly3DCollection(xyz,facecolors=colors,edgecolors=(.10,.14,.18,.65),linewidths=.24))
    for k,func in enumerate([ax.set_xlim,ax.set_ylim,ax.set_zlim]):func(*limits[k])
    from matplotlib.ticker import MaxNLocator
    ax.zaxis.set_major_locator(MaxNLocator(3))
    ax.set_box_aspect(np.diff(limits,axis=1).ravel());ax.view_init(24,-55);ax.set_proj_type('ortho');ax.set_title(title,fontsize=9)
    ax.set_xlabel('x (m)',fontsize=7);ax.set_ylabel('y (m)',fontsize=7);ax.set_zlabel('z (m)',fontsize=7);ax.tick_params(labelsize=6)


def plane_mask(mesh,n,face,owners):
    m=np.asarray(mesh.markers);ids=np.flatnonzero((m>=0)&(m<n));tri=mesh.vertices[mesh.faces[ids]]
    onplane=np.max(abs((tri-face.origin)@face.normal),axis=1)<1e-7
    owned=np.array([owners[int(m[i])]==face.source_key for i in ids]);return ids[onplane&owned]


def images(case,versions,n,owners,scores):
    labels=list(versions);origin=versions[labels[0]].vertices.min(axis=0);allv=versions[labels[0]].vertices-origin
    limits=np.array([allv.min(axis=0),allv.max(axis=0)]).T
    fig=plt.figure(figsize=(5*len(labels),5.2))
    for i,name in enumerate(labels):render3d(fig.add_subplot(1,len(labels),i+1,projection='3d'),versions[name],n,origin,limits,name)
    fig.suptitle(f'Case {case:02d}: identical orthographic camera, true XYZ aspect; terrain included');fig.tight_layout();fig.savefig(OUT/f'case-{case:02d}-whole.png',dpi=170);plt.close(fig)
    old=versions[labels[0]];faces=extract_surface_faces(old,owners);a,_,_,_=arrays(old)
    walls=np.flatnonzero((old.markers>=0)&(old.markers<n))
    problematic=int(walls[np.argmin(triangle_quality(old.vertices,old.faces)[walls])])
    strip=next(f for f in faces if f.kind=='wall' and problematic in f.original_faces)
    facade=max([f for f in faces if f.kind=='wall'],key=lambda f:np.sum(a[f.original_faces].min(axis=1)<10))
    # Show the original problematic corner at the same scale in every variant.
    worst=problematic;center=(old.vertices[old.faces[worst]]-origin).mean(axis=0);radius=2.5
    corner=np.array([[center[0]-radius,center[0]+radius],[center[1]-radius,center[1]+radius],[max(0,center[2]-radius),center[2]+radius]])
    fig=plt.figure(figsize=(5*len(labels),5.2))
    for i,name in enumerate(labels):render3d(fig.add_subplot(1,len(labels),i+1,projection='3d'),versions[name],n,origin,corner,name,zoom=(center,radius))
    fig.suptitle('Problematic corner: 5 m horizontal window, true XYZ aspect; identical magnification');fig.tight_layout();fig.savefig(OUT/f'case-{case:02d}-corner.png',dpi=180);plt.close(fig)
    ring=np.array(facade.polygon.exterior.coords);bounds=facade.polygon.bounds
    width,height=bounds[2]-bounds[0],bounds[3]-bounds[1]
    fig,axes=plt.subplots(len(labels),1,figsize=(max(10,min(18,width)),max(4,len(labels)*height/width*12)),squeeze=False)
    for ax,name in zip(axes[:,0],labels):
        mesh=versions[name];ids=plane_mask(mesh,n,facade,owners);v=(mesh.vertices-facade.origin)@facade.basis.T;p=v[mesh.faces[ids]];angles=arrays(mesh)[0][ids].min(axis=1)
        colors=np.where((angles<19.99999)[:,None],np.array([[.95,.5,.45,.8]]),np.array([[.91,.94,.96,1.]]))
        ax.add_collection(PolyCollection(p,facecolors=colors,edgecolors='#283d48',linewidths=.24));ax.set_xlim(bounds[0]-.1,bounds[2]+.1);ax.set_ylim(bounds[1]-.1,bounds[3]+.1);ax.set_aspect('equal',adjustable='box');ax.set_title(f'{name}: {len(ids)} facade triangles; {sum(angles<19.99999)} below 20°',fontsize=10);ax.set_xlabel('Distance along facade (m)');ax.set_ylabel('Height offset (m)')
    fig.suptitle('Wall plane: true 1:1 aspect. Red: minimum angle below 20°. All boundary vertices retained.');fig.tight_layout();fig.savefig(OUT/f'case-{case:02d}-wall-plane.png',dpi=180);plt.close(fig)
    # The narrowest facade is intrinsically tall. Keep its true aspect in
    # compact columns, and add a square physical-scale detail separately.
    for view in ['wall-strip','wall-detail']:
        fig,axes=plt.subplots(1,len(labels),figsize=(3.3*len(labels),10 if view=='wall-strip' else 4.5),squeeze=False)
        bounds=strip.polygon.bounds
        center=((old.vertices[old.faces[problematic]]-strip.origin)@strip.basis.T).mean(axis=0)
        for ax,name in zip(axes[0],labels):
            mesh=versions[name];ids=plane_mask(mesh,n,strip,owners);v=(mesh.vertices-strip.origin)@strip.basis.T
            ax.add_collection(PolyCollection(v[mesh.faces[ids]],facecolors='#e5eef3',edgecolors='#233d4c',linewidths=.35))
            if view=='wall-strip':
                ax.set_xlim(bounds[0]-.03,bounds[2]+.03);ax.set_ylim(bounds[1]-.03,bounds[3]+.03)
            else:
                ax.set_xlim(center[0]-.5,center[0]+.5);ax.set_ylim(center[1]-.5,center[1]+.5)
            ax.set_aspect('equal',adjustable='box');ax.set_title(name,fontsize=9);ax.set_xlabel('Along facade (m)',fontsize=8);ax.tick_params(labelsize=7)
        axes[0,0].set_ylabel('Height offset (m)')
        fig.suptitle('Worst original sliver facade: true 1:1 aspect'+(' — magnified 1 m square' if view=='wall-detail' else ' — full height'),fontsize=10)
        fig.tight_layout();fig.savefig(OUT/f'case-{case:02d}-{view}.png',dpi=190);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(16,4))
    for ax,part in zip(axes,['walls','roofs','terrain']):
        for name in labels:ax.stairs(scores[name][part]['angle_histogram'],BINS,label=name)
        ax.set_yscale('symlog',linthresh=1);ax.set_xlabel('Minimum angle (degrees)');ax.set_ylabel('Absolute triangle count');ax.set_title(part);ax.axvline(20,color='black',ls=':')
    axes[0].legend(fontsize=7);fig.tight_layout();fig.savefig(OUT/f'case-{case:02d}-angles.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(1,len(labels),figsize=(5*len(labels),5),squeeze=False)
    for ax,name in zip(axes[0],labels):
        mesh=versions[name];t=mesh.vertices[mesh.faces];m=np.asarray(mesh.markers);angle=arrays(mesh)[0].min(axis=1)
        terrain=t[m<0,:,:2]-origin[:2]
        ax.add_collection(PolyCollection(terrain,facecolors='#f1f4ec',edgecolors='#b8c2bc',linewidths=.12))
        for kind,mask,color in [('roof',m>=n,'#126fa0'),('terrain',m<0,'#b02b26')]:
            bad=mask&(angle<19.99999);centers=t[bad].mean(axis=1)[:,:2]-origin[:2]
            ax.scatter(centers[:,0],centers[:,1],s=5,c=color,label=f'{kind} <20° ({sum(bad)})',alpha=.8)
        ax.set_xlim(limits[0]);ax.set_ylim(limits[1]);ax.set_aspect('equal');ax.set_title(name);ax.legend(fontsize=7,loc='upper right');ax.set_xlabel('x offset (m)');ax.set_ylabel('y offset (m)')
    fig.suptitle('Spatial locations of every sub-target roof and terrain triangle (no subsampling)')
    fig.tight_layout();fig.savefig(OUT/f'case-{case:02d}-poor-locations.png',dpi=170);plt.close(fig)
    fig,axes=plt.subplots(1,len(labels),figsize=(5*len(labels),6),squeeze=False)
    sources=list(dict.fromkeys(owners));palette=plt.get_cmap('tab10')
    for ax,name in zip(axes[0],labels):
        mesh=versions[name];m=np.asarray(mesh.markers);t=mesh.vertices[mesh.faces]
        for k,source in enumerate(sources):
            mask=np.array([int(marker)>=n and owners[int(marker)-n]==source if 0<=marker<2*n else False for marker in m])
            xy=t[mask,:,:2]-origin[:2]
            ax.add_collection(PolyCollection(xy,facecolors=[palette(k)],edgecolors=(.15,.2,.25,.35),linewidths=.13,alpha=.7))
            if len(xy):
                domain=shapely.union_all(shapely.polygons(xy))
                for polygon in getattr(domain,'geoms',[domain]):
                    boundary=np.array(polygon.exterior.coords);ax.plot(boundary[:,0],boundary[:,1],color='#9c2025',lw=1.2)
                center=domain.representative_point();ax.text(center.x,center.y,str(k+1),fontsize=9,ha='center',bbox=dict(facecolor='white',alpha=.7,edgecolor='none'))
        ax.set_xlim(limits[0]);ax.set_ylim(limits[1]);ax.set_aspect('equal');ax.set_title(name);ax.set_xlabel('x offset (m)')
    axes[0,0].set_ylabel('y offset (m)')
    fig.suptitle('Source domains and shared interfaces: red outlines, stable ownership colours; original footprint geometry')
    fig.tight_layout();fig.savefig(OUT/f'case-{case:02d}-interfaces.png',dpi=170);plt.close(fig)
    return facade.id


def main():
    sections=[];summary={}
    for case in (1,11):
        oldcat=json.loads((WORK/f'catalogue-integration/case-{case:02d}.json').read_text())['catalogue']['audit']['roof_catalogue'];owners=[r['source_key'] for r in oldcat['patch_markers']];n=len(owners)
        versions={'Original catalogue':load_vtu(WORK/f'catalogue-integration/case-{case:02d}-catalogue.vtu'),'Fan repair':io.load_mesh(WORK/f'catalogue-walls-dtcc_mesher/case-{case:02d}-catalogue.dtcc')}
        records={}
        for backend in ['dtcc_mesher','triangle']:
            record=json.loads((OUT/f'case-{case:02d}-{backend}.json').read_text());records[backend]=record
            if record['audit']['roof_catalogue'].get('surface_refinement',{}).get('accepted'):versions[f'Pilot {backend}']=io.load_mesh(OUT/f'case-{case:02d}-{backend}.dtcc')
        scores={name:score(mesh,n) for name,mesh in versions.items()}
        facade=images(case,versions,n,owners,scores)
        perface={}
        for backend,record in records.items():
            audit=record['audit']['roof_catalogue'].get('surface_refinement',{});mesh=versions.get(f'Pilot {backend}')
            if mesh is None:continue
            ids=next(f.values for f in mesh.fields if f.name=='catalogue_surface_face')
            perface[backend]=[dict(face_id=r['face_id'],kind=r['kind'],source_key=r['source_key'],angles_and_sizes=metrics(mesh,n,ids==r['face_id'],spacing=r['spacing']),constraints=r['angle_constraint'],backend=r['backend'],geometry=r['geometry']) for r in audit.get('faces',[])]
        result=dict(scores=scores,records=records,per_face=perface,problem_facade_id=facade)
        (OUT/f'case-{case:02d}-comparison.json').write_text(json.dumps(result,indent=2));summary[str(case)]={name:row for name,row in scores.items()}
        baseline_runs={}
        for key in ['original','fan']:
            path=OUT/f'isolated-{key}/case-{case:02d}-dtcc_mesher.json'
            if path.exists():baseline_runs[key]=json.loads(path.read_text())
        result['isolated_baseline_resources']={k:dict(seconds=r['seconds'],peak_rss_bytes=r['peak_rss_bytes'],preservation=r['preservation']) for k,r in baseline_runs.items()}
        (OUT/f'case-{case:02d}-comparison.json').write_text(json.dumps(result,indent=2))
        rows=[]
        for name,score_ in scores.items():
            for part in ['walls','roofs','terrain']:
                s=score_[part];rows.append('<tr>'+''.join(f'<td>{html.escape(str(v))}</td>' for v in [name,part,s['triangles'],round(s['minimum_angle']['min'],5),s['poor']['5']['count'],s['poor']['10']['count'],s['poor']['20']['count'],round(s['poor']['20']['area_m2'],4),s['maximum_vertex_incidence'],round((s['adjacent_size_ratio'] or {}).get('p99',0),2)])+'</tr>')
        notes=[]
        source_rows=[]
        for k,source in enumerate(oldcat['sources']):
            source_rows.append('<tr>'+''.join(f'<td>{html.escape(str(value))}</td>' for value in [k+1,source['source_key'],source.get('shape'),source['outcome'],source.get('final_outcome'),'; '.join(source['diagnostics'].get('rejection_reasons',[]))])+'</tr>')
        notes.append('<table><tr><th>Map source</th><th>Stable key</th><th>Provisional shape</th><th>Accepted/rejected</th><th>Final outcome (unchanged)</th><th>Existing rejection reasons</th></tr>'+''.join(source_rows)+'</table>')
        for b,r in records.items():
            a=r['audit']['roof_catalogue'].get('surface_refinement',{});cat=r['audit']['roof_catalogue'];geom=[f['geometry'] for f in a.get('faces',[])]
            notes.append(f'<p><b>{b}</b>: accepted={a.get("accepted")}; {r["seconds"]:.2f} s; process peak {r["peak_rss_bytes"]/1024**2:.1f} MiB; {len(a.get("rounds",[]))} rounds; {a.get("added_boundary_vertices")} boundary / {a.get("added_interior_vertices")} interior additions. Strict topology valid={a.get("topology",{}).get("valid")}; shared-boundary mismatches={a.get("shared_boundary_mismatches")}. Error: {html.escape(str(a.get("error")))}. Preservation {html.escape(str(r["preservation"]))}. Retained {cat.get("retained_fits")} fits ({cat.get("retained_nonflat_fits")} non-flat), {len(cat.get("block_fallbacks",[]))} block fallbacks. Max plane/coverage error {max((g["plane_error_m"] for g in geom),default=0):.3g} m / {max((g["coverage_error_m2"] for g in geom),default=0):.3g} m².</p>')
        notes.append('<p>Isolated baseline resource reruns: '+html.escape(str(result['isolated_baseline_resources']))+'. Saved original and fan meshes match these reruns bitwise. Historical cumulative-process peaks are not compared as per-build allocation measurements.</p>')
        sections.append(f'<section><h2>Case {case:02d}: '+('hip building 14' if case==1 else 'merged [0, 3, 23, 32, 46]')+'</h2>'+''.join(notes)+'<p>Coverage triangulation is DTCC Mesher for all pilot variants; the named face backend performs the actual independent quality refinement.</p><table><tr><th>Variant</th><th>Surface</th><th>Triangles</th><th>Min angle °</th><th>&lt;5°</th><th>&lt;10°</th><th>&lt;20°</th><th>Area &lt;20° m²</th><th>Max incidence</th><th>Adjacent size p99</th></tr>'+''.join(rows)+'</table>'+''.join(f'<a href="case-{case:02d}-{view}.png"><img src="case-{case:02d}-{view}.png"></a>' for view in ['whole','corner','wall-plane','wall-strip','wall-detail','angles','poor-locations','interfaces'])+f'<p><a href="case-{case:02d}-comparison.json">Complete distributions, per-face constraints, geometry, topology and provenance</a></p></section>')
    intro='''<!doctype html><meta charset="utf-8"><title>Shared-edge surface refinement pilot</title><style>body{font:15px system-ui;max-width:1700px;margin:30px auto;color:#20333e;background:#f5f7f8}section{background:white;padding:22px;margin:25px 0}img{width:100%}table{border-collapse:collapse;font-size:12px}td,th{border:1px solid #c7d0d5;padding:7px}code{background:#eee}p{line-height:1.55}</style><h1>Shared-edge constrained refinement — two-case pilot</h1><p>20° face target; nominal 4 m building / 8 m terrain spacing, area caps 8 / 32 m². DTCC Mesher uses maintained constrained circumcentre/encroached-segment refinement with common dyadic corner spacing and shell reuse (protection threshold 20°). Boundary requests propagate through a canonical 3D edge registry. No custom Delaunay/off-centre repair and no proximity weld. Default city meshing is unchanged.</p><p>Budgets: 12 synchronization rounds, 180 s refinement, 100,000 refinement steps per face, 50,000 vertices per face, 300,000 total vertices. Supervised pilot processes: 600 s timeout, 2 GiB peak-memory audit (not an operating-system RSS cap). Failed candidates are retained transactionally as the original complete shell and explicitly marked rejected.</p><p>Neighbour size ratio is √(larger triangle area / smaller triangle area) across a shared edge. q = 4√3 area / sum(edge²). All counts use the complete stored meshes; no subsampling. Runtime measurements below exclude imports and input-cache loading.</p><p>Frozen train data only; no threshold tuning or holdout access. Original and fan baselines remain read-only. Roof domains/planes, footprint conditioning, source ownership and reconstruction decisions must agree before accepting any quality claim. Triangle is optional. Per-triangle DTCC exemptions are exported and counted separately; protection does not establish that every sub-target triangle is mathematically unavoidable.</p>'''
    conclusion="""<section><h2>Assessment and remaining limits</h2>
<p>The DTCC pilot resolves the facade fans: every wall triangle in both cases meets 20°, all non-exempt shell triangles meet 20°, and maximum wall incidence falls from 25 to 10 for the hip building. Adjacent wall size ratios have p99 approximately 2.29 / 2.14 and maxima below 3.17. The wireframes show progressive refinement rather than dense boundary vertices connected to distant sparse vertices.</p>
<p>The hip shell retains 53 sub-target roof triangles (0.0265 m²) and 302 terrain triangles (0.179 m²). The group retains 499 sub-target terrain triangles (6.262 m²). All are backend-exempt; all affected planar faces have an actual input corner below 20°. Original minimum input angles are about 2.57° on the hip roof, 1.42° on its terrain, and 3.46° on group terrain. A fixed corner below 20° cannot be filled entirely by triangles whose angles are at least 20°. Neighbour size ratios above 4 remain at protected regions only (17 roof + 41 terrain edges in case 01; 22 terrain edges in case 11), with maxima up to 7.98. Shield policy also exempts some neighbouring elements: their precise counts and lower angles are not all proved mathematically unavoidable. <a href="remaining-angle-reproducers.json">Reduced planar reproducers and coordinates</a>.</p>
<p>The cost is substantial: 1,378 → 20,661 triangles (15.0×) for the hip shell and 3,184 → 36,428 (11.4×) for the group, largely to retain the original terrain planes and grade short features. The pilot completes within its recorded time, vertex and memory budgets. This does not establish full-city scaling. Isolated baseline timings/peaks are listed above; prior cumulative full-run peaks are not treated as isolated per-build memory.</p>
<p>Optional Triangle refinement remains unsuccessful: both complete component replacements exhaust the 12-round synchronization budget and retain their original valid shells. This is an algorithm/synchronization limitation under the declared budget, not a proof of geometric impossibility and not a successful independent-backend validation. Failed runs and earlier protection experiments are preserved.</p>
<p><b>Recommendation:</b> retain this as an explicitly selected DTCC pilot. Upstream and pin the focused mesher extension, address Triangle synchronization and assess growth before requesting broader integration. No full-city rollout was performed. Default LOD1 and catalogue repair remain unchanged.</p></section>"""
    (OUT/'report.html').write_text(intro+''.join(sections)+conclusion+'<p><a href="validation.json">Validation and baseline hashes</a> · <a href="../../../docs/shared_edge_pilot.md">Implementation and reproduction</a></p>')
    (OUT/'quality-summary.json').write_text(json.dumps(summary,indent=2))

if __name__=='__main__':main()
