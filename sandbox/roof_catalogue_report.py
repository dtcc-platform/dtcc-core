"""Render the catalogue experiment, individual examples, and grouped roofs."""
import argparse
import html
import json
import os
from pathlib import Path
import pickle
import sys

os.environ.setdefault('MPLCONFIGDIR','/private/tmp/roof-catalogue-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
import shapely
from shapely.ops import unary_union

from roof_catalogue import polygon_parts


def triangles(patches):
    output = []
    for p in patches:
        for part in polygon_parts(p['polygon']):
            for tri in shapely.constrained_delaunay_triangles(part).geoms:
                xy = np.array(tri.exterior.coords)[:3,:2]
                z = (xy-p['origin'][:2])@p['coeff'][:2]+p['origin'][2]+p['coeff'][2]
                output.append(np.c_[xy,z])
    return output


def export_obj(path, tris, origin):
    with path.open('w') as f:
        f.write('# Roof faces only. Local origin is recorded in cases-details.json.\n')
        for tri in tris:
            for v in tri-origin:
                f.write('v %.6f %.6f %.6f\n'%tuple(v))
        for i in range(len(tris)):
            f.write(f'f {3*i+1} {3*i+2} {3*i+3}\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    out,src = args.output,args.input
    (out/'figures').mkdir(exist_ok=True)
    (out/'meshes').mkdir(exist_ok=True)
    sys.path.insert(0,str(src.resolve()))
    from evaluate import metrics
    summary = json.loads((out/'summary.json').read_text())
    rows = json.loads((out/'results.json').read_text())
    lookup = {(b['area'],b['index']):b for b in rows}
    base = pickle.load(open(src/'patches.pkl','rb'))
    high = pickle.load(open(src/'sensitivity-high-patches.pkl','rb'))
    cat = pickle.load(open(out/'patches.pkl','rb'))
    models = pickle.load(open(out/'models.pkl','rb'))
    targets = pickle.load(open(src/'targets.pkl','rb'))
    buildings = {a:pickle.load(open(src/a/'buildings.pkl','rb')) for a in ['gothenburg','lund']}
    cases = json.loads((out/'cases.json').read_text())
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    names = ['Points','DTCC','Roofer 0.888','Roofer 0.99','Catalogue']
    methods = ['dtcc','roofer','roofer_high','catalogue']
    details = []
    contact = plt.figure(figsize=(16,4*((len(cases)+3)//4)),layout='constrained')

    def draw(ax,ps,cloud,footprints,origin,bounds,zrange,norm,view):
        lo,hi = zrange
        if cloud is not None:
            if view=='3d':ax.scatter(cloud[:,0]-origin[0],cloud[:,1]-origin[1],cloud[:,2]-origin[2],c=cloud[:,2],norm=norm,cmap='viridis',s=2,alpha=.8,rasterized=True)
            else:ax.scatter(cloud[:,0]-origin[0],cloud[:,1]-origin[1],c=cloud[:,2],norm=norm,cmap='viridis',s=5,alpha=.85,rasterized=True)
        elif ps:
            local, facecolors = [], []
            light = np.array([-.4,-.5,1.]);light /= np.linalg.norm(light)
            for p in ps:
                normal = np.r_[-p['coeff'][:2],1.];normal /= np.linalg.norm(normal)
                shade = .35+.65*max(0.,float(normal@light))
                tris = triangles([p])
                local.extend(t-origin for t in tris)
                facecolors.extend([np.array([.38,.72,.82])*shade]*len(tris))
            if view=='3d':ax.add_collection3d(Poly3DCollection(local,facecolors=facecolors,edgecolors=facecolors,linewidths=.1,antialiased=False))
            else:ax.add_collection(PolyCollection([t[:,:2] for t in local],facecolors=facecolors,edgecolors=facecolors,linewidths=.1,antialiased=False))
            for p in ps:
                for part in polygon_parts(p['polygon']):
                    for ring in [part.exterior,*part.interiors]:
                        xy=np.array(ring.coords);z=(xy-p['origin'][:2])@p['coeff'][:2]+p['origin'][2]+p['coeff'][2]
                        if view=='3d':ax.plot(xy[:,0]-origin[0],xy[:,1]-origin[1],z-origin[2],c='#263d47',lw=.6)
                        else:ax.plot(xy[:,0]-origin[0],xy[:,1]-origin[1],c='#263d47',lw=.6)
        elif view=='3d':ax.text2D(.2,.5,'No accepted roof',transform=ax.transAxes,color='#934e2c')
        else:ax.text(.5,.5,'No accepted roof',transform=ax.transAxes,ha='center',color='#934e2c')
        if view!='3d':
            for footprint in footprints:
                for part in polygon_parts(footprint):
                    for ring in [part.exterior,*part.interiors]:
                        xy=np.asarray(ring.coords)-origin[:2];ax.plot(xy[:,0],xy[:,1],c='#8b9499',lw=.6,ls='--')
        xmin,ymin,xmax,ymax=bounds
        ax.set_xlim(xmin-origin[0],xmax-origin[0]);ax.set_ylim(ymin-origin[1],ymax-origin[1])
        ax.set_xlabel('x (m)',labelpad=0);ax.set_ylabel('y (m)',labelpad=0)
        ax.xaxis.set_major_locator(MaxNLocator(3));ax.yaxis.set_major_locator(MaxNLocator(3));ax.tick_params(labelsize=7,pad=0)
        if view=='3d':
            ax.set_zlim(lo-origin[2],hi-origin[2]);ax.set_zlabel('z (m)',labelpad=0)
            ax.zaxis.set_major_locator(MaxNLocator(3));ax.set_box_aspect((xmax-xmin,ymax-ymin,(hi-lo)*1.7));ax.view_init(50,-65)
        else:ax.set_aspect('equal')

    for number,case in enumerate(cases,1):
        area=case['area'];indices=case.get('indices',[case.get('index')]);keys=[f'{area}/{i}' for i in indices]
        label=f'{area.title()} '+(', '.join('#'+str(i) for i in indices))
        ident=f'case-{number:02d}';footprints=[buildings[area][i].lod0.to_polygon() for i in indices]
        clouds=[np.load(src/area/f'b{i:03d}'/'points.npz') for i in indices]
        cloud=np.vstack([p['all'] for p in clouds]);test=np.vstack([p['test'] for p in clouds])
        footprint=unary_union(footprints);bounds=np.array(footprint.bounds)+[-1,-1,1,1]
        origin=np.r_[np.array(footprint.centroid.coords[0]),min(buildings[area][i].attributes['ground_height'] for i in indices)]
        allpatches=[sum((base[k]['dtcc'] for k in keys),[]),sum((base[k]['roofer'] for k in keys),[]),sum((high.get(k,[]) for k in keys),[]),sum((cat[k] for k in keys),[])]
        zs=[v[2] for ps in allpatches for tri in triangles(ps) for v in tri]
        zmin=min(np.quantile(cloud[:,2],.01),min(zs,default=np.inf))-.3;zmax=max(np.quantile(cloud[:,2],.99),max(zs,default=-np.inf))+.3
        norm=Normalize(zmin,zmax);typelist=[models[k].kind if models[k] else 'none' for k in keys]
        supported=[bool(models[k] and models[k].diagnostics['supported']) for k in keys]
        status='supported' if all(supported) else 'REVIEW: weak or unsupported fit'
        fig=plt.figure(figsize=(18,8),layout='constrained')
        for col,name in enumerate(names):
            ps=allpatches[col-1] if col else None
            stats=metrics(ps,test,footprint) if col else None
            title=name
            if col==4:title+=' · '+('/'.join(typelist) if len(typelist)==1 else f'{sum(supported)}/{len(supported)} supported')
            if stats and 'median_ae' in stats:title+=f"\nmedian / p90: {stats['median_ae']:.2f} / {stats['p90_ae']:.2f} m"
            if col==0:title+=f'\n{len(cloud)} points · {len(test)} held out'
            for row in range(2):
                ax=fig.add_subplot(2,5,row*5+col+1,projection='3d' if row==0 else None)
                draw(ax,ps,cloud if col==0 else None,footprints,origin,bounds,(zmin,zmax),norm,'3d' if row==0 else '2d')
                if row==0:ax.set_title(title,fontsize=10)
            if col:export_obj(out/'meshes'/f'{ident}-{methods[col-1]}.obj',triangles(ps),origin)
        fig.suptitle(f'{number:02d} · {label}\n{case["reason"]} · catalogue {status}',fontsize=14)
        fig.supxlabel('Roof faces only · common coordinates · points coloured by elevation, roofs shaded by plane orientation · vertical exaggeration ×1.7',fontsize=10)
        fig.savefig(out/'figures'/f'{ident}.png',dpi=140);plt.close(fig)
        ax=contact.add_subplot((len(cases)+3)//4,4,number,projection='3d')
        draw(ax,allpatches[-1],None,footprints,origin,bounds,(zmin,zmax),norm,'3d')
        ax.set_xticks([]);ax.set_yticks([]);ax.set_zticks([])
        ax.set_xlabel('');ax.set_ylabel('');ax.set_zlabel('')
        ax.set_title(f'{number:02d} · {label}\n'+('/'.join(typelist) if len(typelist)==1 else 'building group')+(' ✓' if all(supported) else ' · review'),fontsize=9)
        record=dict(number=number,label=label,source_ids=[lookup[area,i]['id'] for i in indices],
                    selection=case['reason'],origin=origin.tolist(),models=typelist,supported=supported,
                    diagnostics=[models[k].diagnostics if models[k] else {} for k in keys],
                    metrics={m:metrics(ps,test,footprint) for m,ps in zip(methods,allpatches)},figure=f'figures/{ident}.png')
        if len(indices)>1:
            cleanfoot=[targets[area][.5][buildings[area][i].id] for i in indices]
            cleanps=sum((models[k].patches(p) if models[k] else [] for k,p in zip(keys,cleanfoot)),[])
            fig=plt.figure(figsize=(10,5),layout='constrained')
            for col,(ps,fps,title) in enumerate([(allpatches[-1],footprints,'Source footprints'),(cleanps,cleanfoot,'Cleaned source cells (0.5 m)')],1):
                ax=fig.add_subplot(1,2,col,projection='3d');draw(ax,ps,None,fps,origin,bounds,(zmin,zmax),norm,'3d');ax.set_title(title)
            fig.suptitle(label+' · same fitted planes on both footprints',fontsize=12)
            fig.savefig(out/'figures'/f'{ident}-transfer.png',dpi=160);plt.close(fig)
            export_obj(out/'meshes'/f'{ident}-catalogue-cleaned.obj',triangles(cleanps),origin)
            record['transfer_figure']=f'figures/{ident}-transfer.png'
        details.append(record)
        print('rendered',ident,label,flush=True)
    contact.suptitle('Catalogue prototype · 13 selected buildings and groups\nReview flags use training data only; roof faces shown without walls',fontsize=16)
    contact.savefig(out/'catalogue-gallery.png',dpi=160);plt.close(contact)
    (out/'cases-details.json').write_text(json.dumps(details,indent=2,default=lambda x:x.item()))

    fig,axs=plt.subplots(1,3,figsize=(14,4),layout='constrained')
    labels=['DTCC','Roofer','Roofer 0.99','Catalogue'];colors=['#287b8b','#be7135','#8362a8','#38764c']
    vals=[summary['all'][m]['roofs'] for m in methods];bars=axs[0].bar(labels,vals,color=colors);axs[0].bar_label(bars,padding=3);axs[0].set_ylim(0,375);axs[0].set_title('Roof outputs / 341 inputs')
    vals=[summary['paired_dtcc'][m]['median_ae'] for m in methods];bars=axs[1].bar(labels,vals,color=colors);axs[1].bar_label(bars,fmt='%.3f',padding=3);axs[1].set_ylim(0,max(vals)*1.3);axs[1].set_ylabel('Vertical error (m)');axs[1].set_title('Typical median error · same 83 buildings')
    vals=[summary['paired_dtcc'][m]['p90_ae'] for m in methods];bars=axs[2].bar(labels,vals,color=colors);axs[2].bar_label(bars,fmt='%.3f',padding=3);axs[2].set_ylim(0,max(vals)*1.3);axs[2].set_title('Typical p90 error · same 83 buildings')
    for ax in axs:ax.tick_params(axis='x',rotation=15,labelsize=9)
    fig.savefig(out/'comparison.png',dpi=180);plt.close(fig)

    def fmt(v):return f'{v:.3f}'
    table='<table><tr><th>Method</th><th>Roof outputs</th><th>Median error, paired 83</th><th>P90 error, paired 83</th></tr>'
    for m,name in zip(methods,labels):
        table+=f'<tr><td>{name}</td><td>{summary["all"][m]["roofs"]}/341</td><td>{fmt(summary["paired_dtcc"][m]["median_ae"])} m</td><td>{fmt(summary["paired_dtcc"][m]["p90_ae"])} m</td></tr>'
    table+='</table>'
    sections=[]
    for case in details:
        stats='<table><tr><th>Method</th><th>Median error</th><th>P90 error</th><th>Held-out point coverage</th></tr>'
        for m,name in zip(methods,labels):
            v=case['metrics'][m];stats+=f'<tr><td>{name}</td><td>{fmt(v["median_ae"]) if "median_ae" in v else "—"}</td><td>{fmt(v["p90_ae"]) if "p90_ae" in v else "—"}</td><td>{100*v.get("point_coverage",0):.1f}%</td></tr>'
        stats+='</table>'
        review='; '.join(sorted(set(reason for d in case['diagnostics'] for reason in d.get('review_reasons',[])))) or 'All fitted members pass the prototype support checks.'
        meshes=' · '.join(f'<a href="meshes/case-{case["number"]:02d}-{m}.obj">{name} OBJ</a>' for m,name in zip(methods,labels))
        sections.append(f'<section id="case-{case["number"]}"><h2>{case["number"]:02d} · {case["label"]}</h2><p>{case["selection"]}. Catalogue: {", ".join(case["models"])}.</p><p class="note">{html.escape(review)}</p><a href="{case["figure"]}"><img loading="lazy" src="{case["figure"]}" alt="Five-method roof comparison for {case["label"]}"></a>{stats}<p>{meshes}</p>'+ (f'<img loading="lazy" src="{case["transfer_figure"]}" alt="Catalogue roofs transferred to cleaned cells">' if 'transfer_figure' in case else '')+'</section>')
    report=f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Roof catalogue prototype</title><style>body{{font:17px/1.6 system-ui,sans-serif;max-width:1400px;margin:36px auto;padding:0 24px;color:#24343d}}h1,h2{{line-height:1.2}}h2{{margin-top:40px}}img{{width:100%;border:1px solid #dbe2e5}}table{{border-collapse:collapse;width:100%}}td,th{{text-align:left;padding:9px;border-bottom:1px solid #dbe2e5}}th{{background:#eef3f5}}.note{{background:#fff5e8;padding:12px;border-left:4px solid #bc7b38}}a{{color:#13697d}}code{{background:#eef3f5;padding:2px 5px}}nav{{display:flex;flex-wrap:wrap;gap:8px 20px}}section{{padding-top:10px;border-top:1px solid #dbe2e5}}</style>
    <h1>Roof catalogue prototype: comparison on the same 341 buildings</h1>
    <p>Flat, shed, symmetric gable and hip models were fitted to exactly the saved training points used by DTCC and standalone Roofer. The catalogue improves some roofs dramatically and reduces large error tails compared with Roofer, but does not beat DTCC on the 83 buildings all methods reconstruct and can evaluate. Irregular footprints and compound roofs remain substantial limitations.</p>
    {table}<p>The errors are medians across buildings of per-building median and 90th-percentile absolute vertical error. Roof outputs include provisional catalogue fits; they are not equivalent to accepted, validated reconstructions. Of 335 catalogue candidates, {summary['supported']} pass training-based support checks. Six inputs have fewer than 12 training points and are skipped.</p>
    <img src="comparison.png" alt="Coverage and paired errors across four methods">
    <h2>Method and validation</h2>
    <p>Model degrees of freedom are flat: 1 (height); shed: 3 (height and two slope components); gable: 3 (orientation, eave height, pitch); hip: 4 (orientation, eave height, side and end pitch). Gable ridges are centred in the oriented bounding rectangle. Hip ridge length follows from the pitches and rectangle. These models impose a shared eave height on that rectangle and are clipped to the actual polygon, retaining courtyards. Concave or irregular boundaries need not have level eaves.</p>
    <p>Candidate orientations come from dominant footprint edges, the bounding rectangle and PCA, with local angle refinement of ±7.5°. Fitting uses soft-L1 residuals and inverse point counts in 1 m cells to reduce density imbalance. Slopes are bounded by 60°. Model selection uses an inner 2 m spatial-block holdout plus a 0.025 m penalty per parameter; 27 small or poorly distributed clouds cannot support this inner split and use penalized training error. The selected model is refitted to all training points. The original outer held-out points are used only for evaluation. Thresholds were not tuned to maximize the outer scores.</p>
    <p>Support checks require optimizer convergence, sufficient inliers on substantial roof faces, training median error ≤0.35 m, training p90 ≤1.5 m, and roof clearance above ground ≥0.25 m. These are provisional diagnostics, not calibrated confidence probabilities. A candidate can cover its whole footprint and still be unsupported. Review flags remain visible in the gallery.</p>
    <p class="note">The observations contain non-ground clutter and are not surveyed roof truth. The original DTCC statistical filter ran before the outer split. The source clouds have no building class; the Roofer adapter relabelled extracted points as class 6. A clean-looking catalogue roof is not proof that it is correct. This experiment neither validates complete solids nor produces a conforming city mesh.</p>
    <h2>What changed compared with Roofer</h2>
    <p>Across all evaluable outputs, the catalogue has 64 buildings meeting the illustrative held-out screen (median ≤0.25 m, p90 ≤1 m and ≥99% point coverage), compared with 58 DTCC, 47 default Roofer and 50 Roofer 0.99. Those counts mix different accepted subsets and should be read alongside the paired table. The prototype's broad coverage comes with many provisional fits; support checks accept only {summary['supported']}.</p>
    <p>Gothenburg #14 is a clear success: the catalogue selects a hip roof with roughly 5 cm median and 10 cm p90 error. On previously inspected Lund #146 and #201 it selects flat approximations and flags both for review. Those are useful counterexamples to adopting the catalogue as a universal replacement. Compound-roof decomposition or additional families would be needed, but additional flexibility should be justified with new validation data.</p>
    <h2>Transfer and grouped buildings</h2>
    <p>Each model is a continuous lower envelope of a few planes. Intersecting their analytic domains with a cleaned source cell provides full coverage without inventing a nearest-face assignment. All 334 candidates with retained cleaned cells pass the 2D coverage/overlap checks at both 0.5 m and 1 m cleaning. One fitted source building is removed by cleaning. These checks do not guarantee bounded extrapolation, good roof truth, or 3D mesh quality.</p>
    <p>The two group examples preserve a separate model per source building and show transfer onto cleaned cells. Shared roof heights can still differ or cross; step walls and common edge noding remain necessary for <code>city_surface_mesh</code>. The prototype returns roof planes and polygon patches, not production Building shells.</p>
    <h2>Thirteen visual cases</h2><p>Eleven individual buildings and two five-building groups were selected using the previously discussed cases and footprint/point-cloud characteristics: typical rectangular footprints, strongly irregular outlines, low point density, courtyards, and connected groups. This is a diagnostic gallery, not a random sample or a claim of balanced roof-type ground truth. Every panel uses common coordinates; points are coloured by elevation and roof faces are shaded by plane orientation; each case has perspective and plan views. Perspective views use vertical exaggeration ×1.7. Points beyond their 1st/99th elevation percentiles may be clipped when model extrema do not expand the view. Coloured surfaces are roof faces only; walls are omitted. Dashed plan outlines mark the input footprints. Group errors are pooled point-level errors and can cover different subsets when a method skips members; point coverage is reported.</p>
    <nav>{''.join(f'<a href="#case-{c["number"]}">{c["number"]:02d} {c["label"]}</a>' for c in details)}</nav>
    <p><a href="catalogue-gallery.png">Open the catalogue-only overview</a> · <a href="summary.json">Aggregate results</a> · <a href="results.json">Per-building results</a> · <a href="cases-details.json">Example IDs and diagnostics</a> · <a href="transfer.json">Transfer checks</a></p>
    {''.join(sections)}
    <h2>Reproduce and inspect</h2><p>The implementation is in <code>sandbox/roof_catalogue.py</code>; run <code>sandbox/roof_catalogue_evaluate.py --input .agent/work/roofer-evaluation --output .agent/work/catalogue-evaluation</code> with the repository Python environment, then <code>sandbox/roof_catalogue_report.py</code> with the same arguments. The report renderer uses Shapely 2.1 constrained triangulation. Original DTCC and Roofer outputs were reused unchanged. Roof model parameters, diagnostics and polygon patches are saved in <code>models.pkl</code> and <code>patches.pkl</code>. Each example has portable OBJ roof meshes in local metres, with the common origin recorded in <code>cases-details.json</code>. Baseline production code is unchanged.</p></html>'''
    (out/'report.html').write_text(report)
    print(out/'report.html',flush=True)


if __name__=='__main__':
    main()
