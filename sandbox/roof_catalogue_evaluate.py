"""Compare the roof catalogue with saved DTCC/Roofer evaluation inputs."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import pickle
import sys
import time

import numpy as np
from shapely.ops import unary_union

from roof_catalogue import fit_catalogue


def select_cases(originals, baseline, targets):
    """Select diagnostic examples from input geometry, before catalogue fitting."""
    by = {(b['area'],b['index']):b for b in baseline}
    cases, used = [], set()

    def add(area, index, reason):
        if (area,index) not in used:
            cases.append(dict(area=area,index=index,reason=reason))
            used.add((area,index))

    for area,index in [('lund',146),('gothenburg',14),('lund',201)]:
        add(area,index,'Previously inspected comparison case')
    for area, buildings in originals.items():
        items = []
        for i,b in enumerate(buildings):
            poly = b.lod0.to_polygon()
            items.append((i,poly.area,poly.area/poly.minimum_rotated_rectangle.area,
                          len(poly.interiors),by[area,i]['train_points']))
        rectangular = sorted([t for t in items if t[2]>.9 and t[4]>60],key=lambda t:t[1])
        add(area,rectangular[len(rectangular)//2][0],'Typical rectangular footprint by median area')
        nonrect = sorted([t for t in items if t[1]>80 and t[4]>60],key=lambda t:t[2])
        add(area,nonrect[0][0],'Strongly irregular footprint')
        sparse = sorted([t for t in items if t[4]>=12],key=lambda t:t[4]/t[1])
        add(area,sparse[0][0],'Sparse points per footprint area')
        holes = sorted([t for t in items if t[3]],key=lambda t:t[1],reverse=True)
        if holes:
            add(area,holes[0][0],'Courtyard / holed footprint')
    largest = max(baseline,key=lambda b:b['area_m2'])
    add(largest['area'],largest['index'],'Largest footprint in the sample')
    for area, buildings in originals.items():
        polygons = [targets[area][.5].get(b.id) for b in buildings]
        neighbors = {i:set() for i in range(len(polygons))}
        for i,p in enumerate(polygons):
            if p is None:
                continue
            for j in range(i):
                if polygons[j] is not None and p.distance(polygons[j])<.002:
                    neighbors[i].add(j)
                    neighbors[j].add(i)
        groups, seen = [], set()
        for i in neighbors:
            if i in seen:
                continue
            group, todo = set(), [i]
            while todo:
                j = todo.pop()
                if j in group:
                    continue
                group.add(j)
                todo.extend(neighbors[j]-group)
            seen |= group
            if 2<=len(group)<=6:
                groups.append(sorted(group))
        group = max(groups,key=lambda g:sum(buildings[i].lod0.to_polygon().area for i in g))
        cases.append(dict(area=area,indices=group,
                          reason='Connected cleaned-footprint group (0.5 m); largest 2–6 member component'))
    return cases


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(args.input.resolve()))
    from evaluate import metrics
    baseline = json.loads((args.input/'results.json').read_text())['buildings']
    originals = {a:pickle.load(open(args.input/a/'buildings.pkl','rb')) for a in ['gothenburg','lund']}
    targets = pickle.load(open(args.input/'targets.pkl','rb'))
    (args.output/'cases.json').write_text(json.dumps(select_cases(originals,baseline,targets),indent=2))
    model_path = Path(__file__).with_name('roof_catalogue.py')
    (args.output/'experiment.json').write_text(json.dumps(dict(
        input=str(args.input.resolve()), model_sha256=hashlib.sha256(model_path.read_bytes()).hexdigest(),
        baseline=json.loads((args.input/'experiment.json').read_text()),
        catalogue='flat, shed, symmetric gable, hip', inner_validation='2 m cells, (2*ix+iy)%5==1',
        model_selection='weighted pseudo-Huber validation error + 0.025 m per parameter',
        outer_holdout_used_for_fitting=False),indent=2))
    high = {(b['area'],b['index']):b for b in json.loads((args.input/'sensitivity-high-results.json').read_text())}
    rows, models, patches, transfer = [], {}, {}, []
    for b in baseline:
        area, index = b['area'], b['index']
        key = f'{area}/{index}'
        original = originals[area][index]
        footprint = original.lod0.to_polygon()
        cloud = np.load(args.input/area/f'b{index:03d}'/'points.npz')
        start = time.perf_counter()
        model = fit_catalogue(footprint, cloud['train'], original.attributes['ground_height'])
        elapsed = time.perf_counter()-start
        models[key] = model
        ps = model.patches(footprint) if model else []
        patches[key] = ps
        cat = dict(seconds=elapsed, model=model.kind if model else None,
                   diagnostics=model.diagnostics if model else {}, **metrics(ps, cloud['test'], footprint))
        if ps:
            union = unary_union([p['polygon'] for p in ps])
            cat['partition_overlap_m2'] = max(0.,sum(p['polygon'].area for p in ps)-union.area)
        rows.append({**b, 'roofer_high':high.get((area,index),{'patches':0}), 'catalogue':cat})
        if model:
            for delta in [.5,1.]:
                target = targets[area][delta].get(b['id'])
                if target is None:
                    continue
                tp = model.patches(target)
                u = unary_union([p['polygon'] for p in tp])
                transfer.append(dict(area=area,index=index,delta=delta,
                                     coverage=u.intersection(target).area/target.area,
                                     overlap=max(0.,sum(p['polygon'].area for p in tp)-u.area),
                                     added_area=target.difference(footprint).area,
                                     supported=model.diagnostics['supported']))
        if index%20==0:
            print(area,index,'model',cat['model'],'seconds',round(elapsed,2),flush=True)
    dump = lambda obj,path:path.write_text(json.dumps(obj,indent=2,default=lambda x:x.item()))
    dump(rows,args.output/'results.json')
    pickle.dump(models,open(args.output/'models.pkl','wb'))
    pickle.dump(patches,open(args.output/'patches.pkl','wb'))
    dump(transfer,args.output/'transfer.json')
    summary = dict(n=len(rows), catalogue_types=dict(Counter(b['catalogue']['model'] for b in rows)),
                   supported=sum(b['catalogue']['diagnostics'].get('supported',False) for b in rows),
                   runtime_seconds=sum(b['catalogue']['seconds'] for b in rows))
    for label, sample in [('all',rows),('paired_dtcc',[b for b in rows if all('median_ae' in b[m] for m in ['dtcc','roofer','catalogue'])]),
                          ('supported_catalogue',[b for b in rows if b['catalogue']['diagnostics'].get('supported')])]:
        eligible = [b for b in sample if b['test_points']>=20]
        methods = {}
        for m in ['dtcc','roofer','roofer_high','catalogue']:
            available = [b for b in eligible if 'median_ae' in b[m]]
            methods[m]=dict(roofs=sum(b[m]['patches']>0 for b in sample), n=len(available),
                            median_ae=float(np.median([b[m]['median_ae'] for b in available])),
                            p90_ae=float(np.median([b[m]['p90_ae'] for b in available])),
                            good_observation_fit=sum(b[m]['median_ae']<=.25 and b[m]['p90_ae']<=1 and b[m]['point_coverage']>=.99 for b in available))
        summary[label]=methods
    summary['transfer'] = {str(d):dict(n=sum(t['delta']==d for t in transfer),
        complete=sum(t['delta']==d and t['coverage']>.99999 and t['overlap']<1e-5 for t in transfer)) for d in [.5,1.]}
    dump(summary,args.output/'summary.json')
    print(json.dumps(summary,indent=2))


if __name__ == '__main__':
    main()
