"""Resume frozen shared geometry with a separate, edge-only additional budget."""
import argparse
from copy import deepcopy
import hashlib
import importlib
import json
from pathlib import Path
import pickle
import resource
import time
from unittest.mock import patch

import numpy as np
from dtcc_core import io
from dtcc_core.builder.meshing.compact_surface import CompactSurfaceOptions
from roof_wall_evaluation import _triangle_signature

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'.agent/work/compact-edge-extension'


def run(case,budget,experiment,*,single_trial=False):
    with (OUT/f'case-{case:02d}-checkpoint.pkl').open('rb') as stream:
        state=pickle.load(stream)
    ci=importlib.import_module('dtcc_core.builder.meshing.compact_improvement')
    cs=importlib.import_module('dtcc_core.builder.meshing.compact_surface')
    improve_shared=ci.improve_shared
    native_improve=ci.improve
    def native_control(*args,**kwargs):
        result=native_improve(*args,**kwargs)
        if single_trial:
            result.segment_trials=()
        return result
    extension={}
    initial=deepcopy(state)
    ceiling=85 if case==1 else 307
    options=CompactSurfaceOptions(max_face_vertices=0,max_source_vertices=4,max_group_vertices=budget)

    def resume(registry,caches,faces,unused):
        caches.clear();caches.update(deepcopy(state['caches']))
        faces[:]=deepcopy(state['faces'])
        registry=deepcopy(state['registry'])
        registry,audit=improve_shared(registry,caches,faces,options,max_building_triangles=ceiling)
        extension.update(audit)
        # Assembly consumes the complete history; the additional-budget audit
        # remains separate from the fixed starting points and their charges.
        total=deepcopy(audit)
        total['optional_vertices']=state['audit']['optional_vertices']+audit['optional_vertices']
        total['source_optional_vertices']={s:n+audit['source_optional_vertices'][s]
                                          for s,n in state['audit']['source_optional_vertices'].items()}
        extension['old_points_bitwise_fixed']=np.array_equal(
            np.asarray(registry.vertices)[:len(initial['registry'].vertices)],initial['registry'].vertices)
        extension['old_points_used']=all(set(initial['caches'][i]['ids'])<=set(c['ids'])
            for i,c in caches.items() if not c['unchanged'])
        extension['no_new_interiors']=all(c['optional']==initial['caches'][i]['optional'] for i,c in caches.items())
        return registry,total

    started=time.perf_counter()
    with patch.object(ci,'improve_shared',resume),patch.object(ci,'improve',native_control):
        mesh,audit=cs.compact_surface(state['reference'],state['patches'],options)
    baseline=io.load_mesh(OUT.parent/f'compact-mesher-interface/case-{case:02d}-shared-dtcc_mesher.dtcc')
    terrain_fixed=np.array_equal(_triangle_signature(mesh,mesh.markers<0),
                                _triangle_signature(baseline,baseline.markers<0))
    assert terrain_fixed and extension['old_points_bitwise_fixed'] and extension['old_points_used'] and extension['no_new_interiors']
    assert np.count_nonzero(mesh.markers>=0)<=ceiling
    if budget==0:
        assert np.array_equal(mesh.vertices,baseline.vertices) and np.array_equal(mesh.faces,baseline.faces)
    import dtcc_mesher._core as native
    record=dict(case=case,budget=budget,experiment=experiment,seconds=time.perf_counter()-started,
        peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        extension=extension,audit=audit,terrain_bitwise_fixed=terrain_fixed,
        checkpoint_sha256=hashlib.sha256((OUT/f'case-{case:02d}-checkpoint.pkl').read_bytes()).hexdigest(),
        binary_sha256=hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest())
    directory=OUT/experiment;directory.mkdir(parents=True,exist_ok=True)
    name=f'case-{case:02d}-extra-{budget}'
    mesh.save(directory/f'{name}.dtcc')
    (directory/f'{name}.json').write_text(json.dumps(record,indent=2))
    print(case,budget,len(extension['optional_vertices']),np.count_nonzero(mesh.markers>=0),
          record['seconds'],extension['termination'],flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--experiment',required=True)
    parser.add_argument('--case',type=int,choices=[1,11])
    parser.add_argument('--single-trial',action='store_true',help='Control: hide returned alternatives, retaining the native local winner')
    args=parser.parse_args()
    for case in ([args.case] if args.case else [1,11]):
        for budget in ([0,2,4] if case==1 else [0,2,4,8]):
            run(case,budget,args.experiment,single_trial=args.single_trial)
