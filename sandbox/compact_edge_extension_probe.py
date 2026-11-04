"""Uncommitted joint diagnostics of every finite native edge candidate."""
import json
import pickle
import subprocess
from types import SimpleNamespace
from compact_edge_extension import OUT
from dtcc_core.builder.meshing import compact_improvement as ci


def probe(case):
    state=pickle.load(open(OUT/f'case-{case:02d}-checkpoint.pkl','rb'))
    registry,caches,faces=(state[x] for x in ('registry','caches','faces'))
    rows=[];fixtures=[]
    for f in faces:
        if f.kind=='terrain':continue
        c=caches[f.id];allowed=ci.permitted_segments(f,c,registry,faces)
        for segment,edge in sorted(allowed.items()):
            text=f"{len(c['points'])} {len(c['triangles'])} {len(c['segments'])} {segment}\n"
            text+='\n'.join(' '.join(map(str,row)) for array in (c['points'],c['triangles'],c['segments']) for row in array)+'\n'
            result=subprocess.run([str(OUT/'edge-probe')],input=text,text=True,capture_output=True,check=True)
            fixtures.append(dict(face=f.id,edge=edge,segment=segment,
                incident={str(i):dict(points=caches[i]['points'].tolist(),triangles=caches[i]['triangles'].tolist(),
                segments=caches[i]['segments'].tolist(),kind=faces[i].kind,source=faces[i].source_key)
                for i in registry.edges[edge].incidents}))
            winner=ci.improve(c['points'],c['triangles'],c['segments'],max_added_points=1,
                              splittable_segments=[segment],allow_interior=False)
            for line in result.stdout.splitlines():
                values=list(map(float,line.split()));t=values[0]
                try:
                    _,record=ci.edge_transaction(registry,caches,faces,f,
                        SimpleNamespace(parent_segment=segment,segment_parameter=t),allowed,backend='dtcc_mesher')
                except (ValueError,RuntimeError) as exc:
                    record=dict(accepted=False,reason='invalid incident retessellation: '+str(exc))
                rows.append(dict(face=f.id,edge=edge,t=t,native_penalty_before=values[1],native_penalty_after=values[2],
                    native_returned_t=winner.segment_parameter if winner.insertion_kind=='segment' else None,
                    joint=record))
    (OUT/f'case-{case:02d}-extension-edge-evidence.json').write_text(json.dumps(dict(fixtures=fixtures,trials=rows),indent=2))
    print(case,len(rows),'joint accepted',sum(r['joint']['accepted'] for r in rows),flush=True)


if __name__=='__main__':
    for case in (1,11):probe(case)
