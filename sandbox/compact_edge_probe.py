"""Joint evidence for native edge candidates that do not produce a returned trial.

The evaluation C driver calls Mesher's private finite generator and primitives;
Core only stages the supplied parent/t on its ordinary shared-edge coordinator.
"""
import json
from pathlib import Path
import subprocess
from types import SimpleNamespace
from unittest.mock import patch

from compact_mesher_evaluation import evaluate,WORK
from dtcc_core.builder.meshing import compact_improvement as ci

OUT=WORK/'compact-mesher-interface'


def main():
    original=ci.improve_shared
    def capture(registry,caches,faces,options):
        rows=[];seen=set();fixtures=[]
        for f in faces:
            if f.kind=='terrain':continue
            c=caches[f.id];allowed=ci.permitted_segments(f,c,registry,faces)
            for segment,edge in sorted(allowed.items()):
                key=(f.id,edge)
                if key in seen:continue
                seen.add(key)
                text=f"{len(c['points'])} {len(c['triangles'])} {len(c['segments'])} {segment}\n"
                text+='\n'.join(' '.join(map(str,row)) for array in (c['points'],c['triangles'],c['segments']) for row in array)+'\n'
                result=subprocess.run([str(OUT/'edge-probe')],input=text,text=True,capture_output=True,check=True)
                incident=sorted(registry.edges[edge].incidents)
                fixtures.append(dict(face_id=f.id,edge_id=edge,segment=segment,
                    incident_faces=incident,face_inputs={str(i):dict(points=caches[i]['points'].tolist(),triangles=caches[i]['triangles'].tolist(),
                        segments=caches[i]['segments'].tolist(),kind=faces[i].kind,source_key=faces[i].source_key) for i in incident}))
                for line in result.stdout.splitlines():
                    values=list(map(float,line.split()));t=values[0]
                    native=dict(parameter=t,penalty_before=values[1],penalty_after=values[2],q_before=values[3],q_after=values[4],
                                poor_5_count=int(values[5]),poor_5_area=values[6],poor_10_count=int(values[7]),poor_10_area=values[8])
                    try:
                        _,record=ci.edge_transaction(registry,caches,faces,f,SimpleNamespace(parent_segment=segment,segment_parameter=t),allowed,backend=options.backend)
                    except (ValueError,RuntimeError) as exc:
                        record=dict(accepted=False,reason='invalid incident retessellation: '+str(exc),
                                    affected_faces=incident,triangle_cost=len(incident))
                    rows.append(dict(proposer=f.id,edge_id=edge,native=native,joint=record))
                    if len(rows)>=256:break
                if len(rows)>=256:break
            if len(rows)>=256:break
        (OUT/'case-01-edge-evidence.json').write_text(json.dumps(dict(fixtures=fixtures,trials=rows,
            limit=256,scope='initial compact hip geometry, native finite edge candidates; diagnostic transactions are never committed'),indent=2))
        print('edge diagnostic trials',len(rows),'joint beneficial',sum(r['joint']['accepted'] for r in rows),flush=True)
        return original(registry,caches,faces,options)
    with patch.object(ci,'improve_shared',capture):
        evaluate(1,'dtcc_mesher',OUT/'edge-probe-run','shared')

if __name__=='__main__':main()
