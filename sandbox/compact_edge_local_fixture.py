"""Find deterministic generic endpoint-star controls, with remote acute components."""
import json
from pathlib import Path
import sys
import numpy as np
import dtcc_mesher as dm
from compact_edge_synthetic import measure


def fixture(index):
    rng=np.random.default_rng(350+index)
    count=10; height=(2.,4.,8.)[index%3]
    x=np.r_[0,np.sort(rng.uniform(.1,19.9,count-1)),20]
    p=np.array([[v,0] for v in x]+[[20,height],[0,height]])
    m=dm.mesh(dm.Domain.from_loops(p),options=dm.MeshingOptions(refine=False,acute_protection='none'))
    top=int(np.flatnonzero(np.all(m.points[m.segments,1]==height,axis=1))[0])
    for i in range(12):
        k=len(m.points);x=100+2*i
        m.points=np.vstack([m.points,[[x,0],[x+1,0],[x,.000001]]])
        m.triangles=np.vstack([m.triangles,[[k,k+1,k+2]]])
        m.segments=np.vstack([m.segments,[[k,k+1],[k+1,k+2],[k+2,k]]])
    return m,top


if __name__=='__main__':
    rows={}
    for i in range(60):
        m,top=fixture(i)
        r=dm.improve(m,max_added_points=1,splittable_segments=[top],allow_interior=False)
        rows[i]=dict(parameter=r.segment_parameter,kind=r.insertion_kind,before=measure(m),after=measure(r.mesh),
                     points=m.points.tolist(),triangles=m.triangles.tolist(),segments=m.segments.tolist(),top=top,
                     evaluations=r.candidate_evaluations)
    Path(sys.argv[1]).write_text(json.dumps(rows,indent=2))
