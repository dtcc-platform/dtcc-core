"""Equal one-point generic controls; select a backend using PYTHONPATH."""
import argparse
import json
from pathlib import Path
import numpy as np
import dtcc_mesher as dm


def fixtures():
    for count in (4,6,10,20):
        for height in (2.,4.,6.,10.):
            for remote in (0,12):
                p=np.array([[x,0] for x in np.linspace(0,12,count+1)]+[[12,height],[0,height]])
                m=dm.mesh(dm.Domain.from_loops(p),options=dm.MeshingOptions(refine=False,acute_protection='none'))
                top=int(np.flatnonzero(np.all(m.points[m.segments,1]==height,axis=1))[0])
                for i in range(remote):
                    k=len(m.points);x=100+2*i
                    m.points=np.vstack([m.points,[[x,0],[x+1,0],[x,.00001]]])
                    m.triangles=np.vstack([m.triangles,[[k,k+1,k+2]]])
                    m.segments=np.vstack([m.segments,[[k,k+1],[k+1,k+2],[k+2,k]]])
                yield f'fan-{count}-{height:g}-remote-{remote}',m,top


def measure(m):
    p=m.points[m.triangles];e=np.roll(p,-1,axis=1)-p
    area=abs(e[:,0,0]*e[:,1,1]-e[:,0,1]*e[:,1,0])/2
    q=4*np.sqrt(3)*area/(e*e).sum(axis=(1,2))
    angle=[]
    for k in range(3):
        u=p[:,(k+1)%3]-p[:,k];v=p[:,(k+2)%3]-p[:,k]
        angle.append(np.degrees(np.arctan2(abs(u[:,0]*v[:,1]-u[:,1]*v[:,0]),(u*v).sum(axis=1))))
    a=np.min(angle,axis=0)
    return dict(triangles=len(p),old_penalty=float(np.maximum(0,.25-q).sum()),
        energy=float(((1-q)**2).sum()),minimum_q=float(q.min()),
        poor={str(k):dict(count=int(sum(a<k-1e-7)),area=float(area[a<k-1e-7].sum())) for k in (5,10,20)})


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('output',type=Path);args=parser.parse_args()
    rows={}
    for name,m,top in fixtures():
        r=dm.improve(m,max_added_points=1,splittable_segments=[top],allow_interior=False)
        rows[name]=dict(before=measure(m),after=measure(r.mesh),kind=r.insertion_kind,
                        parameter=r.segment_parameter,evaluations=r.candidate_evaluations)
    args.output.write_text(json.dumps(rows,indent=2))
