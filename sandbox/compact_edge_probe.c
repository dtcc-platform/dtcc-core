/* Evaluation-only driver for the native bounded candidate set. Not a public API. */
#include "improve.c"

int main(void)
{
    size_t np,nt,ns,parent,i,j,n;dtcc_mesher_mesh input={0};TMMesh base={0};double offset[2];
    Candidate proposals[CANDIDATE_LIMIT];size_t nf=0;
    if(scanf("%zu %zu %zu %zu",&np,&nt,&ns,&parent)!=4) return 2;
    input.num_points=np;input.num_triangles=nt;input.num_segments=ns;
    input.points=calloc(np,sizeof(*input.points));input.triangles=calloc(nt*3,sizeof(*input.triangles));input.segments=calloc(ns,sizeof(*input.segments));
    if(!input.points || !input.triangles || !input.segments) return 3;
    for(i=0;i<np;i++) if(scanf("%lf %lf",&input.points[i].x,&input.points[i].y)!=2) return 2;
    for(i=0;i<nt*3;i++) if(scanf("%u",&input.triangles[i])!=1) return 2;
    for(i=0;i<ns;i++) if(scanf("%u %u",&input.segments[i].a,&input.segments[i].b)!=2) return 2;
    tm_initialize();if(import_mesh(&input,&base,offset)!=TM_OK || parent>=ns) return 4;
    if(flips(&base,&nf)!=TM_OK) return 5;
    n=candidates(&base,&parent,1,0,proposals);
    for(i=0;i<n;i++) {
        TMMesh trial={0};int point=-1;TMStatus status=clone(&base,&trial);size_t flips_count=0;
        if(status!=TM_OK) return 3;
        for(j=0;j<trial.triangle_count;j++) {
            int e=tm_find_edge_in_triangle(&trial.triangles[j],base.segments[parent].v[0],base.segments[parent].v[1]);
            if(e>=0) {status=tm_insert_point_on_edge(&trial,(int)j,e,proposals[i].p,TM_VERTEX_SEGMENT_SPLIT,&point);break;}
        }
        if(status==TM_OK) status=tm_split_segment_registry(&trial,parent,point);
        if(status==TM_OK) status=flips(&trial,&flips_count);
        if(status==TM_OK) {
            dtcc_mesher_shape_score a=score(&base),b=score(&trial);
            printf("%.17g %.17g %.17g %.17g %.17g %zu %.17g %zu %.17g\n",proposals[i].t,
                a.penalty,b.penalty,a.minimum_q,b.minimum_q,b.poor_5_count,b.poor_5_area,b.poor_10_count,b.poor_10_area);
        }
        tm_free_mesh(&trial);
    }
    tm_free_mesh(&base);dtcc_mesher_mesh_free(&input);return 0;
}
