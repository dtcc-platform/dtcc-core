"""Terrain sampling and source-block foundations, independent of roof meshing."""
from collections import defaultdict

import numpy as np
import shapely
from scipy.ndimage import map_coordinates
from shapely import affinity

from ...model import Mesh
from ..geometry_policy import COINCIDENCE
from .roof_interfaces import RoofAssemblyError


def minimum_dem(polygon, raster):
    """Minimum of the clamped bilinear DEM over a polygon, including its holes.

    Bilinear interiors have no strict minimum. Check contained DEM nodes and
    quadratic extrema on each boundary interval between DEM grid lines.
    """
    data = np.asarray(raster.data, dtype=float)
    transform = raster.georef
    xs = transform.c+(np.arange(data.shape[1])+.5)*transform.a
    ys = transform.f+(np.arange(data.shape[0])+.5)*transform.e

    def values(xy):
        xy = np.asarray(xy)
        columns = np.clip((xy[:,0]-transform.c)/transform.a-.5,0,data.shape[1]-1)
        rows = np.clip((xy[:,1]-transform.f)/transform.e-.5,0,data.shape[0]-1)
        return map_coordinates(data,[rows,columns],order=1,mode='nearest',prefilter=False)

    xmin,ymin,xmax,ymax = polygon.bounds
    columns = np.flatnonzero((xs>=xmin)&(xs<=xmax))
    rows = np.flatnonzero((ys>=ymin)&(ys<=ymax))
    minimum = np.inf
    if len(columns) and len(rows):
        x,y = np.meshgrid(xs[columns],ys[rows])
        inside = shapely.intersects_xy(polygon,x,y)
        if np.any(inside): minimum = float(data[np.ix_(rows,columns)][inside].min())
    for piece in getattr(polygon,'geoms',[polygon]):
        for ring in [piece.exterior,*piece.interiors]:
            xy = np.asarray(ring.coords)
            for a,b in zip(xy[:-1],xy[1:]):
                delta = b-a
                cuts = [0.,1.]
                for grid,start,direction in ((xs,a[0],delta[0]),(ys,a[1],delta[1])):
                    if direction:
                        t = (grid-start)/direction
                        cuts.extend(t[(t>0)&(t<1)])
                cuts = np.unique(cuts)
                lo,hi = cuts[:-1],cuts[1:]
                z0 = values(a+lo[:,None]*delta)
                z1 = values(a+hi[:,None]*delta)
                zm = values(a+((lo+hi)/2)[:,None]*delta)
                qa = 2*(z0+z1-2*zm)
                qb = z1-z0-qa
                u = np.divide(-qb,2*qa,out=np.full_like(qa,-1.),where=qa>0)
                interior = (u>0)&(u<1)&(qa>0)
                minimum = min(minimum,float(z0.min()),float(z1.min()))
                if np.any(interior):
                    minimum = min(minimum,float((z0+qb*u+qa*u*u)[interior].min()))
    if not np.isfinite(minimum): raise RoofAssemblyError('Terrain foundation has no finite DEM support')
    return float(minimum)


def terrain_attachment(ground, patches, raster, *, offset, smoothing=0):
    """Sample outside terrain and attach its hole boundaries to block foundations.

    Roof interiors are not triangulated to establish terrain. All boundary
    geometry remains in the authoritative local XY subdivision; native terrain
    sampling alone receives world coordinates.
    """
    from ..geometry_builders import meshes as pipeline

    world = Mesh(vertices=np.array(ground.vertices,copy=True),faces=ground.faces,markers=ground.markers)
    world.vertices[:,:2] += offset
    world = pipeline._snap_ground_mesh_to_raster_bounds(world,raster)
    sampled = pipeline._dtcc_builder.build_terrain_surface_mesh_from_ground_mesh(
        pipeline.mesh_to_builder_mesh(world),pipeline.raster_to_builder_gridfield(raster),0).from_cpp()
    sampled.vertices[:,:2] = ground.vertices[:,:2]
    by_block = defaultdict(list)
    for patch in patches: by_block[patch.block_id].append(patch.polygon)
    blocks = sorted(by_block)
    domains = [shapely.union_all(by_block[block]) for block in blocks]
    bases = {block:minimum_dem(affinity.translate(domain,*offset),raster)
             for block,domain in zip(blocks,domains)}
    frozen = np.zeros(len(sampled.vertices),dtype=bool)
    attached = {}
    if domains:
        tree = shapely.STRtree([domain.boundary for domain in domains])
        pairs = tree.query(shapely.points(sampled.vertices[:,:2]),predicate='dwithin',distance=COINCIDENCE)
        for vertex,index in zip(*pairs):
            block = blocks[int(index)]
            attached.setdefault(int(vertex),set()).add(block)
            bases[block] = min(bases[block],float(sampled.vertices[vertex,2]))
        for vertex,owners in attached.items():
            if len(owners)!=1:
                raise RoofAssemblyError('Distinct cleaned blocks touch at a terrain vertex',block_ids=owners)
            sampled.vertices[vertex,2] = bases[next(iter(owners))]
            frozen[vertex] = True
    sampled.terrain_contact_vertices = np.flatnonzero(frozen)
    if smoothing:
        neighbours = [set() for _ in sampled.vertices]
        for a,b,c in sampled.faces:
            neighbours[a].update((b,c));neighbours[b].update((a,c));neighbours[c].update((a,b))
        for _ in range(smoothing):
            for i,adjacent in enumerate(neighbours):
                if adjacent and not frozen[i]:
                    sampled.vertices[i,2] = sampled.vertices[sorted(adjacent),2].mean()
    return sampled,bases,dict(block_bases_m=bases,terrain_triangles=len(sampled.faces),
        foundation_triangles=0,frozen_vertices=int(frozen.sum()),
        foundation_policy='minimum clamped bilinear DEM over the occupied block and sampled contact',
        xy_unchanged=bool(np.array_equal(sampled.vertices[:,:2],ground.vertices[:,:2])))
