"""Premeshing topology of admitted height sectors at multiway roof junctions."""
from collections import defaultdict
from dataclasses import replace

import numpy as np
import shapely
from shapely.geometry import LineString, Point, box
from shapely.ops import polygonize

from ..geometry_policy import COINCIDENCE, XY_PRECISION, canonical_profile_heights
from ..roofs.partition import polygon_parts
from .roof_interfaces import RoofAssemblyError, node_interfaces


def height_interval_topology(heights):
    """Count exposed vertical sheets between every consecutive incident level.

    For a manifold terrain graph each horizontal slice of the cyclic sector
    star has zero or two solid/empty transitions. Four transitions identify
    disconnected high sectors meeting on the same vertical edge. Equal levels
    are considered together; affine edge crossings are noded before this test.
    """
    levels = sorted(set(heights))
    intervals = []
    for low,high in zip(levels[:-1],levels[1:]):
        occupied = np.asarray(heights)>(low+high)/2
        sheets = int(np.count_nonzero(occupied!=np.roll(occupied,1)))
        intervals.append(dict(low_m=low,high_m=high,vertical_sheets=sheets))
    return dict(manifold=all(i['vertical_sheets']<=2 for i in intervals),intervals=intervals)


def junction_half_widths(directions,radius,clearance):
    """Separate rectangle corners from incident rays with aspect at most sqrt(3).

    Fold all rays into the first quadrant. The largest empty angular interval
    gives the maximin corner direction in [30,60] degrees. Numerical ties use
    the same coordinate uncertainty as the shared geometry, then prefer the
    square and the smaller direction. No triangulation is used in this choice.
    """
    folded=[float(np.arctan2(abs(np.sin(a)),abs(np.cos(a)))) for a in directions]
    low,high=np.pi/6,np.pi/3
    breaks=sorted({low,high,*[a for a in folded if low<a<high]})
    candidates=[low,high,*[(a+b)/2 for a,b in zip(breaks[:-1],breaks[1:])]]
    scores={a:min(abs(a-b) for b in folded) for a in candidates}
    uncertainty=max(64*np.finfo(float).eps,XY_PRECISION/(4*clearance))
    tied=[a for a in candidates if scores[a]>=max(scores.values())-uncertainty]
    nearest=min(abs(a-np.pi/4) for a in tied)
    angle=min(a for a in tied if abs(a-np.pi/4)<=nearest+uncertainty)
    return radius*np.array([np.cos(angle),np.sin(angle)]),dict(
        corner_direction_degrees=float(np.degrees(angle)),
        minimum_ray_corner_separation_degrees=float(np.degrees(scores[angle])))


def _stars(patches, base):
    segments, incident = [], defaultdict(set)
    edges = set()
    for p in patches:
        for ring in [p.polygon.exterior,*p.polygon.interiors]:
            xy=np.asarray(ring.coords)
            for a,b in zip(xy[:-1],xy[1:]):
                edge=tuple(sorted((tuple(a),tuple(b))))
                if edge in edges: continue
                edges.add(edge)
                i=len(segments);segments.append(LineString(edge))
                incident[edge[0]].add(i);incident[edge[1]].add(i)
    if not segments: return [],segments
    tree=shapely.STRtree(segments)
    polygons=[p.polygon for p in patches]
    regions=shapely.STRtree(polygons)
    slopes={i:p.coefficients[:2] for i,p in enumerate(patches)}
    stars=[]
    junctions=sorted(xy for xy,indices in incident.items() if len(indices)>3)
    for xy in junctions:
        origin=np.asarray(xy);point=Point(xy)
        indices=incident[xy]
        clearance=min(segments[i].length for i in indices)/4
        for i in tree.query(point.buffer(4*clearance)):
            if int(i) not in indices:
                distance=point.distance(segments[int(i)])
                if distance>COINCIDENCE:clearance=min(clearance,distance/4)
        for other in junctions:
            if other!=xy:clearance=min(clearance,float(np.linalg.norm(origin-other))/4)
        directions=[]
        for i in indices:
            a,b=np.asarray(segments[i].coords)
            d=b-origin if np.array_equal(a,origin) else a-origin
            directions.append(float(np.arctan2(d[1],d[0])%(2*np.pi)))
        directions.sort()
        sectors=[]
        for k,angle in enumerate(directions):
            following=directions[(k+1)%len(directions)]
            width=(following-angle)%(2*np.pi)
            mid=angle+width/2
            probe=Point(origin+clearance/8*np.array([np.cos(mid),np.sin(mid)]))
            labels=[int(i) for i in regions.query(probe,predicate='covered_by')]
            if len(labels)>1:
                raise RoofAssemblyError('Admitted roof sectors overlap at a junction',block_ids=[p.block_id for p in patches])
            sectors.append(labels[0] if labels else None)
        present=set(sectors)-{None}
        values=[(float(patches[i].height(origin)),i) for i in sorted(present)]
        canonical=canonical_profile_heights(values,slopes)
        heights=[base if i is None else canonical[i] for i in sectors]
        topology=height_interval_topology(heights)
        stars.append(dict(xy=list(xy),valence=len(indices),clearance_m=clearance,
                          directions_rad=directions,sectors=sectors,heights_m=heights,**topology))
    return stars,segments


def regularize_roof_junctions(patches, *, bases, detail):
    """Change only nonmanifold admitted stars, before any mesh is generated.

    One bounded central rectangle separates alternating height
    sectors. Harmless multiway incidence retains its exact geometry. The same
    unchanged supporting planes cover the edited domains; source support,
    occupied extent and resulting interval topology are checked here.
    """
    grouped=defaultdict(list)
    for p in patches:grouped[p.block_id].append(p)
    result=[];audits={}
    for block,parts in grouped.items():
        stars,_=_stars(parts,bases[block])
        polygons=[p.polygon for p in parts]
        union=shapely.union_all(polygons)
        centres=[];changes=[]
        for star in stars:
            if star['manifold']:continue
            radius=min(detail/2,star['clearance_m'])
            if radius<=16*XY_PRECISION:
                raise RoofAssemblyError(f'Nonmanifold roof junction {star["xy"]} cannot be separated within the detail budget',block_ids=[block])
            x,y=star['xy']
            half,construction=junction_half_widths(star['directions_rad'],radius,star['clearance_m'])
            # Quantize the new construction vertices, then remove the GEOS
            # precision model. The occupied outline is immutable: clipping it
            # with a fixed-precision operand would round a second outline copy.
            centre=shapely.set_precision(shapely.set_precision(
                box(x-half[0],y-half[1],x+half[0],y+half[1]),XY_PRECISION),0.)
            touched=[i for i,p in enumerate(polygons) if p.intersection(centre).area>0]
            owner=min(touched,key=lambda i:(-polygons[i].intersection(centre).area,parts[i].source_key,parts[i].patch_id))
            centres.append((centre,owner))
            changes.append(dict(xy=star['xy'],valence=star['valence'],radius_m=radius,
                construction_half_widths_m=half.tolist(),**construction,
                owner=parts[owner].source_key,
                reason='nonmanifold vertical height intervals',height_intervals=star['intervals']))
        edited=parts
        if changes:
            # Derive every intersection once in the shared graph. Independently
            # clipping incident polygons creates competing copies of short edges.
            linework=shapely.union_all([p.boundary for p in polygons])
            for centre,owner in centres:linework=shapely.difference(linework,centre)
            graph=shapely.union_all([linework,union.boundary,*[c.boundary for c,_ in centres]])
            tree=shapely.STRtree(polygons)
            owned=defaultdict(list)
            for piece in polygonize(list(getattr(graph,'geoms',[graph]))):
                probe=piece.representative_point()
                if not union.covers(probe):continue
                label=next((i for centre,i in centres if centre.covers(probe)),None)
                if label is None:
                    labels=[int(i) for i in tree.query(probe,predicate='covered_by')]
                    if len(labels)!=1:
                        raise RoofAssemblyError('Junction subdivision has ambiguous roof ownership',block_ids=[block])
                    label=labels[0]
                owned[label].append(piece)
            if set(owned)!=set(range(len(parts))):
                raise RoofAssemblyError('Junction ownership would erase an admitted roof patch',block_ids=[block])
            edited=[replace(parts[i],polygon=piece) for i,pieces in sorted(owned.items()) for piece in pieces]
            for i,pieces in owned.items():
                for piece in pieces:
                    if min(float(np.min(parts[i].height(np.asarray(r.coords)))) for r in
                           [piece.exterior,*piece.interiors])<=bases[block]+COINCIDENCE:
                        raise RoofAssemblyError('Junction ownership intersects an admitted roof foundation',block_ids=[block])
            edited=node_interfaces(edited,canonical=True)
            after,_=_stars(edited,bases[block])
            invalid=[star for star in after if not star['manifold']]
            if invalid:
                raise RoofAssemblyError(f'Junction ownership has nonmanifold height intervals in block {block}: {invalid}',block_ids=[block])
            if shapely.union_all([p.polygon for p in edited]).symmetric_difference(union).area>max(1e-10,union.area*1e-12):
                raise RoofAssemblyError('Junction ownership changed occupied extent',block_ids=[block])
            final={i:shapely.union_all(pieces) for i,pieces in owned.items()}
            # Measure actual source domains after construction. A central region
            # may extend over terrain, and transfers between patches of the same
            # source do not change that source's occupied area.
            source_before=defaultdict(list)
            source_after=defaultdict(list)
            for patch in parts:source_before[patch.source_key].append(patch.polygon)
            for patch in edited:source_after[patch.source_key].append(patch.polygon)
            source_differences={}
            for source,domains in source_before.items():
                before=shapely.union_all(domains)
                after=shapely.union_all(source_after[source])
                source_differences[source]=(after.difference(before),before.difference(after))
            for change,(centre,owner) in zip(changes,centres):
                touched=[i for i,p in enumerate(polygons) if p.intersection(centre).area>0]
                change['source_changes']={source:dict(
                    added_m2=source_differences[source][0].intersection(centre).area,
                    removed_m2=source_differences[source][1].intersection(centre).area)
                    for source in sorted({parts[i].source_key for i in touched})}
                displacement=max(polygons[i].boundary.hausdorff_distance(final[i].boundary) for i in touched)
                lengths=[];angles=[]
                for i in touched:
                    for piece in polygon_parts(final[i]):
                        for ring in [piece.exterior,*piece.interiors]:
                            points=np.asarray(ring.coords)[:-1]
                            for k,point in enumerate(points):
                                if not centre.buffer(COINCIDENCE).covers(Point(point)) or polygons[i].boundary.distance(Point(point))<=COINCIDENCE:continue
                                a=points[k-1]-point;b=points[(k+1)%len(points)]-point
                                lengths.extend([float(np.linalg.norm(a)),float(np.linalg.norm(b))])
                                angles.append(float(np.degrees(np.arctan2(abs(a[0]*b[1]-a[1]*b[0]),a@b))))
                change.update(maximum_displacement_m=displacement,
                    minimum_new_corner_angle_degrees=min(angles,default=None),
                    minimum_new_corner_edge_m=min(lengths,default=None))
                if displacement>detail/2+2*XY_PRECISION:
                    raise RoofAssemblyError('Junction ownership exceeds the physical displacement budget',block_ids=[block])
        audits[block]=dict(policy='admitted cyclic height intervals; bounded central ownership only when nonmanifold',
            checked_junctions=[dict(s,sector_sources=[None if i is None else parts[i].source_key for i in s['sectors']]) for s in stars],
            changes=changes,displacement_budget_m=detail/2)
        result.extend(edited)
    return result,audits
