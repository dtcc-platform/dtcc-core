"""Bounded source-junction policy for a manifold 2.5D building subdivision."""
from collections import defaultdict

import numpy as np
import shapely
from shapely import affinity
from shapely.geometry import LineString, Point
from shapely.ops import polygonize
from ..geometry_policy import XY_PRECISION


def _boundary_lines(geometry):
    if geometry.geom_type in ('LineString','LinearRing'):
        yield geometry
    else:
        for part in getattr(geometry,'geoms',[]):
            yield from _boundary_lines(part)


def canonical_boundary_graph(graph, domains, audit):
    """Insert the same canonical samples into graph edges and cleaned outlines."""
    from scipy.spatial import cKDTree

    # The labelled cleaner subdivision includes internal source interfaces as
    # well as the occupied outline. Node these canonical grid boundaries once.
    parent_graph = shapely.union_all([domain.boundary for domain in domains.values()])
    segments = []
    for line in _boundary_lines(parent_graph):
        xy = np.asarray(line.coords)
        segments.extend((a,b) for a,b in zip(xy[:-1],xy[1:]))
    lines = [LineString([a,b]) for a,b in segments]
    tree = shapely.STRtree(lines)
    corners = np.unique(np.asarray(segments).reshape(-1,2),axis=0)
    corner_tree = cKDTree(corners)
    graph_lines = [np.asarray(line.coords) for line in _boundary_lines(graph)]
    mapped, samples = {}, {}
    maximum = 0.
    changed = 0
    for xy in np.unique(np.concatenate(graph_lines),axis=0):
        key = tuple(xy)
        distance, corner = corner_tree.query(xy)
        if distance <= 2*XY_PRECISION:
            mapped[key] = corners[corner]
            continue
        candidates = []
        for index in tree.query(Point(xy).buffer(2*XY_PRECISION)):
            a,b = segments[int(index)]
            d = b-a
            t = float(np.clip((xy-a)@d/(d@d),0.,1.))
            projected = a+t*d
            candidates.append((float(np.linalg.norm(projected-xy)),int(index),t))
        selected = min(candidates) if candidates else None
        if selected is not None and selected[0]<=2*XY_PRECISION:
            samples.setdefault(selected[1],[]).append((selected[2],key))
        else:
            mapped[key] = xy
    outlines = []
    for index,(a,b) in enumerate(segments):
        # Alias only indistinguishable samples on the same immutable edge.
        # Inserting these samples in the outline itself avoids nearly coincident
        # full-length outline segments competing with the projected graph.
        previous = a
        for t,key in sorted(samples.get(index,[])):
            projected = a+t*(b-a)
            if np.linalg.norm(projected-previous)>XY_PRECISION:
                outlines.append(LineString([previous,projected]))
                previous = projected
            mapped[key] = previous
        if not np.array_equal(previous,b):
            outlines.append(LineString([previous,b]))
    canonical = []
    boundary = shapely.union_all(lines)
    boundary_band = boundary.buffer(3*XY_PRECISION)
    replaced_segments = 0
    for coordinates in graph_lines:
        points = [mapped[tuple(xy)] for xy in coordinates]
        for a,b in zip(points[:-1],points[1:]):
            if np.array_equal(a,b):
                continue
            segment = LineString([a,b])
            # The outline supplies these edges and every canonical sample.
            # Retaining a second floating-point copy makes spurious thin cells
            # even when its endpoints are mathematically on the same line.
            if boundary_band.covers(segment):
                replaced_segments += 1
            else:
                canonical.append(segment)
        for before,after in zip(coordinates,points):
            distance = float(np.linalg.norm(after-before))
            maximum = max(maximum,distance)
            changed += distance>1e-10
    if audit is not None:
        audit.update(replaced_outline_segments=replaced_segments,projected_boundary_vertices=int(changed),maximum_boundary_projection_m=maximum,
                     boundary_projection_budget_m=3*XY_PRECISION)
    return shapely.union_all(canonical+outlines)


def canonical_source_subdivision(cells, *, domain, source_keys):
    """Normalize source cells against one occupied outline without ownership edits.

    Multiway XY incidence alone does not justify a geometry change. Roof-height
    topology is evaluated separately after roof fitting and domain admission.
    """
    from .construction import polygon_parts
    offset = np.floor(np.asarray(domain.bounds[:2]))
    original_domain = affinity.translate(domain, *-offset)
    normalized_domain = shapely.set_precision(original_domain, XY_PRECISION)
    if not cells:
        return affinity.translate(normalized_domain,*offset), [], dict(changes=[])
    by_label = defaultdict(list)
    for polygon, label in cells:
        by_label[label].append(affinity.translate(polygon, *-offset))
    raw = {label:shapely.union_all(parts) for label, parts in by_label.items()}
    originals = {label:shapely.union_all(polygon_parts(polygon.intersection(original_domain)))
                 for label,polygon in raw.items()}
    regions = {label:shapely.union_all(polygon_parts(shapely.set_precision(polygon,XY_PRECISION).intersection(normalized_domain)))
               for label,polygon in raw.items()}
    excluded = []
    for label in list(regions):
        if regions[label].is_empty:
            if originals[label].area > XY_PRECISION*originals[label].length:
                raise ValueError(f'Numerical domain clipping would erase source {source_keys[label]}')
            excluded.append(dict(source_key=source_keys[label],area_m2=float(originals[label].area),
                                 reason='no two-dimensional interior after numerical domain clipping'))
            del regions[label];del originals[label]
    if any(p.is_empty or not p.is_valid for p in regions.values()):
        raise ValueError('Cleaned source subdivision collapses at numerical overlay precision')
    if not regions:
        raise ValueError("Selected occupied domain has no two-dimensional source support")
    union = shapely.union_all(list(regions.values()))
    original_union = shapely.union_all(list(originals.values()))
    normalized_union = shapely.set_precision(union, 0.)
    normalization_area = original_union.symmetric_difference(normalized_union).area
    normalization_budget = original_union.length*XY_PRECISION + np.pi*XY_PRECISION**2
    displacement = max(p.boundary.hausdorff_distance(regions[label].boundary)
                       for label, p in originals.items())
    if normalization_area > normalization_budget or displacement > 2*XY_PRECISION:
        raise ValueError('Source subdivision normalization exceeds its propagated numerical budget')
    changes = []
    final_union = shapely.union_all(list(regions.values()))
    if final_union.symmetric_difference(union).area > max(1e-10, union.area*1e-12):
        raise ValueError('Cleaned subdivision junction policy changed occupied coverage')
    # Construct the labelled subdivision against the immutable occupied outline.
    # Grid-clipped source edges must not independently redefine that outline.
    source_graph = shapely.union_all([p.boundary for p in regions.values()],grid_size=XY_PRECISION)
    boundary_audit = {}
    graph = canonical_boundary_graph(source_graph,{0:normalized_domain},boundary_audit)
    labels = sorted(regions)
    region_tree = shapely.STRtree([regions[label] for label in labels])
    coordinated = defaultdict(list)
    for polygon in polygonize(list(_boundary_lines(graph))):
        if not normalized_domain.covers(polygon.representative_point()):
            continue
        candidates = region_tree.query(polygon,predicate='intersects')
        overlap = sorted(((regions[labels[int(i)]].intersection(polygon).area,labels[int(i)])
                          for i in candidates),reverse=True)
        if not overlap or overlap[0][0]<=0:
            candidates = region_tree.query(polygon.buffer(3*XY_PRECISION))
            if not len(candidates):
                raise ValueError('Canonical occupied source cell has no label')
            owner = min((regions[labels[int(i)]].distance(polygon),labels[int(i)]) for i in candidates)[1]
        else:
            owner = overlap[0][1]
        coordinated[owner].append(polygon)
    if set(coordinated)!=set(regions):
        raise ValueError('Canonical subdivision would erase a source label')
    regions = {label:shapely.union_all(parts) for label,parts in coordinated.items()}
    final_union = shapely.union_all(list(regions.values()))
    domain_change = final_union.symmetric_difference(normalized_domain).area
    if domain_change>max(1e-10,normalized_domain.area*1e-12):
        raise ValueError('Canonical source cells do not cover the selected occupied domain')
    domain_budget = (original_domain.length+original_union.length)*3*XY_PRECISION
    result = []
    for label, region in sorted(regions.items()):
        for polygon in getattr(region, 'geoms', [region]):
            if polygon.geom_type != 'Polygon':
                raise ValueError('Cleaned subdivision contains a non-polygon source region')
            result.append((affinity.translate(polygon, *offset), label))
    return affinity.translate(normalized_domain,*offset), result, dict(boundary_coordination=boundary_audit, selected_domain_area_change_m2=float(domain_change), selected_domain_area_budget_m2=float(domain_budget), changes=changes, excluded_sources=excluded, occupied_area_change_m2=float(final_union.symmetric_difference(union).area),
                        normalization_area_m2=float(normalization_area),
                        normalization_displacement_m=float(displacement),
                        normalization_area_budget_m2=float(normalization_budget),
                        policy='canonical source subdivision; no planar-valence ownership changes',
                        displacement_budget_m=2*XY_PRECISION)
