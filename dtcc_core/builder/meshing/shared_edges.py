"""Canonical 3D segment identities and monotone shared-boundary discretization."""
from collections import defaultdict
from dataclasses import dataclass, field

import numpy as np


@dataclass
class SharedEdge:
    """Canonical segment endpoints, incident faces and ordered shared split vertices."""
    endpoints: tuple[int, int]
    incidents: set[int] = field(default_factory=set)
    parameters: list[float] = field(default_factory=lambda: [0., 1.])
    vertices: list[int] = field(default_factory=list)


class SharedEdges:
    """Node topologically connected overlapping segments, then share every split.

    Matching uses endpoint identities and 3D collinearity, never XY rounding.
    Original endpoint vertices are immutable. Optional backend proposals within
    numerical tolerance reconcile deterministically along their actual edge.
    """
    def __init__(self, vertices, face_boundaries, *, tolerance=1e-8):
        self.vertices = [np.array(v, dtype=float) for v in vertices]
        self.tolerance = float(tolerance)
        self.edges = []
        self.face_edges = {}
        segments = sorted(set(tuple(sorted(e)) for boundary in face_boundaries.values() for e in boundary))
        incident = defaultdict(list)
        for i, (a, b) in enumerate(segments):
            incident[a].append(i)
            incident[b].append(i)
        # Collinear segment components sharing actual endpoint identities.
        nodes, seen = {}, set()
        for i, (a, b) in enumerate(segments):
            if i in seen:
                continue
            start, end = np.asarray(vertices)[[a,b]]
            delta = end-start
            length2 = float(delta@delta)
            if length2 <= tolerance*tolerance:
                raise ValueError('Original geometric edge is numerically collapsed')
            todo, component, endpoints = [i], [], set()
            seen.add(i)
            while todo:
                j = todo.pop()
                component.append(j)
                endpoints.update(segments[j])
                for vertex in segments[j]:
                    for k in incident[vertex]:
                        if k in seen:
                            continue
                        xy = np.asarray(vertices)[list(segments[k])]-start
                        projection = np.outer(xy@delta/length2, delta)
                        if np.max(np.linalg.norm(xy-projection, axis=1)) <= tolerance:
                            todo.append(k)
                            seen.add(k)
            for j in component:
                x, y = segments[j]
                p, q = np.asarray(vertices)[[x,y]]
                d = q-p
                l2 = float(d@d)
                chain = []
                for vertex in endpoints:
                    t = float((vertices[vertex]-p)@d/l2)
                    residual = np.linalg.norm(vertices[vertex]-(p+t*d))
                    if -1e-12 <= t <= 1+1e-12 and residual <= tolerance:
                        chain.append((t, vertex))
                chain.sort()
                nodes[(x,y)] = [v for _,v in chain]
        lookup = {}
        for face_id, boundary in face_boundaries.items():
            refs = []
            for a, b in boundary:
                key = tuple(sorted((a,b)))
                chain = nodes[key]
                if a != key[0]:
                    chain = chain[::-1]
                for x, y in zip(chain, chain[1:]):
                    edge_key = tuple(sorted((x,y)))
                    if edge_key not in lookup:
                        lookup[edge_key] = len(self.edges)
                        self.edges.append(SharedEdge(edge_key, vertices=list(edge_key)))
                    ei = lookup[edge_key]
                    self.edges[ei].incidents.add(face_id)
                    refs.append((ei, x == edge_key[0]))
            self.face_edges[face_id] = refs

    def chain(self, edge_id, forward=True):
        vertices = self.edges[edge_id].vertices
        return vertices if forward else vertices[::-1]

    def boundary(self, face_id):
        segments = []
        for ei, forward in self.face_edges[face_id]:
            chain = self.chain(ei, forward)
            segments.extend(zip(chain,chain[1:]))
        return segments

    def reconcile(self, requests, *, max_vertices):
        """Apply a batch atomically; return all faces incident on new splits."""
        proposals = defaultdict(list)
        for edge_id, t in requests:
            if not np.isfinite(t) or t < -1e-10 or t > 1+1e-10:
                raise ValueError('Backend boundary request lies outside its canonical edge')
            proposals[int(edge_id)].append(float(np.clip(t,0.,1.)))
        pending = []
        for edge_id in sorted(proposals):
            edge = self.edges[edge_id]
            a,b = [self.vertices[v] for v in edge.endpoints]
            length = np.linalg.norm(b-a)
            accepted = list(edge.parameters)
            for t in sorted(proposals[edge_id]):
                if min(abs(t-u) for u in accepted)*length > self.tolerance:
                    accepted.append(t)
                    pending.append((edge_id,t,a+t*(b-a)))
        if len(self.vertices)+len(pending) > max_vertices:
            raise ValueError(f'Shared-edge vertex budget exceeded ({max_vertices})')
        affected = set()
        for edge_id,t,point in pending:
            edge = self.edges[edge_id]
            vertex = len(self.vertices)
            self.vertices.append(point)
            order = int(np.searchsorted(edge.parameters,t))
            edge.parameters.insert(order,t)
            edge.vertices.insert(order,vertex)
            affected.update(edge.incidents)
        return affected, len(pending)

    def locate(self, edge_id, t):
        edge = self.edges[edge_id]
        k = int(np.argmin(abs(np.asarray(edge.parameters)-t)))
        a,b = [self.vertices[v] for v in edge.endpoints]
        if abs(edge.parameters[k]-t)*np.linalg.norm(b-a) > self.tolerance:
            raise ValueError('Required boundary point is absent from canonical registry')
        return edge.vertices[k]


def seed_corner_shells(registry, faces, *, nominal_spacing, max_vertices):
    """Seed compatible dyadic corner spacing on every incident 3D edge.

    A short feature on one face controls the same corner on all neighbours.
    Equal radial levels prevent incompatible per-face protective shells from
    creating almost coincident optional boundary points. Interior refinement
    remains entirely the maintained backend's responsibility.
    """
    angles = {}
    for face in faces:
        for ring in face.rings:
            for k,vertex in enumerate(ring):
                a=registry.vertices[ring[k-1]]-registry.vertices[vertex]
                b=registry.vertices[ring[(k+1)%len(ring)]]-registry.vertices[vertex]
                angle=float(np.degrees(np.arctan2(np.linalg.norm(np.cross(a,b)),a@b)))
                if angle < 60.:
                    angles[vertex]=min(angle,angles.get(vertex,180.))
    incident=defaultdict(list)
    for ei,edge in enumerate(registry.edges):
        a,b=[registry.vertices[v] for v in edge.endpoints]
        length=float(np.linalg.norm(b-a))
        for vertex in edge.endpoints:incident[vertex].append((ei,length))
    requests=[]
    for vertex,angle in sorted(angles.items()):
        radius=.25*min(nominal_spacing,*[length for _,length in incident[vertex]])
        levels=1
        while levels<6 and angle*2**levels<60.:levels+=1
        for ei,length in incident[vertex]:
            forward=registry.edges[ei].endpoints[0]==vertex
            for level in range(levels):
                t=radius/2**level/length
                requests.append((ei,t if forward else 1-t))
    _,added=registry.reconcile(requests,max_vertices=max_vertices)
    return dict(vertices=added,corners=len(angles),spacing='common dyadic radii from shortest incident feature')
