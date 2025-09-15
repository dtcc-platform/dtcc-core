from dtcc_core.model import Mesh, Bounds
from typing import Tuple, Union, List
from shapely.geometry import Polygon, Point, LineString, box
from shapely.ops import triangulate
import numpy as np


def tile_surface_mesh(
    mesh: Mesh,
    tile_size: Union[float, Tuple[float, float]] = 100.0,
) -> list[Mesh]:
    """
    Tile a surface mesh into smaller meshes of a specified size.

    Args:
        mesh (Mesh): The input surface mesh to be tiled.
        tile_size (float|tuple): The size of each tile (in the same units as the mesh coordinates). if a single float
        is provided, it will be used for both width and height.
        bounds (Bounds): Optional bounds to limit the tiling area.

    Returns:
        list[Mesh]: A list of tiled mesh objects.
    """
    from dtcc_core.builder.meshing import SurfaceMeshClipper

    surface_bounds = mesh.bounds
    surface_tile_bounds = surface_bounds.tiles(tile_size)
    tiled_meshes = []
    clipper = SurfaceMeshClipper(mesh)

    for idx, tile in enumerate(surface_tile_bounds):
        print(f"Tiling mesh: {idx+1}/{len(surface_tile_bounds)}")
        tiled_mesh = clipper.clip_to_bounds(tile)
        tiled_meshes.append(tiled_mesh)
    return tiled_meshes


def _barycentric_coordinates(
    x, y, v0, v1, v2
) -> Union[Tuple[float, float, float], None]:
    """
    Calculate barycentric coordinates of point (x,y) with respect to triangle v0,v1,v2.
    Returns None if the triangle is degenerate.
    """
    denom = (v1[1] - v2[1]) * (v0[0] - v2[0]) + (v2[0] - v1[0]) * (v0[1] - v2[1])

    if abs(denom) < 1e-10:
        return None

    a = ((v1[1] - v2[1]) * (x - v2[0]) + (v2[0] - v1[0]) * (y - v2[1])) / denom
    b = ((v2[1] - v0[1]) * (x - v2[0]) + (v0[0] - v2[0]) * (y - v2[1])) / denom
    c = 1 - a - b

    # Check if point is inside triangle (with some tolerance for boundary points)
    eps = -1e-10
    if a >= eps and b >= eps and c >= eps:
        # Clamp to valid range [0, 1]
        a = max(0, min(1, a))
        b = max(0, min(1, b))
        c = max(0, min(1, c))
        # Normalize to ensure they sum to 1
        total = a + b + c
        return (a / total, b / total, c / total)

    return None


def _interpolate_z_on_edge(pt, p0, p1, p2) -> float:
    point = Point(pt[0], pt[1])
    best_z = p0[2]
    min_dist = float("inf")
    for pa, pb in [(p0, p1), (p1, p2), (p2, p0)]:
        edge = LineString([(pa[0], pa[1]), (pb[0], pb[1])])
        if edge.length < 1e-10:
            continue
        dist = point.distance(edge)
        if dist < min_dist:
            min_dist = dist
            # Project point onto edge to find closest point
            projected = edge.interpolate(edge.project(point))
            # Interpolate z value based on position along the edge

            t = edge.project(projected) / edge.length
            best_z = (1 - t) * pa[2] + t * pb[2]
    return best_z


def _interpolate_z_values(pts, p0, p1, p2) -> List[float]:
    z_values = []
    for pt in pts:
        bary = _barycentric_coordinates(pt[0], pt[1], p0, p1, p2)
        if bary is not None:
            a, b, c = bary
            z = a * p0[2] + b * p1[2] + c * p2[2]
            z_values.append(z)
        else:
            z = _interpolate_z_on_edge(pt, p0, p1, p2)
            z_values.append(z)
    return z_values


def _fan_triangulate(
    vertices: List[Tuple[float, float]]
) -> List[List[Tuple[float, float]]]:
    """
    Triangulate a polygon using a fan triangulation method.
    Assumes the polygon is simple and has at least 3 vertices.
    """
    triangles = []
    if len(vertices) < 3:
        return triangles
    for i in range(1, len(vertices) - 1):
        triangles.append([vertices[0], vertices[i], vertices[i + 1]])
    return triangles


def _compute_normal(v0, v1, v2) -> np.ndarray:
    u = v1 - v0
    v = v2 - v0
    n = np.cross(u, v)
    norm = np.linalg.norm(n)
    if norm < 1e-10:
        return np.array([0.0, 0.0, 0.0])
    return n / norm


def clip_mesh_to_bounds(mesh: Mesh, bounds: Bounds) -> Mesh:
    xmin, ymin, xmax, ymax = bounds.tuple
    clip_box = box(xmin, ymin, xmax, ymax)
    vertices = mesh.vertices

    new_vertices = []
    new_faces = []
    vertex_map = {}  # Map from old vertex index to new vertex index

    def add_vertex(x, y, z):
        # avoid floating point precision issues
        key = (round(x, 6), round(y, 6), round(z, 6))
        if key not in vertex_map:
            vertex_map[key] = len(new_vertices)
            new_vertices.append([x, y, z])
        return vertex_map[key]

    # Create a mask for vertices inside the bounds
    inside_mask = (
        (vertices[:, 0] >= xmin)
        & (vertices[:, 0] <= xmax)
        & (vertices[:, 1] >= ymin)
        & (vertices[:, 1] <= ymax)
    )

    for face in mesh.faces:
        v0, v1, v2 = face[:3]
        p0, p1, p2 = vertices[v0], vertices[v1], vertices[v2]
        base_normal = _compute_normal(p0, p1, p2)
        triangle_2d = Polygon([(p0[0], p0[1]), (p1[0], p1[1]), (p2[0], p2[1])])
        if not triangle_2d.intersects(clip_box):
            continue
        if clip_box.contains(triangle_2d):
            idx0 = add_vertex(p0[0], p0[1], p0[2])
            idx1 = add_vertex(p1[0], p1[1], p1[2])
            idx2 = add_vertex(p2[0], p2[1], p2[2])
            new_faces.append([idx0, idx1, idx2])
            continue

        clipped = triangle_2d.intersection(clip_box)
        if clipped.is_empty or clipped.area == 0:
            continue
        if clipped.geom_type == "Polygon":
            clipped_vertices = list(
                clipped.exterior.coords[:-1]
            )  # Remove duplicate last point
            if len(clipped_vertices) >= 3:
                z_values = _interpolate_z_values(clipped_vertices, p0, p1, p2)
                if len(clipped_vertices) == 3:  # triangle
                    indices = []
                    for (x, y), z in zip(clipped_vertices, z_values):
                        idx = add_vertex(x, y, z)
                        indices.append(idx)
                    if len(indices) != 3:
                        pass
                    clipped_triangle = np.array(new_vertices)[indices]
                    clipped_normal = _compute_normal(
                        clipped_triangle[0],
                        clipped_triangle[1],
                        clipped_triangle[2],
                    )
                    if np.dot(base_normal, clipped_normal) < 0:
                        indices = [indices[0], indices[2], indices[1]]

                    new_faces.append(indices)
                else:  # polygon with > 3 vertices - triangulate
                    triangles = _fan_triangulate(clipped_vertices)
                    z_map = {(x, y): z for (x, y), z in zip(clipped_vertices, z_values)}
                    for tri in triangles:
                        indices = []
                        for x, y in tri:
                            z = z_map[(x, y)]
                            idx = add_vertex(x, y, z)
                            indices.append(idx)
                        if len(indices) != 3:
                            pass
                        clipped_triangle = np.array(new_vertices)[indices]
                        clipped_normal = _compute_normal(
                            clipped_triangle[0],
                            clipped_triangle[1],
                            clipped_triangle[2],
                        )
                        if np.dot(base_normal, clipped_normal) < 0:
                            indices = [indices[0], indices[2], indices[1]]
                        new_faces.append(indices)

    return Mesh(
        vertices=np.array(new_vertices, dtype=np.float64),
        faces=np.array(new_faces, dtype=np.int64),
    )
