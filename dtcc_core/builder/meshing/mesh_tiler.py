import numpy as np
from typing import Tuple, Dict, List, Optional, Set
from shapely.geometry import Polygon, LineString, Point, box
from shapely.ops import triangulate as shapely_triangulate
from shapely.strtree import STRtree
from dataclasses import dataclass


from dtcc_core.model import Mesh, Bounds


@dataclass
class TriangleInfo:
    """Preprocessed triangle information for faster clipping."""

    index: int
    vertices: np.ndarray
    polygon_2d: Optional[Polygon]
    is_vertical: bool
    normal_up: bool
    bbox: Tuple[float, float, float, float]


class SurfaceMeshClipper:
    """
    Optimized clipper for multiple bounding box operations on the same mesh.
    Uses spatial indexing and preprocessing for better performance.
    """

    def __init__(self, mesh: Mesh, chunk_size: int = 1000):
        """
        Initialize the clipper with a mesh and preprocess it.

        Args:
            mesh: Surface mesh to be clipped
        """
        self.mesh = mesh
        self.chunk_size = chunk_size

        # Preprocess all triangles
        self.triangle_infos = self._preprocess_triangles()

        # Build spatial index
        self._build_spatial_index()

        # Precompute triangle adjacency for connected component analysis
        self._build_adjacency()

    def _preprocess_triangles(self) -> List[TriangleInfo]:
        """Preprocess all triangles to avoid repeated calculations."""
        infos = []
        vertices = self.mesh.vertices

        # Vectorized normal computation for all triangles at once
        v0 = vertices[self.mesh.faces[:, 0]]
        v1 = vertices[self.mesh.faces[:, 1]]
        v2 = vertices[self.mesh.faces[:, 2]]

        # Compute all normals vectorized
        edges1 = v1 - v0
        edges2 = v2 - v0
        normals = np.cross(edges1, edges2)
        norms = np.linalg.norm(normals, axis=1)

        # Compute 2D areas vectorized
        areas_2d = 0.5 * np.abs(
            edges1[:, 0] * edges2[:, 1] - edges1[:, 1] * edges2[:, 0]
        )

        for i, face in enumerate(self.mesh.faces):
            tri_verts = vertices[face]

            # Check if vertical
            is_vertical = areas_2d[i] < 1e-10

            # Get normal direction
            normal_up = norms[i] > 0 and normals[i, 2] > 0

            # Create 2D projection if not vertical
            if not is_vertical:
                poly_2d = Polygon(tri_verts[:, :2])
            else:
                poly_2d = None

            # Compute bounding box
            bbox = (
                tri_verts[:, 0].min(),
                tri_verts[:, 1].min(),
                tri_verts[:, 0].max(),
                tri_verts[:, 1].max(),
            )

            infos.append(
                TriangleInfo(
                    index=i,
                    vertices=tri_verts,
                    polygon_2d=poly_2d,
                    is_vertical=is_vertical,
                    normal_up=normal_up,
                    bbox=bbox,
                )
            )

        return infos

    def _build_spatial_index(self):
        """Build an R-tree spatial index for fast triangle lookup."""
        # Create geometry objects for spatial index
        geometries = []
        for info in self.triangle_infos:
            if info.is_vertical:
                # Use line for vertical triangles
                line = LineString(info.vertices[:, :3])
                geometries.append(line)
            else:
                # Use polygon for normal triangles
                geometries.append(info.polygon_2d)

        # Build STRtree (R-tree variant)
        self.spatial_index = STRtree(geometries)
        # self.geometry_to_info = {
        #     info.index: info for geom, info in zip(geometries, self.triangle_infos)
        # }

    def _build_adjacency(self):
        """Build triangle adjacency for connected components."""
        # Create edge to triangle mapping
        edge_to_triangles = {}

        for i, face in enumerate(self.mesh.faces):
            # Sort vertex indices for each edge
            edges = [
                tuple(sorted([face[0], face[1]])),
                tuple(sorted([face[1], face[2]])),
                tuple(sorted([face[2], face[0]])),
            ]

            for edge in edges:
                if edge not in edge_to_triangles:
                    edge_to_triangles[edge] = []
                edge_to_triangles[edge].append(i)

        # Build adjacency list
        self.adjacency = {i: set() for i in range(len(self.mesh.faces))}
        for triangles in edge_to_triangles.values():
            for i in range(len(triangles)):
                for j in range(i + 1, len(triangles)):
                    self.adjacency[triangles[i]].add(triangles[j])
                    self.adjacency[triangles[j]].add(triangles[i])

    def clip_to_bounds(self, bounds: Bounds) -> Mesh:
        """
        Clip the mesh to a single bounding box.

        Args:
            bbox: Tuple of (xmin, ymin, xmax, ymax)

        Returns:
            New Mesh clipped to bbox
        """

        bbox = bounds.tuple
        clip_box = box(*bbox)

        # Query spatial index for potentially intersecting triangles
        potential_triangles = self._query_triangles(clip_box)

        if not potential_triangles:
            return Mesh(vertices=np.array([]), faces=np.array([], dtype=np.int32))

        # Process only relevant triangles
        clipper = SingleBBoxClipper(bbox)

        for info in potential_triangles:
            if info.is_vertical:
                clipper.process_vertical_triangle(info)
            else:
                clipper.process_normal_triangle(info)

        return clipper.get_result()

    def clip_to_multiple_bboxes(
        self, bboxes: List[Tuple[float, float, float, float]]
    ) -> List[Mesh]:
        """
        Clip the mesh to multiple bounding boxes efficiently.

        Args:
            bboxes: List of bounding boxes

        Returns:
            List of clipped meshes
        """
        results = []

        # Process in batches to optimize cache usage
        for bbox in bboxes:
            results.append(self.clip_to_bounds(bbox))

        return results

    def clip_to_grid(
        self, grid_bounds: Tuple[float, float, float, float], nx: int, ny: int
    ) -> Dict[Tuple[int, int], Mesh]:
        """
        Clip mesh to a regular grid of bounding boxes.

        Args:
            grid_bounds: Overall bounds (xmin, ymin, xmax, ymax)
            nx: Number of grid cells in X direction
            ny: Number of grid cells in Y direction

        Returns:
            Dictionary mapping (i, j) grid indices to clipped meshes
        """
        xmin, ymin, xmax, ymax = grid_bounds
        dx = (xmax - xmin) / nx
        dy = (ymax - ymin) / ny

        results = {}

        for i in range(nx):
            for j in range(ny):
                cell_bbox = Bounds(
                    xmin=xmin + i * dx,
                    ymin=ymin + j * dy,
                    xmax=xmin + (i + 1) * dx,
                    ymax=ymin + (j + 1) * dy,
                )

                clipped = self.clip_to_bounds(cell_bbox)
                if len(clipped.vertices) > 0:
                    results[(i, j)] = clipped

        return results

    def _query_triangles(self, clip_box: Polygon) -> List[TriangleInfo]:
        """Query spatial index for triangles that might intersect the clip box."""
        # Query the spatial index
        potential_geoms = self.spatial_index.query(clip_box)

        # Map back to triangle infos
        potential_triangles = []
        for geom_idx in potential_geoms:
            potential_triangles.append(self.triangle_infos[geom_idx])

        return potential_triangles


class SingleBBoxClipper:
    """Handles clipping for a single bounding box."""

    def __init__(self, bbox: Tuple[float, float, float, float]):
        self.bbox = bbox
        self.clip_box = box(*bbox)
        self.vertices = []
        self.faces = []
        self.vertex_map = {}

    def process_normal_triangle(self, info: TriangleInfo):
        """Process a non-vertical triangle."""
        # Clip the triangle
        clipped = info.polygon_2d.intersection(self.clip_box)

        if clipped.is_empty:
            return

        self._process_clipped_geometry(clipped, info.vertices, info.normal_up)

    def process_vertical_triangle(self, info: TriangleInfo):
        """Process a vertical triangle."""
        # For vertical triangles, we need to handle the full 3D geometry
        # Create a 2D line for horizontal clipping
        line_2d = LineString(info.vertices[:, :2])
        
        # Check if the triangle is completely inside the bounding box
        if self.clip_box.contains(line_2d):
            # Triangle is completely inside - add it as-is without modification
            self._add_triangle(info.vertices, info.normal_up)
            return
            
        # Clip the line against the bounding box
        clipped = line_2d.intersection(self.clip_box)

        if clipped.is_empty or clipped.geom_type != "LineString":
            return

        coords = list(clipped.coords)
        if len(coords) < 2:
            return

        # Check if the clipped line is the same as the original
        # (This handles edge cases where intersection returns the original but contains() fails)
        original_coords = list(line_2d.coords)
        if len(coords) == len(original_coords):
            # Compare coordinates with small tolerance
            coords_match = True
            for c1, c2 in zip(coords, original_coords):
                if abs(c1[0] - c2[0]) > 1e-10 or abs(c1[1] - c2[1]) > 1e-10:
                    coords_match = False
                    break
            
            if coords_match:
                # No actual clipping occurred - add original triangle
                self._add_triangle(info.vertices, info.normal_up)
                return

        # Triangle was actually clipped - create a quad for the clipped portion
        # Get the actual Z range from the original vertical triangle
        z_min = info.vertices[:, 2].min()
        z_max = info.vertices[:, 2].max()

        self._create_vertical_quad(coords[0], coords[-1], z_min, z_max, info.normal_up)

    def get_result(self) -> Mesh:
        """Get the final clipped mesh."""
        if not self.vertices:
            return Mesh(vertices=np.array([]), faces=np.array([], dtype=np.int32))

        vertices = np.array(self.vertices)
        faces = np.array(self.faces, dtype=np.int32)

        # Ensure upward normals
        faces = self._ensure_upward_normals(vertices, faces)

        return Mesh(vertices=vertices, faces=faces)

    def _process_clipped_geometry(self, geom, tri: np.ndarray, original_up: bool):
        """Process clipped geometry."""
        if geom.geom_type == "Polygon":
            self._triangulate_polygon(geom, tri, original_up)
        elif geom.geom_type == "MultiPolygon":
            for poly in geom.geoms:
                self._triangulate_polygon(poly, tri, original_up)

    def _triangulate_polygon(self, poly: Polygon, tri: np.ndarray, original_up: bool):
        """Triangulate a clipped polygon."""
        if not poly.is_valid or poly.is_empty:
            return

        vertices_2d = np.array(poly.exterior.coords[:-1])

        if len(vertices_2d) < 3:
            return

        # Interpolate Z values using vectorized barycentric coordinates
        z_values = self._interpolate_z_batch(vertices_2d, tri)
        vertices_3d = np.column_stack([vertices_2d, z_values])

        if len(vertices_3d) == 3:
            self._add_triangle(vertices_3d, original_up)
        else:
            for i in range(1, len(vertices_3d) - 1):
                tri = np.array([vertices_3d[0], vertices_3d[i], vertices_3d[i + 1]])
                self._add_triangle(tri, original_up)

    def _triangulate_complex(
        self, vertices_3d: np.ndarray, poly: Polygon, original_up: bool
    ):
        """Triangulate complex polygon."""
        try:
            triangles = list(shapely_triangulate(poly))
            point_map = {tuple(v[:2]): v for v in vertices_3d}

            for tri in triangles:
                if hasattr(tri, "exterior"):
                    tri_2d = np.array(tri.exterior.coords[:-1])
                    if len(tri_2d) == 3:
                        tri_3d = np.array(
                            [point_map.get(tuple(pt), vertices_3d[0]) for pt in tri_2d]
                        )
                        self._add_triangle(tri_3d, original_up)
        except:
            # Fallback to fan triangulation
            for i in range(1, len(vertices_3d) - 1):
                tri = np.array([vertices_3d[0], vertices_3d[i], vertices_3d[i + 1]])
                self._add_triangle(tri, original_up)

    def _get_z_for_vertical_point(
        self, point_2d: np.ndarray, vertices: np.ndarray
    ) -> float:
        """Get appropriate Z value for a point on a vertical triangle."""
        # For a truly vertical triangle, all vertices share similar X,Y coordinates
        # Return the Z value of the closest vertex
        dists = np.linalg.norm(vertices[:, :2] - point_2d, axis=1)
        closest_idx = np.argmin(dists)
        return vertices[closest_idx, 2]

    def _create_vertical_quad(
        self, p0_2d: Tuple, p1_2d: Tuple, z_min: float, z_max: float, original_up: bool
    ):
        """Create vertical quad for vertical triangle.

        Args:
            p0_2d, p1_2d: The 2D endpoints of the clipped line
            z_min, z_max: The Z range for the vertical surface
            original_up: Whether the original triangle normal pointed up
        """
        dx = p1_2d[0] - p0_2d[0]
        dy = p1_2d[1] - p0_2d[1]
        length = np.sqrt(dx * dx + dy * dy)

        if length < 1e-10:
            # For a point-like vertical triangle, create a small vertical quad
            # Use a small offset in both X and Y directions
            offset = 1e-6
            v0 = np.array([p0_2d[0] - offset, p0_2d[1] - offset, z_min])
            v1 = np.array([p0_2d[0] + offset, p0_2d[1] - offset, z_min])
            v2 = np.array([p0_2d[0] - offset, p0_2d[1] + offset, z_max])
            v3 = np.array([p0_2d[0] + offset, p0_2d[1] + offset, z_max])
        else:
            # Normal case: create a thin quad perpendicular to the line direction
            offset = 1e-6
            perp_x = -dy / length * offset
            perp_y = dx / length * offset

            v0 = np.array([p0_2d[0] - perp_x, p0_2d[1] - perp_y, z_min])
            v1 = np.array([p0_2d[0] + perp_x, p0_2d[1] + perp_y, z_min])
            v2 = np.array([p1_2d[0] - perp_x, p1_2d[1] - perp_y, z_max])
            v3 = np.array([p1_2d[0] + perp_x, p1_2d[1] + perp_y, z_max])

        # Create two triangles to form the quad
        # Ensure correct winding order for upward-facing normals
        self._add_triangle(np.array([v0, v1, v2]), original_up)
        self._add_triangle(np.array([v1, v3, v2]), original_up)

    def _add_triangle(self, tri: np.ndarray, should_face_up: bool):
        """Add triangle with correct normal."""
        indices = [self._add_vertex(v) for v in tri]

        normal = self._compute_normal(tri)
        is_facing_up = normal[2] > 0

        if should_face_up != is_facing_up:
            indices = [indices[0], indices[2], indices[1]]

        self.faces.append(indices)

    def _add_vertex(self, vertex: np.ndarray) -> int:
        """Add vertex, avoiding duplicates."""
        key = tuple(np.round(vertex, 10))

        if key not in self.vertex_map:
            self.vertex_map[key] = len(self.vertices)
            self.vertices.append(vertex)

        return self.vertex_map[key]

    def _interpolate_z_batch(
        self, points_2d: np.ndarray, tri: np.ndarray
    ) -> np.ndarray:
        """Vectorized Z interpolation using barycentric coordinates."""
        v0, v1, v2 = tri[0], tri[1], tri[2]

        # Precompute for barycentric
        v01 = v1[:2] - v0[:2]
        v02 = v2[:2] - v0[:2]
        denom = v01[0] * v02[1] - v01[1] * v02[0]

        if abs(denom) < 1e-10:
            # Degenerate - use edge interpolation
            return np.array([self._interpolate_z_on_line(p, tri) for p in points_2d])

        # Vectorized barycentric computation
        v0p = points_2d - v0[:2]
        u = (v0p[:, 0] * v02[1] - v0p[:, 1] * v02[0]) / denom
        v = (v01[0] * v0p[:, 1] - v01[1] * v0p[:, 0]) / denom
        w = 1 - u - v

        # Check validity and compute Z
        z_values = np.where(
            (u >= -1e-10) & (v >= -1e-10) & (w >= -1e-10),
            w * v0[2] + u * v1[2] + v * v2[2],
            [self._interpolate_z_on_line(p, tri) for p in points_2d],
        )

        return z_values

    def _interpolate_z_on_line(self, point_2d: Tuple, tri: np.ndarray) -> float:
        """Interpolate Z using edge projection."""
        point = Point(point_2d)
        best_z = tri[0, 2]
        min_dist = float("inf")

        for i in range(3):
            j = (i + 1) % 3
            edge = LineString([tri[i, :2], tri[j, :2]])
            dist = point.distance(edge)

            if dist < min_dist:
                min_dist = dist
                if edge.length > 0:
                    t = np.clip(edge.project(point) / edge.length, 0, 1)
                    best_z = (1 - t) * tri[i, 2] + t * tri[j, 2]
                else:
                    best_z = (tri[i, 2] + tri[j, 2]) / 2

        return best_z

    def _compute_normal(self, tri: np.ndarray) -> np.ndarray:
        """Compute triangle normal."""
        v1 = tri[1] - tri[0]
        v2 = tri[2] - tri[0]
        normal = np.cross(v1, v2)
        norm = np.linalg.norm(normal)
        return normal / norm if norm > 0 else normal

    def _ensure_upward_normals(
        self, vertices: np.ndarray, faces: np.ndarray
    ) -> np.ndarray:
        """Ensure all normals point upward."""
        # Vectorized normal check
        v0 = vertices[faces[:, 0]]
        v1 = vertices[faces[:, 1]]
        v2 = vertices[faces[:, 2]]

        normals = np.cross(v1 - v0, v2 - v0)
        should_flip = normals[:, 2] < 0

        # Flip faces where needed
        faces[should_flip] = faces[should_flip][:, [0, 2, 1]]

        return faces
