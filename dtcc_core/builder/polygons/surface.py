import shapely
import shapely.affinity
from shapely.geometry import Polygon, MultiPolygon
from shapely.validation import make_valid

import numpy as np
from scipy.spatial.transform import Rotation as R
from shapely.ops import unary_union
from ...model.geometry import Surface, MultiSurface, PointCloud

from shapely import minimum_rotated_rectangle
from typing import List, Union

from ..polygons.polygons import remove_slivers
from ..logging import info, warning, error, debug


def clean_multisurface(ms: MultiSurface, simplify=1e-2) -> MultiSurface:
    """Clean a MultiSurface by removing empty and small surfaces.

    Args:
        ms (MultiSurface): The MultiSurface to clean.
        simplify (float): The tolerance for simplifying the polygons.

    Returns:
        MultiSurface: The cleaned MultiSurface.
    """
    cleaned_surfaces = []
    for s in ms.surfaces:
        cleaned_surface = clean_surface(s, simplify)
        if cleaned_surface is not None:
            cleaned_surfaces.append(cleaned_surface)
    failed_surfaces = len(ms.surfaces) - len(cleaned_surfaces)
    # if failed_surfaces > 0:
    #     warning(f"Failed to clean {failed_surfaces} surfaces.")
    cleaned_ms = ms.copy(surfaces=cleaned_surfaces)
    return cleaned_ms


def clean_surface(s: Surface, tol=1e-2, smallest_surface=1e-2) -> Surface:
    """Clean a Surface by removing empty and small holes.

    Args:
        s (Surface): The Surface to clean.
        tol (float): The tolerance for simplifying the polygons.
        smallest_surface (float): The smallest surface to keep.

    Returns:
        Surface: The cleaned Surface.
    """
    trans, inv_trans = _transform_to_planar(s)
    if trans is None:
        return None
    surface_poly = _to_polygon(s, trans)
    poly_z = surface_poly.exterior.coords[0][2]
    if not surface_poly.is_valid:
        surface_poly = make_valid(surface_poly)

    if surface_poly.geom_type in ("MultiPolygon", "GeometryCollection"):
        poly = [
            p
            for p in shapely.get_parts(surface_poly)
            if p.geom_type == "Polygon" and p.area > smallest_surface
        ]
        if len(poly) == 0:
            return None
        # return only largest polygon
        surface_poly = max(poly, key=lambda x: x.area)

    surface_poly = remove_slivers(surface_poly, tol)
    surface_poly.simplify(tol, True)
    if (
        surface_poly.is_empty
        or surface_poly.area < smallest_surface
        or not surface_poly.is_valid
    ):
        # warning(f"Failed to clean surface.")
        return None
    else:
        # Add the z-coordinate back to the polygon since we lose them in the simplification
        sp_3d_ext = [(x, y, poly_z) for x, y in surface_poly.exterior.coords]
        sp_3d_holes = []
        for hole in surface_poly.interiors:
            sp_3d_holes.append([(x, y, poly_z) for x, y in hole.coords])
        surface_poly = Polygon(sp_3d_ext, sp_3d_holes)
        return _to_surface(surface_poly, inv_trans)


def subdivide_surface(s: Surface, longest_edge=2) -> MultiSurface:
    """Subdivide a Surface into smaller surfaces.

    Args:
        s (Surface): The Surface to subdivide.
        longest_edge (float): The maximum length of the edges of the new surfaces.

    Returns:
        MultiSurface: The subdivided Surface
    """
    trans, inv_trans = _transform_to_planar(s)
    surface_poly = _to_polygon(s, trans)
    mrr = minimum_rotated_rectangle(surface_poly)

    mrr_coords = np.array(mrr.exterior.coords)[:-1]
    edge_vector = mrr_coords[1] - mrr_coords[0]
    angle = np.arctan2(edge_vector[1], edge_vector[0])
    cetroid = mrr.centroid
    rotated_mrr = shapely.affinity.rotate(mrr, -angle, origin=cetroid, use_radians=True)

    min_x, min_y, max_x, max_y = rotated_mrr.bounds

    tiles = []
    for x in np.arange(min_x, max_x, longest_edge):
        for y in np.arange(min_y, max_y, longest_edge):
            tile = Polygon(
                [
                    (x, y),
                    (x + longest_edge, y),
                    (x + longest_edge, y + longest_edge),
                    (x, y + longest_edge),
                    (x, y),
                ]
            )
            tiles.append(tile)
    mp_tile = MultiPolygon(tiles)
    mp_tile = shapely.affinity.rotate(mp_tile, angle, origin=cetroid, use_radians=True)
    tiles = []
    for t in mp_tile.geoms:
        tiles.append(t.intersection(surface_poly))
    mp_tile_clipped = MultiPolygon(tiles)


def surface_sample_points(s: Surface, spacing=1.0) -> PointCloud:
    """Sample a Surface to a PointCloud.

    Args:
        s (Surface): The Surface to sample.
        spacing (float): spacing between the points.

    Returns:
        PointCloud: The sampled Point
    """
    trans, inv_trans = _transform_to_planar(s)
    if trans is None:
        error(f"Failed to sample surface.")
    surface_poly = _to_polygon(s, trans)
    mrr = minimum_rotated_rectangle(surface_poly)
    mrr_coords = np.array(mrr.exterior.coords)[:-1]

    # Calculate the angle of rotation
    edge_vector = mrr_coords[1] - mrr_coords[0]
    angle = np.arctan2(edge_vector[1], edge_vector[0])

    # Create a rotation matrix for aligning with the main axis
    rotation_matrix = np.array(
        [[np.cos(-angle), -np.sin(-angle)], [np.sin(-angle), np.cos(-angle)]]
    )
    inverse_rotation_matrix = np.linalg.inv(rotation_matrix)

    rotated_rect_coords = np.dot(mrr_coords, rotation_matrix.T)

    # Determine the bounding box of the rotated minimal rectangle
    min_x, min_y = rotated_rect_coords.min(axis=0)
    max_x, max_y = rotated_rect_coords.max(axis=0)

    # Generate a grid of points within the bounding box of the rotated minimal rectangle
    x_range = np.arange(min_x, max_x, spacing)
    y_range = np.arange(min_y, max_y, spacing)
    grid_points = [np.array([x, y]) for x in x_range for y in y_range]

    # Rotate the points back to the original coordinate system
    aligned_points = np.array(
        [np.dot(point, inverse_rotation_matrix.T) for point in grid_points]
    )
    return PointCloud(points=aligned_points)


def union_surfaces(ms: MultiSurface | List[Surface]) -> Surface:
    """Union a list of surfaces into a single Surface. Assumes that all surfaces are co-planar
    and connected. If this is not the case, the result is undefined.

    Args:
        ms (MultiSurface|List[Surface]): The MultiSurface or list of Surfaces to union.
    """

    if isinstance(ms, MultiSurface):
        surfaces = ms.surfaces
    else:
        surfaces = ms

    if len(surfaces) == 0:
        return Surface()
    if len(surfaces) == 1:
        return surfaces[0]

    transform, transform_inv = _transform_to_planar(surfaces[0])
    surface_poly = [_to_polygon(s, transform) for s in surfaces]
    union_poly = unary_union(surface_poly)

    if union_poly.geom_type == "MultiPolygon":
        union_poly = union_poly.convex_hull

    # union_poly = remove_slivers(union_poly, tol=1e-3)
    if not union_poly.is_valid:
        union_poly = make_valid(union_poly)
    union_poly = union_poly.simplify(1e-6, preserve_topology=True)
    merged_surface = _to_surface(union_poly, transform_inv)
    return merged_surface


def _to_polygon(s: Surface, transform) -> Polygon:
    """Convert a Surface to a Polygon."""
    transformed_vertices = np.dot(s.vertices, transform[:3, :3].T) + transform[:3, 3]
    transformed_holes = [
        np.dot(h, transform[:3, :3].T) + transform[:3, 3] for h in s.holes
    ]
    polygon = Polygon(transformed_vertices, transformed_holes)
    return polygon


def _to_surface(p: Polygon, inv_transform) -> Surface:
    """Convert a Polygon to a Surface."""
    transformed_vertices = (
        np.dot(np.array(p.exterior.coords), inv_transform[:3, :3].T)
        + inv_transform[:3, 3]
    )[:-1]

    transformed_holes = [
        (np.dot(np.array(h.coords), inv_transform[:3, :3].T) + inv_transform[:3, 3])[
            :-1
        ]
        for h in p.interiors
    ]
    s = Surface(vertices=transformed_vertices, holes=transformed_holes)
    return s


def _transform_to_planar(s: Surface) -> (np.ndarray, np.ndarray):
    """Transform a Surface to a planar Surface."""
    z_normal = np.array([0, 0, 1])
    normal = s.calculate_normal()
    if np.allclose(normal, 0):
        # warning(f"Surface has zero normal.")
        return None, None
    centroid = s.centroid

    rotation = R.align_vectors(z_normal, normal)[0].as_matrix()
    # Calculate the translation
    translation = -centroid

    transform = np.eye(4)
    transform[:3, :3] = rotation
    transform[:3, 3] = translation

    # Compute the inverse of the transformation matrix
    transform_inv = np.linalg.inv(transform)

    return transform, transform_inv
