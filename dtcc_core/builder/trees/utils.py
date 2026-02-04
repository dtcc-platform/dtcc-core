import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from skimage import morphology, segmentation
import rasterio.features
from shapely.geometry import shape, Polygon
from shapely.validation import make_valid
from affine import Affine
from skimage.measure import perimeter


def crown_radius_from_height(
    height: float,
    a: float = 0.22,
    b: float = 0.87,
    min_radius: float = 0.8,
    max_radius: float = 6.0,
) -> float:
    r = a * pow(height, b)
    return float(np.clip(r, min_radius, max_radius))


def detect_local_maxima(
    chm: np.ndarray, mask: np.ndarray, footprint_size: int
) -> tuple[np.ndarray, np.ndarray]:
    """
    Over-detect candidate tree tops using a small window.
    """
    footprint = morphology.footprint_rectangle((footprint_size, footprint_size))
    local_max = chm == ndi.maximum_filter(
        chm, footprint=footprint, mode="constant", cval=np.nanmin(chm)
    )

    local_max &= mask

    coords = np.column_stack(np.nonzero(local_max))
    heights = chm[local_max]

    order = np.argsort(-heights)
    return coords[order], heights[order]


def detect_candidate_peaks(
    chm: np.ndarray, mask: np.ndarray, footprint_size: int
) -> tuple[np.ndarray, np.ndarray]:
    """
    Detect candidate local maxima using a small fixed window.
    """
    footprint = morphology.square(footprint_size)

    local_max = chm == ndi.maximum_filter(
        chm, footprint=footprint, mode="constant", cval=np.nanmin(chm)
    )

    local_max &= mask

    coords = np.column_stack(np.nonzero(local_max))
    heights = chm[local_max]

    # Sort tallest first (critical for greedy suppression)
    order = np.argsort(-heights)
    return coords[order], heights[order]


def suppress_peaks_height_aware(
    coords: np.ndarray,
    heights: np.ndarray,
    chm: np.ndarray,
    pixel_size: float,
    min_radius: float = 0.6,
    max_relative_drop: float = 0.3,
    max_absolute_drop: float = 8.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Variable-window, height-aware non-maximum suppression.

    A peak suppresses nearby peaks if:
    - within its height-scaled radius AND
    - the neighbor is not sufficiently lower
    """
    if len(coords) == 0:
        return coords, heights, np.array([])

    radius = estimate_radius_from_segment(chm, coords, pixel_size)

    tree = cKDTree(coords)
    suppressed = np.zeros(len(coords), dtype=bool)
    suppressed[radius < min_radius] = True
    keep = []

    for i, (pt, h) in enumerate(zip(coords, heights)):
        if suppressed[i]:
            continue

        keep.append(i)

        radius_m = radius[i]
        radius_px = radius_m / pixel_size

        neighbors = tree.query_ball_point(pt, radius_px)

        for j in neighbors:
            if j <= i or suppressed[j]:
                continue

            dh = h - heights[j]

            # Height-relative + absolute constraints
            if dh < max_absolute_drop and dh < max_relative_drop * h:
                suppressed[j] = True

    keep = np.asarray(keep, dtype=int)
    return coords[keep], heights[keep], radius[keep]


def variable_window_nms(
    coords: np.ndarray, heights: np.ndarray, pixel_size: float, radius_fn
) -> tuple[np.ndarray, np.ndarray]:
    """
    Greedy variable-window non-maximum suppression.
    """
    if len(coords) == 0:
        return coords, heights

    tree = cKDTree(coords)
    suppressed = np.zeros(len(coords), dtype=bool)
    keep_indices = []

    for i, (pt, h) in enumerate(zip(coords, heights)):
        if suppressed[i]:
            continue

        keep_indices.append(i)

        radius_m = radius_fn(h)
        radius_px = radius_m / pixel_size

        neighbors = tree.query_ball_point(pt, radius_px)

        for j in neighbors:
            if j > i:
                suppressed[j] = True

    keep_indices = np.asarray(keep_indices, dtype=int)
    return coords[keep_indices], heights[keep_indices]


def segment_tree_crowns(
    chm: np.ndarray, tree_coords: np.ndarray, min_height: float
) -> np.ndarray:
    """
    Watershed-based crown segmentation from detected tree tops.
    """
    markers = np.zeros(chm.shape, dtype=np.int32)

    for i, (r, c) in enumerate(tree_coords, start=1):
        markers[r, c] = i

    mask = np.isfinite(chm) & (chm >= min_height)

    labels = segmentation.watershed(-chm, markers=markers, mask=mask)

    return labels


def estimate_radius_from_segment(
    chm: np.ndarray,
    coords: np.ndarray,
    cell_size: float,
) -> np.ndarray:
    labels = segment_tree_crowns(chm, coords, 0.001)
    max_label = int(labels.max())
    cell_area = cell_size**2
    area = np.zeros(max_label)
    val, count = np.unique(labels, return_counts=True)
    area[val[1:] - 1] = count[1:]
    area *= cell_area
    area *= np.pi / 2

    return np.sqrt(area / np.pi)


def extract_crown_polygons_from_labels(
    labels: np.ndarray,
    affine_transform: Affine,
    simplify_tolerance: float = 0.5,
    min_area: float = 1.0,
) -> list[Polygon]:
    """
    Extract shapely polygons from watershed-labeled raster.

    Parameters
    ----------
    labels : np.ndarray
        Watershed labeled raster (values 1..N for trees, 0 for background)
    affine_transform : Affine
        Georeferencing transform from Raster.georef
    simplify_tolerance : float, optional
        Polygon simplification tolerance in meters (default: 0.5)
    min_area : float, optional
        Minimum polygon area filter in m² (default: 1.0)

    Returns
    -------
    list[Polygon]
        List of shapely Polygons, indexed by label-1 (label 1 → index 0).
        None entries for invalid/filtered polygons.
    """
    # Initialize result list with correct size
    max_label = int(labels.max())
    polygons = [None] * max_label

    # Extract shapes using rasterio with georeferencing
    shapes_gen = rasterio.features.shapes(
        labels.astype(np.int32), mask=(labels > 0), transform=affine_transform
    )

    for geom_dict, value in shapes_gen:
        if value == 0:
            continue

        # Convert GeoJSON dict to Shapely geometry
        polygon = shape(geom_dict)

        # Validate and fix geometry
        if not polygon.is_valid:
            polygon = make_valid(polygon)

        # Handle MultiPolygon - take largest component
        if polygon.geom_type == "MultiPolygon":
            polygon = max(polygon.geoms, key=lambda p: p.area)

        # Skip if not a polygon after validation
        if polygon.geom_type != "Polygon":
            continue

        # Filter by minimum area
        if polygon.area < min_area:
            continue

        # Simplify polygon
        if simplify_tolerance > 0:
            polygon = polygon.simplify(simplify_tolerance, preserve_topology=True)

        # Store in list indexed by label-1
        label_idx = int(value) - 1
        if 0 <= label_idx < max_label:
            polygons[label_idx] = polygon

    return polygons
