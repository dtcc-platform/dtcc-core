
import numpy as np
from ...model import Mesh


def get_east_boundary_vertices(vertices, xmax=None, tol=1e-3):
    """
    Return indices on the EAST boundary (x ≈ xmax), sorted from SOUTH→NORTH.
    
    Sorting:
        - 2D: by y ascending
        - 3D: by (y, then z) ascending

    Building the EAST wall polygon with OUTWARD normal (+x):
        Let g = indices on ground (z = zmin), r = corresponding indices on roof (z = zmax).
        Use CCW order as seen from +x (looking toward -x):
            [ g (south→north), r (north→south) ]
        i.e., append the roof sequence in REVERSE to close the quad loop.
        This CCW ordering (right-hand rule) makes the facet normal point to +x (outward).

    Parameters
    ----------
    vertices : (N,2) or (N,3) array-like
    xmax : float or None
        If None, computed from vertices[:,0].max()
    tol : float
        Tolerance for boundary membership.

    Returns
    -------
    idx : (K,) ndarray of int
        Vertex indices, sorted south→north.
    """
    V = np.asarray(vertices)
    if xmax is None:
        xmax = V[:, 0].max()
    mask = V[:, 0] >= (xmax - tol)
    idx = np.flatnonzero(mask)
    if idx.size:
        if V.shape[1] >= 3:
            order = np.lexsort((V[idx, 2], V[idx, 1]))  # by y, then z
        else:
            order = np.argsort(V[idx, 1])               # by y
        idx = idx[order]
    return idx


def get_west_boundary_vertices(vertices, ymin=None, ymax=None, xmin=None, tol=1e-3):
    """
    Return indices on the WEST boundary (x ≈ xmin), sorted from NORTH→SOUTH.

    Why NORTH→SOUTH?
        For the WEST wall (outward normal -x), a convenient CCW ordering (seen from -x)
        is: [ g (north→south), r (south→north) ]. This yields outward normal (-x).

    Sorting:
        - 2D: by y descending
        - 3D: by (y descending, then z ascending)

    Parameters
    ----------
    vertices : (N,2) or (N,3) array-like
    xmin : float or None
        If None, computed from vertices[:,0].min()
    tol : float
        Tolerance for boundary membership.

    Returns
    -------
    idx : (K,) ndarray of int
        Vertex indices, sorted north→south.
    """
    V = np.asarray(vertices)
    if xmin is None:
        xmin = V[:, 0].min()
    mask = V[:, 0] <= (xmin + tol)
    idx = np.flatnonzero(mask)
    if idx.size:
        if V.shape[1] >= 3:
            # sort by y DESC, then z ASC
            order = np.lexsort((V[idx, 2], -V[idx, 1]))
        else:
            order = np.argsort(-V[idx, 1])  # y descending
        idx = idx[order]
    return idx


def get_south_boundary_vertices(vertices, ymin=None, tol=1e-3):
    """
    Return indices on the SOUTH boundary (y ≈ ymin), sorted from WEST→EAST.

    Building the SOUTH wall polygon with OUTWARD normal (-y):
        View from -y (outside). CCW ordering:
            [ g (west→east), r (east→west) ]
        i.e., append roof indices in REVERSE to close the quad; normal points to -y.

    Sorting:
        - 2D: by x ascending
        - 3D: by (x ascending, then z ascending)

    Parameters
    ----------
    vertices : (N,2) or (N,3) array-like
    ymin : float or None
        If None, computed from vertices[:,1].min()
    tol : float
        Tolerance for boundary membership.

    Returns
    -------
    idx : (K,) ndarray of int
        Vertex indices, sorted west→east.
    """
    V = np.asarray(vertices)
    if ymin is None:
        ymin = V[:, 1].min()
    mask = V[:, 1] <= (ymin + tol)
    idx = np.flatnonzero(mask)
    if idx.size:
        if V.shape[1] >= 3:
            order = np.lexsort((V[idx, 2], V[idx, 0]))  # by x, then z
        else:
            order = np.argsort(V[idx, 0])               # by x
        idx = idx[order]
    return idx


def get_north_boundary_vertices(vertices, ymax=None, tol=1e-3):
    """
    Return indices on the NORTH boundary (y ≈ ymax), sorted from EAST→WEST.

    Why EAST→WEST?
        For the NORTH wall (outward normal +y), viewing from +y, a CCW ordering is:
            [ g (east→west), r (west→east) ]
        Reversing the roof segment closes the quad and yields outward normal (+y).

    Sorting:
        - 2D: by x descending
        - 3D: by (x descending, then z ascending)

    Parameters
    ----------
    vertices : (N,2) or (N,3) array-like
    ymax : float or None
        If None, computed from vertices[:,1].max()
    tol : float
        Tolerance for boundary membership.

    Returns
    -------
    idx : (K,) ndarray of int
        Vertex indices, sorted east→west.
    """
    V = np.asarray(vertices)
    if ymax is None:
        ymax = V[:, 1].max()
    mask = V[:, 1] >= (ymax - tol)
    idx = np.flatnonzero(mask)
    if idx.size:
        if V.shape[1] >= 3:
            # sort by x DESC, then z ASC
            order = np.lexsort((V[idx, 2], -V[idx, 0]))
        else:
            order = np.argsort(-V[idx, 0])  # x descending
        idx = idx[order]
    return idx



import numpy as np

def compute_boundary_facets(mesh: Mesh, top_height=100.0, tol=1e-3):
    """
    Build five PLC facet polygons (south, east, north, west, top) for a rectangular box
    around the mesh domain, using 4 new top-corner points. Polygons are wound CCW as seen
    from *outside* so their normals point outward (right-hand rule), i.e.:
      - South  -> outward -y
      - East   -> outward +x
      - North  -> outward +y
      - West   -> outward -x
      - Top    -> outward +z

    Returns
    -------
    vertices : (N+4,3) float64
        Original vertices with the 4 added top points appended at the end.
    facets : dict[str, list[int]]
        Indices of the five polygons with outward normals:
        {'south': [...], 'east': [...], 'north': [...], 'west': [...], 'top': [...]}
    """
    V = np.asarray(mesh.vertices, dtype=float)
    xmin, ymin, zmin = np.min(V, axis=0)
    xmax, ymax, zmax = np.max(V, axis=0)

    # 1) Height check / adjust
    domain_h = float(zmax - zmin)
    height = float(top_height)
    if height <= domain_h:
        # Make it clearly taller than the domain to avoid intersecting the terrain/buildings
        height = 1.5 * domain_h if domain_h > 0 else max(1.0, top_height)

    # 2) Grab boundary indices (sorted for correct ground-edge order)
    east_idx  = get_east_boundary_vertices (V, xmax=xmax, tol=tol)   # south->north
    west_idx  = get_west_boundary_vertices (V, xmin=xmin, tol=tol)   # north->south
    south_idx = get_south_boundary_vertices(V, ymin=ymin, tol=tol)   # west->east
    north_idx = get_north_boundary_vertices(V, ymax=ymax, tol=tol)   # east->west

    # 3) Create 4 top-corner points (appended to the vertex array)
    z_top = zmin + height
    top_points = np.array([
        [xmin, ymin, z_top],  # t_sw: south-west  (index = N + 0)
        [xmin, ymax, z_top],  # t_nw: north-west  (index = N + 1)
        [xmax, ymin, z_top],  # t_se: south-east  (index = N + 2)
        [xmax, ymax, z_top],  # t_ne: north-east  (index = N + 3)
    ], dtype=float)

    N0 = V.shape[0]
    t_sw, t_nw, t_se, t_ne = N0 + 0, N0 + 1, N0 + 2, N0 + 3
    V_out = np.vstack([V, top_points])

    # 4) Build five polygons with CCW winding as seen from outside

    # SOUTH wall (y = ymin, outward -y): CCW seen from -y
    #   ground: west->east  then top: east->west (reverse)
    south_poly = list(south_idx) + [t_se, t_sw]

    # EAST wall (x = xmax, outward +x): CCW seen from +x
    #   ground: south->north then top: north->south (reverse)
    east_poly  = list(east_idx)  + [t_ne, t_se]

    # NORTH wall (y = ymax, outward +y): CCW seen from +y
    #   ground: east->west then top: west->east
    north_poly = list(north_idx) + [t_nw, t_ne]

    # WEST wall (x = xmin, outward -x): CCW seen from -x
    #   ground: north->south then top: south->north
    west_poly  = list(west_idx)  + [t_sw, t_nw]

    # TOP cap (z = z_top, outward +z): CCW in XY as seen from above (+z)
    #   CCW rectangle: (xmin,ymin)->(xmax,ymin)->(xmax,ymax)->(xmin,ymax)
    top_poly   = [t_sw, t_se, t_ne, t_nw]

    facets = {
        "south": np.array(south_poly),
        "east":  np.array(east_poly),
        "north": np.array(north_poly),
        "west":  np.array(west_poly),
        "top":   np.array(top_poly),
    }

    return V_out, facets
