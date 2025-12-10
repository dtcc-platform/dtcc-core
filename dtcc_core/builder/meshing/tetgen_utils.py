
import numpy as np
from ...model import Mesh


def get_east_boundary_vertices(vertices, xmax=None, tol=1e-3):
    """
    Return indices on the east boundary (x near ``xmax``), sorted south to north.

    Parameters
    ----------
    vertices : array_like, shape (N, 2) or (N, 3)
        Vertex coordinates.
    xmax : float, optional
        Boundary reference. If ``None``, uses ``vertices[:, 0].max()``.
    tol : float, optional
        Tolerance for boundary membership (default is 1e-3).

    Returns
    -------
    numpy.ndarray
        Vertex indices sorted south to north.

    Notes
    -----
    Sorting:
    - 2D: by ``y`` ascending
    - 3D: by ``(y, then z)`` ascending

    Building the east wall polygon with outward normal ``+x``:
    let ``g`` be indices on ground (``z = zmin``), ``r`` the corresponding roof
    indices (``z = zmax``). Use CCW order as seen from ``+x``:
    ``[ g (south to north), r (north to south) ]``.
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
    Return indices on the west boundary (x near ``xmin``), sorted north to south.

    Parameters
    ----------
    vertices : array_like, shape (N, 2) or (N, 3)
        Vertex coordinates.
    xmin : float, optional
        Boundary reference. If ``None``, uses ``vertices[:, 0].min()``.
    tol : float, optional
        Tolerance for boundary membership (default is 1e-3).

    Returns
    -------
    numpy.ndarray
        Vertex indices sorted north to south.

    Notes
    -----
    Sorting:
    - 2D: by ``y`` descending
    - 3D: by ``(y`` descending, then ``z`` ascending)

    For the west wall (outward normal ``-x``), a CCW ordering seen from ``-x`` is
    ``[ g (north to south), r (south to north) ]``, which yields outward normal
    ``-x``.
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
    Return indices on the south boundary (y near ``ymin``), sorted west to east.

    Parameters
    ----------
    vertices : array_like, shape (N, 2) or (N, 3)
        Vertex coordinates.
    ymin : float, optional
        Boundary reference. If ``None``, uses ``vertices[:, 1].min()``.
    tol : float, optional
        Tolerance for boundary membership (default is 1e-3).

    Returns
    -------
    numpy.ndarray
        Vertex indices sorted west to east.

    Notes
    -----
    Sorting:
    - 2D: by ``x`` ascending
    - 3D: by ``(x, then z)`` ascending

    For the south wall (outward normal ``-y``), viewed from ``-y``, a CCW ordering is
    ``[ g (west to east), r (east to west) ]``; reversing the roof indices closes
    the quad and keeps the normal pointing outward.
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
    Return indices on the north boundary (y near ``ymax``), sorted east to west.

    Parameters
    ----------
    vertices : array_like, shape (N, 2) or (N, 3)
        Vertex coordinates.
    ymax : float, optional
        Boundary reference. If ``None``, uses ``vertices[:, 1].max()``.
    tol : float, optional
        Tolerance for boundary membership (default is 1e-3).

    Returns
    -------
    numpy.ndarray
        Vertex indices sorted east to west.

    Notes
    -----
    Sorting:
    - 2D: by ``x`` descending
    - 3D: by ``(x`` descending, then ``z`` ascending)

    For the north wall (outward normal ``+y``), viewed from ``+y``, a CCW ordering is
    ``[ g (east to west), r (west to east) ]``; reversing the roof segment closes
    the quad and yields outward normal ``+y``.
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
    Build PLC facet polygons for a rectangular box around the mesh domain.

    Polygons are wound CCW as seen from outside so normals point outward
    (right-hand rule):
    - South -> outward ``-y``
    - East  -> outward ``+x``
    - North -> outward ``+y``
    - West  -> outward ``-x``
    - Top   -> outward ``+z``

    Parameters
    ----------
    mesh : Mesh
        Input mesh used to derive the domain bounds.
    top_height : float, optional
        Height of the top cap above ``zmin``; defaults to 100.0.
    tol : float, optional
        Tolerance for boundary membership (default is 1e-3).

    Returns
    -------
    numpy.ndarray
        Vertices with the four added top points appended at the end.
    dict[str, numpy.ndarray]
        Indices of the five polygons with outward normals:
        ``{"south": [...], "east": [...], "north": [...], "west": [...], "top": [...]}``.
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
