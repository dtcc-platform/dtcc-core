from ...model import Mesh, VolumeMesh
import numpy as np
from typing import Any, Dict, Optional, Tuple, Union

from . import tetgen_utils

from ..logging import info,warning

HAS_TETGEN = False
_tetgen_switch_module = None
try:
    import dtcc_tetgen_wrapper as tetwrap

    _tetgen_switch_module = tetwrap.switches
    HAS_TETGEN = True
    info("TetGen is available for volume meshing.")
except ImportError:
    _tetgen_switch_module = None
    warning("TetGen not available. Volume meshing fallback to dtcc base method.")

def is_tetgen_available() -> bool:
    return HAS_TETGEN


def get_default_tetgen_switches() -> Dict[str, Any]:
    """
    Return a fresh copy of TetGen switch defaults if TetGen is available.
    """
    if not HAS_TETGEN or _tetgen_switch_module is None:
        raise RuntimeError("TetGen switch defaults requested but TetGen is not available.")
    return _tetgen_switch_module.tetgen_defaults()


def build_volume_mesh(mesh: Mesh, 
                      build_top_sidewalls: bool=True, 
                      top_height: float=100.0, 
                      return_boundary_faces: bool=True, 
                      switches_params: Optional[Dict[str, Any]] = None, 
                      switches_overrides: Optional[Dict[str, Any]] = None) -> Union[VolumeMesh, Tuple[VolumeMesh, Optional[np.ndarray]]] :
    """
    Build a volume mesh from a surface mesh using TetGen via tetwrap.

    Parameters
    ----------
    mesh : dtcc.VolumeMesh
        The input surface mesh.
    switches_params : dict, optional
        Parameters for TetGen switches.
    switches_overrides : dict, optional
        Overrides for TetGen switches.

    Returns
    -------
    vertices : ndarray
        The vertices of the volume mesh.
    cells : ndarray
        The cells (tetrahedra) of the volume mesh.
    """
    if not isinstance(mesh, Mesh):
        raise TypeError("Input must be a dtcc.Mesh instance.")

    if mesh.faces is None or len(mesh.faces) == 0:
        raise ValueError("Input mesh must have faces defined.")

    if mesh.markers is None or len(mesh.markers) == 0:
        raise ValueError("Input mesh must have face markers defined.")

    b_facets = None
    if build_top_sidewalls:
        new_vertices, boundary_facets = tetgen_utils.compute_boundary_facets(mesh, top_height=top_height)
        mesh = Mesh(vertices=new_vertices, faces=mesh.faces ,markers=mesh.markers)
        b_facets = [facet for facet in boundary_facets.values()]
    
    if b_facets is None:
        raise ValueError(
            "TetGen volume meshing requires boundary facets. "
            "Set build_top_sidewalls=True or provide facets via future extensions."
        )

    # Prepare TetGen switches
    base_switches: Dict[str, Any] = {}
    if _tetgen_switch_module is not None:
        base_switches = get_default_tetgen_switches()
    if switches_params:
        base_switches.update(switches_params)

    # Call tetwrap to build the volume mesh
    tetgen_out: tetwrap.TetwrapIO = tetwrap.tetrahedralize(vertices= mesh.vertices, 
                                             faces= mesh.faces,
                                             face_markers = mesh.markers,
                                             boundary_facets= b_facets, 
                                             switches_params=base_switches,
                                             switches_overrides=switches_overrides,
                                             return_io = True,
                                             return_faces = False,
                                             return_boundary_faces = return_boundary_faces,
                                             return_edges = False,
                                             return_neighbors = False)
    
    vertices = np.asarray(tetgen_out.points)
    cells = np.asarray(tetgen_out.tets)
    volume_mesh = VolumeMesh(vertices=vertices, cells=cells)
    if tetgen_out.boundary_tri_faces is not None:
        boundary_faces = np.asarray(tetgen_out.boundary_tri_faces)
        volume_mesh.boundary_faces = boundary_faces
    if tetgen_out.boundary_tri_markers is not None:
        boundary_markers = np.asarray(tetgen_out.boundary_tri_markers)
        volume_mesh.boundary_markers = boundary_markers

    return volume_mesh
