from typing import Any, Dict, Optional, Tuple, Union

import numpy as np

from ...model import Mesh, VolumeMesh
from . import tetgen_utils

from ..logging import debug, info, warning

HAS_TETGEN = False
_tetgen_switch_module = None
BOUNDARY_FACET_MARKERS: Dict[str, int] = {
    "top": -2,
    "west": -3,
    "east": -4,
    "south": -5,
    "north": -6,
}
try:
    import dtcc_tetgen_wrapper as tetwrap

    _tetgen_switch_module = tetwrap.switches
    HAS_TETGEN = True
    info("TetGen is available for volume meshing.")
except ImportError:
    _tetgen_switch_module = None
    warning("TetGen not available. Volume meshing fallback to dtcc base method.")


def is_tetgen_available() -> bool:
    """
    Indicate whether the TetGen wrapper is importable.

    Returns
    -------
    bool
        ``True`` when ``dtcc_tetgen_wrapper`` is available and switches are
        initialized; otherwise ``False``.
    """
    return HAS_TETGEN


def get_default_tetgen_switches() -> Dict[str, Any]:
    """
    Return a fresh copy of TetGen switch defaults if TetGen is available.
    """
    if not HAS_TETGEN or _tetgen_switch_module is None:
        raise RuntimeError("TetGen switch defaults requested but TetGen is not available.")
    return _tetgen_switch_module.tetgen_defaults()


def build_volume_mesh(
    mesh: Mesh,
    build_top_sidewalls: bool = True,
    top_height: float = 100.0,
    return_boundary_faces: bool = True,
    closure_mesh: Optional[Mesh] = None,
    top_cap_backend: str = "auto",
    top_cap_max_mesh_size: Optional[float] = None,
    top_cap_min_mesh_angle: float = 25.0,
    switches_params: Optional[Dict[str, Any]] = None,
    switches_overrides: Optional[Dict[str, Any]] = None,
) -> Union[VolumeMesh, Tuple[VolumeMesh, Optional[np.ndarray]]]:
    """
    Build a tetrahedral volume mesh from a surface mesh using TetGen.

    Parameters
    ----------
    mesh : Mesh
        Input surface mesh with faces and face markers.
    build_top_sidewalls : bool, optional
        Whether to add top cap and side wall facets around the domain before meshing.
    top_height : float, optional
        Height of the top cap above ``zmin`` when ``build_top_sidewalls`` is True.
    return_boundary_faces : bool, optional
        Request TetGen to return boundary faces; stored on the resulting VolumeMesh.
    closure_mesh : Mesh, optional
        Optional flat ground mesh retained for compatibility with the builder
        API. The live PLC path remeshes the top cap independently from the
        shell's outer boundary ring instead of copying ``closure_mesh``.
    top_cap_backend : str, optional
        2D backend used when re-triangulating the top cap from the shell's
        outer-domain boundary. Defaults to ``"auto"``.
    top_cap_max_mesh_size : float, optional
        Maximum target edge length for the re-triangulated top cap. When ``None``,
        only the prescribed outer boundary vertices constrain the top cap.
    top_cap_min_mesh_angle : float, optional
        Minimum angle target for the re-triangulated top cap.
    switches_params : dict, optional
        Base parameters passed to TetGen switches.
    switches_overrides : dict, optional
        Overrides applied to the TetGen switches after ``switches_params``.

    Notes
    -----
    When boundary faces are requested, the returned boundary markers follow the
    dtcc CFD convention: ``-2`` top, ``-3`` west/xmin, ``-4`` east/xmax,
    ``-5`` south/ymin, ``-6`` north/ymax. Input shell-face markers are
    preserved for the original surface triangles.

    Returns
    -------
    VolumeMesh
        Volume mesh populated with tetrahedra and, when requested, boundary faces and markers.

    Raises
    ------
    TypeError
        If ``mesh`` is not a ``Mesh`` instance.
    ValueError
        If the input mesh lacks faces, face markers, or required boundary facets.
    """
    if not isinstance(mesh, Mesh):
        raise TypeError("Input must be a dtcc.Mesh instance.")

    if mesh.faces is None or len(mesh.faces) == 0:
        raise ValueError("Input mesh must have faces defined.")

    if mesh.markers is None or len(mesh.markers) == 0:
        raise ValueError("Input mesh must have face markers defined.")
    if len(mesh.markers) != len(mesh.faces):
        raise ValueError("Input mesh must have one face marker per face.")

    b_facets = None
    boundary_facet_markers = None
    diagnostic_boundary_facets = None
    named_boundary_facets = None
    if build_top_sidewalls:
        if closure_mesh is not None:
            (
                new_vertices,
                oriented_faces,
                boundary_facets,
                boundary_facet_markers,
                _,
            ) = tetgen_utils.compute_oriented_boundary_plc(
                mesh,
                closure_mesh,
                top_height=top_height,
                top_cap_backend=top_cap_backend,
                top_cap_max_mesh_size=top_cap_max_mesh_size,
                top_cap_min_mesh_angle=top_cap_min_mesh_angle,
            )
            mesh = Mesh(vertices=new_vertices, faces=oriented_faces, markers=mesh.markers)
        else:
            new_vertices, boundary_facets = tetgen_utils.compute_boundary_facets(
                mesh, top_height=top_height
            )
            boundary_facet_markers = {
                name: BOUNDARY_FACET_MARKERS[name] for name in boundary_facets
            }
            mesh = Mesh(vertices=new_vertices, faces=mesh.faces, markers=mesh.markers)
        if isinstance(boundary_facets, dict):
            named_boundary_facets = dict(boundary_facets)
            b_facets = [facet for facet in boundary_facets.values()]
            diagnostic_boundary_facets = named_boundary_facets
        else:
            diagnostic_boundary_facets = {
                f"facet_{i}": facet for i, facet in enumerate(boundary_facets)
            }
            b_facets = list(boundary_facets)

    if not b_facets:
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
    effective_switches = dict(base_switches)
    if switches_overrides:
        effective_switches.update(switches_overrides)

    plc_diagnostics = tetgen_utils.inspect_tetgen_plc(
        mesh.vertices,
        mesh.faces,
        diagnostic_boundary_facets if diagnostic_boundary_facets is not None else b_facets,
    )
    debug(tetgen_utils.format_tetgen_plc_diagnostics(plc_diagnostics))
    if plc_diagnostics.errors:
        summary = "; ".join(plc_diagnostics.errors[:3])
        if len(plc_diagnostics.errors) > 3:
            summary += f"; ... ({len(plc_diagnostics.errors)} total)"
        raise ValueError(f"TetGen PLC precheck failed: {summary}")
    if plc_diagnostics.warnings:
        for msg in plc_diagnostics.warnings:
            warning(f"TetGen PLC precheck: {msg}")
        if effective_switches.get("preserve_surface"):
            warning(
                "TetGen PLC precheck found shell-quality risks while preserve_surface=True; "
                "TetGen will inherit those boundary triangles."
            )
        else:
            warning(
                "TetGen PLC precheck found shell-quality risks while boundary refinement is enabled; "
                "TetGen refinement may be unstable on this PLC."
            )

    # Call tetwrap to build the volume mesh
    tetgen_out: tetwrap.TetwrapIO = tetwrap.tetrahedralize(
        vertices=mesh.vertices,
        faces=mesh.faces,
        face_markers=mesh.markers,
        boundary_facets=named_boundary_facets if named_boundary_facets is not None else b_facets,
        boundary_facet_markers=boundary_facet_markers,
        switches_params=base_switches,
        switches_overrides=switches_overrides,
        return_io=True,
        return_faces=False,
        return_boundary_faces=return_boundary_faces,
        return_edges=False,
        return_neighbors=False,
    )

    vertices = np.asarray(tetgen_out.points)
    cells = np.asarray(tetgen_out.tets)
    if cells.ndim != 2 or cells.shape[1] != 4 or cells.shape[0] == 0:
        raise RuntimeError(
            "TetGen returned no linear tetrahedral cells. "
            f"Got cells shape {cells.shape!r}."
        )
    volume_mesh = VolumeMesh(vertices=vertices, cells=cells)
    if tetgen_out.boundary_tri_faces is not None:
        boundary_faces = np.asarray(tetgen_out.boundary_tri_faces)
        volume_mesh.boundary_faces = boundary_faces
    if tetgen_out.boundary_tri_markers is not None:
        boundary_markers = np.asarray(tetgen_out.boundary_tri_markers)
        volume_mesh.boundary_markers = boundary_markers

    return volume_mesh
