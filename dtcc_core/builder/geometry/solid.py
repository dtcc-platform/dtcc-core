"""Boundary triangulation of native solids using the existing surface mesher."""

from ...model import Mesh, Solid
from ..register import register_model_method
from ..meshing.meshing import _mesh_surface_collection


@register_model_method
def mesh(
    solid: Solid,
    triangle_size=None,
    weld=False,
    snap=0,
    clean=False,
    mesher: str | None = None,
) -> Mesh:
    """Triangulate all exterior and interior boundary surfaces of a Solid.

    Returns a surface Mesh with copied transform and semantic regions reindexed
    onto triangles. The input is unchanged. Shell partition remains on the
    source Solid; the Mesh contains neither shell topology nor volume cells.
    This operation does not certify closure, manifoldness or geometric validity.

    As with semantic MultiSurface meshing, cleaning, welding/snapping, nested
    transforms and attached field transfer are unsupported and fail explicitly.
    ``triangle_size`` and ``mesher`` have the same meaning as on MultiSurface.mesh.
    """
    return _mesh_surface_collection(solid, triangle_size, weld, snap, clean, mesher)
