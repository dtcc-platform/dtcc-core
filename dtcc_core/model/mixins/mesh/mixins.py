from dtcc_core.logging import info, warning, error

from typing import TypeVar, TYPE_CHECKING, Union, List

from ....model.geometry import Bounds

if TYPE_CHECKING:
    from ....model.geometry import Mesh

    T_Mesh = TypeVar("T_Mesh", bound=Mesh)


class MeshProcessingMixin:
    """
    Mixin for processing meshes.
    """

    def merge(self: "T_Mesh", other: "T_Mesh") -> "T_Mesh":
        """
        Merge this mesh with another mesh. This method is non-mutating
        (does not modify data in-place).

        Args:
            other (Mesh): The other mesh to merge with.

        Returns:
            Mesh: The a new mesh object that is the result of merging this mesh with the other mesh.
        """
        from dtcc_core.builder.meshing import merge

        merged_mesh = merge(self, other)

        return merged_mesh

    def snap_vertices(self: "T_Mesh", snap_distance: float = 0.01) -> "T_Mesh":
        """
        Snap the vertices of the mesh to a grid defined by the snap distance. This method is non-mutating
        (does not modify data in-place).

        Args:
            snap_distance (float): The distance to snap the vertices to.

        Returns:
            Mesh: A new mesh object with snapped vertices.
        """
        from dtcc_core.builder.meshing import snap_vertices

        snapped_mesh = snap_vertices(self, snap_distance)

        return snapped_mesh
