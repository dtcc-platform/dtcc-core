from dtcc_core.logging import info, warning, error

from typing import TypeVar, TYPE_CHECKING, Union, List

from ....model.geometry import Bounds

if TYPE_CHECKING:
    from ....model.geometry import Mesh,VolumeMesh

    T_Mesh = TypeVar("T_Mesh", bound=Mesh)
    T_VolumeMesh = TypeVar("T_VolumeMesh", bound=VolumeMesh)

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

    def extrude_to_solid(self: "T_Mesh", extrusion_depth: float = None, 
                        base_z: float = None) -> "T_Mesh":
        """
        Extrude this surface mesh downwards to create a solid mesh suitable for 3D printing.
        This method is non-mutating (does not modify data in-place).

        This function creates a solid by:
        1. Finding the boundary edges of the surface
        2. Extruding the boundary downwards to create side walls
        3. Creating a bottom cap at the specified depth

        Args:
            extrusion_depth (float, optional): How far down to extrude. If None, uses 10% of mesh height.
            base_z (float, optional): Z-coordinate for the bottom. If None, uses min_z - extrusion_depth.

        Returns:
            Mesh: A new solid mesh with extruded sides and bottom cap, suitable for 3D printing.
        """
        from dtcc_core.builder.meshing import extrude_surface_to_solid

        solid_mesh = extrude_surface_to_solid(self, extrusion_depth, base_z)

        return solid_mesh

    def create_printable_solid(self: "T_Mesh", extrusion_depth: float = None,
                              base_z: float = None, minimum_thickness: float = 0.001) -> "T_Mesh":
        """
        Create a solid mesh suitable for 3D printing with additional validation and fixes.
        This method is non-mutating (does not modify data in-place).

        Args:
            extrusion_depth (float, optional): Depth of extrusion
            base_z (float, optional): Base Z coordinate
            minimum_thickness (float): Minimum wall thickness for 3D printing

        Returns:
            Mesh: A new solid mesh optimized for 3D printing
        """
        from dtcc_core.builder.meshing import create_printable_surface_mesh

        printable_mesh = create_printable_surface_mesh(self, extrusion_depth, base_z, minimum_thickness)

        return printable_mesh

    def to_cpp(self: "T_Mesh") -> "_dtcc_builder.Mesh":
        """
        Convert the Mesh to a DTCC builder Mesh.

        Returns
        -------
        _dtcc_builder.Mesh
            A DTCC builder Mesh object.
        """
        from dtcc_core.builder.model_conversion import mesh_to_builder_mesh

        return mesh_to_builder_mesh(self.vertices, self.faces, self.markers)
    

class VolumeMeshProcessingMixin:
    def to_cpp(self: "T_VolumeMesh") -> "_dtcc_builder.VolumeMesh":
        """
        Convert the VolumeMesh to a DTCC builder VolumeMesh.

        Returns
        -------
        _dtcc_builder.VolumeMesh
            A DTCC builder VolumeMesh object.
        """
        from dtcc_core.builder.model_conversion import volume_mesh_to_builder_volume_mesh

        return volume_mesh_to_builder_volume_mesh(self.vertices, self.cells, self.markers)