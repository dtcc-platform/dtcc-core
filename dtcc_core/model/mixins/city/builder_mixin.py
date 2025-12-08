from dtcc_core.logging import info, warning, error


from pathlib import Path
from typing import Union
from typing import TypeVar, TYPE_CHECKING

from ....model.geometry import PointCloud, Mesh 
from ....model.object import GeometryType
from ....model.values import Raster

if TYPE_CHECKING:
    from ....model.object import City

    T_City = TypeVar("T_City", bound=City)


class CityBuilderMixin:

    def build_terrain(
        self: "T_City",
        pc: PointCloud = None,
        cell_size: float = 2.0,
        build_mesh=True,
        max_triangle_size=5.0,
        smoothing=3,
    ) -> "T_City":
        """
        Build terrain for a city using a point cloud.

        Args:
            self (City): The city object to build terrain for.
            pc (PointCloud): The point cloud to use for building the terrain.
            cell_size (float): The size of the cells in the raster (default is 2).

        Returns:
            City: The city object with the terrain added.
        """
        from dtcc_core.builder import build_terrain_raster, build_terrain_mesh
        from ....model.object import Terrain

        if pc is not None and not isinstance(pc, PointCloud):
            raise ValueError("pc must be a PointCloud object")

        if pc is None:
            pc = self.pointcloud
            if pc is None:
                raise ValueError(
                    "No point cloud provided and city has no point cloud geometry\n"
                )

        if len(pc.points) == 0:
            raise ValueError("Point cloud has no points")

        raster = build_terrain_raster(pc, cell_size=cell_size, ground_only=True)

        terrain = Terrain()
        terrain.add_raster(raster)
        if build_mesh:
            mesh = build_terrain_mesh(
                raster,
                max_mesh_size=max_triangle_size,
                smoothing=smoothing,
            )
            terrain.add_mesh(mesh)

        self.add_terrain(terrain)
        return self

    def build_lod1_buildings(
        self: "T_City",
        default_ground_height: float = 0.0,
        min_building_height: float = 2.5,
        always_use_default=False,
        rebuild=True,
        calculate_heights=True,
        building_height_attribute: str = "height",
    ) -> "T_City":
        """
        Build LOD1 buildings for a city.

        Args:
            self (City): The city object to build LOD1 buildings for.
            default_ground_height (float): The default ground height to use if no terrain is available.
            always_use_default (bool): Whether to always use the default ground height or use the ground_height from
            the terrain if available.
            rebuild (bool): Whether to rebuild the LOD1 buildings if they already exist.
            calculate_heights (bool): Whether to calculate building heights from the point cloud or get it from an attribute.
            building_height_attribute (str): The attribute to use for building heights if calculate_heights is False.

        Returns:
            City: The city object with LOD1 buildings added.
        """
        from dtcc_core.builder import (
            building_heights_from_pointcloud,
            build_lod1_buildings,
        )

        if len(self.buildings) == 0:
            raise ValueError(
                "City has no buildings to build LOD1 geometry for\nload building footprints first."
            )

        if not always_use_default and self.terrain.raster is None:
            info("City has no terrain, generating terrain from point cloud")
            self.build_terrain(pc=self.pointcloud)

        if calculate_heights:
            if self.pointcloud is None:
                raise ValueError(
                    "City has no point cloud geometry\nAdd a point cloud to the city to calculate building heights."
                )

            buildings_with_heights = building_heights_from_pointcloud(
                self.buildings,
                self.pointcloud,
                self.terrain.raster,
                statistical_outlier_remover=True,
                roof_outlier_neighbors=5,
                roof_outlier_margin=1.5,
                overwrite=True,
            )
            self.remove_buildings()
            self.add_buildings(buildings_with_heights)
        else:
            for b in self.buildings:
                footprint = b.lod0
                if footprint is None:
                    warning(f"Building {b.id} has no LOD0 geometry.")
                    continue
                if always_use_default:
                    ground_height = default_ground_height
                else:
                    centroid = footprint.centroid
                    if np.isnan(centroid[0]) or np.isnan(centroid[1]):
                        warning(f"Building {b.id} has an invalid centroid.")
                        ground_height = default_ground_height
                    else:
                        ground_height = self.terrain.raster.terrain.get_value(
                            centroid[0], centroid[1]
                        )
                b.attributes["ground_height"] = ground_height
                height = b.attributes.get(building_height_attribute, None)
                if height is None:
                    warning(
                        f"Building {b.id} has no {building_height_attribute} attribute. Using default ground height."
                    )
                    height = default_ground_height
                if height < min_building_height:
                    warning(
                        f"Building {b.id} has a height of {height}, which is less than the minimum building height of {min_building_height}. Setting height to {min_building_height}."
                    )
                    height = min_building_height
                footprint.set_z(ground_height + height)
                b.attributes["height"] = height
        lod1_buildings = build_lod1_buildings(
            self.buildings,
            default_ground_height=default_ground_height,
            always_use_default_ground=always_use_default,
            rebuild=rebuild,
        )
        self.remove_buildings()
        self.add_buildings(lod1_buildings)

        return self

    def build_surface_mesh(
        self: "T_City",
        lod: GeometryType | list[GeometryType] = GeometryType.LOD1 ,
        min_building_detail: float = 0.5,
        min_building_area: float = 15.0,
        merge_buildings: bool = True,
        merge_tolerance: float = 0.5,
        building_mesh_triangle_size: float = 5.0,
        max_mesh_size: float = 10.0,
        min_mesh_angle: float = 25.0,
        merge_meshes: bool = True,
        smoothing: int = 0,
        sort_triangles: bool = False,
        treat_lod0_as_holes: bool = False,
    ) -> Mesh:
        """
            Build a city surface mesh from the buildings and terrain.

             Parameters
        ----------
        `min_building_detail` : float, optional
            The minimum detail of the buildin to resolve, by default 0.5.
        `min_building_area` : float, optional
            The smallest building to include, by default 15.0.
        `merge_buildings` : bool, optional
            merge building footprints, by default True.
        `max_mesh_size` : float, optional
            The maximum size of the mesh, by default 1.0.
        `min_mesh_angle` : float, optional
            The minimum angle of the mesh, by default 30.0.
        `merge_meshes` : bool, optional
            Whether to merge the meshes to a single mesh, by default True.

        `smoothing` : float, optional
            The smoothing of the mesh, by default 0.0.

        Returns
        -------
        `Mesh` : The city surface mesh.
        """
        from dtcc_core.builder import build_city_mesh

        surface_mesh = build_city_mesh(
            self,
            lod=lod,
            min_building_detail=min_building_detail,
            min_building_area=min_building_area,
            merge_buildings=merge_buildings,
            merge_tolerance=merge_tolerance,
            building_mesh_triangle_size=building_mesh_triangle_size,
            max_mesh_size=max_mesh_size,
            min_mesh_angle=min_mesh_angle,
            merge_meshes=merge_meshes,
            smoothing=smoothing,
            sort_triangles=sort_triangles,
            treat_lod0_as_holes=treat_lod0_as_holes,
        )
        return surface_mesh
