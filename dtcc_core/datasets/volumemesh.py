import dtcc_core
from dtcc_core.model import City, Bounds, VolumeMesh
from typing import Literal, Optional, Dict, Any
from pydantic import BaseModel, Field

from .dataset import DatasetDescriptor, DatasetBaseArgs


class VolumeMeshArgs(DatasetBaseArgs):
    max_mesh_size: float = Field(
        25.0, description="Maximum mesh size (h parameter) in meters"
    )
    domain_height: float = Field(
        80.0, description="Height of the computational domain (H parameter) in meters"
    )
    raster_cell_size: float = Field(
        2.0, description="Cell size for terrain raster in meters"
    )
    raster_radius: float = Field(
        3.0, description="Radius for terrain raster interpolation"
    )
    remove_outliers: bool = Field(
        True, description="Whether to remove global outliers from point cloud"
    )
    outlier_threshold: float = Field(
        3.0, description="Threshold for outlier removal (standard deviations)"
    )
    boundary_face_markers: bool = Field(
        True, description="Whether to add boundary face markers to the mesh"
    )
    max_volume: Optional[float] = Field(
        None, description="Maximum tetrahedron volume (defaults to max_mesh_size if not set)"
    )
    tetgen_extra: str = Field(
        " -VV", description="Extra switches to pass to TetGen (e.g., ' -VV' for verbose output)"
    )
    format: Optional[Literal["xdmf", "vtu"]] = Field(
        None, description="Output file format"
    )


class VolumeMeshDataset(DatasetDescriptor):
    name = "volumemesh"
    description = "Generate a tetrahedral volume mesh from point cloud and building data, suitable for CFD/FEM simulations."
    ArgsModel = VolumeMeshArgs

    def build(self, args: VolumeMeshArgs):
        """Build a volume mesh from point cloud and building data.
        
        This method:
        1. Downloads point cloud and building footprints for the given bounds
        2. Removes outliers from the point cloud
        3. Builds a terrain raster
        4. Extracts roof points and computes building heights
        5. Creates a city model with terrain and buildings
        6. Generates a tetrahedral volume mesh suitable for simulations
        
        Args:
            args: Validated arguments containing bounds and meshing parameters
            
        Returns:
            VolumeMesh object or bytes if format is specified
        """
        bounds = self.parse_bounds(args.bounds)
        
        # Download pointcloud and building footprints
        pointcloud = dtcc_core.io.data.download_pointcloud(bounds=bounds)
        buildings = dtcc_core.io.data.download_footprints(bounds=bounds)
        
        # Remove global outliers
        if args.remove_outliers:
            pointcloud = pointcloud.remove_global_outliers(args.outlier_threshold)
        
        # Build terrain raster
        raster = dtcc_core.builder.build_terrain_raster(
            pointcloud, 
            cell_size=args.raster_cell_size, 
            radius=args.raster_radius, 
            ground_only=True
        )
        
        # Extract roof points and compute building heights
        buildings = dtcc_core.builder.extract_roof_points(buildings, pointcloud)
        buildings = dtcc_core.builder.compute_building_heights(
            buildings, raster, overwrite=True
        )
        
        # Create city and add geometries
        city = City()
        city.add_terrain(raster)
        city.add_buildings(buildings, remove_outside_terrain=True)
        
        # Build city volume mesh
        max_vol = args.max_volume if args.max_volume is not None else args.max_mesh_size
        volume_mesh = dtcc_core.builder.build_city_volume_mesh(
            city,
            max_mesh_size=args.max_mesh_size,
            domain_height=args.domain_height,
            boundary_face_markers=args.boundary_face_markers,
            tetgen_switches={
                "max_volume": max_vol,
                "extra": args.tetgen_extra,
            },
        )
        
        # Return as bytes in requested format or as object
        if args.format is not None:
            return self.export_to_bytes(volume_mesh, args.format)
        
        return volume_mesh
