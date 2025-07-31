from ...model import (
    Mesh,
    VolumeMesh,
    Building,
    Terrain,
    City,
    Surface,
    MultiSurface,
    GeometryType,
    PointCloud,
)

from ..model_conversion import (
    create_builder_polygon,
    create_builder_surface,
    create_builder_multisurface,
    builder_mesh_to_mesh,
    builder_volume_mesh_to_volume_mesh,
    mesh_to_builder_mesh,
    create_builder_city,
    raster_to_builder_gridfield,
)

from .. import _dtcc_builder

from ..polygons.polygons import (
    polygon_merger,
    simplify_polygon,
    remove_slivers,
    fix_clearance,
)

from .buildings import (
    extract_roof_points,
    compute_building_heights,
)

from .terrain import (
    build_terrain_mesh,
    build_terrain_raster,
)

from ..building.modify import (
    merge_building_footprints,
    simplify_building_footprints,
    fix_building_footprint_clearance,
)

from ..meshing.convert import mesh_to_raster

from ..logging import debug, info, warning, error


def build_city_mesh(
    city: City,
    lod: GeometryType = GeometryType.LOD1,
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
) -> Mesh:
    """
    Build a mesh from the surfaces of the buildings in the city.

    Parameters
    ----------
    `city` : model.City
        The city to build the mesh from.
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
    `model.Mesh`
    """
    buildings = city.buildings
    if merge_buildings:
        info(f"Merging {len(buildings)} buildings...")
        merged_buildings = merge_building_footprints(
            buildings, lod, min_area=min_building_area
        )
        simplifed_footprints = simplify_building_footprints(
            merged_buildings, min_building_detail / 2, lod=GeometryType.LOD0
        )
        clearance_fix = fix_building_footprint_clearance(
            simplifed_footprints, min_building_detail
        )
        building_footprints = [
            b.get_footprint(GeometryType.LOD0) for b in clearance_fix
        ]
        info(f"After merging: {len(building_footprints)} buildings.")
    else:
        building_footprints = [b.get_footprint(lod) for b in buildings]

    subdomain_resolution = [building_mesh_triangle_size] * len(building_footprints)

    terrain = city.terrain
    if terrain is None:
        raise ValueError("City has no terrain data. Please compute terrain first.")
    terrain_raster = terrain.raster
    terrain_mesh = terrain.mesh
    if terrain_raster is None and terrain_mesh is None:
        raise ValueError("City terrain has no data. Please compute terrain first.")
    if terrain_raster is None and terrain_mesh is not None:
        terrain_raster = mesh_to_raster(terrain_mesh, cell_size=max_mesh_size)
    builder_dem = raster_to_builder_gridfield(terrain_raster)

    builder_surfaces = [
        create_builder_surface(p) for p in building_footprints if p is not None
    ]
    builder_mesh = _dtcc_builder.build_city_surface_mesh(
        builder_surfaces,
        subdomain_resolution,
        builder_dem,
        max_mesh_size,
        min_mesh_angle,
        smoothing,
        merge_meshes,
        sort_triangles,
    )
    if merge_meshes:
        result_mesh = builder_mesh_to_mesh(builder_mesh[0])
    else:
        result_mesh = [builder_mesh_to_mesh(bm) for bm in builder_mesh]
    return result_mesh


def build_volume_mesh(
    city: City,
    lod: GeometryType = GeometryType.LOD1,
    domain_height: float = 100.0,
    max_mesh_size: float = 10.0,
    merge_buildings: bool = True,
    boundary_face_markers: bool = True,
) -> VolumeMesh:
    """
    Build a 3D tetrahedral volume mesh for a city terrain with embedded building volumes.

    This function generates a ground mesh from the city terrain and extrudes building
    footprints to produce a full volume mesh, optionally merging adjacent buildings
    and marking boundary faces.

    Parameters
    ----------
    city : City
        City object containing terrain and building data. Terrain must have either
        a raster or mesh representation.
    lod : GeometryType, optional
        Level of detail for building footprints. Defaults to GeometryType.LOD1.
    domain_height : float, optional
        Height of the mesh domain above the terrain (in the same units as city coordinates).
        Defaults to 100.0.
    max_mesh_size : float, optional
        Maximum allowed mesh element size. Used both for ground mesh cell sizing
        and to cap subdomain resolutions. Defaults to 10.0.
    merge_buildings : bool, optional
        If True, merge adjacent building footprints into larger blocks before meshing.
        Defaults to True.
    boundary_face_markers : bool, optional
        If True, add integer markers to the boundary faces of the volume mesh as a
        post-processing step. Defaults to True. See Notes for marker conventions.

    Returns
    -------
    VolumeMesh
        A VolumeMesh instance representing the 3D tetrahedral mesh of the city domain, including
        building volumes.

    Raises
    ------
    ValueError
        If the city has no terrain data (neither raster nor mesh).
    ValueError
        If the terrain object exists but has no usable raster or mesh data.

    Boundary Face Markers
    ---------------------
    When `boundary_face_markers=True`, integer markers are added as follows (for
    a domain containing N buildings):

    - `0` to `N-1`:  Wall faces of the N buildings  
    - `N` to `2*N-1`:  Roof faces of the N buildings  
    - `-1`:  Ground (terrain) faces  
    - `-2`:  Top faces of the volume domain  
    - `-3`, `-4`, `-5`, `-6`:  The four vertical boundary faces of the domain  

    Notes
    -----
    - Building footprints are extracted at the specified `lod` level (or LOD0 if
      `merge_buildings` is True), optionally merged and simplified using internal
      area/detail thresholds.
    - Subdomain resolution for each building is set to the minimum of its height
      and `max_mesh_size`.
    - Ground mesh is built via the internal DTCC builder, and building surfaces
      are extruded into the volume domain of height `domain_height`.
    - Mesh smoothing and quality parameters (angles, iterations, tolerances, aspect
      ratios) are applied internally.

    Examples
    --------
    >>> mesh = build_volume_mesh(my_city,
    ...                          lod=GeometryType.LOD1,
    ...                          domain_height=150.0,
    ...                          max_mesh_size=5.0,
    ...                          merge_buildings=False,
    ...                          boundary_face_markers=True)
    """

    # FIXME: Where do we set these parameters?
    min_building_area = 10.0
    min_building_detail = 0.5
    min_mesh_angle = 30
    smoother_max_iterations = 5000
    smoothing_relative_tolerance = 0.005
    aspect_ratio_threshold = 10.0
    debug_step = 7
    

    buildings = city.buildings
    if not buildings:
        warning("City has no buildings.")

    if merge_buildings:
        info(f"Merging {len(buildings)} buildings...")
        merged_buildings = merge_building_footprints(
            buildings, GeometryType.LOD0, min_area=min_building_area
        )
        simplifed_footprints = simplify_building_footprints(
            merged_buildings, min_building_detail / 2, lod=GeometryType.LOD0
        )
        clearance_fix = fix_building_footprint_clearance(
            simplifed_footprints, min_building_detail
        )
        
        building_footprints = [
            b.get_footprint(GeometryType.LOD0) for b in clearance_fix
        ]
        info(f"After merging: {len(building_footprints)} buildings.")
    else:
        building_footprints = [b.get_footprint(lod) for b in buildings]
    

     # Set subdomain resolution to half the building height
    subdomain_resolution = [
        min(building.height, max_mesh_size) for building in buildings
    ]

    terrain = city.terrain
    if terrain is None:
        raise ValueError("City has no terrain data. Please compute terrain first.")
    terrain_raster = terrain.raster
    terrain_mesh = terrain.mesh
    if terrain_raster is None and terrain_mesh is None:
        raise ValueError("City terrain has no data. Please compute terrain first.")
    if terrain_raster is None and terrain_mesh is not None:
        terrain_raster = mesh_to_raster(terrain_mesh, cell_size=max_mesh_size)
   
    
    
    # Convert from Python to C++
    _building_polygons = [
        create_builder_polygon(footprint.to_polygon())
        for footprint in building_footprints
        if footprint is not None
    ]
    # FIXME: Pass bounds as argument (not xmin, ymin, xmax, ymax).

    # Build ground mesh
    _ground_mesh = _dtcc_builder.build_ground_mesh(
        _building_polygons,
        subdomain_resolution,
        terrain.bounds.xmin,
        terrain.bounds.ymin,
        terrain.bounds.xmax,
        terrain.bounds.ymax,
        max_mesh_size,
        min_mesh_angle,
        True,
    )

    # FIXME: Should not need to convert from C++ to Python mesh.
    # Convert from Python to C++
    
    _surfaces = [
        create_builder_surface(footprint)
        for footprint in building_footprints
        if footprint is not None
    ]
    
    # Convert from C++ to Python
    # ground_mesh = builder_mesh_to_mesh(_ground_mesh)
    _dem = raster_to_builder_gridfield(terrain.raster)

    # Create volume mesh builder
    volume_mesh_builder = _dtcc_builder.VolumeMeshBuilder(
        _surfaces, _dem, _ground_mesh, domain_height
    )

    # FIXME: How do we handle parameters?

    # Build volume mesh
    _volume_mesh = volume_mesh_builder.build(
        smoother_max_iterations,
        smoothing_relative_tolerance,
        0.0,
        aspect_ratio_threshold,
        debug_step,
    )
    volume_mesh = builder_volume_mesh_to_volume_mesh(_volume_mesh)


    if boundary_face_markers:
        boundary_face_markers = _dtcc_builder.compute_boundary_face_markers(_volume_mesh)
        if boundary_face_markers is not None:
            volume_mesh.boundary_markers = boundary_face_markers
    

    return volume_mesh
