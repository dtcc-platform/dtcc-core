from typing import Any, Dict, Optional

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

_LOD_PRIORITY: dict[GeometryType, int] = {
    GeometryType.LOD0: 0,
    GeometryType.LOD1: 1,
    GeometryType.LOD2: 2,
    GeometryType.LOD3: 3,
}

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

from ..meshing.tetgen import (
    build_volume_mesh as tetgen_build_volume_mesh,
    get_default_tetgen_switches,
    is_tetgen_available,
)

def build_city_mesh(
    city: City,
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
    Build a mesh from the surfaces of the buildings in the city.

    Parameters
    ----------
    `city` : model.City
        The city to build the mesh from.
    `lod` : GeometryType or list of GeometryType, optional
        The level of detail (meshing directive) for each building. If a single value is provided,
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
    `treat_lod0_as_holes` : bool, optional
        When True, building directives resolved to LOD0 are sent to the mesher
        as hole surfaces instead of meshed buildings.

    Returns
    -------
    `model.Mesh`
    """

    def compose_index_map(
        parent_map: list[list[int]], child_map: list[list[int]]
    ) -> list[list[int]]:
        """
        Combine index maps so that the resulting entries always reference
        the original building indices.
        """
        if not child_map:
            return []
        if parent_map is None:
            raise ValueError("Parent index map is undefined for composition.")
        composed: list[list[int]] = []
        for child_indices in child_map:
            combined: list[int] = []
            for idx in child_indices:
                if idx < 0 or idx >= len(parent_map):
                    warning(
                        f"Index map mismatch: child index {idx} outside parent range {len(parent_map)}."
                    )
                    continue
                combined.extend(parent_map[idx])
            if not combined:
                raise ValueError(
                    "Failed to compose index maps: child entry produced no original indices."
                )
            composed.append(combined)
        return composed

    def lod_from_index_map(index_map: list[list[int]]) -> list[GeometryType]:
        """
        Derive a single LoD directive per processed building by reducing the
        directives of the contributing original buildings.
        """
        reduced: list[GeometryType] = []
        for indices in index_map:
            if not indices:
                raise ValueError("Index map entry is empty; cannot determine LoD.")
            subset = [lod[i] for i in indices]
            reduced.append(min(subset, key=lambda x: _LOD_PRIORITY[x]))
        return reduced

    buildings = city.buildings
    
    n_buildings = len(buildings)
    if isinstance(lod, GeometryType):
        # single value -> broadcast to all
        lod = [lod] * n_buildings
    elif isinstance(lod, (list, tuple)):
        if len(lod) != n_buildings:
            raise ValueError(f"lod list length {len(lod)} != number of buildings {n_buildings}")
        if not all(isinstance(x, GeometryType) for x in lod):
            raise TypeError("all elements in lod list must be GeometryType instances")
        lod = list(lod)
    else:
        raise TypeError(
            f"lod must be a single GeometryType or a list/tuple of {n_buildings} GeometryType values, "
            f"got {type(lod).__name__}"
        )
    
    processed_buildings = list(buildings)
    current_index_map: list[list[int]] | None = None
    
    if merge_buildings:
        info(f"Merging {len(buildings)} buildings...")
        merged_buildings, index_map = merge_building_footprints(
            buildings, 
            lod=GeometryType.LOD0,
            max_distance= merge_tolerance, 
            min_area=min_building_area,
            return_index_map= True
        )
        current_index_map = index_map
        processed_buildings = merged_buildings
        
        info(f"Number of buildings after merging: {len(merged_buildings)}")
        simplifed_footprints, simplify_index_map = simplify_building_footprints(
            processed_buildings,
            min_building_detail / 2,
            lod=GeometryType.LOD0,
            return_index_map=True,
        )
        current_index_map = compose_index_map(current_index_map, simplify_index_map)
        info(f"Number of buildings after simplification: {len(simplifed_footprints)}")
        processed_buildings = simplifed_footprints

        clearance_fix, clearance_index_map = fix_building_footprint_clearance(
            processed_buildings,
            min_building_detail,
            return_index_map=True,
        )
        current_index_map = compose_index_map(current_index_map, clearance_index_map)
        info("Number of buildings after clearance fix: ", len(clearance_fix))
        processed_buildings = clearance_fix
        info(f"After merging: {len(processed_buildings)} buildings.")
    
    if merge_buildings:
        target_lods = (
            lod_from_index_map(current_index_map) if current_index_map is not None else []
        )
        building_footprints = [
            b.get_footprint(GeometryType.LOD0) for b in processed_buildings
        ]
    else:
        target_lods = lod
        building_footprints = [
            b.get_footprint(b_lod) for b, b_lod in zip(processed_buildings, target_lods)
        ]
    
    base_resolution = [building_mesh_triangle_size] * len(building_footprints)
    building_surfaces = []
    hole_surfaces = []
    building_resolution = []
    building_lod_switches = []
    default_priority = _LOD_PRIORITY[GeometryType.LOD3]

    for footprint, resolution, lod_value in zip(
        building_footprints, base_resolution, target_lods
    ):
        if footprint is None:
            continue
        builder_surface = create_builder_surface(footprint)
        if treat_lod0_as_holes and lod_value == GeometryType.LOD0:
            hole_surfaces.append(builder_surface)
            continue
        building_surfaces.append(builder_surface)
        building_resolution.append(resolution)
        building_lod_switches.append(_LOD_PRIORITY.get(lod_value, default_priority))

    if not building_surfaces and not hole_surfaces:
        raise ValueError("No valid building footprints available for meshing.")

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

    builder_mesh = _dtcc_builder.build_city_surface_mesh(
        building_surfaces,
        hole_surfaces,
        building_lod_switches,
        building_resolution,
        builder_dem,
        max_mesh_size,
        min_mesh_angle,
        smoothing,
        merge_meshes,
        sort_triangles,
    )
    
    if merge_meshes:
        result_mesh = builder_mesh[0].from_cpp()
    else:
        result_mesh = [bm.from_cpp() for bm in builder_mesh]
    return result_mesh


def build_city_volume_mesh(
    city: City,
    lod: GeometryType = GeometryType.LOD1,
    domain_height: float = 100.0,
    max_mesh_size: float = 10.0,
    merge_buildings: bool = True,
    boundary_face_markers: bool = True,
    tetgen_switches: Optional[Dict[str, Any]] = None,
    tetgen_switch_overrides: Optional[Dict[str, Any]] = None,
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
    tetgen_switches : dict, optional
        Optional TetGen switch parameters passed through to ``dtcc_wrapper_tetgen``.
        Provide keys as defined by ``dtcc_wrapper_tetgen.switches.DEFAULT_TETGEN_PARAMS``.
    tetgen_switch_overrides : dict, optional
        Optional low-level overrides forwarded to ``build_tetgen_switches`` for custom
        text-based switch assembly.

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
    min_mesh_angle = 30.0

    # Fallback dtcc volume meshing parameters
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

    _surfaces = [
        create_builder_surface(footprint)
        for footprint in building_footprints
        if footprint is not None
    ]
    hole_surfaces: list = []
    lod_switch_value = _LOD_PRIORITY.get(lod, _LOD_PRIORITY[GeometryType.LOD3])
    meshing_directives = [lod_switch_value] * len(_surfaces)

    # Convert from C++ to Python
    # ground_mesh = builder_mesh_to_mesh(_ground_mesh)
    _dem = raster_to_builder_gridfield(terrain.raster)

    if is_tetgen_available():

        #FIXME: Where do we set these parameters?
        smoothing = 1
        merge_meshes = True
        sort_triangles = False
        max_edge_radius_ratio = None  # 1.414
        min_dihedral_angle = None  # 30.0
        max_tet_volume = 20.0

        builder_mesh = _dtcc_builder.build_city_surface_mesh(
            _surfaces,
            hole_surfaces,
            meshing_directives,
            subdomain_resolution,
            _dem,
            max_mesh_size,
            min_mesh_angle,
            smoothing,
            merge_meshes,
            sort_triangles,
        )

        surface_mesh = builder_mesh[0].from_cpp()

        if surface_mesh.faces is None or len(surface_mesh.faces) == 0:
            raise ValueError("Surface mesh has no faces. Cannot build volume mesh.")
        if surface_mesh.markers is None or len(surface_mesh.markers) == 0:
            raise ValueError("Surface mesh has no face markers. Cannot build volume mesh.")
        
        switches_params = get_default_tetgen_switches()
        if max_tet_volume is not None:
            switches_params["max_volume"] = max_tet_volume
        if max_edge_radius_ratio is not None or min_dihedral_angle is not None:
            switches_params["quality"] = (
                max_edge_radius_ratio,
                min_dihedral_angle,
            )
        if tetgen_switches:
            switches_params.update(tetgen_switches)

        info("Building volume mesh with TetGen...")
        volume_mesh = tetgen_build_volume_mesh(
            mesh=surface_mesh,
            build_top_sidewalls=True,
            top_height=domain_height,
            switches_params=switches_params,
            switches_overrides=tetgen_switch_overrides,
            return_boundary_faces=boundary_face_markers, # Boundary face markers not implemented but returning boundary faces for now
        )
        return volume_mesh
        
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
        [],
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
    volume_mesh = _volume_mesh.from_cpp()

    if boundary_face_markers:
        boundary_face_markers = _dtcc_builder.compute_boundary_face_markers(
            _volume_mesh
        )
        if boundary_face_markers is not None:
            volume_mesh.boundary_markers = boundary_face_markers

    return volume_mesh
