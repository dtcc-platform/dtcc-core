from typing import Any, Dict, Optional, List, cast
import numpy as np

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
    build_terrain_surface_mesh,
    build_terrain_raster,
)

from ..building.modify import (
    merge_building_footprints,
    simplify_building_footprints,
    fix_building_footprint_clearance,
    clean_building_footprints,
    get_footprint,
)

from ..meshing.convert import mesh_to_raster

from ..logging import debug, info, warning, error

from ..meshing.tetgen import (
    build_volume_mesh as tetgen_build_volume_mesh,
    get_default_tetgen_switches,
    is_tetgen_available,
)

from dtcc_core.common.progress import report_progress


def _preprocess_buildings(
    buildings: List[Building],
    lod: GeometryType,
    merge_buildings: bool,
    merge_tolerance: float,
    min_building_area: float,
    min_building_detail: float,
    max_mesh_size: float,
) -> tuple[list, list[Building], list[float]]:
    """Merge, simplify and validate building footprints for meshing.

    This is the shared preprocessing pipeline used by
    :func:`build_city_flat_mesh` and :func:`build_city_volume_mesh`.

    Parameters
    ----------
    buildings : list[Building]
        Raw buildings from the city.
    lod : GeometryType
        Level-of-Detail used to extract footprints when *merge_buildings* is
        False.
    merge_buildings : bool
        Whether to merge adjacent/overlapping footprints.
    merge_tolerance : float
        Distance tolerance for the merge pass.
    min_building_area : float
        Minimum area threshold — smaller buildings are dropped.
    min_building_detail : float
        Minimum geometric feature size to resolve.
    max_mesh_size : float
        Fallback subdomain resolution when a building has no usable height.

    Returns
    -------
    building_footprints : list
        Validated footprint geometries ready for meshing.
    processed_buildings : list[Building]
        The corresponding Building objects (same order/length).
    subdomain_resolution : list[float]
        Per-building mesh resolution (min of height, max_mesh_size).
    """
    if merge_buildings:
        info(f"Merging {len(buildings)} buildings...")

        step1_buildings = cast(
            List[Building],
            merge_building_footprints(
                buildings,
                lod=GeometryType.LOD0,
                max_distance=merge_tolerance,
                min_area=min_building_area,
                return_index_map=False,
            ),
        )

        step2_buildings = cast(
            List[Building],
            merge_building_footprints(
                step1_buildings,
                GeometryType.LOD0,
                max_distance=0.0,
                min_area=min_building_area,
                return_index_map=False,
            ),
        )

        simplified_footprints = cast(
            List[Building],
            simplify_building_footprints(
                step2_buildings,
                min_building_detail,
                lod=GeometryType.LOD0,
                return_index_map=False,
            ),
        )

        building_footprints = [
            b.get_footprint(GeometryType.LOD0) for b in simplified_footprints
        ]
        processed_buildings = simplified_footprints

        info(f"After merging: {len(building_footprints)} buildings.")
    else:
        building_footprints = [b.get_footprint(lod) for b in buildings]
        processed_buildings = buildings

    # Filter out None and invalid footprints
    valid_indices = []
    for i, fp in enumerate(building_footprints):
        if fp is None:
            continue
        try:
            if hasattr(fp, "is_valid") and not fp.is_valid():
                warning(f"Skipping invalid footprint at index {i}")
                continue
            if hasattr(fp, "is_empty") and fp.is_empty():
                warning(f"Skipping empty footprint at index {i}")
                continue
            if hasattr(fp, "area") and fp.area <= 0:
                warning(f"Skipping zero-area footprint at index {i}")
                continue
            valid_indices.append(i)
        except Exception as e:
            warning(f"Skipping footprint at index {i} due to validation error: {e}")
            continue

    if not valid_indices:
        raise ValueError("No valid building footprints available for meshing.")

    building_footprints = [building_footprints[i] for i in valid_indices]
    processed_buildings = [processed_buildings[i] for i in valid_indices]

    info(f"Using {len(valid_indices)} valid building footprints for meshing.")

    # Set subdomain resolution based on building heights
    subdomain_resolution = []
    for building in processed_buildings:
        try:
            height = building.height
            if height is None or height <= 0:
                height = max_mesh_size
            subdomain_resolution.append(min(height, max_mesh_size))
        except (AttributeError, TypeError):
            subdomain_resolution.append(max_mesh_size)

    return building_footprints, processed_buildings, subdomain_resolution


def build_city_surface_mesh(
    city: City,
    lod: GeometryType | list[GeometryType] = GeometryType.LOD1,
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
    report_mesh_quality: bool = True,
) -> Mesh:
    """
    Build a surface mesh from the surfaces of the buildings in the city.

    Parameters
    ----------
    `city` : model.City
        The city to build the mesh from.
    `lod` : GeometryType or list of GeometryType, optional
        The meshing directive (Level of Detail) to apply to the buildings.
        If a single value is provided, it is applied uniformly to all buildings.
        If a list is provided, it must have the same length as the number of buildings
        in the city, and each entry specifies the directive for the corresponding building.
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
            raise ValueError(
                f"lod list length {len(lod)} != number of buildings {n_buildings}"
            )
        if not all(isinstance(x, GeometryType) for x in lod):
            raise TypeError("all elements in lod list must be GeometryType instances")
        lod = list(lod)
    else:
        raise TypeError(
            f"lod must be a single GeometryType or a list/tuple of {n_buildings} GeometryType values, "
            f"got {type(lod).__name__}"
        )

    current_index_map: list[list[int]] | None = None

    if merge_buildings:
        info(f"Merging {len(buildings)} buildings...")
        merged_buildings, index_map = merge_building_footprints(
            buildings,
            lod=GeometryType.LOD0,
            max_distance=merge_tolerance,
            min_area=min_building_area,
            return_index_map=True,
        )
        current_index_map = index_map

        # city.replace_buildings(merged_buildings)
        # city.save_building_footprints("footprints_merged_mesher.gpkg")

        smallest_hole = max(min_building_detail, min_building_detail**2)
        cleaned_footprints, cleaned_index_map = clean_building_footprints(
            merged_buildings,
            clearance=min_building_detail,
            smallest_hole_area=smallest_hole,
            return_index_map=True,
        )
        current_index_map = compose_index_map(current_index_map, cleaned_index_map)

        # city.replace_buildings(cleaned_footprints)
        # city.save_building_footprints("footprints_merged_cleaned_mesher.gpkg")

        merged_buildings, merged_index_map = merge_building_footprints(
            cleaned_footprints,
            GeometryType.LOD0,
            max_distance=0.0,
            min_area=min_building_area,
            return_index_map=True,
        )
        current_index_map = compose_index_map(current_index_map, merged_index_map)

        # city.replace_buildings(merged_buildings)
        # city.save_building_footprints("footprints_merged_cleaned_merged_mesher.gpkg")

        simplifed_footprints, simplified_index_map = simplify_building_footprints(
            merged_buildings,
            min_building_detail / 2,
            lod=GeometryType.LOD0,
            return_index_map=True,
        )
        current_index_map = compose_index_map(current_index_map, simplified_index_map)

        target_lods = (
            lod_from_index_map(current_index_map)
            if current_index_map is not None
            else []
        )

        building_footprints = [
            b.get_footprint(GeometryType.LOD0) for b in simplifed_footprints
        ]

    else:
        target_lods = lod

        building_footprints = [
            b.get_footprint() for b, b_lod in zip(buildings, target_lods)
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

    if report_mesh_quality:
        from dtcc_core.model.mixins.mesh.quality import triangle_mesh_quality, report_quality
        if merge_meshes:
            q = triangle_mesh_quality(result_mesh.vertices, result_mesh.faces)
            report_quality(q, log_fn=info)
        else:
            for i, m in enumerate(result_mesh):
                q = triangle_mesh_quality(m.vertices, m.faces)
                report_quality(q, log_fn=info)

    return result_mesh


def build_city_flat_mesh(
    city: City,
    lod: GeometryType = GeometryType.LOD1,
    max_mesh_size: float = 10.0,
    min_mesh_angle: float = 25.0,
    merge_buildings: bool = True,
    min_building_detail: float = 0.5,
    min_building_area: float = 15.0,
    merge_tolerance: float = 0.5,
    report_mesh_quality: bool = True,
) -> Mesh:
    """Build a flat 2D triangular mesh of the city with building footprints marked.

    The mesh lies in the z = 0 plane. Triangle edges conform to building
    footprint boundaries and each triangle carries an integer marker:

    * ``-2`` — ground (outside buildings and halos)
    * ``-1`` — halo (triangles that touch a building but are not inside one)
    * ``0, 1, 2, …`` — index of the building whose footprint contains the
      triangle

    Parameters
    ----------
    city : City
        City object containing terrain bounds and building data.
    lod : GeometryType, optional
        Level-of-Detail used when *merge_buildings* is False (default LOD1).
    max_mesh_size : float, optional
        Maximum triangle size (default 10.0).
    min_mesh_angle : float, optional
        Minimum angle quality constraint in degrees (default 25.0).
    merge_buildings : bool, optional
        Merge adjacent/overlapping building footprints (default True).
    min_building_detail : float, optional
        Minimum feature size to resolve in footprints (default 0.5).
    min_building_area : float, optional
        Minimum footprint area; smaller buildings are dropped (default 15.0).
    merge_tolerance : float, optional
        Distance tolerance for merging footprints (default 0.5).

    Returns
    -------
    Mesh
        A flat (z = 0) triangular mesh with per-face markers indicating
        building membership.

    Raises
    ------
    ValueError
        If the city has no terrain data or no valid building footprints.
    """
    # Validate terrain (needed for domain bounds)
    terrain = city.terrain
    if terrain is None:
        raise ValueError("City has no terrain data. Please compute terrain first.")

    buildings = city.buildings
    if not buildings:
        warning("City has no buildings.")

    # Preprocess buildings
    building_footprints, _processed_buildings, subdomain_resolution = (
        _preprocess_buildings(
            buildings,
            lod=lod,
            merge_buildings=merge_buildings,
            merge_tolerance=merge_tolerance,
            min_building_area=min_building_area,
            min_building_detail=min_building_detail,
            max_mesh_size=max_mesh_size,
        )
    )

    report_progress(
        percent=10,
        message=f"Preprocessed {len(building_footprints)} building footprints",
    )

    # Convert footprints to builder polygons
    _building_polygons = [
        create_builder_polygon(footprint.to_polygon())
        for footprint in building_footprints
    ]

    # Call C++ mesher
    report_progress(percent=30, message="Building city flat mesh (C++)...")
    _flat_mesh = _dtcc_builder.build_city_flat_mesh(
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

    flat_mesh = _flat_mesh.from_cpp()

    if report_mesh_quality:
        from dtcc_core.model.mixins.mesh.quality import triangle_mesh_quality, report_quality
        q = triangle_mesh_quality(flat_mesh.vertices, flat_mesh.faces)
        report_quality(q, log_fn=info)

    report_progress(percent=100, message="City flat mesh complete")
    return flat_mesh


def build_city_volume_mesh(
    city: City,
    lod: GeometryType = GeometryType.LOD1,
    domain_height: float = 100.0,
    max_mesh_size: float = 10.0,
    min_mesh_angle: float = 25.0,
    merge_buildings: bool = True,
    min_building_detail: float = 0.5,
    min_building_area: float = 15.0,
    merge_tolerance: float = 0.5,
    smoothing: int = 0,
    boundary_face_markers: bool = True,
    tetgen_switches: Optional[Dict[str, Any]] = None,
    tetgen_switch_overrides: Optional[Dict[str, Any]] = None,
    # Fallback DTCC volume mesher parameters
    smoother_max_iterations: int = 5000,
    smoothing_relative_tolerance: float = 0.005,
    aspect_ratio_threshold: float = 10.0,
    debug_step: int = 7,
    report_mesh_quality: bool = True,
) -> VolumeMesh:
    """
    Build a 3D tetrahedral volume mesh for a city terrain with embedded building volumes.

    This function generates a ground mesh from the city terrain and extrudes building
    footprints to produce a full volume mesh, optionally merging adjacent buildings
    and marking boundary faces.

    Parameters
    ----------
    city : City
        City object containing terrain and building data. The terrain must provide
        either a raster or a mesh representation to support domain surface generation.
    lod : GeometryType, optional
        The meshing directive (Level of Detail) applied to building footprints.
        Defaults to ``GeometryType.LOD1``.
    domain_height : float, optional
        The vertical height of the volume domain above the terrain surface, in the
        same coordinate units as the city. Defaults to 100.0.
    max_mesh_size : float, optional
        Maximum allowed mesh element size. This value governs both the underlying
        ground mesh resolution and the upper bound of element sizing within the
        extruded volume. Defaults to 10.0.
    min_mesh_angle : float, optional
        Minimum allowable mesh angle used as a quality constraint. Defaults to 25.0.
    merge_buildings : bool, optional
        Whether to merge adjacent or overlapping building footprints into larger
        composite blocks prior to meshing. Defaults to True.
    min_building_detail : float, optional
        Minimum geometric feature size to resolve within building footprints.
        Defaults to 0.5.
    min_building_area : float, optional
        Minimum footprint area required for a building to be included in the mesh.
        Buildings below this threshold are omitted. Defaults to 15.0.
    merge_tolerance : float, optional
        Distance tolerance used when merging building footprints. Defaults to 0.5.
    smoothing : int, optional
        Number of mesh-smoothing iterations applied to the terrain and building
        surface meshes prior to volume meshing. Defaults to 0 (no smoothing).
    boundary_face_markers : bool, optional
        If True, annotate boundary faces of the resulting volume mesh with integer
        markers as a post-processing step. This supports downstream workflows such
        as boundary-condition assignment. Defaults to True. See Notes for marker
        conventions.
    tetgen_switches : dict, optional
        Optional high-level TetGen parameter dictionary. Keys must correspond to
        those defined in ``dtcc_wrapper_tetgen.switches.DEFAULT_TETGEN_PARAMS``.
        These values are passed directly to ``dtcc_wrapper_tetgen``.
    tetgen_switch_overrides : dict, optional
        Optional low-level overrides for custom text-based switch assembly via
        ``build_tetgen_switches``. Use this when direct control of TetGen's
        command-string switches is required.
    smoother_max_iterations : int, optional
        Maximum iterations for fallback volume mesh smoother. Defaults to 5000.
    smoothing_relative_tolerance : float, optional
        Relative tolerance for fallback volume mesh smoothing. Defaults to 0.005.
    aspect_ratio_threshold : float, optional
        Aspect ratio threshold for fallback volume mesher. Defaults to 10.0.
    debug_step : int, optional
        Debug step parameter for fallback volume mesher. Defaults to 7.

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
    ValueError
        If no valid building footprints are available after preprocessing.

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
    - TetGen is preferred when available. If not available, falls back to the
      internal DTCC volume mesh builder.

    Examples
    --------
    >>> mesh = build_city_volume_mesh(my_city,
    ...                               lod=GeometryType.LOD1,
    ...                               domain_height=150.0,
    ...                               max_mesh_size=5.0,
    ...                               merge_buildings=False,
    ...                               boundary_face_markers=True)
    """

    # 1. VALIDATE INPUT AND TERRAIN

    buildings = city.buildings
    if not buildings:
        warning("City has no buildings.")

    terrain = city.terrain
    if terrain is None:
        raise ValueError("City has no terrain data. Please compute terrain first.")
    terrain_raster = terrain.raster
    terrain_mesh = terrain.mesh
    if terrain_raster is None and terrain_mesh is None:
        raise ValueError("City terrain has no data. Please compute terrain first.")
    if terrain_raster is None and terrain_mesh is not None:
        terrain_raster = mesh_to_raster(terrain_mesh, cell_size=max_mesh_size)

    # 2. PREPROCESS BUILDINGS
    building_footprints, processed_buildings, subdomain_resolution = (
        _preprocess_buildings(
            buildings,
            lod=lod,
            merge_buildings=merge_buildings,
            merge_tolerance=merge_tolerance,
            min_building_area=min_building_area,
            min_building_detail=min_building_detail,
            max_mesh_size=max_mesh_size,
        )
    )

    report_progress(
        percent=10,
        message=f"Preprocessed {len(building_footprints)} building footprints",
    )

    # 3. prepare builder objects

    _surfaces = [create_builder_surface(footprint) for footprint in building_footprints]
    hole_surfaces: list = []
    lod_switch_value = _LOD_PRIORITY.get(lod, _LOD_PRIORITY[GeometryType.LOD3])
    meshing_directives = [lod_switch_value] * len(_surfaces)
    _dem = raster_to_builder_gridfield(terrain_raster)

    # 4. BUILD VOLUME MESH - TETGEN PATH

    if is_tetgen_available():
        info("Building volume mesh with TetGen...")
        report_progress(percent=30, message="Preparing builder objects...")

        # Validate inputs before calling C++ mesher
        info(f"Number of surfaces: {len(_surfaces)}")
        info(f"Number of subdomain resolutions: {len(subdomain_resolution)}")
        info(
            f"Max mesh size: {max_mesh_size}, Min angle: {min_mesh_angle}, Smoothing: {smoothing}"
        )

        if len(_surfaces) != len(subdomain_resolution):
            raise ValueError(
                f"Mismatch: {len(_surfaces)} surfaces but {len(subdomain_resolution)} resolution values"
            )

        # Build surface mesh
        merge_meshes = True
        sort_triangles = False
        report_progress(percent=40, message="Building surface mesh (C++)...")

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
        report_progress(
            percent=55, message="Surface mesh built, preparing volume mesh..."
        )

        # Validate surface mesh
        if surface_mesh.faces is None or len(surface_mesh.faces) == 0:
            raise ValueError("Surface mesh has no faces. Cannot build volume mesh.")
        if surface_mesh.markers is None or len(surface_mesh.markers) == 0:
            raise ValueError(
                "Surface mesh has no face markers. Cannot build volume mesh."
            )

        # Configure TetGen switches
        switches_params = get_default_tetgen_switches()
        if tetgen_switches:
            switches_params.update(tetgen_switches)

        # Build volume mesh with TetGen
        report_progress(percent=60, message="Running TetGen volume mesher...")
        volume_mesh = tetgen_build_volume_mesh(
            mesh=surface_mesh,
            build_top_sidewalls=True,
            top_height=domain_height,
            switches_params=switches_params,
            switches_overrides=tetgen_switch_overrides,
            return_boundary_faces=boundary_face_markers,
        )
        report_progress(percent=95, message="Volume mesh complete")

        if report_mesh_quality:
            from dtcc_core.model.mixins.mesh.quality import tetrahedron_mesh_quality, report_quality
            q = tetrahedron_mesh_quality(volume_mesh.vertices, volume_mesh.cells)
            report_quality(q, log_fn=info)

        return volume_mesh

    # 5. BUILD VOLUME MESH - FALLBACK DTCC PATH
    info("Building volume mesh with fallback DTCC volume mesher...")
    report_progress(percent=40, message="Building volume mesh (fallback mesher)...")

    # Convert footprints to builder polygons for ground mesh
    _building_polygons = [
        create_builder_polygon(footprint.to_polygon())
        for footprint in building_footprints
    ]

    # Build flat mesh (ground mesh with building markers)
    _ground_mesh = _dtcc_builder.build_city_flat_mesh(
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

    # Create volume mesh builder
    volume_mesh_builder = _dtcc_builder.VolumeMeshBuilder(
        _surfaces, _dem, _ground_mesh, domain_height
    )

    # Build volume mesh
    _volume_mesh = volume_mesh_builder.build(
        smoother_max_iterations,
        smoothing_relative_tolerance,
        0.0,
        aspect_ratio_threshold,
        debug_step,
    )
    volume_mesh = _volume_mesh.from_cpp()
    report_progress(percent=90, message="Volume mesh built, finalizing...")

    # Add boundary face markers if requested
    if boundary_face_markers:
        computed_markers = _dtcc_builder.compute_boundary_face_markers(_volume_mesh)
        if computed_markers is not None:
            volume_mesh.boundary_markers = computed_markers

    if report_mesh_quality:
        from dtcc_core.model.mixins.mesh.quality import tetrahedron_mesh_quality, report_quality
        q = tetrahedron_mesh_quality(volume_mesh.vertices, volume_mesh.cells)
        report_quality(q, log_fn=info)

    return volume_mesh
