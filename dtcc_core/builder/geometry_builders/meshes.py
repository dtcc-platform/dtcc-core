from typing import Any, Dict, Optional, List, Sequence
import numpy as np

import dtcc_core.io
from ...model import (
    Mesh,
    VolumeMesh,
    Building,
    City,
    Surface,
    GeometryType,
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
    raster_to_builder_gridfield,
)

from .. import _dtcc_builder

from ..cleaning import (
    ConditioningOptions,
    condition_building_footprints,
    condition_polygon_coverage,
)

from ..meshing.convert import mesh_to_raster

from ..logging import debug, info, warning, error

from ..meshing.tetgen import (
    build_volume_mesh as tetgen_build_volume_mesh,
    get_default_tetgen_switches,
    is_tetgen_available,
)

from dtcc_core.common.progress import report_progress


def _normalize_lod_values(
    buildings: list[Building],
    lod: GeometryType | Sequence[GeometryType],
) -> list[GeometryType]:
    if isinstance(lod, GeometryType):
        return [lod] * len(buildings)
    if len(lod) != len(buildings):
        raise ValueError(
            f"lod list length {len(lod)} != number of buildings {len(buildings)}"
        )
    if not all(isinstance(value, GeometryType) for value in lod):
        raise TypeError("all elements in lod list must be GeometryType instances")
    return list(lod)


def _extract_meshing_polygon(
    building: Building,
    lod: GeometryType,
):
    geometry = building.flatten_geometry(lod)
    if geometry is None:
        return None, 0.0, None, 0.0

    polygon = geometry.to_polygon(simplify=0.0)
    if polygon is None or polygon.is_empty:
        return None, 0.0, None, 0.0

    roof_z = 0.0
    try:
        roof_z = float(geometry.bounds.zmax)
    except (AttributeError, TypeError):
        roof_z = 0.0

    try:
        height = float(building.height)
    except (AttributeError, TypeError):
        height = None

    if height is not None and height <= 0:
        height = None

    return polygon, float(max(polygon.area, 0.0)), height, roof_z


def _area_weighted_value(
    source_indices: Sequence[int],
    source_areas: Sequence[float],
    values: Sequence[float | None],
    *,
    default: float,
) -> float:
    weighted_values: list[float] = []
    weights: list[float] = []
    for index in source_indices:
        if index < 0 or index >= len(values):
            continue
        value = values[index]
        if value is None:
            continue
        weight = source_areas[index] if source_areas[index] > 0 else 1.0
        weighted_values.append(float(value))
        weights.append(float(weight))
    if not weighted_values:
        return default
    return float(np.average(weighted_values, weights=weights))


def _condition_meshing_footprints(
    buildings: list[Building],
    *,
    lod: GeometryType | list[GeometryType],
    min_building_detail: float,
    min_building_area: float,
    merge_tolerance: float,
    merge_buildings: bool,
    max_mesh_size: float,
) -> tuple[list[Surface], list[list[int]], list[float], dict[str, Any]]:
    if not buildings:
        warning("No buildings to preprocess.")
        return [], [], [], {}

    lod_values = _normalize_lod_values(buildings, lod)
    extracted_polygons = []
    initial_source_map: list[list[int]] = []
    source_areas = [0.0] * len(buildings)
    source_heights: list[float | None] = [None] * len(buildings)
    source_roof_z = [0.0] * len(buildings)

    for index, (building, lod_value) in enumerate(zip(buildings, lod_values)):
        polygon, area, height, roof_z = _extract_meshing_polygon(building, lod_value)
        source_areas[index] = area
        source_heights[index] = height
        source_roof_z[index] = roof_z
        if polygon is None:
            continue
        extracted_polygons.append(polygon)
        initial_source_map.append([index])

    options = ConditioningOptions(
        precision_grid=None,
        min_feature_size=min_building_detail,
        merge_distance=merge_tolerance if merge_buildings else 0.0,
        min_area=min_building_area,
        min_hole_area=min_building_detail**2,
    )

    if isinstance(lod, GeometryType):
        result = condition_building_footprints(
            buildings,
            lod=lod,
            options=options,
        )
    else:
        result = condition_polygon_coverage(
            extracted_polygons,
            source_map=initial_source_map,
            options=options,
        )

    conditioned_surfaces: list[Surface] = []
    subdomain_resolution: list[float] = []

    for polygon, source_indices in zip(result.polygons, result.source_map):
        roof_z = _area_weighted_value(
            source_indices,
            source_areas,
            source_roof_z,
            default=0.0,
        )
        height = _area_weighted_value(
            source_indices,
            source_areas,
            source_heights,
            default=max_mesh_size,
        )

        surface = Surface()
        surface.from_polygon(polygon, roof_z)
        conditioned_surfaces.append(surface)
        subdomain_resolution.append(min(height, max_mesh_size))

    info(
        f"Conditioned {len(buildings)} buildings into {len(conditioned_surfaces)} meshing footprints."
    )
    return (
        conditioned_surfaces,
        result.source_map,
        subdomain_resolution,
        result.diagnostics,
    )


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

    buildings = city.buildings
    lod_values = _normalize_lod_values(buildings, lod)
    building_footprints, source_map, conditioned_resolution, conditioning_diagnostics = (
        _condition_meshing_footprints(
            buildings,
            lod=lod if isinstance(lod, GeometryType) else lod_values,
            min_building_detail=min_building_detail,
            min_building_area=min_building_area,
            merge_tolerance=merge_tolerance,
            merge_buildings=merge_buildings,
            max_mesh_size=max_mesh_size,
        )
    )
    target_lods = [
        min((lod_values[index] for index in indices), key=lambda value: _LOD_PRIORITY[value])
        for indices in source_map
    ]
    base_resolution = [
        min(resolution, building_mesh_triangle_size) for resolution in conditioned_resolution
    ]
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
    debug(f"Surface meshing footprint diagnostics: {conditioning_diagnostics}")

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
        from dtcc_core.model.mixins.mesh.quality import (
            triangle_mesh_quality,
            report_quality,
        )

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

    building_footprints, _source_map, subdomain_resolution, diagnostics = (
        _condition_meshing_footprints(
            buildings,
            lod=lod,
            min_building_detail=min_building_detail,
            min_building_area=min_building_area,
            merge_tolerance=merge_tolerance,
            merge_buildings=merge_buildings,
            max_mesh_size=max_mesh_size,
        )
    )

    if not building_footprints:
        raise ValueError("No valid building footprints available for meshing.")

    report_progress(
        percent=10,
        message=f"Preprocessed {len(building_footprints)} building footprints",
    )
    debug(f"Flat meshing footprint diagnostics: {diagnostics}")

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
        from dtcc_core.model.mixins.mesh.quality import (
            triangle_mesh_quality,
            report_quality,
        )

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
    building_footprints, source_map, subdomain_resolution, diagnostics = (
        _condition_meshing_footprints(
            buildings,
            lod=lod,
            min_building_detail=min_building_detail,
            min_building_area=min_building_area,
            merge_tolerance=merge_tolerance,
            merge_buildings=merge_buildings,
            max_mesh_size=max_mesh_size,
        )
    )
    if not building_footprints:
        raise ValueError("No valid building footprints available for meshing.")

    report_progress(
        percent=10,
        message=f"Preprocessed {len(building_footprints)} building footprints",
    )
    debug(f"Volume meshing footprint diagnostics: {diagnostics}")

    # 3. prepare builder objects

    _surfaces = [create_builder_surface(footprint) for footprint in building_footprints]
    hole_surfaces: list = []
    meshing_directives = [
        _LOD_PRIORITY.get(lod, _LOD_PRIORITY[GeometryType.LOD3])
    ] * len(_surfaces)
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
        try:
            volume_mesh = tetgen_build_volume_mesh(
                mesh=surface_mesh,
                build_top_sidewalls=True,
                top_height=domain_height,
                switches_params=switches_params,
                switches_overrides=tetgen_switch_overrides,
                return_boundary_faces=boundary_face_markers,
            )
        except RuntimeError as exc:
            msg = str(exc)
            # Merged footprints can occasionally create degenerate/overlapping
            # facets that TetGen reports as self-intersections. Retry once
            # without building-merge preprocessing for robustness.
            if merge_buildings and "self-intersections" in msg:
                warning(
                    "TetGen failed with self-intersections after merging buildings; "
                    "retrying once with merge_buildings=False."
                )
                return build_city_volume_mesh(
                    city=city,
                    lod=lod,
                    domain_height=domain_height,
                    max_mesh_size=max_mesh_size,
                    min_mesh_angle=min_mesh_angle,
                    merge_buildings=False,
                    min_building_detail=min_building_detail,
                    min_building_area=min_building_area,
                    merge_tolerance=merge_tolerance,
                    smoothing=smoothing,
                    boundary_face_markers=boundary_face_markers,
                    tetgen_switches=tetgen_switches,
                    tetgen_switch_overrides=tetgen_switch_overrides,
                    smoother_max_iterations=smoother_max_iterations,
                    smoothing_relative_tolerance=smoothing_relative_tolerance,
                    aspect_ratio_threshold=aspect_ratio_threshold,
                    debug_step=debug_step,
                    report_mesh_quality=report_mesh_quality,
                )
            raise
        report_progress(percent=95, message="Volume mesh complete")

        if report_mesh_quality:
            from dtcc_core.model.mixins.mesh.quality import (
                tetrahedron_mesh_quality,
                report_quality,
            )

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
        from dtcc_core.model.mixins.mesh.quality import (
            tetrahedron_mesh_quality,
            report_quality,
        )

        q = tetrahedron_mesh_quality(volume_mesh.vertices, volume_mesh.cells)
        report_quality(q, log_fn=info)

    return volume_mesh
