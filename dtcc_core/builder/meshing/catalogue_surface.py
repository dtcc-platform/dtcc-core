"""Opt-in city orchestration for fitted planar roofs and conforming shells."""
import numpy as np
from shapely.geometry import box
from shapely.ops import unary_union

from ...model import Field, Surface
from ..cleaning.clean import clean_regions, cleaning_record, regions_from_cleaned_city
from ..roofs import RoofFitOptions
from .roof_coverage import local_roof_coverage
from .roof_regions import source_roofs
from .roof_interfaces import RoofAssemblyError, node_interfaces
from .roof_validation import validate_partition, validate_triangulation, validate_shell


def build_catalogue_surface(city, *, lod, min_building_detail, min_building_area,
                            building_mesh_triangle_size, max_mesh_size, min_mesh_angle,
                            merge_meshes, smoothing, sort_triangles, report_mesh_quality,
                            show_footprints, footprint_cleaning_plot_block, mesher,
                            pipeline_mode, stage_audit, roof_options, surface_refinement=None, compact_surface=None,
                            lod2_method='catalogue'):
    # Import here to keep the common orchestration helpers independent of fitting.
    """Fit source roofs and assemble a validated joint terrain/building surface."""
    from ..geometry_builders import meshes as pipeline

    options = RoofFitOptions() if roof_options is None else roof_options
    audit = {} if stage_audit is None else stage_audit
    catalogue = dict(sources=[], blocks=[], assemblies=[], block_fallbacks=[],
                     volume_meshing_supported=bool(merge_meshes))
    audit.update(requested_lod='lod2', lod2_method=lod2_method)
    audit['roof_catalogue'] = catalogue
    attempt = pipeline._start_stage_audit_attempt(
        audit, label='catalogue', retry_reason=None, backend='surface',
        requested_mesher=mesher, max_mesh_size=max_mesh_size,
        min_mesh_angle=min_mesh_angle, domain_height=0., tetgen_switches=None,
        tetgen_switch_overrides=None)
    try:
        max_mesh_size = pipeline._normalize_max_mesh_size(max_mesh_size)
        _, raster = pipeline._require_city_terrain_raster(city, max_mesh_size=max_mesh_size)
        bounds = pipeline._raster_bounds_tuple(raster)
        blocks, rows = _prepare_roof_blocks(
            city, bounds=bounds, lod=lod, options=options, catalogue=catalogue,
            min_building_detail=min_building_detail, min_building_area=min_building_area,
            show_footprints=show_footprints, footprint_cleaning_plot_block=footprint_cleaning_plot_block,
            lod2_method=lod2_method)
        result = _mesh_roof_blocks(
            blocks, rows, raster=raster, bounds=bounds, catalogue=catalogue,
            min_building_detail=min_building_detail, building_mesh_triangle_size=building_mesh_triangle_size,
            max_mesh_size=max_mesh_size, min_mesh_angle=min_mesh_angle, merge_meshes=merge_meshes,
            smoothing=smoothing, sort_triangles=sort_triangles, mesher=mesher,
            compact_surface=compact_surface, surface_refinement=surface_refinement, lod2_method=lod2_method)
        for mesh in [result] if merge_meshes else result:
            if not merge_meshes:
                mesh.fields.append(Field(name='volume_boundary', values=np.array([0.]), association='geometry'))
            mesh.stage_audit = audit
        if report_mesh_quality:
            from ...model.mixins.mesh.quality import triangle_mesh_quality, report_quality
            for mesh in [result] if merge_meshes else result:
                report_quality(triangle_mesh_quality(mesh.vertices, mesh.faces), log_fn=pipeline.info)
        pipeline._mark_stage_audit_success(audit, attempt)
        return result
    except Exception as exc:
        pipeline._mark_stage_audit_failure(attempt, exc)
        raise


def _prepare_roof_blocks(city, *, bounds, lod, options, catalogue, min_building_detail,
                         min_building_area, show_footprints, footprint_cleaning_plot_block,
                         lod2_method, isolated=False):
    """Clean once and fit original source observations for either surface layout."""
    from ..geometry_builders import meshes as pipeline

    domain = box(*bounds)
    record = cleaning_record(city, min_building_detail=min_building_detail, lod=lod, isolated=isolated)
    if record is None:
        cleaned = clean_regions(city.buildings, min_building_detail=min_building_detail,
                                lod=lod, isolated=isolated, show=show_footprints,
                                plot_block=footprint_cleaning_plot_block)
    else:
        cleaned = regions_from_cleaned_city(city, record)
    catalogue['reused_cleaning'] = record is not None
    selected_by_output = {}
    for polygon, output in zip(cleaned.polygons, cleaned.region_outputs):
        if polygon.area < min_building_area:
            continue
        if isolated and not domain.contains(polygon):
            continue
        clipped = polygon if isolated else polygon.intersection(domain)
        if not clipped.is_empty and clipped.area > 0:
            selected_by_output.setdefault(output, []).append(clipped)
    selected = [(unary_union(parts), output) for output, parts in selected_by_output.items()]
    blocks, rows = source_roofs(city, cleaned, selected, options=options,
                               fallback_regions=lambda output, polygon: pipeline._building_level_regions(
                                   cleaned, output, polygon, delta=float(min_building_detail)),
                               reused=record is not None, method=lod2_method)
    catalogue['sources'] = rows
    catalogue['blocks'] = [dict(block_id=b.output, source_keys=b.source_keys,
                                fallback_levels=sorted(set(float(p.coefficients[2]) for p in b.lod1_patches))) for b in blocks]
    return blocks, rows


def _mesh_roof_blocks(blocks, rows, *, raster, bounds, catalogue, min_building_detail,
                      building_mesh_triangle_size, max_mesh_size, min_mesh_angle, merge_meshes,
                      smoothing, sort_triangles, mesher, compact_surface, surface_refinement,
                      lod2_method, base_z=None):
    """Assemble and validate fitted blocks against a terrain surface."""
    from ..geometry_builders import meshes as pipeline

    target = unary_union([b.polygon for b in blocks])
    fallback_blocks = set()
    # Local partition failures have exact block attribution. Each block
    # can be flattened once; a failing flat replacement is an error.
    prepared = {}
    for block in blocks:
        try:
            validate_partition(block.patches, block.polygon)
            prepared[block.output] = node_interfaces(block.patches)
            validate_partition(prepared[block.output], block.polygon)
            # Admit every block at the same coverage precision before the
            # two city assembly attempts; collect all local failures once.
            local_roof_coverage(prepared[block.output], bounds)
        except (RoofAssemblyError, ValueError) as exc:
            fallback_blocks.add(block.output)
            catalogue['block_fallbacks'].append(dict(block_id=block.output,
                source_keys=block.source_keys,
                fitted_shapes=[r['shape'] for r in rows if r.get('source_key') in block.source_keys],
                rejection_stage='partition', reason=str(exc)))
    result = None
    for assembly in range(2):
        row = dict(attempt=assembly + 1, flat_blocks=sorted(fallback_blocks), stage='partition')
        catalogue['assemblies'].append(row)
        try:
            patches = [p for block in blocks for p in
                       (block.flat_patches() if block.output in fallback_blocks else prepared[block.output])]
            row['partition'] = validate_partition(patches, target)
            patches = node_interfaces(patches)
            row['noded_partition'] = validate_partition(patches, target)
            if base_z is not None and compact_surface is not None:
                from .compact_surface import compact_surface as build_compact
                row['stage'] = 'skin_construction'
                result, wall_audit = build_compact(None, patches, compact_surface,
                    base_z=base_z, allow_source_steps=lod2_method == 'ransac')
                row['shell'] = wall_audit['topology']
                wall_audit['shell_accepted'] = True
            else:
                local_patches, patches, ground, local_bounds, offset = local_roof_coverage(patches, bounds)
                compact_patches = list(patches)
                regions = [*ground, *[p.polygon for p in local_patches]]
                size = min(max_mesh_size or float('inf'), building_mesh_triangle_size) if building_mesh_triangle_size > 0 else max_mesh_size
                sizes = {i: size for i in range(len(patches))} if size is not None else None
                # Patches are already a validated subdivision. Never pass them
                # through footprint cleaning or minimum-building-area selection.
                row['stage'] = 'triangulation'
                ground_mesh, backend = pipeline._build_ground_mesh_from_coverage(
                    region_polygons=regions,
                    region_markers=[-2]*len(ground) + list(range(len(patches))),
                    region_points=[np.array(p.representative_point().coords[0]) for p in regions],
                    bounds=local_bounds, max_mesh_size=max_mesh_size, min_mesh_angle=min_mesh_angle,
                    mesher=mesher, sort_triangles=sort_triangles,
                    region_triangle_sizes=sizes, add_halo_markers=False, coordinate_tolerance=1e-8, max_steiner_points=50000)
                row['mesher'] = backend
                row['ground_mesh'] = pipeline._triangle_mesh_audit(ground_mesh)
                row['ground_contract'] = pipeline._triangle_mesh_contract_from_audit(
                    row['ground_mesh'], reference_length=float(min_building_detail),
                    require_markers=True, stage_label='Catalogue ground mesh')
                pipeline._raise_stage_contract_errors('Catalogue ground mesh', row['ground_contract'])
                row['max_steiner_points'] = 50000 if backend == 'triangle' else None
                validate_triangulation(ground_mesh, local_patches)
                ground_mesh.vertices[:, :2] += offset
                surfaces = []
                for patch in patches:
                    surface = Surface()
                    surface.from_polygon(patch.polygon, float(patch.coefficients[2]))
                    surfaces.append(surface)
                owners = {key: i for i, key in enumerate(dict.fromkeys(p.source_key for p in patches))}
                row['stage'] = 'shell_construction'
                result = pipeline._build_city_surface_mesh_from_ground_mesh(
                    ground_mesh=ground_mesh, terrain_raster=raster,
                    building_surfaces=surfaces, meshing_directives=[1]*len(surfaces),
                    smoothing=smoothing, merge_meshes=merge_meshes,
                    min_wall_strip=float(min_building_detail),
                    patch_planes=[p.native_plane() for p in patches],
                    # Native owners enforce roof continuity, not semantic ownership.
                    # RANSAC may have legitimate steps within one source building.
                    patch_owners=(list(range(len(patches))) if lod2_method == 'ransac' else
                                  [owners[p.source_key] for p in patches]))
                welded = result if merge_meshes else pipeline._welded_mesh(result)
                row['stage'] = 'shell_validation'
                row['shell'] = validate_shell(welded, num_patches=len(patches),
                                               closure_audit=pipeline._shell_closure_audit)
                row['shell']['mesh_audit'] = pipeline._triangle_mesh_audit(welded)
                row['shell']['quality_contract'] = pipeline._triangle_mesh_contract_from_audit(
                    row['shell']['mesh_audit'], reference_length=float(min_building_detail),
                    require_markers=True, stage_label='Catalogue surface shell')
                row['shell']['topology_valid'] = row['shell']['valid']
                row['shell']['quality_warning'] = bool(row['shell']['quality_contract']['warnings'])
                pipeline._raise_stage_contract_errors('Catalogue surface shell', row['shell']['quality_contract'])
                if not row['shell']['valid']:
                    raise RoofAssemblyError(f"City shell validation failed: {row['shell']}")
                # Quality refinement is a transaction after a valid shell exists.
                # Its rejection preserves that shell, never reconstruction fits.
                from .wall_faces import refine_wall_faces, separate_wall_mesh
                original_result, original_welded = result, welded
                wall_audit = {}
                try:
                    if compact_surface is not None:
                        from .compact_surface import compact_surface as build_compact
                        candidate, wall_audit = build_compact(welded, compact_patches, compact_surface,
                            **({'allow_source_steps': True} if lod2_method == 'ransac' else {}))
                    elif surface_refinement is None:
                        candidate, wall_audit = refine_wall_faces(
                            welded, [p.source_key for p in patches])
                    else:
                        from .surface_refinement import refine_shared_surface
                        candidate, wall_audit = refine_shared_surface(
                            welded, [p.source_key for p in patches], surface_refinement,
                            building_spacing=size or 4., terrain_spacing=max_mesh_size or 8.)
                        if not wall_audit['accepted']:
                            raise ValueError(wall_audit['error'])
                    candidate_shell = validate_shell(candidate, num_patches=len(patches),
                                                     closure_audit=pipeline._shell_closure_audit)
                    if not candidate_shell['valid']:
                        wall_audit['shell_rejection'] = candidate_shell
                        raise ValueError('Refined walls failed strict shell topology validation')
                    result = candidate if merge_meshes else separate_wall_mesh(candidate, len(patches))
                    welded = candidate
                    row['shell_before_wall_refinement'] = row['shell']
                    row['shell'] = candidate_shell
                    row['shell']['mesh_audit'] = pipeline._triangle_mesh_audit(welded)
                    row['shell']['quality_contract'] = pipeline._triangle_mesh_contract_from_audit(
                        row['shell']['mesh_audit'], reference_length=float(min_building_detail),
                        require_markers=True, stage_label='Catalogue surface shell')
                    pipeline._raise_stage_contract_errors('Catalogue surface shell', row['shell']['quality_contract'])
                    row['shell']['topology_valid'] = row['shell']['valid']
                    row['shell']['quality_warning'] = bool(row['shell']['quality_contract']['warnings'])
                    wall_audit['shell_accepted'] = True
                except (ValueError, RuntimeError) as exc:
                    result, welded = original_result, original_welded
                    row['shell'] = row.get('shell_before_wall_refinement', row['shell'])
                    wall_audit.update(shell_accepted=False, shell_fallback_reason=str(exc))
                if compact_surface is not None and not wall_audit.get('shell_accepted'):
                    raise RoofAssemblyError('Compact surface rejected: '+wall_audit.get('shell_fallback_reason','unknown'))
            catalogue['compact_surface' if compact_surface is not None else
                      ('surface_refinement' if surface_refinement is not None else 'wall_triangulation')] = wall_audit
            catalogue['patch_markers'] = [dict(marker=i, roof_marker=len(patches)+i,
                source_key=p.source_key, block_id=p.block_id, patch_id=p.patch_id,
                origin=p.origin.tolist(), coefficients=p.coefficients.tolist(), area=float(p.polygon.area),
                fit_retained=bool(p.block_id not in fallback_blocks and any(
                    r.get('source_key') == p.source_key and r['outcome']=='accepted_fit' for r in rows)))
                for i, p in enumerate(patches)]
            _summarize_roof_outcomes(catalogue, rows, fallback_blocks)
            row['valid'] = True
            break
        except (RoofAssemblyError, RuntimeError, ValueError) as exc:
            row.update(valid=False, error=str(exc), rejection_stage=row['stage'])
            if assembly == 1:
                raise RoofAssemblyError(
                    f'Catalogue assembly and its one LOD1 block fallback failed. '
                    f'Inspect stage_audit["roof_catalogue"]. Last error: {exc}') from exc
            # Only attributed blocks may fall back; unknown backend failures
            # raise instead of replacing the whole city with LOD1.
            invalid_markers = row.get('shell', {}).get('invalid_patch_markers', [])
            affected = {patches[m].block_id for m in invalid_markers}
            if not affected:
                affected = set(getattr(exc, 'block_ids', ()))
            if not affected:
                raise RoofAssemblyError(
                    f'Catalogue {row["stage"]} failed without local block attribution; '
                    + ('LOD1 fallback failed; ' if fallback_blocks else '') +
                    f'no city-wide fallback is allowed. Inspect stage_audit["roof_catalogue"]: {exc}') from exc
            if affected & fallback_blocks:
                raise RoofAssemblyError(
                    f'LOD1 fallback failed for blocks {sorted(affected & fallback_blocks)}: {exc}') from exc
            fallback_blocks |= affected
            catalogue['block_fallbacks'].extend(dict(block_id=b.output,
                source_keys=b.source_keys,
                fitted_shapes=[r['shape'] for r in rows if r.get('source_key') in b.source_keys],
                rejection_stage=row['rejection_stage'], reason=str(exc)) for b in blocks if b.output in affected)
    return result


def _summarize_roof_outcomes(catalogue, rows, fallback_blocks):
    for source in rows:
        key = source.get('source_key')
        owned = [p for p in catalogue['patch_markers'] if p['source_key']==key]
        source['retained_fit_area'] = sum(p['area'] for p in owned if p['fit_retained'])
        source['fallback_area'] = sum(p['area'] for p in owned if not p['fit_retained'])
        source['final_outcome'] = ('retained_fit' if source['retained_fit_area'] > 0
                                   else 'lod1_fallback')
        source['block_fallbacks'] = [f for f in catalogue['block_fallbacks'] if key in f['source_keys']]
    catalogue['retained_fits'] = sum(r.get('retained_fit_area',0)>0 for r in rows)
    catalogue['retained_nonflat_fits'] = sum(r.get('retained_fit_area',0)>0 and r.get('shape')!='flat' for r in rows)
    for block in catalogue['blocks']:
        active = [p for p in catalogue['patch_markers'] if p['block_id']==block['block_id']]
        block['final_fit_sources'] = sorted(set(p['source_key'] for p in active if p['fit_retained']))
        block['fallback'] = block['block_id'] in fallback_blocks
