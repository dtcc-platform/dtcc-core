"""City orchestration for admitted planar roofs and conforming shared shells."""
import time
import numpy as np
from shapely.geometry import box
from shapely.ops import unary_union

from ...model import Field
from ...common.progress import progress_range, report_progress
from ..cleaning.clean import clean_regions, cleaning_record, regions_from_cleaned_city
from ..roofs import RoofFitOptions
from .roof_coverage import local_roof_coverage
from .roof_regions import source_roofs
from .roof_interfaces import RoofAssemblyError, node_interfaces
from .roof_validation import validate_partition, validate_terrain_coverage, validate_shell


def build_catalogue_surface(city, *, lod, min_building_detail, min_building_area,
                            building_mesh_triangle_size, max_mesh_size, min_mesh_angle,
                            merge_meshes, smoothing, sort_triangles, report_mesh_quality,
                            show_footprints, footprint_cleaning_plot_block, mesher,
                            pipeline_mode, stage_audit, roof_options, compact_surface=None,
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
        with progress_range(0, 20, 'Cleaning footprints and fitting LOD2 roofs...'):
            blocks, rows = _prepare_roof_blocks(
                city, bounds=bounds, lod=lod, options=options, catalogue=catalogue,
                min_building_detail=min_building_detail, min_building_area=min_building_area,
                show_footprints=show_footprints, footprint_cleaning_plot_block=footprint_cleaning_plot_block,
                lod2_method=lod2_method)
        with progress_range(20, 98, 'Meshing LOD2 roof partitions...'):
            result = _mesh_roof_blocks(
                blocks, rows, raster=raster, bounds=bounds, catalogue=catalogue,
                min_building_detail=min_building_detail, building_mesh_triangle_size=building_mesh_triangle_size,
                max_mesh_size=max_mesh_size, min_mesh_angle=min_mesh_angle, merge_meshes=merge_meshes,
                smoothing=smoothing, sort_triangles=sort_triangles, mesher=mesher,
                compact_surface=compact_surface, lod2_method=lod2_method)
        for mesh in [result] if merge_meshes else result:
            if not merge_meshes:
                mesh.fields.append(Field(name='volume_boundary', values=np.array([0.]), association='geometry'))
            mesh.stage_audit = audit
        if report_mesh_quality:
            from ...model.mixins.mesh.quality import triangle_mesh_quality, report_quality
            for mesh in [result] if merge_meshes else result:
                report_quality(triangle_mesh_quality(mesh.vertices, mesh.faces), log_fn=pipeline.info)
        report_progress(percent=100, message='LOD2 surface mesh complete')
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
    clock = time.perf_counter()
    with progress_range(0, 60, 'Cleaning LOD2 footprints...'):
        record = cleaning_record(city, min_building_detail=min_building_detail, lod=lod, isolated=isolated)
        if record is None:
            cleaned = clean_regions(city.buildings, min_building_detail=min_building_detail,
                                    lod=lod, isolated=isolated, show=show_footprints,
                                    plot_block=footprint_cleaning_plot_block)
        else:
            cleaned = regions_from_cleaned_city(city, record)
    catalogue['reused_cleaning'] = record is not None
    catalogue['stage_seconds'] = dict(cleaned_subdivision=time.perf_counter()-clock)
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
    clock = time.perf_counter()
    with progress_range(60, 100, 'Fitting LOD2 roofs...'):
        blocks, rows = source_roofs(city, cleaned, selected, options=options,
                                   reused=record is not None, method=lod2_method, detail=min_building_detail)
    catalogue['sources'] = rows
    catalogue['stage_seconds']['source_fitting_and_transfer'] = time.perf_counter()-clock
    catalogue['blocks'] = [dict(block_id=b.output, source_keys=b.source_keys,
                                fallback_levels=sorted(set(float(p.coefficients[2]) for p in b.lod1_patches)),
                                subdivision=b.subdivision_audit) for b in blocks]
    return blocks, rows


def _mesh_roof_blocks(blocks, rows, *, raster, bounds, catalogue, min_building_detail,
                      building_mesh_triangle_size, max_mesh_size, min_mesh_angle, merge_meshes,
                      smoothing, sort_triangles, mesher, compact_surface,
                      lod2_method, base_z=None):
    """Assemble and validate fitted blocks against a terrain surface."""
    from ..geometry_builders import meshes as pipeline

    from dataclasses import replace
    import time
    from .compact_surface import compact_surface as build_compact, CompactSurfaceOptions
    from .terrain_attachment import terrain_attachment
    from .roof_regions import admit_roofs
    from .wall_faces import separate_wall_mesh

    if compact_surface is None:
        compact_surface = CompactSurfaceOptions(budgeted=False, backend=mesher or 'dtcc_mesher')
    row = dict(attempt=1, stage='roof_partition')
    catalogue['assemblies'].append(row)
    timings = row['stage_seconds'] = {}
    started = time.perf_counter()
    try:
        patches = [p for block in blocks for p in block.patches]
        target = unary_union([b.polygon for b in blocks])
        row['partition'] = validate_partition(patches, target)
        row['planar_subdivision'] = {}
        # All incident domains are noded together in one local numerical frame.
        local, world, ground, local_bounds, offset = local_roof_coverage(patches, bounds,
            height_crossings=False,domains={(b.output,i):p.polygon for b in blocks for i,p in enumerate(b.lod1_patches)},
            audit=row['planar_subdivision'])
        # Terrain sees only the cleaned occupied boundary. Roof ridges and
        # source interfaces enter the shared face graph after roof admission.
        from shapely import affinity
        import shapely
        from ..roofs.partition import polygon_parts
        local_target = affinity.translate(target, *-offset)
        ground = polygon_parts(shapely.difference(box(*local_bounds), local_target))
        timings['roof_partition'] = time.perf_counter()-started
        reference, bases = None, None
        if base_z is None:
            row['stage'] = 'terrain_attachment'
            clock = time.perf_counter()
            regions = ground
            ground_mesh, backend = pipeline._build_ground_mesh_from_coverage(
                region_polygons=regions, region_markers=[-2]*len(ground),
                region_points=[np.array(p.representative_point().coords[0]) for p in regions],
                bounds=local_bounds, max_mesh_size=max_mesh_size, min_mesh_angle=min_mesh_angle,
                mesher=mesher, sort_triangles=sort_triangles,
                region_triangle_sizes=None,
                add_halo_markers=False, coordinate_tolerance=1e-8, max_steiner_points=50000)
            row['terrain_coverage'] = validate_terrain_coverage(ground_mesh,ground)
            reference, bases, row['terrain_attachment'] = terrain_attachment(
                ground_mesh, local, raster, offset=offset, smoothing=smoothing)
            row['mesher'] = backend
            timings['terrain_attachment'] = time.perf_counter()-clock
        else:
            bases = {p.block_id:base_z for p in local}
        row['stage'] = 'roof_admission'
        clock = time.perf_counter()
        local = admit_roofs(local, blocks, rows, bases)
        row['shared_subdivision'] = {}
        local = node_interfaces(local,audit=row['shared_subdivision'],canonical=True)
        timings['roof_admission'] = time.perf_counter()-clock
        row['stage'] = 'roof_junction_geometry'
        clock = time.perf_counter()
        from .roof_junctions import regularize_roof_junctions
        local,junctions = regularize_roof_junctions(local,bases=bases,detail=min_building_detail)
        for block in blocks:
            junction = junctions.get(block.output)
            if junction is not None:
                for change in junction['changes']:
                    change['xy'] = (np.asarray(change['xy'])+offset).tolist()
                for checked in junction['checked_junctions']:
                    checked['xy'] = (np.asarray(checked['xy'])+offset).tolist()
                block.subdivision_audit.update(junction)
        if any(j['changes'] for j in junctions.values()):
            local = node_interfaces(local,canonical=True)
        row['admitted_partition'] = validate_partition(local, affinity.translate(target, *-offset))
        timings['roof_junction_geometry'] = time.perf_counter()-clock
        row['stage'] = 'shared_geometry_and_meshing'
        clock = time.perf_counter()
        with progress_range(35, 95):
            result, wall_audit = build_compact(reference, local, compact_surface,
                base_z=base_z, block_bases=bases, allow_source_steps=lod2_method=='ransac')
        timings['shared_geometry_and_meshing'] = time.perf_counter()-clock
        # Transform only the completed mesh; shared geometry calculations stay local.
        result.vertices[:, :2] += offset
        clock = time.perf_counter()
        world_topology = validate_shell(result,num_patches=len(local),
                                       closure_audit=pipeline._shell_closure_audit,base_z=base_z)
        if not world_topology['valid']:
            raise RoofAssemblyError(f'World-coordinate export topology failed: {world_topology}')
        wall_audit['world_topology'] = world_topology
        timings['world_validation'] = time.perf_counter()-clock
        catalogue['numerical_frame'] = dict(xy_origin=offset.tolist(), units='metres',
                                            xy_precision_m=1e-7, coincidence_m=1e-8)
        catalogue['compact_surface'] = wall_audit
        world = [replace(p, origin=p.origin+offset) for p in local]
        catalogue['patch_markers'] = [dict(marker=i, roof_marker=len(world)+i,
            source_key=p.source_key, block_id=p.block_id, patch_id=p.patch_id,
            origin=p.origin.tolist(), coefficients=p.coefficients.tolist(), area=float(p.polygon.area),
            fit_retained=any(r.get('source_key')==p.source_key and r['outcome']=='accepted_fit' for r in rows))
            for i,p in enumerate(world)]
        _summarize_roof_outcomes(catalogue, rows, set())
        row.update(valid=True, shell=wall_audit['topology'])
        wall_audit['shell_accepted'] = True
        timings['total'] = time.perf_counter()-started
        return result if merge_meshes else separate_wall_mesh(result, len(world))
    except (RoofAssemblyError, RuntimeError, ValueError) as exc:
        row.update(valid=False, error=str(exc), rejection_stage=row['stage'])
        timings['total'] = time.perf_counter()-started
        raise RoofAssemblyError(f'LOD2 {row["stage"]} contract failed: {exc}',
            block_ids=getattr(exc, 'block_ids', ())) from exc


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
