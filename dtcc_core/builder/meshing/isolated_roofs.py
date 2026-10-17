"""Independent fitted building skins above an unconstrained terrain mesh."""
from copy import deepcopy

import numpy as np

from ...model import Bounds, Mesh, Raster
from ..roofs import RoofFitOptions
from .catalogue_surface import _prepare_roof_blocks, _mesh_roof_blocks, _summarize_roof_outcomes


def build_isolated_roofs(city, *, lod, min_building_detail, min_building_area,
                         building_mesh_triangle_size, max_mesh_size, min_mesh_angle,
                         merge_meshes, smoothing, sort_triangles, report_mesh_quality,
                         show_footprints, footprint_cleaning_plot_block, mesher,
                         pipeline_mode, stage_audit, roof_options, surface_refinement=None,
                         compact_surface=None, lod2_method='catalogue'):
    """Fit once, mesh each cleaned output independently, and retain separate indices.

    Compact skins use fitted domains and an explicit base, without ground
    triangulation. The noncompact assembler uses a discarded flat ground ring.
    Buildings may overlap each other and the separately meshed terrain.
    """
    from ..geometry_builders import meshes as pipeline

    audit = {} if stage_audit is None else stage_audit
    catalogue = dict(sources=[], blocks=[], assemblies=[], block_fallbacks=[],
                     building_assemblies=[], patch_markers=[], volume_meshing_supported=False)
    audit.update(requested_lod='lod2', lod2_method=lod2_method, isolated=True)
    audit['roof_catalogue'] = catalogue
    attempt = pipeline._start_stage_audit_attempt(
        audit, label='isolated_lod2', retry_reason=None, backend='surface',
        requested_mesher=mesher, max_mesh_size=max_mesh_size, min_mesh_angle=min_mesh_angle,
        domain_height=0., tetgen_switches=None, tetgen_switch_overrides=None)
    try:
        max_mesh_size = pipeline._normalize_max_mesh_size(max_mesh_size)
        _, raster = pipeline._require_city_terrain_raster(city, max_mesh_size=max_mesh_size)
        terrain, tree, reach, _ = pipeline._isolated_terrain(
            raster, max_mesh_size=max_mesh_size, min_mesh_angle=min_mesh_angle,
            smoothing=smoothing, mesher=mesher)
        blocks, rows = _prepare_roof_blocks(
            city, bounds=pipeline._raster_bounds_tuple(raster), lod=lod,
            options=RoofFitOptions() if roof_options is None else roof_options,
            catalogue=catalogue, min_building_detail=min_building_detail,
            min_building_area=min_building_area, show_footprints=show_footprints,
            footprint_cleaning_plot_block=footprint_cleaning_plot_block,
            lod2_method=lod2_method, isolated=True)
        shells, region_counts = [], []
        for block in blocks:
            base_z = pipeline._lowest_terrain(terrain, tree, reach, block.polygon)
            if base_z is None:
                raise ValueError(f'No terrain support for building {block.output}')
            margin = max(1., 2.*float(min_building_detail))
            xmin, ymin, xmax, ymax = block.polygon.bounds
            bounds = (xmin-margin, ymin-margin, xmax+margin, ymax+margin)
            base = None
            if compact_surface is None:
                base = Raster(data=np.full((2, 2), base_z))
                base.set_bounds(Bounds(*bounds))
            sources = deepcopy([r for r in rows if r.get('source_key') in block.source_keys
                                or r.get('block_id') == block.output])
            local = dict(sources=sources, blocks=[dict(block_id=block.output)],
                         assemblies=[], block_fallbacks=[])
            catalogue['building_assemblies'].append(dict(
                building=block.output, base_z=base_z, roof_catalogue=local))
            reference = _mesh_roof_blocks(
                [block], sources, raster=base, bounds=bounds, catalogue=local,
                min_building_detail=min_building_detail, building_mesh_triangle_size=building_mesh_triangle_size,
                max_mesh_size=max_mesh_size, min_mesh_angle=min_mesh_angle, merge_meshes=True,
                smoothing=0, sort_triangles=sort_triangles, mesher=mesher,
                compact_surface=compact_surface, surface_refinement=surface_refinement,
                lod2_method=lod2_method, base_z=base_z)
            keep = np.asarray(reference.markers) >= 0
            faces = reference.faces[keep]
            used, inverse = np.unique(faces, return_inverse=True)
            shell = Mesh(vertices=reference.vertices[used], faces=inverse.reshape((-1, 3)),
                         markers=np.asarray(reference.markers)[keep])
            shells.append((block.output, shell))
            region_counts.append(len(local['patch_markers']))
            catalogue['block_fallbacks'].extend(local['block_fallbacks'])
            for patch in local['patch_markers']:
                catalogue['patch_markers'].append(dict(patch, marker=len(catalogue['patch_markers'])))
        # Global marker IDs match the shared isolated-shell concatenation rule.
        total = len(catalogue['patch_markers'])
        for patch in catalogue['patch_markers']:
            patch['roof_marker'] = total + patch['marker']
        _summarize_roof_outcomes(catalogue, rows, {f['block_id'] for f in catalogue['block_fallbacks']})
        result = pipeline._combine_isolated_shells(terrain, shells, merge_meshes=merge_meshes,
                                                   region_counts=region_counts)
        for mesh in [result] if merge_meshes else result:
            mesh.stage_audit = audit
            if report_mesh_quality:
                from ...model.mixins.mesh.quality import triangle_mesh_quality, report_quality
                report_quality(triangle_mesh_quality(mesh.vertices, mesh.faces), log_fn=pipeline.info)
        pipeline._mark_stage_audit_success(audit, attempt)
        return result
    except Exception as exc:
        pipeline._mark_stage_audit_failure(attempt, exc)
        raise
