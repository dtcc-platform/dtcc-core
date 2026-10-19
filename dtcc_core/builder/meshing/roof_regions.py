"""Resolve source-labelled cleaned cells and transfer ephemeral fitted roofs."""
from dataclasses import dataclass, field
import time

import numpy as np

from ...model import GeometryType
from ...common.progress import report_progress
from ..cleaning.blocks import source_footprint, source_id
from ..roofs import RoofFitResult, RoofPatch, fit_roof
from ..roofs.partition import polygon_parts


@dataclass
class RoofBlock:
    """Cleaned output polygon with fitted and fallback patches and source identities."""
    output: int
    polygon: object
    patches: list
    lod1_patches: list
    source_keys: list[str]
    subdivision_audit: dict = field(default_factory=dict)



def source_roofs(city, cleaned, selected, *, options, reused=False,
                 method='catalogue', detail=0.5):
    """Use exact cleaner associations; never recover partitions by proximity."""
    fitted, rows, blocks = {}, [], []
    for index, (polygon, output) in enumerate(selected):
        cells = cleaned.output_cells[output] if cleaned.output_cells else []
        from ..cleaning.subdivision import canonical_source_subdivision
        source_keys = {}
        for _, label in cells:
            building_index, part_index = cleaned.labels[label]
            source = city.buildings[building_index]
            if part_index is not None:
                source = source.building_parts[part_index]
            source_keys[label] = str(source_id(source))
        try:
            polygon, cells, subdivision = canonical_source_subdivision(
                cells,domain=polygon,source_keys=source_keys)
        except ValueError as exc:
            from .roof_interfaces import RoofAssemblyError
            raise RoofAssemblyError(f'LOD2 cleaned_subdivision block {output}: {exc}',
                                    block_ids=[output]) from exc
        from shapely.ops import unary_union
        grouped_cells = {}
        for cell, label in cells:
            grouped_cells.setdefault(label, []).append(cell)
        cells = [(unary_union(parts), label) for label, parts in sorted(grouped_cells.items())]
        keys, patches = [], []
        from ..cleaning.levels import group_levels
        original_cells = cleaned.output_cells[output] if cleaned.output_cells else []
        labels = sorted({label for _, label in original_cells})
        levels = dict(zip(labels, group_levels([cleaned.label_roof_z[label] for label in labels],
            [sum(p.area for p,i in original_cells if i==label) for label in labels], detail))) if labels else {}
        fallback = []
        for cell, label in cells:
            key = source_keys[label]
            for piece in polygon_parts(cell):
                fallback.append(RoofPatch(piece, np.asarray(piece.centroid.coords[0]),
                    np.array([0.,0.,levels[label]]), key, len(fallback), output))
        if not cells:
            for piece in polygon_parts(polygon):
                fallback.append(RoofPatch(piece,np.asarray(piece.centroid.coords[0]),
                    np.array([0.,0.,cleaned.output_roof_z[output]]),f'block:{output}:cleaner-fallback',
                    len(fallback),output))
        for cell, label in cells:
            target = cell
            if target.is_empty or target.area <= 0:
                continue
            building_index, part_index = cleaned.labels[label]
            parent = city.buildings[building_index]
            source = parent if part_index is None else parent.building_parts[part_index]
            key = str(source_id(source))
            keys.append(key)
            if label not in fitted:
                footprint = source_footprint(source)
                evidence = ('preserved source footprint' if 'source_footprint' in source.attributes
                            else 'source LOD0 footprint')
                if footprint is None:
                    footprint = cleaned.label_polygons[label]
                    evidence = 'label footprint'
                cloud = source.get_geometry(GeometryType.POINT_CLOUD)
                points = np.asarray(cloud.points) if cloud is not None else np.empty((0, 3))
                try:
                    ground = float(source.attributes.get('ground_height', parent.attributes.get('ground_height', np.nan)))
                except (TypeError, ValueError):
                    ground = float('nan')
                started = time.perf_counter()
                if reused and 'source_footprint' not in source.attributes:
                    result = RoofFitResult(None, dict(accepted=False, rejection_reasons=[
                        'original source footprint unavailable on previously cleaned source']))
                    evidence = 'cleaned cell only; original orientation evidence unavailable'
                elif method == 'ransac':
                    from ..roofs.ransac import fit_ransac
                    result = fit_ransac(footprint, points, float(ground), source_key=key)
                else:
                    result = fit_roof(footprint, points, float(ground), source_key=key, options=options)
                fitted[label] = result
                row = dict(source_key=key, source_id=str(source_id(source)), label=label,
                           building_index=building_index, part_index=part_index,
                           evidence=evidence, ground_height=float(ground) if np.isfinite(ground) else None,
                           points=len(points), shape=result.model.shape if result.model else None,
                           parameters=result.model.parameters.tolist() if result.model else None,
                           fit_seconds=time.perf_counter()-started,
                           outcome='candidate_fit' if result.accepted else 'lod1_fallback',
                           diagnostics=result.diagnostics)
                if part_index is not None and cloud is None:
                    row['diagnostics']['rejection_reasons'].append('source part has no observations; parent pooling is not used')
                rows.append(row)
            result = fitted[label]
            if result.accepted:
                from .roof_validation import validate_partition
                try:
                    transferred = result.model.patches(target, block_id=output)
                    validate_partition(transferred, target)
                except (ValueError, RuntimeError) as exc:
                    result.diagnostics.setdefault('transfer_rejections', []).append(str(exc))
                    for row in rows:
                        if row.get('source_key') == key:
                            row['outcome'] = 'lod1_fallback'
                    patches.extend(p for p in fallback if p.source_key == key
                                   and p.polygon.intersection(target).area > 0)
                else:
                    patches.extend(transferred)
            else:
                patches.extend(p for p in fallback if p.source_key == key
                               and p.polygon.intersection(target).area > 0)
        block = RoofBlock(output, polygon, patches, fallback, list(dict.fromkeys(keys)), subdivision)
        if not cells:
            # The cleaner deliberately discarded the partition. Its single
            # output roof is the only correspondence that remains available.
            block.patches = list(block.lod1_patches)
            rows.append(dict(block_id=output, source_keys=[], outcome='lod1_fallback',
                             diagnostics=dict(accepted=False, rejection_reasons=[
                                 'cleaner retained no source cells; source partition unavailable'])))
        blocks.append(block)
        report_progress(current=index+1, total=len(selected),
                        message=f'Fitted LOD2 roof group {index+1}/{len(selected)}')
    return blocks, rows


def admit_roofs(patches, blocks, rows, bases):
    """Admit transferred candidates against foundations before shell meshing.

    An affine patch attains its extrema on its ring vertices. A candidate that
    extrapolates through its foundation is inadmissible over the cleaned domain;
    the source's declared flat evidence model is selected here, without moving
    either planes or boundaries. Shared-geometry and meshing errors never enter
    this decision.
    """
    from dataclasses import replace
    from .roof_interfaces import RoofAssemblyError

    rejected = {}
    for row in rows:
        if row.get('diagnostics', {}).get('transfer_rejections'):
            rejected.setdefault(row['source_key'], []).append(dict(stage='roof_transfer',
                reason='; '.join(row['diagnostics']['transfer_rejections'])))
    for patch in patches:
        minimum = min(float(np.min(patch.height(np.asarray(r.coords)))) for r in
                      [patch.polygon.exterior, *patch.polygon.interiors])
        if minimum <= bases[patch.block_id] + 1e-8:
            rejected.setdefault(patch.source_key, []).append(dict(block_id=patch.block_id,
                patch_id=patch.patch_id, minimum_roof_m=minimum, foundation_m=bases[patch.block_id],
                reason='transferred affine roof intersects its terrain foundation'))
    levels = {}
    for block in blocks:
        for patch in block.lod1_patches:
            levels.setdefault((patch.block_id,patch.source_key), set()).add(float(patch.coefficients[2]))
    for row in rows:
        key = row.get('source_key')
        if row['outcome'] == 'candidate_fit' or key in rejected:
            if key in rejected:
                row.update(outcome='lod1_fallback', admission=dict(accepted=False,
                    rejections=rejected[key], stage='roof_geometry'))
            else:
                row.update(outcome='accepted_fit', admission=dict(accepted=True, stage='roof_geometry'))
    result = []
    for patch in patches:
        if patch.source_key in rejected:
            choices = levels.get((patch.block_id,patch.source_key), set())
            if len(choices) != 1:
                raise RoofAssemblyError('Roof admission needs one declared flat source model: '
                    +patch.source_key, block_ids=[patch.block_id])
            level = next(iter(choices))
            if level <= bases[patch.block_id]+1e-8:
                raise RoofAssemblyError('Declared source roof is below its terrain foundation: '
                    +patch.source_key, block_ids=[patch.block_id])
            patch = replace(patch, coefficients=np.array([0., 0., level]))
        result.append(patch)
    return result
