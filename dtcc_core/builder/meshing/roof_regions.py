"""Resolve source-labelled cleaned cells and transfer ephemeral fitted roofs."""
from dataclasses import dataclass

import numpy as np

from ...model import GeometryType
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

    def flat_patches(self):
        """The ordinary LOD1 level regions, retaining their source ownership."""
        return self.lod1_patches


def lod1_patches(city, cleaned, output, polygon, regions):
    """Partition ordinary flat fallback levels by the cleaner's source ownership."""
    patches = []
    cells = cleaned.output_cells[output] if cleaned.output_cells else []
    for region, level, _ in regions:
        if cells:
            owned = []
            for cell, label in cells:
                building, part = cleaned.labels[label]
                source = city.buildings[building]
                if part is not None:
                    source = source.building_parts[part]
                owned.append((cell.intersection(region).intersection(polygon), str(source_id(source))))
        else:
            owned = [(region.intersection(polygon), f'block:{output}:cleaner-fallback')]
        for target, key in owned:
            for piece in polygon_parts(target):
                patches.append(RoofPatch(piece, np.array(piece.centroid.coords[0]),
                                         np.array([0., 0., level]), key, len(patches), output))
    return patches


def source_roofs(city, cleaned, selected, *, options, fallback_regions, reused=False,
                 method='catalogue'):
    """Use exact cleaner associations; never recover partitions by proximity."""
    fitted, rows, blocks = {}, [], []
    for polygon, output in selected:
        cells = cleaned.output_cells[output] if cleaned.output_cells else []
        keys, patches = [], []
        regions, level_shift = fallback_regions(output, polygon)
        fallback = lod1_patches(city, cleaned, output, polygon, regions)
        for cell, label in cells:
            target = cell.intersection(polygon)
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
                           outcome='accepted_fit' if result.accepted else 'lod1_fallback',
                           diagnostics=result.diagnostics)
                if part_index is not None and cloud is None:
                    row['diagnostics']['rejection_reasons'].append('source part has no observations; parent pooling is not used')
                rows.append(row)
            result = fitted[label]
            if result.accepted:
                if method == 'ransac':
                    try:
                        patches.extend(result.model.patches(target, block_id=output))
                    except (ValueError, RuntimeError) as exc:
                        # Missing coverage invokes the existing validated block fallback.
                        result.diagnostics.setdefault('transfer_rejections', []).append(str(exc))
                else:
                    patches.extend(result.model.patches(target, block_id=output))
            else:
                patches.extend(p for p in fallback if p.source_key == key
                               and p.polygon.intersection(target).area > 0)
        block = RoofBlock(output, polygon, patches, fallback, list(dict.fromkeys(keys)))
        if not cells:
            # The cleaner deliberately discarded the partition. Its single
            # output roof is the only correspondence that remains available.
            block.patches = block.flat_patches()
            rows.append(dict(block_id=output, source_keys=[], outcome='lod1_fallback',
                             diagnostics=dict(accepted=False, rejection_reasons=[
                                 'cleaner retained no source cells; source partition unavailable'])))
        blocks.append(block)
    return blocks, rows
