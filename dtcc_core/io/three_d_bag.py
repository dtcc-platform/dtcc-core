"""Explicit enrichment of the audited 3DBAG b3_h_dak CityJSON convention.

Source attributes remain evidence on their original owner. See
docs/design/3dbag-attribute-mapping.md for meanings and unsupported conventions.
"""

import math

from .cityjson.cityjson import load as load_cityjson
from ..model import Building, BuildingPart
from ..model._standard_schema import validate_admitted


MAPPING_ID = 'https://github.com/dtcc-platform/dtcc-core/mappings/3dbag/b3_h_dak/1'
NAP_CRS = 'https://www.opengis.net/def/crs/EPSG/0/5709'
_ELEVATIONS = ('b3_h_maaiveld', 'b3_h_dak_min', 'b3_h_dak_50p',
               'b3_h_dak_70p', 'b3_h_dak_max')
_ROOF_CODES = frozenset({'horizontal', 'multiple horizontal', 'slanted'})


def load_3dbag(path, *, extent_policy='validate', validate_schema=True):
    """Load strict CityJSON and map the b3_h_dak attribute convention explicitly.

    Requires source CRS EPSG:7415. Known, non-null elevations become signed metre
    records relative to NAP; physical roof codes retain their 3DBAG namespace.
    Null/missing values stay absent. No scalar height, date or release is inferred.
    Newer b3_h_50p/b3_h_70p/b3_h_min/b3_h_max spellings require a separate mapping.
    Existing destination attributes fail rather than being overwritten or merged.

    validate_schema=False bypasses only DTCC semantic evaluation, never strict
    source admission or this mapping's required source interpretation checks.
    """
    from pyproj import CRS
    from ..datasets.schema import (DatasetContext, DatasetIdentity, DatasetMetadata,
                                   DatasetProvenance, DatasetPresentation, DatasetRequest)

    if type(validate_schema) is not bool:
        raise ValueError('validate_schema must be a boolean')
    # Admit geometry once; evaluate the selected semantics once, after enrichment.
    city = load_cityjson(path, strict=True, extent_policy=extent_policy,
                         validate_schema=False)
    if not city.transform.srs or CRS.from_user_input(city.transform.srs) != CRS(7415):
        raise ValueError('3DBAG b3_h_dak mapping requires source CRS EPSG:7415 (RD New + NAP)')
    seen = 0
    mapped_objects = 0
    stack = [city]
    while stack:
        obj = stack.pop()
        stack.extend(child for group in obj.children.values() for child in group)
        if not isinstance(obj, (Building, BuildingPart)):
            continue
        attrs = obj.attributes
        if any(key in attrs for key in ('b3_h_50p', 'b3_h_70p', 'b3_h_min', 'b3_h_max')):
            raise ValueError(f"objects[{obj.id!r}]: newer 3DBAG elevation names are outside the b3_h_dak mapping")
        elevations = []
        for key in _ELEVATIONS:
            if key not in attrs:
                continue
            seen += 1
            value = attrs[key]
            if value is None:
                continue
            if type(value) not in (int, float) or (type(value) is float and not math.isfinite(value)):
                raise ValueError(f"objects[{obj.id!r}].attributes[{key!r}] must be a finite number or null")
            elevations.append({
                'value': value, 'unit': 'm', 'vertical_reference': NAP_CRS,
                'reference': {'value': key, 'code_space': MAPPING_ID + '/elevation-references'},
                'source': f'{MAPPING_ID}; source attribute {key} on this object',
            })
        additions = {}
        if elevations:
            additions['elevation_measurements'] = elevations
        roof = attrs.get('b3_dak_type')
        if roof is not None and type(roof) is not str:
            raise ValueError(f"objects[{obj.id!r}].attributes['b3_dak_type'] must be a string or null")
        if roof in _ROOF_CODES:
            additions['roof_type'] = {'value': roof, 'code_space': MAPPING_ID + '/roof-types'}
        for key in additions:
            if key in attrs:
                raise ValueError(f"objects[{obj.id!r}].attributes[{key!r}]: 3DBAG mapping destination already exists")
        attrs.update(additions)
        mapped_objects += bool(additions)
    if not seen:
        raise ValueError('No supported 3DBAG b3_h_dak/ground elevation attributes found')
    if validate_schema:
        validate_admitted(city, city.schema_id, city.schema_version)
    if city.dataset_context is None:
        city.dataset_context = DatasetContext(
            identity=DatasetIdentity(name='3dbag-import', title='3DBAG import'),
            metadata=DatasetMetadata(crs=[city.transform.srs]),
            provenance=DatasetProvenance(), presentation=DatasetPresentation(),
            request=DatasetRequest(dataset_name='3dbag-import'),
        )
    city.dataset_context.provenance.processing_steps.append({
        'operation': 'map source attributes to qualified DTCC values',
        'mapping': MAPPING_ID, 'mapped_object_count': mapped_objects,
        'source_release': None,
        'source_attributes': 'preserved on original owners; no height or date inference',
    })
    if not isinstance(path, dict):
        city.dataset_context.provenance.sources.append({'path': str(path)})
    return city
