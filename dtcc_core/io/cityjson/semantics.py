"""CityJSON surface assignments mapped onto generic native semantic regions."""

import numpy as np

from ...model import SemanticRegion
from ...model._standard_schema import SEMANTIC_NAMESPACE
from .attributes import read_attributes, write_attributes

# Allowed unextended Building semantic tokens in CityJSON 2.0.2 section 3.3.
# This is an external-format mapping; the selected DTCC profile may be narrower.
BUILDING_SURFACE_TYPES = frozenset({
    'RoofSurface', 'GroundSurface', 'WallSurface', 'ClosureSurface',
    'OuterCeilingSurface', 'OuterFloorSurface', 'Window', 'Door',
    'InteriorWallSurface', 'CeilingSurface', 'FloorSurface',
})

TRANSPORTATION_SURFACE_TYPES = frozenset({
    'TrafficArea', 'AuxiliaryTrafficArea', 'TransportationMarking', 'TransportationHole',
})
WATER_SURFACE_TYPES = frozenset({'WaterSurface', 'WaterGroundSurface', 'WaterClosureSurface'})
SURFACE_TYPES = BUILDING_SURFACE_TYPES | TRANSPORTATION_SURFACE_TYPES | WATER_SURFACE_TYPES
FEATURE_SURFACE_TYPES = {
    'Building': BUILDING_SURFACE_TYPES, 'BuildingPart': BUILDING_SURFACE_TYPES,
    **{name: TRANSPORTATION_SURFACE_TYPES
       for name in ('Road', 'Railway', 'TransportSquare')},
    'WaterBody': WATER_SURFACE_TYPES,
}


def validate_region_owner(regions, feature_type):
    """The external vocabulary belongs to its feature, even with schema bypass."""
    allowed = FEATURE_SURFACE_TYPES.get(feature_type, frozenset())
    for region in regions:
        if (not region.semantic_type.startswith(SEMANTIC_NAMESPACE) or
                region.semantic_type[len(SEMANTIC_NAMESPACE):] not in allowed):
            raise NotImplementedError(f'{feature_type} semantic regions require its mapped external vocabulary')


def read_regions(semantics, count):
    if not isinstance(semantics, dict) or set(semantics) != {'surfaces', 'values'}:
        raise ValueError("CityJSON semantics requires surfaces and values")
    entities, assignments = semantics['surfaces'], semantics['values']
    if not isinstance(entities, list) or not isinstance(assignments, list) or len(assignments) != count:
        raise ValueError("CityJSON semantic assignments must match the surface count")
    indices = [[] for _ in entities]
    for index, assigned in enumerate(assignments):
        if assigned is None:
            continue
        if type(assigned) is not int or not 0 <= assigned < len(entities):
            raise ValueError("CityJSON semantic index is out of range")
        indices[assigned].append(index)
    regions = []
    for entity, members in zip(entities, indices):
        if not isinstance(entity, dict):
            raise ValueError("CityJSON semantic entity must be an object")
        name = entity.get('type')
        if not isinstance(name, str) or name not in SURFACE_TYPES:
            raise NotImplementedError("Strict CityJSON currently requires unextended semantic type names")
        attrs = read_attributes({key: value for key, value in entity.items()
                                 if key not in ('type', 'parent', 'children')}, name)
        if any(type(value) not in (str, int, float, bool) for value in attrs.values()):
            raise ValueError("CityJSON semantic attributes must be scalar values")
        regions.append(SemanticRegion(SEMANTIC_NAMESPACE + name,
                                      np.array(members, dtype=np.int64), attributes=attrs))
    # CityJSON can supply either direction. Normalize to the single native
    # parent authority; a supplied inverse list must agree with all known edges.
    def reference(value):
        if type(value) is not int or not 0 <= value < len(regions):
            raise ValueError("CityJSON semantic relationship index is out of range")
        return value

    for region, entity in zip(regions, entities):
        if 'parent' in entity:
            region.parent = reference(entity['parent'])
    declared_children = {}
    for index, entity in enumerate(entities):
        if 'children' not in entity:
            continue
        children = entity['children']
        if not isinstance(children, list):
            raise ValueError("CityJSON semantic children must be an array of indices")
        children = [reference(child) for child in children]
        if len(set(children)) != len(children):
            raise ValueError("Duplicate CityJSON semantic child")
        declared_children[index] = set(children)
        for child in children:
            region = regions[child]
            if region.parent is not None and region.parent != index:
                raise ValueError("Conflicting CityJSON semantic parent/children relationships")
            region.parent = index
    actual_children = {}
    for index, region in enumerate(regions):
        if region.parent is not None:
            actual_children.setdefault(region.parent, set()).add(index)
    for index, children in declared_children.items():
        if children != actual_children.get(index, set()):
            raise ValueError("Conflicting CityJSON semantic parent/children relationships")
    # Generic cycle checks run in the final native admission, alongside indices.
    return regions


def write_regions(regions, count):
    """Native admission is performed once by the strict export boundary."""
    entities, assignments = [], [None] * count
    for region in regions:
        if not region.semantic_type.startswith(SEMANTIC_NAMESPACE):
            raise NotImplementedError("Semantic URI has no CityJSON mapping in this profile")
        name = region.semantic_type[len(SEMANTIC_NAMESPACE):]
        if name not in SURFACE_TYPES:
            raise NotImplementedError("Semantic URI has no CityJSON type spelling")
        if region.id is not None:
            raise NotImplementedError("Region IDs require a defined CityJSON extension mapping")
        if {'type', 'parent', 'children'} & region.attributes.keys():
            raise ValueError("Reserved CityJSON semantic attribute name")
        if any(type(value) not in (str, int, float, bool) for value in region.attributes.values()):
            raise ValueError("CityJSON semantic attributes must be scalar values")
        for index in region.indices:
            assignments[int(index)] = len(entities)
        entities.append({'type': name, **write_attributes(region.attributes, name)})
    for index, region in enumerate(regions):
        if region.parent is not None:
            entities[index]['parent'] = region.parent
            entities[region.parent].setdefault('children', []).append(index)
    return {'surfaces': entities, 'values': assignments}
