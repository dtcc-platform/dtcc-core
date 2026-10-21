"""Recognized CityJSON property spellings at the DTCC boundary."""

_BUILDING = {
    'measuredHeight': 'measured_height', 'roofType': 'roof_type',
    'storeysAboveGround': 'storeys_above_ground',
    'storeysBelowGround': 'storeys_below_ground',
}
_VEGETATION = {'crownDiameter': 'crown_diameter', 'trunkDiameter': 'trunk_diameter'}
_NAMES = {'Building': _BUILDING, 'BuildingPart': _BUILDING,
          'SolitaryVegetationObject': _VEGETATION,
          'PlantCover': {'averageHeight': 'average_height'},
          'TrafficArea': {'surfaceMaterial': 'surface_material'},
          'AuxiliaryTrafficArea': {'surfaceMaterial': 'surface_material'},
          'WaterSurface': {'waterLevel': 'water_level'}}


class _AttributeNameError(ValueError):
    """A spelling cannot be translated without losing its identity."""


def _rename(attributes, names):
    if not isinstance(attributes, dict):
        raise ValueError('CityJSON attributes must be an object')
    # Destination-only keys are ambiguous too: preserving one on import and
    # translating it on export would silently change its original spelling.
    for source, target in names.items():
        if target in attributes:
            raise _AttributeNameError(
                f'Attribute {target!r} is reserved by the {source!r} -> {target!r} '
                'mapping; cannot preserve this spelling unambiguously'
            )
    return {names.get(key, key): value for key, value in attributes.items()}


def read_attributes(attributes, feature_type):
    """Rename CityJSON attribute names to DTCC names for a feature type.

    For example, ``measuredHeight`` becomes ``measured_height`` on buildings.
    Names without a mapping are kept.

    Parameters
    ----------
    attributes : dict
        CityJSON attributes.
    feature_type : str
        CityJSON feature type, such as ``"Building"``.

    Returns
    -------
    dict
        Attributes with DTCC names.

    Raises
    ------
    ValueError
        If ``attributes`` is not a dict, or already contains a DTCC name that a
        CityJSON name would be renamed to.
    """
    return _rename(attributes, _NAMES.get(feature_type, {}))


def write_attributes(attributes, feature_type):
    """Rename DTCC attribute names to CityJSON names for a feature type.

    The inverse of :func:`read_attributes`: for example, ``measured_height``
    becomes ``measuredHeight`` on buildings.

    Parameters
    ----------
    attributes : dict
        DTCC attributes.
    feature_type : str
        CityJSON feature type, such as ``"Building"``.

    Returns
    -------
    dict
        Attributes with CityJSON names.

    Raises
    ------
    ValueError
        If ``attributes`` is not a dict, or already contains a CityJSON name
        that a DTCC name would be renamed to.
    """
    return _rename(attributes, {target: source for source, target in _NAMES.get(feature_type, {}).items()})
