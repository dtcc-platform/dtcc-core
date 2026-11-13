"""Private model adapter for closed LOD1 column geometry."""

from collections import defaultdict
from collections.abc import Sequence
from copy import copy, deepcopy

import numpy as np
from shapely.affinity import translate

from ...model import (
    Building,
    BuildingPart,
    City,
    Mesh,
    MultiSurface,
    Solid,
    Surface,
    Transform,
)
from ...model._standard_schema import DEFAULT_VERSION, SCHEMA_ID, SEMANTIC_NAMESPACE
from ..logging import warning
from .buildings import building_from_footprint
from .prisms import (
    _check_metadata,
    _check_parameters,
    _decompose_prisms,
    _extract_faces,
    _frame,
    _geometry_key,
    _parts,
    _PrismError,
    _PrismFace,
    _reject_overlapping_components,
    _volume_budget,
)


def _objects(roots):
    stack, seen, result = list(reversed(roots)), set(), []
    while stack:
        obj = stack.pop()
        if id(obj) in seen:
            raise ValueError("Shared or cyclic object containment is unsupported")
        seen.add(id(obj))
        result.append(obj)
        stack.extend(
            reversed([child for group in obj.children.values() for child in group])
        )
    return result


def _inputs(target):
    if isinstance(target, City):
        sources = list(target.buildings)
        objects = _objects([target])
    elif isinstance(target, Sequence) and not isinstance(target, (str, bytes)):
        sources = list(target)
        if not all(isinstance(source, Building) for source in sources):
            raise TypeError("Expected a City or a sequence of Buildings")
        objects = _objects(sources)
    else:
        raise TypeError("Expected a City or a sequence of Buildings")
    reserved = set()
    for obj in objects:
        if not isinstance(obj.id, str) or not obj.id:
            raise ValueError("Every object requires a nonempty string ID")
        if obj.id in reserved:
            raise ValueError(f"Duplicate object ID {obj.id!r}")
        reserved.add(obj.id)
    referenced = {
        id for obj in objects for targets in obj.relations.values() for id in targets
    }
    return sources, reserved, referenced


def _validate_object(obj, *, part=False):
    expected = BuildingPart if part else Building
    semantic = SEMANTIC_NAMESPACE + expected.__name__
    if type(obj) is not expected or obj.semantic_type not in (None, semantic):
        raise _PrismError(
            "unsupported_metadata",
            "custom object types or semantics require explicit transfer",
        )
    if (
        obj.profile_id is not None
        or obj.profile_version is not None
        or obj.schema_id is not None
        or obj.schema_version is not None
        or obj.dataset_context is not None
        or obj.relations
    ):
        raise _PrismError(
            "unsupported_metadata",
            "object profile, schema, context or relations require explicit transfer",
        )
    if part and obj.attributes:
        raise _PrismError(
            "unsupported_metadata",
            "child-specific attributes cannot be assigned to derived pieces",
        )


def _metre_srs(srs):
    if not srs:
        return
    from pyproj import CRS
    from pyproj.exceptions import CRSError

    try:
        crs = CRS.from_user_input(srs)
    except CRSError as error:
        raise _PrismError(
            "unsupported_frame", "SRS must identify a metre coordinate system"
        ) from error
    if (
        not crs.axis_info
        or not (crs.is_projected or crs.is_engineering)
        or any(
            not np.isclose(axis.unit_conversion_factor, 1.0, rtol=0, atol=1e-12)
            for axis in crs.axis_info
        )
    ):
        raise _PrismError("unsupported_frame", "coordinates must use local metres")


def _check_discarded_footprint(geometry):
    if not isinstance(geometry, (Surface, MultiSurface, Mesh, Solid)):
        raise _PrismError(
            "unsupported_metadata", "non-polygon LOD0 requires explicit transfer"
        )
    _check_metadata(geometry)
    if isinstance(geometry, Mesh) and geometry.markers.size:
        raise _PrismError(
            "unsupported_metadata", "LOD0 mesh markers require explicit transfer"
        )
    for surface in getattr(geometry, "surfaces", ()):
        _check_metadata(surface)


def _source_geometry(source, referenced, inherited_srs):
    records, frames = [], []
    stack = [(source, np.zeros(3), inherited_srs, False)]
    while stack:
        obj, parent_offset, parent_srs, part = stack.pop()
        _validate_object(obj, part=part)
        if obj.id in referenced:
            raise _PrismError(
                "unsupported_metadata", "a relation references an object being replaced"
            )
        offset, srs = _frame(obj)
        if parent_srs and srs and parent_srs != srs:
            raise _PrismError(
                "unsupported_frame", "object hierarchy has incompatible SRS labels"
            )
        srs = srs or parent_srs
        offset += parent_offset
        for record in obj.geometry.values():
            if not part and record.lod == "0" and record.role is None:
                _check_discarded_footprint(record.geometry)
                continue
            if record.lod != "1" or record.role is not None:
                raise _PrismError(
                    "unsupported_metadata",
                    "extra LODs or representation roles cannot be discarded",
                )
            faces, (translation, geometry_srs) = _extract_faces(record.geometry)
            if srs and geometry_srs and srs != geometry_srs:
                raise _PrismError(
                    "unsupported_frame", "object and geometry SRS labels differ"
                )
            frames.append((tuple(offset + translation), geometry_srs or srs))
            for face in faces:
                records.append(_PrismFace(face.rings, len(records)))
        for group in obj.children.values():
            stack.extend((child, offset, srs, True) for child in reversed(group))
    if not records:
        raise _PrismError("missing_geometry", "building has no LOD1 faces")
    translations = {offset for offset, _ in frames}
    labels = {srs for _, srs in frames if srs}
    if len(translations) != 1 or len(labels) > 1:
        raise _PrismError("unsupported_frame", "LOD1 solids must share one local frame")
    srs = next(iter(labels), inherited_srs)
    _metre_srs(srs)
    return records, Transform(srs=srs, affine=np.eye(4)), next(iter(translations))


def _piece_ids(source_id, count, reserved):
    if count == 1:
        return [source_id]
    result, suffix = [], 1
    while len(result) < count:
        candidate = f"{source_id}#{suffix}"
        if candidate not in reserved:
            result.append(candidate)
            reserved.add(candidate)
        suffix += 1
    return result


def _placed_component(component, offset, tolerance, index):
    """Resolve a translation into City coordinates and charge floating error."""
    levels = []
    placement = 0.0
    for roof, polygon in component.levels:
        shifted = translate(polygon, offset[0], offset[1])
        if not shifted.is_valid or shifted.area <= 0:
            raise _PrismError(
                "precision_loss",
                "frame placement invalidates a roof polygon",
                component=index,
            )
        shifted_roof = roof + offset[2]
        for before, after in zip(_parts(polygon), _parts(shifted), strict=True):
            for old, new in zip(
                (before.exterior, *before.interiors),
                (after.exterior, *after.interiors),
                strict=True,
            ):
                old_xy, new_xy = np.asarray(old.coords), np.asarray(new.coords)
                if len(np.unique(new_xy, axis=0)) < len(np.unique(old_xy, axis=0)):
                    raise _PrismError(
                        "precision_loss",
                        "frame placement collapses a nonzero ring edge",
                        component=index,
                    )
                for z, new_z in (
                    (roof, shifted_roof),
                    (component.ground_z, component.ground_z + offset[2]),
                ):
                    original = np.column_stack([old_xy, np.full(len(old_xy), z)])
                    moved = np.column_stack([new_xy, np.full(len(new_xy), new_z)])
                    if (
                        not np.isfinite(moved).all()
                        or np.max(abs(np.spacing(moved))) >= tolerance / 4
                    ):
                        raise _PrismError(
                            "precision_loss",
                            "translated coordinates cannot resolve the requested tolerance",
                            component=index,
                        )
                    placement = max(
                        placement,
                        float(np.linalg.norm(moved - offset - original, axis=1).max()),
                    )
        levels.append((float(shifted_roof), shifted))
    movement = component.diagnostics["max_displacement"] + placement
    if movement > tolerance:
        raise _PrismError(
            "precision_loss",
            "frame placement exceeds the combined displacement budget",
            component=index,
        )
    ground = float(component.ground_z + offset[2])
    extent = max(
        component.base.bounds[2] - component.base.bounds[0],
        component.base.bounds[3] - component.base.bounds[1],
        max(z - component.ground_z for z, _ in component.levels),
        1.0,
    )
    epsilon = np.finfo(float).eps * extent * 64
    budget = component.diagnostics["volume_budget"] + _volume_budget(
        levels, ground, placement, epsilon
    )
    volume = float(sum(p.area * (z - ground) for z, p in levels))
    diagnostics = component.diagnostics
    if (
        abs(volume - diagnostics["source_shell_volume"]) > budget
        or abs(volume - diagnostics["shell_volume"]) > budget
    ):
        raise _PrismError(
            "precision_loss",
            "frame placement exceeds the volume budget",
            component=index,
        )
    diagnostics.update(
        max_displacement=movement,
        frame_placement_displacement=placement,
        prism_volume=volume,
        volume_residual=abs(volume - diagnostics["shell_volume"]),
        source_volume_residual=abs(volume - diagnostics["source_shell_volume"]),
        volume_budget=budget,
    )
    component.ground_z = ground
    component.base = translate(component.base, offset[0], offset[1])
    if not component.base.is_valid or component.base.area <= 0:
        raise _PrismError(
            "precision_loss", "frame placement invalidates the base", component=index
        )
    component.levels = tuple(levels)


def _convert(source, reserved, referenced, inherited_srs, tolerance, angle):
    records, transform, offset = _source_geometry(source, referenced, inherited_srs)
    if not np.isfinite(offset).all():
        raise _PrismError("unsupported_frame", "composed translation must be finite")
    components = _decompose_prisms(
        records, tolerance=tolerance, wall_angle_tolerance=angle
    )
    if any(offset):
        for index, component in enumerate(components):
            _placed_component(component, np.asarray(offset), tolerance, index)
        try:
            _reject_overlapping_components(components)
        except _PrismError as error:
            raise _PrismError(
                "precision_loss",
                "frame placement changes component contact topology",
                component=error.component,
            ) from error
    specifications = [
        (ci, li, component, roof, polygon)
        for ci, component in enumerate(components)
        for li, (roof, geometry) in enumerate(component.levels)
        for polygon in sorted(_parts(geometry), key=_geometry_key)
    ]
    ids = _piece_ids(source.id, len(specifications), reserved)
    pieces = []
    for id, (ci, li, component, roof, polygon) in zip(ids, specifications, strict=True):
        footprint = Surface().from_polygon(polygon, roof)
        footprint.transform = deepcopy(transform)
        piece = building_from_footprint(
            footprint,
            roof - component.ground_z,
            ground_height=component.ground_z,
            id=id,
            attributes=source.attributes,
        )
        piece.transform.srs = transform.srs
        piece.semantic_type = source.semantic_type
        piece.attributes.setdefault("source_id", source.id)
        piece.attributes["prism_decomposition"] = dict(
            version=1, input_id=source.id, component_index=ci, level_index=li
        )
        if any(offset):
            piece.attributes["prism_decomposition"]["input_translation"] = [
                float(value) for value in offset
            ]
        piece.calculate_bounds()
        pieces.append(piece)
    return pieces, [component.diagnostics for component in components]


def _has_lod1(source):
    return any(
        record.lod == "1"
        for obj in _objects([source])
        for record in obj.geometry.values()
    )


def _failure(source, error):
    return dict(
        input_id=source.id,
        reason=error.reason,
        message=str(error),
        component=error.component,
        face_index=error.face_index,
        location=error.location,
    )


def _city_output(target, buildings, report):
    result = copy(target)
    for name, value in vars(target).items():
        if name != "children":
            setattr(result, name, deepcopy(value))
    result.children = defaultdict(
        list,
        {
            key: deepcopy(group)
            for key, group in target.children.items()
            if key is not Building
        },
    )
    result.children[Building] = buildings
    if report["converted"]:
        result.attributes.pop("footprint_cleaning", None)
    result.attributes["footprints_from_lod1"] = report
    result.calculate_bounds()
    return result


def _city_conversion_error(target):
    if isinstance(target, City):
        try:
            _frame(target)
            _metre_srs(target.transform.srs)
            if (
                target.profile_id is not None
                or target.profile_version is not None
                or (target.schema_id, target.schema_version)
                not in ((None, None), (SCHEMA_ID, DEFAULT_VERSION))
                or target.semantic_type not in (None, SEMANTIC_NAMESPACE + "City")
            ):
                raise _PrismError(
                    "unsupported_metadata",
                    "City profile, schema or custom semantics require explicit transfer",
                )
        except _PrismError as error:
            return error
    return None


def _footprints_from_lod1(
    target, *, tolerance, wall_angle_tolerance, overwrite_lod0, on_error
):
    _check_parameters(tolerance, wall_angle_tolerance)
    if not isinstance(overwrite_lod0, bool):
        raise TypeError("overwrite_lod0 must be a bool")
    if on_error not in ("raise", "warn"):
        raise ValueError("on_error must be raise or warn")
    sources, reserved, referenced = _inputs(target)
    inherited_srs = target.transform.srs if isinstance(target, City) else ""
    city_error = _city_conversion_error(target)
    report = dict(
        version=1,
        parameters=dict(
            tolerance=float(tolerance),
            wall_angle_tolerance=float(wall_angle_tolerance),
            overwrite_lod0=overwrite_lod0,
            on_error=on_error,
        ),
        attempted=len(sources),
        converted=0,
        preserved=0,
        output_count=0,
        failures=[],
        residuals=[],
    )
    output = []
    for source in sources:
        try:
            if not _has_lod1(source):
                raise _PrismError("missing_geometry", "building has no LOD1 geometry")
            if not overwrite_lod0 and any(
                record.lod == "0" for record in source.geometry.values()
            ):
                output.append(deepcopy(source))
                report["preserved"] += 1
                continue
            if city_error is not None:
                raise city_error
            pieces, residuals = _convert(
                source,
                reserved,
                referenced,
                inherited_srs,
                tolerance,
                wall_angle_tolerance,
            )
        except _PrismError as error:
            failure = _failure(source, error)
            message = f"Building {source.id!r}: {error}; component={error.component}, face={error.face_index}, location={error.location}"
            if on_error == "raise":
                raise ValueError(message) from error
            warning(message)
            output.append(deepcopy(source))
            report["preserved"] += 1
            report["failures"].append(failure)
        else:
            output.extend(pieces)
            inherited_srs = inherited_srs or pieces[0].lod0.transform.srs
            report["converted"] += 1
            report["residuals"].append(dict(input_id=source.id, components=residuals))
    report["output_count"] = len(output)
    return _city_output(target, output, report) if isinstance(target, City) else output
