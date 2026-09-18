"""Explicit validation of native data against local LinkML schemas and profiles.

Importing this module does not import LinkML. Canonical I/O selects the bundled
standard schema; stricter domain profiles can also be evaluated explicitly.
See docs/design/standard-schema-io.md for the default contract and its scope.
"""

from dataclasses import dataclass
import json
from pathlib import Path

from . import exchange
from .object import Object, Tree, Landuse
from .geometry import Geometry, Bounds, Transform
from .values import Field, Raster
from ._profile_backend import LinkMLProfile, issue

__all__ = ['SemanticProfile', 'ValidationIssue', 'ValidationReport']


@dataclass(frozen=True)
class ValidationIssue:
    """A single problem found when validating data against a profile.

    Attributes
    ----------
    path : str
        Location of the problem in the validated data.
    rule : str
        Identifier of the rule that failed.
    message : str
        Human-readable description.
    """
    path: str
    rule: str
    message: str


@dataclass(frozen=True)
class ValidationReport:
    """Result of validating data against a semantic profile.

    ``valid`` is True when there are no issues.

    Attributes
    ----------
    profile_id : str
        Identifier of the profile used.
    profile_version : str
        Version of the profile used.
    issues : tuple[ValidationIssue, ...]
        Problems found.
    """
    profile_id: str
    profile_version: str
    issues: tuple[ValidationIssue, ...]

    @property
    def valid(self) -> bool:
        return not self.issues


def _project(model, references=None, *, native_types=None, slots=None,
             representation_references=None, representation_issues=None):
    """Project already-admitted native objects and retain paths to their authority."""
    def object_id(value):
        return json.dumps(['object', value])

    records, locations = [], {}

    def semantic_type(value, explicit=None):
        if native_types is None or explicit in slots:
            return explicit
        try:
            return native_types[type(value).__name__]
        except KeyError as exc:
            raise NotImplementedError(f'No schema binding for native {type(value).__name__}') from exc

    def append(id, semantic_type, attributes, relations, base, fields, intrinsic=None, reserved=()):
        intrinsic = intrinsic or {}
        if native_types is not None:
            # Unknown metadata/reference names remain in the native model. Their
            # generic integrity has already been admitted by Core; only declared
            # semantic slots enter the schema's flattened record.
            refs = (references or {}).get(semantic_type, {})
            relations = {key: value for key, value in relations.items() if key in refs}
            attributes = {key: value for key, value in attributes.items()
                          if key in slots[semantic_type] and key not in {'id', 'semantic_type', 'parent', *reserved}}
        for key in sorted(relations.keys() & ({'id', 'semantic_type'} | intrinsic.keys())):
            raise ValueError(f'{base}.relations[{key!r}]: conflicts with a projected native field')
        reserved = {'id', 'semantic_type', 'parent', *reserved} | intrinsic.keys() | relations.keys()
        for key in sorted(attributes.keys() & reserved):
            raise ValueError(f'{base}.attributes[{key!r}]: conflicts with a projected native field')
        records.append({'id': id, 'semantic_type': semantic_type,
                        'attributes': {**attributes, **intrinsic}, 'relations': relations})
        locations[id] = (base, fields)

    def project_field(value, base):
        names = ('name', 'unit', 'description', 'association')
        append(json.dumps(['field', base]), semantic_type(value), {}, {}, base,
               {name: f'{base}.{name}' for name in names},
               {name: getattr(value, name) for name in names})

    def project_geometry(value, base, owner_id, owner_path, region_key):
        # Only standalone roots and field owners need extra geometry records.
        # Plain nested polygons carry no separate semantic metadata to project.
        pending = [(value, base, owner_id, owner_path, region_key)]
        while pending:
            geometry, path, parent_id, parent_path, key = pending.pop()
            if native_types is not None:
                if parent_id is None or getattr(geometry, 'fields', None):
                    geometry_id = json.dumps(['geometry', path])
                    append(geometry_id, semantic_type(geometry), {}, {}, path, {})
                    if parent_id is None:
                        parent_id, parent_path = geometry_id, path
                for index, value in enumerate(getattr(geometry, 'fields', ())):
                    project_field(value, f'{path}.fields[{index}]')
            for index, region in enumerate(getattr(geometry, 'regions', ())):
                region_base = f'{path}.regions[{index}]'
                kind = semantic_type(region, region.semantic_type)
                extra_refs = (references or {}).get(kind, {}).keys() - {'parent', 'host'}
                if extra_refs:
                    raise NotImplementedError(f'{region_base}: region references support only owner parent and host')
                relations = {'parent': parent_id}
                if region.parent is not None:
                    relations['host'] = json.dumps([*key, region.parent])
                fields = {'id': region_base, 'semantic_type': f'{region_base}.semantic_type',
                          'parent': parent_path, 'host': f'{region_base}.parent'}
                append(json.dumps([*key, index]), kind, region.attributes, relations,
                       region_base, fields, reserved=('host',))
            if native_types is not None:
                child_key = 'linestrings' if hasattr(geometry, 'linestrings') else 'surfaces'
                for index, surface in enumerate(getattr(geometry, child_key, ())):
                    if surface.fields or surface.regions:
                        child_path = f'{path}.{child_key}[{index}]'
                        pending.append((surface, child_path, parent_id, parent_path, ['region', child_path]))

    if isinstance(model, Field):
        project_field(model, 'field')
        return records, locations
    if isinstance(model, (Geometry, Raster, Bounds, Transform)):
        project_geometry(model, 'geometry', None, None, ['region', 'geometry'])
        return records, locations

    stack = [(model, None, None)]
    while stack:
        obj, parent, containment_path = stack.pop()
        base = f'objects[{obj.id!r}]'
        kind = semantic_type(obj, obj.semantic_type)
        relations = {name: [object_id(id) for id in ids] for name, ids in obj.relations.items()}
        fields = {name: f'{base}.relations[{name!r}]' for name in relations}
        fields.update(id=f'{base}.id', semantic_type=f'{base}.semantic_type', parent=containment_path or base)
        for slot in (references or {}).get(kind, {}):
            fields.setdefault(slot, f'{base}.relations[{slot!r}]')
            if slot in obj.attributes and (native_types is None or slot != 'parent'):
                raise ValueError(f'{base}.attributes[{slot!r}]: class references belong in Object.relations')
        if parent is not None:
            relations['parent'] = object_id(parent.id)
        intrinsic = {}
        if type(obj) is Tree:
            intrinsic = {'height': obj.height, 'crown_radius': obj.crown_radius}
            fields.update(height=f'{base}.height', crown_radius=f'{base}.crown_radius')
        elif type(obj) is Landuse:
            intrinsic = {'native_landuses': [code.name for code in obj.landuses]}
            fields['native_landuses'] = f'{base}.landuses'
        append(object_id(obj.id), kind, obj.attributes, relations, base, fields, intrinsic)
        for slot, multiple, field, native in (representation_references or {}).get(kind, ()):
            value = obj.attributes.get(slot)
            if multiple:
                values = enumerate(value) if isinstance(value, list) else ()
            else:
                values = [(None, value)]
            for index, record in values:
                target_id = record.get(field) if isinstance(record, dict) else None
                if not isinstance(target_id, str) or not target_id.strip():
                    continue  # LinkML checks missing/malformed record values.
                target = obj.geometry.get(target_id)
                path = (slot, index, field) if multiple else (slot, field)
                if target is None:
                    error = issue(records[-1], path, 'dangling_representation',
                                  f'No representation {target_id!r} on this object')
                elif type(target.geometry).__name__ != native:
                    error = issue(records[-1], path, 'representation_type',
                                  f'Representation {target_id!r} must contain {native}')
                else:
                    continue
                if representation_issues is not None:
                    representation_issues.append(error)
        for key, record in obj.geometry.items():
            project_geometry(record.geometry, f'{base}.geometry[{key!r}].geometry',
                             object_id(obj.id), base, ['region', obj.id, key])
        children = [(child, obj, f'{base}.children[{cls.__name__}][{index}]')
                    for cls, group in obj.children.items() for index, child in enumerate(group)]
        stack.extend(reversed(children))
    return records, locations


class SemanticProfile:
    """Load one self-contained local schema and reuse it for explicit validation.

    Core's LinkML dependencies must be installed in the calling Python environment.
    Loading never follows imports or fetches the profile URI. Edits require a new
    instance. Validation uses this selected schema regardless of stored profile
    labels and never changes the model. Unknown attributes are allowed; named
    relationships must be declared unless the schema permits open semantics.
    """

    def __init__(self, path: str | Path):
        path = Path(path)
        if not path.is_file():
            raise FileNotFoundError(f'Local profile schema not found: {path}')
        try:
            self._backend = LinkMLProfile(path, closed=False, self_contained=True)
        except ModuleNotFoundError as exc:
            raise ImportError(
                'SemanticProfile requires LinkML validation dependencies in this Python environment; '
                'install the current dtcc-core dependencies. See docs/design/standard-schema-io.md'
            ) from exc
        self.profile_id = str(self._backend.view.schema.id)
        self.profile_version = str(self._backend.view.schema.version)

    def validate(self, model: Object) -> ValidationReport:
        """Return semantic issues; raise for invalid/unsupported native state.

        Reuses canonical admission, including numerical checks, before projection.
        A successful report is a snapshot, not continuous enforcement on mutable data.
        Diagnostic ``objects[id]`` identifies a model-scoped object, not a new API.
        Schemas declaring native bindings can also validate Geometry and Field roots.
        """
        if not isinstance(model, Object) and not self._backend.open_semantics:
            raise TypeError('SemanticProfile.validate expects a native Object root')
        exchange.validate(model)
        return self._validate_admitted(model)

    def _validate_admitted(self, model):
        """Semantic stage within an operation that already admitted native state."""
        native_types = self._backend.native_types if self._backend.open_semantics else None
        representation_issues = []
        records, locations = _project(model, self._backend.references,
                                     native_types=native_types, slots=self._backend.slots,
                                     representation_references=self._backend.representation_references,
                                     representation_issues=representation_issues)
        values = [{'id': node['id'], **node['attributes'], **node['relations']} for node in records]
        schema_errors, graph_errors = self._backend.validate(records, values)
        issues = []
        for error in (*schema_errors, *graph_errors, *representation_issues):
            base, fields = locations[error['node_id']]
            slots = error['slots']
            path = base
            if slots:
                path = fields.get(slots[0], f'{base}.attributes[{slots[0]!r}]')
                path += ''.join(f'[{part!r}]' for part in slots[1:])
            issues.append(ValidationIssue(path, error['rule'], error['message']))
        return ValidationReport(self.profile_id, self.profile_version, tuple(issues))
