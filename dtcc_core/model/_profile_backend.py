"""Shared optional semantic-record backend, independent of native geometry.

Also used by the historical process-separated experiments. Optional dependencies
are imported only when a profile is explicitly loaded.
"""


def issue(node, slot, rule, message):
    parts = (slot,) if isinstance(slot, str) and slot else tuple(slot)
    suffix = "." + ".".join(map(str, parts)) if parts else ""
    return {
        "path": f"objects[{node['id']!r}]{suffix}",
        "rule": rule,
        "message": message,
        "node_id": node["id"],
        "slots": parts,
    }


def payload(node):
    """Admit the experiment's JSON projection at the process boundary."""
    if set(node) != {"id", "semantic_type", "attributes", "relations"}:
        raise ValueError("Expected id, semantic_type, attributes and relations")
    for key in ("id", "semantic_type"):
        if not isinstance(node[key], str) or not node[key]:
            raise ValueError(f"Expected nonempty string {key}")
    attrs, refs = node["attributes"], node["relations"]
    if not isinstance(attrs, dict) or not isinstance(refs, dict):
        raise ValueError("Expected attribute and relation dictionaries")
    if "id" in attrs or "id" in refs or attrs.keys() & refs.keys():
        raise ValueError(f"Overlapping projected fields on {node['id']}")
    return {"id": node["id"], **attrs, **refs}


def graph_checks(nodes, values, references):
    """Shared graph algorithm; each backend supplies its own declared target rules.

    Named references are ID lists; projected containment has one scalar parent.
    The LinkML backend reads range/inheritance and our containment annotation from
    SchemaView. The historical baseline supplies its own explicit rules. No domain
    names here.
    """
    errors, index = [], {}
    for node in nodes:
        if node["id"] in index:
            errors.append(
                issue(node, "id", "duplicate_id", "ID is not unique in this model")
            )
        index[node["id"]] = node
    if errors:
        return errors  # Ambiguous IDs must not resolve to an arbitrary occurrence.
    parents = {}
    for node, data in zip(nodes, values):
        for slot, (allowed, containment) in references.get(
            node["semantic_type"], {}
        ).items():
            targets = data.get(slot)
            if isinstance(targets, str):
                targets = [targets]
            if not isinstance(targets, list):
                continue  # Missing/malformed values belong to schema validation.
            for target_id in targets:
                if not isinstance(target_id, str):
                    continue
                target = index.get(target_id)
                if target is None:
                    errors.append(
                        issue(node, slot, "dangling_id", f"No object {target_id!r}")
                    )
                elif target["semantic_type"] not in allowed:
                    errors.append(
                        issue(
                            node,
                            slot,
                            "target_type",
                            f"{target_id!r} is {target['semantic_type']}; expected {sorted(allowed)}",
                        )
                    )
                if containment and target is not None:
                    parents[node["id"]] = target_id
    done = set()
    for start in parents:
        path = set()
        current = start
        while current in parents and current not in done:
            if current in path:
                errors.append(
                    issue(
                        index[current],
                        "parent",
                        "containment_cycle",
                        "Containment contains a cycle",
                    )
                )
                break
            path.add(current)
            current = parents[current]
        done.update(path)
    return errors


class LinkMLProfile:
    def __init__(self, path, *, closed=True, self_contained=False):
        from pathlib import Path

        from linkml.validator import Validator
        from linkml.validator.plugins import JsonschemaValidationPlugin
        from linkml_runtime.utils.schemaview import SchemaView

        if self_contained:
            import yaml
            from linkml_runtime.linkml_model import SchemaDefinition

            try:
                data = yaml.safe_load(Path(path).read_text())
            except yaml.YAMLError as exc:
                raise ValueError(f"Invalid YAML profile {path}: {exc}") from exc
            if not isinstance(data, dict) or data.get("imports"):
                raise ValueError(
                    "Select a self-contained local LinkML schema without imports"
                )
            if not data.get("id") or not data.get("version") or not data.get("classes"):
                raise ValueError("Profile requires id, version and classes")
            self.view = SchemaView(SchemaDefinition(**data))
        else:
            self.view = SchemaView(str(path))
        self.classes = self.view.all_classes()
        # Identifier-free classes are nested values, not nodes in the ID graph.
        self.record_classes = {
            name for name in self.classes if self.view.get_identifier_slot(name) is None
        }
        self.uris = {
            name: self.view.get_uri(cls, expand=True)
            for name, cls in self.classes.items()
        }
        self.class_names = {uri: name for name, uri in self.uris.items()}
        if len(self.class_names) != len(self.classes):
            raise ValueError("Profile class URIs must be unique")
        annotation = self.view.schema.annotations.get("dtcc_open_semantics")
        self.open_semantics = annotation is not None and annotation.value is True
        self.native_types = {}
        self.slots = {}
        for name, cls in self.classes.items():
            annotation = cls.annotations.get("dtcc_native_type")
            if annotation is not None:
                native = str(annotation.value)
                if native in self.native_types:
                    raise ValueError(f"Duplicate native type binding {native!r}")
                self.native_types[native] = self.uris[name]
            self.slots[self.uris[name]] = {
                slot.name for slot in self.view.class_induced_slots(name)
            }
        if self.open_semantics and not self.native_types:
            raise ValueError("Open semantics requires explicit native type bindings")
        # A record may name one representation on its owning Object. The schema
        # supplies the target's existing native binding; this is not an entity
        # reference or a second registry of domain-specific metadata names.
        record_references = {}
        for name in self.classes:
            for slot in self.view.class_induced_slots(name):
                annotation = getattr(
                    slot.annotations, "dtcc_representation_reference", None
                )
                if annotation is None:
                    continue
                target = self.classes.get(str(annotation.value))
                binding = target.annotations.get("dtcc_native_type") if target else None
                if (
                    name not in self.record_classes
                    or slot.range != "string"
                    or slot.multivalued
                    or slot.any_of
                    or binding is None
                ):
                    raise ValueError(
                        f"{name}.{slot.name}: representation reference requires a scalar "
                        "string in an inline record and a target class with a native binding"
                    )
                record_references.setdefault(name, []).append(
                    (slot.name, str(binding.value))
                )
        self.representation_references = {}
        object_class = self.class_names.get(self.native_types.get("Object"))
        object_classes = (
            set(self.view.class_descendants(object_class, reflexive=True))
            if object_class is not None
            else set()
        )
        for name in self.classes:
            for slot in self.view.class_induced_slots(name):
                ranges = [option.range for option in slot.any_of] or [slot.range]
                candidates = {
                    descendant
                    for target in ranges
                    if target in self.classes
                    for descendant in self.view.class_descendants(
                        target, reflexive=True
                    )
                }
                if not candidates.intersection(record_references):
                    continue
                if (
                    name not in object_classes
                    or not (slot.inlined or slot.inlined_as_list)
                    or slot.any_of
                    or slot.range not in record_references
                ):
                    raise ValueError(
                        f"{name}.{slot.name}: representation references require a direct "
                        "inline record slot on an Object-bound class or its descendants"
                    )
                self.representation_references.setdefault(self.uris[name], []).extend(
                    (slot.name, bool(slot.multivalued), field, native)
                    for field, native in record_references[slot.range]
                )
        self.validator = Validator(
            schema=self.view.schema,
            validation_plugins=[JsonschemaValidationPlugin(closed=closed)],
        )
        self.references = {}
        for kind in self.classes:
            refs = {}
            for slot in self.view.class_induced_slots(kind):
                ranges = [option.range for option in slot.any_of] or [slot.range]
                class_ranges = {
                    descendant
                    for target in ranges
                    if target in self.classes
                    for descendant in self.view.class_descendants(
                        target, reflexive=True
                    )
                }
                if class_ranges and (slot.inlined or slot.inlined_as_list):
                    if self_contained and any(
                        target not in self.record_classes for target in class_ranges
                    ):
                        raise ValueError(
                            f"{kind}.{slot.name}: inlined records must be identifier-free"
                        )
                    continue  # The existing JSON Schema plugin checks nested values.
                if self_contained and any(target in self.classes for target in ranges):
                    if not all(target in self.classes for target in ranges):
                        raise ValueError(
                            f"{kind}.{slot.name}: expected non-inlined class ID references"
                        )
                    if kind in self.record_classes or any(
                        target in self.record_classes for target in class_ranges
                    ):
                        raise ValueError(
                            f"{kind}.{slot.name}: records require explicit inlining and cannot contain ID references"
                        )
                if not all(target in self.classes for target in ranges):
                    continue
                allowed = {
                    self.uris[descendant]
                    for target in ranges
                    for descendant in self.view.class_descendants(
                        target, reflexive=True
                    )
                }
                annotation = getattr(slot.annotations, "containment", None)
                containment = annotation is not None and annotation.value is True
                if (
                    self_contained
                    and containment
                    and (slot.name != "parent" or slot.multivalued)
                ):
                    raise ValueError(
                        f"{kind}.{slot.name}: containment must be the scalar parent slot"
                    )
                refs[slot.name] = (allowed, containment)
            self.references[self.uris[kind]] = refs

    def validate(self, nodes, values):
        errors = []
        for node, data in zip(nodes, values):
            kind = self.class_names.get(node["semantic_type"])
            if (
                kind not in self.classes
                or self.classes[kind].abstract
                or kind in self.record_classes
            ):
                errors.append(
                    issue(
                        node,
                        "semantic_type",
                        "unknown_type",
                        f"Unsupported concrete semantic type {node['semantic_type']!r}",
                    )
                )
                continue
            # Open attributes are intentional in the buildings profile, but a
            # relationship must be declared or its target would escape checking.
            for slot in sorted(
                node["relations"].keys() - self.references[node["semantic_type"]].keys()
            ):
                errors.append(
                    issue(
                        node,
                        slot,
                        "unknown_relation",
                        "Relationship is not declared by this profile",
                    )
                )
            missing_paths = set()
            for result in self.validator.validate(data, kind).results:
                source = result.source
                slot = tuple(source.absolute_path)
                if source.validator == "required":
                    for missing in source.validator_value:
                        path = (*slot, missing)
                        if missing not in source.instance and path not in missing_paths:
                            missing_paths.add(path)
                            errors.append(
                                issue(
                                    node,
                                    path,
                                    "required",
                                    f"Missing required value {missing!r}",
                                )
                            )
                else:
                    errors.append(
                        issue(node, slot, str(source.validator), source.message)
                    )
        return errors, graph_checks(nodes, values, self.references)
