"""Compare two optional LinkML profiles with an intentionally frozen Python baseline.

Experiment only: per-object schema validation plus a small ID/reference graph pass.
Run in the isolated LinkML environment; it does not import dtcc_core.
"""

import argparse
from copy import deepcopy
from importlib.metadata import version
import json
from pathlib import Path
import sys



HERE = Path(__file__).resolve().parent
NAMESPACE = "https://example.org/dtcc/"


from linkml_profile import LinkMLProfile, graph_checks, issue, payload


# Deliberately independent duplication for comparison ONLY. This baseline supports
# the original four types and the low-rise height rule, not schema-added classes.
BASE_FIELDS = {
    "City": {"id", "name"},
    "Building": {"id", "parent", "height", "usage"},
    "BuildingPart": {"id", "parent"},
    "Sensor": {"id", "parent", "observes", "phenomenon"},
}
BASE_REFERENCES = {
    "Building": {"parent": ({"City"}, True)},
    "BuildingPart": {"parent": ({"Building", "BuildingPart"}, True)},
    "Sensor": {"parent": ({"City"}, True), "observes": ({"Building"}, False)},
}

BASE_REFERENCES = {
    NAMESPACE + kind: {slot: ({NAMESPACE + target for target in allowed}, containment)
                       for slot, (allowed, containment) in refs.items()}
    for kind, refs in BASE_REFERENCES.items()
}
BASE_TYPE_NAMES = {NAMESPACE + kind: kind for kind in BASE_FIELDS}


def python_baseline(nodes, values, profile):
    errors = []
    for node, data in zip(nodes, values):
        kind = BASE_TYPE_NAMES.get(node["semantic_type"])
        if kind not in BASE_FIELDS:
            errors.append(issue(
                node, "semantic_type", "unknown_type", f"Unsupported concrete type {kind!r}",
            ))
            continue
        required = BASE_FIELDS[kind]
        for slot in sorted(required - data.keys()):
            errors.append(issue(node, slot, "required", "Missing required value"))
        for slot in sorted(data.keys() - required):
            errors.append(issue(node, slot, "additionalProperties", "Unknown attribute/relation"))
        for slot in sorted(required & data.keys()):
            value = data[slot]
            if slot == "height":
                if isinstance(value, bool) or not isinstance(value, (int, float)):
                    errors.append(issue(node, slot, "type", "Expected a number"))
                elif value < 0 or (profile == "low-rise" and value > 10):
                    errors.append(issue(node, slot, "range", "Height outside profile range"))
            elif slot == "observes":
                if not isinstance(value, list) or len(value) != 1 or not isinstance(value[0], str):
                    errors.append(issue(node, slot, "cardinality", "Expected exactly one ID in a list"))
            elif not isinstance(value, str):
                errors.append(issue(node, slot, "type", "Expected a string/scalar ID"))
            elif slot == "usage" and value not in {"residential", "public"}:
                errors.append(issue(node, slot, "enum", "Unknown building usage"))
            elif slot == "phenomenon" and value != "temperature":
                errors.append(issue(node, slot, "const", "Expected temperature"))
    return errors + graph_checks(nodes, values, BASE_REFERENCES)


def cases(example):
    """Focused, named acceptance cases; no Cartesian matrix or synthetic scale test."""
    original = example["objects"]
    yield "valid", original, "city", True, True
    mutations = [
        ("missing_height", 1, "attributes", "height", None),
        ("negative_height", 1, "attributes", "height", -1),
        ("height_wrong_type", 1, "attributes", "height", "12.5"),
        ("invalid_usage", 1, "attributes", "usage", "unknown"),
        ("unknown_attribute", 1, "attributes", "typo", True),
        ("missing_observes", 3, "relations", "observes", None),
        ("too_many_targets", 3, "relations", "observes", ["building-1", "part-1"]),
        ("wrong_target_type", 3, "relations", "observes", ["part-1"]),
        ("dangling_id", 3, "relations", "observes", ["missing"]),
        ("wrong_parent_type", 2, "relations", "parent", "city-1"),
        ("containment_cycle", 2, "relations", "parent", "part-1"),
    ]
    for name, index, section, key, value in mutations:
        nodes = deepcopy(original)
        if value is None:
            del nodes[index][section][key]
        else:
            nodes[index][section][key] = value
        yield name, nodes, "city", False, False
    duplicate = deepcopy(original)
    duplicate.append(deepcopy(original[1]))
    yield "duplicate_id", duplicate, "city", False, False
    unknown = deepcopy(original)
    unknown[3]["semantic_type"] = NAMESPACE + "Unrecognized"
    yield "unknown_type", unknown, "city", False, False
    yield "low_rise_rejects_tall", original, "low-rise", False, False
    low = deepcopy(original)
    low[1]["attributes"]["height"] = 8.0
    yield "low_rise_accepts_short", low, "low-rise", True, True
    yield "schema_added_type", example["extension_objects"], "low-rise", True, False
    extension = deepcopy(example["extension_objects"])
    del extension[3]["attributes"]["accuracy"]
    yield "schema_added_required_attribute", extension, "low-rise", False, False


def read_example(path):
    def reject_constant(value):
        raise ValueError(f"Non-JSON number: {value}")
    example = json.loads(path.read_text(), parse_constant=reject_constant)
    for key in ("objects", "extension_objects"):
        if not isinstance(example.get(key), list) or not example[key]:
            raise ValueError(f"Expected nonempty {key} list")
        for node in example[key]:
            payload(node)
    return example


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("example", type=Path)
    parser.add_argument("--output", type=Path, help="Write the complete comparison evidence")
    parser.add_argument("--case", help="Validate one named case; exit 1 when invalid")
    args = parser.parse_args()
    example = read_example(args.example)
    profiles = {name: LinkMLProfile(HERE / f"{name}.yaml") for name in ("city", "low-rise")}
    results = []
    for name, nodes, profile, expected_linkml, expected_python in cases(example):
        if args.case and args.case != name:
            continue
        before = deepcopy(nodes)
        values = [payload(node) for node in nodes]
        schema_errors, graph_errors = profiles[profile].validate(nodes, values)
        python_errors = python_baseline(nodes, values, profile)
        valid = not (schema_errors or graph_errors)
        python_valid = not python_errors
        assert nodes == before, f"Validation mutated the projection: {name}"
        results.append({
            "case": name, "profile": profile,
            "schema_valid": not schema_errors, "linkml_plus_graph_valid": valid,
            "python_valid": python_valid,
            "expected_outcomes_met": valid == expected_linkml and python_valid == expected_python,
            "schema_errors": schema_errors, "graph_errors": graph_errors,
            "python_errors": python_errors,
        })
        print(
            f"{name}: schema={'PASS' if not schema_errors else 'FAIL'}, "
            f"schema+graph={'PASS' if valid else 'FAIL'}, "
            f"python={'PASS' if python_valid else 'FAIL'}"
        )
    if not results:
        parser.error(f"Unknown case: {args.case}")
    if args.case:
        for error in results[0]["schema_errors"] + results[0]["graph_errors"]:
            print(f"{error['path']}: {error['message']}")
        return 0 if results[0]["linkml_plus_graph_valid"] else 1
    # Reuse both instances in one process and return to A after B.
    values = [payload(node) for node in example["objects"]]
    coexistence = [not any(profiles[name].validate(example["objects"], values))
                   for name in ("city", "low-rise", "city")]
    assert coexistence == [True, False, True], coexistence
    evidence = {
        "versions": {name: version(name) for name in ("linkml", "linkml-runtime", "jsonschema")},
        "profile_ids": {name: {"id": str(p.view.schema.id), "version": p.view.schema.version}
                        for name, p in profiles.items()},
        "coexistence_A_B_A": coexistence,
        "canonical_exchange_observation": example["canonical_exchange_observation"],
        "cases": results,
    }
    if args.output:
        args.output.write_text(json.dumps(evidence, indent=2, allow_nan=False) + "\n")
    if not all(result["expected_outcomes_met"] for result in results):
        print("Unexpected comparison outcome; inspect evidence.", file=sys.stderr)
        return 1
    print(f"All {len(results)} expected outcomes matched; profiles coexist without input mutation.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
